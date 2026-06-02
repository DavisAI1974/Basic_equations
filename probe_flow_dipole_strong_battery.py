"""
STRONG-force two-particle BATTERY -- real CERN open data, many extraction tools.

CONTEXT / honesty up front:
  The dipole tool (windowed-H -> MI -> dMI/d(axis)) gave LOW R^2 on gravity and
  weak.  Greg's "first try, not only try": run a BATTERY of DIFFERENT extraction
  tools in quick succession on REAL strong-force two-particle data and see if ANY
  tool surfaces genuine 2-channel structure.  Each tool = one data point, a metric
  vs an appropriate null/control, and a boolean "hits".  NO new-law claims.

DATA -- what we could actually get:
  There is NO per-jet / dijet / hadron-pair CSV in the CMS open-data EDUCATION
  set (queried the opendata.cern.ch API for dijet/jet/two-jet csv: zero CSV hits;
  only dimuon / dielectron / W-lepton CSVs exist; record/700 MuRun is dimuon too).
  The closest REAL strong-interaction two-particle sample available as a small CSV
  is the broad CMS dimuon set (record/303, dimuon.csv, 100k events) whose mass
  spectrum is dominated by the J/psi and Upsilon -- charmonium / bottomonium, i.e.
  QCD (strong-force) BOUND STATES -- plus the (electroweak) Z.  We use it honestly:
    * channels A,B = leading / subleading muon of the pair (a real correlated pair,
      they share the parent's 4-momentum -- no invented coupling).
    * STRONG-relevant flow axes: pair invariant mass M (J/psi 3.097, Upsilon ~9.46
      are strong-force resonances) AND the femtoscopy momentum difference
      q = |p_A - p_B| toward q->0.
    * CONTROL: event-MIXED pairs (muon A from one event, muon B from another) and
      same-/opposite-charge splits.  Genuine structure must show in real pairs and
      be absent in the mixed/scrambled control.
  Caveat stated in the JSON: muons are leptons; the *pair* is QCD-produced for the
  J/psi/Upsilon component but the per-muon kinematics are not hadronic tracks.  A
  true pion/kaon femtoscopy sample would be better; it is not available as a CSV.

TOOLS (each one block, keep moving):
  1. Baseline dipole, multiple channels (pT / energy / rapidity) x axes (M / q),
     equal-EVENT quantile bins (n_axis_points >> n_params).
  2. Estimator swap: histogram-MI vs KSG kNN-MI; vary window/bin counts.
  3. RAW-signal dipole (regress channel-B raw moment on channel-A raw moment along
     axis) vs entropy-operator dipole.
  4. KNOWN strong-force 2-particle observable: Bose-Einstein / femtoscopy
     C(q)=N_same(q)/N_mixed(q) for LIKE-SIGN pairs.  CRITICAL diagnostic: if this
     real correlation is visible but the dipole tools are not -> "dipole tool blind
     to a real correlation present here."
  5. PySR symbolic regression (if available) / SINDy poly fallback for any governing
     relation between the channel operators along the axis.
  6. Direct coupling references: Pearson/Spearman cross-correlation and a transfer-
     entropy-style binned conditional-MI, real vs mixed control.

Run:  python probe_flow_dipole_strong_battery.py [--canary]
"""
import json
import time
import argparse
import numpy as np
import pandas as pd

DATA = "data/strong/dimuon.csv"
JPSI = "data/strong/dimuon-Jpsi.csv"


# ----------------------------------------------------------------------------- IO
def load_pairs(path=DATA, max_rows=None):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if max_rows:
        df = df.iloc[:max_rows]
    E1, px1, py1, pz1 = df.E1.values, df.px1.values, df.py1.values, df.pz1.values
    E2, px2, py2, pz2 = df.E2.values, df.px2.values, df.py2.values, df.pz2.values
    pt1, pt2 = df.pt1.values, df.pt2.values
    eta1, eta2 = df.eta1.values, df.eta2.values
    Q1, Q2 = df.Q1.values, df.Q2.values
    M = df.M.values
    # rapidity y = 0.5 ln((E+pz)/(E-pz))
    def rap(E, pz):
        num = np.clip(E + pz, 1e-9, None)
        den = np.clip(E - pz, 1e-9, None)
        return 0.5 * np.log(num / den)
    y1, y2 = rap(E1, pz1), rap(E2, pz2)
    # q = |p_A - p_B| (3-momentum difference) -- femtoscopy axis
    q = np.sqrt((px1 - px2) ** 2 + (py1 - py2) ** 2 + (pz1 - pz2) ** 2)
    lead = pt1 >= pt2
    def pick(a1, a2):
        return np.where(lead, a1, a2), np.where(lead, a2, a1)
    pt_l, pt_s = pick(pt1, pt2)
    E_l, E_s = pick(E1, E2)
    y_l, y_s = pick(y1, y2)
    return dict(M=M, q=q, Q1=Q1, Q2=Q2,
                pt_l=pt_l, pt_s=pt_s, E_l=E_l, E_s=E_s, y_l=y_l, y_s=y_s,
                px1=px1, py1=py1, pz1=pz1, E1=E1,
                px2=px2, py2=py2, pz2=pz2, E2=E2,
                pt1=pt1, pt2=pt2)


# ----------------------------------------------------------------- estimators
def entropy_1d(x, bins=12):
    h, e = np.histogram(x, bins=bins, density=True)
    w = e[1] - e[0]
    p = h * w
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def mi_hist(x, y, bins=8):
    c, _, _ = np.histogram2d(x, y, bins=bins)
    pxy = c / c.sum()
    px = pxy.sum(1, keepdims=True)
    py = pxy.sum(0, keepdims=True)
    nz = pxy > 0
    return float(np.sum(pxy[nz] * np.log(pxy[nz] / (px * py)[nz])))


def mi_ksg(x, y, k=4):
    """Kraskov-Stoegbauer-Grassberger kNN MI (estimator 1). Pure numpy."""
    from scipy.spatial import cKDTree
    from scipy.special import digamma
    n = len(x)
    if n < k + 2:
        return 0.0
    xr = (x - x.mean()) / (x.std() + 1e-12)
    yr = (y - y.mean()) / (y.std() + 1e-12)
    z = np.column_stack([xr, yr])
    tree = cKDTree(z)
    # distance to k-th neighbour (Chebyshev)
    d, _ = tree.query(z, k=k + 1, p=np.inf)
    eps = d[:, k]
    tx = cKDTree(xr[:, None])
    ty = cKDTree(yr[:, None])
    nx = np.array([len(tx.query_ball_point([xr[i]], eps[i] - 1e-12, p=np.inf)) - 1
                   for i in range(n)])
    ny = np.array([len(ty.query_ball_point([yr[i]], eps[i] - 1e-12, p=np.inf)) - 1
                   for i in range(n)])
    nx = np.clip(nx, 1, None)
    ny = np.clip(ny, 1, None)
    mi = digamma(k) + digamma(n) - np.mean(digamma(nx + 1) + digamma(ny + 1))
    return float(max(mi, 0.0))


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ssr = float(np.sum((y - yhat) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ssr / sst if sst > 0 else float("nan")
    return beta, r2, ssr


def quantile_axis(axis, chans, ev_per_bin=200, min_bin=80,
                  estimator="hist", h_bins=12, mi_bins=8, k=4):
    """Equal-EVENT quantile bins along `axis`; per bin compute H_a,H_b,MI of the
    two channels.  chans=(a_arr,b_arr).  Returns (axis_center,Ha,Hb,MI)."""
    a, b = chans
    order = np.argsort(axis)
    ax, a, b = axis[order], a[order], b[order]
    n = len(ax)
    nb = max(8, n // ev_per_bin)
    edges = np.linspace(0, n, nb + 1).astype(int)
    Ac, Ha, Hb, MI = [], [], [], []
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        if hi - lo < min_bin:
            continue
        Ac.append(float(np.median(ax[lo:hi])))
        Ha.append(entropy_1d(a[lo:hi], bins=h_bins))
        Hb.append(entropy_1d(b[lo:hi], bins=h_bins))
        if estimator == "ksg":
            MI.append(mi_ksg(a[lo:hi], b[lo:hi], k=k))
        else:
            MI.append(mi_hist(a[lo:hi], b[lo:hi], bins=mi_bins))
    return (np.array(Ac), np.array(Ha), np.array(Hb), np.array(MI))


def fit_dipole(Ac, Ha, Hb, MI):
    """dMI/d(axis) ~ paper dipole library; report full R2, cross-gain, shuffle null."""
    if len(Ac) < 8:
        return None
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    dMI = np.gradient(MIs, Ac)
    sl = slice(1, len(Ac) - 1)
    Ha_, Hb_, MI_, dMI_ = Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha_ ** 2, Hb_ ** 2, Ha_ * Hb_
    one = np.ones_like(Ha_)
    n = len(dMI_)
    Xfull = np.column_stack([Ha2, Hb2, HaHb, Ha_, Hb_, one])
    Xself = np.column_stack([Ha2, Hb2, Ha_, Hb_, one])
    _, r2_full, _ = ols(Xfull, dMI_)
    _, r2_self, _ = ols(Xself, dMI_)
    rng = np.random.default_rng(0)
    nullr2 = []
    for _ in range(40):
        p = rng.permutation(n)
        Xn = np.column_stack([Ha2[p], Hb2[p], HaHb[p], Ha_[p], Hb_[p], one])
        _, r2n, _ = ols(Xn, dMI_)
        nullr2.append(r2n)
    return {"n_axis_points": int(n), "R2_full": float(r2_full),
            "dR2_cross": float(r2_full - r2_self),
            "shuffle_null_R2_mean": float(np.mean(nullr2)),
            "shuffle_null_R2_p95": float(np.percentile(nullr2, 95)),
            "MI_range": [float(MI_.min()), float(MI_.max())]}


# ============================================================= TOOL 1: baseline dipole
def tool1_baseline_dipole(d, sel, ev_per_bin, results):
    chan_sets = {"pt": (d["pt_l"][sel], d["pt_s"][sel]),
                 "energy": (d["E_l"][sel], d["E_s"][sel]),
                 "rapidity": (d["y_l"][sel], d["y_s"][sel])}
    axes = {"mass_M": d["M"][sel], "femto_q": d["q"][sel]}
    best_hit = False
    detail = {}
    for an, ax in axes.items():
        for cn, ch in chan_sets.items():
            Ac, Ha, Hb, MI = quantile_axis(ax, ch, ev_per_bin=ev_per_bin)
            fd = fit_dipole(Ac, Ha, Hb, MI)
            if fd is None:
                continue
            hit = (fd["R2_full"] > 0.5) and (fd["R2_full"] >
                                             fd["shuffle_null_R2_p95"] + 0.1)
            detail[f"{an}|{cn}"] = {"R2_full": fd["R2_full"],
                                    "dR2_cross": fd["dR2_cross"],
                                    "null_p95": fd["shuffle_null_R2_p95"],
                                    "n_axis": fd["n_axis_points"], "hits": hit}
            best_hit = best_hit or hit
    bestkey = max(detail, key=lambda k: detail[k]["R2_full"]) if detail else None
    results.append({
        "tool": "1_baseline_dipole (chan x axis grid)",
        "metric": f"max R2_full = {detail[bestkey]['R2_full']:.3f} @ {bestkey}"
                  if bestkey else "no stable axis",
        "control": f"shuffle-null R2 p95 = {detail[bestkey]['null_p95']:.3f}"
                   if bestkey else "n/a",
        "hits": bool(best_hit), "detail": detail})


# ============================================================= TOOL 2: estimator swap
def tool2_estimator_swap(d, sel, ev_per_bin, results):
    ax = d["M"][sel]
    ch = (d["pt_l"][sel], d["pt_s"][sel])
    rows = {}
    for est in ["hist", "ksg"]:
        for mb in [6, 10]:
            Ac, Ha, Hb, MI = quantile_axis(ax, ch, ev_per_bin=ev_per_bin,
                                           estimator=est, mi_bins=mb, k=4)
            fd = fit_dipole(Ac, Ha, Hb, MI)
            if fd is None:
                continue
            rows[f"{est}_b{mb}"] = {"R2_full": fd["R2_full"],
                                    "null_p95": fd["shuffle_null_R2_p95"],
                                    "hits": fd["R2_full"] >
                                    fd["shuffle_null_R2_p95"] + 0.1}
    any_hit = any(r["hits"] for r in rows.values()) if rows else False
    bk = max(rows, key=lambda k: rows[k]["R2_full"]) if rows else None
    results.append({
        "tool": "2_estimator_swap (hist vs KSG kNN MI, axis=M, chan=pt)",
        "metric": f"max R2_full = {rows[bk]['R2_full']:.3f} @ {bk}" if bk else "n/a",
        "control": f"null p95 = {rows[bk]['null_p95']:.3f}" if bk else "n/a",
        "hits": bool(any_hit), "detail": rows})


# ============================================================= TOOL 3: raw-signal dipole
def tool3_raw_dipole(d, sel, ev_per_bin, results):
    """Raw moments per bin (mean of channel) regressed instead of entropy operators.
    Compares whether RAW per-bin channel statistics carry a stronger axis relation
    than entropy operators do (entropy-dipole R2_full from tool1 mass|pt)."""
    ax = d["M"][sel]
    a, b = d["pt_l"][sel], d["pt_s"][sel]
    order = np.argsort(ax)
    ax2, a2, b2 = ax[order], a[order], b[order]
    n = len(ax2)
    nb = max(8, n // ev_per_bin)
    edges = np.linspace(0, n, nb + 1).astype(int)
    Ac, mA, mB, sA, sB = [], [], [], [], []
    for i in range(nb):
        lo, hi = edges[i], edges[i + 1]
        if hi - lo < 80:
            continue
        Ac.append(float(np.median(ax2[lo:hi])))
        mA.append(float(a2[lo:hi].mean())); mB.append(float(b2[lo:hi].mean()))
        sA.append(float(a2[lo:hi].std())); sB.append(float(b2[lo:hi].std()))
    Ac = np.array(Ac); mA = np.array(mA); mB = np.array(mB)
    sA = np.array(sA); sB = np.array(sB)
    if len(Ac) < 8:
        results.append({"tool": "3_raw_signal_dipole", "metric": "n/a",
                        "control": "n/a", "hits": False})
        return
    # raw cross-channel: does mean-B track mean-A along the mass axis?
    _, r2_raw, _ = ols(np.column_stack([mA, np.ones_like(mA)]), mB)
    # control: shuffle B vs A
    rng = np.random.default_rng(1)
    nullr2 = [ols(np.column_stack([mA, np.ones_like(mA)]),
                  mB[rng.permutation(len(mB))])[1] for _ in range(40)]
    nullp95 = float(np.percentile(nullr2, 95))
    hit = (r2_raw > 0.5) and (r2_raw > nullp95 + 0.1)
    results.append({
        "tool": "3_raw_signal_dipole (mean pt_sub ~ mean pt_lead along M)",
        "metric": f"raw cross R2 = {r2_raw:.3f}",
        "control": f"shuffle-null p95 = {nullp95:.3f}",
        "hits": bool(hit),
        "detail": {"r2_raw_meanB_vs_meanA": float(r2_raw),
                   "null_p95": nullp95}})


# ============================================ TOOL 4: Bose-Einstein / femtoscopy C(q)
def tool4_femtoscopy(d, sel, results, n_mix=None):
    """C(q) = N_same(q) / N_mixed(q) for LIKE-SIGN pairs.
    Same-event = the two muons as recorded.  Mixed = muon A of event i paired with
    muon B of a DIFFERENT random event (event mixing).  A real BE correlation rises
    toward q->0 for same-event; mixed is flat by construction.  CRITICAL diagnostic
    that the kinematics carry a known real 2-particle correlation."""
    Q1, Q2 = d["Q1"][sel], d["Q2"][sel]
    like = (Q1 * Q2) > 0
    # same-event q for like-sign pairs
    q_same = d["q"][sel][like]
    # build mixed (event-mixed) like-sign pairs: muon1 of one event + muon2 of another
    idx = np.where(sel)[0]
    p1 = np.column_stack([d["px1"][idx], d["py1"][idx], d["pz1"][idx]])
    p2 = np.column_stack([d["px2"][idx], d["py2"][idx], d["pz2"][idx]])
    Q1a, Q2a = d["Q1"][idx], d["Q2"][idx]
    rng = np.random.default_rng(42)
    m = len(idx) if n_mix is None else min(n_mix, len(idx))
    perm = rng.permutation(len(idx))[:m]
    perm2 = rng.permutation(len(idx))[:m]
    # mixed pair: muon-1 from event perm[i], muon-2 from event perm2[i]
    qmix = np.sqrt(np.sum((p1[perm] - p2[perm2]) ** 2, axis=1))
    qmix_like = qmix[(Q1a[perm] * Q2a[perm2]) > 0]
    # ratio in q bins (low-q region)
    qmax = float(np.percentile(np.concatenate([q_same, qmix_like]), 60))
    bins = np.linspace(0, max(qmax, 1.0), 13)
    hs, _ = np.histogram(q_same, bins=bins)
    hm, _ = np.histogram(qmix_like, bins=bins)
    hs = hs / max(hs.sum(), 1)
    hm = hm / max(hm.sum(), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.where(hm > 0, hs / hm, np.nan)
    centers = 0.5 * (bins[:-1] + bins[1:])
    # is C(q) enhanced at low q vs high q?  compare first vs last valid thirds
    valid = ~np.isnan(C)
    Cv, cv = C[valid], centers[valid]
    if len(Cv) >= 4:
        lowC = float(np.nanmean(Cv[:max(2, len(Cv) // 3)]))
        highC = float(np.nanmean(Cv[-max(2, len(Cv) // 3):]))
        ratio = lowC / highC if highC > 0 else float("nan")
    else:
        lowC = highC = ratio = float("nan")
    # Visibility diagnostic = is ANY real q-structure present vs the flat event-mixed
    # control?  A genuine same-event/mixed ratio departs from 1 in q (BE peak would be
    # low-q ENHANCEMENT, ratio>1; resonance/back-to-back kinematics give low-q
    # SUPPRESSION, ratio<1).  Either direction means the kinematics carry a real
    # q-correlation the tools could in principle see; report the direction explicitly.
    rms_dev = float(np.sqrt(np.nanmean((Cv - 1.0) ** 2))) if len(Cv) else float("nan")
    direction = ("low-q ENHANCEMENT (BE-like)" if (np.isfinite(ratio) and ratio > 1.0)
                 else "low-q SUPPRESSION (resonance/back-to-back kinematics)")
    hit = bool(np.isfinite(ratio) and (ratio > 1.15 or ratio < 0.87))  # structure either way
    results.append({
        "tool": "4_femtoscopy_C(q)_like-sign (KNOWN 2-particle q-correlation)",
        "metric": f"C(low q)/C(high q) = {ratio:.3f} ({direction}); RMS|C-1|={rms_dev:.3f}",
        "control": "denominator = event-MIXED like-sign pairs (flat by construction)",
        "hits": hit,
        "detail": {"ratio_low_over_high": float(ratio), "direction": direction,
                   "rms_dev_from_flat": rms_dev,
                   "q_centers": [round(x, 3) for x in cv.tolist()],
                   "C_of_q": [round(float(x), 3) for x in Cv.tolist()],
                   "n_same_like": int(len(q_same)),
                   "n_mixed_like": int(len(qmix_like))}})


# ============================================ TOOL 5: PySR / SINDy symbolic relation
def tool5_symbolic(d, sel, ev_per_bin, results, try_pysr=True):
    """Symbolic relation between the per-bin channel operators along M.
    Target: MI ~ f(H_a, H_b).  PySR if available; else SINDy-style poly + report
    best polynomial R2 vs a shuffle null."""
    ax = d["M"][sel]
    ch = (d["pt_l"][sel], d["pt_s"][sel])
    Ac, Ha, Hb, MI = quantile_axis(ax, ch, ev_per_bin=ev_per_bin)
    if len(Ac) < 10:
        results.append({"tool": "5_symbolic", "metric": "n/a (few axis pts)",
                        "control": "n/a", "hits": False})
        return
    # SINDy-style poly library degree<=2 in (H_a, H_b)
    lib = np.column_stack([np.ones_like(Ha), Ha, Hb, Ha ** 2, Hb ** 2, Ha * Hb])
    beta, r2_poly, _ = ols(lib, MI)
    rng = np.random.default_rng(2)
    nullr2 = [ols(lib, MI[rng.permutation(len(MI))])[1] for _ in range(40)]
    nullp95 = float(np.percentile(nullr2, 95))
    detail = {"sindy_poly2_R2": float(r2_poly), "null_p95": nullp95,
              "coeffs[1,Ha,Hb,Ha2,Hb2,HaHb]": [round(float(x), 4) for x in beta]}
    pysr_eq = None
    if try_pysr:
        try:
            from pysr import PySRRegressor
            X = np.column_stack([Ha, Hb])
            model = PySRRegressor(niterations=20, binary_operators=["+", "*", "-"],
                                  unary_operators=["square", "exp"],
                                  maxsize=12, verbosity=0, progress=False,
                                  random_state=0, deterministic=True,
                                  parallelism="serial")
            model.fit(X, MI)
            pysr_eq = str(model.get_best()["equation"])
            yhat = model.predict(X)
            ss_res = float(np.sum((MI - yhat) ** 2))
            ss_tot = float(np.sum((MI - MI.mean()) ** 2))
            detail["pysr_R2"] = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
            detail["pysr_equation"] = pysr_eq
        except Exception as e:
            detail["pysr"] = f"unavailable: {type(e).__name__}: {str(e)[:80]}"
    hit = (r2_poly > 0.6) and (r2_poly > nullp95 + 0.1)
    results.append({
        "tool": "5_symbolic_MI~f(Ha,Hb) (SINDy poly2 + PySR if avail)",
        "metric": f"SINDy poly2 R2 = {r2_poly:.3f}"
                  + (f"; PySR R2 = {detail.get('pysr_R2'):.3f}"
                     if isinstance(detail.get("pysr_R2"), float) else ""),
        "control": f"shuffle-null p95 = {nullp95:.3f}",
        "hits": bool(hit), "detail": detail})


# ============================================ TOOL 6: direct coupling references
def tool6_direct_coupling(d, sel, results):
    """Per-event cross-correlation of the two channels (Pearson/Spearman) and a
    binned conditional-MI ("transfer-entropy-style") gain, REAL pairs vs event-MIXED
    control.  References for whether the two channels are coupled at all."""
    from scipy.stats import pearsonr, spearmanr
    a, b = d["pt_l"][sel], d["pt_s"][sel]
    pear = float(pearsonr(a, b)[0])
    spear = float(spearmanr(a, b)[0])
    # event-mixed control: scramble channel B across events
    rng = np.random.default_rng(3)
    bmix = b[rng.permutation(len(b))]
    pear_mix = float(pearsonr(a, bmix)[0])
    # MI real vs mixed (histogram)
    mi_real = mi_hist(a, b, bins=24)
    mi_mix = mi_hist(a, bmix, bins=24)
    # TE-style: does knowing M (the axis) plus A reduce uncertainty in B beyond A
    # alone?  proxy = MI(B; [A,M]) - MI(B;A) via binned conditional, real vs mixed.
    M = d["M"][sel]
    def cond_gain(a_, b_, m_):
        # I(B; A, M) approx via 3D hist - I(B;A)
        c3, edges = np.histogramdd(np.column_stack([b_, a_, m_]), bins=6)
        p = c3 / c3.sum()
        pb = p.sum(axis=(1, 2), keepdims=True)
        pam = p.sum(axis=0, keepdims=True)
        nz = p > 0
        i_b_am = float(np.sum(p[nz] * np.log(p[nz] / (pb * pam)[nz])))
        return i_b_am - mi_hist(a_, b_, bins=6)
    cg_real = cond_gain(a, b, M)
    cg_mix = cond_gain(a, bmix, M)
    hit = (abs(pear) > 0.1 and abs(pear) > abs(pear_mix) + 0.05) or \
          (mi_real > mi_mix + 0.02)
    results.append({
        "tool": "6_direct_coupling (Pearson/Spearman/MI, real vs event-mixed)",
        "metric": f"Pearson={pear:.3f} Spearman={spear:.3f} MI_real={mi_real:.3f}",
        "control": f"mixed Pearson={pear_mix:.3f} MI_mixed={mi_mix:.3f}",
        "hits": bool(hit),
        "detail": {"pearson": pear, "spearman": spear,
                   "pearson_mixed": pear_mix, "mi_real": mi_real,
                   "mi_mixed": mi_mix, "cond_gain_real": float(cg_real),
                   "cond_gain_mixed": float(cg_mix)}})


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    if args.canary:
        d = load_pairs(DATA, max_rows=8000)
        ev_per_bin = 250
        out_path = "probe_flow_dipole_strong_battery_canary.json"
        try_pysr = False
    else:
        d = load_pairs(DATA)
        ev_per_bin = 200
        out_path = "probe_flow_dipole_strong_battery_results.json"
        try_pysr = True

    sel = np.ones(len(d["M"]), dtype=bool)  # all pairs (signal); controls inside tools
    results = []

    tool1_baseline_dipole(d, sel, ev_per_bin, results)
    tool2_estimator_swap(d, sel, ev_per_bin, results)
    tool3_raw_dipole(d, sel, ev_per_bin, results)
    tool4_femtoscopy(d, sel, results,
                     n_mix=(6000 if args.canary else None))
    tool5_symbolic(d, sel, ev_per_bin, results, try_pysr=try_pysr)
    tool6_direct_coupling(d, sel, results)

    n_hits = sum(1 for r in results if r["hits"])
    # CRITICAL diagnostic logic
    femto = next((r for r in results if r["tool"].startswith("4_")), None)
    dipoles = [r for r in results if r["tool"][0] in ("1", "2", "3")]
    dipole_hit = any(r["hits"] for r in dipoles)
    blind_flag = bool(femto and femto["hits"] and not dipole_hit)

    out = {
        "probe": "STRONG-force two-particle BATTERY (real CMS dimuon open data)",
        "data_source": {
            "url": "https://opendata.cern.ch/record/303/files/dimuon.csv",
            "record": "CMS open data record/303 (events with two muons, 2010)",
            "n_events_used": int(len(d["M"])),
            "honesty_note": "No per-jet/dijet/hadron-pair CSV exists in CMS "
                            "open-data EDUCATION (API query for dijet/jet/two-jet "
                            "csv returned ZERO csv hits). Closest REAL strong-"
                            "interaction two-particle CSV is this broad dimuon set: "
                            "the mass spectrum is dominated by J/psi (3.097) and "
                            "Upsilon (~9.46), which are QCD (strong-force) BOUND "
                            "STATES (charmonium/bottomonium), plus the EW Z. The "
                            "PAIR is QCD-produced for the quarkonium component; the "
                            "per-muon tracks are leptonic, not hadronic. A true "
                            "pion/kaon femtoscopy sample would be preferable but is "
                            "not available as a CSV."},
        "channels": "leading/subleading muon per-object observable (pt/energy/rapidity)",
        "axes": "pair invariant mass M (J/psi, Upsilon, Z resonances) and femtoscopy q=|p_A-p_B|",
        "tools": results,
        "n_tools": len(results),
        "n_hits": n_hits,
        "dipole_tools_hit": dipole_hit,
        "known_femtoscopy_visible": bool(femto["hits"]) if femto else None,
        "CRITICAL_dipole_blind_to_real_correlation": blind_flag,
        "reading_guard": "Each tool = one data point, no verdict. Weight by "
                         "consistency; a lone hit is not promoted. NO new-law claim.",
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    # console table
    print(f"\n=== STRONG BATTERY ({'CANARY' if args.canary else 'FULL'}) "
          f"{out['runtime_s']}s  n_events={len(d['M'])} ===")
    print(f"{'tool':<52}{'hits':<6}metric / control")
    for r in results:
        print(f"{r['tool'][:50]:<52}{str(r['hits']):<6}{r['metric']}")
        print(f"{'':<58}control: {r['control']}")
    print(f"\nn_hits={n_hits}/{len(results)}  dipole_hit={dipole_hit}  "
          f"femto_visible={out['known_femtoscopy_visible']}")
    if blind_flag:
        print("CRITICAL: dipole tool BLIND to a real correlation present here "
              "(femtoscopy C(q) visible, dipole tools not).")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
