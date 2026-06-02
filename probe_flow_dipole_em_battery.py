"""
EM DIPOLE TOOL BATTERY -- real HBT intensity-interferometry data.

"First try, not only try." The windowed-entropy -> MI -> dMI/d(axis) dipole
tool gave LOW R2 on gravity (LIGO) and weak (CMS). Here we run a BATTERY of
DIFFERENT extraction tools in quick succession on the SAME real EM 2-channel
data to see if ANY tool surfaces genuine 2-channel structure the first tool
missed. For EACH tool we report a clear metric vs an appropriate null/control,
plus a boolean "hits" (structure clearly beats control/null).

Data: Zenodo 5113016 "Enhanced Photonic Maxwell's Demon with Correlated Baths"
g2data raw timetags. tab-separated (tagger-channel, clock-cycle), RES=156.25 ps.
  SIGNAL  : splitThermal       A=ch11, B=ch15  (one thermal field split 50/50 ->
            common field -> HBT photon bunching: g2(0) > 1, correlated streams)
  CONTROL : uncorrelatedThermal A=ch15, B=ch16 (dataset's "uncorrelated baths")

EVERY tool is run on BOTH signal and control. A genuine result shows structure
in signal and NOT (or much weaker) in control.

HONESTY: no novelty/new-law claims. Each tool result is ONE data point. No
verdicts. The g2 (tool #4) is the KNOWN real EM 2-channel correlation -- if it
is visible but the dipole tools are not, that is "the dipole tool is blind to a
correlation that demonstrably exists here," reported explicitly.

Run:  python probe_flow_dipole_em_battery.py [--canary]
"""
import json
import time
import argparse
import numpy as np

RES = 156.25e-12  # tagger clock resolution (s)

SOURCES = {
    "split_correlated": ("data/em/splitThermal_rawtags.txt", 11, 15),
    "uncorrelated_control": ("data/em/uncorrelatedThermal_rawtags.txt", 15, 16),
}


# ----------------------------------------------------------------------------- loaders
def load_tags(path, chA, chB):
    """Load raw timetags once, return (tA, tB, T) in SECONDS, t-zeroed."""
    d = np.loadtxt(path)
    ch = d[:, 0].astype(int)
    t = d[:, 1]
    t0 = t.min()
    T = (t.max() - t0) * RES
    tA = (t[ch == chA] - t0) * RES
    tB = (t[ch == chB] - t0) * RES
    return tA, tB, T


def bin_counts(tA, tB, T, dt):
    """Bin two tag streams into count series at bin width dt (s)."""
    nb = int(T / dt)
    a = np.histogram(tA, bins=nb, range=(0, T))[0].astype(float)
    b = np.histogram(tB, bins=nb, range=(0, T))[0].astype(float)
    return a, b


# ----------------------------------------------------------------------------- estimators
def entropy_1d(x, bins=24):
    h, e = np.histogram(x, bins=bins, density=True)
    w = e[1] - e[0]
    p = h * w
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def mi_hist(x, y, bins=14):
    c, _, _ = np.histogram2d(x, y, bins=bins)
    pxy = c / c.sum()
    px = pxy.sum(1, keepdims=True)
    py = pxy.sum(0, keepdims=True)
    nz = pxy > 0
    return float(np.sum(pxy[nz] * np.log(pxy[nz] / (px * py)[nz])))


def mi_ksg(x, y, k=4):
    """Kraskov-Stoegbauer-Grassberger kNN MI estimator (algorithm 1), 1D x,y.
    Pure-numpy, O(n^2). Caller must cap n (<= ~1500)."""
    from scipy.special import digamma
    n = len(x)
    if n < k + 2:
        return 0.0
    # tiny jitter to break count-data ties
    rng = np.random.default_rng(0)
    x = x + rng.normal(0, 1e-6, n)
    y = y + rng.normal(0, 1e-6, n)
    dx = np.abs(x[:, None] - x[None, :])
    dy = np.abs(y[:, None] - y[None, :])
    dz = np.maximum(dx, dy)
    # kth nearest neighbour distance in joint (chebyshev), excluding self
    dz_sort = np.sort(dz, axis=1)
    eps = dz_sort[:, k]  # index k = kth neighbour (0 is self)
    nx = np.array([np.sum(dx[i] < eps[i]) - 1 for i in range(n)])
    ny = np.array([np.sum(dy[i] < eps[i]) - 1 for i in range(n)])
    mi = (digamma(k) + digamma(n)
          - np.mean(digamma(nx + 1) + digamma(ny + 1)))
    return float(max(mi, 0.0))


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2, ss_res


# ============================================================================= TOOL 4
def tool_g2(tA, tB, T, dt, max_lag_bins):
    """The CLASSIC HBT observable: g2(tau) = <I_A(t) I_B(t+tau)> / (<I_A><I_B>).
    Known real 2-channel correlation: peaks at tau=0 for split thermal, flat
    for uncorrelated. Metric = g2(0) (and peak/baseline ratio). Control/null =
    g2 at large |tau| (the uncorrelated baseline, by construction ~1)."""
    a, b = bin_counts(tA, tB, T, dt)
    ma, mb = a.mean(), b.mean()
    if ma <= 0 or mb <= 0:
        return None
    lags = list(range(-max_lag_bins, max_lag_bins + 1))
    g2 = []
    for lag in lags:
        if lag == 0:
            aa, bb = a, b
        elif lag > 0:
            aa, bb = a[:-lag], b[lag:]
        else:
            aa, bb = a[-lag:], b[:lag]
        g2.append(float(np.mean(aa * bb) / (ma * mb)))
    g2 = np.array(g2)
    lags = np.array(lags)
    g2_0 = float(g2[lags == 0][0])
    # baseline = mean of the outer 30% of |tau| (far wings)
    far = np.abs(lags) > 0.7 * max_lag_bins
    base = float(np.mean(g2[far])) if far.any() else 1.0
    peak_ratio = g2_0 / base if base != 0 else float("nan")
    return {
        "dt_s": dt,
        "g2_zero": g2_0,
        "g2_baseline_farwings": base,
        "peak_over_baseline": peak_ratio,
        "g2_curve_lags_us": [float(l * dt * 1e6) for l in lags],
        "g2_curve": [float(v) for v in g2],
    }


# ============================================================================= TOOL 6
def tool_direct_coupling(a, b):
    """Sanity references: zero-lag Pearson cross-correlation + binned-MI; and a
    crude transfer-entropy proxy both directions (1-step, 3-bin discretized).
    Metric = corr / MI / TE. Control comparison handled by running on both
    sources."""
    cc0 = float(np.corrcoef(a, b)[0, 1])
    # lag-1 cross correlations (directional asymmetry hint)
    cc_ab = float(np.corrcoef(a[:-1], b[1:])[0, 1])  # A leads B
    cc_ba = float(np.corrcoef(b[:-1], a[1:])[0, 1])  # B leads A
    mi = mi_hist(a, b, bins=12)

    def te(src, dst, nbins=3):
        # TE(src->dst) = I(dst_{t+1}; src_t | dst_t), discretized
        qs = np.quantile(src, np.linspace(0, 1, nbins + 1))
        qd = np.quantile(dst, np.linspace(0, 1, nbins + 1))
        s = np.clip(np.digitize(src, qs[1:-1]), 0, nbins - 1)
        d = np.clip(np.digitize(dst, qd[1:-1]), 0, nbins - 1)
        dn, dp, sp = d[1:], d[:-1], s[:-1]
        from collections import Counter
        N = len(dn)
        c_jdp_s = Counter(zip(dn, dp, sp))
        c_dp_s = Counter(zip(dp, sp))
        c_jdp = Counter(zip(dn, dp))
        c_dp = Counter(dp)
        te_v = 0.0
        for (i, j, l), c in c_jdp_s.items():
            p_full = c / N
            p1 = c / c_dp_s[(j, l)]                  # P(dn|dp,sp)
            p2 = c_jdp[(i, j)] / c_dp[j]             # P(dn|dp)
            if p1 > 0 and p2 > 0:
                te_v += p_full * np.log(p1 / p2)
        return float(max(te_v, 0.0))

    te_ab = te(a, b)
    te_ba = te(b, a)
    return {
        "pearson_zero_lag": cc0,
        "pearson_lag1_A_leads_B": cc_ab,
        "pearson_lag1_B_leads_A": cc_ba,
        "mi_hist_bits_nats": mi,
        "transfer_entropy_A_to_B": te_ab,
        "transfer_entropy_B_to_A": te_ba,
    }


# ============================================================================= TOOLS 1-3
def build_axis(a, b, lags, win_bins, n_win, h_bins=12, mi_bins=8,
               mi_mode="hist", ksg_cap=1200):
    """Windowed entropy -> MI along flow axis = inter-detector lag tau.
    mi_mode in {hist, ksg}. Returns arrays (Tau, Ha, Hb, MI)."""
    Tau, Ha, Hb, MI = [], [], [], []
    n = len(a)
    for lag in lags:
        if lag == 0:
            aa, bb = a, b
        else:
            aa, bb = a[:-lag], b[lag:]
        m = min(len(aa), len(bb))
        aa, bb = aa[:m], bb[:m]
        W = win_bins
        k = min(n_win, m // W)
        if k < 8:
            continue
        ha_l, hb_l, mi_l = [], [], []
        for i in range(k):
            sa = aa[i * W:(i + 1) * W]
            sb = bb[i * W:(i + 1) * W]
            ha_l.append(entropy_1d(sa, bins=h_bins))
            hb_l.append(entropy_1d(sb, bins=h_bins))
            if mi_mode == "ksg":
                if len(sa) > ksg_cap:
                    idx = np.linspace(0, len(sa) - 1, ksg_cap).astype(int)
                    mi_l.append(mi_ksg(sa[idx], sb[idx]))
                else:
                    mi_l.append(mi_ksg(sa, sb))
            else:
                mi_l.append(mi_hist(sa, sb, bins=mi_bins))
        Tau.append(float(lag)); Ha.append(float(np.mean(ha_l)))
        Hb.append(float(np.mean(hb_l))); MI.append(float(np.mean(mi_l)))
    return (np.array(Tau), np.array(Ha), np.array(Hb), np.array(MI))


def dipole_fit_metrics(Tau, Ha, Hb, MI, n_null=40, seed=0):
    """Fit paper dipole form dMI/dtau ~ c_self*H^2 + c_cross*Ha*Hb + lin + const
    vs shuffle-null and 1-D MI-relaxation. Returns R2_full, R2_relax,
    shuffle_null_R2_mean(+std), dR2_cross. 'hits' = R2_full clearly beats BOTH
    shuffle-null AND relax."""
    if len(Tau) < 8:
        return None
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    dMI = np.gradient(MIs, Tau)
    sl = slice(1, len(Tau) - 1)
    Ha_, Hb_, MI_, dMI_ = Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha_**2, Hb_**2, Ha_ * Hb_
    one = np.ones_like(Ha_)
    n = len(dMI_)
    X_full = np.column_stack([Ha2, Hb2, HaHb, Ha_, Hb_, one])
    X_self = np.column_stack([Ha2, Hb2, Ha_, Hb_, one])
    X_relax = np.column_stack([MI_, one])
    _, r2_full, _ = ols(X_full, dMI_)
    _, r2_self, _ = ols(X_self, dMI_)
    _, r2_relax, _ = ols(X_relax, dMI_)
    rng = np.random.default_rng(seed)
    nullr2 = []
    for _ in range(n_null):
        p = rng.permutation(n)
        Xn = np.column_stack([Ha2[p], Hb2[p], HaHb[p], Ha_[p], Hb_[p], one])
        _, r2n, _ = ols(Xn, dMI_)
        nullr2.append(r2n)
    null_mean, null_std = float(np.mean(nullr2)), float(np.std(nullr2))
    return {
        "n_axis_points": int(n),
        "R2_full": float(r2_full),
        "R2_relax_1d": float(r2_relax),
        "dR2_cross": float(r2_full - r2_self),
        "shuffle_null_R2_mean": null_mean,
        "shuffle_null_R2_std": null_std,
        "MI_range": [float(MI_.min()), float(MI_.max())],
    }


def tool_raw_signal_dipole(a, b, lags, win_bins, n_win, seed=0):
    """RAW-SIGNAL dipole: same flow form but on RAW count means/products per
    window instead of entropy operators. Target = d(cov_AB)/dtau along lag;
    regressors = self <n_A^2>,<n_B^2>, cross <n_A n_B>, linear means.
    Metric = R2 vs shuffle-null."""
    Tau, Sa, Sb, Cab, Ma, Mb = [], [], [], [], [], []
    for lag in lags:
        if lag == 0:
            aa, bb = a, b
        else:
            aa, bb = a[:-lag], b[lag:]
        m = min(len(aa), len(bb)); aa, bb = aa[:m], bb[:m]
        k = min(n_win, m // win_bins)
        if k < 8:
            continue
        sa_l, sb_l, cab_l, ma_l, mb_l = [], [], [], [], []
        for i in range(k):
            x = aa[i * win_bins:(i + 1) * win_bins]
            y = bb[i * win_bins:(i + 1) * win_bins]
            sa_l.append(np.mean(x**2)); sb_l.append(np.mean(y**2))
            cab_l.append(np.mean(x * y)); ma_l.append(np.mean(x)); mb_l.append(np.mean(y))
        Tau.append(float(lag)); Sa.append(np.mean(sa_l)); Sb.append(np.mean(sb_l))
        Cab.append(np.mean(cab_l)); Ma.append(np.mean(ma_l)); Mb.append(np.mean(mb_l))
    if len(Tau) < 8:
        return None
    Tau = np.array(Tau)
    Cab = np.array(Cab)
    # connected covariance as flow quantity
    cov = Cab - np.array(Ma) * np.array(Mb)
    covs = np.convolve(cov, np.ones(3) / 3, mode="same")
    dcov = np.gradient(covs, Tau)
    sl = slice(1, len(Tau) - 1)
    Sa_, Sb_, Cab_, Ma_, Mb_, dcov_ = (np.array(Sa)[sl], np.array(Sb)[sl],
                                       Cab[sl], np.array(Ma)[sl], np.array(Mb)[sl], dcov[sl])
    one = np.ones_like(Sa_)
    n = len(dcov_)
    X = np.column_stack([Sa_, Sb_, Cab_, Ma_, Mb_, one])
    _, r2, _ = ols(X, dcov_)
    rng = np.random.default_rng(seed)
    nullr2 = []
    for _ in range(40):
        p = rng.permutation(n)
        Xn = np.column_stack([Sa_[p], Sb_[p], Cab_[p], Ma_[p], Mb_[p], one])
        _, r2n, _ = ols(Xn, dcov_)
        nullr2.append(r2n)
    return {
        "n_axis_points": int(n),
        "R2_full": float(r2),
        "shuffle_null_R2_mean": float(np.mean(nullr2)),
        "shuffle_null_R2_std": float(np.std(nullr2)),
        "cov_range": [float(cov.min()), float(cov.max())],
    }


# ============================================================================= TOOL 5a (SINDy)
def poly_library(cols, names, degree=2):
    feats = [np.ones_like(cols[0])]
    fnames = ["1"]
    p = len(cols)
    feats += list(cols); fnames += list(names)
    if degree >= 2:
        for i in range(p):
            for j in range(i, p):
                feats.append(cols[i] * cols[j]); fnames.append(f"{names[i]}*{names[j]}")
    return np.column_stack(feats), fnames


def stlsq(X, y, thresh=0.05, iters=10):
    """Sequentially-thresholded least squares (SINDy core)."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for _ in range(iters):
        small = np.abs(beta) < thresh
        if small.all():
            break
        beta[small] = 0.0
        big = ~small
        if not big.any():
            break
        b2, *_ = np.linalg.lstsq(X[:, big], y, rcond=None)
        beta = np.zeros_like(beta)
        beta[big] = b2
    return beta


def tool_sindy(Tau, Ha, Hb, MI, seed=0):
    """SINDy on the operator series: find sparse governing relation for dMI/dtau
    from a deg-2 library of {Ha, Hb}. Metric = R2 of the recovered sparse model
    vs shuffle-null R2."""
    if len(Tau) < 8:
        return None
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    dMI = np.gradient(MIs, Tau)
    sl = slice(1, len(Tau) - 1)
    Ha_, Hb_, dMI_ = Ha[sl], Hb[sl], dMI[sl]
    # standardize columns for fair thresholding
    cols = [Ha_, Hb_]
    names = ["Ha", "Hb"]
    X, fnames = poly_library(cols, names, degree=2)
    Xs = X.copy()
    sc = X.std(0); sc[sc == 0] = 1.0
    Xs = X / sc
    y = (dMI_ - dMI_.mean())
    ss_tot = np.sum((y - y.mean()) ** 2)
    # threshold scaled to coefficient magnitude so SINDy can actually keep terms
    # if they carry signal (small thresh => closer to plain LS; STLSQ still
    # prunes terms that contribute nothing). Entropy regressors are ~tau-flat by
    # construction, so an all-pruned model is itself an honest 'blind' signal.
    beta = stlsq(Xs.copy(), y, thresh=0.01)
    yhat = Xs @ beta
    ss_res = np.sum((y - yhat) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    active = {fnames[i]: float(beta[i] / sc[i]) for i in range(len(fnames)) if beta[i] != 0}
    # full (un-thresholded) LS R2 as upper bound on what this library can do
    bls, r2_ls, _ = ols(Xs, y)
    rng = np.random.default_rng(seed)
    nullr2 = []
    for _ in range(30):
        p = rng.permutation(len(y))
        bn = stlsq((Xs[p]).copy(), y, thresh=0.01)
        yh = Xs[p] @ bn
        nullr2.append(1 - np.sum((y - yh) ** 2) / ss_tot)
    return {
        "n_axis_points": int(len(y)),
        "R2": r2,
        "R2_full_LS_upper_bound": float(r2_ls),
        "active_terms": active,
        "n_active_terms": int(np.sum(beta != 0)),
        "shuffle_null_R2_mean": float(np.mean(nullr2)),
        "shuffle_null_R2_std": float(np.std(nullr2)),
    }


# ============================================================================= TOOL 5b (PySR)
def tool_pysr(Tau, Ha, Hb, MI, niter=20, seed=0):
    """PySR symbolic regression: find ANY governing relation for dMI/dtau from
    (Tau, Ha, Hb, MI). Metric = best-equation R2 vs shuffle-null best R2."""
    if len(Tau) < 8:
        return None
    try:
        from pysr import PySRRegressor
    except Exception as e:
        return {"error": f"pysr unavailable: {e}"}
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    dMI = np.gradient(MIs, Tau)
    sl = slice(1, len(Tau) - 1)
    Xin = np.column_stack([Tau[sl], Ha[sl], Hb[sl], MI[sl]])
    y = dMI[sl]

    def fit_get_r2(Xf, yf):
        m = PySRRegressor(
            niterations=niter, binary_operators=["+", "-", "*", "/"],
            unary_operators=["square", "exp"], maxsize=14,
            populations=8, population_size=30, verbosity=0,
            progress=False, deterministic=True, parallelism="serial",
            random_state=seed, temp_equation_file=True,
        )
        m.fit(Xf, yf)
        yh = m.predict(Xf)
        ss_res = np.sum((yf - yh) ** 2); ss_tot = np.sum((yf - yf.mean()) ** 2)
        r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
        try:
            eq = str(m.get_best()["equation"])
        except Exception:
            eq = "n/a"
        return r2, eq

    r2, eq = fit_get_r2(Xin, y)
    rng = np.random.default_rng(seed + 1)
    p = rng.permutation(len(y))
    r2_null, _ = fit_get_r2(Xin, y[p])
    return {
        "n_axis_points": int(len(y)),
        "R2": r2,
        "best_equation": eq,
        "shuffle_null_R2": r2_null,
    }


# ============================================================================= driver
def run_source(name, path, chA, chB, cfg, run_pysr):
    tA, tB, T = load_tags(path, chA, chB)
    meta = {"n_tags_A": int(len(tA)), "n_tags_B": int(len(tB)),
            "rate_A_hz": float(len(tA) / T), "rate_B_hz": float(len(tB) / T),
            "T_s": float(T)}

    tools = {}

    # ---- TOOL 4: g2(tau) HBT (the KNOWN correlation). dt swept.
    g2_by_dt = {}
    for dt in cfg["g2_dts"]:
        g = tool_g2(tA, tB, T, dt, cfg["g2_max_lag"])
        if g:
            g2_by_dt[f"{dt:.0e}s"] = {kk: g[kk] for kk in
                                      ("g2_zero", "g2_baseline_farwings", "peak_over_baseline")}
    # finest dt as headline
    fine = cfg["g2_dts"][0]
    gfull = tool_g2(tA, tB, T, fine, cfg["g2_max_lag"])
    g2_zero = gfull["g2_zero"] if gfull else float("nan")
    g2_base = gfull["g2_baseline_farwings"] if gfull else float("nan")
    tools["tool4_g2_HBT"] = {
        "metric_name": "g2(0) and peak/baseline ratio",
        "g2_zero": g2_zero, "g2_baseline": g2_base,
        "peak_over_baseline": gfull["peak_over_baseline"] if gfull else float("nan"),
        "by_dt": g2_by_dt,
        "g2_curve_finest": {"lags_us": gfull["g2_curve_lags_us"],
                            "g2": gfull["g2_curve"]} if gfull else None,
        # hit: g2(0) clearly above baseline (>5% bunching above the wing baseline)
        "control_value": g2_base,
        "hits": bool(gfull and (g2_zero - g2_base) > 0.02 and gfull["peak_over_baseline"] > 1.02),
    }

    # ---- prep count series at the dipole bin for tools 1,2,3,5,6
    dt0 = cfg["dipole_dt"]
    a0, b0 = bin_counts(tA, tB, T, dt0)
    lags = cfg["lags"]
    wb, nw = cfg["win_bins"], cfg["n_win"]

    # ---- TOOL 1: baseline dipole, binning-timescale sweep + lag-axis sweep
    t1 = {}
    for dt in cfg["dipole_dt_sweep"]:
        a, b = bin_counts(tA, tB, T, dt)
        Tau, Ha, Hb, MI = build_axis(a, b, lags, wb, nw, mi_mode="hist")
        fit = dipole_fit_metrics(Tau, Ha, Hb, MI)
        if fit:
            t1[f"dt={dt:.0e}s"] = {"R2_full": fit["R2_full"],
                                   "shuffle_null": fit["shuffle_null_R2_mean"],
                                   "R2_relax": fit["R2_relax_1d"],
                                   "dR2_cross": fit["dR2_cross"]}
    # headline: best dt by margin over null
    if t1:
        best_dt = max(t1, key=lambda kk: t1[kk]["R2_full"] - t1[kk]["shuffle_null"])
        h = t1[best_dt]
        tools["tool1_dipole_binsweep"] = {
            "metric_name": "R2_full(dipole) over shuffle-null & 1D-relax, swept over bin dt",
            "best_dt": best_dt, "R2_full": h["R2_full"],
            "control_value": max(h["shuffle_null"], h["R2_relax"]),
            "shuffle_null": h["shuffle_null"], "R2_relax_1d": h["R2_relax"],
            "dR2_cross": h["dR2_cross"], "by_dt": t1,
            "hits": bool(h["R2_full"] - max(h["shuffle_null"], h["R2_relax"]) > 0.25
                         and h["dR2_cross"] > 0.05),
        }
    else:
        tools["tool1_dipole_binsweep"] = {"metric_name": "R2_full dipole",
                                          "note": "no stable axis", "hits": False,
                                          "control_value": None}

    # ---- TOOL 2: KSG kNN MI estimator + window-size variation
    t2 = {}
    for wbk in cfg["ksg_win_bins"]:
        Tau, Ha, Hb, MI = build_axis(a0, b0, lags, wbk, nw, mi_mode="ksg",
                                     ksg_cap=cfg["ksg_cap"])
        fit = dipole_fit_metrics(Tau, Ha, Hb, MI)
        if fit:
            t2[f"win={wbk}"] = {"R2_full": fit["R2_full"],
                                "shuffle_null": fit["shuffle_null_R2_mean"],
                                "R2_relax": fit["R2_relax_1d"],
                                "MI_range": fit["MI_range"]}
    if t2:
        best_w = max(t2, key=lambda kk: t2[kk]["R2_full"] - t2[kk]["shuffle_null"])
        h = t2[best_w]
        tools["tool2_dipole_ksg"] = {
            "metric_name": "R2_full(dipole) with KSG kNN MI, swept window size",
            "best_window": best_w, "R2_full": h["R2_full"],
            "control_value": max(h["shuffle_null"], h["R2_relax"]),
            "shuffle_null": h["shuffle_null"], "R2_relax_1d": h["R2_relax"],
            "MI_range": h["MI_range"], "by_window": t2,
            "hits": bool(h["R2_full"] - max(h["shuffle_null"], h["R2_relax"]) > 0.25),
        }
    else:
        tools["tool2_dipole_ksg"] = {"metric_name": "R2 dipole KSG",
                                     "note": "no stable axis", "hits": False,
                                     "control_value": None}

    # ---- TOOL 3: RAW-SIGNAL dipole (covariance flow on raw counts)
    raw = tool_raw_signal_dipole(a0, b0, lags, wb, nw)
    if raw:
        tools["tool3_raw_signal_dipole"] = {
            "metric_name": "R2 of d(cov_AB)/dtau ~ self+cross+linear on RAW counts, vs shuffle-null",
            "R2_full": raw["R2_full"], "control_value": raw["shuffle_null_R2_mean"],
            "shuffle_null": raw["shuffle_null_R2_mean"], "cov_range": raw["cov_range"],
            "hits": bool(raw["R2_full"] - raw["shuffle_null_R2_mean"] > 0.25),
        }
    else:
        tools["tool3_raw_signal_dipole"] = {"metric_name": "R2 raw dipole",
                                            "note": "no stable axis", "hits": False,
                                            "control_value": None}

    # ---- TOOL 5a: SINDy
    Tau, Ha, Hb, MI = build_axis(a0, b0, lags, wb, nw, mi_mode="hist")
    sindy = tool_sindy(Tau, Ha, Hb, MI)
    if sindy:
        tools["tool5a_sindy"] = {
            "metric_name": "R2 of sparse SINDy model for dMI/dtau over {Ha,Hb} deg-2, vs shuffle-null",
            "R2": sindy["R2"], "R2_full_LS_upper_bound": sindy["R2_full_LS_upper_bound"],
            "control_value": sindy["shuffle_null_R2_mean"],
            "shuffle_null": sindy["shuffle_null_R2_mean"],
            "active_terms": sindy["active_terms"], "n_active_terms": sindy["n_active_terms"],
            "hits": bool(sindy["R2"] - sindy["shuffle_null_R2_mean"] > 0.25),
        }
    else:
        tools["tool5a_sindy"] = {"metric_name": "SINDy R2", "note": "no stable axis",
                                 "hits": False, "control_value": None}

    # ---- TOOL 5b: PySR (heavy; only on full run)
    if run_pysr:
        pysr_res = tool_pysr(Tau, Ha, Hb, MI, niter=cfg["pysr_niter"])
        if pysr_res and "error" not in pysr_res:
            tools["tool5b_pysr"] = {
                "metric_name": "R2 of best PySR equation for dMI/dtau, vs shuffle-null PySR R2",
                "R2": pysr_res["R2"], "control_value": pysr_res["shuffle_null_R2"],
                "shuffle_null": pysr_res["shuffle_null_R2"],
                "best_equation": pysr_res["best_equation"],
                "hits": bool(pysr_res["R2"] - pysr_res["shuffle_null_R2"] > 0.3),
            }
        else:
            tools["tool5b_pysr"] = {"metric_name": "PySR R2",
                                    "note": (pysr_res or {}).get("error", "skipped"),
                                    "hits": False, "control_value": None}

    # ---- TOOL 6: direct coupling references (corr, MI, TE)
    dc = tool_direct_coupling(a0, b0)
    tools["tool6_direct_coupling"] = {
        "metric_name": "zero-lag Pearson corr / hist-MI / transfer-entropy both dirs (cross-source compare)",
        "pearson_zero_lag": dc["pearson_zero_lag"],
        "mi_hist": dc["mi_hist_bits_nats"],
        "te_A_to_B": dc["transfer_entropy_A_to_B"],
        "te_B_to_A": dc["transfer_entropy_B_to_A"],
        "pearson_lag1_A_leads_B": dc["pearson_lag1_A_leads_B"],
        "pearson_lag1_B_leads_A": dc["pearson_lag1_B_leads_A"],
        "control_value": "see uncorrelated_control source",
        # within-source hit flag: nonzero corr / MI. Cross-source comparison done in summary.
        "hits": bool(abs(dc["pearson_zero_lag"]) > 0.02 or dc["mi_hist_bits_nats"] > 0.01),
    }

    return {"meta": meta, "tools": tools}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    if args.canary:
        cfg = {
            "g2_dts": [5e-7, 2e-6],
            "g2_max_lag": 40,
            "dipole_dt": 2e-6,
            "dipole_dt_sweep": [1e-6, 2e-6],
            "lags": [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48],
            "win_bins": 4000, "n_win": 60,
            "ksg_win_bins": [2000], "ksg_cap": 600,
            "pysr_niter": 8,
        }
        out_path = "probe_flow_dipole_em_battery_canary.json"
        run_pysr = False  # too heavy for canary budget
        sources = {"split_correlated": SOURCES["split_correlated"]}
    else:
        cfg = {
            "g2_dts": [2e-7, 5e-7, 1e-6, 2e-6],
            "g2_max_lag": 60,
            "dipole_dt": 2e-6,
            "dipole_dt_sweep": [5e-7, 1e-6, 2e-6, 4e-6],
            "lags": [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 96],
            "win_bins": 4000, "n_win": 200,
            "ksg_win_bins": [1500, 3000], "ksg_cap": 1000,
            "pysr_niter": 20,
        }
        out_path = "probe_flow_dipole_em_battery_results.json"
        run_pysr = True
        sources = SOURCES

    by_source = {}
    for name, (path, chA, chB) in sources.items():
        print(f"[{name}] loading + running battery...", flush=True)
        by_source[name] = run_source(name, path, chA, chB, cfg, run_pysr)

    # cross-source summary table per tool
    summary = {}
    sig = by_source.get("split_correlated", {}).get("tools", {})
    ctl = by_source.get("uncorrelated_control", {}).get("tools", {})
    for tname in sig:
        s = sig[tname]
        c = ctl.get(tname, {})
        summary[tname] = {
            "metric_name": s.get("metric_name"),
            "signal_metric": {k: s[k] for k in s if k not in
                              ("metric_name", "by_dt", "by_window", "g2_curve_finest")},
            "control_metric": {k: c[k] for k in c if k not in
                               ("metric_name", "by_dt", "by_window", "g2_curve_finest")}
            if c else None,
            "signal_hits": s.get("hits"),
            "control_hits": c.get("hits") if c else None,
        }

    result = {
        "probe": "EM dipole TOOL BATTERY -- real HBT intensity interferometry "
                 "(Zenodo 5113016 g2data raw timetags)",
        "purpose": "first-try-not-only-try: run many DIFFERENT extraction tools on "
                   "the SAME real EM 2-channel data; per tool report metric vs "
                   "null/control + hits bool. g2 (tool4) is the KNOWN correlation: "
                   "if it is visible but dipole tools are not, the dipole tool is "
                   "blind to a correlation that demonstrably exists here.",
        "channels": "split: A=ch11,B=ch15 common thermal field (HBT bunching, SIGNAL); "
                    "uncorrelated: A=ch15,B=ch16 dataset's 'uncorrelated baths' (CONTROL)",
        "honesty": "no novelty/new-law claims. each tool = one data point. no verdicts. "
                   "weight by consistency; do not romanticize a lone hit.",
        "config": cfg,
        "by_source": by_source,
        "cross_source_summary": summary,
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    # console table
    print(f"\n=== EM TOOL BATTERY ({'CANARY' if args.canary else 'FULL'}) "
          f"{result['runtime_s']}s ===")
    print(f"{'tool':28s} {'signal_metric':>16s} {'control/null':>16s}  sig_hit ctl_hit")
    for tname in summary:
        s = summary[tname]["signal_metric"]
        # pick a representative scalar metric
        for key in ("R2_full", "R2", "g2_zero", "peak_over_baseline", "pearson_zero_lag"):
            if key in s and isinstance(s[key], (int, float)):
                sm = s[key]; break
        else:
            sm = float("nan")
        cm = summary[tname]["control_metric"]
        cv = s.get("control_value")
        cv = cv if isinstance(cv, (int, float)) else float("nan")
        print(f"{tname:28s} {sm:16.3f} {cv:16.3f}  "
              f"{str(summary[tname]['signal_hits']):>5s}   {str(summary[tname]['control_hits']):>5s}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
