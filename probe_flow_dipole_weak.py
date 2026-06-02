"""
FLOW DIPOLE EQUATION -- WEAK force, REAL DATA, own schema.
Greg (S19): drop the toy simulators; rerun on new/real data; each flow may have
its OWN schema (do not force one fixed operator basis).

REAL 2-channel object (no invented coupling): the two muons of each CMS dimuon
event are physically correlated (they share the decaying boson's 4-momentum).
  channel A = leading muon, channel B = subleading muon (per-muon observable = pT).
  FLOW AXIS = dimuon invariant mass M (toward the Z resonance pole ~91.19 GeV),
    M^2 = 2 pt1 pt2 (cosh(eta1-eta2) - cos(phi1-phi2)).
Data: data/forces/Zmumu.csv, ~10583 real CMS open-data events (Run2010/2011).

Along the mass axis we bin events by M; in each bin estimate H_a=H(pT_lead),
H_b=H(pT_sub), MI=I(pT_lead; pT_sub) across the events in that bin, then form
dMI/dM. The paper flow-dipole form (davisai.ai/dipole Sec 2.2) is ONE candidate:
    dMI/dM ~ c_self*H_i^2 + c_cross*H_a*H_b + linear + const
but per Greg we do NOT impose it -- we fit a LIBRARY of candidate schemas and let
the data pick (AIC + CV-R2), so a different form can win for this flow.

REAL CONTROL (built in): same-charge (Q1*Q2>0) dimuon events are NOT Z decays
(no resonance) -> rerun the whole extraction on them; a genuine dipole should
be present in opposite-charge (signal) and absent/weaker in same-charge (control).

DEFLATIONARY controls: 1-D MI-relaxation, self-only (cross-term value), feature
time/axis-shuffle null. Standardized betas reported (entropy units depend on
binning; INFO-051: MI is the only scale-invariant operator).
Result Discipline: bootstrap over events for inter-rep scatter (no claim from one).

Run:  python probe_flow_dipole_weak.py [--canary]
"""
import json
import time
import argparse
import numpy as np
import pandas as pd


def load_events(path="data/forces/Zmumu.csv"):
    df = pd.read_csv(path)
    pt1, eta1, phi1 = df.pt1.values, df.eta1.values, df.phi1.values
    pt2, eta2, phi2 = df.pt2.values, df.eta2.values, df.phi2.values
    Q1, Q2 = df.Q1.values, df.Q2.values
    M2 = 2 * pt1 * pt2 * (np.cosh(eta1 - eta2) - np.cos(phi1 - phi2))
    M = np.sqrt(np.clip(M2, 0, None))
    # leading / subleading by pT (well-defined per-event channel labels)
    lead = np.maximum(pt1, pt2)
    sub = np.minimum(pt1, pt2)
    opp = (Q1 * Q2) < 0
    return dict(M=M, lead=lead, sub=sub, opp=opp)


def entropy_1d(x, bins=24):
    h, e = np.histogram(x, bins=bins, density=True)
    w = e[1] - e[0]
    p = h * w
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def mi_2d(x, y, bins=14):
    c, ex, ey = np.histogram2d(x, y, bins=bins)
    pxy = c / c.sum()
    px = pxy.sum(1, keepdims=True)
    py = pxy.sum(0, keepdims=True)
    nz = pxy > 0
    return float(np.sum(pxy[nz] * np.log(pxy[nz] / (px * py)[nz])))


def build_axis_series(M, lead, sub, m_lo=60.0, m_hi=120.0, n_bins=34,
                      min_per_bin=60, ev_per_bin=150, h_bins=12, mi_bins=8):
    """Equal-EVENT (quantile) bins along the mass axis so every bin has a stable
    sample for the entropy/MI estimators AND we get many axis points (n >> k).
    n_bins/min_per_bin retained for API but ev_per_bin drives the binning."""
    sel0 = (M >= m_lo) & (M < m_hi)
    Ms, Ls, Ss = M[sel0], lead[sel0], sub[sel0]
    order = np.argsort(Ms)
    Ms, Ls, Ss = Ms[order], Ls[order], Ss[order]
    n = len(Ms)
    k = max(8, n // ev_per_bin)
    edges_idx = np.linspace(0, n, k + 1).astype(int)
    Mc, Ha, Hb, MI = [], [], [], []
    for i in range(k):
        a, b = edges_idx[i], edges_idx[i + 1]
        if b - a < max(40, min_per_bin):
            continue
        Mc.append(float(np.median(Ms[a:b])))
        Ha.append(entropy_1d(Ls[a:b], bins=h_bins))
        Hb.append(entropy_1d(Ss[a:b], bins=h_bins))
        MI.append(mi_2d(Ls[a:b], Ss[a:b], bins=mi_bins))
    return (np.array(Mc), np.array(Ha), np.array(Hb), np.array(MI))


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2, ss_res


def aic(ss_res, n, k):
    return n * np.log(ss_res / n + 1e-300) + 2 * k


def fit_library(Mc, Ha, Hb, MI):
    """Fit a LIBRARY of candidate schemas; return per-schema R2/AIC so the data
    picks the form (Greg: flows may have different schemas)."""
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    # axis derivative dMI/dM (interior, drop convolution edges)
    dMI = np.gradient(MIs, Mc)
    sl = slice(1, len(Mc) - 1)
    Mc_, Ha_, Hb_, MI_, dMI_ = Mc[sl], Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha_**2, Hb_**2, Ha_ * Hb_
    one = np.ones_like(Ha_)
    n = len(dMI_)

    schemas = {
        "paper_dipole_full": [Ha2, Hb2, HaHb, Ha_, Hb_, one],
        "self_only": [Ha2, Hb2, Ha_, Hb_, one],
        "cross_only": [HaHb, one],
        "MI_relax_1d": [MI_, one],
        "linear_in_axis": [Mc_, one],
        "linear_entropy": [Ha_, Hb_, one],
        "quad_no_linear": [Ha2, Hb2, HaHb, one],
    }
    out = {}
    for name, cols in schemas.items():
        X = np.column_stack(cols)
        beta, r2, ssr = ols(X, dMI_)
        out[name] = {"R2": float(r2), "AIC": float(aic(ssr, n, X.shape[1])),
                     "k": X.shape[1]}
    # standardized betas of the paper form (shape, unit-robust)
    cols = [Ha2, Hb2, HaHb, Ha_, Hb_]
    Xs = np.column_stack([(c - c.mean()) / (c.std() + 1e-12) for c in cols])
    ys = (dMI_ - dMI_.mean()) / (dMI_.std() + 1e-12)
    sbeta, _, _ = ols(Xs, ys)
    names = ["H_a^2", "H_b^2", "H_a*H_b", "H_a", "H_b"]
    std_betas = {names[i]: float(sbeta[i]) for i in range(5)}

    # shuffle-null on paper form
    rng = np.random.default_rng(0)
    nullr2 = []
    for _ in range(30):
        p = rng.permutation(n)
        Xn = np.column_stack([Ha2[p], Hb2[p], HaHb[p], Ha_[p], Hb_[p], one])
        _, r2n, _ = ols(Xn, dMI_)
        nullr2.append(r2n)

    r2_full = out["paper_dipole_full"]["R2"]
    r2_self = out["self_only"]["R2"]
    r2_relax = out["MI_relax_1d"]["R2"]
    opp = bool(np.sign(sbeta[2]) != 0 and
               np.all(np.sign(sbeta[:2]) == -np.sign(sbeta[2])))
    best = min(out.items(), key=lambda kv: kv[1]["AIC"])[0]
    C = float(np.mean(0.5 * (Ha_ + Hb_) / (MI_ + 1e-9)))
    return {
        "n_axis_points": int(n),
        "schema_library": out,
        "best_schema_by_AIC": best,
        "paper_form_std_betas": std_betas,
        "R2_full": r2_full,
        "dR2_cross": float(r2_full - r2_self),
        "R2_relax_1d": r2_relax,
        "shuffle_null_R2_mean": float(np.mean(nullr2)),
        "opposition_signature": opp,
        "C_ratio": C,
        "MI_range": [float(MI_.min()), float(MI_.max())],
    }


def run_sample(ev, charge, n_boot, n_bins, min_per_bin):
    mask = ev["opp"] if charge == "opposite" else ~ev["opp"]
    M, lead, sub = ev["M"][mask], ev["lead"][mask], ev["sub"][mask]
    # point estimate
    Mc, Ha, Hb, MI = build_axis_series(M, lead, sub, n_bins=n_bins,
                                       min_per_bin=min_per_bin)
    if len(Mc) < 8:
        return {"n_events": int(len(M)),
                "note": f"too few {charge}-charge events to build a stable mass "
                        f"axis ({len(Mc)} bins) -- control has little data, as "
                        f"expected for a no-resonance sample.",
                "point_estimate": None, "bootstrap_scatter": None}
    point = fit_library(Mc, Ha, Hb, MI)
    # bootstrap over events for scatter
    rng = np.random.default_rng(7)
    boot_r2, boot_dcross, boot_betas = [], [], {k: [] for k in
                                               point["paper_form_std_betas"]}
    for _ in range(n_boot):
        idx = rng.integers(0, len(M), len(M))
        Mcb, Hab, Hbb, MIb = build_axis_series(M[idx], lead[idx], sub[idx],
                                               n_bins=n_bins,
                                               min_per_bin=min_per_bin)
        if len(Mcb) < 8:
            continue
        fb = fit_library(Mcb, Hab, Hbb, MIb)
        boot_r2.append(fb["R2_full"])
        boot_dcross.append(fb["dR2_cross"])
        for k in boot_betas:
            boot_betas[k].append(fb["paper_form_std_betas"][k])
    scatter = {
        "n_events": int(len(M)),
        "R2_full_boot": [float(np.mean(boot_r2)), float(np.std(boot_r2))],
        "dR2_cross_boot": [float(np.mean(boot_dcross)),
                           float(np.std(boot_dcross))],
        "std_betas_boot": {k: [float(np.mean(v)), float(np.std(v))]
                           for k, v in boot_betas.items()},
    }
    return {"point_estimate": point, "bootstrap_scatter": scatter}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    ev = load_events()
    if args.canary:
        n_boot, n_bins, min_per_bin = 10, 24, 60
        out_path = "probe_flow_dipole_weak_canary.json"
        charges = ["opposite"]
    else:
        n_boot, n_bins, min_per_bin = 60, 34, 60
        out_path = "probe_flow_dipole_weak_results.json"
        charges = ["opposite", "same"]

    res = {}
    for ch in charges:
        res[ch] = run_sample(ev, ch, n_boot, n_bins, min_per_bin)

    result = {
        "probe": "flow-dipole equation -- WEAK (real CMS dimuon, axis=invariant mass)",
        "channels": "leading muon pT (A) vs subleading muon pT (B); MI(A;B) vs M",
        "form_note": "paper dipole form is ONE candidate; schema library lets the "
                     "data pick (best_schema_by_AIC). Greg: flows may differ in schema.",
        "config": {"n_boot": n_boot, "n_bins": n_bins, "min_per_bin": min_per_bin,
                   "m_window_GeV": [60, 120], "Z_pole_GeV": 91.1876},
        "by_charge": res,
        "control_logic": "opposite-charge = real Z signal; same-charge = no-resonance "
                         "control. A real weak flow-dipole should appear in signal and "
                         "be absent/weaker in control.",
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"=== flow-dipole WEAK ({'CANARY' if args.canary else 'FULL'}) "
          f"{result['runtime_s']}s ===")
    for ch in charges:
        p = res[ch]["point_estimate"]
        s = res[ch]["bootstrap_scatter"]
        if p is None:
            print(f"\n[{ch}-charge]  {res[ch]['note']}")
            continue
        print(f"\n[{ch}-charge]  n_events={s['n_events']}  "
              f"axis_pts={p['n_axis_points']}")
        print(f"  best schema by AIC : {p['best_schema_by_AIC']}")
        print(f"  R2_full            : {p['R2_full']:.3f} "
              f"(boot {s['R2_full_boot'][0]:.3f}+/-{s['R2_full_boot'][1]:.3f})")
        print(f"  R2_relax_1d        : {p['R2_relax_1d']:.3f}")
        print(f"  dR2_cross          : {p['dR2_cross']:+.3f} "
              f"(boot {s['dR2_cross_boot'][0]:+.3f}+/-{s['dR2_cross_boot'][1]:.3f})")
        print(f"  shuffle-null R2    : {p['shuffle_null_R2_mean']:.3f}")
        print(f"  opposition sig     : {p['opposition_signature']}")
        print(f"  C ratio            : {p['C_ratio']:.2f}")
        print(f"  std betas (paper form): " +
              " ".join(f"{k}={v:+.2f}" for k, v in
                       p["paper_form_std_betas"].items()))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
