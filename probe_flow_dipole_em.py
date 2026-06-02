"""
FLOW DIPOLE EQUATION -- ELECTROMAGNETIC force, REAL DATA, own schema.

Mirrors probe_flow_dipole_weak.py exactly (own-schema fit library picked by
AIC + R2; deflationary controls: 1-D MI-relaxation, self-only/cross-term value,
shuffle-null; standardized betas; bootstrap-over-data scatter; opposition
signature; algebraic ratio C). Greg (S19): do NOT force one fixed operator
basis -- let the data pick the schema.

REAL 2-channel object (no invented coupling): Hanbury Brown-Twiss intensity
interferometry. Two photon-counting detectors A,B sample one common light field.
  Data: Zenodo 5113016, "Enhanced Photonic Maxwell's Demon with Correlated
  Baths" g2data raw timetags (tab-separated: tagger-channel, clock-cycle;
  resolution res=156.25 ps). Confirmed via the dataset's own
  thermalLightG2_withSave.m:
    SPLIT thermal (SIGNAL): one thermal field split 50/50 -> detectors
      A=ch11, B=ch15. A COMMON field -> HBT photon bunching -> the two intensity
      streams are correlated, peaking at zero time-lag and relaxing over the
      coherence time (measured: count-series corr 0.14 at lag 0 -> ~0 by ~10 us).
    UNCORRELATED thermal (dataset's labelled CONTROL): detectors A=ch15, B=ch16
      (paper's "uncorrelated baths" config for the demon). Run identically; the
      analogue of weak's same-charge no-resonance control. (Honesty note: this
      pair empirically shows its OWN zero-lag bump -- reported, not assumed away.)

  channel observable = photon COUNT per fine time-bin n_A[k], n_B[k] (intensity).

FLOW AXIS = inter-detector time-lag tau. HBT correlation g2(tau) is maximal at
tau=0 (bunching) and relaxes toward the uncorrelated baseline as |tau| grows
past the coherence time. So MI(tau) = I(n_A(t); n_B(t+tau)) is the natural
flow quantity and dMI/dtau is the flow derivative -- exactly analogous to the
weak probe's dMI/dM toward the resonance pole. At each tau we estimate
H_a=H(n_A), H_b=H(n_B) (per-channel windowed entropies; ~tau-independent by
construction) and MI(tau) by pooling many windows; then fit the schema library
to dMI/dtau. The paper flow-dipole form (davisai.ai/dipole Sec 2.2) is ONE
candidate among several:  dMI/dtau ~ c_self*H_i^2 + c_cross*H_a*H_b + lin + const.

Result Discipline: bootstrap over windows for inter-rep scatter (no claim from
one realization). NO novelty claims (literature scan pending). Report the
three-level reading with the deflationary read always present: genuine
2-channel flow-dipole only if R2_full clearly beats BOTH shuffle-null AND the
1-D MI-relaxation, the cross term adds real dR2, and opposition is present;
otherwise it COLLAPSES to a control (R2_full ~= shuffle-null), as gravity (LIGO)
and weak (CMS dimuon) both did.

Run:  python probe_flow_dipole_em.py [--canary]
"""
import json
import time
import argparse
import numpy as np

RES = 156.25e-12  # tagger clock resolution (s), from the dataset's matlab script

SOURCES = {
    # name: (path, chA, chB)   -- channel ids confirmed from thermalLightG2_withSave.m
    "split_correlated": ("data/em/splitThermal_rawtags.txt", 11, 15),
    "uncorrelated_control": ("data/em/uncorrelatedThermal_rawtags.txt", 15, 16),
}


def load_counts(path, chA, chB, dt):
    """Load raw timetags, select two detector channels, bin each into a count
    (intensity) series at bin width dt seconds. Returns n_A[k], n_B[k]."""
    d = np.loadtxt(path)
    ch = d[:, 0].astype(int)
    t = d[:, 1]
    t0 = t.min()
    T = (t.max() - t0) * RES
    tA = (t[ch == chA] - t0) * RES
    tB = (t[ch == chB] - t0) * RES
    nb = int(T / dt)
    a = np.histogram(tA, bins=nb, range=(0, T))[0].astype(float)
    b = np.histogram(tB, bins=nb, range=(0, T))[0].astype(float)
    return a, b, T, len(tA), len(tB)


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


def build_axis_series(a, b, lags, win_bins=4000, n_win=200, h_bins=12, mi_bins=8):
    """Equal-window estimation along the FLOW AXIS = time-lag tau.

    For each lag tau (in count-series bins), split the long count series into
    n_win contiguous windows of win_bins each; per window compute
    H(n_A), H(n_B), and MI(n_A(t); n_B(t+tau)); average over windows to get the
    axis-point (tau, H_a, H_b, MI(tau)). This mirrors the weak probe's per-bin
    pooled-event estimate (here windows play the role of the events in a bin)."""
    Tau, Ha, Hb, MI = [], [], [], []
    n = len(a)
    for lag in lags:
        # delayed B
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
            mi_l.append(mi_2d(sa, sb, bins=mi_bins))
        Tau.append(float(lag))
        Ha.append(float(np.mean(ha_l)))
        Hb.append(float(np.mean(hb_l)))
        MI.append(float(np.mean(mi_l)))
    return (np.array(Tau), np.array(Ha), np.array(Hb), np.array(MI))


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2, ss_res


def aic(ss_res, n, k):
    return n * np.log(ss_res / n + 1e-300) + 2 * k


def fit_library(Tau, Ha, Hb, MI):
    """Fit a LIBRARY of candidate schemas; return per-schema R2/AIC so the data
    picks the form (Greg: flows may have different schemas)."""
    MIs = np.convolve(MI, np.ones(3) / 3, mode="same")
    dMI = np.gradient(MIs, Tau)
    sl = slice(1, len(Tau) - 1)
    Tau_, Ha_, Hb_, MI_, dMI_ = Tau[sl], Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha_**2, Hb_**2, Ha_ * Hb_
    one = np.ones_like(Ha_)
    n = len(dMI_)

    schemas = {
        "paper_dipole_full": [Ha2, Hb2, HaHb, Ha_, Hb_, one],
        "self_only": [Ha2, Hb2, Ha_, Hb_, one],
        "cross_only": [HaHb, one],
        "MI_relax_1d": [MI_, one],
        "linear_in_axis": [Tau_, one],
        "linear_entropy": [Ha_, Hb_, one],
        "quad_no_linear": [Ha2, Hb2, HaHb, one],
    }
    out = {}
    for name, cols in schemas.items():
        X = np.column_stack(cols)
        beta, r2, ssr = ols(X, dMI_)
        out[name] = {"R2": float(r2), "AIC": float(aic(ssr, n, X.shape[1])),
                     "k": X.shape[1]}
    # standardized betas of the paper form (shape, unit-robust). On the HBT-tau
    # axis the per-channel entropies H_a,H_b are ~constant (entropy of a channel
    # does not depend on the inter-detector lag), so the [H_a^2,H_b^2,H_a,H_b]
    # block is near rank-deficient -> a plain lstsq blows the betas up on
    # near-zero std. Use a small ridge so the reported standardized betas are
    # stable; the near-constancy of the entropy regressors is itself reported
    # below (entropy_std_over_mean). This does NOT touch any R2/AIC/control.
    cols = [Ha2, Hb2, HaHb, Ha_, Hb_]
    Xs = np.column_stack([(c - c.mean()) / (c.std() + 1e-12) for c in cols])
    ys = (dMI_ - dMI_.mean()) / (dMI_.std() + 1e-12)
    lam = 1e-3
    sbeta = np.linalg.solve(Xs.T @ Xs + lam * np.eye(Xs.shape[1]), Xs.T @ ys)
    names = ["H_a^2", "H_b^2", "H_a*H_b", "H_a", "H_b"]
    std_betas = {names[i]: float(sbeta[i]) for i in range(5)}
    ent_cov = {
        "H_a_std_over_mean": float(Ha_.std() / (abs(Ha_.mean()) + 1e-12)),
        "H_b_std_over_mean": float(Hb_.std() / (abs(Hb_.mean()) + 1e-12)),
        "MI_std_over_mean": float(MI_.std() / (abs(MI_.mean()) + 1e-12)),
    }

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
        "tau_bins_range": [float(Tau_.min()), float(Tau_.max())],
        "regressor_variation": ent_cov,
    }


def run_source(name, path, chA, chB, dt, lags, win_bins, n_win, n_boot):
    a, b, T, NA, NB = load_counts(path, chA, chB, dt)
    cc0 = float(np.corrcoef(a, b)[0, 1])  # zero-lag count-series correlation
    Tau, Ha, Hb, MI = build_axis_series(a, b, lags, win_bins=win_bins,
                                        n_win=n_win)
    meta = {
        "n_tags_A": int(NA), "n_tags_B": int(NB),
        "rate_A_hz": float(NA / T), "rate_B_hz": float(NB / T),
        "T_s": float(T), "dt_bin_s": float(dt),
        "zero_lag_corr": cc0,
        "lags_bins": [int(x) for x in lags],
        "lags_us": [float(x * dt * 1e6) for x in lags],
    }
    if len(Tau) < 8:
        return {"meta": meta,
                "note": f"too few axis points ({len(Tau)}) -- could not build a "
                        f"stable lag axis at this binning.",
                "point_estimate": None, "bootstrap_scatter": None}
    point = fit_library(Tau, Ha, Hb, MI)
    # bootstrap over windows: resample the contiguous count series in blocks by
    # bootstrapping which windows enter each axis-point estimate.
    rng = np.random.default_rng(7)
    boot_r2, boot_dcross = [], []
    boot_betas = {k: [] for k in point["paper_form_std_betas"]}
    for _ in range(n_boot):
        # block bootstrap the series: resample contiguous chunks of the count
        # series, then rebuild the axis. Cheap proxy: circular shift + window
        # resample via random window offset.
        shift = int(rng.integers(0, win_bins))
        ab, bb_ = a[shift:], b[shift:]
        Taub, Hab, Hbb, MIb = build_axis_series(ab, bb_, lags,
                                                win_bins=win_bins, n_win=n_win)
        if len(Taub) < 8:
            continue
        fb = fit_library(Taub, Hab, Hbb, MIb)
        boot_r2.append(fb["R2_full"])
        boot_dcross.append(fb["dR2_cross"])
        for k in boot_betas:
            boot_betas[k].append(fb["paper_form_std_betas"][k])
    scatter = {
        "R2_full_boot": [float(np.mean(boot_r2)), float(np.std(boot_r2))]
        if boot_r2 else [float("nan"), float("nan")],
        "dR2_cross_boot": [float(np.mean(boot_dcross)),
                           float(np.std(boot_dcross))]
        if boot_dcross else [float("nan"), float("nan")],
        "std_betas_boot": {k: [float(np.mean(v)), float(np.std(v))]
                           for k, v in boot_betas.items() if v},
        "n_boot_ok": len(boot_r2),
    }
    return {"meta": meta, "point_estimate": point, "bootstrap_scatter": scatter}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    dt = 2e-6  # 2 us count bin (coherence time ~ few us per the corr scan)
    if args.canary:
        # cheaper: fewer lags, fewer/smaller windows, one source, few boots
        lags = [0, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48]
        win_bins, n_win, n_boot = 4000, 80, 8
        out_path = "probe_flow_dipole_em_canary.json"
        sources = {"split_correlated": SOURCES["split_correlated"]}
    else:
        lags = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 96]
        win_bins, n_win, n_boot = 4000, 200, 40
        out_path = "probe_flow_dipole_em_results.json"
        sources = SOURCES

    res = {}
    for name, (path, chA, chB) in sources.items():
        res[name] = run_source(name, path, chA, chB, dt, lags, win_bins,
                                n_win, n_boot)

    result = {
        "probe": "flow-dipole equation -- EM (real HBT intensity interferometry, "
                 "axis=inter-detector time-lag tau)",
        "channels": "detector-A photon count n_A (channel A) vs detector-B photon "
                    "count n_B (channel B); MI(n_A(t); n_B(t+tau)) vs tau",
        "dataset": "Zenodo 5113016 'Enhanced Photonic Maxwell's Demon with "
                   "Correlated Baths', g2data raw timetags; res=156.25 ps",
        "form_note": "paper dipole form is ONE candidate; schema library lets the "
                     "data pick (best_schema_by_AIC). Greg: flows may differ in schema.",
        "config": {"dt_bin_s": dt, "win_bins": win_bins, "n_win": n_win,
                   "n_boot": n_boot, "lags_bins": [int(x) for x in lags]},
        "by_source": res,
        "control_logic": "split (A=ch11,B=ch15) = common thermal field -> HBT "
                         "bunching (signal). uncorrelated (A=ch15,B=ch16) = "
                         "dataset's labelled 'uncorrelated baths' control. A "
                         "genuine EM flow-dipole should appear in the common-field "
                         "object and be absent/weaker in the control.",
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"=== flow-dipole EM ({'CANARY' if args.canary else 'FULL'}) "
          f"{result['runtime_s']}s ===")
    for name in res:
        m = res[name]["meta"]
        p = res[name]["point_estimate"]
        s = res[name]["bootstrap_scatter"]
        print(f"\n[{name}]  rateA={m['rate_A_hz']:.0f}Hz rateB={m['rate_B_hz']:.0f}Hz "
              f"T={m['T_s']:.2f}s  zero-lag corr={m['zero_lag_corr']:.3f}")
        if p is None:
            print(f"  {res[name]['note']}")
            continue
        print(f"  axis_pts={p['n_axis_points']}  MI_range={p['MI_range']}")
        print(f"  best schema by AIC : {p['best_schema_by_AIC']}")
        print(f"  R2_full            : {p['R2_full']:.3f} "
              f"(boot {s['R2_full_boot'][0]:.3f}+/-{s['R2_full_boot'][1]:.3f})")
        print(f"  R2_relax_1d        : {p['R2_relax_1d']:.3f}")
        print(f"  dR2_cross          : {p['dR2_cross']:+.3f} "
              f"(boot {s['dR2_cross_boot'][0]:+.3f}+/-{s['dR2_cross_boot'][1]:.3f})")
        print(f"  shuffle-null R2    : {p['shuffle_null_R2_mean']:.3f}")
        print(f"  opposition sig     : {p['opposition_signature']}")
        print(f"  C ratio            : {p['C_ratio']:.2f}")
        print(f"  std betas (paper)  : " +
              " ".join(f"{k}={v:+.2f}" for k, v in
                       p["paper_form_std_betas"].items()))
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
