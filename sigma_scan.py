"""
Sigma scan on OU baseline (Session 6, item B).

Question: does the smallest eigenvalue of the operator covariance
scale with sigma in the OU two-stream baseline? If yes, going to
smaller sigma should recover v5's machine-epsilon spectrum
[4e-13, 8e-13, 3e-11, 3e-4, 4e-4, 2e-3].

Hypothesis from Taylor expansion: the algebraic relations hold to
order (dH)^2 where dH = H - <H>. dH scales with sigma (windowed
entropy of OU has variance proportional to sigma at fixed gamma,
window). The variance along the constrained direction scales as
(dH)^4 ~ sigma^4 * log^4(sigma). At sigma=1 my floor was ~5e-5;
predicted at sigma=0.01: 5e-5 * (0.01)^4 = 5e-13. Matches v5
order.

Method: same kbk_pipeline.py OU baseline, vary sigma from 1.0 down
to 0.01 in 5 steps. Same gamma=0.5. Track bottom-3 eigenvalues and
log-log scaling.

Output: sigma_scan_results.json
"""
import numpy as np
import json
import time
import argparse
from kbk_pipeline import (
    simulate_ou_pair, vasicek_entropy, mi_hist_2d,
    extract_v1, kbk_rank_gap, project_234,
)


def build_M(seed, T, dt, window_s, stride_s, n_real, gamma, sigma,
            mi_bins):
    rng = np.random.default_rng(seed)
    window_samples = int(window_s / dt)
    stride_samples = max(1, int(stride_s / dt))
    M_pool = []
    for r in range(n_real):
        a, b = simulate_ou_pair(T, dt, gamma, sigma, rng)
        n = len(a)
        n_w = (n - window_samples) // stride_samples + 1
        if n_w <= 0:
            continue
        M_r = np.zeros((n_w, 6))
        for w in range(n_w):
            s = w * stride_samples
            e = s + window_samples
            aw = a[s:e]
            bw = b[s:e]
            Ha = vasicek_entropy(aw)
            Hb = vasicek_entropy(bw)
            MI = mi_hist_2d(aw, bw, bins=mi_bins)
            M_r[w, 0] = Ha
            M_r[w, 1] = Hb
            M_r[w, 2] = Ha * Ha
            M_r[w, 3] = Hb * Hb
            M_r[w, 4] = Ha * Hb
            M_r[w, 5] = MI
        M_pool.append(M_r)
    return np.concatenate(M_pool, axis=0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=str, default="sigma_scan_results.json")
    args = p.parse_args()

    base_cfg = dict(
        T=100.0, dt=0.05, window_s=40.0, stride_s=2.0,
        gamma=0.5, n_real=80, mi_bins=12,
    )
    sigma_list = [1.0, 0.3, 0.1, 0.03, 0.01]
    seed = 11

    t0 = time.time()
    print(f"[sigma_scan] base_cfg={base_cfg} sigma_list={sigma_list} "
          f"seed={seed}", flush=True)
    runs = []
    for sigma in sigma_list:
        ts = time.time()
        cfg = dict(base_cfg)
        cfg["sigma"] = sigma
        M = build_M(seed=seed, **cfg)
        mu = M.mean(axis=0)
        sd = M.std(axis=0)
        Mc = M - mu
        v_null, S, _, _ = extract_v1(M)
        eigs = (S ** 2) / max(1, len(Mc) - 1)
        kbk = kbk_rank_gap(S)
        sub234, cos_a = project_234(v_null)
        elapsed = time.time() - ts
        r = {
            "sigma": sigma,
            "elapsed_s": elapsed,
            "n_windows": int(M.shape[0]),
            "operator_means": mu.tolist(),
            "operator_stds": sd.tolist(),
            "singular_values": S.tolist(),
            "eigenvalues_op_cov": eigs.tolist(),
            "eig_smallest": float(eigs[-1]),
            "eig_2nd": float(eigs[-2]),
            "eig_3rd": float(eigs[-3]),
            "eig_4th": float(eigs[-4]),
            "kbk_rank_gap": kbk,
            "extract_v1_v_null_234": sub234,
            "cos_to_artifact_+1+1+2": cos_a,
        }
        runs.append(r)
        print(f"  sigma={sigma:.3f}  elapsed={elapsed:.1f}s  "
              f"n_windows={M.shape[0]}", flush=True)
        print(f"    op_means = "
              f"{[f'{m:+.3f}' for m in mu]}", flush=True)
        print(f"    op_stds  = "
              f"{[f'{s:.3e}' for s in sd]}", flush=True)
        print(f"    eigenvalues = "
              f"{[f'{e:.3e}' for e in eigs]}", flush=True)
        print(f"    bottom-4 eigvals = "
              f"[{r['eig_4th']:.3e}, {r['eig_3rd']:.3e}, "
              f"{r['eig_2nd']:.3e}, {r['eig_smallest']:.3e}]",
              flush=True)
        print(f"    cos extract_v1 [2,3,4] to (+,+,+) = "
              f"{r['cos_to_artifact_+1+1+2']:+.4f}", flush=True)

    # log-log scaling fits
    sig = np.array(sigma_list)
    log_sig = np.log(sig)
    bottom_3 = ["eig_smallest", "eig_2nd", "eig_3rd"]
    fits = {}
    for k in bottom_3:
        log_eig = np.log(np.array([r[k] for r in runs]))
        slope, intercept = np.polyfit(log_sig, log_eig, 1)
        fits[k] = {
            "slope_loglog": float(slope),
            "intercept": float(intercept),
            "interpretation": f"{k} ~ sigma^{slope:.2f}",
        }
        print(f"\n  log-log slope: {k} ~ sigma^{slope:.2f} "
              f"(intercept {intercept:.2f}, i.e. coef "
              f"e^intercept = {np.exp(intercept):.2e})", flush=True)

    extrap_001 = {
        k: float(np.exp(fits[k]["intercept"])
                 * (0.001) ** fits[k]["slope_loglog"])
        for k in bottom_3
    }
    print(f"\n  extrapolation to sigma=0.001 (just for context):")
    for k, v in extrap_001.items():
        print(f"    {k} -> {v:.3e}", flush=True)

    results = {
        "base_cfg": base_cfg,
        "sigma_list": sigma_list,
        "seed": seed,
        "runs": runs,
        "scaling_fits_loglog": fits,
        "elapsed_total_s": time.time() - t0,
    }
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[sigma_scan] wrote {args.out} "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
