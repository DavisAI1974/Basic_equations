"""
Window-size scan on Vasicek + histogram-MI baseline.

Question: is the eigenvalue floor of the operator covariance a
sampling-noise floor that scales with the number of samples per
window? If yes, going to larger N_samples_per_window will drive
the floor toward machine epsilon, recovering v5's spectrum
[4e-13, 8e-13, 3e-11, 3e-4, 4e-4, 2e-3] from the same basis.

Holds T fixed (same realization length) and varies dt to control
N_samples_per_window. Tracks the bottom-3 eigenvalues of operator
covariance and the AI Poincare local-PCA rank.

If the eigenvalue floor scales as 1/N_samples_per_window^p for
some p, this isolates a sampling-noise origin for the "three at
machine epsilon" claim.

Output: window_size_scan_results.json
"""
import numpy as np
import json
import time
import argparse
from kbk_pipeline import (
    simulate_ou_pair, vasicek_entropy, mi_hist_2d,
    extract_v1, kbk_rank_gap, project_234,
)
from ai_poincare_rank import local_pca_rank_curve


def build_M(seed, T, dt, window_s, stride_s, n_real,
            gamma, sigma, mi_bins):
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
    M = np.concatenate(M_pool, axis=0)
    return M, window_samples


def run_one(seed, T, dt, window_s, stride_s, n_real, gamma, sigma,
             mi_bins):
    t0 = time.time()
    M, n_per_window = build_M(seed, T, dt, window_s, stride_s,
                                n_real, gamma, sigma, mi_bins)
    mu = M.mean(axis=0)
    sd = M.std(axis=0) + 1e-12
    Mc = M - mu
    v_null, S, _, _ = extract_v1(M)
    eigvals = (S ** 2) / max(1, len(Mc) - 1)
    kbk = kbk_rank_gap(S)
    Xs = (M - mu) / sd
    lpca = local_pca_rank_curve(
        Xs, k_list=[50, 200], var_thresh=0.95, n_anchors=100,
    )
    sub234, cos_a = project_234(v_null)
    return {
        "seed": seed,
        "dt": dt,
        "window_s": window_s,
        "n_samples_per_window": n_per_window,
        "n_windows_total": int(M.shape[0]),
        "elapsed_s": time.time() - t0,
        "eigenvalues_op_cov": eigvals.tolist(),
        "eig_smallest": float(eigvals[-1]),
        "eig_2nd_smallest": float(eigvals[-2]),
        "eig_3rd_smallest": float(eigvals[-3]),
        "eig_4th_smallest": float(eigvals[-4]),
        "kbk_rank_gap": kbk,
        "ai_poincare_local_pca": lpca,
        "extract_v1_234_normalized": sub234,
        "cos_to_artifact_+++": cos_a,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=str,
                   default="window_size_scan_results.json")
    args = p.parse_args()

    # Fix T, window_s, gamma, sigma. Vary dt to change
    # samples-per-window. Keep n_real moderate.
    base_cfg = dict(
        T=100.0, window_s=40.0, stride_s=4.0,
        gamma=0.5, sigma=1.0, mi_bins=12, n_real=80,
    )
    seed = 11
    # dt values targeted at samples-per-window from 200 to 4000
    dt_list = [0.20, 0.10, 0.05, 0.025, 0.01]
    # corresponds to N_samples_per_window = 200, 400, 800, 1600, 4000

    t0 = time.time()
    print(f"[window_size_scan] base_cfg={base_cfg} seed={seed}",
          flush=True)
    runs = []
    for dt in dt_list:
        cfg = dict(base_cfg)
        cfg["dt"] = dt
        ts = time.time()
        r = run_one(seed=seed, **cfg)
        runs.append(r)
        elapsed = time.time() - ts
        print(f"  dt={dt:.3f}  n_per_window={r['n_samples_per_window']}  "
              f"n_windows={r['n_windows_total']}  "
              f"elapsed={elapsed:.1f}s",
              flush=True)
        print(f"    bottom-4 eigenvalues = "
              f"[{r['eig_4th_smallest']:.3e}, "
              f"{r['eig_3rd_smallest']:.3e}, "
              f"{r['eig_2nd_smallest']:.3e}, "
              f"{r['eig_smallest']:.3e}]",
              flush=True)
        print(f"    KBK signal_rank={r['kbk_rank_gap']['k_signal_kbk']}  "
              f"null_rank={r['kbk_rank_gap']['rank_null_subspace_kbk']}  "
              f"gap_ratio={r['kbk_rank_gap']['largest_gap_ratio']:.2e}",
              flush=True)
        for k, v in r["ai_poincare_local_pca"].items():
            print(f"    AI Poincare local-PCA k={k:>3d} "
                  f"mean={v['mean']:.2f} median={v['median']:.1f}",
                  flush=True)
        print(f"    cos extract_v1 [2,3,4] to (+,+,+) = "
              f"{r['cos_to_artifact_+++']:+.4f}", flush=True)

    results = {
        "base_cfg": base_cfg,
        "dt_values": dt_list,
        "seed": seed,
        "runs": runs,
        "elapsed_total_s": time.time() - t0,
    }
    # scaling analysis
    n_w = [r["n_samples_per_window"] for r in runs]
    eig_min = [r["eig_smallest"] for r in runs]
    eig_3rd = [r["eig_3rd_smallest"] for r in runs]
    # fit log(eig) = a + b * log(n_w)
    log_n = np.log(np.array(n_w))
    log_eig_min = np.log(np.array(eig_min))
    log_eig_3rd = np.log(np.array(eig_3rd))
    p_min = np.polyfit(log_n, log_eig_min, 1)
    p_3rd = np.polyfit(log_n, log_eig_3rd, 1)
    results["scaling_fit"] = {
        "eig_smallest_vs_n_samples_per_window": {
            "slope_in_loglog": float(p_min[0]),
            "intercept": float(p_min[1]),
            "interpretation": (
                f"eig_smallest ~ n_per_window^{p_min[0]:.2f}"
            ),
        },
        "eig_3rd_smallest_vs_n_samples_per_window": {
            "slope_in_loglog": float(p_3rd[0]),
            "intercept": float(p_3rd[1]),
            "interpretation": (
                f"eig_3rd_smallest ~ n_per_window^{p_3rd[0]:.2f}"
            ),
        },
    }
    print(f"\n[window_size_scan] log-log slope of smallest eigval "
          f"vs n_samples_per_window: {p_min[0]:.3f}", flush=True)
    print(f"  log-log slope of 3rd-smallest eigval: "
          f"{p_3rd[0]:.3f}", flush=True)
    with open(args.out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[window_size_scan] wrote {args.out} "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
