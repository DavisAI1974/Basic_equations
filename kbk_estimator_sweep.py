"""
Estimator-sensitivity probe on the KBK + AI Poincare stack.

Question: does the eigenvalue floor of the operator covariance
(machine epsilon in v5, ~1e-5 to 1e-4 in our Vasicek run) depend on
the entropy/MI estimator? If yes, the "rank 3" reading of v5 may be
estimator-specific. If no, the eigenvalue magnitude is structural
and the discrepancy is elsewhere (window length, dt, marginal).

Three estimators on the SAME OU baseline data, with the same
operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] and the same
extract_v1 + KBK + AI Poincare pipelines:

  1. Vasicek (sample-spacing). The baseline from kbk_pipeline.py.
     Closed-form deterministic function of sorted samples.

  2. Histogram-KDE (Gaussian KDE bandwidth via Silverman's rule).
     Smoother estimator. Should give smoother H trajectories,
     possibly tighter algebraic relations.

  3. KSG k=4 (Kraskov-Stoegbauer-Grassberger 2004) for MI;
     k-NN-based differential entropy (Singh et al. 2003) for H.
     Nonparametric, captures non-Gaussian structure cleanly.

Compare per estimator:
  - eigenvalue spectrum of operator covariance
  - KBK rank-gap reading (signal_rank, gap_ratio)
  - AI Poincare local-PCA 95%-var rank
  - extract_v1 smallest-SV direction
  - cos to (+1,+1,+2)/sqrt(6) artifact in [2,3,4] projection
  - sparse symbolic null relations

Per "they never stacked" -- the Kaiser-Brunton-Kutz pipeline and AI
Poincare have been applied to phase-space data with various
features but, per the v5 literature scan, never to windowed
differential entropy of coupled species and never compared across
entropy estimator families. This is the comparison.

Output: kbk_estimator_sweep_results.json
"""
import numpy as np
import json
import time
import argparse
from kbk_pipeline import (
    simulate_ou_pair, vasicek_entropy, mi_hist_2d, OP_NAMES,
    extract_v1, kbk_rank_gap, kbk_symbolic, project_234,
)
from ai_poincare_rank import (
    two_nn_dim, local_pca_rank_curve, levina_bickel_mle,
)


def kde_entropy(x, n_grid=256):
    """Gaussian KDE differential entropy with Silverman's rule.

    H = -integral p(x) log p(x) dx, evaluated on a grid covering
    +/- 4 sigma around the data.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    sd = max(np.std(x), 1e-12)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    A = min(sd, iqr / 1.34) if iqr > 0 else sd
    h = 0.9 * A * n ** (-1.0 / 5.0)
    h = max(h, 1e-6)
    lo = x.min() - 4 * h
    hi = x.max() + 4 * h
    grid = np.linspace(lo, hi, n_grid)
    dx = grid[1] - grid[0]
    # density: sum of Gaussians
    diffs = (grid[:, None] - x[None, :]) / h
    dens = np.exp(-0.5 * diffs ** 2).sum(axis=1) / (n * h * np.sqrt(2 * np.pi))
    mask = dens > 1e-300
    H = -float(np.sum(dens[mask] * np.log(dens[mask])) * dx)
    return H


def knn_entropy(x, k=4):
    """k-NN differential entropy (Kozachenko-Leonenko).

    H = log(n-1) - psi(k) + log(c_d) + (d/n) * sum_i log(r_i)
    with d=1, c_d = 2 (1D), and r_i is distance to k-th nearest.
    """
    from scipy.special import digamma
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    n = len(x)
    # 1D: distances = |x_i - x_j|
    sorted_x = np.sort(x.ravel())
    # for each point, distance to k-th nearest via sorted structure
    rks = np.zeros(n)
    for i in range(n):
        d = np.abs(sorted_x - sorted_x[i])
        d_sorted = np.sort(d)
        rks[i] = d_sorted[k]
    rks = np.maximum(rks, 1e-300)
    H = (
        np.log(n - 1) - digamma(k) + np.log(2.0)
        + np.mean(np.log(rks))
    )
    return float(H)


def ksg_mi(x, y, k=4):
    """KSG estimator (Kraskov-Stoegbauer-Grassberger 2004, alg 1).

    For 1D x, y. Uses Chebyshev distance in joint space.
    """
    from scipy.special import digamma
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = len(x)
    # joint Chebyshev distance to k-th neighbor
    # (slow O(n^2) — fine for n ~ few hundred)
    eps = np.zeros(n)
    for i in range(n):
        dx = np.abs(x - x[i])
        dy = np.abs(y - y[i])
        d = np.maximum(dx, dy)
        d_sorted = np.sort(d)
        eps[i] = d_sorted[k]
    # marginal counts strictly within eps_i (exclusive)
    nx = np.zeros(n)
    ny = np.zeros(n)
    for i in range(n):
        nx[i] = np.sum(np.abs(x - x[i]) < eps[i]) - 1
        ny[i] = np.sum(np.abs(y - y[i]) < eps[i]) - 1
    nx = np.maximum(nx, 1)
    ny = np.maximum(ny, 1)
    mi = (
        digamma(k) + digamma(n)
        - np.mean(digamma(nx + 1) + digamma(ny + 1))
    )
    return float(max(0.0, mi))


def compute_operator_matrix_with(
    a, b, window_samples, stride, entropy_fn, mi_fn
):
    n = len(a)
    n_w = (n - window_samples) // stride + 1
    if n_w <= 0:
        return np.zeros((0, 6))
    M = np.zeros((n_w, 6))
    for w in range(n_w):
        s = w * stride
        e = s + window_samples
        aw = a[s:e]
        bw = b[s:e]
        Ha = entropy_fn(aw)
        Hb = entropy_fn(bw)
        MI = mi_fn(aw, bw)
        M[w, 0] = Ha
        M[w, 1] = Hb
        M[w, 2] = Ha * Ha
        M[w, 3] = Hb * Hb
        M[w, 4] = Ha * Hb
        M[w, 5] = MI
    return M


ESTIMATORS = {
    "vasicek_hist": dict(
        entropy_fn=vasicek_entropy,
        mi_fn=lambda a, b: mi_hist_2d(a, b, bins=12),
        desc="Vasicek H + 12-bin histogram MI (baseline from "
             "kbk_pipeline.py)",
    ),
    "kde_hist": dict(
        entropy_fn=kde_entropy,
        mi_fn=lambda a, b: mi_hist_2d(a, b, bins=12),
        desc="Gaussian KDE H (Silverman) + 12-bin histogram MI",
    ),
    "knn_ksg": dict(
        entropy_fn=lambda x: knn_entropy(x, k=4),
        mi_fn=lambda a, b: ksg_mi(a, b, k=4),
        desc="k-NN H (Kozachenko-Leonenko k=4) + KSG MI k=4",
    ),
}


def run_one_estimator(
    name, fns, T, dt, window_s, stride_s, n_real, gamma, sigma, seed
):
    rng = np.random.default_rng(seed)
    window_samples = int(window_s / dt)
    stride_samples = max(1, int(stride_s / dt))
    M_pool = []
    for r in range(n_real):
        a, b = simulate_ou_pair(T, dt, gamma, sigma, rng)
        M_r = compute_operator_matrix_with(
            a, b, window_samples, stride_samples,
            fns["entropy_fn"], fns["mi_fn"],
        )
        if len(M_r) > 0:
            M_pool.append(M_r)
    M = np.concatenate(M_pool, axis=0)
    mu = M.mean(axis=0)
    sd = M.std(axis=0) + 1e-12
    Mc = M - mu

    # SVD + KBK rank gap
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    v_null = Vt[-1].copy()
    j = int(np.argmax(np.abs(v_null)))
    if v_null[j] < 0:
        v_null = -v_null
    kbk = kbk_rank_gap(S)
    eigvals = (S ** 2) / max(1, len(Mc) - 1)
    n_keep = max(1, min(kbk["rank_null_subspace_kbk"], 4))
    syms = kbk_symbolic(Mc, Vt, n_keep)
    sub234, cos_artifact = project_234(v_null)

    # AI Poincare style on z-scored
    Xs = (M - mu) / sd
    d_2nn, _ = two_nn_dim(Xs)
    d_lb, std_lb, _ = levina_bickel_mle(Xs, k=10)
    lpca = local_pca_rank_curve(
        Xs, k_list=[25, 50, 100, 200, 500],
        var_thresh=0.95, n_anchors=200,
    )

    return {
        "estimator": name,
        "desc": fns["desc"],
        "seed": seed,
        "n_windows": int(M.shape[0]),
        "operator_means": mu.tolist(),
        "operator_stds": sd.tolist(),
        "singular_values": S.tolist(),
        "eigenvalues_op_cov": eigvals.tolist(),
        "smallest_eigenvalue": float(eigvals[-1]),
        "eigenvalue_3rd_smallest": float(eigvals[-3])
            if len(eigvals) >= 3 else None,
        "kbk_rank_gap": kbk,
        "extract_v1_v_null_6d": v_null.tolist(),
        "extract_v1_v_null_234_normalized": sub234,
        "cos_extract_v1_to_artifact_+1+1+2": cos_artifact,
        "kbk_symbolic_relations": syms,
        "ai_poincare_two_nn_dim_z": d_2nn,
        "ai_poincare_levina_bickel": {
            "dim": d_lb, "std": std_lb,
        },
        "ai_poincare_local_pca_per_k": lpca,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        # smaller config for speed; KNN/KSG are O(window_size^2) which
        # for window_samples=80 is fine but n_real*window_size matters
        cfg = dict(
            T=20.0, dt=0.05, window_s=4.0, stride_s=0.5,
            n_real=15, gamma=0.5, sigma=1.0,
        )
        out = args.out or "kbk_estimator_sweep_canary.json"
    else:
        # full-ish: same T/window as baseline. Bring n_real down a bit
        # because KNN/KSG are O(n^2) per window. window_samples=800
        # would be too slow for KSG; we keep window=20s with smaller
        # n_per_window to make KSG tractable.
        # Same T as baseline (T=100s, window=40s -> 800 samples per
        # window). That's 800^2 = 640k ops per window per stat, ~OK
        # at modest n_real.
        cfg = dict(
            T=100.0, dt=0.05, window_s=40.0, stride_s=2.0,
            n_real=80, gamma=0.5, sigma=1.0,
        )
        out = args.out or "kbk_estimator_sweep_results.json"

    t0 = time.time()
    seed = args.seed
    print(f"[estimator_sweep] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seed={seed}", flush=True)

    runs = []
    for name, fns in ESTIMATORS.items():
        ts = time.time()
        print(f"  estimator={name}  ({fns['desc']})", flush=True)
        try:
            r = run_one_estimator(name, fns, seed=seed, **cfg)
            elapsed = time.time() - ts
            runs.append(r)
            print(f"    elapsed={elapsed:.1f}s  n_windows={r['n_windows']}",
                  flush=True)
            print(f"    SV spectrum = "
                  f"{[f'{s:.3e}' for s in r['singular_values']]}",
                  flush=True)
            print(f"    eigenvalues = "
                  f"{[f'{e:.3e}' for e in r['eigenvalues_op_cov']]}",
                  flush=True)
            print(f"    smallest eigenvalue = "
                  f"{r['smallest_eigenvalue']:.3e}", flush=True)
            print(f"    KBK signal_rank={r['kbk_rank_gap']['k_signal_kbk']}  "
                  f"null_rank={r['kbk_rank_gap']['rank_null_subspace_kbk']}  "
                  f"gap_ratio={r['kbk_rank_gap']['largest_gap_ratio']:.2e}",
                  flush=True)
            print(f"    cos extract_v1 [2,3,4] to (+,+,+) = "
                  f"{r['cos_extract_v1_to_artifact_+1+1+2']:+.4f}",
                  flush=True)
            print(f"    AI Poincare two-NN(z) = "
                  f"{r['ai_poincare_two_nn_dim_z']:.3f}  "
                  f"LB MLE = {r['ai_poincare_levina_bickel']['dim']:.3f}",
                  flush=True)
            print(f"    AI Poincare local-PCA 95% rank per k:",
                  flush=True)
            for k, v in r["ai_poincare_local_pca_per_k"].items():
                print(f"       k={k:>3d}  mean={v['mean']:.2f}  "
                      f"median={v['median']:.1f}  std={v['std']:.2f}",
                      flush=True)
            for sym in r["kbk_symbolic_relations"][:3]:
                print(f"    null[{sym['null_index']}] sparse: "
                      f"{sym['eq_sparse']}  resid_std="
                      f"{sym['resid_std_sparse']:.3e}", flush=True)
        except Exception as e:
            elapsed = time.time() - ts
            print(f"    FAILED after {elapsed:.1f}s: {type(e).__name__}: "
                  f"{e}", flush=True)
            runs.append({
                "estimator": name, "seed": seed,
                "error": f"{type(e).__name__}: {e}",
            })

    results = {"config": cfg, "seed": seed, "runs": runs,
               "elapsed_total_s": time.time() - t0}
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[estimator_sweep] wrote {out} elapsed={time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
