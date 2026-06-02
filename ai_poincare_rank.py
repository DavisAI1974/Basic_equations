"""
AI Poincare style numerical-rank-of-manifold (g) on extract_v1
baseline operator output.

Stack: Cao-Liu-Tegmark 2021 "AI Poincare" (arXiv:2011.04698) -- the
core observation is that the intrinsic dimension of the data cloud
in operator space tells you 6 - dim conserved quantities (algebraic
relations holding on the data). This is an INDEPENDENT method to
KBK's rank-gap singular-value analysis. If KBK and AI Poincare
agree, the rank reading is robust to method. If they disagree, that
itself is information.

Per "they never stacked" -- Cao-Liu-Tegmark applied local-PCA rank
to phase-space trajectories. We are stacking it onto the windowed-
entropy operator basis.

Three intrinsic-dimension estimators applied to the same baseline
operator matrix that kbk_pipeline.py produced:

  1. Two-NN (Facco et al. 2017). For each point, compute r1, r2.
     mu_i = r2/r1. F(mu) ~ 1 - mu^(-d). Fit d by linear regression
     of -log(1-F(mu)) vs log(mu) over the central quantiles.

  2. Local PCA rank-vs-scale curve. For each anchor point, take its
     k-NN at varying k. Run PCA on the neighborhood. Report number
     of components needed to explain >= 95% of local variance. Mean
     and median over anchors at each k. AI Poincare reads the rank
     as the plateau of this curve.

  3. Levina-Bickel MLE (2005). Maximum-likelihood estimator of
     intrinsic dimension from k-NN distances.

Same baseline as kbk_pipeline.py: OU two-stream, T=100, dt=0.05,
window=40s, stride=1s, gamma=0.5, sigma=1, n_real=500, 3 seeds.

Output: ai_poincare_rank_results.json
"""
import numpy as np
import json
import time
import argparse

from kbk_pipeline import (
    simulate_ou_pair,
    compute_operator_matrix,
    OP_NAMES,
)


def two_nn_dim(X, q_low=0.1, q_high=0.9):
    """Two-NN intrinsic dimension estimator (Facco et al. 2017).

    For each point compute distance to nearest (r1) and second
    nearest (r2). mu = r2/r1 >= 1. CDF of mu: F(mu) = 1 - mu^(-d).
    Linear regression: -log(1 - F(mu)) = d * log(mu), fit on central
    quantiles q_low to q_high to be robust to noise tails.
    """
    n = X.shape[0]
    # subsample for speed if huge
    if n > 5000:
        idx = np.random.default_rng(0).choice(n, 5000, replace=False)
        X = X[idx]
        n = 5000
    # all pairwise distances would be n^2; use a single-pass nearest
    # and second-nearest via sort of each row
    mus = np.zeros(n)
    for i in range(n):
        d = np.linalg.norm(X - X[i], axis=1)
        d_sorted = np.sort(d)
        r1 = d_sorted[1]
        r2 = d_sorted[2]
        if r1 <= 0:
            mus[i] = np.nan
            continue
        mus[i] = r2 / r1
    mus = mus[~np.isnan(mus)]
    mus = mus[mus > 1.0]
    mus_sorted = np.sort(mus)
    F = np.arange(1, len(mus_sorted) + 1) / (len(mus_sorted) + 1.0)
    # fit on central quantiles
    mask = (F >= q_low) & (F <= q_high)
    x = np.log(mus_sorted[mask])
    y = -np.log(1.0 - F[mask])
    # linear regression through origin: y = d*x
    d_est = float(np.dot(x, y) / np.dot(x, x))
    return d_est, len(mus)


def local_pca_rank_curve(X, k_list, var_thresh=0.95, n_anchors=200,
                          rng=None):
    """For each k in k_list, sample n_anchors anchor points, take k-NN
    around each, run PCA, count components for cumulative variance
    >= var_thresh. Return per-k mean and median rank.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    n = X.shape[0]
    anchors = rng.choice(n, size=min(n_anchors, n), replace=False)
    results = {}
    for k in k_list:
        ranks = []
        for a in anchors:
            d = np.linalg.norm(X - X[a], axis=1)
            idx = np.argsort(d)[:k]
            nbrs = X[idx]
            nbrs_c = nbrs - nbrs.mean(axis=0)
            U, S, Vt = np.linalg.svd(nbrs_c, full_matrices=False)
            var = S ** 2 / max(1, len(idx) - 1)
            var_cum = np.cumsum(var) / np.sum(var)
            r = int(np.searchsorted(var_cum, var_thresh) + 1)
            ranks.append(r)
        ranks = np.array(ranks)
        results[int(k)] = {
            "mean": float(ranks.mean()),
            "median": float(np.median(ranks)),
            "std": float(ranks.std()),
            "rank_hist": np.bincount(
                ranks, minlength=X.shape[1] + 1
            ).tolist(),
        }
    return results


def levina_bickel_mle(X, k=10):
    """Levina-Bickel maximum-likelihood intrinsic dimension.

    d_hat(x) = [ (1/(k-1)) sum_{j=1}^{k-1} log( T_k(x) / T_j(x) ) ]^-1
    where T_j(x) is the distance to the j-th nearest neighbor.
    Final d_hat is the mean over anchors.
    """
    n = X.shape[0]
    if n > 3000:
        idx = np.random.default_rng(1).choice(n, 3000, replace=False)
        X = X[idx]
        n = 3000
    d_hats = []
    for i in range(n):
        d = np.linalg.norm(X - X[i], axis=1)
        d_sorted = np.sort(d)[1:k + 1]
        if d_sorted[0] <= 0:
            continue
        log_ratios = np.log(d_sorted[-1] / d_sorted[:-1])
        if len(log_ratios) == 0:
            continue
        denom = np.mean(log_ratios)
        if denom <= 0:
            continue
        d_hats.append(1.0 / denom)
    d_hats = np.array(d_hats)
    return float(d_hats.mean()), float(d_hats.std()), len(d_hats)


def build_operator_matrix(
    seed, T, dt, window_s, stride_s, n_real, gamma, sigma, mi_bins
):
    rng = np.random.default_rng(seed)
    window_samples = int(window_s / dt)
    stride_samples = max(1, int(stride_s / dt))
    M_pool = []
    for r in range(n_real):
        a, b = simulate_ou_pair(T, dt, gamma, sigma, rng)
        M_r = compute_operator_matrix(
            a, b, window_samples, stride_samples, mi_bins=mi_bins
        )
        if len(M_r) > 0:
            M_pool.append(M_r)
    M = np.concatenate(M_pool, axis=0)
    return M


def run_one(seed, **cfg):
    M = build_operator_matrix(seed, **cfg)
    mu = M.mean(axis=0)
    sd = M.std(axis=0) + 1e-12
    Xs = (M - mu) / sd  # z-score so each operator contributes
                       # comparable distance scale
    Xc = M - mu
    out = {"seed": seed, "n_points": int(M.shape[0])}
    # two-NN on z-scored
    d_2nn_z, n_z = two_nn_dim(Xs)
    d_2nn_c, n_c = two_nn_dim(Xc)
    out["two_nn"] = {
        "intrinsic_dim_zscore": d_2nn_z,
        "intrinsic_dim_centered_raw_units": d_2nn_c,
        "n_used": n_z,
    }
    # local PCA at multiple scales (z-scored data)
    k_list = [10, 25, 50, 100, 200, 500]
    out["local_pca"] = {
        "var_thresh": 0.95,
        "per_k": local_pca_rank_curve(Xs, k_list, var_thresh=0.95,
                                       n_anchors=200),
    }
    # Levina-Bickel MLE
    d_lb, std_lb, n_lb = levina_bickel_mle(Xs, k=10)
    out["levina_bickel_k10"] = {
        "intrinsic_dim": d_lb, "std": std_lb, "n_used": n_lb,
    }
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=None)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        cfg = dict(T=20.0, dt=0.05, window_s=4.0, stride_s=0.5,
                   n_real=20, gamma=0.5, sigma=1.0, mi_bins=10)
        seeds = args.seeds if args.seeds else [11]
        out = args.out or "ai_poincare_rank_canary.json"
    else:
        cfg = dict(T=100.0, dt=0.05, window_s=40.0, stride_s=1.0,
                   n_real=500, gamma=0.5, sigma=1.0, mi_bins=12)
        seeds = args.seeds if args.seeds else [11, 22, 33]
        out = args.out or "ai_poincare_rank_results.json"

    t0 = time.time()
    print(f"[ai_poincare] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={seeds}", flush=True)

    runs = []
    for s in seeds:
        ts = time.time()
        r = run_one(s, **cfg)
        runs.append(r)
        elapsed = time.time() - ts
        print(f"  seed={s} n_points={r['n_points']} "
              f"elapsed={elapsed:.1f}s", flush=True)
        print(f"    two-NN dim (z-scored)       = "
              f"{r['two_nn']['intrinsic_dim_zscore']:.3f}", flush=True)
        print(f"    two-NN dim (raw centered)   = "
              f"{r['two_nn']['intrinsic_dim_centered_raw_units']:.3f}",
              flush=True)
        print(f"    Levina-Bickel MLE (k=10)    = "
              f"{r['levina_bickel_k10']['intrinsic_dim']:.3f} "
              f"+/- {r['levina_bickel_k10']['std']:.3f}",
              flush=True)
        print(f"    local-PCA 95%-var rank per k:", flush=True)
        for k, v in r["local_pca"]["per_k"].items():
            print(f"       k={k:>3d}  mean={v['mean']:.2f}  "
                  f"median={v['median']:.1f}  std={v['std']:.2f}",
                  flush=True)

    results = {"config": cfg, "seeds": seeds, "runs": runs,
               "elapsed_total_s": time.time() - t0}
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[ai_poincare] wrote {out} elapsed={time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
