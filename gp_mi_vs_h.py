"""
Session 7 / item (1): GP regression of MI vs H_a (and H_b, and joint)
per domain. Tests v5's "biology MI ~ polynomial(H_a) with R^2 0.66"
against a flexible GP (RBF + WhiteKernel) and checks if H_b is
really absent in biology by comparing GP(H_a) vs GP(H_a, H_b).

Data source: per_domain_kbk.simulate_* and build_ensemble_operator_matrix
(same simulators that produced the Session 6 per-domain claims).

For each domain x seed:
  - Time series of [H_a, H_b, MI] across ensemble at each t (post burn-in)
  - Fit polynomial deg 1, 2, 3 of MI vs H_a (closed-form via lstsq)
  - Fit GP MI vs H_a (RBF + White)
  - Fit GP MI vs H_b
  - Fit GP MI vs (H_a, H_b)  -- joint
  - Out-of-sample R^2 via block 5-fold CV (time series, contiguous blocks)

Readings to look at WITHOUT pre-assigning meaning (Rule on no pre-
assigned outcomes):
  - In-sample vs out-of-sample R^2 for each fit
  - Does GP beat polynomial deg-3 by more than block-CV noise?
  - Does joint (H_a, H_b) beat single H_a?
  - Cross-domain: which domain has the tightest MI-vs-H functional form?

Output: gp_mi_vs_h_results.json + console log
"""
import argparse
import json
import time
import warnings
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF, ConstantKernel, WhiteKernel,
)
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score

from per_domain_kbk import (
    DOMAINS, build_ensemble_operator_matrix,
)


warnings.filterwarnings("ignore", category=UserWarning)


def block_kfold_indices(n, n_splits=5):
    """Contiguous-block k-fold for time series.

    Returns list of (train_idx, test_idx) tuples.
    """
    fold_size = n // n_splits
    folds = []
    for k in range(n_splits):
        lo = k * fold_size
        hi = (k + 1) * fold_size if k < n_splits - 1 else n
        test = np.arange(lo, hi)
        train = np.concatenate([np.arange(0, lo), np.arange(hi, n)])
        folds.append((train, test))
    return folds


def random_kfold_indices(n, n_splits=5, seed=0):
    """Random k-fold (shuffled).

    Note: for time series with temporal correlation, random k-fold
    overestimates true generalization because nearby (correlated)
    samples can land in different folds. We use it here to isolate
    'functional-form' fit from 'stationarity-across-time' confound
    that the block-CV variant carries. Both readings reported.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    fold_size = n // n_splits
    folds = []
    for k in range(n_splits):
        lo = k * fold_size
        hi = (k + 1) * fold_size if k < n_splits - 1 else n
        test = perm[lo:hi]
        train = np.concatenate([perm[:lo], perm[hi:]])
        folds.append((train, test))
    return folds


def poly_fit_eval(X, y, deg, folds_block, folds_rand):
    """Polynomial (in 1 or multi feature) via Vandermonde / monomials,
    closed-form lstsq. Returns in-sample R^2, block-CV R^2, random-CV
    R^2, per-fold values, coefficients.

    X shape (n, d). If d==1: monomials H_a^0..H_a^deg.
    If d==2: full polynomial up to total degree `deg` (lex).
    """
    n, d = X.shape

    def design(Xq):
        if d == 1:
            return np.vstack([Xq[:, 0] ** k for k in range(deg + 1)]).T
        elif d == 2:
            cols = []
            names = []
            for a in range(deg + 1):
                for b in range(deg + 1 - a):
                    cols.append((Xq[:, 0] ** a) * (Xq[:, 1] ** b))
                    names.append(f"x0^{a}*x1^{b}")
            return np.array(cols).T, names
        else:
            raise ValueError("d > 2 not supported here")

    if d == 1:
        Phi = design(X)
        name_list = [f"x0^{k}" for k in range(deg + 1)]
    else:
        Phi, name_list = design(X)

    # in-sample
    coef_full, *_ = np.linalg.lstsq(Phi, y, rcond=None)
    y_pred_full = Phi @ coef_full
    r2_in = r2_score(y, y_pred_full)

    def cv_score(folds):
        r2s = []
        for tr, te in folds:
            coef, *_ = np.linalg.lstsq(Phi[tr], y[tr], rcond=None)
            y_pred = Phi[te] @ coef
            r2s.append(r2_score(y[te], y_pred))
        return r2s

    r2_block = cv_score(folds_block)
    r2_rand = cv_score(folds_rand)

    return {
        "deg": deg,
        "r2_in_sample": float(r2_in),
        "r2_block_mean": float(np.mean(r2_block)),
        "r2_block_std":  float(np.std(r2_block)),
        "r2_block_per_fold": [float(r) for r in r2_block],
        "r2_rand_mean": float(np.mean(r2_rand)),
        "r2_rand_std":  float(np.std(r2_rand)),
        "r2_rand_per_fold": [float(r) for r in r2_rand],
        "coef": [float(c) for c in coef_full],
        "feature_names": name_list,
    }


def gp_fit_eval(X, y, folds_block, folds_rand, normalize_y=True,
                gp_max_n=500):
    """GP with RBF + WhiteKernel kernel, length_scale optimized.

    For n > gp_max_n, uniformly subsamples to gp_max_n before fitting
    (in-sample reported on subsample, OOS evaluates on held-out from
    the subsampled training set). This keeps the GP fit O(gp_max_n^3)
    bounded; polynomial fits in the parent function still use full n.

    Returns in-sample R^2, block-CV R^2, random-CV R^2, fitted
    hyperparameters, and the subsample size actually used.
    """
    d = X.shape[1]
    n = X.shape[0]
    n_used = min(n, gp_max_n)
    if n_used < n:
        # uniformly spaced subsample preserving temporal order
        idx = np.linspace(0, n - 1, n_used).astype(int)
        X = X[idx]; y = y[idx]
        # rebuild folds at the subsampled indices
        folds_block = block_kfold_indices(n_used, n_splits=5)
        # random folds: re-seed deterministically
        folds_rand = random_kfold_indices(n_used, n_splits=5, seed=0)

    # standardize X to make length-scale optimization sane across domains
    Xm = X.mean(axis=0, keepdims=True)
    Xs = X.std(axis=0, keepdims=True) + 1e-12
    Xn = (X - Xm) / Xs

    # build kernel: 1.0 * RBF(ls=1.0) + WhiteKernel
    kernel = (
        ConstantKernel(1.0, (1e-3, 1e3))
        * RBF(length_scale=np.ones(d), length_scale_bounds=(1e-2, 1e2))
        + WhiteKernel(noise_level=1e-2, noise_level_bounds=(1e-6, 1.0))
    )

    # in-sample
    gp = GaussianProcessRegressor(
        kernel=kernel, normalize_y=normalize_y,
        alpha=0.0, n_restarts_optimizer=2, random_state=0,
    )
    gp.fit(Xn, y)
    y_pred_full = gp.predict(Xn)
    r2_in = r2_score(y, y_pred_full)
    kernel_str = str(gp.kernel_)
    lml = float(gp.log_marginal_likelihood_value_)

    def cv_score(folds):
        r2s = []
        for tr, te in folds:
            gp_k = GaussianProcessRegressor(
                kernel=kernel, normalize_y=normalize_y,
                alpha=0.0, n_restarts_optimizer=2, random_state=0,
            )
            gp_k.fit(Xn[tr], y[tr])
            y_pred = gp_k.predict(Xn[te])
            r2s.append(r2_score(y[te], y_pred))
        return r2s

    r2_block = cv_score(folds_block)
    r2_rand = cv_score(folds_rand)

    return {
        "d_input": d,
        "n_used_for_gp": int(n_used),
        "r2_in_sample": float(r2_in),
        "r2_block_mean": float(np.mean(r2_block)),
        "r2_block_std":  float(np.std(r2_block)),
        "r2_block_per_fold": [float(r) for r in r2_block],
        "r2_rand_mean": float(np.mean(r2_rand)),
        "r2_rand_std":  float(np.std(r2_rand)),
        "r2_rand_per_fold": [float(r) for r in r2_rand],
        "kernel_fit": kernel_str,
        "log_marginal_likelihood": lml,
    }


def analyze_domain(name, simfn, seed, cfg, folds_n=5):
    """Per-domain GP and polynomial fits.

    Builds operator matrix, extracts H_a, H_b, MI columns, fits
    polynomial deg 1/2/3 and GP for (H_a), (H_b), (H_a, H_b).
    """
    X1, X2 = simfn(seed=seed, **cfg)
    M, _ = build_ensemble_operator_matrix(X1, X2)
    # columns of M: [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]
    Ha = M[:, 0]
    Hb = M[:, 1]
    MI = M[:, 5]
    n = len(MI)
    folds_block = block_kfold_indices(n, n_splits=folds_n)
    folds_rand = random_kfold_indices(n, n_splits=folds_n, seed=seed)

    X_a = Ha.reshape(-1, 1)
    X_b = Hb.reshape(-1, 1)
    X_ab = np.column_stack([Ha, Hb])

    res = {
        "name": name,
        "n_points": int(n),
        "Ha_stats": {"mean": float(Ha.mean()), "std": float(Ha.std()),
                     "min": float(Ha.min()), "max": float(Ha.max())},
        "Hb_stats": {"mean": float(Hb.mean()), "std": float(Hb.std()),
                     "min": float(Hb.min()), "max": float(Hb.max())},
        "MI_stats": {"mean": float(MI.mean()), "std": float(MI.std()),
                     "min": float(MI.min()), "max": float(MI.max())},
        "corr_MI_Ha": float(np.corrcoef(MI, Ha)[0, 1]),
        "corr_MI_Hb": float(np.corrcoef(MI, Hb)[0, 1]),
        "corr_Ha_Hb": float(np.corrcoef(Ha, Hb)[0, 1]),
    }

    res["poly_MI_vs_Ha"] = {}
    res["poly_MI_vs_Hb"] = {}
    for deg in (1, 2, 3):
        res["poly_MI_vs_Ha"][f"deg{deg}"] = poly_fit_eval(
            X_a, MI, deg, folds_block, folds_rand)
        res["poly_MI_vs_Hb"][f"deg{deg}"] = poly_fit_eval(
            X_b, MI, deg, folds_block, folds_rand)
    res["poly_MI_vs_Ha_Hb_deg3"] = poly_fit_eval(
        X_ab, MI, 3, folds_block, folds_rand)

    res["gp_MI_vs_Ha"]   = gp_fit_eval(X_a,  MI, folds_block, folds_rand)
    res["gp_MI_vs_Hb"]   = gp_fit_eval(X_b,  MI, folds_block, folds_rand)
    res["gp_MI_vs_HaHb"] = gp_fit_eval(X_ab, MI, folds_block, folds_rand)

    return res


def print_row(r):
    def tup(d, k):
        return (d["r2_in_sample"], d[f"r2_{k}_mean"], d[f"r2_{k}_std"])

    print(f"  {r['name']:22s}  n={r['n_points']}", flush=True)
    print(f"    MI stats: mean={r['MI_stats']['mean']:+.3f} "
          f"std={r['MI_stats']['std']:.3f}  "
          f"corr(MI,Ha)={r['corr_MI_Ha']:+.3f}  "
          f"corr(MI,Hb)={r['corr_MI_Hb']:+.3f}  "
          f"corr(Ha,Hb)={r['corr_Ha_Hb']:+.3f}",
          flush=True)
    print(f"    fmt: in / block-CV / random-CV", flush=True)
    for src in ("Ha", "Hb"):
        for deg in (1, 2, 3):
            d = r[f"poly_MI_vs_{src}"][f"deg{deg}"]
            ins, br, _ = tup(d, "block")
            _, rr, _ = tup(d, "rand")
            print(f"    poly{deg}(MI~{src}):   "
                  f"in={ins:+.3f}  block={br:+.3f}  rand={rr:+.3f}",
                  flush=True)
        g = r[f"gp_MI_vs_{src}"]
        ins, br, _ = tup(g, "block")
        _, rr, _ = tup(g, "rand")
        print(f"    GP(MI~{src}):     "
              f"in={ins:+.3f}  block={br:+.3f}  rand={rr:+.3f}",
              flush=True)
    d = r["poly_MI_vs_Ha_Hb_deg3"]
    ins, br, _ = tup(d, "block"); _, rr, _ = tup(d, "rand")
    print(f"    poly3(MI~Ha,Hb): in={ins:+.3f}  block={br:+.3f}  "
          f"rand={rr:+.3f}", flush=True)
    g = r["gp_MI_vs_HaHb"]
    ins, br, _ = tup(g, "block"); _, rr, _ = tup(g, "rand")
    print(f"    GP(MI~Ha,Hb):    in={ins:+.3f}  block={br:+.3f}  "
          f"rand={rr:+.3f}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11, 22])
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        out = args.out or "gp_mi_vs_h_canary.json"
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "gp_mi_vs_h_results.json"

    t0 = time.time()
    print(f"[gp_mi_vs_h] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={args.seeds}", flush=True)

    all_results = {}
    for seed in args.seeds:
        print(f"\n=== seed={seed} ===", flush=True)
        seed_results = []
        for name, simfn in DOMAINS.items():
            ts = time.time()
            r = analyze_domain(name, simfn, seed, cfg)
            r["elapsed_s"] = time.time() - ts
            seed_results.append(r)
            print_row(r)
            print(f"    elapsed: {r['elapsed_s']:.1f}s", flush=True)
        all_results[f"seed{seed}"] = seed_results

    blob = {
        "config": cfg, "seeds": args.seeds,
        "results_by_seed": all_results,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[gp_mi_vs_h] wrote {out}  elapsed={time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
