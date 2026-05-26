"""
SINDy-style sparse symbolic regression on an extended polynomial
operator library, applied to the same OU baseline used by
kbk_pipeline.py.

Stack: Brunton-Proctor-Kutz 2016 SINDy (PNAS 113:3932) and its
sparse-null-direction variant (arXiv:2409.04463, arXiv:2211.05918)
applied to a windowed-differential-entropy feature library of
coupled species. Per the v5 literature scan, no published case of
this exists.

Question: when we extend the operator library beyond [H_a, H_b,
H_a^2, H_b^2, H_a*H_b, MI] to include degree-3 monomials
[H_a^3, H_b^3, H_a^2*H_b, H_a*H_b^2], can sparse regression
discover (a) cleaner symbolic statements of the three known
first-order relations, (b) genuinely new higher-order algebraic
relations, or (c) does it just rediscover the known three at
the same residual?

Method: Sequential Threshold Least Squares (STLS, SINDy core).
Pool centered operator data, run SVD to get null directions,
then for each near-null direction sparsify by iteratively
zeroing small coefficients and refitting on remaining nonzeros.
Compare to KBK's sparsify-by-magnitude approach.

Three sweeps:
  A. Degree-2 library (the 6-op KBK basis, control)
  B. Degree-3 library (10 ops: above + H_a^3, H_b^3, H_a^2*H_b,
     H_a*H_b^2)
  C. Degree-3 library with log/sqrt features added
     [log(MI+eps), sqrt(H_a^2 + H_b^2)] for sanity check on
     whether non-polynomial features help

Output: sindy_symbolic_results.json
"""
import numpy as np
import json
import time
import argparse
from kbk_pipeline import (
    simulate_ou_pair, vasicek_entropy, mi_hist_2d, project_234,
)


def build_M_base(seed, T, dt, window_s, stride_s, n_real,
                  gamma, sigma, mi_bins):
    """Original 6-op basis."""
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


def extend_library(M_base, mode="deg3"):
    """Extend the 6-op matrix to a richer library.

    mode='deg2':  return as is, 6 features.
    mode='deg3':  add H_a^3, H_b^3, H_a^2*H_b, H_a*H_b^2 (10 features).
    mode='deg3+nonpoly': add log(MI+eps), sqrt(H_a^2+H_b^2) (12).
    """
    Ha = M_base[:, 0]
    Hb = M_base[:, 1]
    Ha2 = M_base[:, 2]
    Hb2 = M_base[:, 3]
    HaHb = M_base[:, 4]
    MI = M_base[:, 5]
    names = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]
    cols = [Ha, Hb, Ha2, Hb2, HaHb, MI]
    if mode in ("deg3", "deg3+nonpoly"):
        Ha3 = Ha2 * Ha
        Hb3 = Hb2 * Hb
        Ha2Hb = Ha2 * Hb
        HaHb2 = Ha * Hb2
        cols += [Ha3, Hb3, Ha2Hb, HaHb2]
        names += ["H_a^3", "H_b^3", "H_a^2*H_b", "H_a*H_b^2"]
    if mode == "deg3+nonpoly":
        eps = 1e-8
        logMI = np.log(np.maximum(MI, 0) + eps)
        rH = np.sqrt(Ha2 + Hb2)
        cols += [logMI, rH]
        names += ["log(MI+eps)", "sqrt(H_a^2+H_b^2)"]
    M = np.column_stack(cols)
    return M, names


def stls_null_direction(Mc, v_init, thresh_frac=0.05, max_iter=20):
    """Sequential Threshold Least Squares to sparsify a null
    direction v.

    Iteratively:
      1. zero out entries of v with |v_i| < thresh_frac * max(|v|)
      2. on the remaining support, solve min ||Mc @ v|| s.t.
         ||v||=1 by taking the smallest right singular vector of
         Mc[:, support]
      3. embed back into full v (zeros outside support)
      4. repeat until support stabilizes or max_iter

    Returns sparsified v.
    """
    v = v_init.copy()
    p = len(v)
    prev_support = None
    for it in range(max_iter):
        thresh = thresh_frac * np.max(np.abs(v))
        support = np.where(np.abs(v) >= thresh)[0]
        if len(support) == 0:
            break
        if prev_support is not None and np.array_equal(support, prev_support):
            break
        prev_support = support
        # SVD on support only
        U, S, Vt = np.linalg.svd(Mc[:, support], full_matrices=False)
        v_sub = Vt[-1]
        v_new = np.zeros(p)
        v_new[support] = v_sub
        # sign-canonicalize
        j = int(np.argmax(np.abs(v_new)))
        if v_new[j] < 0:
            v_new = -v_new
        v = v_new
    return v


def find_null_directions(Mc, n_keep, stls_thresh=0.05):
    """Take the smallest-SV directions of Mc, sparsify each via STLS."""
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    out = []
    for k in range(n_keep):
        v_full = Vt[-(k + 1)].copy()
        j = int(np.argmax(np.abs(v_full)))
        if v_full[j] < 0:
            v_full = -v_full
        v_sparse = stls_null_direction(Mc, v_full, thresh_frac=stls_thresh)
        resid_full = Mc @ v_full
        resid_sparse = Mc @ v_sparse
        out.append({
            "null_index": k,
            "sv": float(S[-(k + 1)]),
            "v_full": v_full.tolist(),
            "v_sparse": v_sparse.tolist(),
            "support_size_full": int(np.sum(np.abs(v_full) > 1e-12)),
            "support_size_sparse": int(np.sum(np.abs(v_sparse) > 1e-12)),
            "resid_std_full": float(np.std(resid_full)),
            "resid_std_sparse": float(np.std(resid_sparse)),
            "norm_v_sparse": float(np.linalg.norm(v_sparse)),
        })
    return out


def format_eq(v, names, prec=3):
    terms = [f"{v[i]:+.{prec}f}*{names[i]}"
             for i in range(len(v)) if abs(v[i]) > 1e-12]
    if not terms:
        return "0 ~ 0"
    return " ".join(terms) + " ~ 0"


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
        out = args.out or "sindy_symbolic_canary.json"
    else:
        cfg = dict(T=100.0, dt=0.05, window_s=40.0, stride_s=1.0,
                   n_real=500, gamma=0.5, sigma=1.0, mi_bins=12)
        seeds = args.seeds if args.seeds else [11, 22, 33]
        out = args.out or "sindy_symbolic_results.json"

    modes = ["deg2", "deg3", "deg3+nonpoly"]
    t0 = time.time()
    print(f"[sindy] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={seeds} libs={modes}", flush=True)

    runs = []
    for seed in seeds:
        for mode in modes:
            ts = time.time()
            M_base = build_M_base(seed, **cfg)
            M, names = extend_library(M_base, mode=mode)
            mu = M.mean(axis=0)
            Mc = M - mu
            U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
            # take bottom min(rank_signal, 6) directions for inspection
            n_keep = min(6, M.shape[1] - 1)
            nulls = find_null_directions(Mc, n_keep, stls_thresh=0.05)
            elapsed = time.time() - ts
            run = {
                "seed": seed,
                "library": mode,
                "n_features": len(names),
                "feature_names": names,
                "n_windows": int(M.shape[0]),
                "singular_values": S.tolist(),
                "elapsed_s": elapsed,
                "nulls": nulls,
            }
            runs.append(run)
            print(f"  seed={seed}  lib={mode} ({len(names)} features)  "
                  f"elapsed={elapsed:.1f}s  n_windows={M.shape[0]}",
                  flush=True)
            print(f"    SV spectrum = "
                  f"{[f'{s:.3e}' for s in S]}", flush=True)
            for n in nulls[:4]:
                eq = format_eq(n["v_sparse"], names, prec=3)
                print(f"    null[{n['null_index']}] sv={n['sv']:.3e}  "
                      f"support_size={n['support_size_sparse']}/"
                      f"{len(names)}  "
                      f"resid_std={n['resid_std_sparse']:.3e}",
                      flush=True)
                print(f"       {eq}", flush=True)

    results = {"config": cfg, "seeds": seeds, "libraries": modes,
               "runs": runs, "elapsed_total_s": time.time() - t0}
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[sindy] wrote {out} elapsed={time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
