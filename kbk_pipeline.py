"""
KBK pipeline on extract_v1 operator output.

Stack: Kaiser-Brunton-Kutz 2024 (arXiv:2403.04889) SVD-rank-gap +
symbolic recovery pipeline applied to the extract_v1 protocol
(operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] over windowed
differential entropy of a coupled two-stream input).

This is a stacked methodology experiment per the "they never stacked"
rule. Gomez-Herrero 2015 had the ensemble + windowed-entropy
combinations on coupled circuits but did not assemble them into an
algebraic basis + SVD null extraction. Kaiser-Brunton-Kutz 2024 had
the rank-gap + symbolic recovery pipeline but did not apply it to
windowed differential entropy of coupled species. This script does
both.

What this script does:
  1. Generate baseline OU two-stream input (same protocol as
     extract_v1 baseline).
  2. Slide a window of length W across each realization; per window,
     compute Vasicek differential entropy estimates H_a, H_b and a
     histogram-MI estimate MI(a, b). Stack the operator vector
     [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] per window.
  3. Pool windows across N_REAL realizations into matrix M
     (n_windows x 6). Center columns. Run SVD.
  4. KBK step 1 -- rank-gap analysis. Report singular values and the
     largest gap-ratio between consecutive singular values. The rank
     k of the recovered subspace is the index of the largest ratio.
  5. KBK step 2 -- symbolic recovery. For each near-null right
     singular vector v_i, attempt a symbolic reading: sparsify with
     a small magnitude threshold, normalize, and print as a linear
     relation on the operator basis. Compute residual on the centered
     M.
  6. Compare KBK rank-recovered subspace to the extract_v1 smallest
     singular vector and the (+1,+1,+2)/sqrt(6) protocol artifact
     direction.

Output: kbk_pipeline_results.json
"""
import numpy as np
import json
import time
import argparse
import sys

OP_NAMES = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def vasicek_entropy(x, m=None):
    """Vasicek differential entropy estimator.

    H_hat = (1/n) sum_i log( (X_(i+m) - X_(i-m)) * n / (2m) )
    with edge-symmetric extension. Standard m = round(sqrt(n) + 0.5).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if m is None:
        m = max(1, int(np.round(np.sqrt(n) + 0.5)))
    xs = np.sort(x)
    idx_lo = np.arange(n) - m
    idx_hi = np.arange(n) + m
    idx_lo = np.clip(idx_lo, 0, n - 1)
    idx_hi = np.clip(idx_hi, 0, n - 1)
    diffs = xs[idx_hi] - xs[idx_lo]
    diffs = np.maximum(diffs, 1e-300)
    return float(np.mean(np.log(diffs * n / (2 * m))))


def mi_hist_2d(x, y, bins=12):
    """Histogram MI estimator (nats)."""
    h2, ex, ey = np.histogram2d(x, y, bins=bins, density=True)
    mask = h2 > 0
    if not mask.any():
        return 0.0
    dA = (ex[1] - ex[0]) * (ey[1] - ey[0])
    hxy = -np.sum(h2[mask] * np.log(h2[mask])) * dA
    hx_marg = np.sum(h2, axis=1) * (ey[1] - ey[0])
    hy_marg = np.sum(h2, axis=0) * (ex[1] - ex[0])
    hx_marg = hx_marg[hx_marg > 0]
    hy_marg = hy_marg[hy_marg > 0]
    hx = -np.sum(hx_marg * np.log(hx_marg)) * (ex[1] - ex[0])
    hy = -np.sum(hy_marg * np.log(hy_marg)) * (ey[1] - ey[0])
    return float(max(0.0, hx + hy - hxy))


def simulate_ou_pair(T, dt, gamma, sigma, rng, rho=0.0):
    """OU two-stream a, b. rho is correlation between driving noises."""
    n = int(T / dt)
    a = np.zeros(n)
    b = np.zeros(n)
    if rho != 0.0:
        L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    else:
        L = np.eye(2)
    for t in range(1, n):
        z = rng.standard_normal(2)
        z = L @ z
        a[t] = a[t - 1] + (-gamma * a[t - 1]) * dt + sigma * np.sqrt(dt) * z[0]
        b[t] = b[t - 1] + (-gamma * b[t - 1]) * dt + sigma * np.sqrt(dt) * z[1]
    return a, b


def compute_operator_matrix(
    a, b, window_samples, stride, mi_bins=12
):
    """Slide window across (a, b). Per window compute operator vector.

    Returns array shape (n_windows, 6).
    """
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
        Ha = vasicek_entropy(aw)
        Hb = vasicek_entropy(bw)
        MI = mi_hist_2d(aw, bw, bins=mi_bins)
        M[w, 0] = Ha
        M[w, 1] = Hb
        M[w, 2] = Ha * Ha
        M[w, 3] = Hb * Hb
        M[w, 4] = Ha * Hb
        M[w, 5] = MI
    return M


def extract_v1(M):
    """Center columns, SVD, return smallest right singular vector and
    singular value spectrum.
    """
    mu = M.mean(axis=0)
    Mc = M - mu
    U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
    v_null = Vt[-1]
    # sign-canonicalize: largest |entry| positive
    j = int(np.argmax(np.abs(v_null)))
    if v_null[j] < 0:
        v_null = -v_null
    return v_null, S, mu, Mc


def kbk_rank_gap(S):
    """KBK 2024 rank-gap analysis.

    Compute ratios S[i] / S[i+1] for descending S. Largest ratio
    indicates the rank cutoff: spectrum is rank-k where k is the
    index of the largest ratio (1-indexed from the top).
    """
    S = np.array(S, dtype=float)
    ratios = []
    for i in range(len(S) - 1):
        denom = S[i + 1] if S[i + 1] > 0 else 1e-300
        ratios.append(float(S[i] / denom))
    ratios = np.array(ratios)
    k_signal = int(np.argmax(ratios) + 1)
    rank_null = len(S) - k_signal
    return {
        "singular_values": S.tolist(),
        "ratios_consecutive": ratios.tolist(),
        "k_signal_kbk": k_signal,
        "rank_null_subspace_kbk": rank_null,
        "largest_gap_ratio": float(ratios.max()),
        "largest_gap_index": int(np.argmax(ratios)),
    }


def kbk_symbolic(Mc, Vt, n_null_keep, sparsify_thresh=0.15):
    """For each of the n_null_keep smallest right singular vectors,
    sparsify (zero coefficients below thresh*max(|coef|)), renormalize,
    print as an equation, and compute residual std on Mc @ v.

    Returns list of dicts.
    """
    results = []
    for k in range(n_null_keep):
        v = Vt[-(k + 1)].copy()
        # sign-canonicalize
        j = int(np.argmax(np.abs(v)))
        if v[j] < 0:
            v = -v
        v_full = v.copy()
        thresh = sparsify_thresh * np.max(np.abs(v))
        v_sparse = v.copy()
        v_sparse[np.abs(v_sparse) < thresh] = 0.0
        norm_sparse = np.linalg.norm(v_sparse)
        if norm_sparse > 0:
            v_sparse /= norm_sparse
        # residual std before sparsifying
        resid_full = Mc @ v_full
        resid_sparse = Mc @ v_sparse if norm_sparse > 0 else np.zeros(len(Mc))
        # symbolic eq string
        terms_full = [
            f"{v_full[i]:+.4f}*{OP_NAMES[i]}" for i in range(len(v_full))
        ]
        terms_sparse = [
            f"{v_sparse[i]:+.4f}*{OP_NAMES[i]}"
            for i in range(len(v_sparse))
            if v_sparse[i] != 0.0
        ]
        results.append({
            "null_index": k,
            "v_full": v_full.tolist(),
            "v_sparse": v_sparse.tolist(),
            "eq_full": " ".join(terms_full) + " ~ 0",
            "eq_sparse": (
                " ".join(terms_sparse) + " ~ 0"
                if terms_sparse else "(all-zero after sparsify)"
            ),
            "resid_std_full": float(np.std(resid_full)),
            "resid_std_sparse": float(np.std(resid_sparse)),
            "resid_mean_full": float(np.mean(resid_full)),
            "data_std_per_col": [float(np.std(Mc[:, i])) for i in range(Mc.shape[1])],
        })
    return results


def project_234(v):
    """Project 6D operator vector to coords [2, 3, 4]
    (H_a^2, H_b^2, H_a*H_b), normalize, sign-canonicalize.
    """
    sub = v[[2, 3, 4]].copy()
    n = np.linalg.norm(sub)
    if n == 0:
        return sub.tolist(), None
    sub = sub / n
    j = int(np.argmax(np.abs(sub)))
    if sub[j] < 0:
        sub = -sub
    target = np.array([1.0, 1.0, 2.0]) / np.sqrt(6.0)
    cos = float(np.dot(sub, target))
    return sub.tolist(), cos


def run_one(
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
    v_null, S, mu, Mc = extract_v1(M)
    kbk = kbk_rank_gap(S)
    n_keep = kbk["rank_null_subspace_kbk"]
    n_keep = max(1, min(n_keep, M.shape[1]))
    syms = kbk_symbolic(Mc, np.linalg.svd(Mc, full_matrices=False)[2], n_keep)
    sub234, cos_to_artifact = project_234(v_null)
    return {
        "seed": seed,
        "n_windows_total": int(M.shape[0]),
        "operator_means": mu.tolist(),
        "extract_v1_v_null_6d": v_null.tolist(),
        "extract_v1_v_null_234_normalized": sub234,
        "cos_extract_v1_to_artifact_+1+1+2": cos_to_artifact,
        "kbk_rank_gap": kbk,
        "kbk_symbolic_relations": syms,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true",
                   help="Small run for verification (2 min budget)")
    p.add_argument("--seeds", type=int, nargs="+", default=None)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        cfg = dict(
            T=20.0, dt=0.05, window_s=4.0, stride_s=0.5,
            n_real=20, gamma=0.5, sigma=1.0, mi_bins=10,
        )
        seeds = args.seeds if args.seeds else [11, 22, 33]
        out = args.out or "kbk_pipeline_canary.json"
    else:
        cfg = dict(
            T=100.0, dt=0.05, window_s=40.0, stride_s=1.0,
            n_real=500, gamma=0.5, sigma=1.0, mi_bins=12,
        )
        seeds = args.seeds if args.seeds else [11, 22, 33]
        out = args.out or "kbk_pipeline_results.json"

    t0 = time.time()
    print(f"[kbk_pipeline] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={seeds}", flush=True)

    runs = []
    for s in seeds:
        ts = time.time()
        r = run_one(s, **cfg)
        elapsed = time.time() - ts
        runs.append(r)
        v6 = r["extract_v1_v_null_6d"]
        v3 = r["extract_v1_v_null_234_normalized"]
        kbk = r["kbk_rank_gap"]
        print(f"  seed={s}  elapsed={elapsed:.1f}s  "
              f"n_windows={r['n_windows_total']}", flush=True)
        print(f"    S (singular values) = "
              f"{[f'{s:.3e}' for s in kbk['singular_values']]}",
              flush=True)
        print(f"    KBK rank_null_subspace = "
              f"{kbk['rank_null_subspace_kbk']}  "
              f"largest_gap_ratio = {kbk['largest_gap_ratio']:.3e}  "
              f"at index {kbk['largest_gap_index']}",
              flush=True)
        print(f"    extract_v1 v_null [2,3,4] normalized = "
              f"[{v3[0]:+.4f}, {v3[1]:+.4f}, {v3[2]:+.4f}]  "
              f"cos to (+1,+1,+2)/sqrt6 = "
              f"{r['cos_extract_v1_to_artifact_+1+1+2']:+.4f}",
              flush=True)
        for sym in r["kbk_symbolic_relations"]:
            print(f"    null[{sym['null_index']}] sparse: "
                  f"{sym['eq_sparse']}  resid_std_sparse="
                  f"{sym['resid_std_sparse']:.3e}",
                  flush=True)

    results = {
        "config": cfg,
        "seeds": seeds,
        "runs": runs,
        "elapsed_total_s": time.time() - t0,
    }

    # cross-seed summary
    cos_artifact = [r["cos_extract_v1_to_artifact_+1+1+2"] for r in runs]
    k_signals = [r["kbk_rank_gap"]["k_signal_kbk"] for r in runs]
    rank_nulls = [r["kbk_rank_gap"]["rank_null_subspace_kbk"] for r in runs]
    results["summary"] = {
        "cos_to_artifact_per_seed": cos_artifact,
        "cos_to_artifact_mean": float(np.mean(cos_artifact)),
        "cos_to_artifact_std": float(np.std(cos_artifact)),
        "kbk_k_signal_per_seed": k_signals,
        "kbk_rank_null_per_seed": rank_nulls,
        "kbk_rank_null_consistent": len(set(rank_nulls)) == 1,
    }

    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"[kbk_pipeline] wrote {out} elapsed={time.time()-t0:.1f}s",
          flush=True)


if __name__ == "__main__":
    main()
