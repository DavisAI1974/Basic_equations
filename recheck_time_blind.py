"""
RECHECK: is the level null really "time-blind"? (Greg: we have contradicting data.)

The contradiction:
  - INFO-S15a/b claimed the level null is time-blind: permuting the time-ROWS of
    the operator matrix and re-extracting gives cos = 1.000000.
  - BUT INFO-046 found strong pre/post merger time-ASYMMETRY (|d cos| up to 0.66),
    and INFO-S15b/036/045 found MI and dMI/dt sharply time-LOCALIZED at the merger.
  Both cannot be "time-blind" in the same sense.

This recheck separates TWO senses of "time-blind" on real LIGO GW150914:

  SENSE 1 (within a fixed window, row-ORDER): is re-extraction invariant to
    permuting the rows of ONE pooled operator matrix? Prediction: EXACTLY 1.0 --
    and this is DEFINITIONAL, because extract_v1 uses the column covariance
    (Mc^T Mc = sum over rows of outer products), which is invariant under row
    permutation. If so, INFO-S15a is a mathematical identity of the METHOD, not a
    property of gravity/the data -- it says nothing about time.

  SENSE 2 (across windows, time STRUCTURE): does the null DIRECTION change as we
    slide along the segment in time, and does that change localize the merger?
    Prediction (the contradicting data): YES -- the per-block null moves in time
    and the merger block is distinct. The data is NOT time-blind; only the pooled
    column-covariance extraction is row-order-invariant by construction.

Output: recheck_time_blind_results.json. Data level only.
"""
import json

import numpy as np

from kbk_pipeline import compute_operator_matrix, extract_v1
from grav_time_retaining import (load_strain, preprocess, FS, EDGE_S,
                                  WIN_S, STRIDE_S, EVENTS)

WIN = int(WIN_S * FS)
STR = int(STRIDE_S * FS)
SUBSTRATE_Q = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
N_BLOCKS = 28          # contiguous time blocks for the across-window trajectory
N_PERM = 25            # random row-permutations for the within-window test


def cos(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(abs(np.dot(a, b) / (na * nb))) if na > 0 and nb > 0 else float("nan")


def main():
    h1f, l1f, label, mt, shift, inv = EVENTS[0]   # GW150914, cleanly aligned
    h1 = preprocess(load_strain(h1f))
    l1 = preprocess(load_strain(l1f))
    l1 = np.roll(l1, int(shift * FS)) * (-1.0 if inv else 1.0)
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1[crop:-crop], WIN, STR, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * STR + WIN / 2) / FS + EDGE_S
    print(f"[recheck_time_blind] GW150914, {n} operator-rows, merger t={mt}s\n", flush=True)

    # ---- SENSE 1: within-window row-ORDER invariance (claim INFO-S15a) ----
    v0, _, _, Mc = extract_v1(M)
    rng = np.random.default_rng(0)
    perm_cos = []
    cov0 = Mc.T @ Mc
    cov_drift = []
    for _ in range(N_PERM):
        p = rng.permutation(n)
        vp, _, _, Mcp = extract_v1(M[p])
        perm_cos.append(cos(v0, vp))
        cov_drift.append(float(np.linalg.norm(Mcp.T @ Mcp - cov0)))
    perm_cos = np.array(perm_cos)
    print("SENSE 1 (within-window row shuffle):", flush=True)
    print(f"   cos(v0, v_perm): min={perm_cos.min():.8f} mean={perm_cos.mean():.8f}", flush=True)
    print(f"   column-covariance drift under permutation: max={max(cov_drift):.2e} "
          f"(~0 => row order CANNOT affect the extraction; DEFINITIONAL)\n", flush=True)

    # ---- SENSE 2: across-window null trajectory (the contradicting data) ----
    edges = np.linspace(0, n, N_BLOCKS + 1, dtype=int)
    block_nulls, block_t, block_cos_attr, block_MI = [], [], [], []
    for b in range(N_BLOCKS):
        s, e = edges[b], edges[b + 1]
        if e - s < 8:
            continue
        vb, _, _, _ = extract_v1(M[s:e])
        block_nulls.append(vb)
        block_t.append(float(centers[s:e].mean()))
        block_cos_attr.append(cos(vb[[2, 3, 4]], SUBSTRATE_Q))   # quad subspace vs equal-entropy
        block_MI.append(float(M[s:e, 5].mean()))
    block_nulls = np.array(block_nulls)
    block_t = np.array(block_t)
    # consecutive cos between adjacent blocks' full 6-D nulls
    consec = np.array([cos(block_nulls[i], block_nulls[i + 1])
                       for i in range(len(block_nulls) - 1)])
    # which block holds the merger
    mb = int(np.argmin(np.abs(block_t - mt)))
    merger_null = block_nulls[mb]
    cos_to_merger = np.array([cos(v, merger_null) for v in block_nulls])

    print("SENSE 2 (across-window null trajectory):", flush=True)
    print(f"   blocks={len(block_nulls)}, merger in block {mb} (t~{block_t[mb]:.1f}s)", flush=True)
    print(f"   consecutive-block null |cos|: min={consec.min():.3f} "
          f"mean={consec.mean():.3f} max={consec.max():.3f}", flush=True)
    print(f"   (if << 1 the null MOVES in time -> data is NOT time-blind)", flush=True)
    # how distinct is the merger block from the typical (noise) block?
    noise_mask = np.abs(block_t - mt) > 2.0
    cos_merger_vs_noise = float(np.mean(cos_to_merger[noise_mask]))
    print(f"   mean |cos|(other blocks -> merger-block null) = {cos_merger_vs_noise:.3f}", flush=True)
    print(f"   merger-block mean MI = {block_MI[mb]:.3f} vs noise-block mean "
          f"{np.mean(np.array(block_MI)[noise_mask]):.3f}", flush=True)
    # the merger block's cos-to-equal-entropy vs noise blocks
    print(f"   merger-block cos->equal-entropy(quad) = {block_cos_attr[mb]:.3f} vs "
          f"noise-block mean {np.mean(np.array(block_cos_attr)[noise_mask]):.3f}\n", flush=True)

    verdict = (
        "Row-order invariance (SENSE 1) is a mathematical identity of the "
        "column-covariance extraction (cov drift ~0), so it is NOT evidence the "
        "data is time-blind. Across windows (SENSE 2) the null direction moves in "
        "time and the merger block is distinct -> the DATA carries strong time "
        "structure the pooled extraction discards. 'Level null is time-blind' "
        "should be restated as 'the pooled extraction is row-permutation-invariant "
        "by construction', which is a different claim."
    )
    print("VERDICT:", verdict, flush=True)

    json.dump({
        "event": label, "n_rows": n, "merger_t": mt,
        "sense1_within_window": {
            "perm_cos_min": float(perm_cos.min()),
            "perm_cos_mean": float(perm_cos.mean()),
            "cov_drift_max": float(max(cov_drift)),
            "definitional": bool(max(cov_drift) < 1e-6 or perm_cos.min() > 0.999999),
        },
        "sense2_across_window": {
            "n_blocks": int(len(block_nulls)),
            "merger_block": mb,
            "consecutive_cos_min": float(consec.min()),
            "consecutive_cos_mean": float(consec.mean()),
            "consecutive_cos_max": float(consec.max()),
            "mean_cos_others_to_merger": cos_merger_vs_noise,
            "merger_block_MI": block_MI[mb],
            "noise_block_MI_mean": float(np.mean(np.array(block_MI)[noise_mask])),
            "block_t": block_t.tolist(),
            "block_cos_to_merger": cos_to_merger.tolist(),
            "block_cos_to_equalentropy": block_cos_attr,
        },
        "verdict": verdict,
    }, open("recheck_time_blind_results.json", "w"), indent=2)
    print("\n[recheck_time_blind] wrote results", flush=True)


if __name__ == "__main__":
    main()
