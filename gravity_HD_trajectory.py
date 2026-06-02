"""
H-D trajectory, PINNED: does the information structure go substrate -> departure
-> substrate across the merger? (Operator hunch H-D, recorded NOT graded.)

Fix for gravity_null_trajectory's resolution wall: there the per-block null at
K=16 rows was estimator-noise-dominated away from the merger (INFO-024), so the
"substrate baseline" was unstable. Here:
  - a STABLE pooled-noise reference null v_noise (all windows >2s from merger),
  - a MODERATE sliding window (K=40 rows ~1.4s -> well-estimated null) stepped
    finely (~0.25s), so each window's null is stable AND the ~0.2s merger
    transient still lands inside a window.
Per window track |cos(v, v_noise)| and |cos(v, equal-entropy attractor)|.

H-D shape (data level, no grading): noise windows ~ on substrate; merger window
DEPARTS; post-merger recovers? Honest gate (INFO-046): post-merger at 125ms/
35-350Hz is detector noise, so a "recovery" there is return-to-noise, not
ringdown physics -- reported, not interpreted.

Output: gravity_HD_trajectory_results.json
"""
import json

import numpy as np

from kbk_pipeline import compute_operator_matrix, extract_v1
from grav_time_retaining import load_strain, preprocess, FS, EDGE_S, WIN_S, STRIDE_S

WIN = int(WIN_S * FS)
STR = int(STRIDE_S * FS)
ATTR6 = np.array([0, 0, -1.0, -1.0, 2.0, 0]) / np.sqrt(6.0)
ATTRQ = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
K = 40          # rows per sliding window (~1.4s; stable null)
STEP = 8        # ~0.25s
MERGER_T = 16.4


def cosabs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(abs(np.dot(a, b) / (na * nb))) if na > 0 and nb > 0 else float("nan")


def main():
    h1 = preprocess(load_strain("data/ligo/H-H1_GW150914_32s.hdf5"))
    l1 = preprocess(load_strain("data/ligo/L-L1_GW150914_32s.hdf5"))
    l1 = np.roll(l1, int(0.0069 * FS)) * -1.0
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1[crop:-crop], WIN, STR, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * STR + WIN / 2) / FS + EDGE_S

    # stable POOLED-NOISE reference null (all rows >2s from merger)
    noise_rows = np.abs(centers - MERGER_T) > 2.0
    v_noise, _, _, _ = extract_v1(M[noise_rows])
    print(f"[gravity_HD_trajectory] pooled-noise reference: "
          f"|cos(v_noise, equal-entropy)| = {cosabs(v_noise, ATTR6):.3f} "
          f"(INFO-036: pooled noise sits ON the attractor)\n", flush=True)

    t, cos_noise, cos6, cosQ, mi = [], [], [], [], []
    for s in range(0, n - K + 1, STEP):
        v, _, _, _ = extract_v1(M[s:s + K])
        t.append(float(centers[s:s + K].mean()))
        cos_noise.append(cosabs(v, v_noise))
        cos6.append(cosabs(v, ATTR6))
        cosQ.append(cosabs(v[[2, 3, 4]], ATTRQ))
        mi.append(float(M[s:s + K, 5].mean()))
    t = np.array(t)

    mb = int(np.argmin(np.abs(t - MERGER_T)))
    pre = (t < MERGER_T - 0.6) & (t > MERGER_T - 2.5)
    post = (t > MERGER_T + 0.6) & (t < MERGER_T + 2.5)

    def mm(a, m):
        a = np.array(a)
        return [round(float(a[m].mean()), 3), round(float(a[m].std()), 3)]

    print(f"merger window idx {mb} (t~{t[mb]:.2f}s)", flush=True)
    print(f"  |cos to pooled-noise|:  pre={mm(cos_noise,pre)}  "
          f"MERGER={round(cos_noise[mb],3)}  post={mm(cos_noise,post)}", flush=True)
    print(f"  |cos to equal-entropy|: pre={mm(cos6,pre)}  "
          f"MERGER={round(cos6[mb],3)}  post={mm(cos6,post)}", flush=True)
    print(f"  MI:                     pre={mm(mi,pre)}  "
          f"MERGER={round(mi[mb],3)}  post={mm(mi,post)}", flush=True)

    pre_n, post_n = mm(cos_noise, pre)[0], mm(cos_noise, post)[0]
    dip = pre_n - cos_noise[mb]
    recov = post_n - cos_noise[mb]
    shape = ("substrate->DEPART->return" if (dip > 0.15 and recov > 0.10)
             else "DEPART-no-clear-return" if dip > 0.15
             else "no clear departure")
    print(f"\n  departure at merger (pre - merger) = {dip:+.3f}; "
          f"post-recovery (post - merger) = {recov:+.3f}", flush=True)
    print(f"  SHAPE (data level): {shape}", flush=True)
    print("  [gate: post-merger at 125ms/35-350Hz is detector noise (INFO-046); "
          "a 'return' there is return-to-noise, not ringdown]", flush=True)

    json.dump({
        "event": "GW150914", "K_rows": K, "step": STEP, "merger_t": MERGER_T,
        "v_noise_cos_to_equal_entropy": cosabs(v_noise, ATTR6),
        "merger_idx": mb,
        "cos_to_noise": {"pre": mm(cos_noise, pre), "merger": cos_noise[mb],
                         "post": mm(cos_noise, post)},
        "cos_to_equal_entropy": {"pre": mm(cos6, pre), "merger": cos6[mb],
                                 "post": mm(cos6, post)},
        "MI": {"pre": mm(mi, pre), "merger": mi[mb], "post": mm(mi, post)},
        "departure_dip": dip, "post_recovery": recov, "shape": shape,
        "series": {"t": t.tolist(), "cos_to_noise": cos_noise,
                   "cos_to_equal_entropy": cos6, "cosQ": cosQ, "MI": mi},
        "gate": "post-merger window is detector noise at this band/resolution (INFO-046)",
    }, open("gravity_HD_trajectory_results.json", "w"), indent=2)
    print("\n[gravity_HD_trajectory] wrote results", flush=True)


if __name__ == "__main__":
    main()
