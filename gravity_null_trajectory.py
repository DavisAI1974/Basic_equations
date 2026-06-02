"""
Gravity null-DIRECTION trajectory through inspiral -> merger -> ringdown.

Motivated by the time-blind recheck (recheck_time_blind.py): the POOLED null is
row-permutation-invariant by construction, but the TIME-RESOLVED null moves
strongly in time and localizes the merger. So we can watch gravity's time
structure WITHOUT a rate-dynamics fit (which the merger window can't support):
just extract the null in a sliding local block and track its trajectory.

Tests the SHAPE (data level only; H-D "return to bare substrate" is the
operator's frame, NOT graded here): across a merger, does the information
structure go substrate -> departure -> substrate?

Per sliding block (K operator-rows) we record, as a function of time:
  - cos6 : |cos(v_null, equal-entropy attractor [0,0,-1,-1,2,0]/sqrt6)| in the
           FULL 6-op basis. INFO-036: noise ~0.98 (on), merger ~0.23 (off,
           because MI enters the null at the merger).
  - cosQ : |cos(v_null[H^2 subspace], (-1,-1,2)/sqrt6)| -- the no-MI/quad view.
  - |MI coef| in v_null (index 5): does MI ENTER the null at the merger?
  - mean MI and mean |H_a - H_b| in the block.
  - block-to-block null rotation (1 - consecutive cos): structural transition rate.

Honest gate (INFO-046): at 125 ms / 35-350 Hz the POST-merger window is
dominated by non-stationary detector noise, so any "return" post-merger may be
return-to-noise, not ringdown physics. Reported, not interpreted.

Output: gravity_null_trajectory_results.json
"""
import json

import numpy as np

from kbk_pipeline import compute_operator_matrix, extract_v1
from grav_time_retaining import load_strain, preprocess, FS, EDGE_S, WIN_S, STRIDE_S
from grav_flow_crossevent import align

WIN = int(WIN_S * FS)
STR = int(STRIDE_S * FS)
ATTR6 = np.array([0, 0, -1.0, -1.0, 2.0, 0]) / np.sqrt(6.0)
ATTRQ = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
K = 16          # operator-rows per local block (~0.6 s span)
STEP = 2        # step in rows (~62.5 ms)

# (H1, L1, label, known merger time into 32 s segment, clean?)
EVENTS = [
    ("data/ligo/H-H1_GW150914_32s.hdf5", "data/ligo/L-L1_GW150914_32s.hdf5", "GW150914", 16.4, True),
    ("data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5", "GW170104", 16.6, False),
    ("data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5", "GW151226", 16.6, False),
]


def cosabs(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(abs(np.dot(a, b) / (na * nb))) if na > 0 and nb > 0 else float("nan")


def run(h1f, l1f, label, mt, clean):
    h1 = preprocess(load_strain(h1f))
    l1 = preprocess(load_strain(l1f))
    if clean:
        l1 = np.roll(l1, int(0.0069 * FS)) * -1.0      # GW150914 known shift
        corr = float("nan")
    else:
        l1, lag, sign, corr = align(h1, l1, mt)
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1[crop:-crop], WIN, STR, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * STR + WIN / 2) / FS + EDGE_S

    t, cos6, cosQ, mico, mi, asym, vlist = [], [], [], [], [], [], []
    for s in range(0, n - K + 1, STEP):
        blk = M[s:s + K]
        v, _, _, _ = extract_v1(blk)
        t.append(float(centers[s:s + K].mean()))
        cos6.append(cosabs(v, ATTR6))
        cosQ.append(cosabs(v[[2, 3, 4]], ATTRQ))
        mico.append(float(abs(v[5])))
        mi.append(float(blk[:, 5].mean()))
        asym.append(float(np.abs(blk[:, 0] - blk[:, 1]).mean()))
        vlist.append(v)
    t = np.array(t)
    rot = [1.0 - cosabs(vlist[i], vlist[i + 1]) for i in range(len(vlist) - 1)]

    # merger-localized vs away ("noise") summary
    near = np.abs(t - mt) <= 0.25
    far = np.abs(t - mt) > 2.0
    mi_idx = int(np.argmin(np.abs(t - mt)))

    def mm(a, mask):
        a = np.array(a)
        return [round(float(a[mask].mean()), 3), round(float(a[mask].std()), 3)]

    summary = {
        "align_corr": corr,
        "merger_idx": mi_idx, "merger_t": float(t[mi_idx]),
        "cos6_near_meanstd": mm(cos6, near), "cos6_far_meanstd": mm(cos6, far),
        "cosQ_near_meanstd": mm(cosQ, near), "cosQ_far_meanstd": mm(cosQ, far),
        "MIcoef_near_meanstd": mm(mico, near), "MIcoef_far_meanstd": mm(mico, far),
        "MI_near_meanstd": mm(mi, near), "MI_far_meanstd": mm(mi, far),
        "max_rotation": round(float(np.max(rot)), 3),
        "rotation_peak_t": round(float(t[1:][int(np.argmax(rot))]), 2),
    }
    series = {"t": t.tolist(), "cos6": cos6, "cosQ": cosQ, "MIcoef": mico,
              "MI": mi, "asym": asym}
    return {"label": label, "merger_t_known": mt, "clean": clean,
            "summary": summary, "series": series}


def main():
    print(f"[gravity_null_trajectory] K={K} rows (~0.6s), step={STEP} (~62ms)\n", flush=True)
    out = []
    for e in EVENTS:
        r = run(*e)
        out.append(r)
        s = r["summary"]
        print(f"=== {r['label']} (corr={s['align_corr'] if isinstance(s['align_corr'],str) else round(s['align_corr'],2) if s['align_corr']==s['align_corr'] else 'clean'}) merger t~{s['merger_t']:.1f}s ===", flush=True)
        print(f"   cos6  near={s['cos6_near_meanstd']}  far={s['cos6_far_meanstd']}  "
              f"(INFO-036: near should DIP if MI enters null at merger)", flush=True)
        print(f"   cosQ  near={s['cosQ_near_meanstd']}  far={s['cosQ_far_meanstd']}  "
              f"(quad/no-MI view; near should be HIGH)", flush=True)
        print(f"   |MIcoef| near={s['MIcoef_near_meanstd']}  far={s['MIcoef_far_meanstd']}  "
              f"(does MI ENTER the null at merger?)", flush=True)
        print(f"   MI    near={s['MI_near_meanstd']}  far={s['MI_far_meanstd']}", flush=True)
        print(f"   null rotation peaks {s['max_rotation']} at t={s['rotation_peak_t']}s "
              f"(merger at {r['merger_t_known']}s)\n", flush=True)
    json.dump({"K": K, "step": STEP, "events": out}, open("gravity_null_trajectory_results.json", "w"), indent=2)
    print("[gravity_null_trajectory] wrote results", flush=True)


if __name__ == "__main__":
    main()
