"""
Prebiotic multi-seed replication at Session 7 baseline.

Question carried over from canary: is the k_cat=0.25 off-attractor
point (cos 0.405 at one seed) a real transition-zone signature, a
sample-noise excursion, or something else?

Speaking posture before (Rule C): I think the k_cat=0.25 deviation
might survive multi-seed replication if it's real, or wash out if
it was noise; either way we'll wait on the data.

Design: same Eigen hypercycle as prebiotic_canary.py. Now at
Session 7 baseline (N_ens=600, T=30, dt=0.02), three seeds (11, 22,
33 -- the first two are Session 7's seeds), same k_cat grid.

For each k_cat we report: per-seed v_null6 and cos-to-Session-3
attractor; cross-seed v_null6 cosine matrix (inter-seed consistency);
mean/std MI and functional-family fits.

Output: prebiotic_multiseed.json
"""
from __future__ import annotations
import json
import time
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix
from kbk_pipeline import extract_v1
from prebiotic_canary import (
    simulate_hypercycle, fit_poly_and_exp, cos_to_attractor,
    SESSION3_ATTRACTOR_234,
)


def run_one(k_cat: float, seed: int,
            N_ens: int = 600, T: float = 30.0, dt: float = 0.02) -> dict:
    t0 = time.time()
    XA, XB = simulate_hypercycle(seed, k_cat, N_ens=N_ens, T=T, dt=dt)
    M, _ = build_ensemble_operator_matrix(XA, XB)
    v_null, S, _, _ = extract_v1(M)
    eigs = (S ** 2) / max(1, len(M) - 1)
    cos_attr, sub234 = cos_to_attractor(v_null)
    fits = fit_poly_and_exp(M)
    return {
        "k_cat": k_cat, "seed": seed,
        "v_null6": v_null.tolist(),
        "v_null_234_normalized": sub234,
        "cos_to_session3_attractor": cos_attr,
        "MI_mean": float(M[:, 5].mean()),
        "MI_std": float(M[:, 5].std()),
        "Ha_mean": float(M[:, 0].mean()),
        "Hb_mean": float(M[:, 1].mean()),
        "eigs_top5": eigs[:5].tolist(),
        "eigs_bottom3": eigs[-3:].tolist(),
        "fits": fits,
        "elapsed_s": float(time.time() - t0),
    }


def inter_seed_cosines(v_list: list[list[float]]) -> dict:
    V = np.array(v_list)
    cm = np.zeros((len(V), len(V)))
    for i in range(len(V)):
        for j in range(len(V)):
            ni = np.linalg.norm(V[i])
            nj = np.linalg.norm(V[j])
            if ni > 0 and nj > 0:
                cm[i, j] = abs(float(np.dot(V[i], V[j]) / (ni * nj)))
            else:
                cm[i, j] = float("nan")
    pairs = [cm[0, 1], cm[0, 2], cm[1, 2]]
    return {
        "cos_matrix_abs": cm.tolist(),
        "min_pair_abs_cos": float(np.min(pairs)),
        "mean_pair_abs_cos": float(np.mean(pairs)),
    }


def main():
    t_start = time.time()
    seeds = [11, 22, 33]
    k_grid = [0.0, 0.25, 0.5, 1.0, 2.0]
    out = {
        "speaking_posture_before": (
            "I think the k_cat=0.25 off-attractor point from the canary "
            "might survive replication or wash out at multi-seed. Either "
            "is informative. No verdict in advance."
        ),
        "baseline": {"N_ens": 600, "T": 30.0, "dt": 0.02},
        "seeds": seeds,
        "k_cat_grid": k_grid,
        "session3_attractor_234": SESSION3_ATTRACTOR_234.tolist(),
        "per_k_cat": {},
    }
    for k in k_grid:
        per_seed = []
        for s in seeds:
            print(f"[prebiotic-MS] k_cat={k} seed={s}", flush=True)
            rec = run_one(k, s)
            per_seed.append(rec)
            print(
                f"  done {rec['elapsed_s']:.1f}s  "
                f"MI={rec['MI_mean']:.3f}  "
                f"cos_attr={rec['cos_to_session3_attractor']:.3f}  "
                f"poly1_r2={rec['fits']['poly1']['r2']:.3f}  "
                f"exp_r2={rec['fits']['exp_INFO025']['r2']:.3f}",
                flush=True,
            )
        # Cross-seed analysis at this k_cat
        v_list_6 = [r["v_null6"] for r in per_seed]
        agg_6 = inter_seed_cosines(v_list_6)
        cos_attrs = [r["cos_to_session3_attractor"] for r in per_seed]
        mis = [r["MI_mean"] for r in per_seed]
        out["per_k_cat"][f"{k:.2f}"] = {
            "per_seed": per_seed,
            "cross_seed_v_null6_cosines": agg_6,
            "cos_to_attractor_mean": float(np.mean(cos_attrs)),
            "cos_to_attractor_std":  float(np.std(cos_attrs)),
            "cos_to_attractor_per_seed": cos_attrs,
            "MI_mean_per_seed": mis,
        }
    out["total_elapsed_s"] = float(time.time() - t_start)

    # Summary lines
    print("\n=== summary ===")
    print(f"{'k_cat':>6} {'MI_mean':>9} {'cos_attr':>9} {'cos_attr_std':>13} "
          f"{'v6_inter_min':>13} {'v6_inter_mean':>14}")
    for k in k_grid:
        key = f"{k:.2f}"
        d = out["per_k_cat"][key]
        mi = np.mean(d["MI_mean_per_seed"])
        ca = d["cos_to_attractor_mean"]
        cs = d["cos_to_attractor_std"]
        vm = d["cross_seed_v_null6_cosines"]["min_pair_abs_cos"]
        va = d["cross_seed_v_null6_cosines"]["mean_pair_abs_cos"]
        print(f"{k:>6.2f} {mi:>9.3f} {ca:>9.3f} {cs:>13.3f} "
              f"{vm:>13.3f} {va:>14.3f}")
    print(f"\nTotal: {out['total_elapsed_s']:.1f}s")

    with open("prebiotic_multiseed.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Wrote prebiotic_multiseed.json")


if __name__ == "__main__":
    main()
