"""
Control test: asymmetric hypercycle (different base rates per species,
NO selection / fitness difference, NO mutation).

Disentangles two candidate readings of the quasispecies off-attractor
finding:
  (a/b) selection or information-storage moves systems off attractor
  (c)   channel asymmetry alone moves systems off attractor

If channel asymmetry alone (driven by base_A != base_B with NO
selection mechanism) puts this system off attractor, reading (c) wins.
If it stays on attractor despite asymmetry, reading (a/b) survives.

Speaking posture before (Rule C): I think the asymmetric hypercycle
might stay on attractor at all base-ratio values because the bilinear
mutual catalysis term keeps the channels coupled symmetrically through
each other. Or it might go off attractor because base-rate asymmetry
forces |H_a - H_b| > 0 regardless. I won't pre-commit. We'll wait on
the data.

System: same Eigen hypercycle as prebiotic_canary.py but with
base_A != base_B and fixed k_cat=1.0 (so coupling is active).

  dA/dt = A * (f_A - phi)        f_A = base_A + k_cat * B
  dB/dt = B * (f_B - phi)        f_B = base_B + k_cat * A
  phi    = (A*f_A + B*f_B) / (A + B)

Sweep: base_A in [1.0, 1.1, 1.25, 1.5, 2.0, 3.0] with base_B fixed
at 1.0. Three seeds. Session 7 baseline.

If channel asymmetry alone produces the off-attractor pattern, we
expect cos_to_attractor to drop monotonically from base_ratio=1.0
toward 0 as base_ratio increases.

Output: asymmetric_hypercycle_control.json
"""
from __future__ import annotations
import json
import time
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix
from kbk_pipeline import extract_v1
from prebiotic_canary import (
    fit_poly_and_exp, cos_to_attractor, SESSION3_ATTRACTOR_234,
)


def simulate_asymmetric_hypercycle(
    seed: int,
    base_A: float,
    base_B: float = 1.0,
    k_cat: float = 1.0,
    N_ens: int = 600,
    T: float = 30.0,
    dt: float = 0.02,
    noise: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    A = 0.5 + 0.05 * rng.standard_normal(N_ens)
    B = 0.5 + 0.05 * rng.standard_normal(N_ens)
    A = np.clip(A, 0.05, None)
    B = np.clip(B, 0.05, None)
    XA = np.zeros((steps, N_ens))
    XB = np.zeros((steps, N_ens))
    for t in range(steps):
        fA = base_A + k_cat * B
        fB = base_B + k_cat * A
        tot = A + B
        phi = (A * fA + B * fB) / tot
        dA = A * (fA - phi) * dt
        dB = B * (fB - phi) * dt
        A = A + dA + noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        B = B + dB + noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        A = np.clip(A, 0.01, None)
        B = np.clip(B, 0.01, None)
        XA[t] = A
        XB[t] = B
    return XA, XB


def inter_seed_cosines(v_list: list[list[float]]) -> dict:
    V = np.array(v_list)
    n = len(V)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            ni = np.linalg.norm(V[i])
            nj = np.linalg.norm(V[j])
            if ni > 0 and nj > 0:
                pairs.append(abs(float(np.dot(V[i], V[j]) / (ni * nj))))
    return {
        "min_pair_abs_cos": float(np.min(pairs)) if pairs else float("nan"),
        "mean_pair_abs_cos": float(np.mean(pairs)) if pairs else float("nan"),
    }


def run_one(base_A: float, seed: int) -> dict:
    t0 = time.time()
    XA, XB = simulate_asymmetric_hypercycle(seed, base_A)
    M, _ = build_ensemble_operator_matrix(XA, XB)
    v_null, S, _, _ = extract_v1(M)
    cos_attr, sub234 = cos_to_attractor(v_null)
    fits = fit_poly_and_exp(M)
    return {
        "base_A": base_A, "seed": seed,
        "v_null6": v_null.tolist(),
        "v_null_234_normalized": sub234,
        "cos_to_session3_attractor": cos_attr,
        "MI_mean": float(M[:, 5].mean()),
        "MI_std": float(M[:, 5].std()),
        "Ha_mean": float(M[:, 0].mean()),
        "Hb_mean": float(M[:, 1].mean()),
        "Ha_minus_Hb_mean": float((M[:, 0] - M[:, 1]).mean()),
        "fits": fits,
        "elapsed_s": float(time.time() - t0),
    }


def main():
    t_start = time.time()
    seeds = [11, 22, 33]
    base_A_grid = [1.0, 1.1, 1.25, 1.5, 2.0, 3.0]
    out = {
        "speaking_posture_before": (
            "Two roads. If channel asymmetry alone drives off-attractor "
            "behavior, this control will mirror the quasispecies pattern. "
            "If selection per se is required, this control stays on "
            "attractor at all base_A values. Wait on the data."
        ),
        "baseline": {"N_ens": 600, "T": 30.0, "dt": 0.02,
                     "base_B": 1.0, "k_cat": 1.0},
        "seeds": seeds,
        "base_A_grid": base_A_grid,
        "session3_attractor_234": SESSION3_ATTRACTOR_234.tolist(),
        "per_base_A": {},
    }
    for ba in base_A_grid:
        per_seed = []
        for s in seeds:
            print(f"[asymctrl] base_A={ba} seed={s}", flush=True)
            rec = run_one(ba, s)
            per_seed.append(rec)
            print(
                f"  done {rec['elapsed_s']:.1f}s  "
                f"MI={rec['MI_mean']:.3f}  "
                f"<Ha-Hb>={rec['Ha_minus_Hb_mean']:+.3f}  "
                f"cos_attr={rec['cos_to_session3_attractor']:.3f}",
                flush=True,
            )
        v_list = [r["v_null6"] for r in per_seed]
        agg = inter_seed_cosines(v_list)
        cas = [r["cos_to_session3_attractor"] for r in per_seed]
        mis = [r["MI_mean"] for r in per_seed]
        hms = [r["Ha_minus_Hb_mean"] for r in per_seed]
        out["per_base_A"][f"{ba:.2f}"] = {
            "per_seed": per_seed,
            "cross_seed_v_null6": agg,
            "cos_to_attractor_mean": float(np.mean(cas)),
            "cos_to_attractor_std":  float(np.std(cas)),
            "MI_mean_per_seed": mis,
            "Ha_minus_Hb_mean_per_seed": hms,
        }
    out["total_elapsed_s"] = float(time.time() - t_start)

    print("\n=== summary (asymmetric hypercycle, no selection) ===")
    print(f"{'base_A':>7} {'MI_mean':>9} {'<Ha-Hb>':>9} "
          f"{'cos_attr':>9} {'cos_std':>9} {'v6_min':>9}")
    for ba in base_A_grid:
        key = f"{ba:.2f}"
        d = out["per_base_A"][key]
        mi = float(np.mean(d["MI_mean_per_seed"]))
        hm = float(np.mean(d["Ha_minus_Hb_mean_per_seed"]))
        ca = d["cos_to_attractor_mean"]
        cs = d["cos_to_attractor_std"]
        vm = d["cross_seed_v_null6"]["min_pair_abs_cos"]
        print(f"{ba:>7.2f} {mi:>9.3f} {hm:>+9.3f} "
              f"{ca:>9.3f} {cs:>9.3f} {vm:>9.3f}")
    print(f"\nTotal: {out['total_elapsed_s']:.1f}s")

    with open("asymmetric_hypercycle_control.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Wrote asymmetric_hypercycle_control.json")


if __name__ == "__main__":
    main()
