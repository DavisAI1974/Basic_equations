"""
Eigen quasispecies probe with explicit selection + mutation.

Tests reading (c) from the multi-seed result: does a system with
explicit selection + mutation (not just tunable coupling) move the
operator-space signature off the universal Session 3 attractor?

Speaking posture before (Rule C): I think the quasispecies system
might or might not move off attractor. If selection + mutation is
the missing ingredient for a chemistry-to-biology operator-space
transition, we'll see it. If it's not, we won't. No verdict in
advance.

System: two sequence types A and B (e.g., master vs mutant class),
each replicating with fidelity q (per-site copy fidelity). At each
step:

  Replication:
    A produces A with rate (1-mu) and B with rate mu
    B produces B with rate (1-mu) and A with rate mu
    Both at base rate r times fitness f_A or f_B.
  Selection:
    Total population kept constant by outflow phi.

  dA/dt = ((1-mu)*f_A - phi)*A  +  mu*f_B*B
  dB/dt = ((1-mu)*f_B - phi)*B  +  mu*f_A*A
  phi   = (f_A*A + f_B*B) / (A + B)

Mu (mutation rate) is the swept parameter:
  mu small             -> selection wins, A or B dominates
  mu = mu_critical     -> Eigen error threshold, master loses lock
  mu large             -> mutational drift, both at equilibrium

Three regimes correspond to three different dynamical phases of
information storage. Operating Rule: no pre-assigned meaning to
where each regime lands in operator space.

Per-domain stack on the ensemble: H_a(t), H_b(t), MI(t), 6-op
operator matrix M, KBK null direction v_null6, cos to Session 3
attractor, fits for INFO-025 families.

Sweep: mu in [0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5].
       f_A = 1.5 (selected master), f_B = 1.0 (less fit mutant).
       Three seeds (11, 22, 33) at Session 7 baseline.

Output: quasispecies_probe.json
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


def simulate_quasispecies(
    seed: int,
    mu: float,
    f_A: float = 1.5,
    f_B: float = 1.0,
    N_ens: int = 600,
    T: float = 30.0,
    dt: float = 0.02,
    noise: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """Two-class quasispecies with mutation mu and fitness ratio f_A/f_B.

    Returns (X_A, X_B) shape (steps, N_ens) for ensemble across seeds.
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    # Initial: half-half with small noise to break symmetry per realization
    A = 0.5 + 0.05 * rng.standard_normal(N_ens)
    B = 0.5 + 0.05 * rng.standard_normal(N_ens)
    A = np.clip(A, 0.05, None)
    B = np.clip(B, 0.05, None)
    XA = np.zeros((steps, N_ens))
    XB = np.zeros((steps, N_ens))
    for t in range(steps):
        tot = A + B
        phi = (f_A * A + f_B * B) / tot
        dA = ((1.0 - mu) * f_A * A + mu * f_B * B - phi * A) * dt
        dB = ((1.0 - mu) * f_B * B + mu * f_A * A - phi * B) * dt
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


def run_one(mu: float, seed: int,
            N_ens: int = 600, T: float = 30.0, dt: float = 0.02) -> dict:
    t0 = time.time()
    XA, XB = simulate_quasispecies(seed, mu, N_ens=N_ens, T=T, dt=dt)
    M, _ = build_ensemble_operator_matrix(XA, XB)
    v_null, S, _, _ = extract_v1(M)
    eigs = (S ** 2) / max(1, len(M) - 1)
    cos_attr, sub234 = cos_to_attractor(v_null)
    fits = fit_poly_and_exp(M)
    return {
        "mu": mu, "seed": seed,
        "v_null6": v_null.tolist(),
        "v_null_234_normalized": sub234,
        "cos_to_session3_attractor": cos_attr,
        "MI_mean": float(M[:, 5].mean()),
        "MI_std": float(M[:, 5].std()),
        "Ha_mean": float(M[:, 0].mean()),
        "Hb_mean": float(M[:, 1].mean()),
        "Ha_minus_Hb_mean": float((M[:, 0] - M[:, 1]).mean()),
        "Ha_minus_Hb_std": float((M[:, 0] - M[:, 1]).std()),
        "eigs_bottom3": eigs[-3:].tolist(),
        "fits": fits,
        "elapsed_s": float(time.time() - t0),
    }


def main():
    t_start = time.time()
    seeds = [11, 22, 33]
    mu_grid = [0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5]
    out = {
        "speaking_posture_before": (
            "I think quasispecies may or may not move off attractor "
            "with mutation. The three regimes (selection dominant / "
            "near error threshold / drift dominant) may show in MI "
            "magnitude, in v6 direction, in Ha-Hb asymmetry, or in "
            "none of them. We'll wait on the data."
        ),
        "baseline": {"N_ens": 600, "T": 30.0, "dt": 0.02,
                     "f_A": 1.5, "f_B": 1.0},
        "seeds": seeds,
        "mu_grid": mu_grid,
        "session3_attractor_234": SESSION3_ATTRACTOR_234.tolist(),
        "per_mu": {},
    }
    for mu in mu_grid:
        per_seed = []
        for s in seeds:
            print(f"[qs] mu={mu} seed={s}", flush=True)
            rec = run_one(mu, s)
            per_seed.append(rec)
            print(
                f"  done {rec['elapsed_s']:.1f}s  "
                f"MI={rec['MI_mean']:.3f}  "
                f"|Ha-Hb|={abs(rec['Ha_minus_Hb_mean']):.3f}  "
                f"cos_attr={rec['cos_to_session3_attractor']:.3f}  "
                f"poly1_r2={rec['fits']['poly1']['r2']:.3f}  "
                f"exp_r2={rec['fits']['exp_INFO025']['r2']:.3f}",
                flush=True,
            )
        v_list = [r["v_null6"] for r in per_seed]
        agg = inter_seed_cosines(v_list)
        cas = [r["cos_to_session3_attractor"] for r in per_seed]
        mis = [r["MI_mean"] for r in per_seed]
        hms = [r["Ha_minus_Hb_mean"] for r in per_seed]
        out["per_mu"][f"{mu:.3f}"] = {
            "per_seed": per_seed,
            "cross_seed_v_null6": agg,
            "cos_to_attractor_mean": float(np.mean(cas)),
            "cos_to_attractor_std":  float(np.std(cas)),
            "MI_mean_per_seed": mis,
            "Ha_minus_Hb_mean_per_seed": hms,
        }
    out["total_elapsed_s"] = float(time.time() - t_start)

    print("\n=== summary ===")
    print(f"{'mu':>7} {'MI_mean':>9} {'<|Ha-Hb|>':>11} "
          f"{'cos_attr':>9} {'cos_std':>9} {'v6_min':>9} {'v6_mean':>9}")
    for mu in mu_grid:
        key = f"{mu:.3f}"
        d = out["per_mu"][key]
        mi = float(np.mean(d["MI_mean_per_seed"]))
        hm = float(np.mean([abs(x) for x in d["Ha_minus_Hb_mean_per_seed"]]))
        ca = d["cos_to_attractor_mean"]
        cs = d["cos_to_attractor_std"]
        vm = d["cross_seed_v_null6"]["min_pair_abs_cos"]
        va = d["cross_seed_v_null6"]["mean_pair_abs_cos"]
        print(f"{mu:>7.3f} {mi:>9.3f} {hm:>11.3f} {ca:>9.3f} {cs:>9.3f} "
              f"{vm:>9.3f} {va:>9.3f}")
    print(f"\nTotal: {out['total_elapsed_s']:.1f}s")

    with open("quasispecies_probe.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Wrote quasispecies_probe.json")


if __name__ == "__main__":
    main()
