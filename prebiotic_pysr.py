"""
PySR symbolic regression on prebiotic ensemble-H data.

Two systems, six parameter points each (selected from the full
multi-seed sweep), 2 seeds each. Mirrors Session 7 INFO-025
methodology so results can be placed directly against the four
INFO-025 functional families:
  physics_duffing:   MI = (H_b - H_a)^2 + 0.275
  biology_lotka:     MI = 0.50 * exp(H_a / 2.0)
  chemistry_brussel: MI = 0.73 * H_a + 1.075
  geology_bk:        MI = 0.199 (constant)

For each (system, parameter, seed) we fit PySR with the Session 7
operator set {+, -, *, /, square, cube, exp, log, sqrt} and report
the Pareto front + the sympy "best" equation.

Reading goal:
  - hypercycle: where does the prebiotic-with-coupling sit relative
    to the four INFO-025 families?
  - quasispecies: do the off-attractor mu values (low mu, selection
    regime) and the on-attractor mu values (high mu, drift regime)
    have DIFFERENT best symbolic forms?

Operating Rule C: speaking posture before AND after. No verdict
in advance.

Settings: niter=30, populations=15, maxsize=12, timeout=60s per
fit (Session 7 used 40/20/15/120 - we use lighter for speed at
multi-condition scale). 2 seeds for cross-seed reproducibility
(Session 7 standard).

Output: prebiotic_pysr.json
"""
from __future__ import annotations
import json
import time
import warnings
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix
from prebiotic_canary import simulate_hypercycle
from quasispecies_probe import simulate_quasispecies


warnings.filterwarnings("ignore", category=UserWarning)


def fit_pysr_one(X, y, niter=30, populations=15, maxsize=12, seed=0,
                 timeout_s=60):
    from pysr import PySRRegressor
    model = PySRRegressor(
        niterations=niter,
        populations=populations,
        maxsize=maxsize,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["square", "cube", "exp", "log", "sqrt"],
        elementwise_loss="loss(x, y) = (x - y)^2",
        progress=False,
        random_state=seed,
        deterministic=True,
        parallelism="serial",
        verbosity=0,
        timeout_in_seconds=timeout_s,
        warm_start=False,
    )
    t0 = time.time()
    model.fit(X, y)
    fit_t = time.time() - t0
    eq = model.equations_
    rows = []
    for _, row in eq.iterrows():
        rows.append({
            "complexity": int(row["complexity"]),
            "loss": float(row["loss"]),
            "score": float(row["score"]),
            "equation": str(row["equation"]),
        })
    best_eq = str(model.sympy())
    return {
        "fit_time_s": float(fit_t),
        "pareto_front": rows,
        "best_equation_sympy": best_eq,
    }


def get_M(sim_fn, sim_kwargs, seed, N_ens=600, T=30.0, dt=0.02):
    X1, X2 = sim_fn(seed=seed, N_ens=N_ens, T=T, dt=dt, **sim_kwargs)
    M, _ = build_ensemble_operator_matrix(X1, X2)
    return M


def fit_one_condition(label, sim_fn, sim_kwargs, seeds, max_n=600,
                      niter=30, timeout_s=60):
    """Run PySR per-seed and report cross-seed Pareto fronts."""
    per_seed = []
    for s in seeds:
        M = get_M(sim_fn, sim_kwargs, seed=s)
        Ha = M[:, 0]; Hb = M[:, 1]; MI = M[:, 5]
        n = len(MI)
        if n > max_n:
            idx = np.linspace(0, n - 1, max_n).astype(int)
            Ha = Ha[idx]; Hb = Hb[idx]; MI = MI[idx]
        X = np.column_stack([Ha, Hb])
        print(f"  [pysr] {label}  seed={s}  n={len(MI)}", flush=True)
        r = fit_pysr_one(X, MI, niter=niter, seed=s,
                         timeout_s=timeout_s)
        r["seed"] = s
        r["MI_mean"] = float(MI.mean())
        r["MI_std"] = float(MI.std())
        r["Ha_minus_Hb_mean"] = float((Ha - Hb).mean())
        print(f"    {r['fit_time_s']:.1f}s  best: {r['best_equation_sympy']}",
              flush=True)
        per_seed.append(r)
    return {"label": label, "per_seed": per_seed}


def main():
    t_start = time.time()
    seeds = [11, 22]
    out = {
        "speaking_posture_before": (
            "I think the off-attractor quasispecies mu values may "
            "show a different best symbolic form than the on-attractor "
            "ones. The hypercycle should give a single coherent form "
            "across k_cat or vary smoothly. We'll wait on the data."
        ),
        "baseline": {"N_ens": 600, "T": 30.0, "dt": 0.02,
                     "niter": 30, "populations": 15, "maxsize": 12,
                     "timeout_s": 60},
        "seeds": seeds,
        "info025_reference": {
            "physics_duffing":   "(H_b - H_a)^2 + 0.275",
            "biology_lotka":     "0.5 * exp(H_a / 2)",
            "chemistry_brussel": "0.73 * H_a + 1.075",
            "geology_bk":        "0.199",
        },
        "hypercycle": [],
        "quasispecies": [],
    }

    # Hypercycle conditions (3 representative k_cat values from
    # the multi-seed result -- pre-coupling, active, saturation)
    hc_grid = [0.0, 0.5, 2.0]
    for k in hc_grid:
        label = f"hypercycle_k_cat_{k:.2f}"
        rec = fit_one_condition(
            label,
            simulate_hypercycle,
            sim_kwargs={"k_cat": k},
            seeds=seeds,
        )
        rec["k_cat"] = k
        out["hypercycle"].append(rec)

    # Quasispecies conditions (4 representative mu values, spanning
    # the cos_to_attractor walk from 0 to 1)
    qs_grid = [0.001, 0.05, 0.1, 0.5]
    for mu in qs_grid:
        label = f"quasispecies_mu_{mu:.3f}"
        rec = fit_one_condition(
            label,
            simulate_quasispecies,
            sim_kwargs={"mu": mu},
            seeds=seeds,
        )
        rec["mu"] = mu
        out["quasispecies"].append(rec)

    out["total_elapsed_s"] = float(time.time() - t_start)

    print("\n=== summary ===")
    print("HYPERCYCLE:")
    for rec in out["hypercycle"]:
        print(f"  k_cat={rec['k_cat']:.2f}")
        for s_rec in rec["per_seed"]:
            print(f"    seed={s_rec['seed']}  best: "
                  f"{s_rec['best_equation_sympy']}")
    print("QUASISPECIES:")
    for rec in out["quasispecies"]:
        print(f"  mu={rec['mu']:.3f}")
        for s_rec in rec["per_seed"]:
            print(f"    seed={s_rec['seed']}  best: "
                  f"{s_rec['best_equation_sympy']}")
    print(f"\nTotal: {out['total_elapsed_s']:.1f}s")

    with open("prebiotic_pysr.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Wrote prebiotic_pysr.json")


if __name__ == "__main__":
    main()
