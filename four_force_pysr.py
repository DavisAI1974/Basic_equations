"""
Session 8 item (1) follow-up: PySR symbolic regression on the four
force-field caricatures.

Procedure: same as pysr_symbolic_per_domain.py but with the four-force
simulators from four_force_probe.py. Fit MI = f(H_a, H_b) on the
ensemble-H operator data for EM / weak / strong / gravity, using
extended operator set {+, -, *, /, square, cube, exp, log, sqrt}.

Tests whether the four forces -- whose [2,3,4] null directions are
nearly identical (KBK probe in four_force_probe) -- differ at the
functional-family level the way the four-domain Session 6/7
simulators did (INFO-025).

Speaking posture (Rule C): the KBK result has the four forces
sharing substrate signature, so a strong prior is that they may
share functional family too. But the per-domain Session 6/7 result
showed functional family can differ even when null direction is
similar -- so we wait on the data.

Output: four_force_pysr_results.json
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix
from four_force_probe import FORCES


warnings.filterwarnings("ignore", category=UserWarning)


def fit_pysr_one(X, y, niter=40, populations=20, maxsize=15, seed=0,
                 timeout_seconds=120):
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
        timeout_in_seconds=timeout_seconds,
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
    return {
        "fit_time_s": fit_t,
        "pareto_front": rows,
        "best_equation_sympy": str(model.sympy()),
    }


def analyze_force(name, simfn, seed, cfg, niter, max_n, timeout_s):
    X1, X2 = simfn(seed=seed, **cfg)
    M, _ = build_ensemble_operator_matrix(X1, X2)
    Ha = M[:, 0]; Hb = M[:, 1]; MI = M[:, 5]
    n = len(MI)
    if n > max_n:
        idx = np.linspace(0, n - 1, max_n).astype(int)
        Ha = Ha[idx]; Hb = Hb[idx]; MI = MI[idx]
        n_used = max_n
    else:
        n_used = n
    X = np.column_stack([Ha, Hb])
    out = fit_pysr_one(X, MI, niter=niter, seed=seed,
                       timeout_seconds=timeout_s)
    out["force"] = name
    out["seed"] = seed
    out["n_used"] = int(n_used)
    out["feature_names"] = ["H_a", "H_b"]
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11, 22])
    p.add_argument("--out", type=str, default=None)
    p.add_argument("--niter", type=int, default=40)
    p.add_argument("--max_n", type=int, default=600)
    p.add_argument("--timeout_s", type=int, default=90)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        out = args.out or "four_force_pysr_canary.json"
        args.niter = 15
        args.max_n = 200
        args.timeout_s = 60
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "four_force_pysr_results.json"

    t0 = time.time()
    print(f"[four_force_pysr] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={args.seeds} niter={args.niter} "
          f"max_n={args.max_n} timeout_s={args.timeout_s}", flush=True)

    all_results = {}
    for seed in args.seeds:
        print(f"\n=== seed={seed} ===", flush=True)
        seed_results = []
        for name, simfn in FORCES.items():
            r = analyze_force(name, simfn, seed, cfg,
                              niter=args.niter,
                              max_n=args.max_n,
                              timeout_s=args.timeout_s)
            print(f"\n  {name:10s}  n={r['n_used']}  "
                  f"fit_time={r['fit_time_s']:.1f}s", flush=True)
            print(f"    best (sympy): {r['best_equation_sympy']}",
                  flush=True)
            print(f"    Pareto front (complexity / loss / "
                  f"score / equation):", flush=True)
            for row in r["pareto_front"][:8]:
                print(f"      c={row['complexity']:3d}  "
                      f"loss={row['loss']:.3e}  "
                      f"score={row['score']:+.3f}  "
                      f"{row['equation']}", flush=True)
            seed_results.append(r)
        all_results[f"seed{seed}"] = seed_results

    blob = {
        "config": cfg, "seeds": args.seeds,
        "niter": args.niter, "max_n": args.max_n,
        "timeout_s": args.timeout_s,
        "results_by_seed": all_results,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[four_force_pysr] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
