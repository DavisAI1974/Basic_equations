"""
Session 7 / item (2): PySR symbolic regression on per-domain
ensemble-H data. Tests v5/Session-6's per-domain functional forms
against arbitrary nonlinear shapes that linear-in-features SINDy
cannot find (ratios, exponentials, logs, sqrt, composition).

Per the v5 literature scan, "no SINDy/AI-Feynman/PySR/Eureqa
application to windowed differential entropy of coupled species"
is a clean unstacked target. This script fits PySR to MI = f(H_a, H_b)
for each of the 4 simulated domains.

Why not just SINDy: SINDy needs a feature library specified up front,
so it can only find combinations of features that are pre-listed.
The Session 6 SINDy run found polynomial-in-H content; PySR can
discover ratio / exponential / log structure.

If install fails or runtime explodes, we fall back to gplearn.

Output: pysr_symbolic_per_domain_results.json

Speaking posture (Rule C): I think PySR may find a simple polynomial
or near-linear form similar to what SINDy already found, but I'm
not pre-committing -- if it surfaces a ratio or exponential that
fits better at lower complexity, that's a different signal worth
mapping.
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import (
    DOMAINS, build_ensemble_operator_matrix,
)


warnings.filterwarnings("ignore", category=UserWarning)


def fit_pysr_one(X, y, niter=40, populations=20, maxsize=15, seed=0,
                 timeout_seconds=120):
    """Fit a PySR regressor on (X, y). Returns best-equation summaries
    sorted by complexity, plus the chosen 'best' equation by PySR's
    own scoring."""
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
        # nan/inf-safe defaults; complexity weighting favors simpler
        # eqs when fits are comparable
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
        "fit_time_s": fit_t,
        "pareto_front": rows,
        "best_equation_sympy": best_eq,
    }


def analyze_domain(name, simfn, seed, cfg, niter, max_n,
                   timeout_s):
    """Run PySR on the per-domain (H_a, H_b) -> MI mapping."""
    X1, X2 = simfn(seed=seed, **cfg)
    M, _ = build_ensemble_operator_matrix(X1, X2)
    Ha = M[:, 0]; Hb = M[:, 1]; MI = M[:, 5]
    n = len(MI)
    # uniform subsample to max_n if needed (PySR scales linearly in n
    # per evaluation but evaluations are many; ~500-1000 is plenty)
    if n > max_n:
        idx = np.linspace(0, n - 1, max_n).astype(int)
        Ha = Ha[idx]; Hb = Hb[idx]; MI = MI[idx]
        n_used = max_n
    else:
        n_used = n
    X = np.column_stack([Ha, Hb])
    feat_names = ["H_a", "H_b"]

    out = fit_pysr_one(X, MI, niter=niter, seed=seed,
                       timeout_seconds=timeout_s)
    out["domain"] = name
    out["seed"] = seed
    out["n_used"] = int(n_used)
    out["feature_names"] = feat_names
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
        out = args.out or "pysr_symbolic_per_domain_canary.json"
        args.niter = 15
        args.max_n = 200
        args.timeout_s = 60
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "pysr_symbolic_per_domain_results.json"

    t0 = time.time()
    print(f"[pysr_symbolic] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={args.seeds} niter={args.niter} "
          f"max_n={args.max_n} timeout_s={args.timeout_s}", flush=True)

    all_results = {}
    for seed in args.seeds:
        print(f"\n=== seed={seed} ===", flush=True)
        seed_results = []
        for name, simfn in DOMAINS.items():
            ts = time.time()
            r = analyze_domain(name, simfn, seed, cfg,
                               niter=args.niter,
                               max_n=args.max_n,
                               timeout_s=args.timeout_s)
            print(f"\n  {name:22s}  n={r['n_used']}  "
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
    print(f"\n[pysr_symbolic] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
