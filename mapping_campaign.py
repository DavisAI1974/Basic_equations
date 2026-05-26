"""
Session 8 item (3): Mapping campaign — which dynamical property
predicts which PySR functional family.

Natural follow-up to INFO-025 (Session 7): four reproducible per-
domain MI-vs-H functional families surfaced by PySR on per-domain
ensemble-H data. The next question: WHICH simulator property
selects WHICH family.

Vary one dominant knob per domain (3 values per knob), run PySR on
each. Track whether the INFO-025 family persists, mutates, or
disappears as the knob moves.

Knobs (one per domain to keep grid bounded):
  - physics_duffing:   K  (neighbor coupling)         [0.05, 0.20, 0.50]
  - biology_lotka:     beta (predation rate)          [0.30, 0.50, 0.80]
  - chemistry_brussel: B  (autocatalytic parameter)    [2.0, 3.0, 4.0]
  - geology_bk:        drift_rate                      [0.02, 0.04, 0.08]

Baseline values (matching per_domain_kbk.py originals) sit in the
middle of each range. So the sweep is +/- a factor 2-3 around baseline.

Speaking posture (Rule C): I think the qualitative family (quadratic
difference, exponential, linear, constant) will persist within domain
across the sweep, with coefficients drifting; very large or small
knob values may break the family entirely. But we'll wait on what
the data says.

Output: mapping_campaign_results.json
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix


warnings.filterwarnings("ignore", category=UserWarning)


# Custom parameterized simulators
def simulate_physics_K(seed, N_ens=600, T=30.0, dt=0.02, K=0.2):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        F1 = -x1 - 0.5 * x1 ** 3 + K * (x2 - x1) + 0.1 * np.sin(0.5 * t * dt)
        F2 = -1.2 * x2 - 0.4 * x2 ** 3 + K * (x1 - x2) + 0.08 * np.cos(0.4 * t * dt)
        v1 = v1 + F1 * dt - 0.05 * v1 * dt
        v2 = v2 + F2 * dt - 0.05 * v2 * dt
        x1 = x1 + v1 * dt + 0.02 * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + 0.02 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


def simulate_biology_beta(seed, N_ens=600, T=30.0, dt=0.02, beta=0.5):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    prey = 1.0 + 0.1 * rng.standard_normal(N_ens)
    pred = 0.5 + 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dprey = (1.0 * prey - beta * prey * pred) * dt
        dpred = (0.3 * prey * pred - 0.4 * pred) * dt
        prey = np.clip(prey + dprey + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt),
                        0.01, None)
        pred = np.clip(pred + dpred + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt),
                        0.01, None)
        X1[t] = prey
        X2[t] = pred
    return X1, X2


def simulate_chemistry_B(seed, N_ens=600, T=30.0, dt=0.02, B=3.0):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    X = 1.0 + 0.1 * rng.standard_normal(N_ens)
    Y = 3.0 + 0.1 * rng.standard_normal(N_ens)
    A_ = 1.0
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dX = (A_ - (B + 1) * X + X ** 2 * Y) * dt
        dY = (B * X - X ** 2 * Y) * dt
        X = X + dX + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        Y = Y + dY + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = X
        X2[t] = Y
    return X1, X2


def simulate_geology_drift(seed, N_ens=600, T=30.0, dt=0.02, drift_rate=0.04):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    u1 = 0.1 * rng.standard_normal(N_ens)
    v1 = 0.05 * rng.standard_normal(N_ens)
    u2 = 0.1 * rng.standard_normal(N_ens)
    v2 = 0.05 * rng.standard_normal(N_ens)
    K = 0.5
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        fric1 = -0.6 * v1 / (1.0 + 8.0 * v1 ** 2)
        fric2 = -0.6 * v2 / (1.0 + 8.0 * v2 ** 2)
        F1 = K * (u2 - u1) - 0.6 * (u1 - drift_rate * t * dt) + fric1
        F2 = K * (u1 - u2) - 0.6 * (u2 - drift_rate * t * dt) + fric2
        v1 = v1 + F1 * dt
        v2 = v2 + F2 * dt
        u1 = u1 + v1 * dt + 0.015 * rng.standard_normal(N_ens) * np.sqrt(dt)
        u2 = u2 + v2 * dt + 0.015 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = u1
        X2[t] = u2
    return X1, X2


KNOB_SWEEPS = {
    "physics_duffing":   {"knob": "K",          "values": [0.05, 0.20, 0.50], "simfn": simulate_physics_K,    "baseline": 0.20},
    "biology_lotka":     {"knob": "beta",       "values": [0.30, 0.50, 0.80], "simfn": simulate_biology_beta, "baseline": 0.50},
    "chemistry_brussel": {"knob": "B",          "values": [2.0, 3.0, 4.0],    "simfn": simulate_chemistry_B,  "baseline": 3.0},
    "geology_bk":        {"knob": "drift_rate", "values": [0.02, 0.04, 0.08], "simfn": simulate_geology_drift,"baseline": 0.04},
}


def fit_pysr(X, y, niter=40, populations=20, maxsize=15, seed=0,
             timeout_seconds=90):
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


def run_one(domain, knob, knob_val, simfn, seed, cfg, niter, timeout_s,
             max_n=600):
    kwargs = dict(seed=seed, **cfg)
    kwargs[knob] = knob_val
    X1, X2 = simfn(**kwargs)
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
    pysr_fit = fit_pysr(X, MI, niter=niter, seed=seed,
                         timeout_seconds=timeout_s)
    return {
        "domain": domain, "knob": knob, "knob_val": float(knob_val),
        "seed": seed,
        "n_pysr": int(n_used),
        "MI_mean": float(np.mean(M[:, 5])),
        "MI_std": float(np.std(M[:, 5])),
        "Ha_mean": float(np.mean(M[:, 0])),
        "Ha_std": float(np.std(M[:, 0])),
        "Hb_mean": float(np.mean(M[:, 1])),
        "Hb_std": float(np.std(M[:, 1])),
        "pysr": pysr_fit,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11])
    p.add_argument("--out", type=str, default=None)
    p.add_argument("--niter", type=int, default=40)
    p.add_argument("--timeout_s", type=int, default=90)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        args.niter = 15
        args.timeout_s = 60
        out = args.out or "mapping_campaign_canary.json"
        # canary: physics only, two knob values
        sweeps = {"physics_duffing": {
            **KNOB_SWEEPS["physics_duffing"],
            "values": [0.05, 0.50],
        }}
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "mapping_campaign_results.json"
        sweeps = KNOB_SWEEPS

    t0 = time.time()
    print(f"[mapping_campaign] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={args.seeds} niter={args.niter} "
          f"timeout_s={args.timeout_s}", flush=True)

    results = []
    for seed in args.seeds:
        for domain, spec in sweeps.items():
            knob = spec["knob"]
            for val in spec["values"]:
                label = f"{domain} {knob}={val}"
                print(f"\n--- seed={seed} {label} (baseline={spec['baseline']}) ---",
                      flush=True)
                tc = time.time()
                r = run_one(domain, knob, val, spec["simfn"], seed,
                             cfg, args.niter, args.timeout_s)
                elapsed = time.time() - tc
                r["elapsed_s"] = elapsed
                r["is_baseline"] = bool(val == spec["baseline"])
                results.append(r)
                print(f"  elapsed={elapsed:.1f}s  n_pysr={r['n_pysr']}",
                      flush=True)
                print(f"  H_a mean={r['Ha_mean']:+.3f} std={r['Ha_std']:.3f}  "
                      f"MI mean={r['MI_mean']:+.3f} std={r['MI_std']:.3f}",
                      flush=True)
                print(f"  PySR best (sympy): "
                      f"{r['pysr']['best_equation_sympy']}", flush=True)
                print(f"  PySR Pareto (top 6 low-complexity):", flush=True)
                for row in r["pysr"]["pareto_front"][:6]:
                    print(f"    c={row['complexity']:3d} "
                          f"loss={row['loss']:.3e} "
                          f"score={row['score']:+.3f} "
                          f"{row['equation']}", flush=True)

    blob = {
        "mode": "canary" if args.canary else "full",
        "seeds": args.seeds, "config": cfg,
        "niter": args.niter, "timeout_s": args.timeout_s,
        "sweeps": {k: {"knob": v["knob"], "values": v["values"],
                        "baseline": v["baseline"]}
                    for k, v in sweeps.items()},
        "results": results,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[mapping_campaign] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
