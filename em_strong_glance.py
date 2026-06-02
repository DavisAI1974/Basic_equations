"""
Session 8 follow-up: quick glance at EM vs strong force.

Greg's question at session close: "we really should take a quick
glance at em/sf". Direct comparison of EM and strong force
caricatures.

The Session 8 four-force probe (INFO-027) had EM and strong
showing different KBK + PySR signatures, but the comparison was
confounded: EM's simulator had omega1=1.0, omega2=1.2 (channel-
ASYMMETRIC), strong's had omega=0.5 for both (channel-SYMMETRIC).
INFO-028 said PySR cross-seed reproducibility tracks dynamical
symmetry.

This script swaps the symmetry to disambiguate. Four configs:

  (a) EM baseline      omega1=1.0, omega2=1.2  asymmetric    K=0.30
  (b) EM symmetric     omega1=1.0, omega2=1.0  SYMMETRIC      K=0.30
  (c) strong baseline  omega=0.5 (both)         SYMMETRIC      K_conf=0.80
  (d) strong asymmetric omega1=0.5, omega2=0.7  asymmetric    K_conf=0.80

Three seeds each, KBK + PySR. The clean predictions under INFO-028:
  - (a) EM asymmetric    -> PySR cross-seed family reproduces      (control)
  - (b) EM symmetric     -> PySR breaks symmetry randomly per seed (swap)
  - (c) strong symmetric -> PySR breaks symmetry randomly per seed (control)
  - (d) strong asymmetric -> PySR cross-seed family reproduces     (swap)

If the swaps land where predicted, INFO-028 is a clean reading.
If they don't, something else is going on.

Speaking posture (Rule C): I think the symmetry swap will move
reproducibility with the symmetry, but we'll wait on what the
data says.

Output: em_strong_glance_results.json
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix, analyze_M


warnings.filterwarnings("ignore", category=UserWarning)


def simulate_em_config(seed, N_ens=600, T=30.0, dt=0.02,
                        omega1=1.0, omega2=1.2, K=0.30):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    gamma_rad = 0.02
    sigma_noise = 0.02
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        F1 = -(omega1 ** 2) * x1 + K * x2
        F2 = -(omega2 ** 2) * x2 + K * x1
        v1 = v1 + F1 * dt - gamma_rad * v1 * dt
        v2 = v2 + F2 * dt - gamma_rad * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


def simulate_strong_config(seed, N_ens=600, T=30.0, dt=0.02,
                            omega1=0.5, omega2=0.5,
                            K_lin=0.20, K_conf=0.80):
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    gamma = 0.05
    sigma_noise = 0.02
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dx = x1 - x2
        Fconf = K_lin * dx + K_conf * dx ** 3
        F1 = -(omega1 ** 2) * x1 - Fconf
        F2 = -(omega2 ** 2) * x2 + Fconf
        v1 = v1 + F1 * dt - gamma * v1 * dt
        v2 = v2 + F2 * dt - gamma * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


CONFIGS = [
    {"name": "em_asym_baseline",
     "simfn": simulate_em_config,
     "kwargs": dict(omega1=1.0, omega2=1.2, K=0.30),
     "symmetric": False,
     "description": "EM, omega1=1.0 omega2=1.2 (baseline; asymmetric)"},
    {"name": "em_sym_swap",
     "simfn": simulate_em_config,
     "kwargs": dict(omega1=1.0, omega2=1.0, K=0.30),
     "symmetric": True,
     "description": "EM, omega1=omega2=1.0 (symmetry swap)"},
    {"name": "strong_sym_baseline",
     "simfn": simulate_strong_config,
     "kwargs": dict(omega1=0.5, omega2=0.5, K_lin=0.20, K_conf=0.80),
     "symmetric": True,
     "description": "Strong, omega1=omega2=0.5 (baseline; symmetric)"},
    {"name": "strong_asym_swap",
     "simfn": simulate_strong_config,
     "kwargs": dict(omega1=0.5, omega2=0.7, K_lin=0.20, K_conf=0.80),
     "symmetric": False,
     "description": "Strong, omega1=0.5 omega2=0.7 (symmetry swap)"},
]


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


def cos_between(a, b):
    a = np.asarray(a); b = np.asarray(b)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return None
    return float(np.dot(a, b) / (na * nb))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11, 22, 33])
    p.add_argument("--out", type=str, default=None)
    p.add_argument("--niter", type=int, default=40)
    p.add_argument("--timeout_s", type=int, default=90)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        args.niter = 15
        args.timeout_s = 60
        out = args.out or "em_strong_glance_canary.json"
        seeds = args.seeds[:1]
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "em_strong_glance_results.json"
        seeds = args.seeds

    t0 = time.time()
    print(f"[em_strong_glance] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={seeds} niter={args.niter} "
          f"timeout_s={args.timeout_s}", flush=True)

    results = []
    for cfgrec in CONFIGS:
        for seed in seeds:
            label = f"{cfgrec['name']} seed={seed}"
            print(f"\n=== {label} ({cfgrec['description']}) ===",
                  flush=True)
            tc = time.time()
            X1, X2 = cfgrec["simfn"](seed=seed, **cfg, **cfgrec["kwargs"])
            M, _ = build_ensemble_operator_matrix(X1, X2)
            kbk = analyze_M(M, cfgrec["name"])
            Ha = M[:, 0]; Hb = M[:, 1]; MI = M[:, 5]
            n = len(MI)
            n_pysr = min(600, n)
            if n > n_pysr:
                idx = np.linspace(0, n - 1, n_pysr).astype(int)
                Ha_p = Ha[idx]; Hb_p = Hb[idx]; MI_p = MI[idx]
            else:
                Ha_p = Ha; Hb_p = Hb; MI_p = MI
            X = np.column_stack([Ha_p, Hb_p])
            pysr_fit = fit_pysr(X, MI_p, niter=args.niter, seed=seed,
                                 timeout_seconds=args.timeout_s)
            elapsed = time.time() - tc
            r = {
                "config_name": cfgrec["name"],
                "config_description": cfgrec["description"],
                "config_kwargs": cfgrec["kwargs"],
                "is_symmetric_dynamics": cfgrec["symmetric"],
                "seed": seed,
                "elapsed_s": elapsed,
                "n_points": int(n),
                "MI_mean": float(np.mean(MI)),
                "MI_std": float(np.std(MI)),
                "kbk_eigenvalues": kbk["eigenvalues_op_cov"],
                "kbk_v_null_6d": kbk["extract_v1_v_null_6d"],
                "kbk_v_null_234": kbk["extract_v1_v_null_234"],
                "kbk_cos_to_artifact": kbk["cos_to_artifact_+1+1+2"],
                "kbk_signal_rank": kbk["kbk_rank_gap"]["k_signal_kbk"],
                "kbk_null_rank": kbk["kbk_rank_gap"]["rank_null_subspace_kbk"],
                "pysr": pysr_fit,
            }
            results.append(r)
            print(f"  elapsed={elapsed:.1f}s  MI mean={r['MI_mean']:+.3f} "
                  f"std={r['MI_std']:.3f}",
                  flush=True)
            print(f"  KBK v_null 6D = "
                  f"{[f'{v:+.3f}' for v in r['kbk_v_null_6d']]}",
                  flush=True)
            print(f"  KBK v_null [2,3,4] = "
                  f"{[f'{v:+.3f}' for v in r['kbk_v_null_234']]}  "
                  f"cos to (+1,+1,+2) = {r['kbk_cos_to_artifact']:+.3f}",
                  flush=True)
            print(f"  KBK signal_rank={r['kbk_signal_rank']} "
                  f"null_rank={r['kbk_null_rank']} "
                  f"eig_min={r['kbk_eigenvalues'][-1]:.2e}",
                  flush=True)
            print(f"  PySR best (sympy): "
                  f"{r['pysr']['best_equation_sympy']}", flush=True)
            print(f"  PySR Pareto (top 5 low-complexity):", flush=True)
            for row in r["pysr"]["pareto_front"][:5]:
                print(f"    c={row['complexity']:3d} "
                      f"loss={row['loss']:.3e} "
                      f"score={row['score']:+.3f} "
                      f"{row['equation']}", flush=True)

    # Cross-seed comparison: for each config, compute cosine of v_null
    # between seed pairs.
    cross_seed = {}
    by_config = {}
    for r in results:
        by_config.setdefault(r["config_name"], []).append(r)
    for cname, rs in by_config.items():
        seeds_list = [r["seed"] for r in rs]
        cos_pairs = {}
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                ci = cos_between(rs[i]["kbk_v_null_6d"],
                                  rs[j]["kbk_v_null_6d"])
                cos_pairs[f"seed{rs[i]['seed']}__seed{rs[j]['seed']}"] = ci
        cross_seed[cname] = {
            "seeds": seeds_list,
            "kbk_v_null_cos_pairs": cos_pairs,
            "is_symmetric_dynamics": rs[0]["is_symmetric_dynamics"],
            "pysr_best_per_seed": [
                {"seed": r["seed"],
                 "best_sympy": r["pysr"]["best_equation_sympy"],
                 "low_complexity_eqs": [
                     row["equation"] for row in r["pysr"]["pareto_front"][:5]
                 ]} for r in rs
            ],
        }

    print("\n=== Cross-seed comparison ===", flush=True)
    for cname, info in cross_seed.items():
        sym = "SYM" if info["is_symmetric_dynamics"] else "ASYM"
        print(f"  [{sym}] {cname}:", flush=True)
        for pair, c in info["kbk_v_null_cos_pairs"].items():
            print(f"    KBK v_null 6D cos {pair} = {c:+.4f}", flush=True)
        print(f"    PySR best per seed:", flush=True)
        for entry in info["pysr_best_per_seed"]:
            print(f"      seed {entry['seed']}: {entry['best_sympy']}",
                  flush=True)

    blob = {
        "mode": "canary" if args.canary else "full",
        "seeds": seeds, "config": cfg,
        "niter": args.niter, "timeout_s": args.timeout_s,
        "configs": [{k: v for k, v in c.items() if k != "simfn"}
                     for c in CONFIGS],
        "results": results,
        "cross_seed": cross_seed,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[em_strong_glance] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
