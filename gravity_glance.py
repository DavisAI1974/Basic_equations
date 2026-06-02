"""
Gravity glance -- what is gravity really?

Mirror of em_strong_glance.py for the gravity caricature. Three
configs x 3 seeds each.

  (a) gravity_sym_baseline
      omega1=omega2=0.8, universal E_total coupling, G=0.10
      (Session 8 four-force probe baseline)

  (b) gravity_asym_omegas
      omega1=0.8, omega2=1.0, universal E_total coupling, G=0.10
      (symmetry swap; test if breaking channel symmetry exposes
      gravity's structural form, the way it did for strong in
      INFO-031)

  (c) gravity_nonuniversal
      omega1=0.8, omega2=1.0, coupling to ONLY channel 1's energy
      (E1 = x1^2 + v1^2 instead of E_total), G=0.10
      (tests whether "universality" -- coupling to total energy
      regardless of channel -- has a detectable structural
      signature distinct from non-universal coupling)

KBK + PySR per config-seed. Cross-seed comparison + cross-config
comparison.

Speaking posture (Rule C): I think gravity may produce a third
distinct functional family (different from EM polynomial-in-
difference and strong exponential-in-difference), probably
involving the entropy SUM (H_a + H_b) since the coupling depends
on E_total. The universality test should show a different
signature when coupling is restricted to one channel. But we'll
wait on what the data says.

Output: gravity_glance_results.json
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import build_ensemble_operator_matrix, analyze_M


warnings.filterwarnings("ignore", category=UserWarning)


def simulate_gravity_config(seed, N_ens=600, T=30.0, dt=0.02,
                              omega1=0.8, omega2=0.8, G=0.10,
                              universal=True):
    """Gravity caricature with configurable channel masses and
    optional non-universal coupling.

    universal=True  : coupling depends on E_total = E1 + E2
    universal=False : coupling depends only on E1 = x1^2 + v1^2
                       (asymmetric source-energy)
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    gamma = 0.02
    eps_soft = 0.50
    sigma_noise = 0.02
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        E1 = x1 * x1 + v1 * v1
        E2 = x2 * x2 + v2 * v2
        if universal:
            E_src = E1 + E2
        else:
            E_src = E1
        dx = x1 - x2
        r2 = dx * dx + eps_soft
        Fgrav = G * E_src * dx / r2
        F1 = -(omega1 ** 2) * x1 - Fgrav
        F2 = -(omega2 ** 2) * x2 + Fgrav
        v1 = v1 + F1 * dt - gamma * v1 * dt
        v2 = v2 + F2 * dt - gamma * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


CONFIGS = [
    {"name": "gravity_sym_baseline",
     "kwargs": dict(omega1=0.8, omega2=0.8, G=0.10, universal=True),
     "symmetric": True,
     "universal": True,
     "description": "omega1=omega2=0.8, universal E_total coupling (Session 8 baseline)"},
    {"name": "gravity_asym_omegas",
     "kwargs": dict(omega1=0.8, omega2=1.0, G=0.10, universal=True),
     "symmetric": False,
     "universal": True,
     "description": "omega1=0.8 omega2=1.0, universal E_total coupling (symmetry swap)"},
    {"name": "gravity_nonuniversal",
     "kwargs": dict(omega1=0.8, omega2=1.0, G=0.10, universal=False),
     "symmetric": False,
     "universal": False,
     "description": "omega1=0.8 omega2=1.0, NON-universal (only E1) coupling"},
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
        out = args.out or "gravity_glance_canary.json"
        seeds = args.seeds[:1]
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "gravity_glance_results.json"
        seeds = args.seeds

    t0 = time.time()
    print(f"[gravity_glance] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seeds={seeds} niter={args.niter} "
          f"timeout_s={args.timeout_s}", flush=True)

    results = []
    for cfgrec in CONFIGS:
        for seed in seeds:
            label = f"{cfgrec['name']} seed={seed}"
            print(f"\n=== {label} ({cfgrec['description']}) ===",
                  flush=True)
            tc = time.time()
            X1, X2 = simulate_gravity_config(seed=seed, **cfg,
                                              **cfgrec["kwargs"])
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
                "is_universal_coupling": cfgrec["universal"],
                "seed": seed,
                "elapsed_s": elapsed,
                "n_points": int(n),
                "Ha_mean": float(np.mean(Ha)),
                "Hb_mean": float(np.mean(Hb)),
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
            print(f"  elapsed={elapsed:.1f}s  H_a mean={r['Ha_mean']:+.3f} "
                  f"H_b mean={r['Hb_mean']:+.3f}  "
                  f"MI mean={r['MI_mean']:+.3f} std={r['MI_std']:.3f}",
                  flush=True)
            print(f"  KBK v_null 6D = "
                  f"{[f'{v:+.3f}' for v in r['kbk_v_null_6d']]}",
                  flush=True)
            print(f"  KBK v_null [2,3,4] = "
                  f"{[f'{v:+.3f}' for v in r['kbk_v_null_234']]}  "
                  f"cos to (+1,+1,+2) = {r['kbk_cos_to_artifact']:+.3f}",
                  flush=True)
            print(f"  PySR best (sympy): "
                  f"{r['pysr']['best_equation_sympy']}", flush=True)
            print(f"  PySR Pareto (top 6 low-complexity):", flush=True)
            for row in r["pysr"]["pareto_front"][:6]:
                print(f"    c={row['complexity']:3d} "
                      f"loss={row['loss']:.3e} "
                      f"score={row['score']:+.3f} "
                      f"{row['equation']}", flush=True)

    # Cross-seed + cross-config comparison
    by_config = {}
    for r in results:
        by_config.setdefault(r["config_name"], []).append(r)

    cross_seed = {}
    for cname, rs in by_config.items():
        cos_pairs = {}
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                ci = cos_between(rs[i]["kbk_v_null_6d"],
                                  rs[j]["kbk_v_null_6d"])
                cos_pairs[f"seed{rs[i]['seed']}__seed{rs[j]['seed']}"] = ci
        cross_seed[cname] = {
            "seeds": [r["seed"] for r in rs],
            "kbk_v_null_cos_pairs": cos_pairs,
            "is_symmetric_dynamics": rs[0]["is_symmetric_dynamics"],
            "is_universal_coupling": rs[0]["is_universal_coupling"],
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
        univ = "UNIV" if info["is_universal_coupling"] else "NON-UNIV"
        print(f"  [{sym}/{univ}] {cname}:", flush=True)
        for pair, c in info["kbk_v_null_cos_pairs"].items():
            print(f"    KBK v_null 6D cos {pair} = {c:+.4f}", flush=True)
        print(f"    PySR best per seed:", flush=True)
        for entry in info["pysr_best_per_seed"]:
            print(f"      seed {entry['seed']}: {entry['best_sympy']}",
                  flush=True)

    # Cross-config: compare config means of v_null directions
    print("\n=== Cross-config comparison (mean v_null per config) ===",
          flush=True)
    config_means = {}
    for cname, rs in by_config.items():
        v_null_mean = np.mean(
            [r["kbk_v_null_6d"] for r in rs], axis=0
        )
        config_means[cname] = v_null_mean.tolist()
    cnames = list(config_means.keys())
    for i in range(len(cnames)):
        for j in range(i + 1, len(cnames)):
            c = cos_between(config_means[cnames[i]],
                             config_means[cnames[j]])
            print(f"  cos({cnames[i]} <-> {cnames[j]}) = {c:+.4f}",
                  flush=True)

    blob = {
        "mode": "canary" if args.canary else "full",
        "seeds": seeds, "config": cfg,
        "niter": args.niter, "timeout_s": args.timeout_s,
        "configs": [{k: v for k, v in c.items()} for c in CONFIGS],
        "results": results,
        "cross_seed": cross_seed,
        "config_means_v_null": config_means,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[gravity_glance] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
