"""
Session 8 item (4): Robustness check on INFO-025 biology finding.

INFO-025 (Session 7): biology Lotka-Volterra MI = 0.5 * exp(H_a/2)
at complexity 5, cross-seed coefficient variation <5%. Tested under
baseline conditions T=30, dt=0.02, N_ens=600.

Open questions from Session 7 close (queued #5 in v8 kickoff):
  - Does biology stay exponential at T=10? T=100?
  - Under Gaussian observation noise added to the time series?
  - Under different ensemble sizes (N_ens=200, 1200)?

This script runs PySR (same configuration as Session 7) on
biology data generated under a grid of robustness conditions and
checks whether the exp(H_a/2) family survives.

Speaking posture (Rule C): I think the exponential family will
survive small noise and large T but may degrade at very short T
(insufficient regime exploration) or under heavy observational
noise (signal washed out). But we'll wait on what the data says
and where it points us.

Output: robustness_info025_results.json
"""
import argparse
import json
import time
import warnings
import numpy as np

from per_domain_kbk import simulate_biology, build_ensemble_operator_matrix


warnings.filterwarnings("ignore", category=UserWarning)


def simulate_biology_noisy(seed, N_ens=600, T=30.0, dt=0.02,
                            obs_noise_sigma=0.0):
    """Simulate biology Lotka-Volterra and optionally add Gaussian
    observation noise to the trajectories before computing entropies.
    obs_noise_sigma is in units of the channel std (so 0.1 = 10% of
    channel std as added noise)."""
    X1, X2 = simulate_biology(seed=seed, N_ens=N_ens, T=T, dt=dt)
    if obs_noise_sigma > 0:
        rng = np.random.default_rng(seed + 99991)
        std1 = float(np.std(X1))
        std2 = float(np.std(X2))
        X1 = X1 + obs_noise_sigma * std1 * rng.standard_normal(X1.shape)
        X2 = X2 + obs_noise_sigma * std2 * rng.standard_normal(X2.shape)
    return X1, X2


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


def fit_exp_baseline(Ha, Hb, MI):
    """Fit MI = A * exp(H_a / B) directly and report R^2.
    Anchor against INFO-025 expectation: A ~ 0.5, B ~ 2.0.
    """
    Ha = np.asarray(Ha); MI = np.asarray(MI)
    from scipy.optimize import curve_fit
    def model(h, A, B):
        return A * np.exp(h / B)
    try:
        popt, _ = curve_fit(model, Ha, MI, p0=[0.5, 2.0], maxfev=5000)
        pred = model(Ha, *popt)
        ss_res = float(np.sum((MI - pred) ** 2))
        ss_tot = float(np.sum((MI - MI.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        return {"A": float(popt[0]), "B": float(popt[1]), "R2": r2}
    except Exception as e:
        return {"A": None, "B": None, "R2": None, "error": str(e)}


def fit_linear_baseline(Ha, MI):
    """Fit MI = a * H_a + b directly (the polynomial baseline that
    PySR might prefer in degraded regimes)."""
    Ha = np.asarray(Ha); MI = np.asarray(MI)
    a, b = np.polyfit(Ha, MI, 1)
    pred = a * Ha + b
    ss_res = float(np.sum((MI - pred) ** 2))
    ss_tot = float(np.sum((MI - MI.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"a": float(a), "b": float(b), "R2": r2}


def run_condition(seed, T, N_ens, obs_noise, niter, timeout_s,
                   max_n_pysr=600):
    """Generate biology data under the given condition, build
    ensemble operator matrix, run PySR and the exp/linear baselines.
    """
    t_sim = time.time()
    X1, X2 = simulate_biology_noisy(seed=seed, N_ens=N_ens, T=T,
                                      dt=0.02,
                                      obs_noise_sigma=obs_noise)
    sim_t = time.time() - t_sim
    M, _ = build_ensemble_operator_matrix(X1, X2)
    Ha = M[:, 0]; Hb = M[:, 1]; MI = M[:, 5]
    n = len(MI)
    if n > max_n_pysr:
        idx = np.linspace(0, n - 1, max_n_pysr).astype(int)
        Ha_p = Ha[idx]; Hb_p = Hb[idx]; MI_p = MI[idx]
        n_used = max_n_pysr
    else:
        Ha_p = Ha; Hb_p = Hb; MI_p = MI
        n_used = n
    X = np.column_stack([Ha_p, Hb_p])

    exp_fit = fit_exp_baseline(Ha, Hb, MI)
    lin_fit = fit_linear_baseline(Ha, MI)
    pysr_fit = fit_pysr(X, MI_p, niter=niter, seed=seed,
                         timeout_seconds=timeout_s)
    return {
        "seed": seed, "T": T, "N_ens": N_ens, "obs_noise": obs_noise,
        "n_points": int(n), "n_pysr": int(n_used),
        "sim_time_s": sim_t,
        "Ha_mean": float(np.mean(Ha)), "Ha_std": float(np.std(Ha)),
        "Hb_mean": float(np.mean(Hb)), "Hb_std": float(np.std(Hb)),
        "MI_mean": float(np.mean(MI)), "MI_std": float(np.std(MI)),
        "exp_baseline_A_exp_HaOverB": exp_fit,
        "linear_baseline": lin_fit,
        "pysr": pysr_fit,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11, 22])
    p.add_argument("--out", type=str, default=None)
    p.add_argument("--niter", type=int, default=40)
    p.add_argument("--timeout_s", type=int, default=90)
    args = p.parse_args()

    if args.canary:
        # Minimal grid: T=10 and T=30 only, no noise, baseline N_ens
        conditions = [
            dict(T=10.0, N_ens=200, obs_noise=0.0),
            dict(T=30.0, N_ens=200, obs_noise=0.0),
        ]
        args.niter = 15
        args.timeout_s = 60
        out = args.out or "robustness_info025_canary.json"
        seeds = args.seeds[:1]
    else:
        # Robustness grid
        conditions = [
            # T sweep (N_ens=600 baseline, no noise)
            dict(T=10.0,  N_ens=600,  obs_noise=0.0),
            dict(T=30.0,  N_ens=600,  obs_noise=0.0),   # baseline
            dict(T=100.0, N_ens=600,  obs_noise=0.0),
            # N_ens sweep (T=30 baseline)
            dict(T=30.0,  N_ens=200,  obs_noise=0.0),
            dict(T=30.0,  N_ens=1200, obs_noise=0.0),
            # Noise sweep (T=30, N_ens=600)
            dict(T=30.0,  N_ens=600,  obs_noise=0.1),
            dict(T=30.0,  N_ens=600,  obs_noise=0.5),
        ]
        out = args.out or "robustness_info025_results.json"
        seeds = args.seeds

    t0 = time.time()
    print(f"[robustness_info025] mode={'canary' if args.canary else 'full'} "
          f"seeds={seeds} niter={args.niter} "
          f"timeout_s={args.timeout_s}", flush=True)
    print(f"  conditions: {len(conditions)}", flush=True)

    results = []
    for seed in seeds:
        for cond in conditions:
            label = (f"T={cond['T']:.0f} N_ens={cond['N_ens']} "
                     f"noise={cond['obs_noise']:.2f}")
            print(f"\n--- seed={seed} {label} ---", flush=True)
            tc = time.time()
            r = run_condition(
                seed=seed,
                T=cond["T"], N_ens=cond["N_ens"],
                obs_noise=cond["obs_noise"],
                niter=args.niter, timeout_s=args.timeout_s,
            )
            elapsed = time.time() - tc
            print(f"  elapsed={elapsed:.1f}s  n_points={r['n_points']}  "
                  f"n_pysr={r['n_pysr']}", flush=True)
            print(f"  H_a mean={r['Ha_mean']:+.3f} std={r['Ha_std']:.3f}  "
                  f"MI mean={r['MI_mean']:+.3f} std={r['MI_std']:.3f}",
                  flush=True)
            eb = r["exp_baseline_A_exp_HaOverB"]
            if eb.get("R2") is not None:
                print(f"  exp baseline (A*exp(H_a/B)): "
                      f"A={eb['A']:+.4f} B={eb['B']:+.4f} R^2={eb['R2']:+.4f}",
                      flush=True)
            lb = r["linear_baseline"]
            print(f"  linear baseline (a*H_a+b): "
                  f"a={lb['a']:+.4f} b={lb['b']:+.4f} R^2={lb['R2']:+.4f}",
                  flush=True)
            print(f"  PySR best (sympy): {r['pysr']['best_equation_sympy']}",
                  flush=True)
            print(f"  PySR Pareto (low-complexity rows):", flush=True)
            for row in r["pysr"]["pareto_front"][:6]:
                print(f"    c={row['complexity']:3d} "
                      f"loss={row['loss']:.3e} score={row['score']:+.3f} "
                      f"{row['equation']}", flush=True)
            results.append(r)

    blob = {
        "mode": "canary" if args.canary else "full",
        "seeds": seeds, "conditions": conditions,
        "niter": args.niter, "timeout_s": args.timeout_s,
        "results": results,
        "elapsed_total_s": time.time() - t0,
    }
    with open(out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"\n[robustness_info025] wrote {out}  "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
