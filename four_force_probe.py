"""
Session 8 item (1): Four-force unification probe (KBK stack).

Greg's "real test" from Session 7 close: "can we come up with one
equation? Should they be grouped at all?"

Under the new Operating Rule (literature as conjecture by default,
Session 7 close): only EM+weak unification meets the immense-
replicated-data bar; everything else (GUT couplings at 10^16 GeV,
MSSM, gravity-as-different-category, gravity-as-emergent) is
conjecture and competes on roughly equal footing with the
substrate-vs-expression frame.

Approach: build four 2-channel toy simulators capturing the
characteristic structural feature of each force:

  EM       -- long-range bilinear coupling, weak radiation damping
              (massless mediator caricature: linear in channel
              amplitude, no suppression with separation)
  WEAK     -- Yukawa-suppressed bilinear coupling via massive
              mediator (exp(-M*(x1-x2)^2) factor on linear coupling),
              strong damping (W/Z short-lived)
  STRONG   -- confining nonlinear coupling (force grows with
              separation: K_lin*dx + K_conf*dx^3 in restoring
              direction, flux-tube caricature)
  GRAVITY  -- universal energy-density-mediated attractive coupling
              (all "mass-energy" gravitates; F ~ E_total * dx /
              (dx^2 + eps); softened 1/r^2)

Apply the per-domain stack from Sessions 6/7 to each:
  - ensemble-H operator matrix (H_a, H_b, Ha^2, Hb^2, Ha*Hb, MI)
  - KBK rank gap (Kaiser-Brunton-Kutz 2024 SVD-rank-gap)
  - AI Poincare (Cao-Liu-Tegmark 2021 intrinsic-dim)
  - SINDy deg-3 extended library
  - extract_v1 smallest null, [2,3,4] projection
  - cross-force cosines of null directions

Cross-force comparison: which forces share signature?

The toy simulators are caricatures, NOT predictions of QFT. They
capture qualitative coupling structure (range, mediator mass,
nonlinearity, universality). The mapping from this output to real
physics is interpretive and frame-level, not claim-level.

Speaking posture (Rule C): I think EM and weak may share most of
their signature through the bilinear linear coupling, with the
weak case carrying additional suppression-term content; strong
should show distinctive nonlinear cubic content; gravity should
sit apart through its universal energy-density coupling. But
we'll wait on what the data says and where it points us.

Output: four_force_probe_results.json
"""
import argparse
import json
import time
import numpy as np

from per_domain_kbk import (
    build_ensemble_operator_matrix, analyze_M, cross_domain_compare,
)


def simulate_em(seed, N_ens=600, T=30.0, dt=0.02):
    """EM caricature: long-range linear bilinear coupling, weak damping.

    F_i = -omega_i^2 * x_i + K * x_j  (linear coupling, no
    separation-dependent suppression: massless mediator analog)
    Small radiation-reaction damping gamma_rad.
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    omega1 = 1.0
    omega2 = 1.2
    K = 0.30
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


def simulate_weak(seed, N_ens=600, T=30.0, dt=0.02):
    """Weak caricature: Yukawa-suppressed bilinear coupling, strong damping.

    F_i = -omega^2 * x_i + K * x_j * exp(-M * (x_i - x_j)^2)
    The Gaussian suppression in (x_i - x_j) caricatures the massive
    mediator: coupling strong at small field separation, exponentially
    suppressed at large separation. Strong damping (W/Z short-lived).
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    omega = 1.0
    K = 0.50
    M_med = 5.0
    gamma = 0.30
    sigma_noise = 0.03
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dx = x1 - x2
        suppress = np.exp(-M_med * dx * dx)
        F1 = -(omega ** 2) * x1 + K * x2 * suppress
        F2 = -(omega ** 2) * x2 + K * x1 * suppress
        v1 = v1 + F1 * dt - gamma * v1 * dt
        v2 = v2 + F2 * dt - gamma * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


def simulate_strong(seed, N_ens=600, T=30.0, dt=0.02):
    """Strong caricature: confining nonlinear coupling.

    F_i = -omega^2 * x_i -+ (K_lin * dx + K_conf * dx^3)
    where the linear + cubic term acts AS A RESTORING force on the
    separation -- coupling grows with |dx|, flux-tube caricature.
    Asymptotic freedom: weak at small dx (linear); confining at large
    dx (cubic).
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    omega = 0.5
    K_lin = 0.20
    K_conf = 0.80
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
        F1 = -(omega ** 2) * x1 - Fconf
        F2 = -(omega ** 2) * x2 + Fconf
        v1 = v1 + F1 * dt - gamma * v1 * dt
        v2 = v2 + F2 * dt - gamma * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


def simulate_gravity(seed, N_ens=600, T=30.0, dt=0.02):
    """Gravity caricature: universal energy-mediated attractive coupling.

    F_i = -omega^2 * x_i  -+  G * E_total * dx / (dx^2 + eps_soft)
    where E_total = x1^2 + v1^2 + x2^2 + v2^2 (universal mass-energy
    density analog). Force is attractive (acts toward the other
    channel), softened 1/r^2 in 1D, modulated by total energy. Weak
    damping (gravity is conservative).
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    omega = 0.8
    G = 0.10
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
        E_tot = x1 * x1 + v1 * v1 + x2 * x2 + v2 * v2
        dx = x1 - x2
        r2 = dx * dx + eps_soft
        Fgrav = G * E_tot * dx / r2
        F1 = -(omega ** 2) * x1 - Fgrav
        F2 = -(omega ** 2) * x2 + Fgrav
        v1 = v1 + F1 * dt - gamma * v1 * dt
        v2 = v2 + F2 * dt - gamma * v2 * dt
        x1 = x1 + v1 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        x2 = x2 + v2 * dt + sigma_noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = x1
        X2[t] = x2
    return X1, X2


FORCES = {
    "em":       simulate_em,
    "weak":     simulate_weak,
    "strong":   simulate_strong,
    "gravity":  simulate_gravity,
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        out = args.out or "four_force_probe_canary.json"
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "four_force_probe_results.json"

    t0 = time.time()
    seed = args.seed
    print(f"[four_force_probe] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seed={seed}", flush=True)

    results = []
    for name, simfn in FORCES.items():
        ts = time.time()
        X1, X2 = simfn(seed=seed, **cfg)
        sim_t = time.time() - ts
        M, t_idx = build_ensemble_operator_matrix(X1, X2)
        op_t = time.time() - ts - sim_t
        r = analyze_M(M, name)
        elapsed = time.time() - ts
        r["sim_time_s"] = sim_t
        r["op_time_s"] = op_t
        r["analyze_time_s"] = elapsed - sim_t - op_t
        results.append(r)
        print(f"\n  force={name}  elapsed={elapsed:.1f}s  "
              f"(sim={sim_t:.1f}s op={op_t:.1f}s)  "
              f"n_points={r['n_points']}",
              flush=True)
        print(f"    op_means = "
              f"{[f'{m:+.3f}' for m in r['operator_means']]}",
              flush=True)
        print(f"    SV deg2 = "
              f"{[f'{s:.3e}' for s in r['singular_values_deg2']]}",
              flush=True)
        print(f"    eigvals = "
              f"{[f'{e:.3e}' for e in r['eigenvalues_op_cov']]}",
              flush=True)
        print(f"    KBK signal_rank={r['kbk_rank_gap']['k_signal_kbk']}  "
              f"null_rank={r['kbk_rank_gap']['rank_null_subspace_kbk']}  "
              f"gap_ratio={r['kbk_rank_gap']['largest_gap_ratio']:.2e}",
              flush=True)
        print(f"    extract_v1 v_null 6D = "
              f"{[f'{v:+.4f}' for v in r['extract_v1_v_null_6d']]}",
              flush=True)
        print(f"    [2,3,4] projection = "
              f"{[f'{v:+.4f}' for v in r['extract_v1_v_null_234']]}  "
              f"cos to (+,+,+) = "
              f"{r['cos_to_artifact_+1+1+2']:+.4f}",
              flush=True)
        if r["ai_poincare_two_nn_z"] is not None:
            print(f"    AI Poincare two-NN(z) = "
                  f"{r['ai_poincare_two_nn_z']:.3f}", flush=True)
        if r["ai_poincare_levina_bickel"]["dim"] is not None:
            print(f"    Levina-Bickel dim = "
                  f"{r['ai_poincare_levina_bickel']['dim']:.3f} "
                  f"+/- {r['ai_poincare_levina_bickel']['std']:.3f}",
                  flush=True)
        for k, v in r["ai_poincare_local_pca"].items():
            print(f"    local-PCA k={k} mean={v['mean']:.2f} "
                  f"median={v['median']:.1f} std={v['std']:.2f}",
                  flush=True)
        print(f"    KBK symbolic relations (sparse):", flush=True)
        for sym in r["kbk_symbolic_relations"][:3]:
            print(f"      null[{sym['null_index']}] sparse: "
                  f"{sym['eq_sparse']}  "
                  f"resid_std={sym['resid_std_sparse']:.3e}",
                  flush=True)

    print("\n  cross-force cosines (smallest extract_v1 null):",
          flush=True)
    cd = cross_domain_compare(results)
    for k, v in cd["cos_6d"].items():
        print(f"    6D  cos({k}) = {v:+.4f}", flush=True)
    for k, v in cd["cos_234"].items():
        print(f"    234 cos({k}) = {v:+.4f}", flush=True)

    out_blob = {"config": cfg, "seed": seed,
                "results": results, "cross_force": cd,
                "elapsed_total_s": time.time() - t0}
    with open(out, "w") as f:
        json.dump(out_blob, f, indent=2)
    print(f"\n[four_force_probe] wrote {out} "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
