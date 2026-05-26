"""
KBK + AI Poincare + SINDy on per-domain ensemble-H operator data.

Different procedure from kbk_pipeline.py (which uses single-
trajectory windowed Vasicek H of OU coupled to itself). Here:
ensemble of N_ens trajectories, entropy computed across the
ensemble distribution at each time step. This is the procedure
that produced the v5 per-domain algebraic coefficients:

  Chemistry:  H_a^2 = 0.007 - 0.093*(H_a*H_b)
              + 1.309*(H_a*H_b)^2,  R^2 = 0.943
  Geology:    0.724*(H_a*H_b) - 0.441*H_b^2
              - 0.290*H_a^2 ~ 0,  std/mean = 0.15%
  Biology:    MI ~ polynomial(H_a) with H_b absent, R^2 = 0.66

Per the v5 handoff (Rule D): these per-domain coefficients are
preserved (not affected by Session 5's operator-extraction
artifact). Greg explicitly said they should be subjected to the
same probe discipline. This is that probe.

Four domains, from the 2-system simulators in static_dipole_test.py
re-implemented here to expose the raw ensemble trajectory:

  PHYSICS:    Coupled Duffing oscillators with neighbor coupling
              (slightly anharmonic + nonlinear damping)
  BIOLOGY:    Coupled Lotka-Volterra (prey-predator)
  CHEMISTRY:  Coupled Brusselator (A=1, B=3, two-step ODE)
  GEOLOGY:   Coupled Burridge-Knopoff slip-stick (rate-dependent
             friction, slow drift)

For each domain: at each time step t, compute H_a(t), H_b(t),
MI(t) across the ensemble distribution. Stack into operator
matrix M with rows = time, columns = [H_a, H_b, H_a^2, H_b^2,
H_a*H_b, MI]. Apply KBK rank-gap, AI Poincare local-PCA,
extract_v1 (smallest null), and a sparse SINDy with the deg-3
library.

Cross-domain comparison: cosine between smallest null directions.
Project each to [2,3,4] subspace; cos to (+1,+1,+2)/sqrt(6).
Domain content shows up if the per-domain null directions differ
from each other and from the OU baseline.

Output: per_domain_kbk_results.json
"""
import numpy as np
import json
import time
import argparse
from kbk_pipeline import (
    mi_hist_2d, extract_v1, kbk_rank_gap, kbk_symbolic, project_234,
    OP_NAMES,
)
from ai_poincare_rank import (
    two_nn_dim, local_pca_rank_curve, levina_bickel_mle,
)
from sindy_symbolic import (
    extend_library, find_null_directions, format_eq,
)


def entropy_hist_1d(samples, bins=30):
    """Histogram differential entropy estimator (nats)."""
    h, e = np.histogram(samples, bins=bins, density=True)
    h = h[h > 0]
    dx = e[1] - e[0]
    return -float(np.sum(h * np.log(h)) * dx)


def simulate_physics(seed, N_ens=600, T=30.0, dt=0.02):
    """Coupled Duffing oscillators."""
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    x1 = 0.5 * rng.standard_normal(N_ens)
    v1 = 0.1 * rng.standard_normal(N_ens)
    x2 = 0.5 * rng.standard_normal(N_ens)
    v2 = 0.1 * rng.standard_normal(N_ens)
    K = 0.2
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


def simulate_biology(seed, N_ens=600, T=30.0, dt=0.02):
    """Lotka-Volterra prey-predator ensemble."""
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    prey = 1.0 + 0.1 * rng.standard_normal(N_ens)
    pred = 0.5 + 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dprey = (1.0 * prey - 0.5 * prey * pred) * dt
        dpred = (0.3 * prey * pred - 0.4 * pred) * dt
        prey = np.clip(
            prey + dprey
            + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt),
            0.01, None,
        )
        pred = np.clip(
            pred + dpred
            + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt),
            0.01, None,
        )
        X1[t] = prey
        X2[t] = pred
    return X1, X2


def simulate_chemistry(seed, N_ens=600, T=30.0, dt=0.02):
    """Brusselator."""
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    X = 1.0 + 0.1 * rng.standard_normal(N_ens)
    Y = 3.0 + 0.1 * rng.standard_normal(N_ens)
    A_, B_ = 1.0, 3.0
    X1 = np.zeros((steps, N_ens))
    X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dX = (A_ - (B_ + 1) * X + X ** 2 * Y) * dt
        dY = (B_ * X - X ** 2 * Y) * dt
        X = X + dX + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        Y = Y + dY + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = X
        X2[t] = Y
    return X1, X2


def simulate_geology(seed, N_ens=600, T=30.0, dt=0.02):
    """Burridge-Knopoff slip-stick with rate-dependent friction
    and slow tectonic drift."""
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
        F1 = K * (u2 - u1) - 0.6 * (u1 - 0.04 * t * dt) + fric1
        F2 = K * (u1 - u2) - 0.6 * (u2 - 0.04 * t * dt) + fric2
        v1 = v1 + F1 * dt
        v2 = v2 + F2 * dt
        u1 = u1 + v1 * dt + 0.015 * rng.standard_normal(N_ens) * np.sqrt(dt)
        u2 = u2 + v2 * dt + 0.015 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = u1
        X2[t] = u2
    return X1, X2


DOMAINS = {
    "physics_duffing":    simulate_physics,
    "biology_lotka":      simulate_biology,
    "chemistry_brussel":  simulate_chemistry,
    "geology_bk":         simulate_geology,
}


def build_ensemble_operator_matrix(X1, X2, bins_H=30, bins_MI=18,
                                    burn_in_frac=0.2):
    """Per time-step entropy across the ensemble.

    Returns (M, t_idx) where M is (T_steps_after_burnin, 6).
    """
    T_steps = X1.shape[0]
    t_start = int(burn_in_frac * T_steps)
    rows = []
    t_idx = []
    for t in range(t_start, T_steps):
        Ha = entropy_hist_1d(X1[t], bins=bins_H)
        Hb = entropy_hist_1d(X2[t], bins=bins_H)
        MI = mi_hist_2d(X1[t], X2[t], bins=bins_MI)
        rows.append([Ha, Hb, Ha * Ha, Hb * Hb, Ha * Hb, MI])
        t_idx.append(t)
    return np.array(rows), np.array(t_idx)


def analyze_M(M, name, n_anchors=200):
    """KBK + AI Poincare + sparse symbolic on M."""
    n_pts = M.shape[0]
    mu = M.mean(axis=0)
    sd = M.std(axis=0) + 1e-12
    Mc = M - mu
    v_null, S, _, _ = extract_v1(M)
    eigs = (S ** 2) / max(1, len(Mc) - 1)
    kbk = kbk_rank_gap(S)
    n_keep = max(3, min(kbk["rank_null_subspace_kbk"], 5))
    U, _, Vt = np.linalg.svd(Mc, full_matrices=False)
    syms = kbk_symbolic(Mc, Vt, n_keep)
    sub234, cos_a = project_234(v_null)

    Xs = (M - mu) / sd
    # n_pts may be small for per-domain (~1200 after burn-in). Use
    # k_list that fits.
    k_list = [10, 25, 50, 100]
    k_list = [k for k in k_list if k < n_pts]
    lpca = (
        local_pca_rank_curve(Xs, k_list, var_thresh=0.95,
                              n_anchors=min(n_anchors, n_pts))
        if k_list else {}
    )
    d_2nn, _ = two_nn_dim(Xs) if n_pts >= 5 else (None, 0)
    try:
        d_lb, std_lb, _ = levina_bickel_mle(Xs, k=min(10, n_pts - 2))
    except Exception:
        d_lb, std_lb = None, None

    # extended-library SINDy on the ensemble M (treat as deg-3 lib)
    M_ext, names_ext = extend_library(M, mode="deg3")
    mu_ext = M_ext.mean(axis=0)
    Mc_ext = M_ext - mu_ext
    nulls_ext = find_null_directions(Mc_ext, n_keep=3,
                                       stls_thresh=0.05)

    return {
        "name": name,
        "n_points": int(n_pts),
        "operator_means": mu.tolist(),
        "operator_stds": sd.tolist(),
        "singular_values_deg2": S.tolist(),
        "eigenvalues_op_cov": eigs.tolist(),
        "kbk_rank_gap": kbk,
        "extract_v1_v_null_6d": v_null.tolist(),
        "extract_v1_v_null_234": sub234,
        "cos_to_artifact_+1+1+2": cos_a,
        "kbk_symbolic_relations": syms,
        "ai_poincare_two_nn_z": d_2nn,
        "ai_poincare_levina_bickel": {
            "dim": d_lb, "std": std_lb,
        },
        "ai_poincare_local_pca": lpca,
        "sindy_deg3_feature_names": names_ext,
        "sindy_deg3_nulls": nulls_ext,
    }


def cos_between(a, b):
    a = np.array(a); b = np.array(b)
    na = np.linalg.norm(a); nb = np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return None
    return float(np.dot(a, b) / (na * nb))


def cross_domain_compare(results):
    names = [r["name"] for r in results]
    cos_6d = {}
    cos_234 = {}
    for i, ri in enumerate(results):
        for j, rj in enumerate(results):
            if j <= i:
                continue
            c6 = cos_between(ri["extract_v1_v_null_6d"],
                              rj["extract_v1_v_null_6d"])
            c3 = cos_between(ri["extract_v1_v_null_234"],
                              rj["extract_v1_v_null_234"])
            cos_6d[f"{ri['name']}__{rj['name']}"] = c6
            cos_234[f"{ri['name']}__{rj['name']}"] = c3
    return {"cos_6d": cos_6d, "cos_234": cos_234}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        out = args.out or "per_domain_kbk_canary.json"
    else:
        cfg = dict(N_ens=600, T=30.0, dt=0.02)
        out = args.out or "per_domain_kbk_results.json"

    t0 = time.time()
    seed = args.seed
    print(f"[per_domain_kbk] mode={'canary' if args.canary else 'full'} "
          f"cfg={cfg} seed={seed}", flush=True)

    results = []
    for name, simfn in DOMAINS.items():
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
        print(f"\n  domain={name}  elapsed={elapsed:.1f}s  "
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
        print(f"    SINDy deg-3 nulls (smallest 2):", flush=True)
        for n in r["sindy_deg3_nulls"][:2]:
            eq = format_eq(n["v_sparse"], r["sindy_deg3_feature_names"],
                           prec=3)
            print(f"      null[{n['null_index']}] sv={n['sv']:.3e}  "
                  f"resid_std={n['resid_std_sparse']:.3e}",
                  flush=True)
            print(f"         {eq}", flush=True)

    print("\n  cross-domain cosines (smallest extract_v1 null):",
          flush=True)
    cd = cross_domain_compare(results)
    for k, v in cd["cos_6d"].items():
        print(f"    6D  cos({k}) = {v:+.4f}", flush=True)
    for k, v in cd["cos_234"].items():
        print(f"    234 cos({k}) = {v:+.4f}", flush=True)

    out_blob = {"config": cfg, "seed": seed,
                "results": results, "cross_domain": cd,
                "elapsed_total_s": time.time() - t0}
    with open(out, "w") as f:
        json.dump(out_blob, f, indent=2)
    print(f"\n[per_domain_kbk] wrote {out} "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
