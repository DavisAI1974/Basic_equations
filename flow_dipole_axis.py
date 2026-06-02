"""
MEASURE the flow dipole's axis -- the object we never placed in operator space.

Greg (S15-eve): "we never measured what axis the 'flow' dipole is on!" Correct.
INFO-039 only structurally IDENTIFIED the paper's flow form as our operator
family; INFO-040 decomposed the LEVEL null; today's time_flow_probe put rates in
as passive null members. None of these measured the flow dipole itself.

The flow dipole (info-dipole paper, davisai.ai/dipole):
    d(MI)/dt ~ c_a*H_a + c_b*H_b + c_self_a*H_a^2 + c_self_b*H_b^2 + c_cross*H_a*H_b + const
with the OPPOSITION SIGNATURE: c_self and c_cross of opposing sign.

Here dMI/dt is the regression TARGET (left-hand side), not a passive column. We
regress it on the self+cross+linear basis along the operator trajectory and read:
  (1) R^2 -- does the flow form fit at all (is there a flow dipole)?
  (2) opposition sign -- do the self-terms (H_a^2,H_b^2) oppose the cross-term
      (H_a*H_b), as the paper claims?
  (3) AXIS -- the normalized quadratic coefficient vector (c_self_a, c_self_b,
      c_cross). cos to the equal-entropy substrate (-1,-1,+2)/sqrt6, and cos to
      the LEVEL-null quadratic axis we extract the usual way. Is the flow axis
      the substrate, the level-null, or a NEW direction we've never seen?

Across 4 simulated domains (3 seeds) + real gravity (GW150914 LIGO).

Caveats kept load-bearing (Result Discipline):
  - dMI/dt is a finite-difference rate on autocorrelated time series; R^2 and
    significance are soft (INFO-026). We lean on the AXIS + sign, not R^2 alone.
  - Per-domain regression is OLS with intercept on raw operators (no standardize)
    so the coefficient vector is the flow dipole in native operator units; we
    report the quadratic-subspace DIRECTION (scale-free) for axis comparison.
  - Gravity = 1 clean event (alignment) -> ISOLATED, not located.
No pre-assigned meaning; we read what the axis actually is.
"""
import json
import time

import numpy as np

from per_domain_kbk import (
    simulate_physics, simulate_biology, simulate_chemistry, simulate_geology,
    build_ensemble_operator_matrix,
)
from kbk_pipeline import extract_v1, compute_operator_matrix
from grav_time_retaining import load_strain, preprocess, FS, EDGE_S, WIN_S, STRIDE_S

SUBSTRATE_Q = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)   # equal-entropy quad axis
SIM_DOMAINS = {
    "physics_duffing": simulate_physics,
    "biology_lotka": simulate_biology,
    "chemistry_brussel": simulate_chemistry,
    "geology_bk": simulate_geology,
}


def cos(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else float("nan")


def measure_flow_axis(M, dt_row):
    """M columns = [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]. Regress dMI/dt on
    [H_a,H_b,H_a^2,H_b^2,H_a*H_b]+const along the row (time) trajectory."""
    Ha, Hb, Ha2, Hb2, HaHb, MI = (M[:, i] for i in range(6))
    dMI = np.gradient(MI) / dt_row                       # flow target
    X = np.column_stack([Ha, Hb, Ha2, Hb2, HaHb, np.ones_like(Ha)])
    beta, *_ = np.linalg.lstsq(X, dMI, rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((dMI - pred) ** 2))
    ss_tot = float(np.sum((dMI - dMI.mean()) ** 2)) + 1e-18
    r2 = 1.0 - ss_res / ss_tot
    c_a, c_b, c_sa, c_sb, c_cross, _ = beta

    quad = np.array([c_sa, c_sb, c_cross])               # flow dipole quad axis
    # opposition: self-terms vs cross-term sign (paper: opposite)
    self_sign = np.sign(c_sa + c_sb)
    opposition = bool(self_sign != 0 and np.sign(c_cross) == -self_sign)

    # level-null (the usual measured axis) quadratic part for comparison
    v_lvl, _, _, _ = extract_v1(M)
    lvl_quad = v_lvl[2:5]

    return {
        "flow_R2": round(r2, 3),
        "flow_coefs": {"H_a": round(float(c_a), 4), "H_b": round(float(c_b), 4),
                       "H_a^2": round(float(c_sa), 4), "H_b^2": round(float(c_sb), 4),
                       "H_a*H_b": round(float(c_cross), 4)},
        "opposition_self_vs_cross": opposition,
        "flow_quad_axis": [round(float(x), 4) for x in quad / (np.linalg.norm(quad) + 1e-18)],
        "cos_flow_to_substrate": round(abs(cos(quad, SUBSTRATE_Q)), 3),
        "cos_flow_to_level_null_quad": round(abs(cos(quad, lvl_quad)), 3),
        "level_null_quad_axis": [round(float(x), 4) for x in
                                 lvl_quad / (np.linalg.norm(lvl_quad) + 1e-18)],
    }


def run_sim(seeds=(11, 22, 33)):
    cfg = dict(N_ens=400, T=20.0, dt=0.02)
    res = {}
    for name, simfn in SIM_DOMAINS.items():
        per_seed = []
        for s in seeds:
            X1, X2 = simfn(seed=s, **cfg)
            M, _ = build_ensemble_operator_matrix(X1, X2)
            r = measure_flow_axis(M, dt_row=cfg["dt"])
            r["seed"] = s
            per_seed.append(r)
        # cross-seed axis stability
        axes = np.array([p["flow_quad_axis"] for p in per_seed])
        cseed = [abs(cos(axes[0], axes[i])) for i in range(1, len(axes))]
        res[name] = {
            "per_seed": per_seed,
            "flow_R2_mean": round(float(np.mean([p["flow_R2"] for p in per_seed])), 3),
            "cos_flow_to_substrate_mean": round(
                float(np.mean([p["cos_flow_to_substrate"] for p in per_seed])), 3),
            "cos_flow_to_level_null_mean": round(
                float(np.mean([p["cos_flow_to_level_null_quad"] for p in per_seed])), 3),
            "opposition_fraction": round(
                float(np.mean([p["opposition_self_vs_cross"] for p in per_seed])), 2),
            "cross_seed_axis_cos": [round(c, 3) for c in cseed],
        }
    return res, cfg


def run_gravity():
    h1 = preprocess(load_strain("data/ligo/H-H1_GW150914_32s.hdf5"))
    l1 = preprocess(load_strain("data/ligo/L-L1_GW150914_32s.hdf5"))
    l1 = -np.roll(l1, int(round(0.0069 * FS)))
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1[crop:-crop],
                                int(WIN_S * FS), int(STRIDE_S * FS), mi_bins=16)
    return measure_flow_axis(M, dt_row=STRIDE_S)


def main():
    t0 = time.time()
    print("[flow_dipole_axis] measuring the flow dipole d(MI)/dt axis\n", flush=True)
    sim, cfg = run_sim()
    print(f"SIM domains (cfg={cfg}, 3 seeds):", flush=True)
    for name, r in sim.items():
        print(f"\n  {name}:", flush=True)
        print(f"    flow R^2 mean       = {r['flow_R2_mean']}", flush=True)
        print(f"    opposition fraction = {r['opposition_fraction']}  "
              f"(paper: self vs cross opposite sign)", flush=True)
        print(f"    cos(flow -> equal-entropy substrate) = "
              f"{r['cos_flow_to_substrate_mean']}", flush=True)
        print(f"    cos(flow -> our LEVEL-null axis)     = "
              f"{r['cos_flow_to_level_null_mean']}", flush=True)
        print(f"    cross-seed flow-axis cos = {r['cross_seed_axis_cos']}", flush=True)
        print(f"    seed11 flow quad axis = {r['per_seed'][0]['flow_quad_axis']} "
              f"(coefs {r['per_seed'][0]['flow_coefs']})", flush=True)

    grav = run_gravity()
    print(f"\nGRAVITY (GW150914, real LIGO, 1 clean event -> ISOLATED):", flush=True)
    print(f"    flow R^2            = {grav['flow_R2']}", flush=True)
    print(f"    opposition          = {grav['opposition_self_vs_cross']}", flush=True)
    print(f"    flow quad axis      = {grav['flow_quad_axis']}", flush=True)
    print(f"    cos(flow -> substrate)        = {grav['cos_flow_to_substrate']}", flush=True)
    print(f"    cos(flow -> level-null quad)  = {grav['cos_flow_to_level_null_quad']}", flush=True)
    print(f"    flow coefs          = {grav['flow_coefs']}", flush=True)

    json.dump({"sim_config": cfg, "sim": sim, "gravity_GW150914": grav,
               "substrate_quad_axis": SUBSTRATE_Q.tolist(),
               "elapsed_s": round(time.time() - t0, 1)},
              open("flow_dipole_axis_results.json", "w"), indent=2)
    print(f"\n[flow_dipole_axis] wrote flow_dipole_axis_results.json "
          f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
