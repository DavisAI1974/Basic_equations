"""
S13 #2 -- the "how-coupled" knob test for CHEMISTRY's residual.

The S12 5-seed decomposition (s12_coupling_decomposition.py / INFO-040) found
chemistry's null[0] is 0.83 equal-entropy + a STABLE 0.167 residual whose
direction is +0.54*H_a +0.54*H_b +0.32*H_a^2 -0.55*H_b^2 ~ 0 (cross-seed
|cos|=0.9996). That is NOT the MI axis (MI-frac = 0): chemistry's coupling is
"partial / residual", a different kind of object from biology's MI-coupling.

This is the chemistry analogue of s12_biology_coupling.py (which showed biology's
MI-slope tracks the Lotka-Volterra interaction g). Question (kickoff #2): does
chemistry's residual track the Brusselator B parameter the way biology's MI-slope
tracks g? Sweep B (baseline 3.0) across the Hopf bifurcation (B_crit = 1+A^2 = 2,
A=1; B<2 stable fixed point, B>2 limit cycle), re-extract null[0], and watch
whether (a) the residual FRACTION and (b) the residual DIRECTION (its |cos| to
the B=3 baseline relation) move systematically with B.

Speaking posture (Rule C): I THINK the residual may track B, since B is the
Brusselator's bifurcation control and INFO-035 already flagged chemistry's
substrate as B-sensitive near the Hopf threshold. But I want to see what the
data says and where it points -- no verdict in advance, and B=2 sits EXACTLY at
the Hopf threshold so it gets inspected on its own terms (no tent-widening).

Reuses build_ensemble_operator_matrix + analyze_M from per_domain_kbk (same
extraction as INFO-023/INFO-040). Basis [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI].
Outputs s13_chemistry_residual_results.json.
"""
import sys, json, time
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix, analyze_M

CANARY = "--canary" in sys.argv
B_GRID = [3.0, 2.0] if CANARY else [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
SEEDS = [11] if CANARY else [11, 22, 33]
N_ENS, T, DT = 600, 30.0, 0.02
A_ = 1.0

# decomposition axes (same as s12_coupling_decomposition.py)
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
BASIS = [d_lin, d_quad, d_mi]
# B=3 baseline residual relation from INFO-040
res_base = np.array([0.54, 0.54, 0.32, -0.55, 0.0, 0.0])
res_base = res_base / np.linalg.norm(res_base)


def simulate_chemistry_B(seed, B, N_ens=600, T=30.0, dt=0.02):
    """Brusselator with the B control parameter as the knob (A fixed at 1).
    B=3.0 reproduces simulate_chemistry exactly; B_crit = 1+A^2 = 2."""
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    X = 1.0 + 0.1 * rng.standard_normal(N_ens)
    Y = B / A_ + 0.1 * rng.standard_normal(N_ens)   # fixed point Y*=B/A
    X1 = np.zeros((steps, N_ens)); X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dX = (A_ - (B + 1) * X + X ** 2 * Y) * dt
        dY = (B * X - X ** 2 * Y) * dt
        X = X + dX + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        Y = Y + dY + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt)
        X1[t] = X; X2[t] = Y
    return X1, X2


def decompose(v):
    """Project null[0] onto eq-entropy / MI / residual; return fractions,
    residual direction, and |cos| of residual dir to the B=3 baseline relation."""
    v = np.array(v) / np.linalg.norm(v)
    a = (v @ d_lin) ** 2; b = (v @ d_quad) ** 2; c = (v @ d_mi) ** 2
    eq_frac = a + b; mi_frac = c; res_frac = max(0.0, 1 - a - b - c)
    rv = v.copy()
    for bb in BASIS:
        rv = rv - (v @ bb) * bb
    n = np.linalg.norm(rv)
    if n > 1e-6:
        rdir = rv / n
        res_cos_base = float(abs(rdir @ res_base))
    else:
        rdir = np.zeros(6); res_cos_base = None
    return eq_frac, mi_frac, res_frac, rdir, res_cos_base


t0 = time.time()
print("=" * 96)
print(f"S13 #2 chemistry residual knob test {'[CANARY]' if CANARY else ''}: "
      f"B in {B_GRID}, seeds {SEEDS}  (Hopf B_crit=2)")
print("=" * 96)
print(f"{'B':>5s} {'eqEnt':>14s} {'MI':>12s} {'residual':>14s} "
      f"{'res|cos|->B3base':>16s} {'res-dir xseed':>13s}")
results = {}
for B in B_GRID:
    eqs, mis, ress, coss, rdirs = [], [], [], [], []
    per_seed = []
    for s in SEEDS:
        X1, X2 = simulate_chemistry_B(s, B, N_ENS, T, DT)
        M, _ = build_ensemble_operator_matrix(X1, X2)
        r = analyze_M(M, f"chem_B{B}")
        v = r["extract_v1_v_null_6d"]
        eq, mi, res, rdir, rcos = decompose(v)
        eqs.append(eq); mis.append(mi); ress.append(res); rdirs.append(rdir)
        if rcos is not None:
            coss.append(rcos)
        per_seed.append(dict(seed=s, eq_frac=eq, mi_frac=mi, res_frac=res,
                             res_cos_to_B3base=rcos,
                             null_rank=r["kbk_rank_gap"]["rank_null_subspace_kbk"],
                             mean_MI=float(r["operator_means"][5])))
    eqs, mis, ress = map(np.array, (eqs, mis, ress))
    # cross-seed residual-direction stability (sign-aligned)
    rd = np.array(rdirs); ref = rd[0]
    rd_al = np.array([x * np.sign(x @ ref) for x in rd])
    rcs = [abs(rd[i] @ rd[j]) for i in range(len(rd)) for j in range(i + 1, len(rd))]
    res_xseed = float(np.mean(rcs)) if rcs else 1.0
    res_mean_dir = rd_al.mean(0)
    res_mean_dir = res_mean_dir / (np.linalg.norm(res_mean_dir) + 1e-12)
    OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]
    res_rel = " ".join(f"{res_mean_dir[i]:+.2f}*{OPS[i]}"
                       for i in range(6) if abs(res_mean_dir[i]) > 0.15) + " ~ 0"
    cmean = float(np.mean(coss)) if coss else None
    results[f"B={B}"] = dict(
        per_seed=per_seed,
        eq_frac=[float(eqs.mean()), float(eqs.std())],
        mi_frac=[float(mis.mean()), float(mis.std())],
        res_frac=[float(ress.mean()), float(ress.std())],
        res_cos_to_B3base_mean=cmean,
        res_dir_xseed_abscos=res_xseed,
        res_relation=res_rel,
    )
    cstr = f"{cmean:.3f}" if cmean is not None else " n/a"
    print(f"{B:5.1f} {eqs.mean():6.3f}+/-{eqs.std():4.3f} "
          f"{mis.mean():5.3f}+/-{mis.std():4.3f} "
          f"{ress.mean():6.3f}+/-{ress.std():4.3f} "
          f"{cstr:>16s} {res_xseed:13.3f}")

print("\nresidual relation per B (does the direction rotate with B?):")
for B in B_GRID:
    print(f"  B={B:4.1f}  {results[f'B={B}']['res_relation']}")

out = dict(meta=dict(canary=CANARY, N_ens=N_ENS, T=T, dt=DT, seeds=SEEDS, A=A_,
                     B_crit_hopf=1 + A_ ** 2, baseline_B=3.0,
                     elapsed_s=round(time.time() - t0, 1)),
           B_grid=B_GRID, baseline_residual_relation_B3=res_base.tolist(),
           results=results)
out["reading"] = (
    "Chemistry residual knob test: sweep Brusselator B across the Hopf "
    "threshold (B_crit=2). Compare to biology, whose MI-slope tracks the "
    "Lotka-Volterra g monotonically (s12_biology_coupling). See res_frac vs B "
    "and res_cos_to_B3base vs B for whether chemistry's residual is a "
    "B-controlled coupling-strength readout or a B-stable structural relation."
)
json.dump(out, open("s13_chemistry_residual_canary.json" if CANARY
                    else "s13_chemistry_residual_results.json", "w"), indent=2)
print(f"\nWrote results ({time.time()-t0:.0f}s)")
