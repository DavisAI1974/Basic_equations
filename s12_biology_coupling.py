"""
S12 #3 -- the "how-coupled" test for biology.

#1 (5-seed gate) showed biology's null[0] is 0.906 +/- 0.007 MI-coupling
(MI ~= 0.28*H_a), reproducible. #3 asks the substantive question: is that MI
participation actually the DYNAMICAL COUPLING, or incidental? Test: scale the
Lotka-Volterra interaction terms (prey*pred) by a factor g (g=1 baseline;
g=0 fully decouples prey and predator), re-extract null[0], and watch whether
(a) the MI-coupling FRACTION and (b) the slope in MI ~= slope*H_a track g.

Prediction (Rule C speaking posture): I THINK lowering g toward 0 will
collapse the MI-coupling fraction (decoupled species -> no inter-channel MI ->
biology should fall back toward equal-entropy/degenerate, like physics), and
that raising g will hold or strengthen it. If instead MI persists at g=0, that
is informative too (the 'coupling' reading would be partly shared-driver
artifact, not interaction). We wait on the data and where it points; no verdict
in advance, no pre-assigned meaning.

Reuses build_ensemble_operator_matrix + analyze_M from per_domain_kbk (same
extraction as INFO-023). Basis [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI].
Outputs s12_biology_coupling_results.json.
"""
import sys, json, time
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix, analyze_M

CANARY = "--canary" in sys.argv
G_GRID = [1.0, 0.0] if CANARY else [0.0, 0.25, 0.5, 1.0, 1.5, 2.0]
SEEDS = [11] if CANARY else [11, 22, 33]
N_ENS, T, DT = 600, 30.0, 0.02

d_mi = np.array([0, 0, 0, 0, 0, 1.0])


def simulate_biology_coupled(seed, g, N_ens=600, T=30.0, dt=0.02):
    """Lotka-Volterra with interaction (prey*pred) terms scaled by g.
    g=1 reproduces simulate_biology exactly; g=0 fully decouples the species."""
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    prey = 1.0 + 0.1 * rng.standard_normal(N_ens)
    pred = 0.5 + 0.1 * rng.standard_normal(N_ens)
    X1 = np.zeros((steps, N_ens)); X2 = np.zeros((steps, N_ens))
    for t in range(steps):
        dprey = (1.0 * prey - g * 0.5 * prey * pred) * dt
        dpred = (g * 0.3 * prey * pred - 0.4 * pred) * dt
        prey = np.clip(prey + dprey + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt), 0.01, None)
        pred = np.clip(pred + dpred + 0.01 * rng.standard_normal(N_ens) * np.sqrt(dt), 0.01, None)
        X1[t] = prey; X2[t] = pred
    return X1, X2


def readout(v):
    v = np.array(v) / np.linalg.norm(v)
    mi_frac = float((v @ d_mi) ** 2)              # MI-coupling fraction of null
    # slope in MI ~= slope*H_a from the relation v0*H_a + v5*MI (+...) ~ 0
    slope = float(-v[0] / v[5]) if abs(v[5]) > 0.1 else None
    return mi_frac, slope


t0 = time.time()
print("=" * 92)
print(f"S12 #3 biology coupling sweep {'[CANARY]' if CANARY else ''}: "
      f"g in {G_GRID}, seeds {SEEDS}")
print("=" * 92)
results = {}
for g in G_GRID:
    per_seed = []
    for s in SEEDS:
        X1, X2 = simulate_biology_coupled(s, g, N_ENS, T, DT)
        M, _ = build_ensemble_operator_matrix(X1, X2)
        r = analyze_M(M, f"biology_g{g}")
        v = r["extract_v1_v_null_6d"]
        mi_frac, slope = readout(v)
        mi_mean = r["operator_means"][5]          # mean MI level (absolute coupling)
        per_seed.append(dict(seed=s, mi_frac=mi_frac, slope=slope,
                             mi_mean=float(mi_mean),
                             null_rank=r["kbk_rank_gap"]["rank_null_subspace_kbk"]))
    mf = np.array([p["mi_frac"] for p in per_seed])
    sl = np.array([p["slope"] for p in per_seed if p["slope"] is not None])
    mm = np.array([p["mi_mean"] for p in per_seed])
    results[f"g={g}"] = dict(per_seed=per_seed,
                             mi_frac_mean=float(mf.mean()), mi_frac_std=float(mf.std()),
                             slope_mean=(float(sl.mean()) if len(sl) else None),
                             slope_std=(float(sl.std()) if len(sl) else None),
                             mi_mean_mean=float(mm.mean()))
    sl_s = f"{sl.mean():+.3f}+/-{sl.std():.3f}" if len(sl) else "n/a (MI absent)"
    print(f"  g={g:4.2f}  MI-frac={mf.mean():.3f}+/-{mf.std():.3f}  "
          f"slope(MI~s*H_a)={sl_s}  meanMI={mm.mean():.3f}")

json.dump(dict(meta=dict(canary=CANARY, N_ens=N_ENS, T=T, dt=DT, seeds=SEEDS,
                         elapsed_s=round(time.time() - t0, 1)),
               g_grid=G_GRID, results=results),
          open("s12_biology_coupling_canary.json" if CANARY
               else "s12_biology_coupling_results.json", "w"), indent=2)
print(f"\nWrote results ({time.time()-t0:.0f}s)")
