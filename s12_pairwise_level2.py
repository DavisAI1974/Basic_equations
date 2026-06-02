"""
S12 pairwise Level-2 probe -- "couple two sciences at a time, look for the
coupling dipole between PAIRS" (the right shape after INFO-009 showed no
UNIVERSAL-across-4 Level-2 dipole, and INFO-040 showed coupling is per-domain).

For each of the 6 ordered science pairs (X,Y) we run BOTH dynamical systems
(the exact per_domain_kbk dynamics) and add a scale-free diffusive coupling g
between their PRIMARY channels: each primary gets +g*(z_other - z_own), where
z = per-step ensemble z-score (dimensionless, O(1), so g is comparable across
heterogeneous systems). Channel A = X's primary, channel B = Y's primary. We
extract the 6-op null[0] (basis [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI]) with the SAME
machinery as INFO-023/040 and read:
  - MI-coupling fraction of null[0]  (does coupling two DIFFERENT sciences
    create the MI-participation signature, like biology's intra-domain coupling)
  - the coupled-null DIRECTION       (is it the SAME for all 6 pairs => a
    universal cross-science coupling dipole, or pair-specific => heterogeneous,
    consistent with INFO-040/009)

g sweep {0 (uncoupled control), 0.5}; seeds {11,22,33}.

CAVEAT (frame can't grade itself): the diffusive coupling is a TOY, methods-
level coupling (same family as level2's diffusion/migration), chosen for
neutrality. This probes whether the EXTRACTION detects and distinguishes
cross-science coupling, NOT a claim that physics and biology physically couple.
Speaking posture: I THINK coupling will create MI-participation (any coupling
creates MI) and that the coupled-null DIRECTIONS will be pair-specific (no
universal Level-2 dipole), but we wait on the data.

Output s12_pairwise_level2_results.json.
"""
import sys, json, time, itertools
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix, analyze_M

CANARY = "--canary" in sys.argv
G_LIST = [0.0, 0.5]
SEEDS = [11] if CANARY else [11, 22, 33]
N_ENS = 200 if CANARY else 600
T, DT = 30.0, 0.02
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def z(p):
    return (p - p.mean()) / (p.std() + 1e-9)


# Each science: init(rng,N)->state; step(state,c,rng,t,N)-> (state, primary).
# c is the external coupling term added to the primary's update. Dynamics match
# per_domain_kbk exactly when c=0.
def phys_init(rng, N):
    return dict(x1=0.5*rng.standard_normal(N), v1=0.1*rng.standard_normal(N),
                x2=0.5*rng.standard_normal(N), v2=0.1*rng.standard_normal(N))
def phys_step(s, c, rng, t, N):
    K = 0.2
    F1 = -s['x1'] - 0.5*s['x1']**3 + K*(s['x2']-s['x1']) + 0.1*np.sin(0.5*t*DT) + c
    F2 = -1.2*s['x2'] - 0.4*s['x2']**3 + K*(s['x1']-s['x2']) + 0.08*np.cos(0.4*t*DT)
    s['v1'] += F1*DT - 0.05*s['v1']*DT; s['v2'] += F2*DT - 0.05*s['v2']*DT
    s['x1'] += s['v1']*DT + 0.02*rng.standard_normal(N)*np.sqrt(DT)
    s['x2'] += s['v2']*DT + 0.02*rng.standard_normal(N)*np.sqrt(DT)
    return s, s['x1']

def bio_init(rng, N):
    return dict(prey=1.0+0.1*rng.standard_normal(N), pred=0.5+0.1*rng.standard_normal(N))
def bio_step(s, c, rng, t, N):
    dprey = (1.0*s['prey'] - 0.5*s['prey']*s['pred'])*DT + c*DT
    dpred = (0.3*s['prey']*s['pred'] - 0.4*s['pred'])*DT
    s['prey'] = np.clip(s['prey'] + dprey + 0.01*rng.standard_normal(N)*np.sqrt(DT), 0.01, None)
    s['pred'] = np.clip(s['pred'] + dpred + 0.01*rng.standard_normal(N)*np.sqrt(DT), 0.01, None)
    return s, s['prey']

def chem_init(rng, N):
    return dict(X=1.0+0.1*rng.standard_normal(N), Y=3.0+0.1*rng.standard_normal(N))
def chem_step(s, c, rng, t, N):
    A_, B_ = 1.0, 3.0
    dX = (A_ - (B_+1)*s['X'] + s['X']**2*s['Y'])*DT + c*DT
    dY = (B_*s['X'] - s['X']**2*s['Y'])*DT
    s['X'] = s['X'] + dX + 0.01*rng.standard_normal(N)*np.sqrt(DT)
    s['Y'] = s['Y'] + dY + 0.01*rng.standard_normal(N)*np.sqrt(DT)
    return s, s['X']

def geo_init(rng, N):
    return dict(u1=0.1*rng.standard_normal(N), v1=0.05*rng.standard_normal(N),
                u2=0.1*rng.standard_normal(N), v2=0.05*rng.standard_normal(N))
def geo_step(s, c, rng, t, N):
    K = 0.5
    fr1 = -0.6*s['v1']/(1.0+8.0*s['v1']**2); fr2 = -0.6*s['v2']/(1.0+8.0*s['v2']**2)
    F1 = K*(s['u2']-s['u1']) - 0.6*(s['u1']-0.04*t*DT) + fr1 + c
    F2 = K*(s['u1']-s['u2']) - 0.6*(s['u2']-0.04*t*DT) + fr2
    s['v1'] += F1*DT; s['v2'] += F2*DT
    s['u1'] += s['v1']*DT + 0.015*rng.standard_normal(N)*np.sqrt(DT)
    s['u2'] += s['v2']*DT + 0.015*rng.standard_normal(N)*np.sqrt(DT)
    return s, s['u1']

SCI = {"physics": (phys_init, phys_step), "biology": (bio_init, bio_step),
       "chemistry": (chem_init, chem_step), "geology": (geo_init, geo_step)}
ORDER = ["physics", "biology", "chemistry", "geology"]
PAIRS = list(itertools.combinations(ORDER, 2)) if not CANARY else [("physics", "biology")]


def run_pair(nameX, nameY, g, seed):
    rng = np.random.default_rng(seed)
    iX, sX = SCI[nameX]; iY, sY = SCI[nameY]
    stX = iX(rng, N_ENS); stY = iY(rng, N_ENS)
    steps = int(T/DT)
    A = np.zeros((steps, N_ENS)); B = np.zeros((steps, N_ENS))
    pX = stX[list(stX)[0]]; pY = stY[list(stY)[0]]  # current primaries
    for t in range(steps):
        cX = g*(z(pY) - z(pX)); cY = g*(z(pX) - z(pY))
        stX, pX = sX(stX, cX, rng, t, N_ENS)
        stY, pY = sY(stY, cY, rng, t, N_ENS)
        A[t] = pX; B[t] = pY
    M, _ = build_ensemble_operator_matrix(A, B)
    r = analyze_M(M, f"{nameX}-{nameY}")
    v = np.array(r["extract_v1_v_null_6d"]); v = v/np.linalg.norm(v)
    return v, float((v @ d_mi)**2), float(r["operator_means"][5])


t0 = time.time()
print("="*92)
print(f"S12 pairwise Level-2 {'[CANARY]' if CANARY else ''}: cross-science coupling, "
      f"g={G_LIST}, seeds={SEEDS}")
print("="*92)
results = {}
for X, Y in PAIRS:
    pr = {}
    for g in G_LIST:
        vs, mfs, mis = [], [], []
        for s in SEEDS:
            v, mf, mi = run_pair(X, Y, g, s)
            vs.append(v); mfs.append(mf); mis.append(mi)
        vs = np.array(vs)
        # cross-seed direction stability
        cs = [abs(vs[i]@vs[j]) for i in range(len(vs)) for j in range(i+1, len(vs))]
        pr[f"g={g}"] = dict(mi_frac_mean=float(np.mean(mfs)), mi_frac_std=float(np.std(mfs)),
                            mean_MI=float(np.mean(mis)),
                            dir_stability=float(np.mean(cs)) if cs else None,
                            null_mean=[round(float(x), 3) for x in vs.mean(0)])
    results[f"{X}-{Y}"] = pr
    print(f"  {X[:4]}-{Y[:4]:9s} g0 MIfrac={pr['g=0.0']['mi_frac_mean']:.3f} "
          f"-> g0.5 MIfrac={pr['g=0.5']['mi_frac_mean']:.3f}  "
          f"meanMI {pr['g=0.0']['mean_MI']:.3f}->{pr['g=0.5']['mean_MI']:.3f}  "
          f"dirStab={pr['g=0.5']['dir_stability']}")

# Universal vs pair-specific: pairwise cos among the 6 coupled (g=0.5) nulls.
if not CANARY:
    keys = list(results)
    nv = {k: np.array(results[k]["g=0.5"]["null_mean"]) for k in keys}
    print("\n  cross-PAIR coupled-null |cos| (universal Level-2 dipole? => all ~1):")
    cross = {}
    for a, b in itertools.combinations(keys, 2):
        c = abs(nv[a]@nv[b]/(np.linalg.norm(nv[a])*np.linalg.norm(nv[b])+1e-9))
        cross[f"{a}__{b}"] = round(float(c), 3)
    vals = list(cross.values())
    print(f"    mean|cos|={np.mean(vals):.3f}  min={np.min(vals):.3f}  max={np.max(vals):.3f}")
    results["_cross_pair_cos"] = cross
    results["_cross_pair_summary"] = dict(mean=float(np.mean(vals)),
                                          min=float(np.min(vals)), max=float(np.max(vals)))

json.dump(dict(meta=dict(canary=CANARY, N_ens=N_ENS, T=T, dt=DT, seeds=SEEDS,
                         g_list=G_LIST, elapsed_s=round(time.time()-t0, 1)),
               results=results),
          open("s12_pairwise_level2_canary.json" if CANARY
               else "s12_pairwise_level2_results.json", "w"), indent=2)
print(f"\nWrote results ({time.time()-t0:.0f}s)")
