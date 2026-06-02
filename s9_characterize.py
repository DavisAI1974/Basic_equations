"""
S9 characterization: map the substrate-flip threshold across all four
force caricatures, and test whether adding a universal-energy (E_total)
coupling term LOWERS the threshold (i.e. whether energy-mediated
coupling PROMOTES the flip).

Baseline asymmetry sweeps use the ORIGINAL verified simulators
(simulate_em_config, simulate_strong_config, simulate_gravity_config).

Two minimal, clearly-labelled experimental variants (needed for the
question; force laws otherwise byte-identical to the originals):
  - em_plus_etotal / strong_plus_etotal: original force law + gravity's
    E_total term  Fadd = G_add * E_total * dx / (dx^2 + eps_soft).
  - weak_asym: original weak force law + an omega2 knob (original
    hardcodes omega=1.0 for both channels, so it cannot be swept).

Extraction is the ORIGINAL raw-covariance KBK (build_ensemble_operator_
matrix + extract_v1). 3 seeds; mean +/- std reported (Result Discipline:
no operator-space coordinate without >=3 seeds and inter-seed scatter).
"""
import sys, json, time
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix
from kbk_pipeline import extract_v1
from em_strong_glance import simulate_em_config, simulate_strong_config
from gravity_glance import simulate_gravity_config

ATTR3 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22, 33)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)
EPS_SOFT = 0.50

def _etotal(x1, v1, x2, v2, G_add):
    E = x1*x1 + v1*v1 + x2*x2 + v2*v2
    dx = x1 - x2
    return G_add * E * dx / (dx*dx + EPS_SOFT)

def em_plus_etotal(seed, N_ens, T, dt, omega1=1.0, omega2=1.2, K=0.30, G_add=0.0):
    """Original EM force law + optional E_total term."""
    rng = np.random.default_rng(seed); steps = int(T/dt)
    gamma_rad = 0.02; sigma_noise = 0.02
    x1=0.5*rng.standard_normal(N_ens); v1=0.1*rng.standard_normal(N_ens)
    x2=0.5*rng.standard_normal(N_ens); v2=0.1*rng.standard_normal(N_ens)
    X1=np.zeros((steps,N_ens)); X2=np.zeros((steps,N_ens))
    for t in range(steps):
        Fa = _etotal(x1,v1,x2,v2,G_add)
        F1 = -(omega1**2)*x1 + K*x2 - Fa
        F2 = -(omega2**2)*x2 + K*x1 + Fa
        v1=v1+F1*dt-gamma_rad*v1*dt; v2=v2+F2*dt-gamma_rad*v2*dt
        x1=x1+v1*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        x2=x2+v2*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        X1[t]=x1; X2[t]=x2
    return X1,X2

def strong_plus_etotal(seed, N_ens, T, dt, omega1=0.5, omega2=0.5,
                       K_lin=0.20, K_conf=0.80, G_add=0.0):
    """Original strong force law + optional E_total term."""
    rng = np.random.default_rng(seed); steps = int(T/dt)
    gamma = 0.05; sigma_noise = 0.02
    x1=0.5*rng.standard_normal(N_ens); v1=0.1*rng.standard_normal(N_ens)
    x2=0.5*rng.standard_normal(N_ens); v2=0.1*rng.standard_normal(N_ens)
    X1=np.zeros((steps,N_ens)); X2=np.zeros((steps,N_ens))
    for t in range(steps):
        dx=x1-x2; Fconf=K_lin*dx+K_conf*dx**3
        Fa = _etotal(x1,v1,x2,v2,G_add)
        F1 = -(omega1**2)*x1 - Fconf - Fa
        F2 = -(omega2**2)*x2 + Fconf + Fa
        v1=v1+F1*dt-gamma*v1*dt; v2=v2+F2*dt-gamma*v2*dt
        x1=x1+v1*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        x2=x2+v2*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        X1[t]=x1; X2[t]=x2
    return X1,X2

def weak_asym(seed, N_ens, T, dt, omega1=1.0, omega2=1.0, K=0.50, M_med=5.0):
    """Original weak force law + an omega2 knob (original hardcodes 1.0)."""
    rng = np.random.default_rng(seed); steps = int(T/dt)
    gamma = 0.30; sigma_noise = 0.03
    x1=0.5*rng.standard_normal(N_ens); v1=0.1*rng.standard_normal(N_ens)
    x2=0.5*rng.standard_normal(N_ens); v2=0.1*rng.standard_normal(N_ens)
    X1=np.zeros((steps,N_ens)); X2=np.zeros((steps,N_ens))
    for t in range(steps):
        dx=x1-x2; suppress=np.exp(-M_med*dx*dx)
        F1 = -(omega1**2)*x1 + K*x2*suppress
        F2 = -(omega2**2)*x2 + K*x1*suppress
        v1=v1+F1*dt-gamma*v1*dt; v2=v2+F2*dt-gamma*v2*dt
        x1=x1+v1*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        x2=x2+v2*dt+sigma_noise*rng.standard_normal(N_ens)*np.sqrt(dt)
        X1[t]=x1; X2[t]=x2
    return X1,X2

def sig(X1, X2):
    M,_ = build_ensemble_operator_matrix(X1, X2)
    v,S,_,_ = extract_v1(M)
    MI = M[:,5]; sub = v[2:5]; sub = sub/(np.linalg.norm(sub)+1e-30)
    return abs(v[5]), float(MI.std()), float(abs(sub@ATTR3))

def run(name, simfn, **skw):
    mc=[]; ms=[]; ac=[]
    for s in seeds:
        a,b,c = sig(*simfn(seed=s, **cfg, **skw)); mc.append(a); ms.append(b); ac.append(c)
    mc=np.array(mc)
    flip = "FLIP" if mc.mean()>0.5 else "    "
    print(f"{name:30s} |MI|={mc.mean():.3f}+/-{mc.std():.3f}  MIstd={np.mean(ms):.4f}  attr={np.mean(ac):.3f}  {flip}")
    return dict(name=name, mi_coef_mean=float(mc.mean()), mi_coef_std=float(mc.std()),
                mi_coef_seeds=[float(x) for x in mc], mi_std=float(np.mean(ms)),
                attractor_cos=float(np.mean(ac)))

t0=time.time(); R={}
print("="*92); print(f"S9 CHARACTERIZE {'[CANARY]' if CANARY else '[FULL]'} seeds={seeds}"); print("="*92)

print("\n--- weak asymmetry sweep (new omega2 knob) ---")
for o2 in [1.0, 1.3, 1.7, 2.0, 2.5, 3.0]:
    R[f"weak_1.0/{o2}"] = run(f"weak 1.0/{o2}", weak_asym, omega1=1.0, omega2=o2)

print("\n--- STRONG + E_total: does energy-coupling lower the (infinite) threshold? ---")
for G in [0.0, 0.05, 0.15, 0.40]:
    print(f"  G_add={G}:")
    for o2 in [0.7, 1.0, 1.5, 2.0]:
        R[f"strong_0.5/{o2}_G{G}"] = run(f"  strong 0.5/{o2} G={G}", strong_plus_etotal, omega1=0.5, omega2=o2, G_add=G)

print("\n--- EM + E_total: does it shift EM's ~1.3x threshold? ---")
for G in [0.0, 0.15]:
    print(f"  G_add={G}:")
    for o2 in [1.1, 1.2, 1.3, 1.4]:
        R[f"em_1.0/{o2}_G{G}"] = run(f"  EM 1.0/{o2} G={G}", em_plus_etotal, omega1=1.0, omega2=o2, G_add=G)

out=dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg, eps_soft=EPS_SOFT,
                   elapsed_s=round(time.time()-t0,1)), results=R)
fn="s9_characterize_canary.json" if CANARY else "s9_characterize_results.json"
json.dump(out, open(fn,"w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
