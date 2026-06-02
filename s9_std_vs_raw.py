"""
S9 follow-up A: standardized-vs-raw covariance fork on the flip.

The INFO-033 substrate flip was extracted with RAW centered-covariance
SVD (extract_v1): the null is the smallest-ABSOLUTE-variance operator
combination, so an operator whose variance collapses (MI under
asymmetry) dominates it. Question (INFO-024 fork): how much of the flip
is that procedure choice vs a structural linear dependency?

Two extractions on the SAME operator matrix M (cols [Ha,Hb,Ha^2,Hb^2,
Ha*Hb,MI]):
  RAW  : null of centered M            (constancy detector; original)
  STD  : null of per-column standardized M  (correlation/linear-
         dependency detector -- removes absolute-scale, so a merely
         small-variance MI is NOT favored).

If the MI-dominant flip is a raw-scale effect, STD null should NOT be
MI-dominant (it should stay on the channel-balance direction). If MI is
structurally locked as a linear relation, STD null stays MI-dominant.

Original sims + original raw extract_v1 reused; STD is the only addition.
3 seeds, mean +/- std.
"""
import sys, json, time
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix
from kbk_pipeline import extract_v1
from em_strong_glance import simulate_em_config, simulate_strong_config
from gravity_glance import simulate_gravity_config
from four_force_probe import simulate_weak

ATTR3 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22, 33)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)

def extract_std(M):
    mu = M.mean(0); sd = M.std(0) + 1e-12
    Z = (M - mu) / sd
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    v = Vt[-1]
    j = int(np.argmax(np.abs(v)))
    if v[j] < 0:
        v = -v
    return v, S

def both(X1, X2):
    M, _ = build_ensemble_operator_matrix(X1, X2)
    vr, _, _, _ = extract_v1(M)            # raw (original)
    vs, _ = extract_std(M)                 # standardized
    subr = vr[2:5] / (np.linalg.norm(vr[2:5]) + 1e-30)
    subs = vs[2:5] / (np.linalg.norm(vs[2:5]) + 1e-30)
    return (abs(vr[5]), abs(vs[5]), float(abs(subr@ATTR3)), float(abs(subs@ATTR3)))

def run(name, simfn, **skw):
    rmi=[]; smi=[]; rat=[]; sat=[]
    for s in seeds:
        a,b,c,d = both(*simfn(seed=s, **cfg, **skw))
        rmi.append(a); smi.append(b); rat.append(c); sat.append(d)
    rmi=np.array(rmi); smi=np.array(smi)
    print(f"{name:24s}  RAW|MI|={rmi.mean():.3f}+/-{rmi.std():.3f} (attr {np.mean(rat):.2f})   "
          f"STD|MI|={smi.mean():.3f}+/-{smi.std():.3f} (attr {np.mean(sat):.2f})")
    return dict(name=name, raw_mi=float(rmi.mean()), raw_mi_std=float(rmi.std()),
                std_mi=float(smi.mean()), std_mi_std=float(smi.std()),
                raw_attr=float(np.mean(rat)), std_attr=float(np.mean(sat)))

t0=time.time(); R={}
print("="*100); print(f"S9 STD-vs-RAW {'[CANARY]' if CANARY else '[FULL]'} seeds={seeds}"); print("="*100)
print("(RAW = original constancy-detector; STD = standardized correlation-detector)\n")

print("--- gravity sweep ---")
for o2 in [0.8, 0.9, 1.0, 1.2]:
    R[f"grav_0.8/{o2}"] = run(f"gravity 0.8/{o2}", simulate_gravity_config, omega1=0.8, omega2=o2, universal=True)
print("--- EM sweep ---")
for o2 in [1.0, 1.3, 1.5, 2.0]:
    R[f"em_1.0/{o2}"] = run(f"EM 1.0/{o2}", simulate_em_config, omega1=1.0, omega2=o2)
print("--- strong sweep ---")
for o2 in [0.5, 1.0, 1.5, 2.0]:
    R[f"strong_0.5/{o2}"] = run(f"strong 0.5/{o2}", simulate_strong_config, omega1=0.5, omega2=o2)
print("--- weak baseline ---")
R["weak_sym"] = run("weak 1.0/1.0", simulate_weak)

out=dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg,
                   elapsed_s=round(time.time()-t0,1)), results=R)
fn="s9_std_vs_raw_canary.json" if CANARY else "s9_std_vs_raw_results.json"
json.dump(out, open(fn,"w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
