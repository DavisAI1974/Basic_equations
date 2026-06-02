"""
S9 double-check (uses ORIGINAL Session 8 code, no reimplementation):

Imports the original simulators + original KBK extraction and asks one
question the Session 8 glances did NOT test: INFO-031 concluded "EM and
strong keep their substrate; only gravity (INFO-032) flips to MI-dominant
under asymmetry." But each glance tested only ONE asymmetry value
(EM 1.0/1.2; strong 0.5/0.7). Is the EM/strong non-flip robust, or does
it just reflect a mild asymmetry? Sweep omega2 and watch MI_coef / MI_std.

Original code used verbatim:
  build_ensemble_operator_matrix, extract_v1  (per_domain_kbk/kbk_pipeline)
  simulate_em_config, simulate_strong_config  (em_strong_glance)
  simulate_gravity_config                     (gravity_glance)
  simulate_weak                               (four_force_probe)

No pre-assigned meaning. Speaking posture: I think pushing asymmetry
MIGHT flip EM/strong too (making the flip a generic asymmetry effect,
not gravity-specific) OR might not (confirming gravity-specificity) --
we wait on the data.
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

def sig(X1, X2):
    M, _ = build_ensemble_operator_matrix(X1, X2)
    v, S, _, _ = extract_v1(M)
    MI = M[:, 5]
    sub = v[2:5]; sub = sub / (np.linalg.norm(sub) + 1e-30)
    return abs(v[5]), float(MI.std()), float(MI.mean()), float(abs(sub @ ATTR3)), v

def run(name, simfn, **skw):
    mc, ms, mm, ac, vv = [], [], [], [], []
    for s in seeds:
        X1, X2 = simfn(seed=s, **cfg, **skw)
        a, b, c, d, v = sig(X1, X2)
        mc.append(a); ms.append(b); mm.append(c); ac.append(d); vv.append(v)
    v = np.mean(vv, axis=0)
    flip = "FLIP" if np.mean(mc) > 0.5 else "    "
    print(f"{name:34s} |MI|={np.mean(mc):.3f}  MIstd={np.mean(ms):.4f}  "
          f"MImean={np.mean(mm):.3f}  attr={np.mean(ac):.3f}  {flip}")
    return dict(name=name, mi_coef=float(np.mean(mc)), mi_std=float(np.mean(ms)),
                mi_mean=float(np.mean(mm)), attractor_cos=float(np.mean(ac)),
                v_null=v.tolist())

t0 = time.time(); R = {}
print("="*86); print(f"S9 DOUBLE-CHECK (original code) {'[CANARY]' if CANARY else '[FULL]'} seeds={seeds}"); print("="*86)

print("\n--- reproduce Session 8 baselines (sanity) ---")
R["grav_sym_0.8/0.8"]  = run("gravity sym 0.8/0.8 (baseline)",  simulate_gravity_config, omega1=0.8, omega2=0.8, universal=True)
R["grav_asym_0.8/1.0"] = run("gravity asym 0.8/1.0 (INFO-032)", simulate_gravity_config, omega1=0.8, omega2=1.0, universal=True)
R["em_asym_1.0/1.2"]   = run("EM asym 1.0/1.2 (INFO-031)",      simulate_em_config, omega1=1.0, omega2=1.2)
R["strong_asym_0.5/0.7"]= run("strong asym 0.5/0.7 (INFO-031)", simulate_strong_config, omega1=0.5, omega2=0.7)

print("\n--- EM asymmetry sweep (does EM flip if pushed?) ---")
for o2 in [1.0, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0]:
    R[f"em_1.0/{o2}"] = run(f"EM 1.0/{o2}", simulate_em_config, omega1=1.0, omega2=o2)

print("\n--- strong asymmetry sweep ---")
for o2 in [0.5, 0.7, 1.0, 1.3, 1.5, 2.0]:
    R[f"strong_0.5/{o2}"] = run(f"strong 0.5/{o2}", simulate_strong_config, omega1=0.5, omega2=o2)

print("\n--- gravity asymmetry sweep (how mild an asym still flips?) ---")
for o2 in [0.8, 0.85, 0.9, 0.95, 1.0]:
    R[f"grav_0.8/{o2}"] = run(f"gravity 0.8/{o2}", simulate_gravity_config, omega1=0.8, omega2=o2, universal=True)

print("\n--- weak baseline (original omega hardcoded 1.0/1.0) ---")
R["weak_baseline"] = run("weak 1.0/1.0 (sym)", simulate_weak)

out = dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg,
                     elapsed_s=round(time.time()-t0,1)), results=R)
fn = "s9_doublecheck_flip_canary.json" if CANARY else "s9_doublecheck_flip_results.json"
json.dump(out, open(fn,"w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
