"""
S9 follow-up B' -- REAL PySR adjudication of leg (ii) (now that PySR
works in-container). Settles the open question follow-up B (curve-fit
fallback) could not: does INFO-031's strong = exp((H_b - H_a)) - const
reproduce with the actual tool, and is the EM/strong contrast regime-
stable across asymmetry?

Reuses the ORIGINAL fit_pysr (same operator set, loss, maxsize as the
Session 8 glances) imported from gravity_glance -- not a reimplementation.
MI vs (H_a, H_b), 2 seeds per config (INFO-028 cross-seed check).
"""
import sys, json, time
import numpy as np
from per_domain_kbk import build_ensemble_operator_matrix
from em_strong_glance import simulate_em_config, simulate_strong_config
from gravity_glance import simulate_gravity_config, fit_pysr
from four_force_probe import simulate_weak

CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)
NITER = 12 if CANARY else 30
TIMEOUT = 40 if CANARY else 60

CONFIGS = [
    ("EM_sym_1.0/1.0",    simulate_em_config,     dict(omega1=1.0, omega2=1.0)),
    ("EM_asym_1.0/1.2",   simulate_em_config,     dict(omega1=1.0, omega2=1.2)),
    ("EM_asym_1.0/1.5",   simulate_em_config,     dict(omega1=1.0, omega2=1.5)),
    ("strong_sym_0.5/0.5",simulate_strong_config, dict(omega1=0.5, omega2=0.5)),
    ("strong_asym_0.5/0.7",simulate_strong_config,dict(omega1=0.5, omega2=0.7)),
    ("strong_asym_0.5/1.5",simulate_strong_config,dict(omega1=0.5, omega2=1.5)),
]
if not CANARY:
    CONFIGS += [
        ("weak_sym_1.0/1.0",   simulate_weak,           dict()),
        ("gravity_sym_0.8/0.8",simulate_gravity_config, dict(omega1=0.8, omega2=0.8, universal=True)),
        ("gravity_asym_0.8/1.0",simulate_gravity_config,dict(omega1=0.8, omega2=1.0, universal=True)),
    ]

t0=time.time(); R={}
print("="*92); print(f"S9 PySR ADJUDICATION (real PySR) {'[CANARY]' if CANARY else '[FULL]'} seeds={seeds}"); print("="*92)
print("INFO-031 claim: EM -> (H_a-H_b)^2 + const ; strong -> exp((H_b-H_a) - const)\n")

for name, simfn, skw in CONFIGS:
    per_seed=[]
    for s in seeds:
        M,_ = build_ensemble_operator_matrix(*simfn(seed=s, **cfg, **skw))
        Ha, Hb, MI = M[:,0], M[:,1], M[:,5]
        n = len(MI); idx = np.linspace(0, n-1, min(600,n)).astype(int)
        X = np.column_stack([Ha[idx], Hb[idx]])
        fit = fit_pysr(X, MI[idx], niter=NITER, seed=s, timeout_seconds=TIMEOUT)
        low = [r["equation"] for r in fit["pareto_front"][:5]]
        per_seed.append(dict(seed=s, best=fit["best_equation_sympy"], low_complexity=low))
        print(f"{name:22s} seed{s}: best = {fit['best_equation_sympy']}")
        for r in fit["pareto_front"][:5]:
            print(f"      c={r['complexity']:2d} loss={r['loss']:.2e}  {r['equation']}")
    R[name]=dict(config=skw, per_seed=per_seed)
    print()

out=dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg, niter=NITER,
                   timeout_s=TIMEOUT, elapsed_s=round(time.time()-t0,1),
                   note="real PySR via original fit_pysr"), results=R)
fn="s9_pysr_adjudicate_canary.json" if CANARY else "s9_pysr_adjudicate_results.json"
json.dump(out, open(fn,"w"), indent=2)
print(f"Wrote {fn} ({out['meta']['elapsed_s']}s)")
