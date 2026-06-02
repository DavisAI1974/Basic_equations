"""
S9 follow-up B: re-examine leg (ii) -- the EM/strong EXPRESSION-level
functional-family contrast (INFO-031), which Session 9 flagged as
untrusted because (a) it came from a single config and (b) INFO-030
already showed functional families are regime-conditional.

PySR is unavailable (no Julia). This is a FIXED-LIBRARY CURVE-FIT
FALLBACK, NOT PySR: it cannot discover novel forms, but it can test
whether the SPECIFIC families INFO-031 reported actually best-fit, and
whether the EM-vs-strong contrast survives an asymmetry sweep.

For MI as a function of (H_a, H_b) it fits a fixed library and reports,
per config, the simplest family within 0.02 R^2 of the best (a low-
complexity Pareto pick). Families include the INFO-031 / INFO-025 forms:
  quad_diff : a*(H_a-H_b)^2 + c     (EM / physics Duffing family)
  exp_diff  : a*exp(b*(H_b-H_a)) + c (strong family)
  exp_Ha    : a*exp(H_a/b) + c       (biology family)
  + const, linear, single-channel quadratic baselines.

Original sims + original operator matrix reused. 3 seeds; per-seed best
family reported so cross-seed stability (INFO-028) is visible.
"""
import sys, json, time, warnings
import numpy as np
from scipy.optimize import curve_fit
from per_domain_kbk import build_ensemble_operator_matrix
from em_strong_glance import simulate_em_config, simulate_strong_config
from gravity_glance import simulate_gravity_config
from four_force_probe import simulate_weak
warnings.filterwarnings("ignore")

CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22, 33)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)

# family: (name, n_params, func(P, Ha, Hb), p0)
FAMILIES = [
    ("const",     1, lambda Ha, Hb, c: c*np.ones_like(Ha), [0.3]),
    ("lin_Ha",    2, lambda Ha, Hb, a, c: a*Ha + c, [0.5, 0.0]),
    ("lin_Hb",    2, lambda Ha, Hb, a, c: a*Hb + c, [0.5, 0.0]),
    ("lin_both",  3, lambda Ha, Hb, a, b, c: a*Ha + b*Hb + c, [0.3, 0.3, 0.0]),
    ("quad_Ha",   2, lambda Ha, Hb, a, c: a*Ha**2 + c, [0.2, 0.0]),
    ("quad_diff", 2, lambda Ha, Hb, a, c: a*(Ha-Hb)**2 + c, [0.3, 0.2]),
    ("exp_diff",  3, lambda Ha, Hb, a, b, c: a*np.exp(np.clip(b*(Hb-Ha), -20, 20)) + c, [0.2, 1.0, 0.0]),
    ("exp_Ha",    3, lambda Ha, Hb, a, b, c: a*np.exp(np.clip(Ha/b, -20, 20)) + c, [0.5, 2.0, 0.0]),
]

def r2(y, yhat):
    ss = np.sum((y - y.mean())**2)
    return 1.0 - np.sum((y - yhat)**2)/ss if ss > 0 else 0.0

def fit_all(Ha, Hb, MI):
    out = []
    for name, npar, f, p0 in FAMILIES:
        try:
            ff = (lambda HaHb, *p, _f=f: _f(HaHb[0], HaHb[1], *p))
            popt, _ = curve_fit(ff, (Ha, Hb), MI, p0=p0, maxfev=20000)
            yhat = ff((Ha, Hb), *popt)
            out.append((name, npar, float(r2(MI, yhat)), [float(x) for x in popt]))
        except Exception:
            out.append((name, npar, -9.99, None))
    return out

def pareto_pick(fits, tol=0.02):
    best = max(f[2] for f in fits)
    cands = [f for f in fits if f[2] >= best - tol]
    cands.sort(key=lambda f: (f[1], -f[2]))   # simplest, then best
    return cands[0], best

def run(name, simfn, **skw):
    picks = []
    for s in seeds:
        M, _ = build_ensemble_operator_matrix(*simfn(seed=s, **cfg, **skw))
        Ha, Hb, MI = M[:,0], M[:,1], M[:,5]
        fits = fit_all(Ha, Hb, MI)
        pick, best = pareto_pick(fits)
        picks.append((pick[0], pick[2], best))
    fam_seeds = [p[0] for p in picks]
    consensus = max(set(fam_seeds), key=fam_seeds.count)
    stable = "STABLE" if fam_seeds.count(consensus) == len(fam_seeds) else "SPLIT "
    r2m = np.mean([p[1] for p in picks])
    print(f"{name:24s} -> {consensus:10s} R^2={r2m:.3f}  {stable}  per-seed={fam_seeds}")
    return dict(name=name, consensus_family=consensus, r2_mean=float(r2m),
                per_seed_family=fam_seeds, stable=(stable.strip()=="STABLE"))

t0=time.time(); R={}
print("="*100); print(f"S9 EXPRESSION REFIT (curve-fit fallback, NOT PySR) {'[CANARY]' if CANARY else '[FULL]'} seeds={seeds}"); print("="*100)
print("INFO-031 claim: EM -> quad_diff (H_a-H_b)^2+c ; strong -> exp_diff exp((H_b-H_a))+c\n")

print("--- INFO-031 original configs ---")
R["em_asym_1.0/1.2"]   = run("EM asym 1.0/1.2", simulate_em_config, omega1=1.0, omega2=1.2)
R["em_sym_1.0/1.0"]    = run("EM sym 1.0/1.0", simulate_em_config, omega1=1.0, omega2=1.0)
R["strong_sym_0.5/0.5"]= run("strong sym 0.5/0.5", simulate_strong_config, omega1=0.5, omega2=0.5)
R["strong_asym_0.5/0.7"]=run("strong asym 0.5/0.7", simulate_strong_config, omega1=0.5, omega2=0.7)

print("\n--- EM across asymmetry (is the family regime-stable?) ---")
for o2 in [1.0, 1.2, 1.5, 2.0]:
    R[f"em_1.0/{o2}"] = run(f"EM 1.0/{o2}", simulate_em_config, omega1=1.0, omega2=o2)
print("--- strong across asymmetry ---")
for o2 in [0.5, 0.7, 1.0, 1.5]:
    R[f"strong_0.5/{o2}"] = run(f"strong 0.5/{o2}", simulate_strong_config, omega1=0.5, omega2=o2)
print("--- gravity / weak for reference ---")
R["gravity_sym"]  = run("gravity 0.8/0.8", simulate_gravity_config, omega1=0.8, omega2=0.8, universal=True)
R["gravity_asym"] = run("gravity 0.8/1.0", simulate_gravity_config, omega1=0.8, omega2=1.0, universal=True)
R["weak_sym"]     = run("weak 1.0/1.0", simulate_weak)

out=dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg, method="scipy curve_fit fixed library (PySR unavailable)",
                   elapsed_s=round(time.time()-t0,1)), results=R)
fn="s9_expression_refit_canary.json" if CANARY else "s9_expression_refit_results.json"
json.dump(out, open(fn,"w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
