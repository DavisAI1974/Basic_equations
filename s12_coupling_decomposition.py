"""
S12 #1 + #2 -- "where do the per-domain constraints live, and are the dipoles
coupled?" Decomposes each domain's null[0] across 5 seeds into three orthogonal
pieces and identifies the chemistry residual. Companion to s12_biology_coupling
(#3, the dynamical-knob test). Greg (S12): "find where the constraints live and
find how the dipoles are coupled which i think they are."

Three orthogonal probe axes in basis [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI]:
  - equal-entropy identity:  d_lin=(H_a-H_b)/sqrt2  and  d_quad=(-H_a^2-H_b^2
    +2 H_a*H_b)/sqrt6  (the geometric H_a~=H_b bookkeeping, INFO-036)
  - MI/coupling axis:        d_mi = MI
  - residual:                everything orthogonal to those (domain-specific)

Reads per_domain_kbk_results_seed{11,22,33,44,55}.json (extract_v1 null[0]).
Outputs s12_coupling_decomposition.json.
"""
import json
import numpy as np

SEEDS = [11, 22, 33, 44, 55]
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]
DOMS = ["physics_duffing", "biology_lotka", "chemistry_brussel", "geology_bk"]
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
BASIS = [d_lin, d_quad, d_mi]

vecs = {dom: [] for dom in DOMS}
for s in SEEDS:
    dd = json.load(open(f"per_domain_kbk_results_seed{s}.json"))
    for r in dd["results"]:
        v = np.array(r["extract_v1_v_null_6d"])
        vecs[r["name"]].append(v / np.linalg.norm(v))

out = {"seeds": SEEDS, "basis": OPS, "per_domain": {}}
print("=" * 92)
print(f"S12 per-domain null decomposition ({len(SEEDS)} seeds)")
print("=" * 92)
print(f"{'domain':10s} {'eqEntropy':>14s} {'MI/coupling':>15s} {'residual':>14s} {'dir|cos|':>9s}")
for dom in DOMS:
    V = vecs[dom]
    eq, mi, res = [], [], []
    res_dirs = []
    for v in V:
        a = (v @ d_lin) ** 2; b = (v @ d_quad) ** 2; c = (v @ d_mi) ** 2
        eq.append(a + b); mi.append(c); res.append(max(0.0, 1 - a - b - c))
        rv = v.copy()
        for bb in BASIS:
            rv = rv - (v @ bb) * bb
        n = np.linalg.norm(rv)
        if n > 1e-6:
            res_dirs.append(rv / n)
    eq, mi, res = map(np.array, (eq, mi, res))
    cs = [abs(V[i] @ V[j]) for i in range(len(V)) for j in range(i + 1, len(V))]
    # residual direction stability + mean (sign-aligned)
    rd = np.array(res_dirs)
    ref = rd[0]
    rd_al = np.array([x * np.sign(x @ ref) for x in rd])
    rcs = [abs(rd[i] @ rd[j]) for i in range(len(rd)) for j in range(i + 1, len(rd))]
    res_mean_dir = rd_al.mean(0)
    res_relation = " ".join(f"{res_mean_dir[i]:+.2f}*{OPS[i]}"
                            for i in range(6) if abs(res_mean_dir[i]) > 0.15) + " ~ 0"
    out["per_domain"][dom.split("_")[0]] = dict(
        eqEntropy_frac=[round(eq.mean(), 4), round(eq.std(), 4)],
        MI_coupling_frac=[round(mi.mean(), 4), round(mi.std(), 4)],
        residual_frac=[round(res.mean(), 4), round(res.std(), 4)],
        null_dir_stability_mean_abscos=round(float(np.mean(cs)), 4),
        residual_dir_stability_mean_abscos=round(float(np.mean(rcs)), 4),
        residual_relation=res_relation,
    )
    print(f"{dom.split('_')[0]:10s} {eq.mean():6.3f}+/-{eq.std():4.3f} "
          f"{mi.mean():7.3f}+/-{mi.std():4.3f} {res.mean():6.3f}+/-{res.std():4.3f} "
          f"{np.mean(cs):9.3f}")

out["reading"] = (
    "Where each per-domain constraint LIVES, 5 seeds with scatter: physics + "
    "geology are ~pure equal-entropy bookkeeping (no MI), reproducibly. BIOLOGY "
    "is 0.906+/-0.007 MI-COUPLING (MI~=0.28*H_a) -- genuine dynamical coupling, "
    "confirmed by the g-knob sweep (s12_biology_coupling: g=0 -> MI-frac 0.006, "
    "g>0 -> 0.81-0.97; slope tracks coupling strength). CHEMISTRY carries a "
    "stable 0.167+/-0.020 residual = 0.54*(H_a+H_b)+0.32*H_a^2-0.55*H_b^2 "
    "(total-entropy vs asymmetric quadratic, cross-seed |cos|=1.000) on top of "
    "0.83 equal-entropy. So 'the dipoles are coupled' is TRUE for biology, "
    "PARTIAL (residual, not MI) for chemistry, and FALSE (pure bookkeeping) for "
    "physics/geology. Refines INFO-039: opposition-as-artifact only for the "
    "equal-entropy domains; biology's coupling is real and dynamical."
)
print("\nresidual relations:")
for dom in DOMS:
    rec = out["per_domain"][dom.split("_")[0]]
    print(f"  {dom.split('_')[0]:10s} {rec['residual_relation']}  "
          f"(res dir |cos|={rec['residual_dir_stability_mean_abscos']})")
json.dump(out, open("s12_coupling_decomposition.json", "w"), indent=2)
print("\nWrote s12_coupling_decomposition.json")
