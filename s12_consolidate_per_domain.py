"""
S12 -- consolidate the FOUR per-domain governing equations (physics, biology,
chemistry, geology) into ONE OD-ingest JSON, built from in-repo result files
(per Greg, Session 12: "look in the basics repo for the 4 individual
equations and all relevant information ... construct to json. don't worry
about markets").

Two INDEPENDENT, REPRODUCIBLE per-domain signatures are pulled directly from
the committed result JSONs -- nothing hand-typed:

  (1) NULL DIRECTION (INFO-023, Session 6): the operator-space null[0] of the
      6-op basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI], from
      per_domain_kbk_results_seed{11,22}.json (KBK rank-gap + extract_v1).
  (2) MI-vs-H FUNCTIONAL FAMILY (INFO-025, Session 7): the PySR Pareto pick
      per domain, from pysr_symbolic_per_domain_results.json (x0=H_a, x1=H_b).

These are the per-domain EXPRESSIONS. The shared SUBSTRATE is the flow dipole:
the Level 1 windowed-null whose attractor (-1,-1,+2)/sqrt(6) in
(H_a^2,H_b^2,H_a*H_b) is -- as Sessions 5/10/11 established -- the
equal-marginal-entropy identity H_a~=H_b (a geometric statistics fact), NOT a
coupling. "Flow" (differential dH/dt = f) vs "algebraic" (f(H)=const) dipole
distinction is from static_dipole_test.py.

Discipline tags carried verbatim: every per-domain entry is a LOCATED finding
with cross-seed scatter; the substrate attractor reading is the deflationary
(equal-entropy) one confirmed on real LIGO data (INFO-036/038). Frames stay
separate from data. Provenance recorded per field.

Output: od_per_domain_equations.json
"""
import json
import numpy as np

OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]
DOM_ORDER = ["physics_duffing", "biology_lotka", "chemistry_brussel", "geology_bk"]
DOM_LABEL = {"physics_duffing": "physics", "biology_lotka": "biology",
             "chemistry_brussel": "chemistry", "geology_bk": "geology"}
SIMULATOR = {"physics_duffing": "Duffing oscillator",
             "biology_lotka": "Lotka-Volterra",
             "chemistry_brussel": "Brusselator",
             "geology_bk": "Burridge-Knopoff"}


def vnull_to_terms(v, tol=0.05):
    """Human-readable null relation: keep operator terms with |coef|>tol."""
    parts = [f"{v[i]:+.3f}*{OPS[i]}" for i in range(len(OPS)) if abs(v[i]) > tol]
    return " ".join(parts) + " ~ 0"


# ---- (1) NULL DIRECTIONS from the KBK per-domain JSONs (both seeds) -------
kbk = {s: json.load(open(f"per_domain_kbk_results_seed{s}.json")) for s in (11, 22)}
null_by_dom = {}
for dom in DOM_ORDER:
    seeds = {}
    for s in (11, 22):
        r = next(x for x in kbk[s]["results"] if x["name"] == dom)
        seeds[f"seed{s}"] = dict(
            v_null_6d=[round(float(x), 4) for x in r["extract_v1_v_null_6d"]],
            rank_null_subspace=r["kbk_rank_gap"]["rank_null_subspace_kbk"],
            min_eig_op_cov=float(r["eigenvalues_op_cov"][-1]),
        )
    v11 = np.array(seeds["seed11"]["v_null_6d"])
    v22 = np.array(seeds["seed22"]["v_null_6d"])
    cos = float(abs(v11 @ v22) / (np.linalg.norm(v11) * np.linalg.norm(v22)))
    null_by_dom[dom] = dict(
        seeds=seeds,
        cross_seed_cos=round(cos, 4),
        readable_relation_seed11=vnull_to_terms(v11),
    )

# ---- (2) MI-vs-H FUNCTIONAL FAMILY from the PySR Pareto fronts ------------
# The documented INFO-025 family pick per domain (the Pareto entry that names
# the reproducible class), pulled by matching the equation string.
pysr = json.load(open("pysr_symbolic_per_domain_results.json"))
FAMILY_PICK = {  # (domain): (complexity to pick, class label)
    "physics_duffing":   (6, "symmetric quadratic in difference: (H_b-H_a)^2 + const"),
    "biology_lotka":     (5, "exponential in H_a, H_b absent: ~exp(H_a/2)"),
    "chemistry_brussel": (5, "linear in H_a: a*H_a + b"),
    "geology_bk":        (1, "decoupled constant"),
}
fam_by_dom = {}
for dom in DOM_ORDER:
    cpick, label = FAMILY_PICK[dom]
    seeds = {}
    for s in ("seed11", "seed22"):
        lst = pysr["results_by_seed"][s]
        domobj = next(x for x in lst if x["domain"] == dom)
        ent = next((e for e in domobj["pareto_front"] if e["complexity"] == cpick), None)
        if ent is None:  # fall back to lowest-loss <= cpick
            ent = min((e for e in domobj["pareto_front"] if e["complexity"] <= cpick),
                      key=lambda e: e["loss"])
        seeds[s] = dict(complexity=ent["complexity"], loss=ent["loss"],
                        equation_x0Ha_x1Hb=ent["equation"])
    fam_by_dom[dom] = dict(family_class=label, seeds=seeds)

# ---- assemble per-domain records -----------------------------------------
domains = {}
for dom in DOM_ORDER:
    domains[DOM_LABEL[dom]] = dict(
        simulator=SIMULATOR[dom],
        null_direction=null_by_dom[dom],          # INFO-023 signature
        mi_vs_h_functional_family=fam_by_dom[dom],  # INFO-025 signature
        discipline_status="LOCATED (cross-seed scatter reported)",
    )

# ---- shared substrate: the FLOW DIPOLE -----------------------------------
flow_dipole = dict(
    name="flow dipole (Level 1 windowed-null, differential form)",
    operator_basis=OPS,
    substrate_attractor_dir_in_Ha2_Hb2_HaHb="(-1, -1, +2)/sqrt(6)",
    algebraic_identity="-(H_a - H_b)^2 ~ 0  (equivalently H_a^2 + H_b^2 ~ 2*H_a*H_b)",
    reading=("the attractor is the EQUAL-MARGINAL-ENTROPY identity H_a ~= H_b, "
             "a geometric statistics fact, NOT a coupling (Sessions 5/10/11; "
             "confirmed on real LIGO noise where whitened detector noise lands "
             "ON it and the chirp leaves it while MI -- a separate operator -- "
             "spikes). Deflationary read is the supported one."),
    flow_vs_algebraic=("FLOW dipole = differential (dH/dt = f). ALGEBRAIC dipole "
                       "= instantaneous constraint f(H)=const, e.g. markets "
                       "H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2. See "
                       "static_dipole_test.py."),
    ledger=["INFO-012", "INFO-014 (interp retracted)", "INFO-022 (rank-3 null subspace)",
            "INFO-036/038 (equal-entropy on real LIGO data)"],
)

# ---- preserved early-session ALGEBRAIC per-domain equations (Rule D) ------
# From SESSION_HANDOFF_2026-05-25.md, preserved through Session 5; algebraic
# (H_a^2 vs H_a*H_b) fits, distinct procedure from the null/PySR above.
algebraic_preserved = dict(
    chemistry_quadratic=dict(
        equation="H_a^2 = 0.007 - 0.093*(H_a*H_b) + 1.309*(H_a*H_b)^2",
        R2=0.943, note="algebraic dipole, chemistry"),
    geology_rank3_constraint=dict(
        equation="0.724*(H_a*H_b) - 0.441*H_b^2 - 0.290*H_a^2 ~ 0",
        residual_std_over_mean="0.15%", note="tight linear coupling H_a~=H_b"),
    biology_MI_poly=dict(
        equation="MI ~ polynomial(H_a)  [INCOMPLETE per Rule D: true family is "
                 "exponential ~exp(H_a/2), see INFO-025]",
        R2=0.66),
    nonlinear_deepdive_Ha2_vs_HaHb_bestcubic_R2=dict(
        physics=0.874, biology=0.448, chemistry=0.948, geology=0.994,
        note="best-cubic fit H_a^2 vs H_a*H_b per domain (May-25 deep dive); "
             "recorded from session figure, not a current repo JSON"),
)

out = dict(
    title="OD per-domain governing equations (4 sciences) -- consolidated",
    built_by="s12_consolidate_per_domain.py",
    source_files=["per_domain_kbk_results_seed11.json",
                  "per_domain_kbk_results_seed22.json",
                  "pysr_symbolic_per_domain_results.json",
                  "static_dipole_test.py (flow vs algebraic)",
                  "SESSION_HANDOFF_2026-05-25.md (preserved algebraic eqns)"],
    config=pysr["config"],
    operator_basis=OPS,
    note_x_mapping="In PySR equations x0=H_a, x1=H_b.",
    shared_substrate_flow_dipole=flow_dipole,
    per_domain=domains,
    preserved_algebraic_equations=algebraic_preserved,
    frame=("Per-domain differentiation has TWO independent reproducible "
           "signatures (null direction INFO-023; functional family INFO-025) "
           "on a shared equal-entropy substrate. Still a FRAME, not a claim. "
           "'don't worry about markets' (Greg, S12); markets algebraic dipole "
           "form retained in CLAUDE.md placeholder only."),
)
json.dump(out, open("od_per_domain_equations.json", "w"), indent=2)

# ---- console summary ------------------------------------------------------
print("=" * 92)
print("OD per-domain governing equations -- consolidated from in-repo JSONs")
print("=" * 92)
print(f"basis = {OPS}   (PySR x0=H_a, x1=H_b)")
print("\nSHARED SUBSTRATE (flow dipole): attractor (-1,-1,+2)/sqrt6 = "
      "equal-entropy identity -(H_a-H_b)^2 ~ 0\n")
for lab, rec in domains.items():
    nd = rec["null_direction"]
    fam = rec["mi_vs_h_functional_family"]
    print(f"-- {lab.upper():9s} ({rec['simulator']})")
    print(f"   null[0]: {nd['readable_relation_seed11']}")
    print(f"            cross-seed cos {nd['cross_seed_cos']}, "
          f"rank_null {nd['seeds']['seed11']['rank_null_subspace']}, "
          f"min_eig {nd['seeds']['seed11']['min_eig_op_cov']:.1e}")
    print(f"   MI-vs-H: {fam['family_class']}")
    print(f"            seed11 [c{fam['seeds']['seed11']['complexity']}] "
          f"{fam['seeds']['seed11']['equation_x0Ha_x1Hb']}")
    print(f"            seed22 [c{fam['seeds']['seed22']['complexity']}] "
          f"{fam['seeds']['seed22']['equation_x0Ha_x1Hb']}")
print("\nWrote od_per_domain_equations.json")
