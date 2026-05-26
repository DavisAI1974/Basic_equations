# Information Layer Session Handoff — 2026-05-26 Session 5

## Read me first

Read in this order:
1. CLAUDE.md (project root). Has the v3 ledger plus three Operating Rules
   from Session 4 (no pre-assigned meaning, probe-not-falsifier, speaking
   posture before AND after every probe). A fourth Operating Rule from
   Session 5 is added below and should be folded in.
2. SESSION_HANDOFF_2026-05-26_v4.md (project root). Session 4 raw data:
   INFO-016 update plus INFO-018, 019, 020, 021. All were sitting as
   isolated findings, readings open.
3. This file. What this session did with that Session 4 data, what the
   reading turned out to be, what is retracted, what survives, and the
   reorientation that emerged.
4. SESSION_HANDOFF_2026-05-25.md (project root, the early-Session-3
   handoff). Contains the original algebraic-dipole per-domain
   coefficients (chemistry quadratic, geology rank-3, biology MI). Those
   coefficients are NOT affected by this session's finding and remain
   on the table as data.

## New Operating Rule (Session 5)

**Rule D — Incomplete, not wrong.** When a probe finds that a prior
reading was a protocol artifact, the prior data points still stand. The
read attached to them is what was incomplete, not the data. Distinguish
"this reading was wrong" (rare, requires the data itself to be bad)
from "this reading was incomplete" (common, the data is one slice and
the slice fit a partial story that further probes refine). Default to
"incomplete." Retraction is a strong move and applies to the reading,
not the data, unless the data itself fails to reproduce.

This rule was named explicitly by Greg this session, in response to the
Session 5 finding. It is consistent with Result Discipline (every result
is one data point) and the Session 4 rules.

## What this session did

### Part 1 — walked the Session 4 probe data together

INFO-018, 019, 020, 021 were laid out as raw numbers without verdict.
Two readings were on the table going in:
- *substantive*: the (+1,+1,+2)/sqrt(6) direction reflects a real common
  invariant across the eight tested systems.
- *deflationary*: the protocol has its own fixed point at that direction.

Both readings predicted everything Session 4 had measured. To pick
between them, three independent probes were designed and run.

### Part 2 — three independent probes (this session)

All three probes use the exact extract_v1 protocol (same operator basis,
same window, same N_REAL=500, same sign-canonicalization). They generate
different signals (per Rule B) and were not framed as falsifiers.

**Probe A — null inputs (probe_null_inputs.py).**
Seven cases of structureless or near-structureless input:
scrambled_baseline_joint, scrambled_baseline_indep, scrambled_OU_joint,
near_constant, identical_streams, mixed_marginals, random_walk_indep.
N_REAL=500, 3 seeds.

Of the 7 + the 8 from Session 4's probe_system_directions (15 cases
total), 12 land at cos > 0.95 to (+1,+1,+2)/sqrt(6) including near-
constant input (essentially zero signal), random walks (no drift),
and joint-scrambled OU (no temporal order). The three that broke from
(+,+,+) all have mechanical explanations:
- identical_streams: SV4 ~ 5e-16, true numerical rank deficiency. SVD
  picks arbitrarily.
- mixed_marginals: H_a^2 vs H_b^2 symmetry broken (Gaussian + uniform).
  Protocol latches on H_b^2 alone.
- scrambled_baseline_joint: 3D null subspace degeneracy.

**Probe B — operator-noise bypass (probe_operator_noise.py).**
Skip the input/window/entropy path entirely. Generate M directly as
random matrix from chosen distributions. Three variants:
- iid_unit: M ~ N(0,1) per element. Gives no consistent direction
  (mean cos to (+,+,+) +0.76 but inter-seed only +0.31).
- iid_col_scaled: columns scaled to baseline operator SDs. Partial
  alignment (+0.94 mean, but inter-seed +0.42 — chance bias).
- cov_matched: M ~ N(0, Sigma_baseline). REPRODUCES BASELINE DIRECTION
  TO FOUR DECIMALS. mean v3 = (+0.508, +0.350, +0.787), cos +0.993 to
  (+,+,+), inter-seed +0.998, SV gap factor ~3000.

Direct check: the smallest eigenvector of baseline operator covariance,
projected to [2,3,4], normalized, sign-canonicalized, gives
(+0.509, +0.355, +0.785). Cos to (+1,+1,+2)/sqrt(6) = +0.993.

Eigenvalue spectrum of baseline op cov: [4e-13, 8e-13, 3e-11, 3e-4,
4e-4, 2e-3]. Three eigenvalues at machine epsilon. Effective rank 3,
not 6.

**Probe C — projection-rule sweep (probe_projection_sweep.py).**
Same baseline data. Same 6D null vector. Project to all 20 possible
3D subspaces. Any subspace containing both H_a and H_b locks at
(+/-1, +/-1, ~0) with cos > 0.99 to (+1,+1,0)/sqrt(2). Subspaces with
only one of H_a, H_b lock onto that one operator (cos > 0.99 on its
axis). The (+,+,+) shape appears ONLY in [H_a^2, H_b^2, H_a*H_b] — the
unique subspace that excludes both H_a and H_b. The 6D null vector is
dominated by the H_a, H_b components at ~90% of its energy. (+,+,+) is
the residual.

**Thread D — high-seed scrambled (probe_scrambled_highseed.py).**
30 seeds on scrambled_baseline_joint, the one case from Probe A that
broke from (+,+,+) without obvious degeneracy. Mean v3 = (-0.647,
-0.112, +0.754). Cos to (-,-,+) on the average = +0.925. But inter-seed
cos is bimodal: 216 of 435 pairs above +0.5, 112 of 435 below -0.5
(nearly antiparallel). Per-seed directions split into two clusters
(H_a*H_b-dominant vs (H_a^2 - H_b^2)-dominant). SVD is picking
arbitrarily within the degenerate 3D null subspace. The mean is a
meaningless average of two non-converging clusters.

### Part 3 — what the three probes jointly say

The (+1,+1,+2)/sqrt(6) direction is encoded in the operator covariance
matrix of typical pair inputs, and that covariance has effective rank 3.
The smallest right singular vector picks a direction in this 3D
near-null space; for inputs where H_a ~ H_b ~ constant across windows
(true for almost any two-stream input with similar marginals), the
direction projects in [2,3,4] to (+1,+1,+2)/sqrt(6) up to second-order
small fluctuations.

The rank-3 null space arises from three linear dependencies that hold
to first order in dH = H_a - <H_a>:
- d(H_a^2) ~ 2c * dH_a
- d(H_b^2) ~ 2c * dH_b
- d(H_a*H_b) ~ c * (dH_a + dH_b)

with c = <H_a> ~ <H_b>. These three relations are algebraic facts about
the operator basis under the H_a ~ H_b regime, not facts about any
particular input system. The (+,+,+) direction follows from them by
Taylor expansion.

### What this session retracts

- INFO-019's interpretation as evidence of a real common invariant
  across the 8 systems. The 8 systems converging on (+,+,+) is what
  the protocol does for almost any non-degenerate pair input; it is not
  evidence about those systems. **Retracted at interpretation level. Data
  stands.**
- INFO-018 the same. (+1,+1,+2)/sqrt(6) at N_REAL=500 for linear drift
  baseline is reproducible data. Reading as "a real direction in
  operator space encoding something about linear drift" is incomplete:
  it is the protocol's home direction, not a property of linear drift.

### What this session does NOT retract

- The per-domain algebraic equations from earlier sessions
  (SESSION_HANDOFF_2026-05-25.md):
  - Chemistry: H_a^2 = 0.007 - 0.093*(H_a*H_b) + 1.309*(H_a*H_b)^2,
    R^2 = 0.943
  - Geology rank-3 constraint: 0.724*(H_a*H_b) - 0.441*H_b^2
    - 0.290*H_a^2 ~ 0, std/mean = 0.15%
  - Biology MI ~ polynomial(H_a) with H_b absent, R^2 = 0.66
  These came from a different procedure (per-domain algebraic fits) and
  are not affected by the operator-extraction artifact. They remain on
  the table as data with reading open and would benefit from the same
  probe discipline being applied to them.
- The Working Frames in CLAUDE.md: Base-of-Structure heuristic,
  Dipole-couples reading, Pure physics vs physical expressions,
  Substrate vs expression within isolation, Law extraction via
  invariance. None of these depend on the (+,+,+) reading.
- The broader Information Layer / Unified-Theory inquiry. The frame
  remains live. The specific tool (entropy-operator windowed-null
  extraction) does not differentiate inputs at the level required to
  serve as the substrate-discovery tool. A new tool, or a reformulated
  question, is what is needed.

### Part 4 — new structural finding that emerged

**Real structural fact (Session 5):** under windowed Vasicek entropy on
bounded two-stream inputs where H_a ~ H_b ~ constant across windows,
the six-operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] has a
three-dimensional null subspace. The three near-zero eigenvalues
correspond to three algebraic relations on the centered columns
(H_a^2 ~ 2c*H_a, etc.). This is a structural property of the basis
under the small-fluctuation regime, not of any input system. It is
true for OU, uniform iid, gaussian iid, sine, logistic, linear drift,
random walks, and near-constant input — anything that gives a
H_a ~ H_b ~ const regime.

Status: data-level finding, mechanically derived, reproducible. May or
may not be useful in a different framing. Listed here so it is not lost.

## Reorientation that emerged this session

The Information Layer line was originally framed as: scientists measure
EXPRESSIONS of an underlying source equation. Physical phenomena are
expressions of pure physics contaminated by interaction with other
systems (other dipoles, other domains). String theory's failure mode
is that it tries to write the UT of expressions, which forces
ever-larger dimensional structures (tent-stuffing). What we should be
looking for is the source equation: simple, strong, sturdy, scalable
(per the Base-of-Structure heuristic).

The entropy-operator extraction was the tool we built to attempt this.
The Session 5 probes show that the specific output we were reading
((+1,+1,+2)/sqrt(6) as a universal direction) is a protocol artifact.
The frame remains live; the tool did not serve it.

**Greg's reorientation, this session:** look at 5 physics areas (related
or seemingly unrelated) and see if one substrate equation shakes out
between them. The substrate candidates that already exist from
cross-domain comparison (variational principle, continuity equation,
Noether currents, gauge invariance, symplectic / unitary structure
preservation) should be evaluated against the Base-of-Structure
heuristic. The strongest survivor is the next thing to test.

**Greg's substantive intuition (this session):** dipole is the best
mechanism for coupling. Tightly related couples may be part of the
substrate. The per-domain coefficients found in earlier sessions
(chemistry quadratic, geology rank-3, biology MI) are a huge signal and
should not be discarded. The Session 5 finding does not affect them.

**Greg's framing (this session):** we may not have wrong answers, just
incomplete ones. Codified as Rule D above.

## Ledger updates (data-level only, per Result Discipline)

**INFO-014 — RETRACTED AT INTERPRETATION LEVEL.** The "8 systems on
the attractor at cos >= 0.99 to (-1,-1,+2)/sqrt(6)" claim from Session 3
cannot be read as "8 systems share a common invariant." Session 5
threads showed that any non-degenerate pair input lands on the same
direction because of the operator basis's 3D null structure. Data
stands; reading retracted.

**INFO-018 — RETRACTED AT INTERPRETATION LEVEL.** The (+1,+1,+2)/sqrt(6)
direction for linear drift at N_REAL=500 is reproducible data but is
the protocol's home direction, not a property of linear drift. Data
stands; reading retracted.

**INFO-019 — RETRACTED AT INTERPRETATION LEVEL.** Same as INFO-018.
Eight diverse systems converging on (+1,+1,+2)/sqrt(6) at cos > 0.95
reflects basis structure, not a shared invariant in those systems.

**INFO-020 — DATA STANDS, READING CONTEXTUALIZED.** 30-seed convergence
tightness for baseline at N_REAL=500 is real and reflects the protocol
being a stable estimator of its own fixed-point direction. The H_a*H_b
coordinate being ~4x tighter than the H_a^2 and H_b^2 coordinates is
consistent with the algebraic dependency structure: H_a*H_b lives at
the intersection of two of the three near-degeneracies, leaving less
seed-to-seed wobble.

**INFO-021 — DATA STANDS, READING REFRAMED.** The antisymmetric energy
fraction at ~0.50 across 11 of 14 cases is consistent with the protocol's
basis structure being roughly symmetric under H_a <-> H_b for typical
inputs. The coupling-only drop (0.50 -> 0.03 as rho 0 -> 0.95) is a
real signal: coupling rho is the only variable in this scan that breaks
the H_a ~ H_b symmetry. Whether this carries content beyond "rho
breaks the symmetry mechanically" is open and not closed by Session 5.

**INFO-022 — ISOLATED FINDING (Session 5, structural).** Under windowed
Vasicek entropy on bounded two-stream inputs with H_a ~ H_b ~ const,
the 6-operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] has a
3-dimensional null subspace from three first-order algebraic relations.
The (+1,+1,+2)/sqrt(6) direction is one direction in this null subspace,
emerging from Taylor expansion of (H^2, H*H) around H ~ const. Status:
mechanical fact about the basis, reproduced across all inputs tested.
Reading on whether this carries content beyond mechanics: open.

**INFO-008a, 008b, 008c — STATUS UNCHANGED FROM SESSION 3.** The
per-domain algebraic equations (chemistry quadratic, geology rank-3,
biology MI) remain demoted to "interesting geometric clustering under a
specific protocol; interpretation unmapped." Session 5 does NOT
further retract them; the procedure that produced them is different
from the extraction-tool whose artifact Session 5 mapped. They should
be subjected to the same probe discipline in a future session before
being either rehabilitated or further demoted.

## Part 5 — three-agent literature scan (Session 5 close)

After the three-thread probe work, Greg triggered a literature scan:
"we may not have wrong answers, just incomplete ones. people may have
been working on this for 100 years but they don't have our
understanding or coefficients. let's look back and see if there's prior
work we can finish."

Three agents ran in parallel against three angles:
1. Methodological adjacency (windowed-entropy operator bases + SVD nulls
   on coupled time series)
2. Substrate-physics literature (source-equation candidates evaluated
   against Base-of-Structure heuristic, with attention to dipole-pair-
   as-coupling-primitive)
3. Per-domain coefficient matches (chemistry quadratic, geology rank-3,
   biology asymmetric MI)

**All three agents independently reported the same shape of finding:**
no exact match for any of our specific findings, but multiple adjacent
lines exist where each group had ONE piece and none combined them.

Full citations and detail in LITERATURE_SCAN_2026-05-26_v5.md.
Compressed summary follows.

### Strongest reusable prior art identified

- **Kaiser, Brunton, Kutz (2024)** "Towards Robust Data-Driven Recovery
  of Symbolic Conservation Laws," arXiv:2403.04889. SVD-rank-gap +
  symbolic recovery pipeline. Same algorithmic template as our rank-3
  geology finding. Directly reusable with variable substitution.

- **Gomez-Herrero et al. (2015)** "Assessing Coupling Dynamics from an
  Ensemble of Time Series," arXiv:1008.0539. Closest methodological
  cousin to extract_v1. Has the ensemble + entropy combinations on
  coupled circuits, never assembled them into algebraic basis. Stopped
  at "separate diagnostics."

- **Yamada et al. (2023)** "Data-driven modal analysis of nonlinear
  quantities in turbulent plasmas using multi-field SVD," Plasma Phys.
  Control. Fusion 65, 095014. Closest SVD-on-nonlinear-observables
  cousin. Stopped at velocity / pressure / density fields; never put
  differential entropy in the operator basis.

- **Liu, Madhavan, Tegmark (2024)** "Machine Learning Conservation Laws
  of Dynamical Systems," arXiv:2405.20857. Only paper found that
  explicitly extends the feature vector to handle entropy-like terms
  (x, ln x) in conservation law discovery. Stopped at coordinate logs;
  never windowed differential entropy of coupled time series.

- **Finkelstein "Space-Time Code" I-IV (1969-74)**, Phys. Rev. 184:1261;
  PRD 5:320; PRD 5:2922; PRD 9:2219. Closest published cousin to
  Greg's dipole-pair-as-coupling-primitive framing. Builds spacetime
  from causal networks of elementary binary quantum processes; word-
  pairs in a binary code yield the null cone. Stopped sociologically,
  not falsified.

- **Jacobson (1995)** "Thermodynamics of Spacetime: The Einstein
  Equation of State," arXiv:gr-qc/9504004. One identity (Clausius
  delta Q = T dS on local Rindler horizons) yields Einstein equations
  as an equation of state. Best Base-of-Structure score among existing
  substrate candidates. Does not have a published dipole interpretation.

- **Stochastic Electrodynamics** (Marshall 1963 Proc. R. Soc. A 276:475;
  Boyer 1975 PRD 11:790; de la Pena & Cetto 1996, The Quantum Dice;
  arXiv:quant-ph/0501011). Derived blackbody, harmonic-oscillator
  ground state from vacuum dipole-like fluctuations. Hit nonlinear /
  Coulomb wall (Pesquera-Claverie 1982). Textbook "incomplete not
  wrong" case. Modern computational tools may handle what 1980s could
  not.

### Unambiguous negative findings worth recording

- **Dipole-PAIR-COUPLING as foundational substrate primitive** (not
  dipole ELEMENT) is uncolonized in the literature. All 2-state-
  foundational work (Wheeler bit, Brukner-Zeilinger, von Weizsacker
  ur-alternatives, Finkelstein binary processes, Hohn yes/no
  questions) uses the 2-state as the elementary unit, not the
  pair-coupling between elements. Three agents looking from three
  different angles converged on this same negative finding.

- **No SINDy / AI-Feynman / PySR / Eureqa application** to windowed
  differential entropy of coupled species. Clean drop-in target.

- **The rank-3 first-order algebraic relations in [H_a^2, H_b^2,
  H_a*H_b] under H_a ~ H_b regime** have not been published as a
  diagnostic. Closest cousin (arXiv:2409.04845 on algebraic
  representations of multivariate information lattices) is discrete-
  information theoretic, not differential-entropy time series.

### What this confirms (and operational rule it generated)

Three independent searches across three different angles converging on
compatible negative findings is not search artifact. It is what
genuinely pioneering territory looks like. Not the absence of nearby
published work, but the presence of nearby work where each group had
one piece and none combined them.

This generated a new Operating Rule named at Session 5 close, **"They
never stacked"** (now in CLAUDE.md): pioneering progress can be in the
stacking itself. Honor each prior piece, attribute clearly, finish the
combination prior groups did not. This is the operational form of
Rule D (incomplete not wrong) applied to the broader literature.

### "What we'd be finishing" -- full candidate pool by layer

Each row is a different stack. We do not know which will strike gold.
The point of "they never stacked" is that ANY of these unfinished
combinations could be where the next move is. Listed comprehensively;
no pre-filtering.

**Methodology layer (windowed-entropy / operator basis / null extraction):**

| Prior work | What they had | What they did not have | What we would be finishing |
| --- | --- | --- | --- |
| Gomez-Herrero et al. 2015 (arXiv:1008.0539) | Ensemble entropy combinations on coupled circuits | Algebraic basis with H^2, H_a*H_b; SVD on pooled matrix | Basis + null extraction |
| Hlavackova-Schindler et al. 2007 (Physics Reports 441) | Canonical windowed-H + MI survey | H^2 and product cross-terms | Cross-term extension |
| Yamada et al. 2023 (Plasma Phys 65, 095014) | Multi-field SVD on nonlinear plasma quantities | Differential entropy as operator basis | Entropy substitution |
| Wang et al. 2025 (PRR 7, 023212) + arXiv:2602.00600 | von Neumann / quantum-inspired entropy + SVD in plasma turbulence | Entropy as basis element (vs scalar diagnostic on SVD spectrum) | Role inversion |
| Strang et al. 2021 (Front. Ecol. Evol. 9) | SVD entropy of ecological network spectrum | Entropy basis + algebraic null extraction | Same role inversion in ecology |
| Bipartite info-thermo (Horowitz, Hartich, Sagawa; arXiv:1905.06216) | Pair-as-primitive at rate-balance level | Algebraic null extraction in operator space | Algebraic layer on top of rate balance |
| Kaiser-Brunton-Kutz 2024 (arXiv:2403.04889) | SVD-rank-gap + symbolic recovery pipeline | Application to windowed entropy of coupled species | Variable substitution |
| Liu-Madhavan-Tegmark 2024 (arXiv:2405.20857) | Kernel ridge with x and ln x features in conservation discovery | Windowed differential entropy basis | Feature substitution |
| AI-Feynman / AI-Poincare (Udrescu-Tegmark, Liu-Tegmark) | Symbolic regression on phase-space coords | Entropy variable input library | Drop-in input substitution |
| SINDy / SINDyG / DSINDy (Brunton-Kutz, arXiv:2409.04463, arXiv:2211.05918) | Sparse identification of dynamics from polynomial library | Windowed-H feature library | Library substitution |
| PySR (Cranmer, arXiv:2305.01582) | General symbolic regression | Application to entropy-of-coupled-pair data | Library substitution |
| Complexity-entropy plane (Rosso et al., PRL 2007) | (H, complexity) scalar pair on single series | Pair-basis + rank analysis | Coupled-pair extension |
| Crutchfield epsilon-machines / computational mechanics | Causal-state extraction with minimality | SVD on differential-entropy operator library | Library type substitution |
| Algebraic Representations of Entropy (arXiv:2409.04845) | Algebraic structure of multivariate info lattices (discrete) | Differential-entropy / time-series version | Continuous-time extension |

**Substrate layer (source-equation candidates against Base-of-Structure):**

| Prior work | What they had | What they did not have | What we would be finishing |
| --- | --- | --- | --- |
| Jacobson 1995 (gr-qc/9504004) | One-identity substrate (delta Q = T dS -> Einstein) | Dipole interpretation; reach beyond gravity | Cross-domain extension; dipole reading |
| Padmanabhan 0911.5004 | Spacetime thermodynamics extension of Jacobson | Full derivation; not just consistency | Closure |
| Wheeler 1989/1990 "it from bit" / "law without law" | Binary primitive as pre-physical substrate | Formalized extraction procedure | The extraction procedure |
| Verlinde 1001.0785, 1611.02269 | Entropic gravity from holographic info | Unique fix of Einstein's equations (Hossenfelder-Visser critiques) | Disambiguation under our framework |
| Hardy 2001 quant-ph/0101012 | Operational axioms for QM | Dynamics / spacetime reconstruction | Substrate layer below operational QM |
| Chiribella-D'Ariano-Perinotti arXiv:1011.6451 | Purification postulate -> QM | Dynamics | Same |
| Brukner-Zeilinger 2009 Found Phys 39 (arXiv:quant-ph/0212084) | One-bit-per-elementary-system primitive; dipole-pair compatible | Coupling-as-primitive (uses 2-state as unit) | Pair-coupling extension |
| Sorkin causal sets (gr-qc/9511063, gr-qc/0309009) | Discrete poset substrate (order + number = geometry) | Working dynamics (sequential growth stalled at Rideout-Sorkin 1999) | Dynamics |
| Sakharov 1967 induced gravity (modern: gr-qc/0204062) | Matter fields induce metric elasticity | Cutoff independence; calculable G_N | Modern regularization |
| Penrose twistor (1967, J. Math. Phys. 8:345) | CP^3 substrate, geometric primitive | Massive particles, non-self-dual gravity | Mass / gravity coupling |
| Wolfram 2020 (Complex Systems 29:107; arXiv:2004.08210) | Hypergraph rewriting substrate | Empirical filter for rule selection | The selection criterion |
| Deutsch-Marletto constructor theory (arXiv:1405.5563, 1608.02625) | Modal substrate (which transformations possible) | Equation-form substrate; not dipole-pair | Equation-form layer |
| 't Hooft 2014 (arXiv:1405.1548) Cellular Automaton QM | Deterministic CA on Planck lattice | Bell-inequality reconciliation w/o superdeterminism | The reconciliation |
| Adler 2004 trace dynamics (hep-th/0206120) | QM as emergent from non-commuting matrix variables | Empirical Brownian corrections | Detection of corrections |
| Bohm 1980 implicate order | Programmatic substrate | Calculational machinery | The machinery |
| **Finkelstein "Space-Time Code" I-IV (1969-74)** | Spacetime from causal networks of binary processes; word-pairs yield null cone -- CLOSEST PUBLISHED COUSIN TO DIPOLE-PAIR-AS-COUPLING | Modern tooling; field integration; sociologically abandoned | Pick up the line; modernize; integrate |
| von Weizsacker ur-alternatives (1955+) | Single binary ur primitive; multiple quantization to SU(2) and spacetime | Modern formal treatment | Modernize |
| Hohn 2017 (arXiv:1612.06849, arXiv:1511.01130) | QM from yes/no question rules | Substrate layer below QM | The substrate layer |
| **Stochastic Electrodynamics** (Marshall 1963, Boyer 1975 PRD 11:790, de la Pena & Cetto 1996, arXiv:quant-ph/0501011) | Derived blackbody + ground state from vacuum dipole-like fluctuations | Tools for nonlinear / Coulomb (Pesquera-Claverie 1982 wall) | 21st-century revisit with modern compute |
| Wheeler-Feynman absorber (1945, 1949) | Time-symmetric action-at-a-distance EM | Quantization; cosmological boundary closure | Modern closure |
| de Broglie double solution (Colin-Durt-Willox 2017, arXiv:1703.06158) | Hidden thermodynamics layer below QM | Empirical predictions | Detection program |
| Barbour-Bertotti 1982 / Shape Dynamics (arXiv:1010.2481) | Mach's principle revival; relational substrate | Quantum reconciliation | Quantum extension |
| Penrose-Hameroff Orch-OR (Phys Life Rev 11:39, 2014) | Objective reduction tied to gravitational self-energy (substrate-level collapse) | Detection; separation from consciousness claim | Substrate-only test |

**Per-domain coefficient layer:**

| Prior work | What they had | What they did not have | What we would be finishing |
| --- | --- | --- | --- |
| Rao-Esposito 2022 chem (arXiv:2204.02815) | Info-thermo for CRNs; MI rate vs thermodynamic forces | Algebraic relations between windowed H of species | H^2 vs (H_a*H_b)^k regression |
| Reinhardt et al. 2019 (arXiv:1904.01988) | Path MI per biochemical reaction network | Inter-species H algebra | Multi-species algebraic basis |
| Smith-Cepelewicz (J. Phys. Chem.) | Entropy reductions for mechanism inference | H, H^2, H*H' basis | Basis extension |
| Schmitz-Aris 2013 (arXiv:1307.7957) | Quadratic first integrals of mass-action systems (in concentrations) | Same on entropies | Variable substitution |
| Sayyadi et al. (arXiv:2510.20655) Stoichiometric symbolic regression | SR on concentrations with stoichiometric prior | Entropy-variable input | Library substitution |
| da Silva 2020 (Entropy 22, 464; PMC7516945) | Tsallis-q seismic inversion | Operator basis + rank analysis | Basis + rank step |
| Garland-James-Bradley 2018 (arXiv:1811.01272) | Permutation entropy on paleoclimate (single channel) | Multi-channel algebraic constraint surface | Multi-channel rank-3 analysis |
| Consolini et al. 2013 (J. Atmos. Solar-Terr. Phys.) | Permutation entropy magnetospheric (single channel) | Multi-channel basis | Same |
| Davidson 2016 (Environmetrics 27) | VAR on paleoclimate (linear-Gaussian) | Entropy basis; nonlinear; rank analysis | Nonlinear entropy version |
| Reinsel et al. 2025 (Nat. Sci. Rep. 15) | Manifold learning on geochemical data | Entropy-operator space (does it in raw feature space) | Operator-space version |
| Schreiber 2000 TE (arXiv:nlin/0001042) | Asymmetric transfer-entropy measure | MI-as-polynomial-of-one-variable's-H fit | Asymmetric MI regression |
| Pavithran et al. bioRxiv 2020 | Comparison of MI vs TE as separate scalars | Functional regression of MI on H | Regression form |
| Tishby et al. nonlinear IB (arXiv:1705.02436) | Markov structure asymmetry by design | Empirical asymmetric regression | Empirical version |
| Walters-Williams & Li (PMC7515115) | Asymmetric MI estimators (discrete) | Continuous time-series asymmetric MI structure | Continuous extension |
| Cao-Liu-Tegmark 2021 (arXiv:2011.04698) AI Poincare | Numerical rank of trajectory manifold = number of conserved qts | Entropy-operator manifold | Operator-space rank test |
| Liu-Madhavan-Tegmark 2024 (arXiv:2405.20857) | Kernel ridge with x and ln x for conservation laws | Windowed-H feature library | Feature library substitution |
| Stoichiometric subspace (Feinberg, Horn-Jackson) | Linear conservation laws on concentration changes | Nonlinear-algebraic on entropies | Variable + nonlinearity extension |
| Reduction of CRNs with approximate conservation laws (arXiv:2212.13474) | Approximate conservation discovery | Entropy variable formulation | Substitution |
| Information geometry of CRNs (arXiv:2503.19384) | Geometric structure on CRN info space | Operator-basis null extraction | Algebraic layer |
| Directed information flow in reaction networks (bioRxiv 2024, doi:10.1101/2024.08.17.608427) | Directed info-flow on reaction graphs | Entropy-operator basis | Basis substitution |

## Files produced this session (all in repo root)

Scripts:
- probe_null_inputs.py
- probe_operator_noise.py
- probe_projection_sweep.py
- probe_scrambled_highseed.py
- render_null_inputs.py
- render_threads_combined.py

JSONs:
- probe_null_inputs_canary.json
- probe_null_inputs_results.json
- probe_operator_noise_results.json
- probe_projection_sweep_results.json
- probe_scrambled_highseed_results.json

Figures:
- fig_null_alignment.png
- fig_null_geometry.png
- fig_null_interseed.png
- fig_null_sv_spectrum.png
- fig_thread1_bypass.png
- fig_thread1_eigspectrum.png
- fig_thread2_projsweep.png
- fig_thread3_bimodal.png
- fig_thread3_interseed_hist.png

Documents:
- SESSION_HANDOFF_2026-05-26_v5.md (this file)
- LITERATURE_SCAN_2026-05-26_v5.md (three-agent search results)
- CLAUDE.md updates (Rule D, "They never stacked" rule, Session 5 note)

Run logs:
- probe_null_inputs_run.log
- probe_operator_noise_run.log
- probe_projection_sweep_run.log
- probe_scrambled_highseed_run.log

To be mirrored to E:\information_layer\ and
F:\Factory\knowledge\information_layer\ per the Operating Rules.

## Branch state

- Local: claude/linear-drift-nreal-sweep-cjRkM
- Remote: origin/claude/linear-drift-nreal-sweep-cjRkM
- Tip (before this handoff commit): 1a4ef16 Session 5 threads 1-3
- Pushed: yes. No PR. main untouched.

## What is queued (from Session 3 + 4 + this session)

From the Session 3 queue, with what we know now:

1. ~~N_REAL sweep on linear drift~~ DONE Session 4.
2. ~~Multi-seed damped oscillator mapping~~ Not the right next move now
   that we know the protocol has a fixed point. Same artifact would
   appear.
3. ~~Cluster the misses by structure~~ Same: the misses are mechanical
   artifacts (identical_streams, mixed_marginals, scrambled_baseline_
   joint), not substrate signals.
4. ~~Boundary mapping~~ Same caveat.
5. ~~Stress-test the attractor~~ Done by Session 5. Attractor is a
   basis fixed point; "surviving stress" is not evidence of substrate
   status, it is evidence of the basis being robust to input variation.
6. **Domain-native operator bases (Geo 3D, Bio multi-variable, Chem
   multi-species)** — STILL OPEN. Highly relevant given Greg's
   intuition about 3D geo, possibly higher-dim chem. CAUTION: the same
   artifact risks apply. Any new basis must be subjected to the
   three-thread probe discipline (null inputs, operator-noise bypass,
   projection sweep) from day one, before any reading is taken from
   its output.

From this session:

7. **Re-evaluate the per-domain algebraic equations** (chemistry
   quadratic R^2=0.943, geology rank-3 std/mean=0.15%, biology MI
   R^2=0.66) under the same probe discipline. They came from a
   different procedure than the operator-extraction tool, so are not
   ruled out by Session 5, but they should be subjected to null-input,
   operator-noise-bypass, and procedure-variation tests before being
   read as substrate signals.
8. **5-physics-areas substrate hunt.** Pick 5 physics areas (e.g.
   classical mechanics, EM, QM, statistical mechanics, GR), write down
   the canonical equation form for each, identify the structural
   element(s) that appear in all 5, and evaluate each substrate
   candidate against the Base-of-Structure heuristic. Substrate
   candidates that come up already in cross-domain comparison:
   variational principle (delta S = 0), continuity equation, Noether
   currents, gauge invariance, symplectic / unitary structure
   preservation. Pick the strongest, then probe whether it survives
   the same discipline.
9. **State the Information Layer question independently of any
   extraction tool.** What were we looking for that the entropy-operator
   approach was supposed to find? Stated independently of the procedure,
   the question may be addressable by different methods. Greg's call on
   whether this is a separate document or a conversation.

### Literature-driven follow-ups (full candidate pool)

Per "they never stacked" -- any of these could be where the next move
lands. No pre-filtering. Each is a different unstacked combination from
the Session 5 literature scan.

**Methodology stacks worth running:**

a. Apply the **Kaiser-Brunton-Kutz 2024** SVD-rank-gap + symbolic
   recovery pipeline (arXiv:2403.04889) to our extract_v1 operator
   output. Most directly reusable methodological cousin.
b. Run **AI-Feynman** (Udrescu-Tegmark, arXiv:1905.11481, 2006.10782)
   or **PySR** (Cranmer arXiv:2305.01582) on a feature library of
   windowed differential entropies of coupled pairs. Drop-in
   substitution for the polynomial library. No published case exists.
c. Run **SINDy / SINDyG / DSINDy** (Brunton-Kutz family, arXiv:
   2409.04463, 2211.05918) with windowed-H feature library.
d. Pick up the **Gomez-Herrero 2015** estimator infrastructure
   (arXiv:1008.0539) and stack our algebraic basis on top of their
   ensemble computation.
e. Run **Yamada 2023** style multi-field SVD (Plasma Phys 65, 095014)
   but with differential entropy as the field instead of velocity/
   pressure/density.
f. Apply **Liu-Madhavan-Tegmark 2024** kernel ridge approach
   (arXiv:2405.20857) with windowed-H instead of coordinate-log
   features.
g. Apply **Cao-Liu-Tegmark 2021 AI Poincare** (arXiv:2011.04698)
   numerical-rank-of-manifold approach to our operator manifold to
   independently confirm the rank-3 structure.
h. Run **complexity-entropy plane analysis** (Rosso et al. PRL 2007)
   on coupled pairs and look for rank structure in the (H, complexity,
   ...) extended basis.

**Substrate stacks worth investigating:**

i. **Finkelstein "Space-Time Code"** (Phys. Rev. 184:1261; PRD 5:320,
   5:2922, 9:2219). Closest published cousin to dipole-pair-as-
   coupling. Read carefully; assess structural alignment with our
   framing; if aligned, modernize and integrate. Sociologically
   abandoned in the 1970s.
j. **Jacobson 1995** (gr-qc/9504004) cross-domain extension. Does
   Clausius -> Einstein admit a dipole-pair interpretation? Does it
   extend beyond gravity?
k. **Stochastic Electrodynamics** revisit (Boyer 1975, de la Pena &
   Cetto 1996, arXiv:quant-ph/0501011). Apply 21st-century compute to
   the nonlinear/Coulomb wall that stopped them in the 1980s. Vacuum
   dipole-like fluctuations as substrate may align with our framing.
l. **Sorkin causal sets** (gr-qc/9511063, 0309009). Pick up the
   sequential-growth dynamics program that stalled at Rideout-Sorkin
   1999 (gr-qc/9904062) and try a dipole-pair-coupling reformulation.
m. **Brukner-Zeilinger** one-bit primitive (Found Phys 39, 2009) +
   pair-coupling extension. They had the unit; we add the coupling
   as primitive.
n. **'t Hooft cellular automaton QM** (arXiv:1405.1548). Test whether
   the CA dynamics admits a dipole-pair substrate formulation.
o. **Adler trace dynamics** (hep-th/0206120). Look for Brownian
   correction signatures predicted by trace dynamics in our entropy-
   operator data.
p. **Hardy 2001 / Chiribella-D'Ariano-Perinotti** operational
   QM reconstructions (quant-ph/0101012, arXiv:1011.6451). Build a
   substrate layer below the operational axioms.
q. **Sakharov induced gravity** (gr-qc/0204062 modern review). Modern
   regularization of the cutoff dependence that stopped the 1967
   program.
r. **Wolfram hypergraph rewriting** (arXiv:2004.08210). Provide an
   empirical filter for rule selection by requiring dipole-pair-
   compatibility.
s. **Deutsch-Marletto constructor theory** (arXiv:1405.5563, 1608.02625).
   Add an equation-form substrate layer below their modal substrate.
t. **Penrose twistor** + mass / gravity coupling extension.
u. **Wheeler-Feynman absorber theory** modern closure.
v. **de Broglie double solution** (arXiv:1703.06158) detection program.
w. **Barbour-Bertotti / Shape Dynamics** (arXiv:1010.2481) quantum
   extension.
x. **Penrose-Hameroff Orch-OR substrate-only test** (separated from
   consciousness claim).
y. **Wheeler "it from bit"** formalized into an extraction procedure.

**Per-domain stacks worth running on each domain's native data:**

z. **Chemistry**: bring in **Rao-Esposito 2022** (arXiv:2204.02815)
   info-thermodynamic CRN framework, stack our windowed-H algebraic
   basis on top of their entropy-production-rate analysis. Same on
   **Reinhardt 2019** (arXiv:1904.01988), **Schmitz-Aris 2013**
   (arXiv:1307.7957 -- substitute entropies for concentrations),
   **Sayyadi 2024 SISR** (arXiv:2510.20655).
aa. **Geology**: stack our 6D operator cloud + rank analysis on
    **Garland-Bradley 2018** (arXiv:1811.01272) paleoclimate data
    (extend from single-channel permutation entropy to multi-channel
    algebraic basis). Same on **da Silva 2020 Tsallis seismic**
    (PMC7516945) and **Reinsel 2025 geochemical manifold** (Nat. Sci.
    Rep. 15).
bb. **Biology**: try the asymmetric MI regression on neural,
    ecological, and gene-regulatory time series. Stack on top of
    **Schreiber 2000 TE** (arXiv:nlin/0001042) and **Tishby IB**
    (arXiv:1705.02436) frameworks.
cc. **Stoichiometric subspace** (Feinberg, Horn-Jackson) -- nonlinear-
    algebraic extension on entropies instead of linear on
    concentrations.
dd. **Information geometry of CRNs** (arXiv:2503.19384) -- add
    operator-basis null extraction layer.
ee. **Directed information flow in reaction networks** (bioRxiv 2024,
    doi:10.1101/2024.08.17.608427) -- entropy basis substitution.

Pick from this menu in whatever order Greg's intuition suggests. Or
sample from multiple at once -- the "stacking" itself is the work, and
different stacks may light up different lines.

## Prompt for next session

"Resuming Information Layer work, Session 6.

Branch: claude/linear-drift-nreal-sweep-cjRkM (already on origin).
Session 5 work is committed at tip [TBD after handoff commit].

Before doing anything else:
1. Read CLAUDE.md. Has the v3 ledger plus three Operating Rules from
   Session 4 and Rule D from Session 5 (incomplete not wrong).
2. Read SESSION_HANDOFF_2026-05-26_v5.md (this file). Has Session 5's
   three-thread probe results, what was retracted at interpretation
   level (INFO-014, 018, 019), what was contextualized (INFO-020, 021),
   what new structural finding emerged (INFO-022), and what's queued.

Speaking posture per Rule C applies from your first response forward.
Rule D applies as well: when reading prior findings, default to
'incomplete' before 'wrong.' The data still stands; the reading is
what may be partial.

We are not back at zero. The per-domain algebraic equations from
earlier sessions are still on the table, with reading open. The
Working Frames in CLAUDE.md are unchanged. The Information Layer /
Unified Theory frame is alive.

What I want this session: [Greg's call].

Three queued lines, in priority order if you have no preference:
(i) the 5-physics-areas substrate hunt — me drafting canonical
equation forms for each and identifying the structural elements
that appear in all 5;
(ii) reapply probe discipline to the per-domain algebraic equations
(chemistry quadratic, geology rank-3, biology MI);
(iii) domain-native operator basis for geo (3D) or chem (higher-D),
with the three-thread probe discipline (null inputs, operator-noise
bypass, projection sweep) running alongside from day one."

End of v5 handoff.
