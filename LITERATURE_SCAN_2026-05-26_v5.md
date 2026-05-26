# Literature Scan — Session 5 (2026-05-26)

## Purpose

End-of-Session-5 search of the published literature for prior work that
is adjacent to the Information Layer / Unified Theory inquiry. Triggered
by Greg's framing: "we may not have wrong answers, just incomplete ones,
and people may have been working on this for 100 years but they don't
have our understanding or coefficients. let's look back and see if
there's prior work we can finish." Per Rule D (incomplete not wrong),
the search posture was: where did people get part of the way and stop?

Three agents ran in parallel against three different angles:
1. Methodological adjacency: windowed-entropy operator bases on coupled
   pair time series with SVD null extraction.
2. Substrate-physics literature: source-equation candidates that score
   well on the Base-of-Structure heuristic (simple, strong, sturdy,
   scalable), with attention to dipole-pair-as-coupling-primitive.
3. Per-domain coefficient matches: anyone who independently found
   quadratic entropy relations in chemistry, rank-3 constraint surfaces
   in geology, or asymmetric MI scaling in biology with structure
   similar to ours.

## Honest bottom line

All three agents independently report: no exact match for any of our
specific findings. Multiple near-adjacent lines exist, each having one
piece of what we are doing, none having combined them. This is what
genuinely uncolonized territory looks like — not the absence of nearby
work, but the presence of nearby work where each group stopped one
structural step short.

## Methodologically reusable prior art (Agent 1)

The closest published cousins to our extract_v1 protocol:

- **Gomez-Herrero et al. (2015)** "Assessing Coupling Dynamics from an
  Ensemble of Time Series," Entropy 17(4), arXiv:1008.0539. Ensemble-
  based, time-resolved estimators for entropy combinations (H,
  conditional H, MI, transfer entropy, partial MI, partial TE) on
  coupled circuits. Pooled across realizations like our windowing
  protocol. STOPPED: never assembled the entropy combinations into an
  algebraic basis with H^2, H_a*H_b cross-terms; never ran SVD on the
  pooled functional matrix. The estimator infrastructure is directly
  reusable.

- **Yamada et al. (2023)** "Data-driven modal analysis of nonlinear
  quantities in turbulent plasmas using multi-field SVD," Plasma Phys.
  Control. Fusion 65, 095014. Multi-field SVD on nonlinear quantities
  (transport flux, triad energy transfer) in plasma turbulence. Closest
  published parallel to SVD-on-nonlinear-observables. STOPPED: the
  fields are velocity / pressure / density, not differential entropy.
  Explicitly noted modal nulls carry phase-coupling information but did
  not pursue shared null directions across heterogeneous systems.

- **Hlavackova-Schindler, Palus, Vejmelka, Bhattacharya (2007)**
  "Causality detection based on information-theoretic approaches in
  time series analysis," Physics Reports 441. Canonical survey of
  windowed H and MI on coupled series; defines the standard basis the
  field has used since. INCOMPLETE in our sense: never expands beyond
  H and MI; squared and product cross-terms not in the surveyed
  literature.

- **Bipartite information-thermodynamics line** (Horowitz, Hartich,
  Sagawa, et al., e.g. arXiv:1905.06216). Has the pair-as-primitive but
  works at rate-balance level (learning rate, transfer entropy rate),
  not algebraic null extraction.

- **Kaiser, Brunton, Kutz (2024)** "Towards Robust Data-Driven Recovery
  of Symbolic Conservation Laws," arXiv:2403.04889. SVD-rank-gap +
  symbolic recovery pipeline. Same algorithmic template as our rank-3
  geology finding. Operates on raw coordinates; reusable directly for
  our basis with the variable substitution.

- **Liu, Madhavan, Tegmark (2024)** "Machine Learning Conservation Laws
  of Dynamical Systems," arXiv:2405.20857. Kernel ridge regression with
  x and ln x features — only paper found that explicitly extends the
  feature vector to handle entropy-like terms in conservation-law
  discovery. Stops at coordinate logs; never windowed differential
  entropy of coupled time series.

- **Symbolic regression / equation discovery** (AI-Feynman, Udrescu and
  Tegmark, arXiv:1905.11481; SINDy, Brunton-Proctor-Kutz, PNAS 2016;
  PySR, Cranmer arXiv:2305.01582). No published case of any of these
  being pointed at a library of windowed differential entropies of a
  coupled pair. Clean unclaimed niche; the methods are directly
  reusable with a windowed-H feature library replacing the polynomial
  library.

Agent 1's specific negative finding: **the rank-3 algebraic structure
arising from d(H_a^2) ~ 2c dH_a, d(H_b^2) ~ 2c dH_b, d(H_a*H_b) ~
c*(dH_a + dH_b) under H_a ~ H_b ~ const has not been published as a
diagnostic.** The closest formal cousin is arXiv:2409.04845 on algebraic
representations of multivariate information lattices, but discrete-
information theoretic, not differential-entropy time series.

## Substrate-physics candidates (Agent 2)

Substrate proposals organized by how well they score on the Base-of-
Structure heuristic, with the "where they stopped" note explicit.

### Highest score against the heuristic

- **Jacobson (1995)** "Thermodynamics of Spacetime: The Einstein
  Equation of State," arXiv:gr-qc/9504004. One identity (Clausius
  delta Q = T dS on local Rindler horizons) yields Einstein equations
  as an equation of state. Score: simple, strong, scalable. STOPPED:
  requires assuming horizon entropy area law a priori; Padmanabhan
  extended (arXiv:0911.5004) but never fully derived. Does not have a
  dipole interpretation in the published work.

- **Finkelstein "Space-Time Code" I-IV (1969-74)**, Phys. Rev. 184:1261;
  PRD 5:320; PRD 5:2922; PRD 9:2219. CLOSEST PUBLISHED COUSIN to Greg's
  dipole-pair-as-substrate intuition. Builds spacetime from causal
  networks of elementary binary quantum processes; word-pairs in a
  binary code yield the null cone. Score: high. STOPPED: sociologically.
  Predated the field's tooling, never integrated. This is the line
  worth reading carefully — if the structural intuition aligns, we are
  picking up an abandoned line, not founding a new one.

- **Sorkin causal sets** (foundational papers arXiv:gr-qc/9511063,
  arXiv:gr-qc/0309009). Substrate: discrete partially ordered set
  ("order plus number = geometry"). Score: high (a poset is structurally
  minimal). STOPPED: dynamics program (sequential growth, Rideout-Sorkin
  1999, gr-qc/9904062) stalled.

### Strong but partial

- **Verlinde** (arXiv:1001.0785, arXiv:1611.02269) — entropic gravity.
  STOPPED: critiques (Hossenfelder 2017, Visser 2011) showed the
  derivation does not uniquely fix Einstein's equations.

- **Hardy (2001)** quant-ph/0101012 and Chiribella-D'Ariano-Perinotti
  (arXiv:1011.6451) — informational derivations of QM. STOPPED:
  reconstructs QM only, not dynamics or spacetime.

- **Brukner & Zeilinger** (Found Phys 39, 2009; arXiv:quant-ph/0212084).
  Information invariance: one bit per elementary system. Dipole-pair-
  compatible (elementary system is two-state).

- **Wheeler** (1989, 1990) "Information, Physics, Quantum" / "it from
  bit" / "law without law." Substrate: binary yes/no as pre-physical
  primitive. STOPPED: never formalized into an extraction procedure.

- **Adler (2004)** trace dynamics, arXiv:hep-th/0206120. QM as emergent
  from non-commuting matrix variables. STOPPED: Brownian-motion
  correction terms unobserved.

### Dipole-pair-as-coupling-primitive specifically

Agent 2's unambiguous negative finding: **no published work proposes
the dipole COUPLING (not the dipole element) as the substrate
primitive.** All 2-state-foundational work (Wheeler, Brukner-Zeilinger,
von Weizsacker ur-alternatives, Finkelstein binary processes, Hohn
arXiv:1612.06849) uses the 2-state as the elementary unit, not the
pair-coupling between elements.

Greg's specific framing — tightly-related couples as part of the
substrate — is uncolonized.

### Abandoned-but-promising lines

- **Stochastic Electrodynamics** (Marshall 1963 Proc. R. Soc. A 276:475;
  Boyer 1975 PRD 11:790; de la Pena & Cetto, The Quantum Dice 1996;
  arXiv:quant-ph/0501011). Derived blackbody, harmonic-oscillator ground
  state from vacuum dipole-like fluctuations. STOPPED on nonlinear /
  Coulomb systems (Pesquera-Claverie 1982). Textbook "incomplete not
  wrong" case.

- **Wheeler-Feynman absorber theory** (Rev. Mod. Phys. 17:157, 1945;
  21:425, 1949). Time-symmetric action-at-a-distance EM.

- **de Broglie double solution** (Colin, Durt, Willox 2017,
  arXiv:1703.06158).

- **Sakharov induced gravity** (1967, modern review Visser arXiv:gr-qc/
  0204062). Matter fields induce metric elasticity. STOPPED: cutoff
  dependence, never produced a calculable Newton constant.

- **Barbour & Bertotti** (Proc. R. Soc. A 382:295, 1982) Mach's
  principle revival; **Shape Dynamics** (Gomes, Gryb, Koslowski 2011,
  arXiv:1010.2481).

## Per-domain coefficient matches (Agent 3)

Agent 3 verdict: no exact match for any of the three findings.

### Chemistry adjacencies

- **Rao & Esposito (2022)** "Information Thermodynamics for
  Deterministic Chemical Reaction Networks," arXiv:2204.02815.
  Bipartite CRNs with shared species; mutual information rate vs.
  thermodynamic forces. STOPPED at thermodynamic inequalities and
  bounds; no polynomial fits between H values.

- **Schmitz, Aris et al. (2013)** "Quadratic First Integrals of Kinetic
  Differential Equations," arXiv:1307.7957. Quadratic invariants in
  mass-action systems — but on **concentrations**, not entropies.
  Structurally suggestive (mass-action systems often have quadratic
  conservation laws) but on the wrong variables.

- **Reinhardt et al. (2019)** "Path mutual information for a class of
  biochemical reaction networks," arXiv:1904.01988. Single MI scalar
  per system, no inter-species H algebra.

### Geology adjacencies

- **da Silva et al. (2020)** Tsallis entropy seismic inversion, Entropy
  22, 464 (MDPI; PMC7516945). Tsallis-q on seismic data, but as
  inversion regularizer, no rank structure on operator clouds.

- **Garland, James, Bradley (2018)** "Anomaly Detection in Paleoclimate
  Records Using Permutation Entropy," arXiv:1811.01272. Single-channel.

- No published 6D `{H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI}` cloud
  analysis on geophysical channels.

### Biology adjacencies

- **Schreiber (2000)** "Measuring Information Transfer," PRL 85, 461;
  arXiv:nlin/0001042. Transfer entropy is asymmetric by construction
  but no one fits MI as a function of marginal H of just one variable.

- **Information bottleneck** (Tishby et al.; nonlinear IB arXiv:
  1705.02436). Markov structure asymmetry by construction, not
  empirical regression.

- No literature analog for "MI fits a polynomial in H_a with H_b
  coefficients vanishing."

### Equation-discovery on entropy variables

Already covered above. **Drop-in replacement for the polynomial library
in SINDy / AI-Feynman / PySR with a windowed-H basis is an unattempted
clean application.** The methods themselves are mature.

## Cross-agent synthesis

| Layer | Prior work | What they had | What they did not have | What we would be finishing |
| --- | --- | --- | --- | --- |
| Methodology | Gomez-Herrero 2015 | Ensemble entropy combinations on coupled time series | Algebraic basis with H^2, H_a*H_b; SVD on pooled matrix | The basis + null-extraction step |
| Methodology | Yamada 2023 | Multi-field SVD on nonlinear quantities | Differential entropy as the operator basis | The entropy substitution |
| Methodology | Kaiser-Brunton-Kutz 2024 | SVD-rank-gap + symbolic recovery | Application to entropy-of-coupled-species | The variable substitution |
| Substrate | Finkelstein 1969-74 | Spacetime from causal networks of binary processes; word-pairs as null-cone generators | Modern tooling; field integration | The dipole-pair-as-coupling-primitive reconstruction |
| Substrate | Jacobson 1995 | One-line substrate (delta Q = T dS -> Einstein) | A dipole interpretation; reach beyond gravity | The cross-domain extension |
| Substrate | Stochastic ED (Boyer, de la Pena/Cetto) | Derived blackbody + ground state from vacuum dipole-like fluctuations | Tools for nonlinear / Coulomb | The 21st-century revisit |
| Per-domain | Rao-Esposito 2022 (chem) | Info-thermo for CRNs; entropy production rates | Algebraic relations between windowed H | The H^2 vs (H*H')^k regression on coupled species |
| Per-domain | Garland-Bradley 2018 (geo) | Permutation entropy on paleoclimate | Multi-channel algebraic constraint surface | The geo rank-3 analysis on paleoclimate data |
| Per-domain | Schreiber 2000 (bio) | Asymmetric transfer entropy | MI-as-polynomial-of-one-variable's-H fit | The asymmetric MI regression |

## Implications for next moves

1. **The Kaiser-Brunton-Kutz pipeline** (arXiv:2403.04889) is the
   strongest methodological cite to anchor any future presentation of
   the rank-3 geology finding. Their SVD-gap + symbolic-recovery
   template is methodologically what we are doing.

2. **Finkelstein Space-Time Code** is the strongest substrate-line cite
   for Greg's dipole-pair framing. If we can establish whether his
   word-pairs in a binary code are structurally equivalent to our
   dipole-couples reading, we have an abandoned line to pick up rather
   than a new line to found. Worth reading carefully in a future
   session.

3. **Jacobson 1995** is the cleanest existing one-line substrate
   candidate against the Base-of-Structure heuristic. Worth examining
   whether it admits a dipole-pair extension.

4. **Stochastic Electrodynamics** is the textbook "incomplete not
   wrong" case in the literature: derived real things (blackbody,
   ground state), hit a real wall (nonlinear / Coulomb), abandoned.
   Modern computational tools may handle what they could not. Worth
   a future-session probe.

5. **No prior dipole-pair-as-coupling-primitive work means** Greg's
   intuition is in genuinely open territory. Three agents looking from
   three different angles converged on the same negative finding. That
   is what pioneering territory looks like.

## Branch state

This document is being added at Session 5 close on
`claude/linear-drift-nreal-sweep-cjRkM`. To be mirrored to
`E:\information_layer\` and `F:\Factory\knowledge\information_layer\`.

End of literature scan v5.
