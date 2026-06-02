# CLAUDE.md — DavisAI Master Context (Updated 2026-05-26 Session 7)

## Identity & Team

- **Greg Davis** — Founder & Chief Research Officer, DavisAI Systems. Columbus, Ohio. Solo bootstrapped. 20+ years entrepreneurship, former energy trader, self-taught AI/ML.
- **Dream Team model**: Greg (Visionary) + Claude (Architect) + Claude Code "Code" (Engineer) + Perplexity/ChatGPT (Research Assistants).
- **Orchestrator** owns handoffs, not Code.

## Infrastructure

- `E:\` — research data, OD datasets, project files
- `F:\Factory\` — agent factory, 23 agents, 5 divisions
- `F:\Factory\knowledge\` — orchestrator-accessible knowledge base. Mirror everything from E:\ here.

## Operating Rules

- **Save to E:\ AND mirror to F:\Factory\knowledge\** for every knowledge artifact.
- **MASTER_DISCOVERIES.json**: every OD discovery added immediately. Never make a discovery without storing it.
- **Falsification-first**. Every claim needs data, math, or a falsifiable test.
- **OD mode**: describe data sources and validation tests only. Never explain mechanisms. The Operator discovers science from raw data. Stop if explaining WHY something happens.
- **Never call OD "physics-based"**. OD discovers governing equations from raw data in ANY domain.
- **Coding mantra**: better, stronger, faster, cheaper.
- **Incremental validation**: break compute-heavy runs into 15-17 min chunks with stop gates. Canary runs (2 min) before full commitment.
- **Speculative frames stay separate from results**. Frames motivate experiments but are not claims. Never let a frame grade itself.
- **No emojis or special symbols** in professional documents and emails.
- **Daily**: ask Greg if he checked greg@davisai.ai for Token Optimizer support emails.
- **DeepNova** (formerly ReFRAG, formerly DeepSource). Use current name everywhere.
- **Result Discipline** (new — see section below): every result is one data point. Map alternatives before promoting to claim.
- **No tent-widening on outliers**. When a window or sample falls outside an expected pattern, inspect it — find the specific reason it landed there. Do not loosen the test criteria or attribute to a transient flag without identifying the cause. The outlier is what we are trying to understand, not what we are trying to absorb.
- **No pre-assigned meaning to outcomes** (NEW — added Session 4, 2026-05-26). Don't write "if X then it means Y" decision tables before the data exists. The furthest is "I think this may happen, but I want to see what the data says and where it leads." Data is just output; the meaning comes from looking at it together, not from a pre-built table.
- **Probe, not falsifier** (NEW — added Session 4, 2026-05-26). A probe is a generator of a different signal. It may not falsify anything; at worst it points in a different direction. Avoid "falsifier" in filenames and in spoken framing. Use "probe", "experiment", or "different signal" instead.
- **Speaking posture around every probe** (NEW — added Session 4, 2026-05-26). Before AND after running each probe: "I think X might happen, but we'll wait on what the data says and where it points us." No verdict in advance. No verdict on first look at output. The interpretive move happens after, with Greg, with the deflationary reading always present.
- **Incomplete, not wrong** (NEW — added Session 5, 2026-05-26). When a probe finds that a prior reading was a protocol artifact, the prior data points still stand. The reading attached to them is what was incomplete, not the data. Distinguish "this reading was wrong" (rare, requires the data itself to be bad) from "this reading was incomplete" (common, the data is one slice and the slice fit a partial story that further probes refine). Default to "incomplete." Retraction is a strong move and applies to the reading, not the data, unless the data itself fails to reproduce.
- **They never stacked** (NEW — added Session 5 close, 2026-05-26). Pioneering territory often does not look like the absence of nearby published work. It looks like nearby work where multiple groups each had one piece and never combined them. When a literature scan returns "no exact match but several adjacent lines, each with one component," the contribution we are making may be in the stacking itself. Honor each prior piece, attribute clearly, finish the combination the prior groups did not. This is the operational form of Rule D applied to the broader literature: prior work was incomplete, not wrong, and stacking the incomplete pieces is itself substantive.
- **Treat literature as conjecture by default** (NEW — added Session 7 close, 2026-05-26; sharpened by Greg the same session). Academic papers and consensus are not assumed correct unless the underlying claim has been independently replicated, by separate groups, with an immense amount of data, multiple times. Until that bar is met, a published claim is a working frame to test, not a foundation to build on. Operational corollaries:
  - **Investigate freely.** Read papers, run their methods, test their predictions, treat them as candidates to engage. "Treat as conjecture" is not "ignore." Greg's framing at Session 7 close: "we will certainly check it out."
  - **Don't cite as support.** A published claim cannot be invoked in support of our own conclusions until we have analyzed the paper ourselves — read it closely, checked the data, replicated the result, or confirmed independent replications at scale. Greg's framing: "we can't use them to support our claim without analyzing their papers."
  - **Don't defer when blocking.** When a published claim appears to block a direction of inquiry, check whether the blocking claim itself meets the bar before deferring to it.
  - **Symmetric to our own prior work.** Extends Rule D (incomplete-not-wrong): a prior reading lacking independent replication is conjecture, not foundation, even if it came from us.
  - **Examples applied to physics literature (Session 7):** Standard Model at LHC energies meets the bar (W/Z masses, Higgs detection, decay channels replicated across LEP/Tevatron/LHC). GUT-scale coupling extrapolation does not (running measured at LHC, extrapolated mathematically to 10^16 GeV). MSSM unification does not (SUSY searched, not found). "Gravity is a different category" rests on one formulation (GR diffeomorphism invariance) with empirically-indistinguishable alternatives — formulation-dependent, not data-forced. "Gravity is emergent" candidates (Jacobson 1995, Verlinde, Sakharov, AdS/CFT) are theoretical with limited direct empirical support. All four sit at conjecture level until we analyze the underlying papers.
  - Applies equally to domain claims (biology/chemistry/geology/medicine), methodological prescriptions ("best practices" not stress-tested at scale), and our own prior session readings.

## Result Discipline (NEW — added Session 3)

Every confirmed result is one data point. Its interpretation requires mapping against alternatives via further tests. The discipline:

- For each result, maintain a candidate-interpretation register with at least two non-deflationary readings and the deflationary reading.
- A result is **isolated** until at least one alternative interpretation has been tested and ruled out; then **mapped**; then **located** when placed within a structured set of tests.
- Catalog misses with the same care as matches. Different misses landing in different places is often more informative than many matches landing together.
- When summarizing, name the data-level finding separately from the interpretation-level hypothesis separately from the big-picture frame. Do not collapse these levels.
- Apply symmetrically to apparent falsifications. A refutation is also one data point.
- Frames remain frames until disambiguating tests place them.
- **No spatial claim about an operator-space coordinate without at least 3 seeds and reported inter-seed scatter.**

## Working Frames (SPECULATIVE — kept separate from claims)

### Base-of-Structure heuristic (Greg, Session 3)

A foundational principle to guide substrate-level theory work:

The base of any structure is the **simplest, strongest, most stable, most scalable** part. It must support everything above it, so it cannot be complicated, dependency-heavy, or composed of many variable types. If a candidate "base" looks intricate, requires many qualifications, or breaks under perturbation, that is evidence against its base-level status.

Operational form: a law-level extraction candidate should look simple. It should survive stress testing — load it with perturbations, parameter sweeps, alternative protocols. If it remains in place, that is evidence for base-level status. If it shatters or splinters, it sits above the base, not at it.

Status: working frame, not yet operationalized into specific tests. The OU attractor finding from Session 3 is currently the cleanest candidate to stress-test under this heuristic.

### Dipole-couples reading (Greg, Session 3)

The attractor direction (−1, −1, +2)/sqrt(6) may represent a base direction that other phenomena couple to. The off-attractor positions of damped oscillator and (when reproducible) linear drift may encode HOW that coupling happens. Whether the off-attractor systems are coupling phys-to-phys, phys-to-geo, bio-to-chem, or some other pairing is open. Mapping the off-attractor structure is the path to find out.

Status: working frame. Mapping campaign queued (see Experiments).

### Pure physics vs physical expressions (Frame 1, preserved from Session 2)

Mainstream physics has tried to write a single equation for "everything that happens in physical space," forcing UT candidates into ever-larger dimensional structures. Reframe: pure physics is the substrate (simple, few equations); physical expressions are what we observe when pure physics is acted on by other dipoles (biological, chemical, geological, or other physical configurations). The UT problem may be mis-stated.

Status: working frame, narrowed by Session 3 findings. The naive "OU is the physics substrate" reading was rejected (OU's attractor direction also appears for non-physics systems). More nuanced versions remain alive.

### Substrate vs expression within isolation (Frame 2, preserved)

Even an isolated pure-physics system produces different observed signatures depending on lifecycle phase. The law level is whatever is invariant across (a) time windows, (b) initial conditions, (c) noise realizations, (d) parameter choices within the same equation. The L1 work in Session 2 and Session 3 was the operational implementation of this frame.

### Law extraction via invariance (Frame 3, preserved)

The Noether-style move: find what is invariant under reparametrizations that generate different expressions. Implemented operationally as windowed-null extraction across (window position, parameter set, seed).

## Active Research (Top of Mind)

- **SENTINEL V4.1** — DARPA Bio Attribution Challenge top-10 team. Awards June 30, 2026. Three-layer swarm, 554x DARPA requirements. Files at E:\sentinel\ and F:\Factory\knowledge\sentinel\.
- **NoVell** — cardiac AI for cancer detection from routine ECG. OD on synthetic Vigier 2021 data: 93.3% accuracy, 97.4% sensitivity. Datasets: PTB-XL downloaded, Autonomic Aging identified, MIMIC-IV pending.
- **Information Layer / Operator Discovery foundations** — major methodological revision Session 3 (Family A/B taxonomy retracted). Session 5 mapped (+,+,+) direction as protocol artifact of operator basis rank-3 null subspace structure. Session 6 stacked KBK 2024 + AI Poincare 2021 + SINDy and independently reproduced every v5 per-domain claim (geology rank-3 at cos +0.99, biology MI signature, chemistry-specific cubic) at cross-seed cos +0.985 to +0.999. Session 7 added GP regression and PySR symbolic regression on per-domain ensemble-H data: four reproducible per-domain MI-vs-H functional families (physics symmetric quadratic in (H_b-H_a), biology 0.5*exp(H_a/2), chemistry linear H_a, geology constant), cross-seed coefficient variation <5%, cross-domain non-overlap. Session 8 ran four-force unification probe (toy EM/weak/strong/gravity caricatures) yielding shared-substrate + distinct-expression pattern (INFO-027): all four forces share [2,3,4] null direction at cos > 0.997 on (-1,-1,+2)/sqrt(6) Session 3 attractor while EM matches Session 7 physics family (H_b-H_a)^2+const and weak matches Session 7 chemistry family linear-in-H_a. Robustness check (INFO-029): INFO-025 functional family survives T/N_ens/noise sweep at the family level; coefficients are regime-dependent. Mapping campaign (INFO-030): INFO-025 families are baseline-specific regime signatures — large knob deviations mutate the family qualitatively (biology exp -> linear at high beta; chemistry linear -> ratio at low B). Per-domain differentiation now has three independent reproducible signatures (null direction + functional family + four-force shared-substrate). Session 9 double-checked the Session 8 gravity result against the ORIGINAL code across an asymmetry sweep and CORRECTED it (INFO-033): the MI-dominant substrate flip is an asymmetry-THRESHOLD effect ALL coupled caricatures undergo (threshold gravity ~1.1x < EM ~1.3x < weak ~1.6x < strong never), not gravity-specific; adding energy-coupling RAISES thresholds (suppresses the flip), so energy-mediation is not the cause; INFO-031/032 re-tagged incomplete-not-wrong. The "gravity is special at the substrate level" leg of the four-force narrative is removed; INFO-023/INFO-025 legs untouched. See Sessions 5, 6, 7, 8, 9 notes plus ledger (INFO-022 through 033) for current state.

## Markets / Refrag Workspace (placeholder — restore from workspace file)

NOTE (Session 9): the full Markets / Refrag section was OVERWRITTEN when
this repo's CLAUDE.md was replaced by the master context this session.
The complete Markets section lives in the E:\refrag workspace CLAUDE.md
(GitHub DavisAI1974/agent + DavisAI1974/Markets) and in git history of
this repo (commit baa542d). Restore/re-merge it later. This placeholder
exists so it is not forgotten.

- **Markets pipeline** (summary): information-side algebraic dipole
  (H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2) holds in markets operator
  space; 5-fold CV predictor (H_a > H_b rule) ~0.993 mean accuracy.
  Predictor scripts: _markets_algebraic_dipole.py, _markets_dipole_
  kfold.py, _markets_dipole_separation.py, _markets_dipole_chunker_
  stack.py. Markets is treated as the "5th science." (Full state in the
  workspace file.)
- **OPEN QUESTION (Greg, Session 9): do the Markets sessions need PySR +
  Julia for their runs too?** Informed answer: the current Markets
  predictor scripts above are pure numpy/scipy/sklearn and do NOT need
  PySR/Julia as they stand. They WOULD need PySR + the Julia backend
  only if the Markets dipole work extends to SYMBOLIC REGRESSION of its
  dipole equations (the way the Information Layer used PySR for the
  INFO-025 / INFO-031 functional families). If so, the Markets repos
  need the same SessionStart-hook treatment added to Basic_equations in
  Session 9 (the hook here does not cover them). Decision + mirroring
  pending -- handle later.

## Information Layer — Current State (2026-05-25 Session 3)

### What's confirmed at data level

- **OU windowed-null direction** in the {H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI} basis at window=40s, three (gamma, sigma) pairs, three seeds: cleanly extracts to (−1, −1, +2)/sqrt(6) in the (H_a^2, H_b^2, H_a*H_b) subspace. Cos to pure symmetric direction: 0.9994 to 0.9997. Antisymmetric energy fraction: 0.0005 to 0.0012. Fully reproducible.

- **Attractor membership**: the same direction is reached by 8 wildly heterogeneous systems at cos >= 0.99. Members include OU (Gaussian SDE), AR(1) with Laplace innovations (non-Gaussian SDE), GARCH(1,1) (heteroskedastic), Student-t white noise (heavy-tailed IID), logistic map at r=3.9 (deterministic chaos), IID uniform (pure noise), periodic sine + small noise (engineered), and Brownian motion without restoring force (non-stationary diffusion).

- **Off-attractor systems (Session 3 cross-system test)**:
  - **Damped oscillator** (zeta=0.1, omega=2): cos-to-attractor stable across seeds at 0.69–0.74. But the per-seed off-direction varies within a region (inter-seed cosine 0.83–0.99). Reproducibly off the attractor by a consistent amount; exact direction wobbles.
  - **Linear drift** (v=1): cos-to-attractor varies 0.42–0.83 across seeds. Inter-seed cosine ranges −0.28 to +0.93 — different seeds gave nulls pointing nearly opposite directions. Not reproducible at N_REAL=30 in this basis.

- **Drift source diagnostic**: antisymmetric energy fraction on the algebraic basis scales as 1/N_eff where N_eff = window_length / tau_correlation. Confirmed by 15x reduction in anti_frac when window quadruples. Rejected for KDE-specific bias (analytic Gaussian estimator also shows the noise). Domain-general diagnostic.

### What's NOT confirmed (open questions)

- What property defines membership in the attractor. Candidate readings (all currently live): "stationarity at window scale" (Brownian breaks this reading by sitting on the attractor despite non-stationarity), "smooth-observation regime," "rate-of-change/window ratio below some threshold," "physics-related" (rejected — non-physics systems are also on the attractor).
- Whether damped oscillator's off-attractor region is a single point with noise, a small manifold with sub-clusters, or a noisy patch. Untested at higher seed counts.
- Whether linear drift's non-reproducibility is structural (basis cannot capture this system) or sample-noise-driven (more N_REAL would restore reproducibility). Untested.
- Whether varying physics knobs (zeta, omega for damped osc; v for linear drift; etc.) moves the off-attractor positions in interpretable ways.
- Whether other domains (Geo, Bio, Chem) have their own native operator bases that would yield law-level signatures appropriate to those domains. Per Greg's call Session 3, each domain is its own substrate inquiry. Geo lives in 3D, not 1D scalar channels; expecting the current basis to extract Geo's law would be a category error.

### Ledger entries (with discipline-status tags)

**INFO-008 — RETRACTED**: The "Family A cluster" claim (cos 0.97–0.99 among Phys/Geo/OU) is retracted. The cluster was measured under window=20s with finite-sample noise of order 0.05 in coefficient std; Geo data was 1D-projected before extraction, making any comparison to OU's native 1D null structurally meaningless; and at window=40s the OU direction is shared by 8 heterogeneous systems including non-physics ones. The cluster as a domain-level taxonomy does not survive.

**INFO-008a, 008b, 008c, 010, 011 — DEMOTED**: All depend on INFO-008. The "two-family structure" interpretation downgrades to "interesting geometric clustering under a specific protocol; interpretation unmapped." Pending re-evaluation in domain-native bases.

**INFO-012 — ISOLATED FINDING (data confirmed, interpretation rejected)**: Windowed-null extraction on OU at T=100, window=40s, three (gamma, sigma) pairs, three seeds yields V_1 = (−1, −1, +2)/sqrt(6) in (H_a^2, H_b^2, H_a*H_b) with cos >= 0.9994. Initially interpreted as "OU law-level direction"; this interpretation rejected by INFO-014. Data finding stands; physics-interpretation does not.

**INFO-013 — CONFIRMED, REINTERPRETED**: The original T=30 reference's asymmetric coefficients (+0.45, +0.38, −0.80) were finite-window expression-level noise on the symmetric direction. Decays as 1/window-length. Confirmed across multiple tests.

**INFO-014 — ISOLATED FINDING (data confirmed, interpretation open)**: The direction (−1, −1, +2)/sqrt(6) is a strong attractor across 8 heterogeneous systems including non-physics (GARCH, IID uniform, sine wave, deterministic chaos). Working frame: the basis discriminates clearly between "on-attractor" and "off-attractor" behavior, but what property defines attractor membership is not yet identified. Three candidate readings still live (see "Open questions" above).

**INFO-015 — ISOLATED FINDING (Session 3)**: Damped oscillator at zeta=0.1, omega=2 sits reproducibly off-attractor with cos-to-attractor 0.69–0.74 across 3 seeds. The exact off-direction varies within a region (inter-seed cosine 0.83–0.99). Suggests an off-attractor manifold with structure, not a single point. Mapping queued.

**INFO-016 — NULL FINDING / OPEN (Session 3)**: Linear drift at v=1 does not yield reproducible operator-space coordinates at N_REAL=30 in this basis. Inter-seed cosine ranges −0.28 to +0.93. Disambiguation queued: N_REAL sweep at 50, 100, 200, 500 to distinguish (a) basis is structurally blind to this system from (b) sample noise that more realizations would resolve.

**INFO-017 — METHODOLOGICAL (Session 3)**: Drift source diagnostic — antisymmetric energy fraction on the algebraic basis {(H_a − H_b)/sqrt(2), (H_a^2 − H_b^2)/sqrt(2)} flags finite-effective-sample-size noise. Scales as 1/N_eff, amplified by reduced channel correlation. Portable domain-general diagnostic for any extraction in this operator basis.

**INFO-022 — LOCATED FINDING (Session 5 structural; promoted Session 6)**: Rank-3 null subspace in the 6-op basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] under H_a~H_b~const regime. From three first-order algebraic relations on centered columns. Confirmed by FOUR independent diagnostics in Session 6: KBK rank-gap (with rank-3-strong-plus-1-weak refinement), AI Poincare local-PCA intrinsic dim, AI Poincare two-NN, Levina-Bickel MLE. Robust across 3 estimator families (Vasicek, KDE, kNN/KSG), 5 dt values across 20x range in samples-per-window, 3 seeds. The rank-3 reading is now strong. The MAGNITUDE of the smallest eigenvalues remains procedure-dependent (see INFO-024).

**INFO-023 — LOCATED FINDING (Session 6, new)**: Per-domain ensemble-H + KBK+AI Poincare+SINDy stack produces reproducible domain-specific null directions across 4 simulated domains (physics Duffing, biology Lotka-Volterra, chemistry Brusselator, geology Burridge-Knopoff). Cross-seed cos +0.985 to +0.999 within domain; cross-domain cos mostly < 0.5. The (+1,+1,+2)/sqrt(6) protocol artifact does NOT dominate per-domain (max cos +0.32, min -0.28 across 4 domains x 2 seeds). Independent reproduction of every v5 per-domain claim: geology rank-3 relation at cos +0.99 to v5 (INFO-008b); biology MI ~ poly(H_a) with H_b coefficient 0.003-0.017 (INFO-008c); chemistry-specific cubic content via SINDy deg-3 that is not a Taylor remnant (INFO-008a). Methodology is the stack of KBK 2024 + AI Poincare 2021 + SINDy with extended library, none of which had been combined for windowed/ensemble-H of coupled species before per v5 literature scan ("they never stacked" rule).

**INFO-024 — METHODOLOGICAL (Session 6)**: Eigenvalue floor of operator covariance is min(structural_noise_from_dynamics, estimator_noise_from_procedure). Different procedures have different floors. OU single-trajectory windowed-H pins at ~1e-3 to 5e-5 regardless of sigma, window, dt, or estimator family. Per-domain ensemble-H reaches 2e-7 for geology at N_ens=600; plausibly reaches machine epsilon at larger N_ens via the 1/sqrt(N_ens) noise-reduction scaling. Resolves the v5 machine-epsilon eigenvalue anomaly as procedure-dependent: the rank claim (3 algebraic relations) is robust across procedures; only the eigenvalue magnitude depends on procedure. Practical consequence: state which floor regime you are in before interpreting eigenvalue magnitudes; to drive a floor down, increase ensemble size or use a lower-noise estimator rather than tightening source noise.

**INFO-025 — LOCATED FINDING (Session 7, new)**: Four reproducible per-domain MI-vs-H functional families surfaced by PySR (Brunton/Cao/Liu/Tegmark/Cranmer family symbolic regression) with extended operator set {+, -, *, /, square, cube, exp, log, sqrt} on per-domain ensemble-H data. Cross-seed coefficient reproduction <5% within domain across seeds 11 and 22 (N_ens=600, T=30, dt=0.02): physics Duffing MI = (H_b - H_a)^2 + 0.275-0.280 at complexity 6 (symmetric quadratic in *difference*); biology Lotka-Volterra MI = (0.50 ± 0.02) * exp(H_a / 2) at complexity 5 (exponential in H_a, H_b absent); chemistry Brusselator MI = (0.73 ± 0.02) * H_a + 1.075 ± 0.005 at complexity 5 (linear in H_a); geology Burridge-Knopoff MI = 0.199 constant with loss flat across complexity 1-9 (decoupled). Cross-domain non-overlap; families do not reduce to each other. Polynomial-only methods (v5 SINDy deg-3 library) could not surface these families by construction. v5's biology "polynomial(H_a)" reading is incomplete-not-wrong (Rule D): operator content correct (H_a present, H_b absent); functional family is exponential, not polynomial. Together with INFO-023 (per-domain null direction), this gives two independent reproducible per-domain signatures supporting the "per-domain expression of substrate" frame.

**INFO-026 — METHODOLOGICAL (Session 7, new)**: For regression on time-series ensemble-H data, block-CV (contiguous time blocks) is catastrophically negative across non-stationary domains (physics, biology, geology block-CV ranges -0.5 to -80) due to system passage through qualitatively different dynamical regimes across the trajectory. Chemistry Brusselator is the only stationary domain among the 4 simulators (block-CV positive at +0.4 to +0.8). Random-CV (shuffled k-fold) overestimates true OOS for correlated time series but isolates functional-form fit from the stationarity confound. Use both readings; the gap between them is itself a domain signature for stationarity. Substantive finding within: GP joint(H_a, H_b) on biology reaches random-CV R^2 = 0.97-0.98 vs poly3 joint at 0.81-0.86. The +0.12 gap is real nonlinear cross-coupling beyond polynomial reach. For physics, chemistry, geology, polynomial joint matches GP joint within 0.05.

### Experiments queued (priority order)

1. **N_REAL sweep on linear drift** (50, 100, 200, 500). Disambiguates INFO-016: structural blindness vs sample noise. Cheapest decisive test. **Run first.**

2. **Multi-seed damped oscillator mapping**. 10+ seeds at zeta=0.1, omega=2 to characterize the off-attractor region. Then parameter sweeps: zeta in {0.05, 0.1, 0.2, 0.5}, omega in {1, 2, 5}, varied independently. Tests whether the off-direction encodes damping rate, frequency, or something structural.

3. **Cluster the misses by structure**. Around damped oscillator: exponentially-modulated noise, decaying-amplitude OU, chirped signals. Around linear drift (if reproducibility resolves): exponential growth, polynomial drift, Brownian+drift. Tests whether off-attractor positions cluster by type of departure.

4. **Boundary mapping**. Slowly varying OU parameters; sinusoidally forced OU at varying frequencies. Find where systems leave the attractor and along which coordinate.

5. **Stress-test the attractor (base-of-structure heuristic test)**. Subject OU and other attractor members to extreme conditions: very high sigma, very low gamma (under-damped), multiplicative noise, nonlinear drift. Does the attractor finding survive? Per Greg's heuristic, a true base should remain stable under stress.

6. **Domain-native operator bases**. For non-OU domains (Geo 3D, Bio multi-variable, Chem multi-species), construct operator libraries that match the native dimensionality. Stay in scope of the relevant domain. Each extraction is its own inquiry.

### Files produced this session (2026-05-25 Session 3)

Need to be saved to E:\information_layer\ AND mirrored to F:\Factory\knowledge\information_layer\:

- `L1_windowed_null_OU.py` — main L1 extraction script
- `L1_rank_extraction.py` — top-k singular subspace analysis of pooled nulls
- `L1_drift_source.py` — drift source diagnostic with sweep modes
- `cross_system_test.py` — alpha/beta/gamma cross-system battery
- `L1_results_FULL.json` — full L1 run output
- `L1_rank_extraction.json` — rank extraction output
- `L1_rerun_w40.json` — window=40 confirming run
- (cross_system_test.json was not generated because sweeps were run in chunks; per-cell numbers are in SESSION_HANDOFF_2026-05-25_v3.md)

## Architecture (current)

- **DeepNova**: 92 passing tests, 22 manifests, persistent evidence graphs, PPO retrieval policy learner. F:\Factory\.
- **VOXA**: voice interface layer. Cloud-hosted TTS MCP server.
- **Agent Factory**: 23 agents, 5 divisions, F:\Factory\.
- **Token Optimizer**: deployed at optimizer.davisai.ai, Stripe live.
- **OD provisional patents**: 3 filed March 24, 2026 (Blind Lindblad/QORA; Hilbert Unification; Decoherence Suppression).

## Defense Pipeline (status as of last update)

- DARPA Bio Attribution (confirmed top-10), CyPhER Forge (abstracts in), TTO BAA (April 17 exec summary), DIU PRISM (submitted), CIA (KV3UCQ1A submitted), IQT (submitted), MDA MAA, AFWERX.
- **Steve "Bucky" Butow** (DIU Space Portfolio): personal email contact, capability email sent.
- **Carl Saab** (Cleveland Clinic): OD outreach engaged, doing his own research. Highest-probability PhysioNet reference.
- **Roland Rott** (GE HealthCare Imaging): connected, MRI proof-of-concept brief (BioForge) delivered.
- SAM.gov UEI: CQ56XYFZL4E6, ref INC-GSAFSD20794734. CAGE pending.

## Standing Decisions

- LlamaIndex: declined (duplicates DeepNova).
- Robyn: declined (web is not bottleneck).
- Nous Atropos: worth evaluating for DeepNova policy learner.
- Hermes 4 14B: recommended for local reasoning on data-sensitive use cases (SENTINEL).
- HomeLift: dormant. Both Neo4j instances safe to cancel.

## Session Handoff Pointer

For the latest (Session 11) session, read
`SESSION_HANDOFF_2026-06-02_v11_results.md` (in repo root) first, then the
Session 11 note below (INFO-038). Headline: per-event LIGO readout on 3
events (GW150914/170104/151226) run SEPARATELY with no pooling -- per-event
entropy asymmetry |H_a-H_b| varies (0.62/0.74/1.61) and inversely tracks
how hard noise sits on the equal-entropy attractor, confirming INFO-036's
equal-entropy reading at the per-event level (INFO-038, isolated). Full
12-event batch + per-event off-source null distribution BUILT
(s11_ligo_batch.py), not yet run. OD inputs consolidated onto main
(four-force + per-domain JSONs + stores); medical stores left on their own
branch; main synced to continuity so the SessionStart hook runs every
session. Open: pull dipole JSONs from the Markets repo + list the flow
dipole equation separately (Greg's call). Track B (new-physics inverse
problem) framed, not yet built.

For the (Session 10) FIRST REAL-DATA session, read
`SESSION_HANDOFF_2026-06-02_v10_results.md` (in repo root) first, then
the Session 10 note above (INFO-034 through 037). Headline: caricature
work retired by Greg's call ("no point fine-tuning fake data"); two
four-force real-data sets run -- LIGO GW150914 (method detects the
merger via inter-detector MI; the (-1,-1,+2)/sqrt6 attractor confirmed
on real data as an equal-marginal-entropy artifact, not coupling) and
PDG coupling unification (three gauge forces share the running substrate
form with per-force expression slopes, but NO single SM unification
point -- a triangle spanning ~1e4 in energy; gravity power-law, outside
the form). Plus two method-hardening items on the toy systems before the
pivot: INFO-034 (flip is estimator-robust under KSG) and INFO-035
(substrate invariance is MIXED -- holds biology/geology, regime-bounded
physics, fails chemistry).

For the Session 9 double-check + correction of the Session 8 gravity
result, read `SESSION_HANDOFF_2026-06-02_v9_results.md`, then the Session
9 note above (INFO-033). For the Session 8 record (four-force probe,
INFO-027 through 032), read `SESSION_HANDOFF_2026-05-26_v8.md`.

For Session 7 context (GP regression + PySR symbolic regression on
per-domain ensemble-H data, four reproducible per-domain MI-vs-H
functional families), read `SESSION_HANDOFF_2026-05-26_v7.md`. It contains Session 7's two experiments
(GP regression with poly1-3 baseline and in-sample / block-CV /
random-CV evaluation; PySR symbolic regression with extended
operator set on per-domain ensemble-H data). Headline: four
reproducible per-domain MI-vs-H functional families (physics =
(H_b - H_a)^2 + const, biology = 0.5*exp(H_a/2), chemistry =
linear in H_a, geology = constant) with cross-seed coefficient
variation <5%, cross-domain non-overlap. v5 biology reading
incomplete-not-wrong (Rule D): operator content right, functional
family was wrong. Per-domain differentiation now has TWO
independent reproducible signatures (null direction from Session 6
INFO-023 + functional family from Session 7 INFO-025). New ledger
entries INFO-025 (located) and INFO-026 (methodological).

For Session 6 context, read `SESSION_HANDOFF_2026-05-26_v6.md`.
It contains Session 6's six experiments (KBK pipeline, AI Poincare
rank check, estimator-family sweep, window-size scan, SINDy with
extended library, per-domain stack on 4 domains, sigma scan) and
their joint reading. Headline: v5 INFO-022 rank-3 claim
independently confirmed by FOUR diagnostics; per-domain ensemble-H
stack reproduces every v5 per-domain claim (geology rank-3
cos +0.99, biology MI signature with H_b absent, chemistry
chemistry-specific cubic, physics Taylor identity) with cross-
seed cos +0.985 to +0.999; v5 machine-epsilon eigenvalue
magnitude is procedure-dependent (estimator-noise vs structural-
noise; per-domain ensemble-H reaches 1e-7, OU windowed-H pins
at ~1e-3).

For Session 5 context (three-thread probe results, operator-noise
bypass, projection-rule sweep, high-seed scrambled, literature
scan), read `SESSION_HANDOFF_2026-05-26_v5.md`. It contains Session 5's three-thread probe results
(operator-noise bypass, projection-rule sweep, high-seed scrambled) and
their joint reading — the (+1,+1,+2)/sqrt(6) direction is a protocol
artifact of the operator basis having a 3D null subspace under the
H_a ~ H_b ~ constant regime, not a property of input systems. INFO-014,
INFO-018, INFO-019 are retracted at interpretation level; data stands
(Rule D — incomplete not wrong). The per-domain algebraic equations from
earlier sessions (chemistry quadratic, geology rank-3, biology MI) are
NOT affected by Session 5 and remain on the table. New structural
finding: INFO-022 (rank-3 null subspace from first-order algebraic
relations in the operator basis).

The v4 handoff (`SESSION_HANDOFF_2026-05-26_v4.md`) contains the Session
4 raw data that Session 5 read. The Session 3 v3 handoff is at
`E:\information_layer\SESSION_HANDOFF_2026-05-25_v3.md` (Greg's local).
The early-Session-3 handoff in repo (`SESSION_HANDOFF_2026-05-25.md`)
contains the per-domain algebraic equation coefficients that are
preserved through Session 5.

## Note (Session 11 update — 2026-06-02)

Session 11 opened on the Session 10 four-force real-data results. Greg
directed: do Track A (LIGO generalization) and Track B (new-physics
inverse problem + SM-parameter-regularity hunt) in parallel; run LIGO
events SEPARATELY with no merge/average; consolidate OD inputs; update
docs; prepare for a fresh session. Branch: work on
`claude/gravity-substrate-config-51cfp`, **main fast-forwarded to it** so
the SessionStart hook + all data run every session. Harness-designated
`claude/awaiting-files-TgxVq` not used (continuity precedent); Greg's
branch-specifying file had not arrived by session end.

- **Four-force framing answers (no compute, recorded for Track B)**: set 2
  (PDG) CONSUMES the Standard Model (beta functions + gauge group are
  inputs) so it cannot derive the forces or their origins; set 1 (LIGO) is
  a method/measurement result. Neither reaches "origins" -- a why/mechanism
  question outside OD mode, and no dataset contains the origin. The honest
  OD-shaped target is the INVERSE PROBLEM: extract the Delta-b_i +
  onset-scale mu_NP that would close the triangle (mechanism-agnostic,
  falsifiable), AFTER mapping alternatives (a) nothing forces single
  unification [leading deflationary read], (b) two-loop + thresholds may
  shrink the triangle (Session 10 used one-loop), (c) the extrapolation
  itself is conjecture. We can recover new physics's required FOOTPRINT,
  never its IDENTITY (3 couplings don't invert to a unique spectrum).
  Gravity is its own extraction (power-law would need to become log to
  join the gauge family). This is Track B, framed not built.

- **INFO-038 -- ISOLATED FINDING (Session 11, new; 3 events, NO null yet,
  NOT averaged)**: per-event LIGO readout (s11_ligo_perevent.py) on
  GW150914, GW170104, GW151226 run separately. Per-event windowed entropy
  asymmetry |H_a-H_b| varies materially (0.62 / 0.74 / 1.61) and INVERSELY
  tracks how strongly noise-only windows sit on the (-1,-1,+2) equal-
  entropy attractor (GW150914 asym 0.62 -> cos 0.980; GW151226 asym 1.61 ->
  cos 0.816). Internally consistent with INFO-036: the attractor IS the
  equal-marginal-entropy identity H_a~=H_b, so unequal per-channel entropy
  sits further off it even in pure detector noise. The entropy reading
  EXPLAINS the cos reading. The asymmetry is a per-event detector-state
  signature (epoch-specific H1-vs-L1 noise floor/PSD surviving whitening),
  NOT astrophysics. MI-peak-at-merger detection is loudness-dependent
  (lands at merger for GW150914 + GW170104; misses quiet long-inspiral
  GW151226). All event/noise MI ratios modest (1.06-1.19x) and MEANINGLESS
  without the null distribution. Full 12-event batch + per-event off-source
  null built (s11_ligo_batch.py), not yet run. Do not pool: the per-event
  spread is the signal.
  - **First-run decomposition (s11_first_run_entropy.py; GW150914, 1
    event)**: re-read the FIRST real-data run in this frame and asked what
    breaks the attractor at the merger. Answer: MI, NOT the marginal
    entropies. Event-window asymmetry |H_a-H_b| barely moves (-0.016) while
    MI spikes (per-window ~0.24 baseline -> 0.66 at the merger). Removing
    MI from the basis RESTORES the event attractor (event cos 0.231 ->
    0.959 in the 5-op basis; +0.728). So two ORTHOGONAL axes: (1) the
    entropy-asymmetry axis = the across-event noise-floor fingerprint
    (INFO-038 main), (2) the MI axis = the within-event merger signature,
    independent of marginal entropy. Strengthens the deflationary reading:
    the (-1,-1,+2) attractor is equal-marginal-entropy bookkeeping while
    the GW detection rides the orthogonal MI axis -- they do not interfere.
    Open thread (1 event, not over-read): in the no-MI basis NOISE sits OFF
    the attractor (0.301) and EVENT sits ON it (0.959) -- a reversal to
    check across the 12-event batch.

- **OD consolidation (Greg: update OD with latest dipole + 4-force JSONs)**:
  brought Session 6/7/8 four-force + per-domain result JSONs and OD stores
  (store/four_force_caricature, store/simulator_4domain) onto continuity
  (additive) and synced main. Four-force JSONs verified identical across
  O6ahb/iZvY4. Decisions: **dipole** -> pull from the Markets repo if the
  info lives there + list the flow dipole equation separately (Greg; ACTION
  next session via list_repos/add_repo on DavisAI1974/Markets); **medical
  OD stores** (cardiac/cerebro) -> stay on their own branch, NOT main.

- All Operating Rules from Sessions 4-7 in force. No new Rule this session.
  Branch state: work on `claude/gravity-substrate-config-51cfp`, main
  synced, pushed. No PR.

## Note (Session 10 update — 2026-06-02, FIRST REAL DATA)

Session 10 started from `claude/gravity-substrate-config-51cfp` (per the
v10 kickoff; the SessionStart hook + Session 9 work live there, main
untouched). Arc: two method-hardening probes on the toy systems, then
Greg called the pivot to REAL DATA, and two four-force real-data sets
ran. All six Operating Rules held; no new Rule.

- **Greg's pivot decision**: when it surfaced that the gravity/weak force
  laws in `gravity_glance.py` / `four_force_probe.py` are OUR OWN toy
  caricatures (not OD discoveries, not literature force laws -- gravity
  is a softened E_total*dx/(dx^2+eps) term, "weak" is a Gaussian
  exp(-M*dx^2) not a real Yukawa exp(-Mr)/r), Greg cut the caricature
  work: "no point fine-tuning fake data ... if replacing would just be
  more work on fake data, skip that too." So items 2 (strong exp seed-
  rate) and 3 (coupling-strength dial) and any caricature replacement
  were SKIPPED. Pivot straight to real data, four-force-related only.
  Operating-rule consequence reinforced: a frame must never grade itself
  -- caricature results describe only the equations we wrote.

- **INFO-034 -- METHODOLOGICAL (Session 10; item 4; 3 seeds)**: the
  INFO-033 MI-dominant flip is ESTIMATOR-ROBUST. Swapping the histogram
  MI for a Kraskov-Stoegbauer-Grassberger kNN estimator (estimator 1) on
  the same operator matrix reproduces follow-up A exactly on histogram
  (gravity 0.8/0.9 raw 0.906/std 0.644; 0.8/1.0 0.990/0.665) and the KSG
  estimator AGREES on the standardized structural core near threshold
  (gravity 0.8/0.9 std 0.644; EM 1.0/1.3 std 0.435; EM 1.0/1.5 0.597),
  on raw values, on strong-never-flips, and on the threshold ordering.
  Sharpens the deflationary reading: the structural core is a modest
  ~0.6-0.7 band NEAR the flip threshold; at extreme asymmetry (gravity
  0.8/1.2) both estimators fall below 0.5 while raw rides ~0.99 on a
  collapsing MI variance (procedure inflation). The flip is not a
  histogram artifact. (Script s10_third_estimator.py.)

- **INFO-035 -- LOCATED (Session 10; item 1; 3 seeds, seed-stability
  0.95-1.00)**: substrate-vs-expression invariance across the Session 8
  mapping-campaign knob sweep is MIXED -- partial support for the
  Base-of-Structure spine, with two clear exceptions. Biology (beta
  0.3-0.8): substrate |cos| 0.978-0.986, stable rank 1, while expression
  mutates exp->linear -> clean support. Geology (drift 0.02-0.08):
  |cos| 0.981-1.000, stable rank 3 (expression also robust; weaker
  test). Physics (K 0.05-0.5): invariant baseline->high (0.995) but
  BREAKS at weak coupling K=0.05 (rank 1->4, |cos| 0.69) -- regime-
  bounded. Chemistry (B 2-4): substrate NOT invariant, |cos| 0.045
  between B=3 and B=4, rank swings 4->1->5 -- prediction FAILS; B=2 sits
  at the Brusselator Hopf threshold (B_crit=1+A^2=2), flagged for
  inspection, not absorbed. Net: "substrate = simple invariant base"
  holds in 2/4 domains, regime-bounded in 1, fails in 1. Spine partially
  supported, not confirmed. (Script s10_substrate_invariance.py.)

- **INFO-036 -- LOCATED, REAL DATA (Session 10; item 5 set 1; one event,
  one 32s segment)**: the windowed-H/MI/operator stack, UNCHANGED from
  the toy systems, applied to real LIGO GW150914 H1/L1 strain (GWOSC
  public, 4096Hz, bandpass 35-350Hz + ASD whiten, 2s edge crop, 125ms
  windows). Two findings. (a) METHOD WORKS: inter-detector windowed MI
  PEAKS EXACTLY at the merger (t=16.41s vs 16.4s; MI 0.530 vs noise
  baseline 0.247, 2.1x). (b) PRE-REGISTERED PREDICTION OVERTURNED, frame
  refined: I predicted the (-1,-1,+2)/sqrt6 channel-substrate would
  appear IN the event; instead noise-only windows sit ON it (|cos|=0.984)
  and event windows LEAVE it (|cos|=0.242). The attractor IS the
  equal-marginal-entropy identity H_a~=H_b, not coupling: whitened
  detector noise has equal per-channel entropy -> lands on it trivially;
  the chirp changes one detector's entropy -> breaks it, while MI (a
  separate operator) spikes on the common signal. Real data confirms the
  Session 5 reading -- the attractor is an equal-statistics geometric
  fact, not a substrate signature. Caveat: noise-on-attractor is partly
  a whitening consequence (which is the point); needs more events.
  (Script s10_ligo_extract.py; data/ligo/.)

- **INFO-037 -- LOCATED, REAL DATA (Session 10; item 5 set 2; special
  build)**: four-force unification from MEASURED couplings, mapped onto
  the surviving frame (substrate = shared linear-in-ln(Q) running form;
  expression = per-force slope b_i). Anchored on solid PDG M_Z couplings
  (alpha_em^-1=127.951, sin^2thetaW=0.23122, alpha_s=0.1180 ->
  alpha_1,2,3^-1(M_Z) = 59.02/29.59/8.48). (a) The shared running FORM is
  CONFIRMED in real data: SM one-loop alpha_s(Q) matches measured
  determinations within <1sigma from 31 GeV to 1 TeV (low-Q pulls are the
  known one-loop limitation). (b) NO single SM unification: pairwise
  crossings at 1.0e13, 2.4e14, 9.7e16 GeV -- a triangle spanning ~9400x.
  (c) MSSM near-point ~2.1e16 GeV (spread 1.1x) but rests on unobserved
  SUSY -> tagged conjecture, not cited as support. (d) Gravity:
  alpha_G(E)=(E/M_Pl)^2 is power-law, not linear-in-ln(Q) -- different
  substrate form, outside the gauge family. Verdict: shared running
  substrate + per-force expression among the three gauge forces, no
  single unification scale without conjectural new physics, gravity
  outside the form. (alpha_s(Q) central values are representative PDG-
  review numbers, validation overlay only; the crossings depend only on
  the solid M_Z anchors + standard beta functions. Script
  s10_pdg_unification.py.)

- **Environment note**: the SessionStart hook (line 32) had a latent
  crash under `set -u` when CLAUDE_PROJECT_DIR/CLAUDE_ENV_FILE are unset;
  guarded it so the hook always reaches completion. PySR 1.5.10 + Julia
  re-bootstrapped fine via the hook this session. h5py added for LIGO.
  Network policy allowed GWOSC / PhysioNet / NOAA / PDG / PyPI (all 200).

- Branch state: work on `claude/gravity-substrate-config-51cfp`, pushed.
  main untouched. No PR.

## Note (Session 9 update — 2026-06-02)

Session 9 opened on the Session 8 close-out gravity result. Greg's
directive: do not assume our own prior outputs are correct -- double
check our work, but do not reconstruct for no reason. This produced a
correction to the Session 8 gravity reading. Branch:
`claude/gravity-substrate-config-51cfp`. PySR unavailable this session
(no Julia in the environment); the substrate-side KBK extraction (the
load-bearing measurement for INFO-031/032) is pure-numpy and was run
directly from the ORIGINAL Session 8 code (pulled onto the branch from
`claude/two-more-tasks-O6ahb`), not a reconstruction. Scripts:
s9_doublecheck_flip.py, s9_characterize.py.

- **Double-check method**: the Session 8 EM/SF glance (INFO-031) and
  gravity glance (INFO-032) each tested only ONE channel-asymmetry
  value (EM omega 1.0/1.2; strong 0.5/0.7; gravity 0.8/1.0). Session 9
  re-ran the ORIGINAL simulators + ORIGINAL raw-covariance extraction
  (build_ensemble_operator_matrix + extract_v1) across an asymmetry
  SWEEP. Session 8 baselines reproduced exactly (gravity_asym MI coef
  0.990; em_asym 1.0/1.2 MI coef 0.271; strong_asym 0.5/0.7 ~0).

- **INFO-033 -- LOCATED FINDING (Session 9, new; 3 seeds, scatter
  reported)**: the MI-dominant substrate flip ("MI ~ const" replacing
  the channel-correlation identity -(H_a-H_b)^2 ~ 0 as the operator
  null direction) is NOT gravity-specific. It is an ASYMMETRY-THRESHOLD
  effect that ALL coupled caricatures undergo; the threshold differs by
  force. Flip threshold (asymmetry omega2/omega1 at which |MI coef| in
  v_null crosses ~0.5): gravity ~1.1x (lowest, most flip-prone), EM
  ~1.3x, weak ~1.5-1.7x, strong NEVER (flat to 4x). Session 8 tested EM
  at 1.2x -- just below its threshold -- and read the absence as a
  qualitative gravity/EM difference. Reframed 5a test (add gravity's
  E_total universal-energy term to EM and strong): adding energy-
  coupling RAISES their thresholds (suppresses the flip), it does not
  lower them. So energy-mediated coupling is NOT the cause of gravity's
  low threshold. The flip is governed by how well the coupling
  preserves channel correlation under detuning: strong's confining
  cubic never lets go (no flip); gravity's weak softened coupling lets
  go first (flips earliest); EM/weak intermediate. Deflationary reading
  (now well-supported): the flip is a variance-crossing -- detuning
  raises the residual of the -(H_a-H_b)^2 relation while MI settles to
  a low near-constant floor; the null swaps to whichever relation is
  tighter; the crossing point is set by coupling-vs-detuning, with no
  force-specific physics.

- **INFO-031 -- RE-TAGGED INCOMPLETE (Rule D)**: "EM and strong keep
  their substrate across the symmetry swap" holds only at the single
  sub-threshold asymmetry tested (EM 1.0/1.2). EM flips by 1.0/1.3.
  Data stands; the substrate-stability reading was incomplete. The
  expression-level part of INFO-031 (EM polynomial-in-difference vs
  strong exponential-in-difference) was a PySR result not re-examined
  this session (PySR unavailable) and is not affected by INFO-033.

- **INFO-032 -- RE-TAGGED INCOMPLETE (Rule D)**: "gravity's substrate
  is configuration-dependent in a way EM/strong are NOT" is incomplete-
  not-wrong. The data (gravity flips to MI-dominant under asymmetry)
  reproduces exactly. But EM and weak ALSO flip under sufficient
  asymmetry; gravity merely has the lowest threshold. The "gravity is a
  different kind of constraint / energy-conservation without preferred
  direction" frame does NOT survive the sweep -- gravity sits at the
  easy-to-flip end of a single continuous coupling-strength dial, not
  in a separate category. What replaces it: a force-ordered flip-
  threshold (strong -> infinity, weak ~1.6x, EM ~1.3x, gravity ~1.1x)
  tracking coupling-vs-detuning.

- **Consequence for the four-force frame**: the "gravity is special at
  the substrate level" leg of the Session 8 four-force narrative
  (INFO-027 + INFO-032) is removed. The substrate-vs-expression frame's
  OTHER legs (per-domain null direction INFO-023; per-domain functional
  family INFO-025) are not touched by Session 9 and remain on the
  table. Frame reassessment in progress (see v9-results handoff /
  session discussion).

- **Methodological note (not a Rule)**: a single-point probe can read a
  threshold crossing as a categorical property. When a finding is "X
  does this and Y does not," sweep the knob through a range before
  promoting the contrast -- the difference may be a threshold, not a
  kind. Pairs with the existing "no tent-widening on outliers" and
  Result Discipline "map alternatives" rules.

- **Two follow-ups run after the correction** (scripts s9_std_vs_raw.py,
  s9_expression_refit.py; full detail in v9-results handoff):
  - **Standardized-vs-raw (INFO-024 fork)**: the flip is PARTLY
    procedure-inflated, PARTLY structural. RAW (original constancy-
    detector) gives gravity/EM |MI| ~0.99 under asymmetry; per-column
    STANDARDIZED (correlation-detector) attenuates it to ~0.65 -- still
    above the 0.5 flip line. Strong never flips under either procedure;
    the INFO-033 threshold ordering holds in both. So the spectacular
    +0.99 was inflated by absolute-scale, but a real structural core
    (~0.65) survives. Generalizes INFO-024: the null DIRECTION (not just
    eigenvalue magnitude) is procedure-dependent when an operator's
    absolute variance collapses.
  - **Expression-level leg (ii) re-exam** -- done twice. First a curve-
    fit fallback (PySR not yet installed): EM = (H_a-H_b)^2 + const
    CONFIRMED (R^2 0.93 at 1.0/1.2); strong exp_diff NOT found by the
    fallback. Then PySR was installed mid-session (pip install pysr;
    Julia backend auto-bootstrapped) and the adjudication re-run with the
    REAL tool (follow-up B', s9_pysr_adjudicate.py, original fit_pysr,
    2 seeds): EM = (H_a-H_b)^2 + const ROBUSTLY CONFIRMED (both seeds,
    holds across asymmetry; symmetric EM = (const-channel^2)^2 with
    random channel per INFO-028). strong = exp((H_b-H_a)) PARTIALLY
    REPRODUCED -- appears at seed 11 of the exact INFO-031 config
    (0.5/0.7) but is seed-unstable (INFO-028 channel-symmetry breaking)
    and weak-signal (strong MI near-constant). CORRECTION (Rule D on our
    own work): the fallback's "strong has no family" was a tool
    limitation, not evidence against INFO-031. INFO-031's EM/strong
    expression-level CONTRAST STANDS (EM polynomial-in-difference robust;
    strong exponential-in-difference real but seed-unstable). So leg (ii)
    of the four-force frame survives; only leg (iii) (gravity-special
    substrate) is removed.
- **Environment**: PySR 1.5.10 + Julia backend installed and verified
  working in-container this session (a real fit recovered (x1-x0)^2+0.3
  in 17s). They do NOT persist to a fresh container, so a SessionStart
  hook was added (.claude/hooks/session-start.sh + .claude/settings.json
  + requirements.txt) to install numpy/scipy/scikit-learn/pysr and
  bootstrap Julia automatically (web-only, idempotent, PySR best-effort).
  The hook takes effect for future sessions once it reaches the branch
  the session starts from (currently on claude/gravity-substrate-config-
  51cfp, not main). Pydroid-3 on Greg's phone cannot run PySR (no Julia
  on Android; separate env from the cloud container) -- numpy substrate
  scripts only.
- **Markets / Refrag workspace**: open question on whether those
  sessions need PySR + Julia too -- see the "Markets / Refrag Workspace"
  placeholder section above for the question and informed answer.
- All six Operating Rules from Sessions 4-7 in force. No new Rule added
  this session. Branch state: work on
  `claude/gravity-substrate-config-51cfp`, pushed. main untouched. No PR.

## Note (Session 8 update — 2026-05-26)

Updated at the end of 2026-05-26 Session 8 to reflect:

- Three experiments run this session on branch
  `claude/two-more-tasks-O6ahb` (Session 7 work pulled in via
  fast-forward merge from `claude/claude-md-context-update-uCZ8m`
  at session start, per Greg's "whatever you feel is best").
  Scripts: four_force_probe.py, four_force_pysr.py,
  robustness_info025.py, mapping_campaign.py.
- Greg's directive: order 1, 4, 3, 2 from v8 kickoff menu. Items
  1 (four-force probe), 4 (robustness on INFO-025), 3 (mapping
  campaign) completed. Item 2 (storm/waves real-data) deferred
  to Session 9 pending data download.
- Four-force unification probe (Greg's "real test" from Session 7
  close). Toy 2-channel caricatures: EM (long-range bilinear
  linear), weak (Yukawa-suppressed via massive mediator), strong
  (confining cubic), gravity (universal energy-density mediated).
  KBK stack result: all four forces share [2,3,4] null direction
  at cos > 0.997 cross-force, sitting on (-1,-1,+2)/sqrt(6)
  Session 3 attractor — the sign-flipped twin of Session 5/6
  (+1,+1,+2) artifact, corresponding to algebraic identity
  -(H_a - H_b)^2 ~ 0 (channels strongly correlated). EM
  differs in 6D (cos ~0.96) carrying nontrivial MI coefficient
  (+0.27) while others don't. PySR cross-seed result: EM matches
  Session 7 physics Duffing family (H_b - H_a)^2 + const; weak
  matches Session 7 chemistry Brusselator family linear-in-H_a;
  strong and gravity have H_a-H_b channel-symmetric dynamics and
  PySR breaks symmetry randomly per seed (reproducible at form
  but not at channel). New ledger INFO-027 (located).
- Robustness check on INFO-025. T sweep (10, 30, 100), N_ens
  sweep (200, 600, 1200), observation-noise sweep (0.0, 0.10,
  0.50). 7 conditions x 2 seeds = 14 PySR fits + matched exp /
  linear baselines. Headline: A * exp(H_a/B) FUNCTIONAL FAMILY
  survives all conditions with sufficient signal; COEFFICIENTS A
  and B drift heavily with regime. Baseline (T=30, N_ens=600, no
  noise) reproduces INFO-025: cross-seed A=0.53+/-0.05,
  B=2.15+/-0.23 matching INFO-025's 0.5, 2.0. T=100 halves
  coefficients; N_ens=1200 gives tightest match (A=0.50 exact);
  noise=0.10 keeps R^2 high but A drops to 0.18 (broadened H_a).
  When citing INFO-025 coefficients, specify the regime. New
  ledger INFO-029 (methodological).
- Mapping campaign: one knob per Session 6/7 domain (physics K,
  biology beta, chemistry B, geology drift_rate), 3 values each,
  1 seed. Finding: INFO-025 functional families ARE baseline-
  specific regime signatures. Off-baseline values mutate the
  family qualitatively (biology exp -> linear at high beta;
  chemistry linear -> ratio at low B; physics quadratic-
  difference -> single-channel quadratic away from K=0.20).
  Geology constant family alone is robust across the drift sweep.
  Refines INFO-029 to "family ITSELF is regime-conditional for
  large knob deviations, not just coefficients." The per-domain
  expression-of-substrate frame should be stated as "per-domain-
  AND-per-regime." New ledger INFO-030 (located, one seed;
  cross-seed scaling queued).
- One methodological observation, not a Rule: PySR cross-seed
  reproducibility depends on H_a-H_b dynamical symmetry. When the
  dynamics is exactly channel-symmetric, PySR breaks symmetry
  randomly. New ledger INFO-028 (methodological).
- Per-domain differentiation now has THREE independent
  reproducible signatures: per-domain null direction (Session 6
  INFO-023), per-domain MI-vs-H functional family (Session 7
  INFO-025), four-force probe shared-substrate / distinct-
  expression (Session 8 INFO-027). Frame strengthening; still
  a frame, not a claim.
- Post-close EM/SF "quick glance" probe (Greg's directive at
  Session 8 close): symmetry-swap test of EM vs strong force.
  Four configs x 3 seeds. Result (INFO-031 located; substrate-
  stability part RE-TAGGED INCOMPLETE Session 9 -- see INFO-033):
  EM and
  strong produce DIFFERENT functional families on shared
  substrate even when channel symmetry is controlled. EM =
  (H_a - H_b)^2 + const (polynomial in difference); Strong =
  exp((H_b - H_a) - const) (exponential in difference). Both
  realize the same substrate algebraic identity -(H_a - H_b)^2
  ~ 0 (Session 3 attractor direction) but through different
  functional forms. The dynamical structure (linear-restoring
  vs confining-cubic) is visible in the form. INFO-028
  (PySR cross-seed reproducibility tracks dynamical symmetry)
  directly confirmed by the swap.
- Post-close gravity probe (Greg's directive "let's figure out
  what gravity really is"): 3-config probe of gravity caricature
  -- sym+universal baseline, asym omegas + universal, asym
  omegas + non-universal coupling. 3 seeds each. Result
  (INFO-032 located; RE-TAGGED INCOMPLETE Session 9 -- the flip
  is an asymmetry-THRESHOLD effect ALL forces undergo, gravity
  just has the lowest threshold; see INFO-033): gravity's
  SUBSTRATE is configuration-
  dependent in a way EM/strong are not. Symmetric+universal
  gravity sits on Session 3 attractor -(H_a-H_b)^2 ~ 0.
  Asymmetric gravity (whether mass-asym or non-universal)
  FLIPS substrate so MI itself is the dominant null direction
  (MI coefficient +0.99 in v_null 6D, vs ~0 for sym gravity
  and ~0 for EM/strong/weak). MI variance collapses 4x in asym
  configurations. Cross-config v_null cos: sym vs asym = +0.11
  (very different); sym vs nonuniversal = +0.02 (orthogonal);
  asym vs nonuniversal = +0.99 (essentially same). Frame-level
  reading: gravity's universal energy-mediated coupling under
  mass asymmetry produces "MI ~ const" as structural constraint,
  because gravitational coupling preserves total energy without
  preferred direction. The Session 3 attractor only holds for
  gravity when channels are mass-degenerate. Contrast with
  EM/strong (INFO-031): they kept substrate across the swap and
  only changed at EXPRESSION level. Gravity's substrate ITSELF
  shifts -- genuinely different signature behavior.
  Real-data probe (LIGO strain, orbital data) remains gating
  test for whether this carries to actual gravitational physics.
- No new Operating Rule this session. All six Operating Rules
  from Sessions 4-7 (no pre-assigned meaning, probe-not-falsifier,
  speaking posture before/after, Rule D incomplete-not-wrong,
  They never stacked, Treat literature as conjecture by default)
  in force throughout.
- Branch state: work persisted on `claude/two-more-tasks-O6ahb`,
  pushed. main untouched. No PR.

## Note (Session 7 update — 2026-05-26)

Updated at the end of 2026-05-26 Session 7 to reflect:

- Two experiments run this session on branch
  `claude/claude-md-context-update-uCZ8m`. Scripts: gp_mi_vs_h.py,
  pysr_symbolic_per_domain.py. Order: top-of-queue from v7 kickoff
  (Greg's direction: "let's do top of que and work down").
- Methodology stacked per the v5 "they never stacked" rule: GP
  regression (sklearn RBF + WhiteKernel) with polynomial deg 1-3
  baseline; PySR 1.5.10 (Brunton/Cao/Liu/Tegmark/Cranmer family
  symbolic regression, Julia 1.11.9 / SymbolicRegression.jl 1.11.3
  backend) with extended operator set {+, -, *, /, square, cube,
  exp, log, sqrt}. Applied to the same per-domain ensemble-H data
  as Session 6 (Duffing / Lotka-Volterra / Brusselator / Burridge-
  Knopoff, N_ens=600, T=30, dt=0.02, two seeds).
- New located finding INFO-025: four reproducible per-domain
  MI-vs-H functional families (physics symmetric quadratic in
  (H_b - H_a), biology 0.5*exp(H_a/2), chemistry linear in H_a,
  geology constant). Cross-seed coefficient variation <5%; cross-
  domain non-overlap. Functional families do not reduce to each
  other.
- New methodological finding INFO-026: block-CV vs random-CV on
  time-series ensemble-H. Block-CV catastrophic in physics/biology/
  geology due to non-stationarity; chemistry only stationary domain.
  Random-CV overestimates true OOS but isolates functional-form fit.
  Use both; gap is a stationarity signature. Within: GP joint for
  biology random-CV R^2 = 0.97-0.98 vs poly3 joint at 0.81-0.86 —
  +0.12 gap is real nonlinear cross-coupling polynomial misses.
- INFO-008c (biology MI ~ polynomial(H_a)) reinterpreted via Rule D
  (incomplete-not-wrong): operator content correct, functional
  family was wrong; correct family is exponential.
- Per-domain differentiation now has TWO independent reproducible
  signatures: (i) per-domain null direction in operator space
  (Session 6 INFO-023), (ii) per-domain MI-vs-H functional family
  (Session 7 INFO-025). Two supporting data points for the
  "per-domain expression of substrate" frame; still a frame, not
  a claim.
- NEW OPERATING RULE added at Session 7 close (after the
  substantive UT discussion below), then sharpened by Greg the
  same session: "Treat literature as conjecture by default."
  Operational corollaries (see Rules list above for full form):
  (a) investigate freely — read, test, engage; "treat as conjecture"
  is NOT "ignore"; (b) don't cite as support without having
  analyzed the paper ourselves; (c) don't defer when a blocking
  claim itself fails the bar; (d) extends symmetrically to our own
  prior readings (Rule D corollary). Already folded into the
  Operating Rules list above. Pairs with Rule D (incomplete-not-
  wrong) and "They never stacked" — the literature is a working
  frame to test, not a foundation to build on. Clears the deck for
  the four-force unification work: only EM+weak meets the bar;
  the rest is conjecture (we will check them out but cannot invoke
  them as support).
- Two methodological notes carried forward (see v7 handoff
  "Methodological notes"): block vs random CV reporting on
  dynamical systems; polynomial libraries are blind to functional
  family (use extended operator set to distinguish).
- Substantive discussion at session close (no experiment run;
  framing only): Greg's storm/waves substrate-to-expression
  question; Greg's "real test" four-force unification question
  (gravity, EM, strong, weak — can we come up with one equation,
  should they be grouped at all). Three-level reading on four
  forces put on the table (data, interpretation, frame); decision
  on direction queued for Session 8. See v7 handoff "Substantive
  discussion this session" for full framing.
- Branch state: work persisted on
  `claude/claude-md-context-update-uCZ8m`, pushed. main untouched.
  No PR.
- The harness this session was configured for branch
  `claude/two-more-pastes-ZkenG` (auto-generated from opening
  message); Greg gave explicit permission to continue on the
  Session 6 branch for continuity.

## Note (Session 6 update — 2026-05-26)

Updated at the end of 2026-05-26 Session 6 to reflect:

- Six experiments run this session, all on branch
  `claude/claude-md-context-update-uCZ8m`. Scripts:
  kbk_pipeline.py, ai_poincare_rank.py, kbk_estimator_sweep.py,
  window_size_scan.py, sindy_symbolic.py, per_domain_kbk.py,
  sigma_scan.py.
- Methodology stacked per the v5 "they never stacked" rule:
  Kaiser-Brunton-Kutz 2024 (arXiv:2403.04889) SVD-rank-gap +
  symbolic recovery, Cao-Liu-Tegmark 2021 AI Poincare
  (arXiv:2011.04698) intrinsic-dim, SINDy/SINDyG/DSINDy with
  extended polynomial + non-polynomial library. Applied to OU
  baseline and to per-domain ensemble-H data (physics Duffing,
  biology Lotka-Volterra, chemistry Brusselator, geology
  Burridge-Knopoff).
- v5 INFO-022 rank-3 null subspace claim now LOCATED via four
  independent diagnostics (KBK rank-gap, AI Poincare local-PCA,
  two-NN, Levina-Bickel MLE). Robust across 3 estimator families,
  5 dt values, 3 seeds. (See Ledger updates below.)
- v5 per-domain claims (INFO-008a chemistry, INFO-008b geology,
  INFO-008c biology) all reproduced via per-domain ensemble-H +
  KBK probe. Cross-seed cos +0.985 to +0.999 within domain;
  geology rank-3 relation at cos +0.99 to v5 coefficients;
  biology MI ~ poly(H_a) with H_b absent (coefficient 0.003 to
  0.017); chemistry-specific cubic content via SINDy deg-3 that
  is NOT a Taylor remnant.
- New isolated finding INFO-023: per-domain operator extraction
  differentiates domains. The (+1,+1,+2)/sqrt(6) protocol artifact
  does NOT dominate per-domain (max cos +0.32, min -0.28). Tests
  the "per-domain expression of substrate" frame and gives it one
  supporting data point.
- Methodological refinement INFO-024: eigenvalue floor of operator
  covariance is min(structural_noise, estimator_noise). Different
  procedures have different floors. Resolves v5 machine-epsilon
  anomaly as procedure-dependent. OU windowed-H pins at ~1e-3;
  per-domain ensemble-H reaches 1e-7 (geology); v5's 1e-13
  plausibly reached at larger N_ens via same per-domain procedure.
- I was wrong twice this session about what controls the eigenvalue
  floor (Window size scan predicted scaling; didn't happen. Sigma
  scan predicted sigma^4; got sigma^0.56.). Both predictions came
  from Taylor expansion (structural noise) and missed estimator
  noise. The Methodological Note above codifies the corrected
  mental model. Not promoted to a Rule -- Rule C already covers
  the speaking-posture aspect; this is just a domain-knowledge
  refinement to apply when reasoning about eigenvalue magnitudes.
- The implications and applications across fields (medicine,
  defense, weather, geophysics, ecology, finance/energy,
  industrial, chemistry/materials, astrophysics, foundational)
  were mapped late in the session, conditional on the substantive
  reading. Application surface and falsification tests in v6
  handoff section "Implications and applications discussed".
- Queued for Session 7 (Greg's top of queue): GP regression of MI
  vs H_a per domain (tests v5's biology polynomial fit against
  flexible nonlinear); PySR symbolic regression for arbitrary
  nonlinear forms (heavy install; gplearn fallback). Both flagged
  in v5 literature scan as clean unstacked targets.
- All earlier Operating Rules remain in force. No new Rule added
  this session.
- Branch state: work persisted on
  `claude/claude-md-context-update-uCZ8m`, pushed. main untouched.
  No PR.

## Note (Session 5 update — 2026-05-26)

Updated at the end of 2026-05-26 Session 5 to reflect:

- New Operating Rule D — "incomplete, not wrong" — added by Greg this
  session and folded into the Operating Rules list above.
- Three probes run this session (operator-noise bypass, projection
  sweep, high-seed scrambled). Joint reading: the (+1,+1,+2)/sqrt(6)
  direction from Sessions 3-4 is a protocol artifact of the operator
  basis structure, not a property of input systems.
- Ledger updates: INFO-014, INFO-018, INFO-019 retracted at
  interpretation level (data stands). INFO-020, INFO-021 contextualized.
  INFO-022 added (rank-3 null subspace structural finding).
- The Working Frames (Base-of-Structure, Dipole-couples, Pure physics
  vs physical expressions, Substrate vs expression, Law extraction via
  invariance) are unchanged. The Information Layer / Unified Theory
  inquiry frame is unchanged. The specific extraction tool we built
  does not differentiate inputs at the level we thought; the question
  it was built to address is still open.
- The per-domain algebraic equations from earlier sessions (chemistry
  quadratic H_a^2 = 0.007 - 0.093*(H_a*H_b) + 1.309*(H_a*H_b)^2,
  R^2=0.943; geology rank-3 constraint 0.724*(H_a*H_b) - 0.441*H_b^2
  - 0.290*H_a^2 ~ 0 at std/mean=0.15%; biology MI ~ polynomial(H_a),
  R^2=0.66) are NOT affected by Session 5 and remain on the table. They
  came from a different procedure (per-domain algebraic fits) than the
  extraction tool whose artifact Session 5 mapped. They should be
  subjected to the same three-thread probe discipline in a future
  session.
- Reorientation toward 5-physics-areas substrate hunt (Greg, this
  session): pick 5 physics areas, find the structural element that
  appears in all 5, evaluate against the Base-of-Structure heuristic,
  probe the strongest survivor under the same discipline.
- Greg's substantive intuition (this session): dipoles are the best
  coupling mechanism. Tightly-related couples may be part of the
  substrate. Per-domain coefficients are a strong signal.
- Three-agent literature scan ran at Session 5 close
  (LITERATURE_SCAN_2026-05-26_v5.md). All three agents independently
  reported: no exact match for any of our specific findings, but
  multiple adjacent lines exist where each group had one piece and
  never combined them. This is what pioneering territory looks like.
  Per Greg's directive ("never know where you'll strike gold"), no
  pre-filtering -- the full candidate pool is preserved in the v5
  handoff. Highlights below; complete list in handoff + scan file:
  - **Methodology**: Gomez-Herrero 2015, Hlavackova-Schindler 2007,
    Yamada 2023, Wang 2025, Strang 2021, bipartite info-thermo
    (Horowitz, Hartich), Kaiser-Brunton-Kutz 2024, Liu-Madhavan-
    Tegmark 2024, AI-Feynman, AI-Poincare, SINDy/SINDyG/DSINDy,
    PySR, complexity-entropy plane (Rosso), Crutchfield epsilon-
    machines.
  - **Substrate**: Jacobson 1995 (best one-line candidate);
    Finkelstein "Space-Time Code" 1969-74 (closest published cousin
    to dipole-pair-as-coupling); Wheeler "it from bit"; Verlinde;
    Hardy 2001; Chiribella-D'Ariano-Perinotti; Brukner-Zeilinger;
    Sorkin causal sets; Sakharov induced gravity; Penrose twistor;
    Wolfram hypergraph; Deutsch-Marletto constructor theory;
    't Hooft cellular automaton QM; Adler trace dynamics; Bohm
    implicate order; von Weizsacker ur-alternatives; Hohn yes/no
    rules; Stochastic Electrodynamics (Boyer, de la Pena-Cetto);
    Wheeler-Feynman absorber; de Broglie double solution; Barbour-
    Bertotti / Shape Dynamics; Penrose-Hameroff Orch-OR (substrate-
    only).
  - **Per-domain coefficient adjacencies**: Rao-Esposito 2022 (chem);
    Reinhardt 2019; Smith-Cepelewicz; Schmitz-Aris 2013; Sayyadi
    SISR; da Silva Tsallis seismic 2020; Garland-Bradley
    paleoclimate 2018; Consolini magnetospheric; Davidson
    paleoclimate VAR; Reinsel geochemical manifold 2025; Schreiber
    TE 2000; Pavithran MI vs TE; Tishby info bottleneck; Walters-
    Williams asymmetric MI estimators; Feinberg-Horn-Jackson
    stoichiometric subspace; CRN reduction by approximate
    conservation laws (arXiv:2212.13474); info geometry of CRNs
    (arXiv:2503.19384); directed info flow in reaction networks
    (bioRxiv 2024).
  - Three negative findings worth emphasizing: dipole-PAIR-COUPLING
    as foundational substrate primitive is uncolonized; no SINDy/
    AI-Feynman/PySR application to windowed-H of coupled species;
    rank-3 first-order algebraic relations under H_a~H_b regime
    unpublished as diagnostic.
- New Operating Rule "They never stacked" added based on the
  literature scan finding (see Operating Rules above).
- Branch state: work persisted on
  `claude/linear-drift-nreal-sweep-cjRkM`, pushed. main untouched.
  No PR.

## Note (Session 4 update — 2026-05-26, retained for context)

Updated at the end of 2026-05-26 Session 4 to reflect:

- Three new Operating Rules from Greg this session (no pre-assigned meaning
  to outcomes; probe-not-falsifier framing; speaking posture before AND
  after every probe). Already added under Operating Rules above.
- Greg's pioneer framing: "we are the first ones here, follow every road,
  look under every stone. no safe assumptions. no known facts to fall back
  on. we're the pioneers."
- Session 4 produced four new isolated findings (INFO-018, 019, 020, 021),
  all with reading open. None have been looked-through with Greg yet.
  Raw data is in `SESSION_HANDOFF_2026-05-26_v4.md` and the four JSON
  result files in repo root.
- INFO-016 disambiguation data point exists; reading happens with Greg
  in the next session.
- Branch state: work persisted on `claude/linear-drift-nreal-sweep-cjRkM`,
  pushed. `main` untouched. No PR.

Earlier note (Session 3, retained for context):

This CLAUDE.md was updated at the end of 2026-05-25 Session 3 to reflect:

- Retraction of the Family A/B taxonomy (INFO-008 and dependents)
- Addition of the Result Discipline rule
- Addition of the Base-of-Structure heuristic as a working frame
- Addition of the Dipole-Couples reading as a working frame
- Updated ledger with proper status tags (isolated finding / mapped finding / located finding / null finding / retracted)
- New experiments queue starting with linear drift N_REAL sweep
- Greg's call that each domain is its own substrate inquiry — no more pooled cross-domain claims without domain-native bases

The Information Layer line of work is at a methodologically clarified but interpretively narrower point than at the end of Session 2. Multiple frames are alive; none has been promoted to claim.
