# Information Layer Session Handoff -- 2026-05-26 Session 6

## Read me first

Read in this order:
1. CLAUDE.md (project root). Now has Session 6 note added at the
   top alongside Sessions 3, 4, 5 notes. All five Operating Rules
   from Sessions 4-5 remain in force (no pre-assigned meaning,
   probe-not-falsifier, speaking posture before and after, Rule D
   incomplete-not-wrong, They never stacked). One methodological
   note added this session, not promoted to a Rule.
2. This file. What Session 6 did, what is confirmed, what is
   refined, what is queued.
3. SESSION_HANDOFF_2026-05-26_v5.md (project root). Session 5
   record. Provides the full setup that Session 6 built on
   (literature scan candidates, three-thread probe pattern,
   per-domain coefficients).
4. SESSION_HANDOFF_2026-05-26_v4.md (project root). Session 4 raw
   data. Optional unless reviewing the chain.
5. SESSION_HANDOFF_2026-05-25.md (project root). Earlier per-domain
   algebraic coefficients (chemistry quadratic, geology rank-3,
   biology MI) -- now independently reproduced by Session 6 (see
   below).

## What Session 6 did

Greg's directive at start of session: "let's do your order".
The order proposed was the four-item menu from the v5 handoff
filtered to "strongest experimental moves" (not literature reads).
Greg added "A then B" at the per-domain juncture. Six experiments
ran. All six pushed to branch
`claude/claude-md-context-update-uCZ8m`.

### Experiment 1 -- KBK pipeline on OU baseline (kbk_pipeline.py)

Stack: Kaiser-Brunton-Kutz 2024 (arXiv:2403.04889) SVD-rank-gap +
symbolic recovery applied to the extract_v1 protocol on the same
OU two-stream baseline that produced v5's (+,+,+) finding.

Result (3 seeds, full T=100 / window=40 / N_REAL=500):
- Singular value spectrum: [80.8, 68.3, 6.81, 3.80, 1.68, 1.22]
- KBK rank-gap: signal_rank=2, null_rank=4, gap_ratio ~10 between
  SV2/SV3 across all seeds
- Smallest extract_v1 null direction (= KBK null[0]) is dominated
  by H_a/H_b (~90% energy); its [2,3,4] projection has cos +0.998
  to (+1,+1,+2)/sqrt(6) (independent confirmation of v5 Probe C)
- KBK null[2] (third-smallest SV) is the clean quadratic-only
  relation -0.42*H_a^2 - 0.43*H_b^2 + 0.80*H_a*H_b ~ 0 with
  residual 0.026, coefficients reproducing to +/-0.01 across seeds.
  This IS the v5 (-1,-1,+2)/sqrt(6) attractor direction in pure
  [2,3,4] form
- KBK null[3] (fourth-smallest) is MI ~ 0 alone (finite-sample
  remnant for uncoupled OU)

Adds vs v5: published-method anchor (citation path), clean symbolic
content per null direction, refinement from "rank 3" to
"rank 4 with internal structure: 2 H-mixed + 1 quadratic-pure
+ 1 MI-pure".

### Experiment 2 -- AI Poincare rank check on OU (ai_poincare_rank.py)

Stack: Cao-Liu-Tegmark 2021 "AI Poincare" (arXiv:2011.04698)
intrinsic-dimension-of-manifold approach. Three estimator
families: two-NN (Facco 2017), Levina-Bickel MLE 2005, local-PCA
rank-vs-scale curve at 95% variance threshold.

Result (3 seeds, same baseline):
- two-NN dim (z-scored):    2.82 - 2.90
- Levina-Bickel MLE k=10:   3.28 - 3.34 +/- ~1.3
- local-PCA 95%-var rank:   3.00 cleanly across k=25 to k=500
                            (median 3, std 0 at most scales)

Intrinsic dim = 3, robust to estimator x neighborhood x seed.
6 ambient - 3 intrinsic = 3 algebraic relations. **Independent
confirmation of v5 INFO-022's 3-relation count.**

Reconciles with experiment 1: KBK rank=4 reading and AI Poincare
rank=3 reading agree once you account for what each measures.
The 4 KBK null directions split as 3 STRONG algebraic relations
(residual 0.007 to 0.026) and 1 WEAK MI~0 constraint (residual
~0.05); AI Poincare at 95% variance treats the weak MI direction
as free, giving intrinsic dim 3.

### Experiment 3 -- Estimator-family sweep + window-size scan (kbk_estimator_sweep.py, window_size_scan.py)

Two follow-up probes on the same OU baseline.

Estimator sweep (canary, single seed, 3 estimator families):
- Vasicek + 12-bin hist MI:  smallest eigval 9.9e-4
- Gaussian KDE H (Silverman) + 12-bin hist MI: smallest eigval 1.1e-3
- k-NN H (k=4) + KSG MI (k=4): smallest eigval 1.2e-3
- All three: KBK signal_rank=2, null_rank=4; AI Poincare local-PCA
  = 3.00 cleanly

Result: v5's machine-epsilon eigenvalues (4e-13 to 3e-11) are NOT
reproduced by any of three estimator families. Eigenvalue floor
is ~1e-3 regardless of estimator. Estimator family is NOT the
cause of the v5/my-run discrepancy.

Window size scan (5 dt values from 0.2 down to 0.01,
N_samples_per_window = 200 to 4000):
- Smallest eigval log-log slope vs N_samples_per_window: -0.058
  (flat)
- 3rd-smallest log-log slope: +0.035 (flat)
- KBK rank stays signal=2/null=4 across all dt
- AI Poincare local-PCA stays 3.00 cleanly across all dt

Result: eigenvalue floor is STRUCTURAL at ~1e-3 to 5e-5, NOT
sampling-noise. Going to denser sampling does not recover v5's
spectrum. Sampling resolution is also NOT the cause.

### Experiment 4 -- SINDy with extended library (sindy_symbolic.py)

Stack: SINDy / SINDyG / DSINDy (Brunton-Kutz family) sparse-null-
direction method with three feature libraries on the same OU
baseline. Per v5 literature scan: "no SINDy/AI-Feynman/PySR/Eureqa
application to windowed differential entropy of coupled species"
exists. This is that stack.

Result (3 seeds, full):
- deg2 (6 features, KBK basis):  smallest SV 1.18-1.27, resid 7e-3
- deg3 (10 features, + H_a^3, H_b^3, H_a^2*H_b, H_a*H_b^2):
  smallest SV 0.087-0.097 (13x smaller), resid 5e-4 (13x smaller)
- deg3+nonpoly (12 features, + log(MI+eps), sqrt(H_a^2+H_b^2)):
  smallest SV 7e-3 (another 13x), resid 4e-5 (another 13x).
  Coefficients reproduce to +/-0.01 across seeds.

Pattern: each library step drops residual ~13x. Reading at the
substantive level: each step captures the next Taylor term in
dH = H - <H>. Set up the sigma-scan hypothesis: smaller sigma ->
smaller dH -> tighter first-order relations -> machine-epsilon
recoverable. (Tested in Experiment 6.)

### Experiment 5 -- Per-domain KBK+AI Poincare+SINDy (per_domain_kbk.py) -- item A

Per v5 handoff Rule D explicitly: per-domain coefficients
(chemistry quadratic, geology rank-3, biology MI) were preserved
through Session 5 and flagged for the same probe discipline. This
is that probe.

Procedure: ensemble-H operator extraction (N_ens=600 at each t)
on 4 domains (physics Duffing, biology Lotka-Volterra, chemistry
Brusselator, geology Burridge-Knopoff). Same KBK+AI Poincare
+SINDy stack. Two seeds (11, 22), full config T=30s dt=0.02
N_ens=600.

Cross-seed reproducibility of smallest extract_v1 null
(cos seed11 vs seed22, per domain):
  Physics:    +0.999
  Biology:    +0.999
  Chemistry:  +0.997
  Geology:    +0.985

Cross-domain cosines stay low (mostly < 0.5, occasionally up to
+0.81). Each domain has its own signature. The (+1,+1,+2)/sqrt(6)
protocol artifact does NOT dominate per-domain (max cos to
(+,+,+) across all 4 domains x 2 seeds: +0.32; minimum -0.28).

**Reproduces every v5 per-domain claim:**

1. **Geology rank-3 collapse** (INFO-008b): KBK signal_rank=3,
   null_rank=3, gap_ratio 14-17 (largest of 4 domains). AI Poincare
   local-PCA at k=100 = 2.98 +/- 0.12 -> intrinsic dim = 3.
   v5's "geology collapses to 3D" confirmed at both rank-gap and
   intrinsic-dim level.

2. **Geology rank-3 relation coefficients** (INFO-008b): the v5
   relation 0.724*H_a*H_b - 0.441*H_b^2 - 0.290*H_a^2 ~ 0 appears
   as null[1] in both seeds (null[0] is the even-tighter linear
   coupling H_a ~ H_b). Coefficients:
     mine seed11: -0.39*H_a^2 - 0.43*H_b^2 + 0.75*H_a*H_b
     mine seed22: -0.41*H_a^2 - 0.42*H_b^2 + 0.81*H_a*H_b
     v5:          -0.29*H_a^2 - 0.44*H_b^2 + 0.72*H_a*H_b
   cos to v5: +0.99 both seeds.

3. **Biology MI ~ polynomial(H_a) with H_b absent** (INFO-008c):
   smallest null in both seeds is -0.27*H_a + 0.96*MI ~ 0; H_b
   coefficient 0.003-0.017 (essentially absent). Exactly v5's
   "H_b absent" finding. AI Poincare intrinsic dim for biology
   = 1-2 (lower than other domains), consistent with MI as tight
   constraint.

4. **Chemistry chemistry-specific signature** (INFO-008a): KBK
   null[0] is -0.3*H_a + 0.78*H_b - 0.43*H_b^2 + 0.33*H_a*H_b
   (different from all other domains; H_b dominates because the
   ensemble means H_a != H_b for Brusselator). SINDy deg-3 finds
   genuine cubic content: H_b^3 ~ -0.22, H_a*H_b^2 ~ +0.25 in
   both seeds. The chemistry-specific cubic terms are NOT a
   Taylor remnant (different magnitude ratio than physics deg-3).

5. **Physics shows the Taylor quadratic identity** (-1,-1,+2)/sqrt(6)
   as tightest constraint. Geology shows the same identity as
   null[1] (secondary). Other domains do not show it cleanly --
   Brusselator's H_a != H_b means breaks the Taylor regime.

### Experiment 6 -- Sigma scan on OU (sigma_scan.py) -- item B

Hypothesis from Experiment 4: smaller sigma -> smaller dH ->
tighter first-order relations -> recover v5's machine-epsilon
spectrum on OU baseline. Predicted scaling: smallest eigval ~
sigma^4 * log^4(sigma).

Result (sigma from 1.0 to 0.01, 5 values):
- smallest eigval log-log slope vs sigma: +0.56 (NOT +4)
- 3rd-smallest log-log slope: +0.04 (FLAT across 5 orders of
  magnitude in sigma)
- op_std for H_a, H_b stays at 0.14 across all sigma (Vasicek
  estimator variance is sample-count limited, not sigma limited)
- Extrapolation to sigma=0.001 gives 3.7e-6, not 1e-13

Hypothesis falsified. Why: the Vasicek estimator variance is set
by sample-count per window, not by underlying sigma. dH is set by
estimator noise, not by physical dH ~ sigma. So algebraic relations
cannot tighten beyond the estimator-noise floor.

**Resolution of v5 machine-epsilon anomaly (joint with item A):**
- v5's 1e-13 floor likely came from per-domain ensemble-H, not OU
  windowed-H. Ensemble averaging reduces estimator noise as
  1/sqrt(N_ens). Geology at N_ens=600 already reached 2e-7
  (Experiment 5). v5 at N_ens 10^4 - 10^5 plausibly reaches 1e-13.
- OU windowed-H is a different procedure with a higher noise
  floor. No choice of sigma, window, dt, or estimator family within
  windowed-H paradigm recovers v5's spectrum.
- The RANK claim (3 algebraic relations) is robust to both
  procedures. Only the eigenvalue MAGNITUDE depends on procedure.

## Methodological note (not promoted to Rule)

I was wrong twice this session about what controls the eigenvalue
floor of the operator covariance:
- First (Experiment 3, window_size_scan): predicted floor would
  scale with N_samples_per_window. It didn't.
- Second (Experiment 6, sigma_scan): predicted floor would scale
  as sigma^4. It scaled as sigma^0.56 (weak; estimator-noise
  dominated).

Both predictions came from Taylor expansion of (dH)^2 -> structural
noise model. Both missed the estimator-noise contribution. The
correct mental model:

  eigval_floor = min(structural_noise_from_dynamics,
                     estimator_noise_from_procedure)

For per-domain ensemble-H (low estimator noise via 1/sqrt(N_ens))
structural noise dominates and can be very small. For OU single-
trajectory windowed-H (high estimator noise from single-window
Vasicek), estimator noise dominates and pins the floor at ~1e-3
regardless of sigma or window.

Practical consequence: when probing for tighter algebraic
relations, the path forward is increasing ensemble size (or using
a lower-noise estimator) NOT decreasing the source noise. This is
a refinement to the v5 understanding: "tighter relations need
tighter noise" is too coarse; "tighter relations need tighter
estimator noise" is more accurate.

This is consistent with Rule C (the data corrected the prediction)
and Rule D (the v5 reading was incomplete-not-wrong about its OWN
procedure; the eigenvalue-magnitude reading does not generalize
across procedures, but the rank reading does).

## Ledger updates

**INFO-022 (Session 5 structural) -- LOCATED FINDING (was isolated)**:
Rank-3 null subspace under H_a~H_b regime in the 6-op basis is
now confirmed by FOUR independent diagnostics: KBK rank-gap (with
the rank-3-strong-plus-1-weak refinement), AI Poincare local-PCA
intrinsic dim, AI Poincare two-NN, Levina-Bickel MLE. Robust to
3 estimator families, 5 dt values across 20x range in samples-
per-window, 3 seeds. The rank-3 reading itself is now strong; the
*magnitude* of the smallest eigenvalues remains procedure-
dependent (see Methodological Note above).

**INFO-008a (chemistry quadratic) -- DATA HOLDS, READING
CONFIRMED via per-domain probe**: The chemistry-specific
quadratic content from v5 is reproduced via KBK null[0] (mixed
H_a-H_b-quadratic shape distinct from other domains) and SINDy
deg-3 with H_b^3 and H_a*H_b^2 coefficients reproducing across
seeds. NOT a Taylor remnant (different magnitude ratio than
physics or geology deg-3). Status: located finding. The exact
v5 formula H_a^2 = polynomial(H_a*H_b) is a specific functional
fit on top of the same data; the KBK procedure gives the linear
algebraic shape of the constraint manifold, the SINDy deg-3
gives the cubic shape. Both consistent at data level.

**INFO-008b (geology rank-3 relation) -- LOCATED FINDING (was
demoted)**: Independently reproduced at cos +0.99 in both seeds
via per-domain ensemble-H + KBK. AI Poincare confirms intrinsic
dim = 3. Eigenvalues reach 2e-7 (approaching machine epsilon
direction predicted by v5). Status: substantive at data level
and at rank-claim level. Reading rehabilitation: the v5 "geology
rank-3" claim was demoted at Session 3 (depending on INFO-008);
Session 6 per-domain probe rehabilitates it.

**INFO-008c (biology MI ~ polynomial(H_a), H_b absent) -- DATA
HOLDS, READING CONFIRMED**: Smallest null direction in both seeds
is -0.27*H_a + 0.96*MI ~ 0 with H_b coefficient 0.003-0.017.
Exactly v5's "H_b absent" finding. Status: located finding. The
functional form (polynomial vs linear) is testable in next session
via GP regression of MI on H_a per domain.

**INFO-023 -- LOCATED FINDING (Session 6, new)**: Per-domain
ensemble-H + KBK+AI Poincare+SINDy stack produces reproducible
domain-specific null directions. Cross-seed cos +0.985 to +0.999
within domain; cross-domain cos mostly < 0.5. The
(+1,+1,+2)/sqrt(6) protocol artifact does NOT dominate per-domain
(max cos +0.32, min -0.28). Tests the "per-domain expression of
substrate" frame; gives it one supporting data point. Methodology
is the stack of KBK 2024 + AI Poincare 2021 + SINDy with extended
library, none of which had been combined for windowed/ensemble-
H of coupled species before (per v5 literature scan).

**INFO-024 -- METHODOLOGICAL (Session 6)**: Eigenvalue floor of
operator covariance is min(structural_noise, estimator_noise).
Different procedures have different floors. OU windowed-H pins at
~1e-3 to 5e-5 regardless of sigma, window, dt, estimator family;
per-domain ensemble-H reaches 2e-7 (geology), plausibly extends
to machine epsilon at larger N_ens. Resolves the v5 machine-
epsilon anomaly as procedure-dependent. Implications for any
future operator-cloud analysis: state which floor regime you are
in before interpreting eigenvalue magnitudes.

## What this session did NOT do (preserved for later)

- GP regression of MI vs H_a per domain (queued; tests v5's
  "MI ~ polynomial(H_a) with R^2 0.66" against flexible nonlinear
  fit; quick, uses sklearn)
- PySR symbolic regression for arbitrary nonlinear forms (queued;
  tests whether per-domain null directions have simpler ratio/
  exponential representations beyond polynomial; heavy install,
  Julia runtime)
- Stress-test per-domain stack on noisy / missing data / shorter
  T (robustness check)
- Apply methodology to real data (NoVell ECG, SENTINEL biomarkers,
  any PhysioNet or financial time series)
- Five-physics-areas substrate hunt (Greg's reorientation from
  Session 5 still queued)
- Investigate the chemistry SINDy-deg-3 cubic terms further -- are
  they capturing reaction-order content per v5's
  "first-order vs second-order chemistry test" queued from
  Session 3?
- Cross-domain comparison of the Taylor identity strength: it is
  null[0] in physics, null[1] in geology, weaker in chemistry,
  weakest in biology. What predicts where it lands? Likely
  H_a-H_b mean separation, but untested.

## Files produced this session (all in repo root)

Scripts:
- kbk_pipeline.py
- ai_poincare_rank.py
- kbk_estimator_sweep.py
- window_size_scan.py
- sindy_symbolic.py
- per_domain_kbk.py
- sigma_scan.py

JSONs:
- kbk_pipeline_canary.json, kbk_pipeline_results.json
- ai_poincare_rank_canary.json, ai_poincare_rank_results.json
- kbk_estimator_sweep_canary.json
- window_size_scan_results.json
- sindy_symbolic_canary.json, sindy_symbolic_results.json
- per_domain_kbk_canary.json,
  per_domain_kbk_results_seed11.json,
  per_domain_kbk_results_seed22.json
- sigma_scan_results.json

Run logs:
- kbk_pipeline_canary.log, kbk_pipeline_run.log
- ai_poincare_rank_canary.log, ai_poincare_rank_run.log
- kbk_estimator_sweep_canary.log
- window_size_scan.log
- sindy_symbolic_canary.log, sindy_symbolic_run.log
- per_domain_kbk_canary.log,
  per_domain_kbk_run_seed11.log,
  per_domain_kbk_run_seed22.log
- sigma_scan.log

Docs:
- SESSION_HANDOFF_2026-05-26_v5.md (Session 5 record, persisted
  from earlier in the same Session-6 chat)
- SESSION_HANDOFF_2026-05-26_v6.md (this file)
- CLAUDE.md updates (Session 6 note)

To be mirrored to E:\information_layer\ and
F:\Factory\knowledge\information_layer\ per the Operating Rules.

## Implications and applications discussed (substantive reading)

Greg asked late in the session: if the substantive reading holds,
what are the implications across fields? Conditional on the per-
domain finding NOT being a deflationary "fancier statistics"
result, the core claim is: any system with two coupled measurable
subsystems has a reproducible domain-specific algebraic signature
in entropy-operator space, usable for regime classification,
anomaly detection, early warning, and substrate inference.

Application surface mapped (with falsification tests per
application -- see chat for specifics):
- Medicine: NoVell extension (per-patient null direction across
  ECG leads); sepsis early warning; anesthesia depth; diabetes
  phenotyping; seizure/sleep/cognitive decline; drug-organ
  coupling. Carl Saab Cleveland Clinic angle for PhysioNet
  cohort testing.
- Defense: SENTINEL extension (per-agent null in coupled
  biomarker streams); ELINT/SIGINT emitter signatures; sonar/
  radar discrimination; cyber intrusion detection; C4ISR fusion;
  DIU PRISM spacecraft anomaly detection.
- Weather: hurricane rapid-intensification precursors; ENSO state
  classification; pollution event signatures; tornado warning.
- Geophysics: earthquake precursors; volcanic forecasting;
  reservoir monitoring. Squarely in OD provisional patent
  territory.
- Ecology: regime-shift early warning; reef bleaching; crop yield.
- Finance/energy: regime classification in coupled asset returns;
  energy market signatures (Greg's home turf); grid stability.
- Industrial: predictive maintenance; process fault detection;
  battery health; nuclear reactor monitoring.
- Chemistry/materials: reaction network identification; drug
  discovery scaffold classification; catalysis screening.
- Astrophysics: stellar classification; pulsar timing; multi-
  messenger transient signatures.
- Foundational: per-domain null IS the "per-domain expression"
  frame from v5. Five-physics-areas substrate hunt uses this
  methodology as one of the tools.

All conditional on the substantive reading holding. Each
application has its own falsification test (signature does or
does not cluster by class / shift before event / discriminate
regimes on holdout data). Listed in chat in more detail.

## Branch state

- Local: claude/claude-md-context-update-uCZ8m
- Remote: origin/claude/claude-md-context-update-uCZ8m
- Tip (Session 6 close): committed at end of v6 handoff write
- Pushed: yes. No PR. main untouched.

## Prompt for next session

See NEW_SESSION_KICKOFF_v7.md for the drop-in.

End of v6 handoff.
