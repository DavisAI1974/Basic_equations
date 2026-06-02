# Information Layer Session Handoff -- 2026-06-02 Session 9 (results)

## Read me first

1. CLAUDE.md "Note (Session 9 update — 2026-06-02)" -- the ledger
   update (INFO-033; INFO-031/032 re-tagged incomplete-not-wrong).
2. This file -- Session 9 record and the four-force frame reassessment.
3. SESSION_HANDOFF_2026-05-26_v8.md -- Session 8 record being corrected.
4. NEW_SESSION_KICKOFF_v9.md -- the kickoff that opened this session.

Branch: `claude/gravity-substrate-config-51cfp`. PySR unavailable (no
Julia). All substrate-side work ran the ORIGINAL Session 8 KBK code
(pulled from `claude/two-more-tasks-O6ahb`), not a reconstruction.

## What Session 9 did

Greg's directive: do not assume our own prior outputs are correct;
double-check, but do not reconstruct for no reason. Target: the Session
8 close-out gravity result (INFO-032, "gravity's substrate is
configuration-dependent in a way EM/strong are not").

### The flaw found

INFO-031 and INFO-032 each tested ONE channel-asymmetry value (EM
omega 1.0/1.2; strong 0.5/0.7; gravity 0.8/1.0). A single point cannot
distinguish a categorical property ("gravity flips, EM doesn't") from a
threshold crossing ("everything flips, at different asymmetries").

### Experiment 1 -- double-check (s9_doublecheck_flip.py)

Original simulators + original raw-covariance extraction
(build_ensemble_operator_matrix + extract_v1), asymmetry SWEEP, 3 seeds.

Session 8 baselines reproduced exactly:
  gravity_asym 0.8/1.0  |MI coef| = 0.990   (INFO-032 said +0.99)
  em_asym 1.0/1.2       |MI coef| = 0.271   (INFO-031: no flip)
  strong_asym 0.5/0.7   |MI coef| = 0.008   (INFO-031: no flip)

Flip thresholds (asymmetry omega2/omega1 at which |MI coef| in v_null
crosses ~0.5, i.e. substrate switches from -(H_a-H_b)^2 ~ 0 to MI ~ const):
  gravity : ~1.1x   (0.8/0.9 -> 0.906)        lowest / most flip-prone
  EM      : ~1.3x   (1.0/1.2 0.27 -> 1.0/1.3 0.570)
  weak    : ~1.6x   (1.0/1.3 0.025 -> 1.0/1.7 0.857)
  strong  : never   (flat < 0.06 out to 0.5/2.0 = 4x)
EM sweep is monotonic and clean; canary (1 seed) and full (3 seed)
agree on flip/no-flip at every point.

### Experiment 2 -- characterization + reframed 5a (s9_characterize.py)

3 seeds, mean +/- std reported. Two minimal experimental variants
(force laws byte-identical to originals, clearly labelled): weak with an
omega2 knob (original hardcodes 1.0); EM/strong with gravity's E_total
universal-energy term added (Fadd = G_add * E_total * dx/(dx^2+eps)).

Reframed 5a question: does energy-mediated coupling CAUSE the easy flip?
  - strong + E_total (G_add 0.05, 0.15, 0.40): stays flip-resistant at
    ALL G; adding energy-coupling does not make strong flip.
  - EM + E_total (G_add 0.15): RAISES EM's threshold (1.0/1.3 flips at
    G=0, no longer flips at G=0.15; flip moves to ~1.4x).
  => Energy-coupling SUPPRESSES the flip (raises threshold). It is NOT
     the cause of gravity's low threshold. Opposite of the 5a guess.

### Joint reading (three levels)

A. Data: the MI-dominant flip is an asymmetry-threshold effect all
   coupled caricatures undergo; threshold ordering strong(inf) > weak
   (~1.6x) > EM (~1.3x) > gravity (~1.1x). Adding energy-coupling raises
   thresholds. All reproduced at 3 seeds with tight scatter.

B. Mechanism (deflationary, now well-supported): the flip is a
   variance-crossing. Channel asymmetry detunes the two oscillators ->
   their entropies decorrelate -> the -(H_a-H_b)^2 relation's residual
   rises, while MI settles to a low near-constant floor. The raw-
   covariance null direction picks whichever relation is tighter; at the
   threshold they swap. Coupling strength resists (keeps channels
   correlated): strong's confining cubic never lets go; gravity's weak
   softened coupling lets go first. No force-specific physics.

C. Frame: "gravity is special at the substrate level" (INFO-032) does
   not survive. Gravity sits at the easy-to-flip end of one continuous
   coupling-strength dial, not in a separate category.

## Ledger

- INFO-033 LOCATED (Session 9): asymmetry-threshold flip; see CLAUDE.md.
- INFO-031 re-tagged INCOMPLETE (substrate-stability part; expression-
  level PySR part not re-examined, PySR unavailable).
- INFO-032 re-tagged INCOMPLETE (data reproduces; "gravity-specific /
  different-kind-of-constraint" reading removed).

## Four-force frame reassessment

What the Session 8 four-force narrative had:
  (i)   INFO-027: all four share the [2,3,4] substrate null on the
        Session 3 attractor (channel-correlation identity).
  (ii)  INFO-031: EM vs strong differ at EXPRESSION level (functional
        family) on shared substrate.
  (iii) INFO-032: gravity differs at SUBSTRATE level (flips to MI).

Session 9 removes leg (iii). Re-reading the other two:
  - Leg (i) is now MORE coherent, not less: "all four share the same
    substrate identity (-(H_a-H_b)^2 ~ 0) when channels are symmetric"
    is exactly what INFO-033 predicts -- the shared substrate is the
    symmetric-channel regime, and asymmetry moves any of them off it by
    a threshold set by coupling strength. The shared-substrate finding
    is intact and better explained.
  - Leg (ii) (expression-level family differences) is a PySR result and
    was NOT touched by Session 9 (PySR unavailable). It needs re-
    examination under the same "sweep the knob, don't trust one point"
    discipline before being trusted -- INFO-030 already showed the
    functional families are regime-conditional, so the EM/strong
    expression contrast may also be regime-specific. Flagged.

Net: the four-force frame collapses from "three distinct signatures
incl. a gravity-special substrate" to "one shared substrate (symmetric
regime) + force-dependent off-substrate thresholds set by coupling
strength." Simpler, and consistent with the Base-of-Structure heuristic
(the base should be simple and shared; the differences are how each
force leaves the base under perturbation). Still a frame, not a claim.

## Session 9 follow-ups (done) -- standardized-vs-raw + expression refit

### Follow-up A -- standardized-vs-raw covariance fork (s9_std_vs_raw.py)

Same operator matrix, two extractions: RAW centered SVD (original;
constancy detector) vs per-column STANDARDIZED SVD (correlation
detector; removes absolute scale so a merely small-variance MI is not
favored). 3 seeds.

  config            RAW |MI|   STD |MI|
  gravity 0.8/0.9   0.906      0.644
  gravity 0.8/1.0   0.990      0.665
  EM 1.0/1.5        0.923      0.600
  EM 1.0/2.0        0.996      0.626
  strong (all)      ~0.02      ~0.02

Reading: standardization ATTENUATES but does NOT erase the flip. The
spectacular +0.99 under RAW is partly inflated by the absolute-scale
procedure (MI picked because its variance collapses), but a real
structural core (~0.6-0.65, still above the 0.5 flip line) survives the
scale-controlled procedure. Strong never flips under either procedure.
The INFO-033 threshold ordering holds in both. So the flip is PARTLY
procedure-inflated, PARTLY structural -- both true. (Note: the 1-seed
canary underestimated STD at ~0.14; the 3-seed full gives ~0.65. Trust
the full run.) Generalizes INFO-024: the null DIRECTION, not just the
eigenvalue magnitude, is procedure-dependent when an operator's absolute
variance collapses.

### Follow-up B -- expression-family refit, leg (ii) re-exam (s9_expression_refit.py)

PySR unavailable -> FIXED-LIBRARY scipy curve_fit fallback (NOT PySR;
cannot discover novel forms, but can test whether INFO-031's specific
families best-fit and whether the contrast survives an asymmetry sweep).
Library includes quad_diff (H_a-H_b)^2+c [EM/INFO-031], exp_diff
a*exp(b*(H_b-H_a))+c [strong/INFO-031], exp_Ha [biology/INFO-025], plus
const/linear/single-channel-quadratic baselines. 3 seeds.

  - EM = quad_diff (H_a-H_b)^2 + const: CONFIRMED. EM 1.0/1.2 R^2=0.929,
    identical family all 3 seeds. Holds as the label across asymmetry
    (R^2 0.93 -> 0.85 -> 0.56 for 1.2/1.5/2.0). Symmetric EM (1.0/1.0)
    has NO clean family (best quad_Ha R^2=0.18) -- consistent with
    INFO-028 (symmetric dynamics -> no reproducible family).
  - strong = exp_diff: NOT CONFIRMED. Strong's MI is poorly fit by every
    family in the library (R^2 <= 0.35 across all configs); best picks
    are unstable linear forms (lin_both/lin_Hb/lin_Ha, split across
    seeds at 0.5/1.0). exp_diff never wins despite being in the library.

Reading: INFO-031's EM/strong expression-level CONTRAST is half-
confirmed. The EM half (quad_diff) is real and reproduces robustly in an
independent method. The strong half (exp_diff) does NOT reproduce; strong
may simply have no clean MI-vs-(H_a,H_b) functional family. Caveat: a
fixed-library fallback cannot fully adjudicate -- a PySR re-run is needed
to settle whether strong's exponential was a PySR-specific overfit or a
real form the fallback's optimizer missed. INFO-031 expression-level part
re-tagged: EM family confirmed; strong family unverified / likely
regime-or-method-specific. Leg (ii) is NOT a clean "two distinct
families" result; it is "EM has a family, strong does not clearly."

## Queued for Session 10

1. [DONE this session, follow-up B] Leg (ii) re-exam via curve-fit
   fallback: EM family confirmed, strong family not reproduced. STILL
   NEEDS a PySR re-run to settle the strong=exp_diff question once Julia
   is available -- the fallback cannot adjudicate a form it failed to
   fit.
2. [DONE this session, follow-up A] Standardized-vs-raw fork: flip is
   partly procedure-inflated, partly structural (STD core ~0.65). Could
   extend with a third estimator (kNN/KSG MI) to triangulate.
3. Remaining v9 kickoff items not yet run: 7 (KBK across mapping-campaign
   knob values -- is substrate more invariant than expression?), and the
   real-data probes (storm/waves, PDG running, LIGO, PhysioNet) once
   data is staged and network access confirmed.
4. Coupling-strength dial as the organizing variable: quantify each
   caricature's effective channel-binding stiffness and check it
   predicts the flip threshold ordering.

## Files produced this session (repo root)

Scripts:
- s9_doublecheck_flip.py    (asymmetry sweep, original code)
- s9_characterize.py        (threshold map + reframed 5a E_total test)
Result JSONs:
- s9_doublecheck_flip_canary.json, s9_doublecheck_flip_results.json
- s9_characterize_canary.json, s9_characterize_results.json
Original Session 8 modules pulled onto this branch (verbatim, for the
double-check): kbk_pipeline.py, per_domain_kbk.py, ai_poincare_rank.py,
sindy_symbolic.py, gravity_glance.py, em_strong_glance.py,
four_force_probe.py.

Dependencies installed: numpy 2.4.6, scipy 1.17.1, scikit-learn 1.8.0.
PySR NOT available (no Julia). Substrate-side KBK is pure-numpy and
unaffected; expression-level (PySR) work deferred.

## Branch state

- Local + remote: claude/gravity-substrate-config-51cfp. Pushed.
- main untouched. No PR.

End of v9 results handoff.
