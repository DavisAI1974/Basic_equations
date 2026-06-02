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
  - Leg (ii) (expression-level family differences) WAS re-examined this
    session -- first with a curve-fit fallback (follow-up B, which under-
    reported strong), then with REAL PySR (follow-up B'). Verdict: the
    EM/strong contrast STANDS. EM = (H_a-H_b)^2 + const is robust across
    seeds and across the asymmetry sweep (as a family; coefficient drifts
    with regime per INFO-030). strong = exp((H_b-H_a)) reproduces at one
    seed but is seed-unstable (INFO-028 channel-symmetry breaking) and
    weak-signal (strong's MI is near-constant). So leg (ii) is real, with
    a seed-instability caveat on the strong half.

Net: the four-force frame goes from "three distinct signatures incl. a
gravity-special SUBSTRATE" to "one shared substrate (symmetric regime) +
force-dependent off-substrate flip-thresholds set by coupling strength
(leg iii reinterpreted) + genuine EXPRESSION-level family differences
(leg ii survives real-PySR re-exam)." The gravity-special-substrate
claim is gone; the shared-substrate (leg i) and expression-difference
(leg ii) legs stand. Simpler on the substrate side, and consistent with
the Base-of-Structure heuristic
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

### Follow-up B' -- REAL PySR adjudication (s9_pysr_adjudicate.py)

After follow-ups A/B, PySR was installed in-container (pip install pysr;
the Julia backend auto-bootstrapped -- network policy allowed it). This
let leg (ii) be settled with the ACTUAL tool (original fit_pysr config
reused), which follow-up B's curve-fit fallback could not do. 2 seeds,
N_ens=600, T=30, niter=30.

  - EM = (H_a - H_b)^2 + const: ROBUSTLY CONFIRMED. EM_asym 1.0/1.2 ->
    (x0-x1)^2 + 0.20 BOTH seeds; holds at 1.0/1.5 with shrinking
    coefficient; EM_sym -> (const - channel^2)^2 with channel randomly
    assigned (INFO-028). The fallback (B) and real PySR (B') agree on EM.
  - strong = exp((H_b - H_a)) - const: PARTIALLY REPRODUCED. At strong's
    exact INFO-031 config (0.5/0.7), seed 11 -> exp(-x0 + x1)*0.24 =
    exp(H_b - H_a)*c, which IS the INFO-031 form. Seed 22 -> a different
    shape (x0*(x0-x1)+c, with exp(square(x0)) also on its Pareto):
    SEED-UNSTABLE, exactly as INFO-028 predicts for channel-symmetric
    strong dynamics. So the exponential is in strong's family but not a
    stable per-seed realization. Losses are tiny (~5e-3) and the Pareto
    is flat -- strong's MI is nearly constant, so the family fits a weak
    residual signal.

CORRECTION to follow-up B (Rule D, applied to our own work): B's "strong
has no clean family / exp_diff not reproduced" reading was INCOMPLETE --
it was a limitation of the fixed-library curve_fit (exp_diff was in the
library but the optimizer did not surface it), NOT evidence against
INFO-031. Real PySR surfaces the exponential at one seed. INFO-031's
EM/strong expression-level CONTRAST therefore STANDS: EM polynomial-in-
difference (robust), strong exponential-in-difference (real but weak-
signal and seed-unstable per INFO-028). Leg (ii) of the four-force frame
is in better shape than the follow-up-B reassessment suggested.

Reference configs (real PySR): gravity_asym MI ~ 0.21 near-constant
(both seeds, loss ~7e-4) -- confirms the MI->const floor under asymmetry;
gravity_sym -> (const - H^2) forms; weak_sym -> linear in H_a ~ 0.30.

## Queued for Session 10

1. [DONE this session, follow-ups B + B'] Leg (ii) settled with REAL
   PySR: EM family robustly confirmed; strong exponential reproduced at
   one seed (seed-unstable per INFO-028). INFO-031 contrast stands.
   Open extension: more seeds on strong_asym to quantify how often the
   exponential vs alternative forms wins (the INFO-028 symmetry-breaking
   rate).
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
- s9_std_vs_raw.py          (follow-up A: standardized-vs-raw fork)
- s9_expression_refit.py    (follow-up B: curve-fit fallback family)
- s9_pysr_adjudicate.py     (follow-up B': REAL PySR family adjudication)
- .claude/hooks/session-start.sh, .claude/settings.json, requirements.txt
                            (PySR/numeric-stack persistence hook)
Result JSONs:
- s9_doublecheck_flip_canary.json, s9_doublecheck_flip_results.json
- s9_characterize_canary.json, s9_characterize_results.json
- s9_std_vs_raw_canary.json, s9_std_vs_raw_results.json
- s9_expression_refit_canary.json, s9_expression_refit_results.json
- s9_pysr_adjudicate_canary.json, s9_pysr_adjudicate_results.json
Original Session 8 modules pulled onto this branch (verbatim, for the
double-check): kbk_pipeline.py, per_domain_kbk.py, ai_poincare_rank.py,
sindy_symbolic.py, gravity_glance.py, em_strong_glance.py,
four_force_probe.py.

Dependencies: numpy 2.4.6, scipy 1.17.1, scikit-learn 1.8.0. PySR 1.5.10
was INSTALLED mid-session (pip install pysr; Julia backend auto-
bootstrapped -- the network policy allowed the Julia download). A real
PySR fit runs end-to-end (verified: recovered (x1-x0)^2+0.3 in 17s after
precompile). This unblocked follow-up B'.

PERSISTENCE: PySR/Julia do NOT survive into a fresh container. Added a
SessionStart hook (.claude/hooks/session-start.sh + .claude/settings.json
+ requirements.txt) that installs numpy/scipy/scikit-learn/pysr and
bootstraps Julia at container start (web/remote only; idempotent; PySR
best-effort so a blocked Julia download cannot break session start).
NOTE: the hook only takes effect for future sessions once it is on the
branch the session starts from -- it is on claude/gravity-substrate-
config-51cfp, NOT yet on main. Merge to main (or start the next session
from this branch) for it to run automatically.

Pydroid-3 (Python on Android, on Greg's phone) was raised as an option
but CANNOT run PySR: Julia does not run on Android, and Pydroid is a
separate environment from the cloud container where this work executes.
It could run the pure-numpy substrate scripts only.

## Branch state

- Local + remote: claude/gravity-substrate-config-51cfp. Pushed.
- main untouched. No PR.

End of v9 results handoff.
