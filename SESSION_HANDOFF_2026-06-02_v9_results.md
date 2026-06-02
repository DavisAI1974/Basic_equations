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

## Queued for Session 10

1. Re-examine leg (ii) (EM/strong/weak/gravity expression-level
   functional families) across asymmetry + knob sweeps, once PySR (or a
   symbolic-regression fallback) is available. Test whether the family
   contrast is regime-specific like INFO-030.
2. Standardized-vs-raw covariance fork (INFO-024): re-run the flip
   sweeps with per-column standardization to quantify how much of the
   flip is the raw-covariance procedure vs structural. Decides physical
   weight of the flip.
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
