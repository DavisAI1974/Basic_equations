Resuming Information Layer work, Session 10.

================================================================
BRANCH: claude/gravity-substrate-config-51cfp
START THIS SESSION FROM THAT BRANCH (not main).
Two reasons: (1) the Session 9 work (correction + PySR + hook) lives
there, main is untouched; (2) the SessionStart hook that installs
PySR + Julia is on that branch only -- start elsewhere and PySR will
not be available until reinstalled. No PR has been opened.
================================================================

Before doing anything else, read in this order:
1. CLAUDE.md (project root) -- now has the "Note (Session 9 update --
   2026-06-02)" at the top of the notes, plus the new "Markets / Refrag
   Workspace (placeholder)" section. Six Operating Rules from Sessions
   4-7 remain in force (no pre-assigned meaning; probe-not-falsifier;
   speaking posture before/after; Rule D incomplete-not-wrong; They
   never stacked; Treat literature as conjecture by default). NO new
   Rule in Session 9.
2. SESSION_HANDOFF_2026-06-02_v9_results.md -- full Session 9 record
   (double-check, characterization, std-vs-raw, expression refit, real-
   PySR adjudication, four-force frame reassessment, queue).
3. SESSION_HANDOFF_2026-05-26_v8.md -- the Session 8 record Session 9
   corrected.

Speaking posture (Rule C) applies from the first response: "I think X
might happen, but we'll wait on what the data says and where it points
us." No verdict in advance. Rule D and "literature as conjecture"
especially load-bearing for the real-data probes queued below.

------------------------------------------------------------
WHAT SESSION 9 ESTABLISHED (we are NOT at zero)
------------------------------------------------------------
Greg's directive: do not assume our own prior outputs are correct;
double-check, but do not reconstruct for no reason. Result -- a
correction to the Session 8 close-out gravity result, using the
ORIGINAL Session 8 code (pulled onto this branch), not a rebuild.

- INFO-033 (LOCATED, 3 seeds): the MI-dominant substrate flip is NOT
  gravity-specific. It is an asymmetry-THRESHOLD effect ALL coupled
  caricatures undergo. Flip threshold (omega2/omega1 at |MI coef| ~0.5):
  gravity ~1.1x < EM ~1.3x < weak ~1.6x < strong NEVER (flat to 4x).
  Session 8 tested EM at 1.2x (just below its threshold) and read the
  absence as a gravity/EM category difference. Adding gravity's E_total
  term to EM/strong RAISES their thresholds (suppresses the flip) -- so
  energy-mediation is NOT the cause. The flip tracks coupling-vs-
  detuning: strong's confining cubic never lets go; gravity's weak
  softened coupling lets go first.
- INFO-031 / INFO-032 RE-TAGGED INCOMPLETE (Rule D): data reproduces
  exactly; the "gravity is special at the substrate level" reading is
  removed. Gravity sits at the easy-to-flip end of a coupling-strength
  dial, not a separate category.
- Follow-up A (std-vs-raw, INFO-024 fork): the flip is PARTLY procedure-
  inflated, PARTLY structural. RAW gives |MI|~0.99; per-column
  standardized attenuates to ~0.65 (still past the 0.5 flip line).
  Strong never flips under either procedure.
- Follow-up B then B' (curve-fit fallback, then REAL PySR once
  installed): EM = (H_a-H_b)^2 + const ROBUSTLY confirmed; strong =
  exp((H_b-H_a)) reproduced at one seed but seed-unstable (INFO-028) and
  weak-signal. The fallback under-reported strong (tool limitation, Rule
  D on our own work). INFO-031 EM/strong expression contrast STANDS.

FOUR-FORCE FRAME after Session 9:
  leg (i) shared substrate (symmetric regime) -- intact, better
    explained (the shared substrate IS the symmetric-channel regime).
  leg (ii) expression-level family differences (EM polynomial vs strong
    exponential) -- survives real-PySR re-exam, with INFO-028 seed-
    instability caveat on strong.
  leg (iii) gravity-special SUBSTRATE -- REMOVED.
Net: simpler frame -- one shared substrate + coupling-strength flip-
thresholds + genuine expression-level differences. Still a frame, not a
claim. Nothing has touched real physics data yet.

Methodological note added (not a Rule): a single-point probe can read a
threshold crossing as a categorical property. Sweep the knob before
promoting an "X does this, Y does not" contrast.

------------------------------------------------------------
ENVIRONMENT
------------------------------------------------------------
- PySR 1.5.10 + Julia backend now install automatically via the
  SessionStart hook (.claude/hooks/session-start.sh). Verified working
  in Session 9 (real fit recovered (x1-x0)^2+0.3 in 17s). If the hook
  ran, PySR is ready; if PySR is missing, run the hook manually:
  CLAUDE_CODE_REMOTE=true ./.claude/hooks/session-start.sh
- Pure-numpy substrate scripts (KBK extraction, flip sweeps) never need
  PySR. PySR is only for expression-level functional-family work.
- Original Session 8 modules are on this branch (kbk_pipeline.py,
  per_domain_kbk.py, gravity_glance.py, em_strong_glance.py,
  four_force_probe.py, etc.) -- reuse them, do not reconstruct.

------------------------------------------------------------
MENU FOR SESSION 10 (Greg's call -- pick from here or anywhere)
------------------------------------------------------------
(1) KBK across mapping-campaign knob values (v9 kickoff item 7, still
    open): does the SUBSTRATE-level signature drift with knob value, or
    is only the expression-level family regime-dependent? Tests whether
    substrate is more invariant than expression. Pure-numpy, runs now.
(2) Strong exponential seed-rate: more seeds on strong_asym 0.5/0.7 to
    quantify how often exp((H_b-H_a)) vs alternative forms wins (the
    INFO-028 symmetry-breaking rate). PySR, runs now.
(3) Coupling-strength dial: quantify each caricature's effective
    channel-binding stiffness and check it predicts the INFO-033 flip-
    threshold ordering. Pure-numpy.
(4) Third-estimator triangulation on the flip (kNN/KSG MI alongside
    histogram) to further separate procedure vs structure (extends
    follow-up A).
(5) REAL-DATA probes (the gating tests; need data staged + network
    checked first): storm/waves (NOAA buoy / HRRR), four-force PDG
    coupling-running, LIGO O3/O4 strain, PhysioNet ECG. Each is a fresh
    stack on previously-unstacked data ("they never stacked").
(6) Five-physics-areas substrate hunt (Greg's Session 5 reorientation,
    still open).
(7) Markets / Refrag: decide whether those sessions need PySR+Julia and,
    if extending to symbolic regression of the dipole equations, mirror
    the SessionStart hook to the Markets repos (DavisAI1974/agent +
    Markets). See the Markets placeholder section in CLAUDE.md. The full
    Markets section needs restoring from the E:\refrag workspace file /
    git history commit baa542d.

Hold all six Operating Rules through Session 10. We are still pioneers;
the literature-as-conjecture Rule is load-bearing for the real-data
probes -- engage published coupling-running / emergent-gravity claims on
roughly equal footing with the substrate-vs-expression frame.
