# SESSION KICKOFF — Gravity / Flow-Dipole / Time thread (continues 2026-06-02 Session 15)

## Attach to the new chat
1. `CLAUDE.md` (master, Session 14)
2. `CLAUDE_session_note_2026-06-02.md` (this thread's single-day note — read first)
3. this file

(Deeper probe detail if needed: `SESSION_HANDOFF_2026-06-02_FULL.md` +
`..._gravity_flow_time_EXPLORATORY.md`.)

## Branch + env
- Work on the branch given at session start. Prior work is on
  `claude/file-upload-session-context-pOsLp`; pull the 5 probe scripts + result
  JSONs + `data/ligo/` from there if not present.
- `h5py` needed for the LIGO probes (`pip install h5py`; add to requirements).
  SessionStart hook bootstraps numpy/scipy/sklearn/PySR/Julia.
- No PR unless asked.

## State in one breath
Exploration of the info-dipole "flow" form vs gravity/time. The level null is
provably time-blind; the measured flow-dipole axis sits on the equal-entropy
substrate in 4 sim domains (opposition = the equal-entropy identity, INFO-039 by
measurement). A gravity flow-axis displacement off the substrate appears on
GW150914 but is loosely pinned, equals its own detector noise, and does NOT
reproduce across 3 events (cos→substrate 0.78/0.64/0.99). No claim either way.
Operator hunches H-A..H-G logged in the session note (recorded, not graded).

## First actions (in order)
1. **PROBE 1 — dMI/dt rate-DYNAMICS (run first).** The level extraction is
   time-blind, so gravity's measured time-coupling could only surface in the RATE
   dynamics. Build dMI/dt (+ dH/dt) and test for LAGGED/temporal structure (does
   dMI/dt depend on past operator values; autocorrelation; a rate-form with memory),
   gravity (real LIGO event-window) vs sim domains, within-construction. Decision
   gate first: is ~32 event-window samples enough for a rate-dynamics fit, or is
   finer windowing needed? Scope: data-level rate structure, NOT mechanism.
2. **PROBE 2 — construction control.** Measure a non-gravity REAL detector-pair
   object (EM/HBT, INFO-048) with the P3 flow-axis method; on-substrate (~0.9) or
   off (~0.6–0.8)? Separates "gravity off-substrate" from "real-detector objects
   off-substrate." Needs HBT data.
3. **PROBE 3 — louder events.** Fetch GW170814/GW190521 (4096 s / 4096 Hz hdf5
   archives are GWOSC-API-available, ~134 MB each, slice 32 s), re-run
   `grav_flow_crossevent.py` for ≥3 clean events.

Remaining queued (session note §8): biology joint time-shuffle; inspiral→ringdown
trajectory; scattered-coords structure (H-G); INFO-039 ≥3-seed promotion; gravity
vs gauge forces within-construction.

## Standing constraints
- Frames never grade themselves; no pre-assigned meaning; speaking posture
  before+after; Rule D; no tent-widening; ≥3 seeds (or real-data analogue =
  agreement across events/observables); never discard data — test outliers, locate
  them. NO synthesis across construction types (real detector-pair vs particle-pair
  vs sim ensemble) — compare within type.
- Keep this thread's results at DATA level until the operator chooses to interpret.

## Note
This gravity/flow/time line is PARALLEL to the master v15 "how are the forces
derived" 3-piece program (weak Z-propagator first); that work is unchanged and
independent. Run whichever the operator picks.
