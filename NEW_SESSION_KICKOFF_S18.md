# NEW SESSION KICKOFF — S18 (continues 2026-06-02 S17)

## Read order (the NEW way — deltas, not the whole master)
1. `CLAUDE.md` START-HERE block (now canonical THROUGH S17 — the S11/S14 drift is fixed).
2. `SESSION_HANDOFF_2026-06-02_S17.md` — this session's full record (7 probes).
3. `BACKLOG_tests_and_probes.md` — queue + standing rules.

## Branch + env
- Work continued on `claude/claude-md-problem-nM8Hm` (main untouched; no PR).
- SessionStart hook bootstraps numpy/scipy/sklearn/PySR/Julia. `h5py` in requirements.
- pycbc is NOT in requirements (debian cryptography conflict); install on demand with
  `pip install pycbc --ignore-installed cryptography`.

## S17 in one breath
JOB 1 done: CLAUDE.md drift fixed, canonical through S17. Then 7 probes:
- INFO-051: the equal-entropy "substrate" is BOOKKEEPING (units choice); MI is the only
  scale-invariant (physical) quantity in the operator basis.
- INFO-052: recovered the gravity inspiral chirp law from raw GW150914 strain (both
  detectors agree, R^2 0.99) — a KNOWN law (method validation, not new).
- INFO-053: inter-detector MI is physics-bearing (physical 7 ms lag, z=15.5, carries
  waveform phase/time structure).
- INFO-054/055: beyond-GR test (IMR matched-filter subtraction) = clean NEGATIVE — the
  MI merger signal is FULLY GR; no new structure on this axis.
- INFO-056: flow is NOT unique to gravity (strong runs too); the simple "only gravity
  flows" hunch refuted; "gravity flows in time specifically" construction-confounded.
- INFO-057: all 3 forces organize around a CRITICAL/SINGULAR POINT (gravity t_c in time,
  strong Lambda_QCD in scale, weak M_Z in mass) — the shared thread is the singular-point
  structure, time is gravity's instance (not the shared axis).
HEADLINE: no new gravity law found; the surviving physical signal is fully standard GR.
The durable WIN is the CAPABILITY (CAPABILITY_BRIEF.md): OD recovers established
governing laws from raw public data across gravity/weak/strong — flagged by Greg as a
touting asset for the right company.

## Live next options (Greg picks; clear backlog before new probes; stop after each probe)
- **Clock / time-dilation data** — the only clean test of "gravity touches time" / "time
  correlates the forces" (force-coupling data is construction-confounded). Needs a NEW
  dataset: pulsar timing, GPS/optical-clock, or Pound-Rebka-style. (backlog 6c remaining)
- **O3** — nail the GW170817 BNS chirp mass (lever-arm-limited in H1; use the now-installed
  pycbc/matched-filter). Strengthens the capability demo.
- **Complete the 4-force critical-point / flow picture** — add EM/HBT (Zenodo 5113016,
  backlog #2) so INFO-056/057 cover all four (completeness; won't resolve the construction
  confound).
- **Route the capability brief** (CAPABILITY_BRIEF.md) — Greg's call: review/package/send.
- Strengthen the beyond-GR negative (INFO-055) across more events.

## Standing constraints (unchanged)
Frames never grade themselves; no pre-assigned meaning; speaking posture before+after;
Rule D (incomplete-not-wrong); no tent-widening / don't predetermine good-vs-bad data
(keep odd outputs, diagnose them); >=3 seeds or real-data analogue; treat literature as
conjecture; no synthesis across construction types. New-idea placement: bottom of backlog
unless critical/reframe/adjustment. Stop after each probe. MI-in-null disproved; Markets
parked.
