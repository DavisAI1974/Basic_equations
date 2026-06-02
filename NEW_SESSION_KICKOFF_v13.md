# NEW SESSION KICKOFF v13 — opens Session 13 (continues 2026-06-02 Session 12)

## Attach to the new chat
1. `CLAUDE.md` (master context, updated through Session 12)
2. `SESSION_HANDOFF_2026-06-02_v12_results.md` (read this first)
3. this file

## Branch + environment
- **Work on the branch given at session start** (Session 12 used
  `claude/file-attachment-hold-DjGSW`; Greg's standing rule: use the
  session-start branch, keep `main` synced, disregard any stale branch name in
  older docs). Confirm at start.
- SessionStart hook installs numpy/scipy/sklearn/pysr + **h5py** (added S12) and
  bootstraps Julia. ~29G disk. GWOSC / PDG / PyPI reachable. `data/ligo_bulk/`
  (~1.9GB) is gitignored.
- No PR unless Greg asks.

## State in one breath
The information-dipole work now has a paper connection and a coupling map.
Per-domain (physics/biology/chemistry/geology) governing equations are
consolidated (`od_per_domain_equations.json`). The flow-dipole paper
(davisai.ai/dipole) IS the operator family the windowed-null extraction works
on (INFO-039). Coupling is PER-DOMAIN and heterogeneous: biology is genuinely
MI-coupled (knob-confirmed, INFO-040), chemistry has a stable residual,
physics/geology are pure equal-entropy; generic cross-science coupling does NOT
reproduce biology's law-like signature (INFO-041); domains split self-pole (3)
vs cross-pole (biology only), complementary not oppositional (INFO-043). The
four-force line: LIGO 12-event null done, Track B inverse problem done, SM-
parameter-regularity hunt done (real structure -- lepton Koide, GST, QLC,
Wolfenstein -- but all conjecture, INFO-042). Force<->equation has NO real-data
bridge yet (the EM<->physics hit is caricature-contaminated).

## First actions (Greg's pick; he flagged the dipole/force thread to dive into)
1. **HEADLINE — principled force-operator-space (the force<->equation dive).**
   DECISION GATE FIRST (Result Discipline, do NOT build before this): is there a
   REAL dataset that gives a gauge force a 2-channel entropy/MI object WITHOUT
   inventing the coupling? Candidates: collider event / cross-section / decay-
   rate distributions vs energy as the two channels (gravity already has its
   real object via LIGO). If yes -> put each force in the per-domain 6-op basis
   and test whether EM really resembles the physics equation (vs caricature
   artifact) and whether the coupling-type mapping (gravity <-> equal-entropy
   self-domains) survives. If no real dataset -> force<->equation stays a frame,
   say so, move on. Map alternatives before committing compute.
2. **Coupling-thread knob tests** (cheap, high-value): does chemistry's residual
   `0.54(H_a+H_b)+0.32H_a^2-0.55H_b^2` track the Brusselator B parameter the way
   biology's MI-slope tracks the Lotka-Volterra coupling g (INFO-040)? Map
   biology's slope-vs-g as a clean coupling-strength readout.
3. **LIGO loose end**: re-run the s11_first_run_entropy no-MI-basis decomposition
   per event across the full 12 -- does the noise-OFF / event-ON attractor
   reversal (seen on GW150914) hold or flip across events?

## Standing constraints
- Do NOT overwrite the Markets section placeholder in CLAUDE.md; ADD to it (the
  S12 dipole subsection with all equations + INFO-039..043 already lives there).
- Medical OD stores (cardiac/cerebro) stay on their own branch, not main.
- All Operating Rules in force: no pre-assigned meaning, probe-not-falsifier,
  speaking posture before+after, Rule D incomplete-not-wrong, they-never-stacked,
  treat-literature-as-conjecture. Report data / interpretation / frame
  separately. No tent-widening on outliers. >=3 seeds + scatter for any
  operator-space spatial claim. Frames never grade themselves (the caricature
  retirement and the INFO-039 deflationary caveat are both this rule in action).

## Ledger pointer (Session 12 new entries)
INFO-039 (dipole-paper connection, MAPPED/deflationary), INFO-040 (per-domain
coupling: biology genuinely coupled, LOCATED), INFO-041 (pairwise Level-2:
generic coupling != law-like coupling; no universal Level-2 dipole), INFO-042
(SM parameter-regularity hunt + force<->equation answer), INFO-043 (cross-domain
self/cross balance is complementary, not oppositional). Track A 12-event LIGO
null complete. Track B inverse-problem complete.
