# NEW SESSION KICKOFF v14 — opens Session 14 (continues 2026-06-02 Session 13)

## Attach to the new chat
1. `CLAUDE.md` (master context, updated through Session 13 / INFO-049)
2. `SESSION_HANDOFF_2026-06-02_v13_results.md` (read this first)
3. this file

## Branch + environment
- **Work on the branch given at session start** (Session 13 used
  `claude/file-upload-memory-LKCoB`; keep `main` synced; disregard any stale
  branch name in older docs). Confirm at start.
- SessionStart hook installs numpy/scipy/sklearn/pysr + **h5py** and bootstraps
  Julia. Additional this session: `pip install libarchive-c` (reads RAR via
  system libarchive.so.13). For the big build you will likely add **uproot** +
  **awkward** (pure-Python ROOT reading -- try this BEFORE any VM).
- Network reaches GWOSC / CERN Open Data (opendata.cern.ch) / Zenodo / PDG /
  HEPData / PyPI (all 200 in S13). `data/ligo_bulk/` and `data/forces/em/` are
  gitignored (large, re-fetchable); `data/ligo_M/*.npz` and `data/forces/Zmumu.csv`
  are kept.
- No PR unless Greg asks.

## State in one breath
The force<->equation dive has a validated method and three real gauge-force
objects. ONE model-free pipeline (windowed entropy H_a,H_b + mutual information
MI -> conserved "null" direction) now runs on simulated sciences AND real physics
data. Tonight's real-data result: GRAVITY (LIGO), WEAK (CMS Z->mumu, INFO-047),
and EM (HBT + g2(0), INFO-048/049) all sit on the equal-entropy SELF-POLE -- MI
stays OUT of the null -- so the gauge forces look like the bookkeeping domains
(physics/geology), NOT the coupling pole. The caricature EM<->physics / weak<->
chemistry matches were self-generated and do NOT reproduce on real data.
"MI-in-the-null = genuine dynamical coupling" appears ONLY in simulated biology
(knob-confirmed) and -- via its predictive algebraic dipole -- markets. Chemistry
has a weaker residual coupling (knob-confirmed via Brusselator B, INFO-044). The
grand "info substrate of physics" frame was PRUNED, not advanced; the load-bearing
survivors are the coupling discriminator + per-domain differentiation.

## First action (Greg's pick: "let's do the big one next") -- STRONG FORCE
**Build the strong-force real-data operator object** -- the last within-type force
and the 4-for-4 test of the self-pole frame. If strong ALSO lands self-pole ->
clean 4/4 real-data statement (gauge forces are bookkeeping-side). If it does NOT
-> the most interesting outcome of the whole thread. (No verdict in advance --
Rule C speaking posture.)

DECISION GATE FIRST (Result Discipline -- map before committing TB/compute):
- Data exists: CERN Open Data heavy-ion EVENT-LEVEL -- ALICE Pb-Pb ESD (record
  1102 etc., run 139038 ~1.2 TiB, CC0) or CMS Run-1 PbPb RECO/AOD (the dihadron
  "ridge" substrate). The inter-hadron (two-pion femtoscopy / Bose-Einstein)
  correlation is the intrinsic strong-interaction signal -- not invented.
- HEPData femtoscopy records are SUMMARY-ONLY (precomputed C(q) curves) -> CANNOT
  build per-channel entropies; sanity-check only, not a build.
- The lift is TB data + ROOT tooling. **KEY first step (cheap, do before any TB
  download):** `pip install uproot awkward`, fetch ONE open-data heavy-ion file,
  and test whether uproot can read event-level charged-hadron track branches
  (pt, eta, phi). 
    - If YES -> build two-identical-hadron channels and run the SAME 6-op
      extraction on a SMALL subset (a handful of files, NOT the full TB), binned
      by an energy axis (pair relative-momentum q, or centrality/multiplicity).
    - If the reco objects need the full experiment framework (CMSSW / AliPhysics
      VM) -> that is the genuinely heavy path; report the cost and check with Greg
      before committing, or look for a skimmed/derived event-level dataset first.
- Construction = same family as weak (event-ensemble binned by an energy axis;
  channels = two identical hadrons). Within-type comparison is to WEAK. Report
  data / interpretation / frame separately; the self-pole reading is "MI NOT in
  the null"; deviations are the interesting part.

## Standing constraints (Greg, S13)
- **NO synthesis across objects** -- synthesizing isn't real data, it earns us
  nothing. Each real-data object is its own probe. Clean comparisons are
  WITHIN-construction-type (strong vs weak; EM vs gravity).
- **Don't worry about Markets** (repo pull dropped). BUT keep the CLAUDE.md
  master context updated in lockstep with the work ("update the claude section
  when the rest is updated") -- including the Markets SECTION when dipole/force
  findings touch it.
- Keep the CLAUDE.md header line current (header-currency Operating Rule, S13):
  bump the date + session number whenever the master context is updated.
- All Operating Rules in force: no pre-assigned meaning, probe-not-falsifier,
  speaking posture before+after, Rule D incomplete-not-wrong, they-never-stacked,
  treat-literature-as-conjecture, no tent-widening on outliers, >=3 seeds (or the
  real-data analogue: agreement across observables/binnings) for any operator-
  space spatial claim, frames never grade themselves.

## Ledger pointer (Session 13 new entries)
INFO-044 (chemistry residual tracks Brusselator B, knob-confirmed, LOCATED),
INFO-045 (LIGO MI-merger-departure generalizes 10-11/11; full reversal does not,
LOCATED), INFO-046 (peri-event before/after = detector non-stationarity, not
source physics, LOCATED), INFO-047 (WEAK real-data self-pole, MI not in null,
robust across 5 observables, LOCATED/REAL), INFO-048 (EM HBT self-pole, attractor
= channel symmetry, LOCATED/REAL), INFO-049 (EM g2(0) bunching real but stays out
of the null; coherence-scale MI-in-null = INFO-024 sparse-count artifact,
LOCATED/REAL). Force-operator-space gate PASSED 3/4; strong force is the open
within-type extension (this session).

## Commercial thread (recorded S13, for when Greg wants to act on it)
Core asset = a model-free COUPLING ENGINE (coupled-vs-coincidental, strength,
structural-vs-incidental, sampling-artifact flag) + template-free event detection.
Top fits: MEDICAL (physiological-coupling biomarker / ICU early-warning -> NoVell),
DEFENSE (multi-sensor corroboration / spoof detection -> SENTINEL), FINANCE/ENERGY
(structural-vs-spurious correlation / regime detection). Patentable core = the
MI-in-null discriminator + strength readout + sparse-data artifact detector.
Research-stage; needs per-domain validation; do not oversell "force unification".
