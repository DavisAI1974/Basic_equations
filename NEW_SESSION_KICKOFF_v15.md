# NEW SESSION KICKOFF v15 — opens Session 15 (continues 2026-06-02 Session 14)

## Attach to the new chat
1. `CLAUDE.md` (master context, updated through Session 14 / INFO-050)
2. `SESSION_HANDOFF_2026-06-02_v14_results.md` (read this first)
3. this file

## Branch + environment
- **Work on the branch given at session start** (confirm at start). S14 used
  `claude/upload-to-memory-1g6eY`. Keep `main` synced if/when Greg asks; disregard
  stale branch names in older docs.
- SessionStart hook installs numpy/scipy/sklearn/pysr + h5py and bootstraps Julia.
  This session also needs: `pip install uproot awkward aiohttp requests` (already
  in requirements.txt) only if re-touching the strong/ATLAS data.
- **Network reality (S14)**: `eospublic.cern.ch` is TLS-BLOCKED this environment
  (proxy can't verify CERN's CA); xrootd:1094 blocked. **Use the bypass**:
  `https://opendata.cern.ch/eos/opendata/<path>` streams CERN open-data files
  (200 + range support). GWOSC / Zenodo / PDG / PyPI reachable. data/strong/*.root
  and large raw downloads are gitignored (re-fetchable).
- No PR unless Greg asks.

## State in one breath
The force<->equation dive's DATA-COLLECTION phase is COMPLETE: all FOUR real
gauge-force operator objects are built and all sit on the equal-entropy SELF-POLE
-- gravity (LIGO), weak (CMS Z->mumu), EM (HBT+g2), strong (ATLAS Pb-Pb
femtoscopy, INFO-050 this session). MI is active but NEVER a null constraint for
any real force; the gauge forces look like the bookkeeping domains
(physics/geology). "MI-in-null = genuine coupling" remains ONLY in simulated
BIOLOGY (knob-confirmed) and -- via its predictive algebraic dipole -- MARKETS.
The grand "info substrate of physics" frame stays pruned; the load-bearing
survivors are the coupling discriminator + per-domain differentiation. (Re-run
the discriminator when models improve -- Markets section flag.)

## First actions (Greg's S14 decision): the "how are the forces DERIVED" 3-piece program
**Do ALL THREE pieces, EASIEST ONE FIRST.** Honest boundary (state it every time):
"what a force IS" (ontology/mechanism/origin) is NOT a data question and stays a
frame; we recover the data-level governing FORM only.

- **PIECE 1a (DO FIRST -- cheapest, data in hand) -- WEAK Z-propagator from raw
  dimuon data.** data/forces/Zmumu.csv (CMS record 545; if not on the branch,
  re-fetch via the opendata.cern.ch bypass). Build the dimuon invariant-mass
  spectrum; recover the Breit-Wigner resonance WITHOUT assuming it -- (i) direct
  fit for M_Z, Gamma_Z; (ii) PySR free-form symbolic regression on the peak to see
  if it independently lands on the resonance/propagator form. Validate recovered
  M_Z, Gamma_Z vs PDG (91.19 GeV, 2.495 GeV). Deliverable: "the weak neutral
  current's data-level propagator shape, derived from raw data." Scope: NOT the
  mechanism.
- **PIECE 1b -- GRAVITY chirp law from the cached LIGO inspiral.** data/ligo_M/ +
  GWOSC strain (the bypass / GWOSC). Extract instantaneous GW frequency f(t)
  through a loud inspiral (e.g. GW150914), and let OD recover the chirp law
  df/dt ~ f^(11/3) (the quadrupole radiation form) without assuming it; check the
  exponent + the chirp-mass relation. Deliverable: "the gravitational radiation
  governing form from raw data."
- **PIECE 3 (+ PIECE 2 folded in) -- unification footprint.** Extend Track B
  (INFO-037, s12_track_b_inverse.py): the new-physics FOOTPRINT (Delta-b_i, onset
  mu_NP) that closes the running-coupling triangle; map alternatives (no-closure
  [deflationary], two-loop, extrapolation-is-conjecture) FIRST. Piece 2 (a
  held-out-PREDICTIVE relation among the ~26 PDG parameters; INFO-042 found
  structure but nothing derived) folds in here. FOOTPRINT recoverable, IDENTITY
  never.

**RUN THE MARKETS ANALOGUE IN LOCKSTEP** (Greg S14: "a lot of it transfers"). The
Piece-1 method (re-derive a governing equation from raw data by symbolic
regression) IS how the markets dipole was found. Apply 1:1 -- Piece 1 = re-derive
the markets dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2 from raw market data;
Piece 2 = where (a,b,c) come from / a predictive relation that survives held-out
data; Piece 3 = markets<->biology unification (the two coupling-side domains).
(Markets-repo JSON pull still blocked -- no list_repos/add_repo; keep the Markets
section updated meanwhile.)

## Standing constraints (Greg, S13-S14)
- **NO synthesis across construction types** -- each real-data object is its own
  probe; clean comparisons are WITHIN type (strong vs weak; EM vs gravity).
- **Keep the CLAUDE.md master context updated in lockstep** with the work,
  INCLUDING the Markets section. Keep the header line current (header-currency
  Operating Rule): bump date + session number whenever the master context updates.
- **Re-run the MI/coupling discriminator when models improve** (Markets re-run flag).
- All Operating Rules in force: no pre-assigned meaning, probe-not-falsifier,
  speaking posture before+after, Rule D incomplete-not-wrong, they-never-stacked,
  treat-literature-as-conjecture (Koide/Cabibbo etc. = conjecture until WE derive
  them), no tent-widening on outliers, >=3 seeds (or the real-data analogue:
  agreement across observables/binnings), frames never grade themselves.

## Ledger pointer (Session 14 new entry)
INFO-050 (STRONG force real-data self-pole, ATLAS DAOD_HION14 Pb-Pb femtoscopy,
MI not in null, robust across charge/observable/seed + low-q + BE-presence checks;
LOCATED/REAL). Completes 4/4 real gauge forces on the self-pole.

## Commercial thread (recorded S13, unchanged)
Core asset = a model-free COUPLING ENGINE (coupled-vs-coincidental, strength,
structural-vs-incidental, sampling-artifact flag) + template-free event detection.
Top fits: MEDICAL (physiological-coupling biomarker / ICU early-warning -> NoVell),
DEFENSE (multi-sensor corroboration / spoof detection -> SENTINEL), FINANCE/ENERGY
(structural-vs-spurious correlation / regime detection -- markets is the one real
non-bio coupler). Patentable core = the MI-in-null discriminator + strength
readout + sparse-data artifact detector. Research-stage; needs per-domain
validation; do not oversell "force unification".
