# SESSION HANDOFF v13 — 2026-06-02 (Session 13): principled force-operator-space (real gauge-force objects), chemistry-residual + LIGO follow-ups

Read this first, then the Session 13 note in CLAUDE.md (INFO-044 through INFO-049).
Branch: `claude/file-upload-memory-LKCoB` (the session-start branch; main synced;
the v12 scripts/result JSONs were pulled onto it from
`claude/file-attachment-hold-DjGSW` so the runs had their inputs). No PR.

## What happened this session

Greg uploaded the S12 CLAUDE.md to memory, added a header-currency Operating Rule
(the title line had drifted to "Session 7"), then directed: do the three v13
queued items, the two cheap ones first, then the headline. Then "do everything
but the big heavy build" (skip strong force) and "no synthesis across objects --
that isn't real data, each is its own real-data probe." Markets-repo pull dropped
("don't worry about markets; just keep the CLAUDE.md section updated").

### Two cheap follow-ups
1. **INFO-044 -- chemistry residual knob test** (`s13_chemistry_residual.py`).
   Sweep Brusselator B across the Hopf bifurcation (B_crit=2). In the oscillatory
   regime (B>=3) the null[0] residual FRACTION rises monotonically (0.17 -> 0.34
   -> 0.40) while the residual DIRECTION stays pinned to the INFO-040 baseline
   relation (|cos| 0.97-0.99, xseed 0.999+). Chemistry's coupling-strength readout
   = residual fraction of a FIXED non-MI relation (vs biology's MI-slope). Flips
   at the Hopf threshold B=2 (critical point, inspected not absorbed); below
   threshold noise-limited. Refines INFO-040; per-domain heterogeneous coupling
   strengthened (biology MI-coupled + chemistry residual-coupled, both knob-
   confirmed; physics/geology pure equal-entropy).
2. **INFO-045 -- LIGO no-MI decomposition across all 12 events**
   (`s13_ligo_nomi_batch.py`). The event window leaves the 6-op equal-entropy
   attractor in 11/11; removing MI restores it in 10/11 -> the MI-driven merger
   departure (S11 GW150914 thread) GENERALIZES. The full noise-OFF/event-ON
   reversal does NOT (6/11; noise-side is event-dependent). The 3 events whose
   noise also leaves the attractor are exactly the high |H_a-H_b| events (confirms
   INFO-038 inverse relation). S11 open thread resolved.
3. **INFO-046 -- peri-event before/after trajectory** (`s13_ligo_trajectory.py`;
   per-event operator matrices cached to `data/ligo_M/`). Greg's "does anything
   right before/after the merger tell us something?" (a) MI excess is a SPIKE
   confined to the merger window (+/-0.125s); NO inspiral precursor, NO ringdown
   tail -- even the long-inspiral BNS GW170817 (peak only 1.05x far). (b) the
   no-MI attractor IS time-asymmetric pre vs post BUT event-specific in sign and
   lives in off-merger NOISE windows -> DEFLATIONARY: non-stationary detector
   noise floor, NOT source inspiral/ringdown. "Before differs from after" = the
   instrument's noise state, not source physics.

### HEADLINE -- principled force-operator-space (the force<->equation dive)
Decision gate (S12): is there a REAL dataset giving a gauge force a 2-channel
entropy/MI object WITHOUT inventing the coupling? A 3-agent open-data scan
answered: **PASS for 3 of 4 forces.**
- Gravity: already have LIGO (INFO-036/038; data cached `data/ligo_M/`).
- Weak: CMS Open Data dimuon (record 545 Zmumu_Run2011A) -- per-event 2-lepton
  4-vectors, ~MB CSV, direct fetch. Z propagator -> intrinsic correlation.
- EM: Zenodo 5113016 HBT raw two-channel photon timetags (~21MB .rar; read with
  pip `libarchive-c` via system libarchive.so.13). Bose bunching intrinsic.
- Strong: ALICE/CMS heavy-ion femtoscopy event-level data EXISTS but TB-scale +
  ROOT/VM -> DEFERRED by Greg (the "big heavy build", next session).

STRUCTURAL CAVEAT (a finding): the four objects are TWO construction types --
detector-pair time series (gravity + EM) vs particle-pair event ensembles binned
by energy (weak + strong). Clean tests are WITHIN-type. No synthesis across types
(Greg's rule).

- **INFO-047 (WEAK, real CMS Z->mumu, 10227 events; s13_force_dimuon.py +
  s13_force_dimuon_observables.py)**. Channels = mu+/mu- (by charge, symmetric);
  energy axis = dimuon mass M, equal-count bins. Weak object sits on the
  equal-marginal-entropy SELF-POLE (H_a~=H_b, MI NOT in null = 0.000) -- like
  physics/geology, NOT biology. ROBUST across pt/E/signed-pz/signed-eta/CS-cos
  (5 observables = real-data analogue of cross-seed); even cos(theta*) (muons
  anti-correlated, MI=2.6) keeps MI out of the null (INFO-041 mechanism). The
  caricature "weak~chemistry linear-in-H_a" hit does NOT reproduce on real data.
- **INFO-048 (EM, real HBT, Zenodo 5113016; s13_force_hbt.py)**. Same construction
  type as LIGO. MI does NOT enter the null (max 0.031); equal-entropy attractor
  membership tracks the marginal-entropy SYMMETRY of the two channels, not the
  force (asym split detector sits off; symmetric pair sits closer). Confirms
  INFO-036 "attractor = equal-entropy geometry" on a THIRD construction.
- **INFO-049 (EM g2(0) closure; s13_force_hbt_g2.py)**. Direct g2(tau): g2(0)=
  1.83/1.92, coherence ~2us, decays to 1.00 by +-400us -> HBT bunching CONFIRMED
  real. Well-sampled (dt=50us): MI substantial (tracks g2-1), OUT of null, self-
  pole (reproduces INFO-048). Sparse coherence-scale (dt=1-2us, count<1) shows a
  spurious "+1.00*MI~0" null -- a SPARSE-COUNT ESTIMATOR-FLOOR ARTIFACT (INFO-024
  on real data: MI std 0.0004 vs H_a std 0.114), NOT coupling. Genuine bunching
  does not enter the null. Biology remains the ONLY genuine MI-in-null coupling.

### Frame (held as frame, no synthesis)
Across gravity + weak + EM (real data), the gauge forces look like the SELF-POLE/
bookkeeping domains (physics/geology) -- MI active, never a null constraint;
equal-entropy = pure channel symmetry. The caricature EM<->physics / weak<->
chemistry hits do NOT reproduce on real data. "MI-in-the-null = law-like
coupling" appears so far ONLY in simulated biology (and markets, via its
predictive algebraic dipole, sits on the same coupling side). The grand "info
substrate of physics" frame was PRUNED tonight (attractor-as-law, gravity-special
leg, caricature matches all removed), not advanced -- what survives is load-
bearing: a validated coupling discriminator + per-domain differentiation.

### Plain-English synthesis + commercial discussion (recorded for continuity)
Greg asked for a plain-English readout and commercial applications. Core asset =
a MODEL-FREE COUPLING ENGINE: from 2+ data streams, answers (1) genuinely coupled
vs coincidental, (2) how strongly, (3) structural vs incidental, (4) is an
apparent signal a sampling artifact; plus template-free event detection (the
LIGO MI-spike). Top commercial fits (grounded, research-stage IP): MEDICAL/
digital-health (physiological-coupling biomarker + ICU early-warning; best fit --
biology is the one real coupler, slots into NoVell), DEFENSE (multi-sensor
corroboration / spoof detection / template-free detection -- SENTINEL),
FINANCE/ENERGY (structural-vs-spurious correlation, regime/contagion detection --
markets is the one real-coupling non-bio domain). Adjacent: industrial predictive
maintenance (the sparse-data artifact detector is a differentiator), pharma/
systems-biology (coupling-strength + network inference), climate/geophysics,
scientific instrumentation (HBT/LIGO common-signal extraction, g2 from timetags).
Patentable core = the MI-in-null discriminator + strength readout + sparse-data
artifact detector. Honest: method IP, needs per-domain validation; do NOT oversell
"unifying the forces" (tonight pushed against the strong version).

## State of repo / runs / data
- Committed + pushed on `claude/file-upload-memory-LKCoB`:
  CLAUDE.md (through INFO-049, header S13 + header-currency Operating Rule),
  the v12 scripts/JSONs pulled over, requirements.txt (+h5py), and S13:
  s13_chemistry_residual.py(+results+canary), s13_ligo_nomi_batch.py(+results+
  canary), s13_ligo_trajectory.py(+results+canary), s13_force_dimuon.py(+results),
  s13_force_dimuon_observables.py(+results), s13_force_hbt.py(+results),
  s13_force_hbt_g2.py(+results), this handoff, the v12 handoff + v13 kickoff.
- Data: `data/ligo_M/*.npz` (cached per-event operator matrices, small, kept);
  `data/forces/Zmumu.csv` (~1MB, kept); `data/ligo_bulk/` and `data/forces/em/`
  are gitignored (large raw downloads, re-fetchable).
- Env: h5py installed + in requirements; pip `libarchive-c` reads the RAR via
  system libarchive.so.13; network reaches GWOSC / CERN Open Data / Zenodo / PDG
  / HEPData / PyPI (all 200). PySR+Julia bootstrapped by the SessionStart hook.

## NEXT SESSION -- the big heavy build (Greg: "let's do the big one next")
**STRONG FORCE real-data operator object** -- the last within-type force, and the
4-for-4 test of the self-pole frame. See the v14 kickoff for the scoped plan.
Headline question: if strong ALSO lands self-pole, that's a clean 4/4 real-data
statement that the gauge forces are bookkeeping-side; if it does NOT, that is the
most interesting outcome of the whole thread.
- Data: CERN Open Data heavy-ion event-level -- ALICE Pb-Pb ESD (record 1102 etc.,
  run 139038 ~1.2TiB, CC0) or CMS Run-1 PbPb RECO/AOD (the dihadron "ridge"
  substrate). HEPData femtoscopy = summary C(q) only (CANNOT build per-channel
  entropies; sanity-check only).
- The lift = TB data + ROOT tooling. KEY scoping move (test FIRST, before any TB
  download): try pure-Python `uproot` on ONE open-data file to read event-level
  charged-hadron track lists (pt,eta,phi). If uproot reads the track branches,
  build two-pion channels + the same 6-op M-binned extraction on a SMALL subset
  (a handful of files, not the full TB). If reco objects need the full experiment
  framework, fall back to the VM/container (heavy) or a skimmed derived dataset.
  Map this before committing compute/disk (Result Discipline decision gate).
- Construction: same as weak (event-ensemble binned by an energy axis -- e.g.
  pair relative-momentum q or centrality/multiplicity -- channels = two identical
  hadrons). Within-type comparison is to weak (the other particle-pair object).
  NO synthesis across construction types.

## Other open threads
- Markets-repo JSON pull: DROPPED per Greg ("don't worry about markets"); the
  markets algebraic dipole FORM + its coupling-side placement are in CLAUDE.md.
  (Still blocked anyway: no list_repos/add_repo, GitHub locked to basic_equations.)
- EM g2(0) is closed (INFO-049); no further EM needed unless a new angle appears.

All Operating Rules in force (incl. the new header-currency rule + Greg's
no-synthesis directive this session). No new formal Rule. No PR.
