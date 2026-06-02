# SESSION HANDOFF v14 — 2026-06-02 (Session 14): STRONG force real-data object (4/4 self-pole) + "how are the forces derived" 3-piece roadmap

Read this first, then the Session 14 note in CLAUDE.md (INFO-050 + the queued
"what ARE the 4 forces / how derived" thread). Branch:
`claude/upload-to-memory-1g6eY` (the three S13 master files were uploaded to
memory at session start; main NOT synced -- the S13 force scripts live on
`claude/file-upload-memory-LKCoB`, this branch carries the S14 work). No PR.

## What happened this session

Greg uploaded the S13 CLAUDE.md + v13 handoff + v14 kickoff to memory, had me
pull the S12/S13 notes into context (the auto-loaded block was stale at S11),
then: "let's do the heavy build" -- the STRONG force, the last within-type gauge
force and the 4-for-4 test of the self-pole frame. After it landed, Greg asked
(1) to flag a Markets re-run when models improve, and (2) "what would we actually
need to figure out what the 4 forces actually are and where/how are they derived
-- even one part." Decision at close: do ALL THREE derivation pieces, easiest
first, in a NEW session; update all 3 continuity files; mirror the roadmap into
the Markets section ("a lot of it transfers").

## THE HEADLINE -- strong force real-data operator object (INFO-050)

Decision gate (Result Discipline, mapped before any TB download):
- **Network**: this session's egress proxy CANNOT reach `eospublic.cern.ch` (the
  EOS host for ALL CMS/ATLAS open-data files) -- HTTPS 503 (proxy can't verify
  CERN's TLS cert chain), xrootd:1094 blocked. DIFFERS from S13. **Reusable
  workaround**: `https://opendata.cern.ch/eos/opendata/<path>` streams the same
  files (200 + HTTP range support), so uproot partial remote reads work.
- **CMS HI RECO REJECTED** (records 14010/14011/14014, 19.3 TB, ~2.5-3.9 GB/file):
  uproot opens the Events tree (2380 branches), sees hiSelectedTracks /
  hiGlobalPrimTracks with momentum_.fCoordinates leaves, but returns 0-length for
  EVERY track member (momentum/chi2_/ndof_/charge_) while `.present`=True --
  a systematic failure to reconstruct the vector<reco::Track> member-wise counts
  (AsGroup / UnknownInterpretation EDProduct). Needs CMSSW (heavy VM). Confirmed
  not-empty-events via .present + mid-file sampling.
- **ATLAS DAOD_HION14 PASSED** (record 80036, 2015 Pb-Pb "Open Data for Research",
  1913 files / 4.4 TB, CC0, /eos/opendata/atlas/rucio/data15_hi/): pure-Python
  uproot reads CollectionTree; InDetTrackParticlesAuxDyn.{phi,theta,qOverP,HITight}
  are flat readable arrays. A real central PbPb event = 3644 tracks, pt 0.5-2.8
  GeV, eta +-2.3, read in 0.3s. NO VM. pt=sin(theta)/|qOverP|, eta=-ln tan(theta/2),
  p=1/|qOverP|.

Build (s14_force_strong.py), SAME family as WEAK (s13_force_dimuon):
- CHANNELS = two identical same-charge hadrons of a pair, A/B assigned AT RANDOM
  (symmetric labeling; identical bosons have no distinguishing charge -- the
  strong analogue of weak's mu+/mu- symmetry). Observable = hadron pT (eta robust).
- ENERGY AXIS = pair relative momentum q_inv (femtoscopy scale, BE peak at q->0),
  equal-count bins. Per q-bin: 6-op [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI] -> extract_v1.
- 1500 events / 5 files; cuts HITight, pt in [0.2,2.0] GeV, |eta|<2.5, MAX_TRK=150.

FINDING: strong sits on the equal-entropy SELF-POLE. MI-coupling = 0.000,
equal-entropy ~ 1.000, ROBUST across charge (++/--) x observable (pT/eta) x 3
random-A/B seeds (real-data analogue of >=3-seed). Two no-tent-widening checks
(s14_strong_lowq_becheck.py):
  (1) BE correlation CONFIRMED PRESENT -- C(q)=same-event/mixed-event for
      same-charge pairs rises at low q, C(q<0.1)~1.10 (lowest bins 1.19-1.23).
      MI-not-in-null is NOT the absence of a signal (the INFO-049 lesson).
  (2) Restricting to the low-q femtoscopy window q<0.4 GeV (down to q~0.075, where
      BE peaks): self-pole SURVIVES, MI-coupling still 0.000 all configs/seeds.
The genuine BE coupling creates MI but it stays an ACTIVE high-variance variable,
never a low-variance null constraint -- exactly like EM/HBT (same Bose statistics,
INFO-048/049) and consistent with INFO-041. CAVEATS (honest): MAX_TRK=150 cap
(central events have thousands -- compute bound; no Coulomb/purity correction, so
modest BE enhancement); 1500 of 4.4 TB events (small subset, but null stable
across files/seeds/observables/q-windows); one construction (femtoscopy pair).

### Force-frame (held as frame, NO synthesis per Greg)
**4/4 real gauge forces -- gravity (LIGO INFO-036/038), weak (CMS Z->mumu
INFO-047), EM (HBT+g2 INFO-048/049), strong (ATLAS femtoscopy INFO-050) -- all
sit on the equal-entropy SELF-POLE; MI active but NEVER a null constraint.** The
gauge forces look like the bookkeeping domains (physics/geology), NOT the coupling
pole. "MI-in-null = law-like coupling" remains confined to simulated BIOLOGY
(knob-confirmed, INFO-040) and -- via its predictive algebraic dipole -- MARKETS.
The clean 4/4 statement the v14 kickoff named as the goal is achieved. Strong did
NOT give the "most interesting" outcome (MI in the null); it confirmed the
self-pole. Compared WITHIN type to weak (both particle-pair event ensembles);
gravity/EM are the other (detector-pair time-series) type.

## Two recorded items (Greg's two asks)

1. **MARKETS re-run flag** (Markets section, S14 block): RE-RUN THE WHOLE
   MI/COUPLING DISCRIMINATOR ONCE WE HAVE BETTER MODELS. The coupling-vs-self-pole
   split rests on toy ODE simulators + the ~0.993 markets predictor; as the models
   get higher-fidelity (real biology/market data, more channels) the discriminator
   must be re-run end-to-end -- does MI still enter the null only for
   biology+markets, do the 4 forces stay self-pole? Load-bearing claim of the
   whole arc; first thing to re-check when models improve (Rule D: one slice).

2. **"What ARE the 4 forces / how derived" roadmap** (queued thread + Markets
   transfer). HONEST BOUNDARY: "what a force IS" (ontology/mechanism/origin) is
   NOT a data question -- OD extracts governing equations, not mechanism, and no
   dataset contains the origin. That stays a frame. But "where/how derived"
   decomposes into 3 data-shaped falsifiable pieces:
   - PIECE 1: re-derive each force's GOVERNING EQUATION from raw data without
     assuming it (OD mantra). gravity/EM have classical force laws (symbolic
     regression precedent: Lemos-Cranmer 2022 rediscovered Newton from real
     ephemerides). weak/strong have no classical law -> recover the propagator /
     running-coupling / resonance shape from measured distributions. Data IN HAND:
     weak Z Breit-Wigner from Zmumu.csv; strong alpha_s/femtoscopy R from data/strong/.
   - PIECE 2: where the COUPLING STRENGTHS come from -- a PREDICTIVE relation among
     the ~26 SM parameters (held-out-predictive). INFO-042 found structure (Koide,
     Cabibbo) but nothing derived. PDG set in hand. High-risk (likely no clean rule).
   - PIECE 3: are the 4 ONE thing -- Track B inverse problem (INFO-037, partly
     built): the new-physics FOOTPRINT (Delta-b_i, mu_NP) that closes the running
     triangle. FOOTPRINT recoverable, IDENTITY never.
   TRANSFERS TO MARKETS (Greg): the Piece-1 method IS how the markets dipole was
   found; the 3 pieces apply to markets one-to-one (re-derive the dipole eqn /
   where coefficients come from / markets<->biology unification). Run in lockstep.

## State of repo / runs / data
- Committed + pushed on `claude/upload-to-memory-1g6eY`:
  the 3 master files (CLAUDE.md S14, v13 handoff, v14 kickoff) uploaded to memory;
  CLAUDE.md (header S14, INFO-050, S14 note, Markets re-run flag, forces roadmap +
  markets transfer, handoff pointer); requirements.txt (+uproot/awkward/aiohttp/
  requests/h5py); s14_force_strong.py(+results+canary),
  s14_strong_lowq_becheck.py(+results); this handoff; the v15 kickoff.
- Data: data/strong/*.root (5 ATLAS DAOD files, ~1.2 GB) are GITIGNORED
  (re-fetchable via the opendata.cern.ch streaming bypass); data/forces/Zmumu.csv
  and data/ligo_M/ are needed for the S15 pieces (Zmumu.csv lives on
  claude/file-upload-memory-LKCoB / record 545 -- re-fetch via the bypass if not
  on this branch).
- Env: uproot 5.7.4 / awkward 2.9.0 / aiohttp installed (in requirements). xrootd
  client installed but xrootd:1094 is blocked this session -- use the
  opendata.cern.ch HTTP-streaming bypass. eospublic.cern.ch is TLS-blocked this
  session (was reachable S13 -- environment-level network-policy difference).
  PySR + Julia bootstrapped by the SessionStart hook.

## NEXT SESSION (Greg S14): the 3-piece derivation program, EASIEST FIRST
See NEW_SESSION_KICKOFF_v15.md. Order: (1) weak Z-propagator from Zmumu.csv
(cheapest, pure-fit + PySR Breit-Wigner recover, check M_Z/Gamma_Z vs PDG); (2)
gravity chirp law from cached LIGO inspiral (f(t) -> df/dt ~ f^(11/3)); (3)
unification footprint (extend Track B / s12_track_b_inverse.py; Piece 2 folds in).
Run the MARKETS analogue in lockstep. Scope honesty on every piece: we recover the
data-level governing FORM, not the mechanism / "what the force is."

All Operating Rules in force (header-currency + Greg's no-synthesis directive).
No new formal Rule this session. No PR.
