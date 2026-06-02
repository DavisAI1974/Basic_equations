# Information Layer Session Handoff -- 2026-06-02 Session 10 (results)

## Read me first

1. CLAUDE.md "Note (Session 10 update — 2026-06-02, FIRST REAL DATA)" --
   the ledger update (INFO-034 through 037).
2. This file -- Session 10 record: two method-hardening probes, Greg's
   pivot to real data, two four-force real-data sets.
3. SESSION_HANDOFF_2026-06-02_v9_results.md -- Session 9 (the
   correction that set up this session).
4. NEW_SESSION_KICKOFF_v10.md -- the kickoff that opened this session.

Branch: `claude/gravity-substrate-config-51cfp`. main untouched. No PR.

## Arc of the session

Started on the v10 menu. Greg's framing question: "does real data take
care of these other problems? if not get everything working properly
then real data." Answer given: no -- the toy method-hardening items
(estimator robustness, substrate invariance) are PREREQUISITES for
interpreting any real-data result; real data cannot validate the method.
So we began the "get it working properly" sequence -- then it surfaced
that the four-force laws are our own caricatures, and Greg pivoted the
whole session to real data.

## What ran (in order)

### Item 4 -- third-estimator triangulation (INFO-034)

`s10_third_estimator.py`. Follow-up A left the flip "partly procedure-
inflated, partly structural" using ONE histogram MI estimator for both
raw and standardized SVD. This swapped in a genuinely different MI
estimator family (Kraskov-Stoegbauer-Grassberger kNN, estimator 1) on
the same operator matrix, entropy columns held identical.

Result (3 seeds, N_ens=600, T=30): the flip is ESTIMATOR-ROBUST.
- Histogram reproduces follow-up A exactly (gravity 0.8/0.9 raw 0.906/
  std 0.644; 0.8/1.0 0.990/0.665; EM 1.0/1.5 0.923/0.600).
- KSG agrees on the standardized structural core near threshold
  (gravity 0.8/0.9 std 0.644; EM 1.0/1.3 0.435; EM 1.0/1.5 0.597), on
  raw values, on strong-never-flips, and on the threshold ordering.
- Sharpened reading: structural core is a modest ~0.6-0.7 band NEAR the
  threshold; at extreme asymmetry both estimators drop below 0.5 while
  raw rides ~0.99 on collapsing MI variance (procedure inflation).

### Item 1 -- substrate invariance across the knob sweep (INFO-035)

`s10_substrate_invariance.py`. Reused the ORIGINAL parametrized
simulators + KNOB_SWEEPS from `mapping_campaign.py` (pulled verbatim
from the Session 8 branch). For each domain x knob value x seed,
extracted the SUBSTRATE null direction and measured cross-knob |cos|
(invariance) vs the known expression-family drift (INFO-030).

Result (3 seeds; seed-stability 0.95-1.00, so drifts are structural) --
MIXED, partial support for Base-of-Structure with two exceptions:
- biology (beta 0.3-0.8): |cos| 0.978-0.986, stable rank 1, expression
  mutates -> clean support.
- geology (drift 0.02-0.08): |cos| 0.981-1.000, stable rank 3 (weaker
  test, expression also robust).
- physics (K 0.05-0.5): invariant baseline->high (0.995) but breaks at
  weak coupling K=0.05 (rank 1->4, |cos| 0.69) -- regime-bounded.
- chemistry (B 2-4): NOT invariant, |cos| 0.045 between B=3 and B=4,
  rank swings 4->1->5 -- FAILS. B=2 = Brusselator Hopf threshold
  (B_crit=1+A^2=2), flagged for inspection, not absorbed.

### THE PIVOT -- caricatures retired

Greg asked where the gravity/weak laws come from. Answer: our own toy
caricatures (Session 8), not OD and not literature. Gravity =
G*E_total*dx/(dx^2+eps_soft) (softened, energy-mediated -- not Newton,
not GR). "Weak" = K*x_j*exp(-M*dx^2) -- a GAUSSIAN, not a Yukawa
exp(-Mr)/r. EM = plain bilinear; strong = linear+cubic confinement. All
hand-built. Every four-force result (INFO-027/031/032/033, the flip
ordering, EM-poly-vs-strong-exp) is therefore a statement about EQUATIONS
WE WROTE, not about real forces -- the frame had been partly grading
itself. Greg's call: "no point fine-tuning fake data ... if replacing
would just be more work on fake data, skip that too." -> items 2, 3, and
caricature replacement SKIPPED; pivot to real data, four-force-related.

### Item 5 set 1 -- LIGO GW150914 (INFO-036, REAL DATA)

`s10_ligo_extract.py`, `data/ligo/`. GWOSC public 32s/4096Hz H1+L1
strain. Bandpass 35-350Hz + ASD whiten, 2s edge crop (the canary's
spurious peak-MI at a segment edge flagged the filter-edge artifact),
125ms windows / 31ms stride, L1 aligned (-7ms, sign-flipped). Reused
`kbk_pipeline.compute_operator_matrix` verbatim.

- METHOD WORKS: inter-detector windowed MI peaks EXACTLY at the merger
  (t=16.41s vs 16.4s; MI 0.530 vs noise baseline 0.247, 2.1x).
- PRE-REGISTERED PREDICTION OVERTURNED: predicted the channel-substrate
  appears IN the event; instead noise-only windows sit ON the
  (-1,-1,+2)/sqrt6 attractor (|cos|=0.984) and event windows LEAVE it
  (|cos|=0.242). The attractor IS the equal-marginal-entropy identity
  H_a~=H_b: whitened noise has equal per-channel entropy -> trivially on
  it; the chirp changes one detector's entropy -> off it, while MI (a
  separate operator) spikes. Confirms the Session 5 deflationary reading
  on REAL data. Caveat: one event; noise-on-attractor partly a whitening
  consequence (the point).

### Item 5 set 2 -- PDG four-force unification (INFO-037, REAL DATA, special build)

`s10_pdg_unification.py`. The unification question in native form
(couplings vs energy, not a time series -> purpose-built). Frame map:
substrate = shared linear-in-ln(Q) running form; expression = per-force
slope b_i. Anchored on solid PDG M_Z couplings.

- Shared running FORM confirmed in data: SM one-loop alpha_s(Q) matches
  measured determinations <1sigma from 31 GeV to 1 TeV (low-Q pulls =
  known one-loop limitation).
- NO single SM unification: crossings at 1.0e13 / 2.4e14 / 9.7e16 GeV --
  a triangle spanning ~9400x.
- MSSM near-point ~2.1e16 GeV (spread 1.1x) but unobserved SUSY ->
  conjecture, not support.
- Gravity: alpha_G(E)=(E/M_Pl)^2 power-law, outside the log-running
  family.
- Verdict: shared running substrate + per-force expression among the 3
  gauge forces; no single unification without conjectural new physics;
  gravity outside the form.

## Joint reading (three levels, deflationary present)

- Data: (4) flip is estimator-robust; (1) substrate invariance is mixed
  (2/4 clean, 1 regime-bounded, 1 failure); (5a) our stack detects a real
  GW merger via MI and the "attractor" is equal-entropy on real noise;
  (5b) three gauge couplings share a running form, do not meet at a
  point, gravity is outside the form.
- Interpretation: the "shared substrate + per-force expression" pattern
  appears in BOTH real-data sets in a deflated form -- in LIGO the
  "substrate" is an equal-statistics fact, in PDG the substrate is a real
  shared running form but does NOT force unification. The strong version
  of the frame (a deep common substrate that forces the forces together)
  is not supported by either real-data set.
- Frame: still a frame. What changed: it now has its first real-data
  contact, and the contact is deflationary on the substrate side and
  confirmatory on the method side (the extraction works on real data).

## Queued for Session 11

1. LIGO: more events (GW151226, GW170817 BNS, O3/O4) and a proper
   noise-vs-signal null distribution -- is the MI-at-merger peak
   significant against many noise segments? Is the equal-entropy-attractor
   reading event-independent?
2. LIGO without whitening / with different whitening -- separate the
   equal-entropy-attractor result from the preprocessing (the open caveat).
3. PDG: two-loop running + flavor thresholds to tighten the alpha_s(Q)
   validation chi2 (one-loop gave 2.5, low-Q dominated); pull the REAL
   PDG alpha_s(Q) table programmatically rather than the curated overlay.
4. The chemistry substrate failure (INFO-035): inspect the B=2 Brusselator
   Hopf-threshold region directly -- is the rank swing a bifurcation
   crossing? (find the specific reason; do not widen the tent.)
5. Open: storm/waves and ECG real-data sets remain (NOT four-force; Greg
   scoped this session to four-force-related data only).

## Files produced this session (repo root)

Scripts:
- s10_third_estimator.py        (item 4: KSG triangulation)
- s10_substrate_invariance.py   (item 1: substrate vs expression invariance)
- s10_ligo_extract.py           (item 5 set 1: LIGO real data)
- s10_pdg_unification.py        (item 5 set 2: PDG unification, special build)
- mapping_campaign.py           (pulled verbatim from Session 8 branch, for item 1)

Result JSONs:
- s10_third_estimator_canary.json, s10_third_estimator_results.json
- s10_substrate_invariance_canary.json, s10_substrate_invariance_results.json
- s10_ligo_extract_canary.json, s10_ligo_extract_results.json
- s10_pdg_unification_results.json

Data:
- data/ligo/H-H1_GW150914_32s.hdf5, data/ligo/L-L1_GW150914_32s.hdf5
  (GWOSC public strain)

Env: .claude/hooks/session-start.sh hardened (guarded the
CLAUDE_PROJECT_DIR/CLAUDE_ENV_FILE line). h5py 3.16 added. PySR 1.5.10 +
Julia re-bootstrapped via the hook.

## Branch state

- Local + remote: claude/gravity-substrate-config-51cfp. Pushed.
- main untouched. No PR.

End of v10 results handoff.
