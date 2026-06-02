# SESSION HANDOFF — 2026-06-02 (S17) — CLAUDE.md drift fix + construction-vs-nature probe

Branch `claude/claude-md-problem-nM8Hm` (main untouched; no PR). All Operating Rules
in force. Two things this session: JOB 1 (the CLAUDE.md drift) then backlog #1 (the
construction-vs-nature probe Greg flagged "clean, decisive").

## JOB 1 — CLAUDE.md drift FIXED

The repo `CLAUDE.md` body was canonical only through S11 with a header stale at "S7";
the true master was S14 (with S15/S16 in handoffs). This branch was also cut from
`main` at the S12 kickoff, so it was missing all S13-S16 work.

- Fast-forwarded this branch to absorb the S16 work branch
  (`claude/claude-md-strategy-pwfZr`, 10 commits) — pulled in `CLAUDE_master_through_S14.md`,
  the S15 note, the S14-S15 combined handoff, the S16 handoff, the backlog, the S17
  kickoff, plus all S13-S16 scripts + result JSONs.
- Verified the S14 master is a STRICT SUPERSET of the S11 body before overwriting
  (diffed section headers + ledger entries: nothing dropped; all 44 INFO entries
  preserved). The S14 master already carried S13/S14/S15 notes.
- Rebuilt `CLAUDE.md` from the S14 master, folded in a Session 16 note + (this
  session) a Session 17 note, re-attached the START-HERE workflow block, fixed the
  header, and updated the Session Handoff Pointer.
- Kept `CLAUDE_master_through_S14.md` + the S15/S16 handoffs in-repo as historical
  artifacts. The "keep the header current" Operating Rule (added S14) is what stops
  the drift recurring.

## PROBE (backlog #1) — construction-vs-nature of the equal-entropy substrate

THE question (Greg, S16): everything we measure sits ON the equal-entropy region; is
equal-marginal-entropy a FORCED property of physical 2-channel observables (nature),
or our construction choice (bookkeeping)?

THE LEVER (`probe_construction_vs_nature.py`): scale ONE channel by a constant s,
`b -> s*b`. This is (a) MI-INVARIANT — `mi_hist_2d` uses adaptive (data-range) bin
edges, so scaling rescales the edges proportionally and the 2D histogram is
IDENTICAL, MI(a, s*b) = MI(a, b) exactly; the coupling/physics is untouched — and
(b) a pure UNITS choice — differential entropy is not scale-invariant,
`H(s*b) = H(b) + ln|s|`. So scaling injects pure construction asymmetry into the
marginal entropies while leaving the genuine coupling invariant. Reuses
`kbk_pipeline` (vasicek H + hist MI + extract_v1 + project_234) verbatim.

Design: 4 heterogeneous sim generators (OU correlated, AR(1)-Laplace, coupled
logistic chaos, sine+noise) x 5 seeds [11,22,33,44,55] x scale sweep
[0.1,0.2,0.5,1,2,5,10]; PLUS 3 real LIGO segments scaled post-whitening
(GW150914 H1/L1 aligned per s10; two independent noise-only segments V1/V2). Metric:
mean|H_a-H_b| (asymmetry, the INFO-038 knob), cos(v_null_234,(1,1,2)/sqrt6) [=
project_234, the s10/s11 substrate metric], cos to (-1,-1,2)/sqrt6 (the -(H_a-H_b)^2
identity), |MI coef in null|, mean MI.

### Results (data level)

1. **MI is EXACTLY scale-invariant.** MI_cv (coeff of var of mean MI across the
   whole scale sweep) = 1e-16 to 9e-16 (machine precision) for ALL 7 systems. The
   genuine coupling is representation-independent — a units choice changes no physics.

2. **Marginal-entropy asymmetry is a pure units knob.** asym_range (max-min of
   mean|H_a-H_b| over the scale sweep) = 2.1 to 3.8 for 6 of 7 systems, tracking
   `|H_a - H_b0 - ln s|` exactly: symmetric-construction sims (baseline H_a~=H_b)
   give a monotone `~ln s` rise; LIGO segments give the predicted V-shape with the
   minimum where `ln s` matches the baseline gap (V1 asym 0.78->0.32->1.52, min near
   s~3; GW150914 min near s~2). You can place any system on or off the
   equal-marginal-entropy point by choosing units for one channel.

3. **The substrate cos metric is NOT scale-invariant** (cos112_std over the sweep
   0.07-0.27; cosident_std 0.14-0.28) and is fragile/segment-specific. At s=1 the
   independent noise segments V1/V2 do NOT reproduce s10's GW150914-noise cos 0.984
   (they sit at 0.12 / 0.38). NOTE on the s10 "0.984": that null is DOMINATED by the
   LINEAR H_a,H_b terms (null6d [+0.569,+0.774,...]); project_234 keeps only the
   small quadratic residual, whose projection onto (1,1,2) happens to be ~0.98 for
   that particular near-constant-entropy (whitened) segment. It is an estimator-noise
   residual, not a robust geometric fact. cos generally peaks near the
   asymmetry-minimum and degrades as you rescale away.

### Reading (Result Discipline: data / interpretation / frame, alternatives mapped)

- DATA finding (INFO-051, located): MI is exactly scale-invariant; marginal-entropy
  asymmetry is set by a per-channel units choice (`H(sX)=H(X)+ln s`); the substrate
  cos metric is representation-dependent and fragile.
- INTERPRETATION (deflationary, supported): "everything falls on the equal-entropy
  substrate" is BOOKKEEPING. Systems sit there because we build/normalize channels
  to comparable scales (sims: equal noise amplitude; LIGO: whitening to unit
  variance). A physics-preserving units change (MI fixed exactly) dissolves the
  clustering. The genuinely physical, scale-invariant content in the basis is MI.
- ALTERNATIVES mapped + ruled against: (1) a hidden scale-invariant substrate
  signature — the only scale-invariant quantity in the basis is MI, which is not the
  entropy-quadratic "substrate", so no; (2) nature enforces equal scales for real
  channel pairs — for LIGO that is instrument-engineering + whitening (a construction
  choice), not a law. Both reduce to construction.
- This RESOLVES the S16 backlog-#1 open question on the construction side and is
  internally consistent with INFO-036 (attractor = equal-marginal-entropy identity)
  and INFO-038 (asymmetry moves objects off): asymmetry moves objects off BECAUSE it
  is a units knob.

### Distinct regime kept as data (don't predetermine good/bad -- Greg, S17)

`logistic_chaos` is the lone system with small asym_range (0.23). Cause: coupled
logistic maps (r=3.9, eps=0.15) ANTI-synchronize (corr -0.995) -- the two channels
become near-perfect mirror images. In that LOCKED regime the units-knob barely moves
the asymmetry (ln s shift +0.076 vs expected +2.303 at s=10), unlike every
independent system. Two things are jointly true and both are data: (1) the windowed
Vasicek estimator hits its limit on a near-deterministic signal (entropy ~ -620,
dominated by near-zero-diff windows), and (2) strongly-coupled / locked channels
genuinely respond differently to rescaling than independent ones. This is NOT garbage
to exclude -- it is a real "what happens at near-perfect coupling" data point that may
deserve its own probe (e.g. a different entropy estimator, or treating perfect
coupling as a regime of interest). Logged, not predetermined as bad.

### Caveats / open
- The quadratic-null cos metric (project_234) is estimator-noise-sensitive
  (INFO-024) and segment-specific; do not over-read individual cos values.
- Speaking posture (after): the data supports the bookkeeping reading decisively at
  the marginal-entropy level (MI exactly invariant; asymmetry = units knob). The
  fragile substrate-cos level is consistent with it but not load-bearing.

## PROBE 2 (backlog #3) — gravity chirp law recovered from real strain (INFO-052)

GRAVITY-FORWARD by design (the pivot from S17 PROBE 1's deflationary edge result):
recover a governing LAW on the scale-invariant content, the gravity analogue of the
S16 WF (Z propagator) + SF (QCD running) wins. Newtonian inspiral predicts
u(t)=f(t)^(-8/3) is LINEAR in t with slope k = -(256/5) pi^(8/3) (G M_c/c^3)^(5/3),
so a clean linear u-vs-t + slope -> chirp mass = the law recovered. Method
(`probe_gravity_chirp_ridge.py`): own FFT-based Morlet CWT (scipy.signal.cwt removed
in 1.17) -> scalogram; ridge tracked BACKWARD from the merger column with continuity
constraint (search band +/-25% log-index, gap tolerance 12) so it follows the
descending inspiral arc instead of jumping to noise; fit u vs t on the rising
inspiral. Replaces the S16 per-sample whitened-Hilbert that gave R^2 0.001.

Results (data level):
- **GW150914 (BBH): CLEAN recovery, both detectors agree.** H1 R^2 0.995, M_c 38.40
  Msun, ridge f 34->155 Hz, t_c 16.424 s; L1 R^2 0.987, M_c 38.18 Msun, f 49->174 Hz,
  t_c 16.419 s (true merger 16.40 s). Catalog detector-frame chirp mass ~31 Msun, so
  recovered M_c is ~24% high -- the expected bias of a Newtonian-only f^(-8/3) fit on
  the late/strong-field arc (PN + whitening shape the ridge). The f^(-8/3)-linear-in-t
  LAW is confirmed (R^2 0.99) and the two independent detectors agree to <1% on M_c.
  DECISIVELY beats S16 (Hilbert R^2 0.001 -> ridge R^2 0.99), exactly as backlog #3
  predicted.
- **GW170817 (BNS): honest partial.** Merger from file metadata (GPSstart 1187008867
  + known merger GPS 1187008882.43 = 15.43 s; envelope detection had failed). The
  f^(-8/3) LAW FORM fits this second, physically-distinct system (R^2 0.879,
  t_c 15.430 s), but absolute M_c is biased (15.67 vs catalog ~1.20 Msun) because the
  visible H1 arc is a narrow low-frequency band (51->76 Hz, short lever arm); the
  faint high-frequency sweep to merger is data/method-limited in this 32 s H1 segment
  (BNS SNR is mostly in L1, which carries the known glitch). Kept as data: the law
  FORM generalizes to neutron stars; absolute-mass recovery needs the full sweep
  (longer data / louder detector / matched filter), not manufactured.

Reading (Result Discipline): DATA -- the Newtonian inspiral law u=f^(-8/3) linear in
t is recovered from raw GW150914 strain on both detectors (R^2 0.99), M_c ~38 vs
catalog ~31; the law FORM also fits GW170817 (R^2 0.88) with a lever-arm-limited mass.
INTERPRETATION -- this is a positive gravity governing-law recovery on the scale-
invariant content (not the bookkeeping substrate), the gravity analogue of WF/SF
Piece-1. CAVEAT -- Newtonian-only fit => ~20-25% M_c bias on the strong-field arc is
expected; GW170817 absolute mass is not constrained by the H1 arc here. No new
Operating Rule. Speaking posture (after): GW150914 is a clean positive; GW170817 is
a partial that confirms the form, not the mass -- recorded as such.

ADJUSTMENT NOTE (within-probe, per the new rule): backlog #3's S16 diagnosis
("needs Q-transform/matched-filter, not naive Hilbert") is CONFIRMED -- the CWT ridge
was the fix. Merger-seeding from metadata (not envelope-peak / max-power) was needed
for the faint BNS.

## PROBE 3 ("MI next" -- Greg) — inter-detector MI is physics-bearing (INFO-053)

Characterize the ONE genuinely physical footing left after PROBE 1 (the equal-entropy
substrate is bookkeeping; MI is the only scale-invariant operator). THE question: is
the MI-at-merger (INFO-036/045) genuine GRAVITATIONAL information carrying structure
BEYOND loudness, or just "two detectors saw the same loud thing"? Four tests on
GW150914 H1/L1 (`probe_mi_merger_axis.py`), same hist-MI estimator. (MI is invariant
under L1's sign inversion, so only the inter-detector LAG matters.)

- **T1 MI(t) through merger**: peak MI 0.455 at t 16.366 s (merger 16.42 s) vs
  baseline 0.205 (2.2x). Reproduces INFO-036.
- **T2 LAG SCAN (decisive geometry test)**: merger-window MI peaks at lag
  **+7.0 ms == the physical H1-L1 light-travel delay** for GW150914; the off-merger
  noise window is FLAT (max 0.256). A loudness coincidence is lag-independent -- this
  tracks the gravitational GEOMETRY.
- **T3 TIME-SLIDE NULL**: merger MI 0.534 vs null (10 large non-physical lags) mean
  0.225 sd 0.020 => **z = 15.5**. Highly significant; not coincidence.
- **T4 PHASE-SCRAMBLE (beyond loudness)**: phase-randomize L1 (keep amplitude
  spectrum, destroy waveform phase) -> MI 0.240 +/- 0.029 (50 realizations), ~= the
  null (0.225) and baseline (0.205). So the ENTIRE MI excess above baseline is carried
  by waveform PHASE/TIME structure, not amplitude/loudness.

Reading (Result Discipline): DATA -- inter-detector MI at the merger peaks at the
physical 7 ms lag, is 15.5 sigma above the time-slide null, and its excess is entirely
destroyed by phase-scrambling. INTERPRETATION -- the MI operator is PHYSICS-BEARING:
genuine common gravitational information carrying the waveform's phase/time structure,
NOT loudness coincidence or bookkeeping. This validates MI as the real physical
footing (consistent with INFO-051: MI is the only scale-invariant quantity). NOT new
physics -- LIGO routinely uses inter-detector consistency + the ~7 ms delay; our
contribution is that OUR information-theoretic MI operator captures it AND that the MI
signal is time-structure-bearing. FRAME CONTACT (ungraded, frames don't grade
themselves): MI -- the surviving physical quantity -- carries the gravitational TIME
structure; this is the empirical contact point for the "gravity/flow couples to time"
hunch (H-A/H-B/H-C), but it stays a frame until a probe separates "MI tracks the
chirp's time structure (already in GR)" from "MI carries time-coupling beyond GR."
No new Operating Rule.

## Credibility / outreach asset (Greg flagged HUGE, S17)

Greg: "this is huge and something we would want to tout to the right company ...
proves our OD machinery can pull a real gravity governing law out of raw data with
both detectors agreeing. That's credibility." Created a standing **Capability
Demonstrations** section in CLAUDE.md collecting the three raw-data governing-law
recoveries (GRAVITY chirp INFO-052 + MI INFO-053; WEAK Z propagator S16; STRONG QCD
running S16). Sellable claim: a domain-agnostic engine that recovers ESTABLISHED laws
from RAW PUBLIC DATA with no physics assumptions, across physics' three hardest force
domains -> credibility that it can find governing laws where they are UNKNOWN.
Honest: these RECOVER known laws (checkable against ground truth), not new physics --
which is the point. ACTION pending (Greg): package a one-page capability brief.

## Gravity thread — footholds scaffold (start of framework accounting, Greg S17)

Honest accounting of where the gravity work stands. This is NOT an equation yet; it is
the ledger of what may stack toward one. Did we find a NEW law? NO -- INFO-052
recovered a KNOWN law (Newtonian chirp); its value is method validation.

CONFIRMED FOOTHOLDS (physical, build on these):
- F1 (INFO-053): inter-detector MI at merger = genuine gravitational common info --
  physical 7 ms lag, z 15.5, carries waveform phase/time structure.
- F2 (INFO-052): inspiral chirp governing law f^(-8/3) ~ (t_c - t) recovered from raw
  strain, both detectors agree (KNOWN law; validates the OD method on real gravity).
- F3 (INFO-051): MI is the only scale-invariant (representation-independent) quantity
  in the operator basis.

CLEARED (bookkeeping/artifact -- do NOT build on):
- C1 (INFO-051): the equal-entropy "substrate" everything clustered on = units/
  construction bookkeeping.
- C2 (S16): post-merger "return to substrate" = return to detector NOISE, not ringdown.
- C3 (Greg S16): MI-in-null coupling discriminator DISPROVED.

OPEN (toward a possible gravity equation):
- O1: does MI carry MORE than the GR chirp -- any time/flow structure beyond what the
  inspiral already encodes? (the new-physics question; the only road to a NEW law)
- O2: the time-coupling thesis (gravity CONTAINS a time term) -- design on MI/dynamics
  (the scale-invariant content), NOT the entropy substrate.
- O3: GW170817 absolute chirp mass (lever-arm-limited in H1; needs full sweep / L1 /
  matched filter).

## PROBE 4 (O1) — does MI carry structure beyond the recovered chirp? (INFO-054)

The only road to a NEW gravity law: is the inter-detector MI more than the common GR
chirp? HONEST SCOPE up front: a rigorous beyond-GR test needs IMR matched-filter
templates (pycbc/lalsuite), not in our self-contained stack. This tests "beyond the
recovered INSPIRAL law (INFO-052)": fit + subtract the Newtonian chirp model from each
GW150914 detector over the late inspiral (where the MI peak lives, T1 ~16.366 s < t_c)
and recompute inter-detector MI of the residuals vs a time-slide null
(`probe_mi_beyond_chirp.py`).

Result (data): the Newtonian model removed only 51% (H1) / 38% (L1) of window variance
(crude; covers only t<t_c). Inter-detector MI peak at the physical lag: FULL 0.624 at
+7.5 ms (z=10.0 vs null 0.290+/-0.034); RESIDUAL 0.325 at +7.5 ms (z=8.7 vs null
0.199+/-0.014). So subtracting the recovered inspiral cuts MI ~1.9x but a SIGNIFICANT
residual remains at the correct physical lag.

Reading (Result Discipline): INCONCLUSIVE for "beyond GR" -- NOT a new-law signal. The
residual is most plausibly un-modeled GR (merger + ringdown + higher PN, all omitted by
the crude Newtonian inspiral model that removed <55% of variance and stops at t_c), not
new physics. DECISIVE methodological output: separating "beyond GR" from "un-modeled
GR" REQUIRES IMR matched-filter templates to remove the FULL GR waveform; with
self-contained tools O1 cannot go further. Establishes the test method (residual
inter-detector MI vs time-slide null) and pins the exact requirement. No new law. (Bug
caught + fixed mid-probe: the time-slide null first applied large lags to the short
extracted window -> empty -> z~1e11; fixed to pull L1 null windows from far times in the
full array. Don't-predetermine: kept and diagnosed.) Speaking posture (after): I
predicted reduce-but-not-zero (merger/ringdown remain); that is what happened; the
honest verdict is "need templates," not "found new structure."

Scaffold update: O1 -> ATTEMPTED, inconclusive, requires IMR templates (new backlog
item). The new-law question on the MI axis is now a DEFINED next step, not open-ended.

## PROBE 5 (O1-real) — beyond-GR via IMR templates: clean NEGATIVE (INFO-055)

The rigorous version O1/INFO-054 said was needed. pycbc 2.11.0 installed
(`--ignore-installed cryptography` to dodge a debian uninstall conflict; add to
requirements if persisting). Generate IMRPhenomD (m1=36, m2=29 Msun), whiten it with
each detector's ASD, fine time-shift search + 2-quadrature lstsq fit to the whitened
GW150914 data over the merger window, subtract the maximum-likelihood GR waveform,
recompute inter-detector MI of residuals vs time-slide null (`probe_mi_beyond_GR_imr.py`).

Result (data): IMR template removed 67.7% (H1) / 53.3% (L1) of merger-window variance.
MI at the physical +7 ms lag: FULL 0.536 (z=18.0 vs null 0.227+/-0.017) -> RESIDUAL
0.226 == null 0.209+/-0.019. The residual PEAK (0.263) drifts off to -18.5 ms (search-
edge) at z=2.9 -- consistent with the null, and NOT at the physical lag.

Reading (Result Discipline): LOCATED, NEGATIVE. The inter-detector MI merger signal is
FULLY accounted for by the GR waveform -- removing it collapses the physical-lag MI to
chance. NO detectable common structure beyond GR on the MI axis for GW150914 at this
sensitivity. This is a real, valuable NEGATIVE (falsification-first): we built the
rigorous new-physics test and the answer is "it's GR." CAVEATS: single event; template
masses fixed (not per-detector refit), so subtraction removed ~55-68% of variance --
but the decisive metric is that the PHYSICAL-LAG MI drops to null and the residual peak
is off-lag/insignificant (the leftover variance is uncorrelated detector noise, which
carries no inter-detector MI -- exactly why residual MI -> null). Speaking posture
(after): I expected the full template to drop residual MI toward null; it did; verdict =
no structure BEYOND GR on this axis.

CORRECTION (Rule D -- Greg caught S17): an earlier line here said this "also closes O2's
premise (the time-coupling hunch)." That was incomplete-not-wrong and is RETRACTED. O2
(gravity-couples-to-time) was NEVER a beyond-GR question: gravity/time coupling is
INTRINSIC to GR (the chirp IS time-evolution; time dilation; proper time). So finding
the signal is "all GR" is CONSISTENT with the flow/time hunch, not against it -- indeed
the gravity signal's entire informational content BEING its time structure is a (frame-
level) data point consistent with the hunch. INFO-055 closes ONLY "beyond-GR structure
on the MI axis." O2 is reframed as a FRAMEWORK/REPRESENTATION question (below) and stays
OPEN.

Scaffold update: O1 -> RESOLVED NEGATIVE (beyond-GR only). O2 -> OPEN, REFRAMED: not a
new-physics question but "can OD express gravity's GR-real time/flow coupling as a flow
dipole term, and does that term distinguish gravity from the gauge forces (H-C: only
gravity touches time)?" -- testable by comparing a flow/time operator's response on
gravity (chirp, has time-evolution) vs the gauge-force objects (no chirp-like
time-evolution). The honest status: NO new gravity LAW found; the surviving physical
signal is fully standard GR (good -- that is the falsifiable negative). A BEYOND-GR
new-law search would need a different observable/axis, or many
events, not this one.

## Files this session
- `CLAUDE.md` (rebuilt canonical through S17), `SESSION_HANDOFF_2026-06-02_S17.md`,
  `BACKLOG_tests_and_probes.md` (#1 marked DONE).
- `probe_construction_vs_nature.py` + `probe_construction_vs_nature_results.json`
  (+ `_canary.json`).
- `probe_gravity_chirp_ridge.py` + `probe_gravity_chirp_ridge_results.json`
  (+ `_canary.json`).
- `probe_mi_merger_axis.py` + `probe_mi_merger_axis_results.json` (+ `_canary.json`).
- `probe_mi_beyond_chirp.py` + `probe_mi_beyond_chirp_results.json`.
- `probe_mi_beyond_GR_imr.py` + `probe_mi_beyond_GR_imr_results.json` (needs pycbc).
- h5py confirmed in `requirements.txt` (installed at runtime; the fresh container's
  SessionStart hook had not installed it this run).

## Next (backlog, Greg picks order; clear backlog before new probes)
Remaining strong candidates: #2 construction control gravity-LIGO vs EM-HBT flow
(needs Zenodo 5113016 .rar + S13 HBT loader); #3 proper chirp recovery
(Q-transform/matched-filter); #7 SF femtoscopy-R; plus housekeeping #11 (OD /
MASTER_DISCOVERIES) and the now-mostly-done #12 (master merge).
