# BACKLOG — tests & probes (as of 2026-06-02, S16)

## STANDING RULE (Greg, S16)
**Complete the backlogged probes below BEFORE starting any new probe we think of.**
New ideas get appended here first; we clear the backlog before opening new lines.
A probe that didn't give a 100% hit still carries a piece — keep its output, don't discard.

**Run-order is NOT fixed (Greg, S16): the next session decides which PROBE to run
first.** The numbering below is just an index, and the HIGH/MED/LOW tags are rough
suggestions — not a prescribed sequence.

**New-idea placement (Greg, S17): new probe ideas go to the BOTTOM of the list,
UNLESS they are (a) critical, (b) a reframe of an existing probe, or (c) an
adjustment of a previous probe — those can slot in place.**

**Don't predetermine good/bad (Greg, S17): crazy data might still be valuable data.
Do not label an output "pathology / bad / exclude" on first look. Diagnose it and
keep it as a data point; an odd regime may deserve its own probe.**

**Stop after each probe (Greg, S17): pause after each probe so Greg can steer the
next one. Do not chain probes without a check-in.**

**EXCEPTION — JOB 1 (Greg, S16): fix the CLAUDE.md drift FIRST, before any probe.**
The repo CLAUDE.md body is canonical only through S11; the true master is S14
(preserved as `CLAUDE_master_through_S14.md`); S15/S16 are in handoffs. Fold S12-S16
into one canonical, current CLAUDE.md (safe procedure in the CLAUDE.md START-HERE
block). This was backlog #12 — promoted to JOB 1. Do it, THEN pick a probe.

---

## NEXT PROBE (Greg, S18) -- FLOW DIPOLE EQUATION on the flows (TOP of queue)
Reframe/continuation of the O2/flow thread (slots in, not bottom). S18 confirmed FLOW is a
substrate with TIME as gravity's expression, generalizing to all 4 forces + chemistry
(INFO-062/063/064). Greg's next ask: figure out the flow DIPOLE EQUATION on these flows.
Each flow is currently a single characteristic vs an axis; a flow DIPOLE needs TWO coupled
channels in the info-dipole paper's form:
  dMI/dt ~ sum_i c_self,i*H_i^2 + sum_{i<j} c_cross,ij*H_i*H_j + linear   (opposition signature)
PLAN: pick 2 channels per flow -- chemistry Brusselator = native 2 species x,y (START HERE,
cleanest); gravity = H1/L1 detectors (inter-detector MI already in INFO-053); strong/EM/weak
= two observables. Extract the flow-dipole equation per system, then test substrate/expression
AT THE EQUATION LEVEL: shared dipole FORM (substrate) + per-domain COEFFICIENTS (expressions).
Connect to the info-dipole paper (davisai.ai/dipole) + the Markets flow dipole. Guard against
the INFO-051 bookkeeping artifact (MI is the only scale-invariant operator). Greg's epistemic
rule: a non-conforming system is one data point / maybe wrong observable, NOT a falsification.

## GRAVITY thread (primary focus)

1. **Substrate construction-vs-nature test** -- DONE (S17, INFO-051). Answer:
   BOOKKEEPING. Scaling one channel (b->s*b) is exactly MI-invariant (MI_cv ~1e-16
   across 7 systems) yet sets the marginal-entropy asymmetry by a pure units choice
   (H(sX)=H(X)+ln s; asym_range 2-4 on 6/7 systems). So equal-marginal-entropy is a
   construction choice, not nature; the only scale-invariant (genuinely physical)
   quantity in the basis is MI. The substrate cos metric is itself representation-
   dependent + fragile (independent noise segments don't reproduce s10's 0.984).
   See `SESSION_HANDOFF_2026-06-02_S17.md`, `probe_construction_vs_nature.py`.
   CONSEQUENCE for the gravity thread: the equal-entropy "substrate" the gravity
   story leaned on is a coordinate artifact -- pivot gravity-forward probes onto the
   scale-invariant content (MI-merger axis + governing laws), not the substrate.

2. **Construction control: gravity-LIGO vs EM-HBT flow axis** (HIGH; decisive for
   the "only gravity's flow is unbound/off-substrate" hunch).
   Measure a non-gravity REAL detector-pair (EM/HBT, INFO-048) flow axis the same
   way (P3 method). On-substrate (~0.9) => gravity special; off (~0.6-0.8) =>
   real-detector/noise effect. DATA: Zenodo 5113016 (demonDataPublic.rar, 20.8MB)
   + unrar + port the S13 HBT loader (s13_force_hbt.py, on claude/file-upload-memory-LKCoB).

3. **Proper chirp recovery (gravity Piece-1, redo)** -- DONE (S17, INFO-052) for
   GW150914; PARTIAL for GW170817. A Morlet-CWT ridge tracked backward from the merger
   recovers u=f^(-8/3) linear in t: GW150914 both detectors R^2 0.99, M_c ~38 vs
   catalog ~31 (Newtonian-late-inspiral ~24% high), decisively beating the S16 Hilbert
   (R^2 0.001) -- so the S16 limiter was the METHOD, confirmed. GW170817 confirms the
   LAW FORM (R^2 0.88) but absolute M_c is biased (15.7 vs 1.20): the visible H1 arc is
   a narrow 51->76 Hz band (short lever arm). See `probe_gravity_chirp_ridge.py`.
   REMAINING (new sub-item, bottom-of-thread): GW170817 absolute chirp mass needs the
   full sweep to merger -- louder detector (L1, has the glitch), longer data, or
   matched-filter template tracking. Append as a follow-up, not a blocker.
   -- DONE (S18, O3, INFO-061): GW170817 BNS chirp mass NAILED. Fetched GWOSC 4096s/4096Hz
   H1+L1, TaylorF2 chirp-mass template bank + matched filter, L1 glitch gated with pycbc
   .gate(): detector-frame M_c = 1.200 Msun (network, 0.19% from catalog 1.1977), H1 1.200 /
   L1 1.195 independently agree, net SNR 14.8 (L1 null 9.7 sigma). Removes the S17 ridge
   absolute-mass bias. See `probe_o3_gw170817_chirpmass.py`, `SESSION_HANDOFF_2026-06-02_S18.md`.

4. **Louder-event flow reproduction** (BACKBURNER per Greg, S16).
   Fetch GW170814/GW190521 strain; re-run per-event flow axis; test whether
   gravity's off-substrate flow displacement is consistent (magnitude steady /
   direction wandering = "shapeless but real" vs both wandering = noise).

5. **Scatter structure (H-G)** (LOW; depends on #4).
   Does per-event gravity coordinate spread (cos->substrate 0.78/0.64/0.99) track
   SNR / total mass / detector-noise state / alignment lag? Needs several clean events.

6. **Gravity rate-dynamics (PROBE 1 gravity half)** (LOW; gated).
   dMI/dt lagged structure on real gravity. Merger window has only ~8 independent
   samples (INFO-024 wall) — needs the long inspiral at finer/lower band, not louder
   mergers. Sim half already done (chem>phys>bio memory; geology null).

6b. **MI beyond GR -- IMR-template residual test** -- DONE (S17, INFO-055), clean
   NEGATIVE. pycbc 2.11.0 IMRPhenomD matched-filter subtraction on GW150914: removing
   the full GR waveform collapses the physical-lag inter-detector MI to the null (0.54
   z=18 -> 0.23 == null; residual peak off-lag, z=2.9). The MI merger signal is FULLY
   GR -- NO beyond-GR structure on this axis. (Does NOT close O2: gravity/time coupling
   is intrinsic to GR, so "all GR" is consistent with the time-coupling hunch, not
   against it -- Rule D correction, Greg S17. O2 stays open, reframed as a framework
   question -- see below.)
   See `probe_mi_beyond_GR_imr.py`. INSTALL NOTE: pycbc needs
   `pip install pycbc --ignore-installed cryptography` (debian cryptography uninstall
   conflict); NOT added to requirements.txt to avoid breaking the SessionStart hook.
   REMAINING (bottom, optional): repeat across many events / per-detector mass refit to
   strengthen the negative; but for GW150914 the new-law road on the MI axis is closed.

6c. **O2 (reframed) -- gravity's time/flow coupling as a framework term** -- PARTIALLY
   DONE (S17, INFO-056). First pass via a flow operator |Spearman(characteristic, axis)|
   on the 3 recovered governing relations: gravity chirp |rho|=1.0 (FLOW, axis time),
   strong alpha_s(lnQ) |rho|=1.0 (FLOW, axis scale), weak Breit-Wigner |rho|=0.55
   (resonance, no flow). VERDICT: flow is NOT unique to gravity (strong RUNS) -- the
   simple "only gravity flows" hunch is REFUTED. Surviving narrow form ("gravity flows
   in TIME specifically") is construction-confounded and NOT decided by this data. See
   `probe_flow_operator.py`. REMAINING (the only clean way forward on H-C, bottom):
   gravity-time-specialness is fundamentally a TIME-DILATION/CLOCK claim -- test it with
   clock/time-dilation data (e.g. GPS/optical-clock/Pound-Rebka-style or pulsar timing),
   NOT force-coupling objects. Add EM/HBT to the flow comparison only as a completeness
   check (won't resolve the construction confound).
   -- REMAINING NOW DONE (S18): clock/time-dilation data run, gravity-couples-to-time
   recovered POSITIVELY three ways. GPS term-retaining route (INFO-059) recovered the
   time-dilation coefficient -2/c^2 from raw RINEX observations on eccentric Galileo sats
   at k/truth 1.02-1.04, z=380 (precise-product route was a definitional NULL, INFO-058);
   pulsar route (INFO-060, raw Arecibo TOAs PSR B1913+16) recovered orbital decay dP_b/dt
   (ratio 1.005 to GR) + Einstein-delay gamma (0.014%). See `SESSION_HANDOFF_2026-06-02_S18.md`,
   `probe_6c_gps_positive.py`, `probe_pulsar_time.py`, `probe_gravity_time_dilation.py`.
   STILL OPEN (gravity-SPECIAL-vs-gauge-forces in touching time): needs the multi-force
   comparison (EM/HBT, backlog #2) -- a completeness check, won't resolve the construction
   confound. The narrow "gravity dilates clocks" claim is now POSITIVELY supported by data.

## 4-FORCE thread

7. **SF femtoscopy-R Piece-1** (MED; deferred — heavy fetch).
   Recover the femtoscopy source radius R from C(q) on ATLAS heavy-ion data
   (data/strong/, gitignored; opendata.cern.ch streaming bypass). The alpha_s(Q)
   running Piece-1 is DONE (asymptotic freedom recovered from data).

8. **Forces-within-construction flow comparison** (MED; depends on #2).
   Once real gauge-force flow objects exist (HBT for EM, etc.), measure each
   force's flow axis (P3 method) and compare gravity directly (H-A/H-C: is gravity
   the one that sits differently?). Respect no-synthesis-across-construction-type.

## FOUNDATIONAL / hunch-linked

9. **INFO-039 promotion** (MED). dMI/dt as target, >=3 seeds on the off-attractor
   residual; separate opposition-beyond-equal-entropy from the equal-entropy identity.

10. **Biology JOINT time-shuffle** (LOW; off-gravity; tests H-B "flow=time").
    Same permutation both channels; does the coupling survive (instantaneous) or
    die (temporal)? NOTE: MI-in-null coupling reading is DISPROVED per Greg (S16) —
    re-scope before running.

## HOUSEKEEPING

11. **OD run / MASTER_DISCOVERIES.json** (confirm + do).
    Pin what the planned "OD run" was; most likely = store this session's
    discoveries (Z-propagator, QCD running, time-blind clarification) per the
    "every OD discovery added immediately" rule. MASTER_DISCOVERIES.json is not in
    this repo (lives on E:\).

12. **CLAUDE.md master merge** (the parked restructure).
    Repo CLAUDE.md is canonical only through S11; the master is S15; this is S16.
    Fold S12-S16 deltas into the repo CLAUDE.md so it is canonical-and-current,
    keep detail in handoff files, then only upload deltas going forward
    (the "new way" — see CLAUDE.md START-HERE block).

---

## STATUS NOTES (what's already settled — do NOT relitigate)
- MI-in-null coupling discriminator: DISPROVED (Greg, S16). Do not build on it.
- Equal-entropy attractor = equal-marginal-entropy geometric identity (settled S5+S10).
  The OPEN question is #1 (construction vs nature), not whether the identity holds.
- Time-blind: RESOLVED (S16). The pooled extraction is row-permutation-invariant by
  construction (definitional); the TIME-RESOLVED null carries strong time structure.
- WF + SF Piece-1: DONE (Z propagator M_Z 99.5%; QCD asymptotic freedom from data).
