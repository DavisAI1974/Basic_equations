# SESSION HANDOFF — 2026-06-03 (S19) — flow-dipole equation: wrong tool, and the OTHER dipole solves for time

Branch `claude/kickoff-claude-handoff-0vfgz` (main untouched; no PR). All Operating Rules
in force; no new Rule this session.

## Opening housekeeping (file continuity)
- The session branch was cut from a STALE point (main at S12 kickoff; CLAUDE.md header
  "Session 7", no S13-S18 work present). Fast-forwarded it onto the full S18 work branch
  `claude/kickoff-handoff-sequence-4jk6N` (38 commits, clean ancestor — no divergence),
  pushed. All four uploaded continuity docs verified byte-identical to the repo S18
  versions. Continuity restored before any probe.

## Task: the FLOW DIPOLE EQUATION (S18 backlog #1), then Greg's live steers
Confirmed the paper flow form from davisai.ai/dipole §2.2 verbatim:
  dMI_total/dt ~ sum_i c_self,i*H_i^2 + sum_{i<j} c_cross,ij*H_i*H_j + linear
plus the algebraic ratio C = H_self/H_cross (H_self = internal Shannon entropy, H_cross =
MI), opposition = self-coeffs vs cross-coeffs opposite sign.

### Epistemic recalibration (Greg, load-bearing this session)
A run of corrections that govern how the results below are read:
- Don't be determinative about one test; each result is one data point.
- Not all data points are equal: **4 consistent tests are a stronger signal than 1
  outlier**; a lone dissenter is more likely a weak/wrong candidate than hidden gold —
  do not romanticize the outlier, but inspect it (no tent-widening).
- A candidate that only shows up in 1 of N systems is, by that inconsistency, a POOR
  candidate (Base-of-Structure: a real base shows up consistently).
- "First try, not only try": a uniform low result across systems is at least as much
  evidence the TOOL is wrong as that the structure is absent.
- Category discipline: chemistry is NOT a 5th force — the four forces are gravity/EM/
  weak/strong; sciences/domains are a separate inquiry. Don't blend them.

## ARC OF RESULTS

### 1. The flow dipole is BLIND on real data (and the toy result was an artifact)
- **Toy simulators retired** (Greg: "why are we using toy brusselator? ... untrustworthy";
  S10 "no fake data"). `probe_flow_dipole_4domain.py` + `_brusselator.py` marked SUPERSEDED
  (audit trail kept). Their apparent shared opposition signature did NOT survive real data.
- **Real-data flow dipole, first tool**: gravity (LIGO H1/L1) R2_full 0.022 ~= shuffle-null
  0.005; weak (CMS dimuon, two muons, axis=invariant mass, 66 quantile bins) R2 0.044 ~=
  shuffle-null 0.065. Flat. (`probe_flow_dipole_gravity.py`, `probe_flow_dipole_weak.py`.)

### 2. Tool-batteries: the flatness is TOOL-BLINDNESS, not absence of structure
Ran 6-7 different extraction tools per force on the SAME real data, each with a KNOWN-
correlation diagnostic that demonstrably exists in the data.
- **STRONG** (`probe_flow_dipole_strong_battery.py`, real CMS two-particle/dimuon-quarkonium):
  1 of 6 tools hit — the KNOWN femtoscopy C(q) (low-q correlation vs event-mixed control,
  ratio 0.80). ALL entropy/dipole/symbolic/direct-coupling tools blind (R2 <= 0.09 at/below
  null). Data caveat: no hadron-pair/dijet CSV in CMS open data -> dimuon used.
- **EM** (`probe_flow_dipole_em_battery.py`, real HBT photon timetags, Zenodo 5113016):
  the known g2(0)=1.85 HBT bunching FIRES; the RAW-count covariance dipole HITS (R2 0.49 vs
  null 0.21) and direct coupling hits (Pearson 0.142, ~250 sigma); but the entropy->MI->
  dMI/dt operators are BLIND (R2 0.11-0.13 ~= null). MECHANISM (the key diagnosis): the
  windowed marginal entropies H_a,H_b are essentially independent of the inter-channel TIME
  LAG, so coupling/time info never enters the entropy operators — it lives in the RAW
  covariance. Caveat: the dataset's "uncorrelated" control also shows g2~1.9 (not a clean
  no-correlation control); within-script time-shift/shuffle nulls are the load-bearing ones.
- READING: on the two forces with a known real correlation to check against, the dipole/
  entropy tool is provably BLIND to a correlation that exists. "Wrong tool" over "no structure."

### 3. Greg's pivot: don't discard the positive TIME findings; try the OTHER dipole to solve for time
Two dipole forms beyond the differential/flow one (precedent in `static_dipole_test.py`:
DNA gave R2=0 differential but preserved the dipole as an ALGEBRAIC constraint):
  (A) RAW cross-covariance dipole (over lag, for time-series channels).
  (B) STATIC ALGEBRAIC dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2.

- **GRAVITY** (`probe_time_dipole_gravity.py`, LIGO H1/L1 GW150914):
  (A) raw cross-cov dipole SOLVES FOR TIME — recovers the inter-detector light-travel lag
      7.32 ms (|physical| ~6.9 ms), peak|cc| 0.64 vs off-source null 0.17+/-0.03, z=14.2;
      the ENTROPY-MI dipole on the SAME data is lag-FLAT (flatness 0.019) = blind to time.
  (B) static algebraic dipole R2 0.93 (event) vs the flow dipole's 0.022.
- **GRAVITY tautology-killing null** (`probe_time_dipole_gravity_null.py`): circular-shift of
  H_b preserves H_b smoothness + the shared-H_a factor, kills only the instantaneous
  H_a<->H_b pairing. EVENT R2 0.93 vs shift-null 0.80 -> excess +0.125, z=2.7, p=0.000,
  corr(Ha,Hb)=+0.73 (GENUINE event-specific joint structure = the common GW signal coupling
  both detectors' entropies). NOISE R2 0.71 vs shift-null 0.72 -> excess ~0, z=-0.7,
  corr -0.06 (PURE shared-H_a tautology). So the algebraic dipole carries a real (modest)
  event signal surviving a proper null; the noise apparent-dipole is tautology.
- **GPS** (`probe_time_dipole_gpspulsar.py`): decision-gate PASS (genuine 2-channel: A =
  cleaned RINEX pseudorange relativistic residual, B = independent broadcast-element geometry
  driver c*F*e*sqrt(a)*sin(E) — NOT the S18 regression relabeled). Raw cross-cov dipole:
  eccentric Galileo E18 lag-0 |cc|~0.99 z=37, E14 z=22; near-circular control G15 |cc|~0.44
  lag-12 (no lag-0) — relativistic coupling recovered, control fails. Static algebraic + shift
  null: E18 R2 0.989, excess +0.571, z=5.7, corr(Ha,Hb) 0.97; control G15 excess +0.374,
  z=1.5, corr 0.17 — signal beats tautology null far more than control.
- **PULSAR DECLINED — category stretch** (the chemistry lesson applied): PSR B1913+16's
  dP_b/dt and gamma are SCALAR parameters of a global PINT WLS fit; no co-sampled 2nd
  PHYSICAL channel, the secular dP_b/dt is a cumulative parabola not an instantaneous A<->B
  coupling; a synthesized channel-B would make any fit a model tautology. Clean negative, no
  pulsar dipole run, no fabrication.

## NET READING (three levels, deflationary brakes on)
- DATA: the differential/flow entropy dipole is flat on every real force; tool-batteries show
  this is blindness (known correlations C(q), g2 fire; raw-covariance dipole hits). The
  RAW-covariance and STATIC-ALGEBRAIC dipoles travel to TWO independent gravity-time systems
  (LIGO 7ms z=14; GPS -2/c^2 coupling z=37) and the algebraic form survives a tautology-
  killing null (LIGO event excess +0.125 z=2.7; GPS E18 +0.571 z=5.7), where the noise window
  and the controls do not.
- INTERPRETATION: the flow dipole was the WRONG TOOL — it discards the inter-channel covariance
  where time/coupling lives. The positive time findings were not lost; we were measuring them
  with the wrong dipole. "Time is something else" reads as: time is carried by the raw/algebraic
  dipole, not the flow dipole.
- FRAME (ungraded): these are RECOVERIES OF KNOWN physics (7 ms light-travel time; GR -2/c^2) —
  a positive CONTROL that the dipole tool travels, NOT new physics. Stats are modest where data
  is thin (GPS one-station/one-day/~10 windows). The toy "4 forces share the dipole" claim is
  retired. Whether the algebraic dipole is a real universal law vs a positive-control recovery
  is NOT settled — that is what the in-flight replications test.

## FOLLOW-UPS (a) + (b) -- BOTH COMPLETED
- (a) STRENGTHEN GPS (`probe_time_dipole_gps_strengthen.py`) -- DONE. Multi-station/multi-day/
  dual-frequency. Eccentric-Galileo algebraic excess-over-tautology-null = mean +0.464+/-0.278
  (z~3.0) across 12 sat-replicas / 2 stations (BRUX, ALGO; ONSA/MATE/GRAZ never see the GREAT
  eccentric sats -- real GNSS visibility) / 2-3 days; +0.572+/-0.134 on the clean days = the S18
  BRUX/001 anchor (+0.571) REPRODUCED as a multi-station/day mean. Raw cross-cov |cc0| 0.787+/-
  0.327, z~22 across 13 replicas. k/truth +1.14+/-0.30 (recovers GR -2/c^2). Dual-frequency P3 did
  NOT rescue the 17-33 m circular-sat controls -- they STAY clean low-excess controls (good
  control behavior; the tiny signal is below the residual/receiver-clock floor even iono-free).
  Honest outlier KEPT (no tent-widening): day-003 degrades E18 at BOTH stations identically under
  C1 and P3 -> a real station-day condition (receiver-clock/geometry), not a frequency artifact;
  E14 still recovers. NET: the GPS gravity-time dipole is a robust AGGREGATE, not a one-day fluke.
- (b) DOES THE OTHER DIPOLE TRAVEL TO THE GAUGE FORCES (`probe_time_dipole_forces.py`) -- DONE,
  mostly NO. Static algebraic dipole + circular-shift tautology null: EM excess +0.085 z=2.0
  (and its 'uncorrelated' control shows a LARGER excess +0.294 -> contaminated control, inverts
  the logic), strong +0.057 z=1.4 -- both COLLAPSE toward tautology like the gravity NOISE window;
  weak borderline +0.114 z=2.5 but fragile (52 near-resonance bins, no usable same-charge control).
  The RAW cross-cov-over-lag form DOES fire on EM (z=420 lag-0 = known HBT bunching) but applies
  only to time-series channels (dimuon ensembles have no lag). The gravity analog (corr rises where
  the correlation is strong) does NOT reproduce on the forces. NET: the algebraic dipole's genuine
  excess is confined to the 2 gravity-TIME systems (LIGO event, GPS eccentric) -> reads GRAVITY-
  TIME-SPECIFIC, consistent with the H-C hunch -- with brakes (weak fragile, thin samples,
  recoveries of known physics).

CONSOLIDATED S19 PICTURE: flow/entropy dipole = wrong tool (blind, INFO-065); the OTHER dipole
(raw cross-cov + static algebraic, tautology-null-survived) carries the gravity-TIME coupling on
TWO independent systems (LIGO + GPS, now multi-station-aggregated) and does NOT travel to the
gauge forces -> gravity-time-specific. All RECOVERIES of known physics (positive control the tool
travels), not new physics.

## FILES THIS SESSION
- Real-data flow dipole: probe_flow_dipole_gravity.py, probe_flow_dipole_weak.py (+ results).
- Tool-batteries: probe_flow_dipole_strong_battery.py, probe_flow_dipole_em_battery.py (+ results).
- OTHER dipole / solve-for-time: probe_time_dipole_gravity.py, probe_time_dipole_gravity_null.py,
  probe_time_dipole_gpspulsar.py (+ results/canaries).
- Deprecated (audit trail): probe_flow_dipole_4domain.py, probe_flow_dipole_brusselator.py.
- In flight: probe_time_dipole_gps_strengthen.py, probe_time_dipole_forces.py.
- Raw data (gitignored): data/em/ (HBT), data/strong/ (dimuon CSV), data/gps/, data/chem/ (BZ,
  fetched then parked — chemistry is off-category for the force inquiry).

## NEXT (Greg picks; stop after each probe)
- Land (a)+(b); read whether the algebraic dipole TRAVELS to the gauge forces or is gravity-time-
  specific (decides "universal dipole" vs "gravity-time positive control").
- If it travels: a proper cross-system comparison of the algebraic-dipole coefficients/excess
  (substrate form vs per-force expression) — the original decisive test, now on the RIGHT tool.
- Chemistry/sciences inquiry is a SEPARATE line (BZ data parked in data/chem/), not the forces.
