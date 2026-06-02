# SESSION HANDOFF — 2026-06-02 (S18) — gravity-couples-to-time (backlog 6c)

Branch `claude/kickoff-handoff-sequence-4jk6N` (main untouched; no PR). All Operating
Rules in force; no new Rule this session.

## Opening housekeeping (file continuity)
- Folded the S12-S17 deltas into the canonical `CLAUDE.md` (header bumped to S18, body
  canonical through S18) and brought the full S17 work branch
  (`claude/claude-md-problem-nM8Hm`) onto this branch via merge — all probe scripts +
  result JSONs + intermediate artifacts (S14 master, S15/S16 handoffs, LIGO/CMS data).
  The 5 uploaded continuity docs were byte-identical to that branch's versions.
- Added `data/gps/`, `data/pulsar/`, `data/ligo_bulk/` to `.gitignore` (reproducible raw
  downloads; results JSON committed).

## Backlog 6c — does gravity couple to TIME in a recoverable governing law?
The H-C hunch. INFO-056 (S17) showed force-coupling data is construction-confounded
(gravity is our only time-series observable), so the clean test is clock / time-dilation
data. Greg steered: run 6c first, then O3; for 6c run the GPS-positive and pulsar routes
IN PARALLEL. Answer: YES — three independent recoveries.

### PROBE 6c-A — GPS precise-product route: definitional NULL (INFO-058)
`probe_gravity_time_dilation.py` on CODE final products (`data/gps/COD_20230010000_05M_ORB.SP3`,
GPS+Galileo, 5-min orbit+clock, 2023 DOY 001). Target: the relativistic clock offset
`dt_rel = -2(r.v)/c^2 = -2(r*rdot)/c^2` (the identity `r.v = r*rdot` is frame-invariant, so
it comes from SP3 position magnitudes alone). Joint per-sat fit `clock(t) = poly3(t) +
k*(r*rdot)`; ground truth `k = -2/c^2 = -2.2256e-17 s^2/m^2`.
- RESULT (57 sats): recovered k = 0.4% of truth; median residual 0.10 ns vs orbit-predicted
  term 24-385 ns (eccentric Galileo GREAT E14/E18 ~275-385 ns). NOT buried in noise (a cubic
  cannot absorb ~2 cycles/day; 0.05 ns << 275 ns).
- READING: the GR time-dilation term is MODELED OUT of IGS/CODE precise clock products
  (standard processing applies it in the observation model), so recovery on the cleaned
  product is a definitional NULL. Clean null with a clear reason; kept as data
  (don't-predetermine). Motivated the term-retaining route.

### PROBE 6c-B — GPS term-retaining route: POSITIVE recovery of -2/c^2 (INFO-059)
`probe_6c_gps_positive.py` (+ helper `_route2_obs.py`). Route 2 (raw observations): raw
RINEX C1 pseudorange (IGS station BRUX, Brussels, 30 s, 2023-001) - geometric range (SP3
orbit + known station ECEF + Sagnac + travel-time) - precise SP3 sat clock (term removed) -
per-epoch receiver clock (estimated from near-circular Galileo sats) leaves the relativistic
modulation in the residual. Fit `cleaned_resid = poly2(t) + k*[-c*e*sqrt(a)*sin(E)]`; the
`e*sqrt(a)*sin(E)` regressor is computed INDEPENDENTLY from broadcast Keplerian elements, so
high R^2 is genuine cross-source agreement (non-circular).
- RESULT (eccentric Galileo GREAT sats, e~0.162, ~275 ns / ~200 m signal): E18 k/truth
  +1.038, R^2 0.978, corr -0.987, z=380 sigma vs scramble null; E14 k/truth +1.019, R^2
  0.995, corr -0.997. Recovered F = -2*sqrt(mu)/c^2 directly.
- TWO DEAD ROUTES mapped as data (route1_null in JSON): precise-product NULL (= INFO-058);
  broadcast-differencing NULL — precise SP3 clock and broadcast af0/af1/af2 polynomial agree
  to ~3 ns over the day while the relativistic term is ~760 ns ptp, so it is in NEITHER (why
  raw observations were necessary).
- CAVEATS: clean only on the 2 eccentric Galileo sats; the high-e GPS sats (e~0.015-0.025,
  signal 17-33 m) scatter (k/truth -2.3 to +1.9) because single-frequency C1 carries
  uncorrected ionosphere comparable to their tiny signal — diagnosed, not absorbed; a
  dual-frequency ionosphere-free combination would clean them. One station/one day -> LOCATED;
  replicate across stations/days before promotion. Sign is a pseudorange convention folded
  into the regressor; the physics is magnitude + form + null separation.

### PROBE 6c-C — Pulsar route: independent confirmation (INFO-060)
`probe_pulsar_time.py`. PSR B1913+16 (Hulse-Taylor), 9261 raw Arecibo TOAs 1981-2012
(Weisberg & Huang 2016, Zenodo doi:10.5281/zenodo.54764; RAW topocentric TOAs, NOT a
digitized figure). PINT (validated TEMPO successor) does only the standard clock/DE405-
barycenter/DM/Keplerian reductions (no dP_b/dt or gamma assumption); relativistic params
recovered by weighted least squares. GNSS-convention-independent.
- EXP1 orbital-decay detection (PBDOT-fixed delta-chi^2): PBDOT=0 wRMS 292 us (chi2 3.15e8)
  -> PBDOT free wRMS 26.98 us (chi2 2.69e6); Delta-chi2 = 3.13e8 (~17684 sigma). PBDOT-fixed-
  at-GR (27.03 us) reproduces the free fit -> data GR-consistent.
- EXP2 direct dP_b/dt: recovered -2.4151e-12 vs GR -2.40263e-12 (ratio 1.005; the ~0.5% is
  the known galactic-acceleration term -0.025e-12 required for GR-consistency), bracketed by
  observed -2.42297e-12.
- EXP3 Einstein delay gamma (gravitational redshift + 2nd-order Doppler, the most direct
  "gravity slows the clock" parameter): recovered 4.30737e-3 s vs published 4.30675e-3
  (ratio 1.00014, 0.014%).
- CAVEATS (Rule D, honest): PINT does the metrology (not from scratch); 9261 vs published
  9257 TOAs (tolerant re-parse of legacy fixed-column files); per-session JUMPs, high-order
  spin (F3..F10), Shapiro M2/SINI, red noise NOT all fit -> postfit wRMS ~27 us (paper 16.3)
  and FORMAL ERRORS ARE OPTIMISTIC; the large sigma_from_GR (91.6) is an artifact — the POINT
  ESTIMATES (ratios 1.005, 1.0001) are the robust result. DTHETA dropped (mis-scaled in .in).
  Published GR values treated as conjecture-to-check, not cited as support. An exp2 1e12
  units-reporting bug (PINT PBDOT.value is already dimensionless s/s; old *1e-12 was double-
  scaling) was caught and fixed; root cause verified by setting -2.4151e-12 -> .value
  -2.4151e-12.

### 6c joint reading (Result Discipline)
- DATA: three independent positive recoveries of gravity-time governing relations — GPS
  time-dilation coefficient -2/c^2 (z=380), pulsar dP_b/dt (ratio 1.005, ~17684 sigma), pulsar
  gamma (0.014%) — plus a clean definitional null (precise-product GPS) with a clear reason.
- INTERPRETATION: gravity's coupling to TIME is a recoverable governing law, via two distinct
  mechanisms (clock-rate dilation: GPS coeff + pulsar gamma; gravity altering orbital timing
  via GW emission: pulsar dP_b/dt), across two fully independent data domains and independent
  of GNSS processing conventions. This is the positive answer to the H-C hunch the S17
  force-coupling work could not reach (construction-confounded). NOT new physics — these
  RECOVER known GR laws (checkable against ground truth); that is the capability point.
- FRAME (ungraded): consistent with "gravity couples to time" as a real, recoverable law; the
  frame did not grade itself — the data carried it. The deeper "is gravity SPECIAL vs the
  gauge forces in touching time" remains a separate question (gauge forces have no clock-rate
  law, but a fair multi-force comparison is future work; EM/HBT not done).

## Capability brief
`CAPABILITY_BRIEF.md` + the CLAUDE.md "Capability Demonstrations" section updated to add the
gravity time-dilation recovery (GPS coeff + pulsar gamma + pulsar dP_b/dt) alongside the S17
inspiral-chirp recovery. Gravity now demonstrates TWO independent raw-data governing-law
recoveries.

## O3 — GW170817 BNS chirp mass: NAILED (INFO-061)
`probe_o3_gw170817_chirpmass.py`. S17 got the inspiral law FORM (u=f^-8/3 linear, R^2 0.88)
but biased the absolute mass (15.7 vs catalog ~1.20), lever-arm-limited in the 32 s H1 band.
O3 fetched GWOSC GWTC-1 GW170817 4096s/4096Hz H1+L1 (data/ligo_bulk/, gitignored), ran a
TaylorF2 chirp-mass template bank + Welch-PSD whitening + matched filter, gated the L1 glitch
(t_merger -1.05 s) with pycbc `.gate()`, took the peak complex-SNR within +/-0.1 s of the
merger GPS.
- RESULT: detector-frame M_c = 1.200 Msun (network), 0.19% from catalog 1.1977; H1 1.200
  (SNR 11.1, +0.014 s) and L1 1.195 (SNR 9.8, -0.008 s) INDEPENDENTLY agree; net SNR 14.8
  sharply peaked (drops to ~7.5 by +/-0.015 Msun). Null: L1 on 9.8 vs off-source 9.7 sigma
  (clean); H1 5.6 sigma (noisier H1-band off-source tail -- the S17 ridge ugliness -- but the
  on-source peak lands at the right M_c AND time).
- READING: matched filtering REMOVES the S17 ridge absolute-mass bias -> 0.19% recovery; the
  BNS chirp mass is a genuine gravity-domain governing-law (inspiral phasing) recovery on a
  second, physically-distinct event.
- CAVEATS (limit SNR, NOT M_c): TaylorF2 inspiral-only, no spin, q=1, f_final 1024 Hz, single
  PSD -> net SNR 14.8 << catalog ~32; inspiral phasing still fixes M_c. Source-frame ~1.187
  reached after z~0.0099 correction (not applied). Catalog = comparison target only. V1 not
  run (SNR ~2). Env: pycbc 2.11.0 silently shadows scipy with a broken 1.16.3 -> force-
  reinstall --no-deps scipy; use pycbc `.gate()` not a hand-rolled taper (the latter rings the
  matched filter).

Gravity now has THREE independent raw-data governing-law recoveries: inspiral chirp
(INFO-052), time dilation (INFO-059 GPS + INFO-060 pulsar), and BNS chirp mass (INFO-061).

## FLOW pivot (Greg) — flow=substrate, time=one expression (INFO-062/063/064)
Greg: "flow can have many expressions and time is one of them, just like quantum physics is
one expression of physics." Tested as substrate-vs-expression. Greg's epistemic rule
(load-bearing): a non-conforming system is ONE data point / possibly the wrong observable,
NOT a falsification of the frame or of that system.

- **INFO-062 (`probe_flow_substrate_expression.py`)**: SUBSTRATE = a characteristic running
  monotonically and DIVERGING toward a critical point. gravity chirp + strong alpha_s both
  flow (|rho|=1.000); WEAK Z (a real critical point at M_Z) does NOT (|rho|=0.054, peaks) =
  control beats "monotonic=trivial". Critical point real (gravity R^2 scan locks to true
  t_c). EXPRESSION: axis distinct (time/scale/mass); form open here.
- **INFO-063 (`probe_flow_battery.py`) -- ALL 4 FORCES**: 5 flows -- gravity (time), strong
  (scale), EM alpha_em running (scale, OPPOSITE sign), weak alpha_2 running (scale),
  Brusselator->Hopf critical slowing (CONTROL-PARAMETER axis, new domain, q=1.000 R^2=1.000).
  Controls weak-Z + driven-oscillator correctly non-flow. WEAK is flow (running) AND non-flow
  (Z resonance) -> "non-flow" was an OBSERVABLE CHOICE (Greg's wrong-observable point in
  data). EM was previously omitted because its program object was HBT (a correlation, not a
  flow); its flow is the alpha_em running. MISS kept as data: GW170817 H1 ridge (wrong
  observable). FORM distinct where measurable: gravity 3/8 vs Brusselator 1.0.
- **INFO-064 (`probe_flow_form_strong.py`)**: strong's divergence FORM (power vs log)
  UNDETERMINED -- the critical region Q->Lambda is nonperturbative/inaccessible (duality
  breaks <0.84 GeV; the measured coupling FREEZES to alpha_s(0)~0.76, so the Landau-pole
  divergence is a SCHEME ARTIFACT). |dR^2|<0.02 even at Q/Lambda~5.6. Legitimate undetermined,
  not a falsification.

NET: FLOW substrate confirmed + generalizes to all 4 forces + chemistry; TIME is gravity's
expression. Axis distinct; form distinct where measurable; axis UNIFICATION still not
data-forced (construction-confounded). No new Operating Rule.

## QUEUED NEXT PROBE (Greg) — the FLOW DIPOLE EQUATION on these flows
Each flow is currently a single characteristic vs an axis. A flow DIPOLE needs TWO coupled
channels, in the info-dipole paper's form:
  dMI/dt ~ sum_i c_self,i * H_i^2 + sum_{i<j} c_cross,ij * H_i*H_j + linear   (opposition signature)
Plan: pick the 2 channels per flow -- gravity = H1/L1 detectors (inter-detector MI already in
INFO-053); chemistry Brusselator = its native 2 species x,y (cleanest, start here); strong/
EM/weak = two observables. Extract the flow-dipole equation per system, then test substrate/
expression AT THE EQUATION LEVEL: a shared dipole FORM (substrate) with per-domain
COEFFICIENTS (expressions). Connect to the info-dipole paper (davisai.ai/dipole) + the
Markets flow dipole. This closes the loop on Greg's "time = flow dipole". IDEAL FIRST PROBE
for the fresh session.

## Files this session
- `probe_gravity_time_dilation.py` + `_results.json` + `_canary.json` (6c-A, INFO-058).
- `probe_6c_gps_positive.py` + `_route2_obs.py` + `_results.json` + `_canary.json` (6c-B,
  INFO-059).
- `probe_pulsar_time.py` + `_results.json` + `_run.log` (6c-C, INFO-060).
- `probe_o3_gw170817_chirpmass.py` + `_results.json` + `_canary.json` (O3, INFO-061).
- `probe_flow_substrate_expression.py` + `_results.json` (INFO-062).
- `probe_flow_battery.py` + `_results.json` (INFO-063, all 4 forces).
- `probe_flow_form_strong.py` + `_results.json` (INFO-064, strong form undetermined).
- `CLAUDE.md` (header S18, ledger INFO-058..064, Capability Demonstrations + Session 18
  note + handoff pointer), `CAPABILITY_BRIEF.md`, `BACKLOG_tests_and_probes.md`.
- Raw data under gitignored `data/gps/`, `data/pulsar/`, `data/strong/` (+ `data/ligo_bulk/`).

## Env notes (do not persist across containers)
- Pulsar route: `pip install pint-pulsar pdfminer.six`, `pip install --force-reinstall cffi`
  (broken container cryptography binding), certifi CA-bundle patch for PINT clock downloads.
- O3 route: `pip install pycbc --ignore-installed cryptography`.
- numpy 2.4 removed `ndarray.ptp()` -> use `np.ptp(...)`.

## Next (backlog; Greg picks; clear backlog before new probes; stop after each probe)
- **FIRST PROBE for the fresh session (Greg): the FLOW DIPOLE EQUATION on these flows** --
  see the "QUEUED NEXT PROBE" section above. 2-channel dMI/dt dipole form per flow; start
  with the Brusselator (native 2 species) + gravity H1/L1; test substrate(form)/expression
  (coefficients); connect to the info-dipole paper + Markets flow dipole.
- 6c follow-ups: replicate the GPS positive recovery across stations/days + dual-frequency
  ionosphere-free combination to recover the circular-GPS sats.
- Remaining backlog: #2 EM/HBT construction control (Zenodo 5113016), #7 SF femtoscopy-R,
  #9 INFO-039 promotion, housekeeping #11 (MASTER_DISCOVERIES).
