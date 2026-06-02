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

## O3 (in progress) — GW170817 BNS chirp mass
Launched a probe to nail the GW170817 chirp mass (S17 got the law form R^2 0.88 but mass
biased 15.7 vs catalog ~1.20, lever-arm-limited in the 32 s H1 band) via pycbc matched-filter
on longer/L1 data. Deliverables `probe_o3_gw170817_chirpmass.*`; INFO-061 reserved. Fold the
result into the ledger + capability brief when it lands.

## Files this session
- `probe_gravity_time_dilation.py` + `_results.json` + `_canary.json` (6c-A, INFO-058).
- `probe_6c_gps_positive.py` + `_route2_obs.py` + `_results.json` + `_canary.json` (6c-B,
  INFO-059).
- `probe_pulsar_time.py` + `_results.json` + `_run.log` (6c-C, INFO-060).
- `CLAUDE.md` (header S18, ledger INFO-058/059/060, Capability Demonstrations + Session 18
  note + handoff pointer), `CAPABILITY_BRIEF.md`, `BACKLOG_tests_and_probes.md` (6c marked).
- Raw data under gitignored `data/gps/`, `data/pulsar/` (+ `data/ligo_bulk/` for O3).

## Env notes (do not persist across containers)
- Pulsar route: `pip install pint-pulsar pdfminer.six`, `pip install --force-reinstall cffi`
  (broken container cryptography binding), certifi CA-bundle patch for PINT clock downloads.
- O3 route: `pip install pycbc --ignore-installed cryptography`.
- numpy 2.4 removed `ndarray.ptp()` -> use `np.ptp(...)`.

## Next (backlog; Greg picks; clear backlog before new probes; stop after each probe)
- O3 result (in flight) -> fold + possibly promote.
- 6c follow-ups: replicate the GPS positive recovery across stations/days + dual-frequency
  ionosphere-free combination to recover the circular-GPS sats; GW170817 absolute mass (O3).
- Remaining backlog: #2 EM/HBT construction control (Zenodo 5113016), #7 SF femtoscopy-R,
  #9 INFO-039 promotion, housekeeping #11 (MASTER_DISCOVERIES).
