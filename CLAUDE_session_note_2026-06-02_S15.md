## Note (Session 15 — 2026-06-02) — Gravity / Flow-Dipole / Time / Black-Hole exploration

Single-session CLAUDE note, in master format, ready to merge. Branch
`claude/file-upload-session-context-pOsLp`; main untouched; no PR. Opened by
installing the Session 14 master `CLAUDE.md` + v14 handoff + v15 kickoff into
memory, then a free-form exploration. Findings below are stated at DATA level only
(operator's request: no frame-grading this session — fresh eyes to judge). All
Operating Rules in force.

### Operator hunches raised (conjectures, recorded WITHOUT evaluation)
- **H-A** gravity = mixture of "pure physics" + "flow" dipole; affects time/flow.
- **H-B** the flow dipole expresses itself AS time.
- **H-C** EM/strong/weak don't directly affect time, only gravity does ("something
  in the gravity equation cancels the time equation out").
- **H-D** a black hole turning objects into info = zeroing the dipole equations back
  to the bare substrate/dipole level (no-hair analogue: settled BH = mass/charge/
  spin only).
- **H-E** if flow = equal-entropy substrate, that's why nothing (phys/bio/chem/geo)
  exists in a black hole.
- **H-F** (clarified, one consistent model) existence gradient: high cos→substrate
  (~0.9) = exists; distance away (gravity ~0.6) = toward non-existence; substrate =
  existence pole, distance = the gradient.
- **H-G** (newest, unprobed) the SCATTER in per-event gravity coordinates may itself
  carry structure.
Anchored physics (re-verify, treat as conjecture): gravitational time dilation is
measured/replicated (GPS, optical clocks, Pound-Rebka); static EM/strong/weak
potentials don't directly dilate a clock; other forces touch time only via energy
sourcing gravity.

### Probes run + data-level findings (no interpretation)
- **INFO-S15a (P1, `time_flow_probe.py`)** — sim physics+biology, 3 seeds. Level
  null is exactly time-ORDER invariant (cos 1.000000; partly definitional, row
  covariance). With rate ops [dH/dt, dMI/dt] added, null weight is 99.2–99.7% on
  LEVELS, 0.3–0.8% on RATES. Constraint is a static level relation.
- **INFO-S15b (P2, `grav_time_retaining.py`)** — real LIGO, 3 events (GW150914
  cleanly aligned). Order-shuffle cos 1.000000 all events (time-blind on real data).
  GW150914: MI peaks AT merger (16.41 s, kurtosis 18.25); dMI/dt peaks 16.47 s; |MI
  coef in null| 0.001 (self-pole). Reproduces INFO-036/045.
- **INFO-S15c (P3, `flow_dipole_axis.py`)** — FIRST measurement of the flow-dipole
  axis (d(MI)/dt regressed on levels). Opposition signature present fraction 1.0 all
  4 sim domains; flow axis cos→equal-entropy substrate 0.943–0.995; flow R^2 weak
  0.003–0.13. Biology flow axis (0.969→substrate) ≠ its coupling axis (0.588→MI).
  Gravity whole-segment: R^2 0.024, cos→substrate 0.608. Measures by regression what
  INFO-039 identified structurally (opposition ≈ equal-entropy identity).
- **INFO-S15d (P4, `grav_flow_reproduce.py`)** — gravity full/event/noise +
  bootstrap. GW150914 event window R^2 0.438, cos→substrate 0.668 (noise 0.780);
  bootstrap [0.31, 0.95] (loose); off-axis component also in noise. Other 2 events
  alignment failed.
- **INFO-S15e (P5, `grav_flow_crossevent.py`)** — fixed alignment on known merger,
  3 events. Event-window cos→substrate 0.781 / 0.639 / 0.994; GW150914 event ≈ its
  noise (0.781 vs 0.792); cross-event axes scatter (GW150914–GW170104 |cos| 0.106).
  Off-substrate flow axis does NOT reproduce across events; consistent with INFO-038/
  046 detector-state. Per-event scatter unexplained (H-G).

Method-sensitivity audit (GW150914 flow cos→substrate): 0.608 (whole seg, R^2 0.024)
→ 0.668 (event, peak-MI, R^2 0.44) → 0.781 (event, known merger, R^2 0.46; = noise).
Value not yet pinned; cross-event reproduction not established (1 clean event; other
two weak detections, MI_event/noise ≈1.07–1.08).

### Fits vs does-not-fit preexisting data
FITS: level-null time-blindness; flow axis ≈ equal-entropy substrate + opposition =
identity (measures INFO-039); gravity self-pole + MI-peaks-at-merger (INFO-036/045);
off-substrate value present in noise + event-varying (INFO-038/046).
NEW/does-not-fit a prior frame: biology flow axis ≠ coupling axis (new distinction);
gravity event-window off-substrate displacement loosely determined, = noise on clean
event, non-reproducing across events; per-event coordinate scatter unexplained.

### Queued probes (PROBE 1 first)
1. **(FIRST) dMI/dt rate-DYNAMICS** — test lagged/temporal rate structure the
   time-blind level null discards; gravity vs sim within-construction. Gate: is
   event-window sampling enough for a rate-dynamics fit?
2. **Construction-type control** — measure a non-gravity REAL detector-pair object
   (EM/HBT, INFO-048) identically; on-substrate (~0.9) or off (~0.6–0.8)? Needs HBT.
3. **Louder-event reproduction** — GW170814/GW190521 (4096 s hdf5 archives,
   ~134 MB, sliceable); re-run P5 for ≥3 clean events.
4. **Biology JOINT time-shuffle** — same permutation both channels; MI≈0.28*H_a
   survives (instantaneous, H-B unsupported) or dies (temporal, H-B supported)?
5. **Inspiral→merger→ringdown trajectory** — cos-to-attractor(t) shape
   substrate→off→substrate (H-D)? Gate: INFO-046 post-merger is detector-noise at
   125 ms/35–350 Hz; needs finer/lower band.
6. **Scattered-coordinates structure (H-G)** — does per-event spread track SNR/mass/
   noise-state/alignment-lag? Needs ≥ several clean events (depends on #3).
7. **INFO-039 promotion** — dMI/dt as target, ≥3 seeds on off-attractor residual;
   separate opposition-beyond-equal-entropy from the identity.
8. **Gravity vs gauge forces within-construction** — once weak/EM/strong real
   objects on a branch, P3-method flow axis each; compare gravity directly (H-A/H-C).

### Files (branch `claude/file-upload-session-context-pOsLp`)
Scripts: time_flow_probe.py, grav_time_retaining.py, flow_dipole_axis.py,
grav_flow_reproduce.py, grav_flow_crossevent.py. Results: matching *_results.json
(+ time_flow_probe_canary.json). Data: data/ligo/ (H1/L1 GW150914 + GW170104 +
GW151226, 32 s/4096 Hz hdf5). Env: h5py pip-installed for the LIGO probes (add to
requirements.txt if persisting). Companion docs: SESSION_HANDOFF_2026-06-02_FULL.md,
SESSION_HANDOFF_2026-06-02_gravity_flow_time_EXPLORATORY.md.
