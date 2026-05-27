Resuming Information Layer work, Session 11.

Branch: `claude/unused-tools-review-r3CJV` (Session 10 closed here; on origin).
main untouched. No PR.

Session 10 opened the prebiotic-chemistry-to-biology line. Ran the per-domain stack on three new simulator systems (hypercycle, quasispecies, asymmetric hypercycle) and on PySR cross-seed test. Produced four candidate INFO entries (INFO-037 through INFO-040). Surfaced a new working frame at session close.

The new frame is why this session is starting fresh: **"Time as expression of dipole flow"** (Greg, Session 10). Frame stays speculative until dedicated probes run, but the operational connection to Session 10 data is real and worth handling cleanly.

Before doing anything else:

1. Read **CLAUDE.md** (project root). Session 10 note at top alongside Sessions 3-9 notes. Read all seven Operating Rules (no pre-assigned meaning, probe-not-falsifier, speaking posture before/after, Rule D incomplete-not-wrong, They never stacked, Treat literature as conjecture). Read the "Time as expression of dipole flow" working frame in the Working Frames section. Hold all of this through this session. **NO new Operating Rule was added Session 10.**

2. Read **SESSION_HANDOFF_2026-05-27_v10.md** (project root). Full Session 10 record: five experiments, INFO-037 through INFO-040, the substantive close where the new frame surfaced.

3. (Optional) Read **SESSION_HANDOFF_2026-05-26_v9.md** for the OD platform + real-data probe context. Session 10 did not touch the platform; the simulator line and the real-data disease-detection line run in parallel.

Speaking posture from Rule C applies from your first response: "I think X might happen, but we'll wait on what the data says and where it points us." No verdict in advance. Rule D applies (incomplete-not-wrong). Rule "literature as conjecture" applies. "They never stacked" applies.

## What we know going in

The Session 3 universal attractor identity is `-(H_a - H_b)^2 ≈ 0` — channel-symmetric, no preferred direction. Every system tested across Sessions 3-9 (8+ heterogeneous dynamical systems, real cardiac data, real cerebrovascular data, simulator four-domain, four-force toy) lands on this attractor at cos ≥ 0.99 EXCEPT — discovered Session 10 — channel-asymmetric prebiotic systems.

Session 10 produced two reproducibly-off-attractor systems:
- Quasispecies with selection (smooth walk in cos_attr from 0.056 to 1.000 across mu)
- Asymmetric hypercycle WITHOUT selection (apparently sharp transition at base_A ~ 1.5-2.0)

Reading: channel asymmetry per se — not selection — is sufficient for off-attractor behavior. But the manner of the transition may differ between selection-driven (smooth) and base-rate-driven (sharp) asymmetry.

Greg's new frame says: read on/off attractor as no-time-arrow / time-arrow-exists. The off-attractor regime is where the dipole has charge separation, where flow has a direction. The frame is consistent with the data but the data does not prove it. Two probes (top of queue) test the frame specifically.

## Top of queue for Session 11

Greg's call to make. Default order is by frame priority + cost; pick or redirect.

**(1) Time-reversal probe on quasispecies.** Take a quasispecies trajectory at low mu (off-attractor, e.g., mu=0.001) and at high mu (on-attractor, e.g., mu=0.5). Reverse the time series in place. Recompute the operator matrix from the reversed series. Compute the new v_null6 and cos_to_attractor.

Prediction from the frame (held as prediction, not as expectation): if the off-attractor signature IS time-direction-coded, the reversed low-mu system lands on a predictably-different operator-space direction — ideally the sign-flipped twin on the asymmetric subspace (H_b - H_a swaps sign). The on-attractor high-mu system, being symmetric, should be invariant under reversal.

Operating Rule "no pre-assigned meaning": this is the prediction the frame WOULD make; we run the probe and look at what comes out. If the prediction holds, frame promotes to candidate. If it doesn't, frame is incomplete-or-wrong and we map why.

Cost: ~15 min including write + run. Code already has simulate_quasispecies and build_ensemble_operator_matrix. Add a time-reversed variant of the operator matrix builder. Three seeds at each mu, mu in [0.001, 0.05, 0.5] (off, transition, on).

**(2) Entropy-production rate correlation.** Compute dS/dt along each trajectory in the existing quasispecies + asymmetric hypercycle sweeps. Quantitative test of the link between dipole charge separation (|H_a - H_b|, off-attractor distance) and irreversibility (dS/dt).

Operating posture: if the frame is right, off-attractor distance should correlate with dS/dt across both sweeps. If it doesn't correlate, the frame's quantitative link is weaker than the qualitative one.

Cost: ~10 min. Read existing JSONs, recompute dS/dt from the time-series, plot/report correlations.

**(3) Densify asymmetric-hypercycle base_A grid** between 1.5 and 2.0. Characterize whether the cos_to_attractor transition is genuinely sharp (threshold/phase transition) or whether the apparent sharpness is grid-sampling artifact.

Cost: ~10 min. Same script as Exp 4, finer grid (e.g., 1.5, 1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 2.0).

**(4) Quasispecies with asymmetric base rates** (compose selection + base-rate asymmetry). Vary f_A AND base_A independently. Tests whether the two off-attractor mechanisms compose linearly, cancel, or stay distinct.

Cost: ~15-20 min.

**(5) PDG four-force real-data probe** (universe-origin side, deferred from Sessions 8/9). Pull PDG running of g1, g2, g3 with energy. Apply per-domain stack. Tests whether real EM/weak/strong behave like the Session 8 toy caricatures.

Cost: 30-45 min including data download.

## Also queued (Session 9-10 carryovers)

(6) Compartmentalization probe (protocell boundary autocatalytic system).
(7) Tools Greg attaches in next session (carried from v10 kickoff queue).
(8) AF cohort full fetch + plug into store (carried from v10).
(9) Active-learning hook in od_learning_store (carried).
(10) Online classifier swap (SGDClassifier with partial_fit).
(11) Cross-problem transfer (cardiac warm-start to sepsis).
(12) Storm/waves real-data probe (NOAA buoy / HRRR atmospheric).
(13) MIMIC credentialing (durable infrastructure).
(14) Five-physics-areas substrate hunt (Session 5 reorientation, still open).
(15) Complexity-entropy plane (Rosso 2007).
(16) Yamada 2023 multi-field SVD with entropy as field.
(17) Gomez-Herrero 2015 ensemble infrastructure.
(18) Jacobson 1995 cross-domain extension.
(19) Stochastic Electrodynamics 21st-century revisit.

## Files to know about (Session 10 outputs in repo root)

Stack (copied from `claude/two-more-tasks-iZvY4` Session 6-7):
- kbk_pipeline.py        Vasicek-H + KBK rank-gap + projection
- ai_poincare_rank.py    AI Poincare local-PCA / two-NN / MLE
- sindy_symbolic.py      SINDy with extended library
- per_domain_kbk.py      Four-domain simulators + ensemble-H builder
- pysr_symbolic_per_domain.py  PySR runner template

Session 10 simulators + probes:
- prebiotic_canary.py              Hypercycle canary (1 seed, k_cat sweep)
- prebiotic_multiseed.py           Hypercycle 3 seeds at Session 7 baseline
- quasispecies_probe.py            Eigen quasispecies (mu sweep) — first off-attractor system
- asymmetric_hypercycle_control.py Asymmetric hypercycle, NO selection (base_A sweep)
- prebiotic_pysr.py                PySR cross-seed test

Result JSONs / logs paired with each.

## Platform code (NOT pulled in Session 10)

The OD learning platform (od_features.py, od_learning_store.py, od_ingest_patient.py, od_predict_patient.py, od_learn_cerebro.py, store/) lives on `claude/two-more-tasks-iZvY4`. Session 10 did not pull these in because the simulator probes did not need them. If Session 11 work requires platform features, pull them in via `git checkout origin/claude/two-more-tasks-iZvY4 -- <files>` as Session 10 did for the stack scripts.

## Dependency notes

`numpy`, `scipy`, `scikit-learn` installed. `pysr` installed (Julia 1.x backend via juliapkg). First PySR fit per session is slow (Julia JIT); subsequent fits ~5s each.

## What's at stake operationally

If the time-reversal probe confirms the frame, Session 10's INFO-037 through INFO-039 promote from "channel asymmetry signature" to "time-arrow signature." That changes the reading of every off-attractor result we've ever seen — including the cardiac AF outlier from Session 9 (which sat off the cardiac substrate cluster) and the damped oscillator from Session 3 (which sat off the universal attractor with cos 0.69-0.74).

If the probe doesn't confirm, the frame is incomplete-or-wrong (Rule D), the off-attractor signature remains a "channel asymmetry signature" — still substantive but without the arrow-of-time interpretation, and we stay on the queued mechanism-distinction work (densification + selection-vs-base-rate composition).

Either way, we run the probe and look at what the data says.

We are still pioneers. All seven Operating Rules in force. Hold them through Session 11.
