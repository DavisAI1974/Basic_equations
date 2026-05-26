# Information Layer Session Handoff — 2026-05-26 Session 4

## Read me first

Read in this order:
1. CLAUDE.md (project root) — has the v3 ledger plus Session 4 additions and
   the three new operating rules from this session.
2. This file — what happened this session, raw numbers, what's queued.
3. The Session 3 handoff at `E:\information_layer\SESSION_HANDOFF_2026-05-25_v3.md`
   (Greg's local) for the full Session 3 context that this session builds on.

Three operating rules were added to CLAUDE.md this session (Greg, 2026-05-26).
They apply to all future probe work:

- **Rule A** — Don't pre-assign meaning to potential outcomes before the data
  exists. Furthest is "I think this may happen, but I want to see what the data
  says and where it leads."
- **Rule B** — A probe is a generator of a different signal, not a falsifier.
  May not falsify anything; at worst points in a different direction.
  Avoid "falsifier" naming.
- **Rule C** — Speaking posture before AND after every probe: "I think X might
  happen, but we'll wait on what the data says and where it points us." No
  verdict in advance. No verdict on first look at output.

Greg's framing this session: we are the first ones here. Follow every road,
look under every stone. No safe assumptions. No known facts to fall back on.

## What this session did

### Part 1 — Session 4 baseline (N_REAL sweep on linear drift)

Greg ran the N_REAL sweep on linear drift per the top-priority Session 3
queued experiment. Files: `extract_v1_linear_drift.py`,
`run_NREAL_sweep.py`, `analyze_spectrum.py`,
`L1_linear_drift_NREAL_sweep.json`. All in repo root.

Protocol: extract_v1 — Vasicek 1D entropy, 6D operator basis
`[H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]`, w=40s, T=100, dt=0.1, 50% window
overlap (4 windows per realization), pool across N_REAL, center,
smallest right singular vector = null direction, project to indices
`[2, 3, 4]` for the 3D subspace.

System: `x_a(t) = v*t + eps_a`, `x_b(t) = v*t + eps_b`, v=1, sigma=0.1,
independent noise.

Sweep: N_REAL in {50, 100, 200, 500}, 3 seeds each. Twelve extractions.

Raw measurement:

- Inter-seed cos in 3D subspace: 0.79 at N=50, 0.85 at N=100, 0.985 at
  N=200, 0.99 at N=500. Monotonically non-decreasing.
- N_REAL=500 mean v3: (+0.466, +0.370, +0.804).
  Cos to (+1,+1,+2)/sqrt(6) = +0.9975. Cos to (-1,-1,+2)/sqrt(6) = +0.31.
- SV spectrum at all 12 cells: 3 large SVs of order 10^0, 3 small SVs of
  order 10^-4 to 10^-5, SV3/SV4 gap factor ~3000.
- The N_REAL=30 inter-seed spread reported in Session 3 (INFO-016) is not
  present at N >= 200.

### Part 2 — Three independent probes (this session)

Three probe scripts written and run this session. All use the exact
extract_v1 protocol. Each produces a different signal. None has a
verdict attached. Files: `probe_system_directions.py`,
`probe_highseed_convergence.py`, `probe_antisymmetric_diagnostic.py`.
JSONs alongside.

**Probe 1 — system directions (8 systems, 3 seeds, N_REAL=500):**

System                | mean v3                       | cos (-,-,+) | cos (+,+,+) | inter-seed cos
----------------------|-------------------------------|-------------|-------------|---------------
OU (gamma=1, sigma=1) | (+0.412, +0.373, +0.831)      | +0.36       | +0.999      | +0.994
uniform[-1,+1]        | (+0.218, +0.304, +0.928)      | +0.55       | +0.970      | +0.997
gauss sigma=0.1       | (+0.443, +0.307, +0.842)      | +0.38       | +0.994      | +0.998
gauss sigma=1.0       | (+0.469, +0.325, +0.821)      | +0.35       | +0.995      | +0.999
student-t df=3        | (+0.350, +0.426, +0.834)      | +0.36       | +0.998      | +0.998
laplace               | (+0.383, +0.406, +0.830)      | +0.36       | +1.000      | +0.995
sine + small noise    | (+0.485, +0.260, +0.835)      | +0.38       | +0.986      | +0.998
logistic r=3.9        | (+0.167, +0.282, +0.945)      | +0.59       | +0.955      | +0.976

All eight had cos > 0.95 to (+1,+1,+2)/sqrt(6) and cos < 0.59 to
(-1,-1,+2)/sqrt(6). Inter-seed cos > 0.97 for all eight.

Full 6D null vectors stored in `probe_system_directions_results.json`.

This may speak to INFO-014 (the 8-system attractor claim). It also may
not — Session 3's `cross_system_test.py` (the canonical INFO-014 source)
is not in this repo. Side-by-side protocol verification would be a
separate probe.

**Probe 2 — high-seed convergence (30 seeds, N_REAL=500):**

Case        | mean v3                       | per-coord SD                | cos (+,+,+) | cos (-,-,+) | inter-seed cos mean/min
------------|-------------------------------|----------------------------|-------------|-------------|------------------------
BASELINE v=1| (+0.424, +0.397, +0.814)      | (0.058, 0.050, 0.016)      | +0.9998     | +0.329      | +0.994 / +0.930
GAUSS_V0    | (+0.390, +0.352, +0.851)      | (0.066, 0.073, 0.015)      | +0.998      | +0.392      | +0.990 / +0.941

The H_a^2 and H_b^2 coordinates have ~6-7% SD across seeds. The H_a*H_b
coordinate has ~1.5% SD. Both cases still hit cos +0.998 to
(+1,+1,+2)/sqrt(6) at 30 seeds.

**Probe 3 — INFO-017 antisymmetric energy fraction:**

Case                              | anti_frac (mean of 3 seeds)
----------------------------------|-----------------------------
linear_drift_baseline_v1          | 0.495
linear_drift_gauss_v0             | 0.501
linear_drift_asym_v1v2            | 0.493
linear_drift_rho_0p3              | 0.364
linear_drift_rho_0p7              | 0.172
linear_drift_rho_0p95             | 0.032
OU_gamma1_sigma1                  | 0.501
uniform_pm1                       | 0.495
gauss_small_0p1                   | 0.501
gauss_large_1p0                   | 0.501
student_t_df3                     | 0.500
laplace_scale1                    | 0.490
sine_w0p5_phi0p7                  | 0.536
logistic_r3p9                     | 0.507

Everything except the coupled linear-drift cases sits in 0.49-0.54.
The coupling sweep is the only thing that moves anti_frac, monotonically
down from ~0.50 to ~0.03 as rho goes from 0 to 0.95.

### Part 3 — earlier suite at N_REAL=500 (reference)

Before the three probes, a six-case suite was run to confirm the protocol
reconstruction matched Session 4 exactly. File: `falsifier_suite_INFO018.py`
(name kept — already committed under that name; new work uses "probe").
Results in `falsifier_suite_INFO018_results.json`.

Cases: BASELINE (v=v_b=1), v=v_b=0 (no drift), v_a=1/v_b=2 (asymmetric),
rho=0.30, 0.70, 0.95 (coupling sweep). All N_REAL=500, 3 seeds.

BASELINE matches Session 4 to 4 decimals (mean v3 +0.466,+0.370,+0.804;
cos to (+1,+1,+2)/sqrt(6) = +0.9975). Protocol reconstruction is right.

## Files produced this session (all in repo root)

From Greg (Session 4 baseline):
- `extract_v1_linear_drift.py`
- `run_NREAL_sweep.py`
- `analyze_spectrum.py`
- `L1_linear_drift_NREAL_sweep.json`

From this session:
- `falsifier_suite_INFO018.py`
- `falsifier_suite_INFO018_results.json`
- `probe_system_directions.py`
- `probe_system_directions_results.json`
- `probe_highseed_convergence.py`
- `probe_highseed_convergence_results.json`
- `probe_antisymmetric_diagnostic.py`
- `probe_antisymmetric_diagnostic_results.json`

All to be mirrored to `E:\information_layer\` and
`F:\Factory\knowledge\information_layer\` per the Operating Rules.

## Branch state

- Local branch: `claude/linear-drift-nreal-sweep-cjRkM`
- Remote: `origin/claude/linear-drift-nreal-sweep-cjRkM`
- Tip commit: `8a28479` Session 4 probes + 3 new Operating Rules in CLAUDE.md
- Pushed: yes. No PR opened. `main` untouched.

## Ledger updates (data-level only, per Result Discipline)

**INFO-016 — DATA UPDATED**: linear drift extraction inter-seed cos rises
monotonically with N_REAL: 0.79 (N=50), 0.85 (N=100), 0.985 (N=200),
0.99 (N=500). The N=30 spread reported in Session 3 is not present at
N >= 200. Open question from Session 3 about "structural blindness vs
sample noise" — the data point is now on the table; reading happens
with Greg in the next session.

**INFO-018 — ISOLATED FINDING (data; reading open)**: at N_REAL=500
the linear drift extraction in 3D subspace lands at mean v3
(+0.466, +0.370, +0.804), cos +0.9975 to (+1,+1,+2)/sqrt(6) and cos +0.31
to (-1,-1,+2)/sqrt(6). 6D null vector dominated by H_a, H_b components
(~90% of energy); the 3D subspace projection is the ~10% residual. MI
component ~10^-7.

**INFO-019 — DATA, READING OPEN (Session 4)**: at N_REAL=500 with three
seeds, eight diverse systems (OU, uniform, gauss x2, student-t, laplace,
sine, logistic) all land at cos > 0.95 to (+1,+1,+2)/sqrt(6) and cos
< 0.59 to (-1,-1,+2)/sqrt(6) in the 3D subspace projection of this
protocol. Inter-seed cos > 0.97 for all eight. The Session 3 INFO-014
claim that eight (different) systems land at cos >= 0.99 to
(-1,-1,+2)/sqrt(6) is not reproduced in this measurement — but the
canonical Session 3 `cross_system_test.py` was not run here, so the two
results cannot be directly compared without protocol alignment.

**INFO-020 — DATA, READING OPEN (Session 4)**: at N_REAL=500 with 30
seeds (10x the original three), the linear drift BASELINE mean v3 is
(+0.424, +0.397, +0.814) (cos +0.9998 to (+1,+1,+2)/sqrt(6)); the v=0
Gaussian case lands at (+0.390, +0.352, +0.851) (cos +0.998). Per-seed
SDs are (0.058, 0.050, 0.016) and (0.066, 0.073, 0.015) — H_a*H_b
coordinate is ~4x tighter than the H_a^2 and H_b^2 coordinates.

**INFO-021 — DATA, READING OPEN (Session 4)**: INFO-017 antisymmetric
energy fraction measured across 14 cases. All eight system-direction
probes and three of the six suite cases (baseline / v=0 / asymmetric)
sit at anti_frac in [0.49, 0.54]. The three coupled linear-drift cases
move anti_frac monotonically down: 0.364 (rho=0.30), 0.172 (rho=0.70),
0.032 (rho=0.95). Coupling is the only variable in this scan that moves
anti_frac off the ~0.50 region.

No tag promotions in this session. All four new data entries are
isolated findings pending look-through with Greg.

## What's queued (from Session 3 plus this session)

From the Session 3 queue, with what we did this session:

1. ~~N_REAL sweep on linear drift~~ DONE — Greg ran it; data above.
2. Multi-seed damped oscillator mapping. NOT YET RUN.
3. Cluster the misses by structure. NOT YET RUN.
4. Boundary mapping. NOT YET RUN.
5. Stress-test the attractor (base-of-structure heuristic test). NOT YET RUN.
6. Domain-native operator bases (Geo 3D, Bio multi-variable, Chem multi-species).
   NOT YET RUN.

Open from this session (no priority assigned — Greg's call):

- Side-by-side protocol verification against Session 3's canonical
  `cross_system_test.py` (the INFO-014 source). Not in this repo. Would
  need Greg to upload or point to it. The INFO-019 / INFO-014 tension
  cannot be read until this is done.
- Look-through of the Session 4 probe data with Greg, against the
  candidate-interpretation register, with the deflationary reading
  always present per Result Discipline.

## Prompt for next session

"Resuming Information Layer work. Read CLAUDE.md first (it has the three
new Operating Rules from Session 4, plus the Session 3 ledger). Then read
SESSION_HANDOFF_2026-05-26_v4.md. Open question on the table: look at
INFO-018, 019, 020, 021 raw data together. No verdict in advance. No
verdict on first look. Use the speaking posture from Rule C. We are
the first ones here — follow every road, look under every stone. No
safe assumptions."
