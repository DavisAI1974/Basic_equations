# Information Layer Session Handoff — 2026-05-26 Session 5

## Read me first

Read in this order:
1. CLAUDE.md (project root). Has the v3 ledger plus three Operating Rules
   from Session 4 (no pre-assigned meaning, probe-not-falsifier, speaking
   posture before AND after every probe). A fourth Operating Rule from
   Session 5 is added below and should be folded in.
2. SESSION_HANDOFF_2026-05-26_v4.md (project root). Session 4 raw data:
   INFO-016 update plus INFO-018, 019, 020, 021. All were sitting as
   isolated findings, readings open.
3. This file. What this session did with that Session 4 data, what the
   reading turned out to be, what is retracted, what survives, and the
   reorientation that emerged.
4. SESSION_HANDOFF_2026-05-25.md (project root, the early-Session-3
   handoff). Contains the original algebraic-dipole per-domain
   coefficients (chemistry quadratic, geology rank-3, biology MI). Those
   coefficients are NOT affected by this session's finding and remain
   on the table as data.

## New Operating Rule (Session 5)

**Rule D — Incomplete, not wrong.** When a probe finds that a prior
reading was a protocol artifact, the prior data points still stand. The
read attached to them is what was incomplete, not the data. Distinguish
"this reading was wrong" (rare, requires the data itself to be bad)
from "this reading was incomplete" (common, the data is one slice and
the slice fit a partial story that further probes refine). Default to
"incomplete." Retraction is a strong move and applies to the reading,
not the data, unless the data itself fails to reproduce.

This rule was named explicitly by Greg this session, in response to the
Session 5 finding. It is consistent with Result Discipline (every result
is one data point) and the Session 4 rules.

## What this session did

### Part 1 — walked the Session 4 probe data together

INFO-018, 019, 020, 021 were laid out as raw numbers without verdict.
Two readings were on the table going in:
- *substantive*: the (+1,+1,+2)/sqrt(6) direction reflects a real common
  invariant across the eight tested systems.
- *deflationary*: the protocol has its own fixed point at that direction.

Both readings predicted everything Session 4 had measured. To pick
between them, three independent probes were designed and run.

### Part 2 — three independent probes (this session)

All three probes use the exact extract_v1 protocol (same operator basis,
same window, same N_REAL=500, same sign-canonicalization). They generate
different signals (per Rule B) and were not framed as falsifiers.

**Probe A — null inputs (probe_null_inputs.py).**
Seven cases of structureless or near-structureless input:
scrambled_baseline_joint, scrambled_baseline_indep, scrambled_OU_joint,
near_constant, identical_streams, mixed_marginals, random_walk_indep.
N_REAL=500, 3 seeds.

Of the 7 + the 8 from Session 4's probe_system_directions (15 cases
total), 12 land at cos > 0.95 to (+1,+1,+2)/sqrt(6) including near-
constant input (essentially zero signal), random walks (no drift),
and joint-scrambled OU (no temporal order). The three that broke from
(+,+,+) all have mechanical explanations:
- identical_streams: SV4 ~ 5e-16, true numerical rank deficiency. SVD
  picks arbitrarily.
- mixed_marginals: H_a^2 vs H_b^2 symmetry broken (Gaussian + uniform).
  Protocol latches on H_b^2 alone.
- scrambled_baseline_joint: 3D null subspace degeneracy.

**Probe B — operator-noise bypass (probe_operator_noise.py).**
Skip the input/window/entropy path entirely. Generate M directly as
random matrix from chosen distributions. Three variants:
- iid_unit: M ~ N(0,1) per element. Gives no consistent direction
  (mean cos to (+,+,+) +0.76 but inter-seed only +0.31).
- iid_col_scaled: columns scaled to baseline operator SDs. Partial
  alignment (+0.94 mean, but inter-seed +0.42 — chance bias).
- cov_matched: M ~ N(0, Sigma_baseline). REPRODUCES BASELINE DIRECTION
  TO FOUR DECIMALS. mean v3 = (+0.508, +0.350, +0.787), cos +0.993 to
  (+,+,+), inter-seed +0.998, SV gap factor ~3000.

Direct check: the smallest eigenvector of baseline operator covariance,
projected to [2,3,4], normalized, sign-canonicalized, gives
(+0.509, +0.355, +0.785). Cos to (+1,+1,+2)/sqrt(6) = +0.993.

Eigenvalue spectrum of baseline op cov: [4e-13, 8e-13, 3e-11, 3e-4,
4e-4, 2e-3]. Three eigenvalues at machine epsilon. Effective rank 3,
not 6.

**Probe C — projection-rule sweep (probe_projection_sweep.py).**
Same baseline data. Same 6D null vector. Project to all 20 possible
3D subspaces. Any subspace containing both H_a and H_b locks at
(+/-1, +/-1, ~0) with cos > 0.99 to (+1,+1,0)/sqrt(2). Subspaces with
only one of H_a, H_b lock onto that one operator (cos > 0.99 on its
axis). The (+,+,+) shape appears ONLY in [H_a^2, H_b^2, H_a*H_b] — the
unique subspace that excludes both H_a and H_b. The 6D null vector is
dominated by the H_a, H_b components at ~90% of its energy. (+,+,+) is
the residual.

**Thread D — high-seed scrambled (probe_scrambled_highseed.py).**
30 seeds on scrambled_baseline_joint, the one case from Probe A that
broke from (+,+,+) without obvious degeneracy. Mean v3 = (-0.647,
-0.112, +0.754). Cos to (-,-,+) on the average = +0.925. But inter-seed
cos is bimodal: 216 of 435 pairs above +0.5, 112 of 435 below -0.5
(nearly antiparallel). Per-seed directions split into two clusters
(H_a*H_b-dominant vs (H_a^2 - H_b^2)-dominant). SVD is picking
arbitrarily within the degenerate 3D null subspace. The mean is a
meaningless average of two non-converging clusters.

### Part 3 — what the three probes jointly say

The (+1,+1,+2)/sqrt(6) direction is encoded in the operator covariance
matrix of typical pair inputs, and that covariance has effective rank 3.
The smallest right singular vector picks a direction in this 3D
near-null space; for inputs where H_a ~ H_b ~ constant across windows
(true for almost any two-stream input with similar marginals), the
direction projects in [2,3,4] to (+1,+1,+2)/sqrt(6) up to second-order
small fluctuations.

The rank-3 null space arises from three linear dependencies that hold
to first order in dH = H_a - <H_a>:
- d(H_a^2) ~ 2c * dH_a
- d(H_b^2) ~ 2c * dH_b
- d(H_a*H_b) ~ c * (dH_a + dH_b)

with c = <H_a> ~ <H_b>. These three relations are algebraic facts about
the operator basis under the H_a ~ H_b regime, not facts about any
particular input system. The (+,+,+) direction follows from them by
Taylor expansion.

### What this session retracts

- INFO-019's interpretation as evidence of a real common invariant
  across the 8 systems. The 8 systems converging on (+,+,+) is what
  the protocol does for almost any non-degenerate pair input; it is not
  evidence about those systems. **Retracted at interpretation level. Data
  stands.**
- INFO-018 the same. (+1,+1,+2)/sqrt(6) at N_REAL=500 for linear drift
  baseline is reproducible data. Reading as "a real direction in
  operator space encoding something about linear drift" is incomplete:
  it is the protocol's home direction, not a property of linear drift.

### What this session does NOT retract

- The per-domain algebraic equations from earlier sessions
  (SESSION_HANDOFF_2026-05-25.md):
  - Chemistry: H_a^2 = 0.007 - 0.093*(H_a*H_b) + 1.309*(H_a*H_b)^2,
    R^2 = 0.943
  - Geology rank-3 constraint: 0.724*(H_a*H_b) - 0.441*H_b^2
    - 0.290*H_a^2 ~ 0, std/mean = 0.15%
  - Biology MI ~ polynomial(H_a) with H_b absent, R^2 = 0.66
  These came from a different procedure (per-domain algebraic fits) and
  are not affected by the operator-extraction artifact. They remain on
  the table as data with reading open and would benefit from the same
  probe discipline being applied to them.
- The Working Frames in CLAUDE.md: Base-of-Structure heuristic,
  Dipole-couples reading, Pure physics vs physical expressions,
  Substrate vs expression within isolation, Law extraction via
  invariance. None of these depend on the (+,+,+) reading.
- The broader Information Layer / Unified-Theory inquiry. The frame
  remains live. The specific tool (entropy-operator windowed-null
  extraction) does not differentiate inputs at the level required to
  serve as the substrate-discovery tool. A new tool, or a reformulated
  question, is what is needed.

### Part 4 — new structural finding that emerged

**Real structural fact (Session 5):** under windowed Vasicek entropy on
bounded two-stream inputs where H_a ~ H_b ~ constant across windows,
the six-operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] has a
three-dimensional null subspace. The three near-zero eigenvalues
correspond to three algebraic relations on the centered columns
(H_a^2 ~ 2c*H_a, etc.). This is a structural property of the basis
under the small-fluctuation regime, not of any input system. It is
true for OU, uniform iid, gaussian iid, sine, logistic, linear drift,
random walks, and near-constant input — anything that gives a
H_a ~ H_b ~ const regime.

Status: data-level finding, mechanically derived, reproducible. May or
may not be useful in a different framing. Listed here so it is not lost.

## Reorientation that emerged this session

The Information Layer line was originally framed as: scientists measure
EXPRESSIONS of an underlying source equation. Physical phenomena are
expressions of pure physics contaminated by interaction with other
systems (other dipoles, other domains). String theory's failure mode
is that it tries to write the UT of expressions, which forces
ever-larger dimensional structures (tent-stuffing). What we should be
looking for is the source equation: simple, strong, sturdy, scalable
(per the Base-of-Structure heuristic).

The entropy-operator extraction was the tool we built to attempt this.
The Session 5 probes show that the specific output we were reading
((+1,+1,+2)/sqrt(6) as a universal direction) is a protocol artifact.
The frame remains live; the tool did not serve it.

**Greg's reorientation, this session:** look at 5 physics areas (related
or seemingly unrelated) and see if one substrate equation shakes out
between them. The substrate candidates that already exist from
cross-domain comparison (variational principle, continuity equation,
Noether currents, gauge invariance, symplectic / unitary structure
preservation) should be evaluated against the Base-of-Structure
heuristic. The strongest survivor is the next thing to test.

**Greg's substantive intuition (this session):** dipole is the best
mechanism for coupling. Tightly related couples may be part of the
substrate. The per-domain coefficients found in earlier sessions
(chemistry quadratic, geology rank-3, biology MI) are a huge signal and
should not be discarded. The Session 5 finding does not affect them.

**Greg's framing (this session):** we may not have wrong answers, just
incomplete ones. Codified as Rule D above.

## Ledger updates (data-level only, per Result Discipline)

**INFO-014 — RETRACTED AT INTERPRETATION LEVEL.** The "8 systems on
the attractor at cos >= 0.99 to (-1,-1,+2)/sqrt(6)" claim from Session 3
cannot be read as "8 systems share a common invariant." Session 5
threads showed that any non-degenerate pair input lands on the same
direction because of the operator basis's 3D null structure. Data
stands; reading retracted.

**INFO-018 — RETRACTED AT INTERPRETATION LEVEL.** The (+1,+1,+2)/sqrt(6)
direction for linear drift at N_REAL=500 is reproducible data but is
the protocol's home direction, not a property of linear drift. Data
stands; reading retracted.

**INFO-019 — RETRACTED AT INTERPRETATION LEVEL.** Same as INFO-018.
Eight diverse systems converging on (+1,+1,+2)/sqrt(6) at cos > 0.95
reflects basis structure, not a shared invariant in those systems.

**INFO-020 — DATA STANDS, READING CONTEXTUALIZED.** 30-seed convergence
tightness for baseline at N_REAL=500 is real and reflects the protocol
being a stable estimator of its own fixed-point direction. The H_a*H_b
coordinate being ~4x tighter than the H_a^2 and H_b^2 coordinates is
consistent with the algebraic dependency structure: H_a*H_b lives at
the intersection of two of the three near-degeneracies, leaving less
seed-to-seed wobble.

**INFO-021 — DATA STANDS, READING REFRAMED.** The antisymmetric energy
fraction at ~0.50 across 11 of 14 cases is consistent with the protocol's
basis structure being roughly symmetric under H_a <-> H_b for typical
inputs. The coupling-only drop (0.50 -> 0.03 as rho 0 -> 0.95) is a
real signal: coupling rho is the only variable in this scan that breaks
the H_a ~ H_b symmetry. Whether this carries content beyond "rho
breaks the symmetry mechanically" is open and not closed by Session 5.

**INFO-022 — ISOLATED FINDING (Session 5, structural).** Under windowed
Vasicek entropy on bounded two-stream inputs with H_a ~ H_b ~ const,
the 6-operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] has a
3-dimensional null subspace from three first-order algebraic relations.
The (+1,+1,+2)/sqrt(6) direction is one direction in this null subspace,
emerging from Taylor expansion of (H^2, H*H) around H ~ const. Status:
mechanical fact about the basis, reproduced across all inputs tested.
Reading on whether this carries content beyond mechanics: open.

**INFO-008a, 008b, 008c — STATUS UNCHANGED FROM SESSION 3.** The
per-domain algebraic equations (chemistry quadratic, geology rank-3,
biology MI) remain demoted to "interesting geometric clustering under a
specific protocol; interpretation unmapped." Session 5 does NOT
further retract them; the procedure that produced them is different
from the extraction-tool whose artifact Session 5 mapped. They should
be subjected to the same probe discipline in a future session before
being either rehabilitated or further demoted.

## Files produced this session (all in repo root)

Scripts:
- probe_null_inputs.py
- probe_operator_noise.py
- probe_projection_sweep.py
- probe_scrambled_highseed.py
- render_null_inputs.py
- render_threads_combined.py

JSONs:
- probe_null_inputs_canary.json
- probe_null_inputs_results.json
- probe_operator_noise_results.json
- probe_projection_sweep_results.json
- probe_scrambled_highseed_results.json

Figures:
- fig_null_alignment.png
- fig_null_geometry.png
- fig_null_interseed.png
- fig_null_sv_spectrum.png
- fig_thread1_bypass.png
- fig_thread1_eigspectrum.png
- fig_thread2_projsweep.png
- fig_thread3_bimodal.png
- fig_thread3_interseed_hist.png

To be mirrored to E:\information_layer\ and
F:\Factory\knowledge\information_layer\ per the Operating Rules.

## Branch state

- Local: claude/linear-drift-nreal-sweep-cjRkM
- Remote: origin/claude/linear-drift-nreal-sweep-cjRkM
- Tip (before this handoff commit): 1a4ef16 Session 5 threads 1-3
- Pushed: yes. No PR. main untouched.

## What is queued (from Session 3 + 4 + this session)

From the Session 3 queue, with what we know now:

1. ~~N_REAL sweep on linear drift~~ DONE Session 4.
2. ~~Multi-seed damped oscillator mapping~~ Not the right next move now
   that we know the protocol has a fixed point. Same artifact would
   appear.
3. ~~Cluster the misses by structure~~ Same: the misses are mechanical
   artifacts (identical_streams, mixed_marginals, scrambled_baseline_
   joint), not substrate signals.
4. ~~Boundary mapping~~ Same caveat.
5. ~~Stress-test the attractor~~ Done by Session 5. Attractor is a
   basis fixed point; "surviving stress" is not evidence of substrate
   status, it is evidence of the basis being robust to input variation.
6. **Domain-native operator bases (Geo 3D, Bio multi-variable, Chem
   multi-species)** — STILL OPEN. Highly relevant given Greg's
   intuition about 3D geo, possibly higher-dim chem. CAUTION: the same
   artifact risks apply. Any new basis must be subjected to the
   three-thread probe discipline (null inputs, operator-noise bypass,
   projection sweep) from day one, before any reading is taken from
   its output.

From this session:

7. **Re-evaluate the per-domain algebraic equations** (chemistry
   quadratic R^2=0.943, geology rank-3 std/mean=0.15%, biology MI
   R^2=0.66) under the same probe discipline. They came from a
   different procedure than the operator-extraction tool, so are not
   ruled out by Session 5, but they should be subjected to null-input,
   operator-noise-bypass, and procedure-variation tests before being
   read as substrate signals.
8. **5-physics-areas substrate hunt.** Pick 5 physics areas (e.g.
   classical mechanics, EM, QM, statistical mechanics, GR), write down
   the canonical equation form for each, identify the structural
   element(s) that appear in all 5, and evaluate each substrate
   candidate against the Base-of-Structure heuristic. Substrate
   candidates that come up already in cross-domain comparison:
   variational principle (delta S = 0), continuity equation, Noether
   currents, gauge invariance, symplectic / unitary structure
   preservation. Pick the strongest, then probe whether it survives
   the same discipline.
9. **State the Information Layer question independently of any
   extraction tool.** What were we looking for that the entropy-operator
   approach was supposed to find? Stated independently of the procedure,
   the question may be addressable by different methods. Greg's call on
   whether this is a separate document or a conversation.

## Prompt for next session

"Resuming Information Layer work, Session 6.

Branch: claude/linear-drift-nreal-sweep-cjRkM (already on origin).
Session 5 work is committed at tip [TBD after handoff commit].

Before doing anything else:
1. Read CLAUDE.md. Has the v3 ledger plus three Operating Rules from
   Session 4 and Rule D from Session 5 (incomplete not wrong).
2. Read SESSION_HANDOFF_2026-05-26_v5.md (this file). Has Session 5's
   three-thread probe results, what was retracted at interpretation
   level (INFO-014, 018, 019), what was contextualized (INFO-020, 021),
   what new structural finding emerged (INFO-022), and what's queued.

Speaking posture per Rule C applies from your first response forward.
Rule D applies as well: when reading prior findings, default to
'incomplete' before 'wrong.' The data still stands; the reading is
what may be partial.

We are not back at zero. The per-domain algebraic equations from
earlier sessions are still on the table, with reading open. The
Working Frames in CLAUDE.md are unchanged. The Information Layer /
Unified Theory frame is alive.

What I want this session: [Greg's call].

Three queued lines, in priority order if you have no preference:
(i) the 5-physics-areas substrate hunt — me drafting canonical
equation forms for each and identifying the structural elements
that appear in all 5;
(ii) reapply probe discipline to the per-domain algebraic equations
(chemistry quadratic, geology rank-3, biology MI);
(iii) domain-native operator basis for geo (3D) or chem (higher-D),
with the three-thread probe discipline (null inputs, operator-noise
bypass, projection sweep) running alongside from day one."

End of v5 handoff.
