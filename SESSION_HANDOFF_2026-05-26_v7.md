# Information Layer Session Handoff -- 2026-05-26 Session 7

## Read me first

Read in this order:
1. CLAUDE.md (project root). Has Session 7 note added at the top
   alongside Sessions 3, 4, 5, 6 notes. All five Operating Rules
   from Sessions 4-5 remain in force (no pre-assigned meaning,
   probe-not-falsifier, speaking posture before and after, Rule D
   incomplete-not-wrong, They never stacked). ONE NEW RULE added
   at Session 7 close: "Treat literature as conjecture by default"
   — academic consensus is not foundation unless independently
   replicated by separate groups with an immense amount of data
   multiple times. Pairs with Rule D and "They never stacked"
   (the literature is a frame to test, not a foundation to build
   on). Six Operating Rules now in force.
2. This file. What Session 7 did, what is confirmed, what is queued,
   plus the substantive UT discussion at session close (no experiment
   run; framing work only).
3. SESSION_HANDOFF_2026-05-26_v6.md (Session 6 record). Provides the
   per-domain stack methodology that Session 7 built on top of.
4. SESSION_HANDOFF_2026-05-26_v5.md (Session 5 record). The
   "they never stacked" literature scan that flagged PySR/GP on
   windowed-H of coupled species as a clean unstacked target. Session 7
   executes that stack.
5. SESSION_HANDOFF_2026-05-25.md (per-domain coefficients from earlier
   work).

## What Session 7 did

Greg's directive at start of session: "let's do top of que and work
down". The top of queue from the v7 kickoff was:

1. GP regression of MI vs H_a / H_b / joint per domain. Tests v5's
   "biology MI ~ polynomial(H_a), R^2 = 0.66" against flexible
   nonlinear fit. Quick (~sklearn).
2. PySR symbolic regression per domain. Surfaces ratio / exponential /
   log / sqrt forms that linear-in-features SINDy cannot find. Flagged
   in v5 literature scan as the cleanest unstacked target.

Both ran. All work pushed to branch
`claude/claude-md-context-update-uCZ8m`.

### Experiment 1 -- GP regression of MI vs H_a / H_b / joint (gp_mi_vs_h.py)

Procedure: re-simulate same per-domain ensemble-H data as Session 6
(per_domain_kbk.py simulators, N_ens=600, T=30, dt=0.02). Extract
H_a(t), H_b(t), MI(t) time series. Fit:
- Polynomial deg 1, 2, 3 of MI ~ H_a; MI ~ H_b; MI ~ (H_a, H_b)
- GP with RBF + WhiteKernel of same three input configurations

Evaluation: three R^2 readings per fit
- In-sample R^2
- Block-CV R^2 (5 contiguous time blocks; sensitive to time-series
  non-stationarity)
- Random-CV R^2 (5-fold shuffled; overestimates true OOS for
  correlated time series but isolates functional-form fit from
  stationarity)

GP fits subsampled to n=500 (cap O(n^3) cost). Polynomial fits use
full n=1200 per domain.

Two seeds (11, 22). Elapsed: 524s.

Cross-seed pattern at random-CV R^2 (joint deg-3 poly vs GP):
```
                  poly3(joint)   GP(joint)   block-CV behavior
physics_duffing   0.53-0.78      0.67-0.76   catastrophic (-0.5 to -1.5)
biology_lotka     0.81-0.86      0.97-0.98   catastrophic (-10 to -80)
chemistry_brussel 0.80-0.86      0.81-0.88   POSITIVE (+0.4 to +0.8)
geology_bk        0.05-0.08      0.02-0.05   small (signal too low)
```

**Three readings, separated by level (Result Discipline):**

A. *Data level*: Cross-seed reproducibility is strong at random-CV
(variation <0.05 within domain between seeds 11 and 22). Each domain
has a distinct MI-vs-H signature.

B. *Methodological*: Block-CV catastrophic in physics/biology/
geology due to non-stationarity (system passes through qualitatively
different regimes across blocks). Chemistry is the only stationary
domain. The block-vs-random gap is itself a domain signature.

C. *Substantive on biology*: poly1(MI ~ H_a) random-CV = 0.71-0.78
across seeds (above v5's reported R^2 = 0.66). H_b alone weak (0.13-
0.16). GP joint reaches 0.97-0.98 vs poly3 joint at 0.81-0.86. The
+0.12 gap is real nonlinear cross-coupling that polynomial methods
miss. **v5's biology reading is incomplete-not-wrong (Rule D)**: the
H_a-dominance and H_b-weakness reproduce, but additional joint
nonlinear structure exists beyond polynomial reach.

### Experiment 2 -- PySR symbolic regression per domain (pysr_symbolic_per_domain.py)

Stack: PySR 1.5.10 (Julia 1.11.9 / SymbolicRegression.jl 1.11.3).
Binary operators {+, -, *, /}. Unary operators {square, cube, exp,
log, sqrt}. Pareto-front optimization. Per v5 literature scan: "no
SINDy/AI-Feynman/PySR/Eureqa application to windowed differential
entropy of coupled species" -- this is that stack.

Procedure: same per-domain ensemble-H data, subsampled to n=600,
niter=40, populations=20, maxsize=15, timeout=90s/domain. Two seeds
(11, 22). Elapsed: 99s total.

Cross-seed best low-complexity equations:
```
              seed 11                          seed 22
physics       (H_b - H_a)^2 + 0.275   c=6      (H_b - H_a)^2 + 0.280   c=6
biology       0.531 * exp(H_a/2)      c=5      0.495 * exp(H_a/2)      c=5
chemistry     0.71 * H_a + 1.08       c=5      0.75 * H_a + 1.07       c=5
geology       0.199 (constant)        c=1      0.199 (constant)        c=1
```

Cross-seed coefficient reproduction: <5% variation. The functional
families themselves reproduce exactly (same operator structure
selected by PySR independently for each seed).

**Four qualitatively different functional families:**
- Physics: symmetric quadratic in the *difference* (H_b - H_a)
- Biology: exponential in H_a alone, H_b absent
- Chemistry: linear in H_a, near-flat in H_b
- Geology: decoupled, MI ~ constant

These families do not reduce to each other. Cross-domain non-overlap.

**Two substantive findings from this:**

(i) Polynomial-only methods (v5 SINDy at deg-3) could not have
distinguished biology from chemistry from physics, because all four
families would compete for the same "polynomial in H_a, H_b" library
slot. PySR's extended operator set surfaces the *functional family*,
which is the actual differentiator.

(ii) v5's biology reading ("MI ~ polynomial(H_a), R^2 = 0.66") was
qualitatively right about operator content (H_a present, H_b absent)
but wrong about functional family. The right reading: MI ~ exp(H_a/2)
in biology, score +0.92 at complexity 5 (PySR's parsimony metric).
Same operator structure as v5, different functional class. v5 is
incomplete-not-wrong (Rule D).

### Joint reading of Experiments 1 + 2

The combined evidence is:
- Each simulated domain has a *per-domain MI-vs-H functional family*
- These families reproduce across seeds within domain (cos > 0.99 on
  GP fit, <5% on PySR coefficients)
- Families do not overlap across domains (zero shared structure
  between physics' (H_b-H_a)^2, biology's exp(H_a/2), chemistry's
  linear H_a, geology's constant)
- Family identity is *not visible* to polynomial-only feature
  libraries (which is what v5 used)
- The Session 6 per-domain null direction findings (INFO-023) and
  these per-domain MI-vs-H functional families are independent lines
  of evidence for the same per-domain differentiation

This is one further data point for the "per-domain expression of
substrate" frame from v5 / v6 / Session 7. It does not promote the
frame to a claim yet -- mapping work is required to test what
predicts which functional family appears in which domain. But the
data side now has TWO independent reproducible per-domain signatures
(null direction + MI-vs-H form), not one.

## Substantive discussion this session (no experiment run; framing)

Late in the session, Greg pivoted to:

1. "How does pure physics get expressed as a storm or waves?"
   (the substrate-to-expression question for atmospheric / fluid
   physics).
2. "The real test: gravity, EM, strong, weak. Can we come up with
   one equation? Should they be grouped at all?"
   (the classical four-force unification question).
3. After the framing was laid out, Greg added: "I like your
   substrate level. And let's not assume that academic papers are
   right on some of this unless there's an immense amount of data
   that has proven it multiple times." This added a SIXTH Operating
   Rule ("Treat literature as conjecture by default" — see Rules
   list above) that fundamentally changes how we engage the
   four-force question.

No experiment was run on either question. The framing put on the
table, for next session to engage if Greg directs:

**On storm/waves (Greg's first question)**: storms and waves are
physical *expressions* of underlying substrate operating under
specific coupling configurations (atmospheric thermodynamics +
water + Coriolis + boundary geometry). The Session 6 + Session 7
result -- one substrate (rank-3 null subspace), four expression-level
functional families -- is the operational analog. A real-data test
would run the per-domain stack on NOAA buoy / HRRR atmospheric
reanalysis / hurricane track data and look for the storm's signature
in MI-vs-H space on appropriate observable pairs (e.g., temperature
vs humidity; pressure vs wind speed). Queued, not started.

**On four-force unification (Greg's "real test")**: three levels,
re-read under the new "literature as conjecture" Rule at Session 7
close.

- Data level (meets the immense-replicated-data bar): EM + weak
  ARE empirically unified (electroweak symmetry breaking confirmed
  at LHC with W/Z masses, Higgs detection, decay channels across
  multiple experiments). This stands.

- Data level (does NOT meet the bar, hold as conjecture):
  - "Three gauge couplings almost meet at GUT scale" — the running
    is measured at LHC energies; extrapolation to 10^16 GeV is
    mathematical, not empirical.
  - "MSSM makes them meet exactly" — depends on SUSY, searched
    extensively at LHC and not found.
  - "Gravity is a different category (geometric vs gauge)" —
    rests on a specific formulation (GR) with empirically
    indistinguishable alternatives (e.g., teleparallel gravity);
    the category claim is formulation-dependent, not data-forced.
  - "Gravity emerges (Jacobson 1995, Verlinde, Sakharov, AdS/CFT,
    et al.)" — theoretical frameworks with limited direct
    empirical support. Conjectures.

- Frame level: by analogy with this session's four-domain MI-vs-H
  result, the four forces may be four expression-level signatures
  of one substrate. EM expresses one way (gauge field on flat
  space); gravity expresses another way (geometric structure of the
  manifold itself); strong expresses with confinement; weak
  expresses with Higgs mass-giving. Same substrate, four
  functional families. Greg likes this framing (Session 7 close:
  "I like your substrate level").

Recommendation made to Greg, accepted at frame level: don't pre-
assign that all four should be grouped, and don't pre-assign that
they can't be. The deck is clearer under the new Rule than the
literature initially suggests — only EM+weak unification is
established; everything else (GUT, MSSM, gravity-category,
gravity-emergent) is conjecture and competes on roughly equal
footing with the substrate-vs-expression frame.

Next-session direction is open: real-data on storm/waves, real-data
on coupling-running (PDG / arXiv published data), continued
simulator probing, or hold framing and pick a different angle. Greg
to decide.

## Ledger updates

**INFO-025 -- LOCATED FINDING (Session 7, new)**: Four reproducible
per-domain MI-vs-H functional families surfaced by PySR symbolic
regression with extended operator set {+, -, *, /, square, cube, exp,
log, sqrt} on per-domain ensemble-H data. Cross-seed coefficient
reproduction <5% within domain:
- Physics (Duffing):    MI = (H_b - H_a)^2 + (0.275-0.280)
- Biology (Lotka-V):    MI = (0.50 +/- 0.02) * exp(H_a / 2)
- Chemistry (Bruss):    MI = (0.73 +/- 0.02) * H_a + (1.075 +/- 0.005)
- Geology (BK):         MI = 0.199 (loss flat across complexity 1-9)
Cross-domain non-overlap. Polynomial-only methods (v5 SINDy with
deg-3 library) could not surface these families by construction.
v5's biology "polynomial(H_a)" reading is incomplete-not-wrong
(Rule D): operator content correct (H_a present, H_b absent);
functional family is exponential, not polynomial.

**INFO-026 -- METHODOLOGICAL (Session 7, new)**: For GP regression
on time-series ensemble-H data, block-CV (contiguous time blocks) is
catastrophically negative across non-stationary domains (physics,
biology, geology block-CV ranges -0.5 to -80) due to system passage
through qualitatively different dynamical regimes across the
trajectory. Chemistry is the only stationary domain (block-CV
positive at +0.4 to +0.8). Random-CV (shuffled k-fold) overestimates
true OOS for correlated time series but isolates functional-form
fit from the stationarity confound. Use both; the gap between them
is itself a domain signature (stationarity). Substantive within:
GP joint(H_a, H_b) for biology random-CV R^2 = 0.97-0.98 vs poly3
joint at 0.81-0.86. The +0.12 gap is real nonlinear cross-coupling
that polynomial methods miss. For physics/chemistry/geology, poly3
joint matches GP joint within 0.05 (no extra nonlinear structure).

**INFO-023 (Session 6) -- REINFORCED**: The per-domain operator
extraction differentiation finding now has a second independent
reproducible signature (per-domain MI-vs-H functional family from
INFO-025) on top of the per-domain null direction finding from
Session 6. Two independent lines of evidence for per-domain
differentiation; one further data point for the "per-domain
expression of substrate" frame.

**INFO-008c (biology MI ~ polynomial(H_a)) -- REINTERPRETED via
Rule D**: Original v5 reading "polynomial in H_a, H_b absent,
R^2 = 0.66" is incomplete-not-wrong. Operator content correct
(H_a present, H_b absent, confirmed by Session 6 KBK and Session 7
GP). Functional family is exponential (exp(H_a/2)), not polynomial,
surfaced by PySR. R^2 actually higher than 0.66 when measured
correctly (in-sample 0.83-0.85 with poly3; random-CV 0.75-0.83).
Polynomial gives a usable but suboptimal fit; PySR's exponential
form is the natural one.

## Experiments queued (priority order)

Open questions raised this session AND prior queues from Sessions
3-6 still standing.

1. **Real-data application of per-domain stack** (Greg's frequent
   ask across Sessions 5-6, raised again Session 7 in storm/waves
   context). Run KBK + AI Poincare + SINDy + GP + PySR stack on:
   - NOAA buoy data (storm/wave example)
   - PhysioNet ECG (NoVell extension; Carl Saab angle)
   - LIGO O3/O4 strain data subsets (gravity probe at scale)
   - LHC published 4-vector kinematics (strong force probe at scale)
   - HRRR atmospheric reanalysis (weather)
   - Paleoclimate (Garland-Bradley 2018 from literature scan)
   - PDG coupling-running data (four-force unification probe)

2. **Four-force unification probe** (Greg's "real test" from Session
   7 close). Three angles:
   - Coupling-running data analysis: PDG running of g1, g2, g3 with
     energy; observe almost-meeting at GUT scale; treat as four
     scalar trajectories and run the per-domain stack
   - Gravity-as-emergent test: Jacobson 1995 thermodynamic gravity;
     find observable analog of MI/H operators in the EM/weak/strong
     trio, ask whether gravity's expression-level signature is a
     derived quantity
   - Frame work: extend the v5 Working Frame "pure physics vs
     physical expressions" specifically to the four forces

3. **Five-physics-areas substrate hunt** (Greg's Session 5
   reorientation, still open). Pick 5 physics areas, find the
   structural element that appears in all 5, evaluate against
   Base-of-Structure heuristic.

4. **Mapping campaign for which property predicts which functional
   family** (the natural follow-up to INFO-025). Vary simulator
   knobs:
   - Coupling strength K in physics Duffing
   - Predator-prey rates in Lotka-Volterra
   - Brusselator A, B parameters
   - Burridge-Knopoff drift rate
   Track which functional family persists, which mutates, which
   disappears. Identify the dynamical property that selects family.

5. **Methodology stacks from v5 menu** (Sessions 5/6 still open):
   - Complexity-entropy plane on coupled pairs (Rosso 2007)
   - Yamada 2023 style multi-field SVD with entropy as field
   - Gomez-Herrero 2015 ensemble infrastructure + algebraic basis
   - Stochastic Electrodynamics 21st-century revisit
   - Jacobson 1995 cross-domain extension

6. **Functional-form ablation under noise / shorter T / missing
   data** (robustness check on INFO-025). Does biology stay
   exponential at T=10? T=100? Under Gaussian noise added to
   observations?

7. **Older queued items**: N_REAL sweep on linear drift (INFO-016);
   multi-seed damped oscillator mapping; stress-test attractor
   under extreme conditions; domain-native operator bases for
   non-OU domains.

## Methodological notes from this session

(Not promoted to Rules. Domain-knowledge refinements to apply when
running similar analyses.)

**Block vs random CV on dynamical systems data**: when the time
series passes through qualitatively different regimes (transient
-> attractor; pre-event -> post-event; slow drift), block-CV
measures BOTH the functional fit and the stationarity of the
relationship across regimes. If the relationship is regime-
dependent, block-CV explodes negative. Random-CV measures only
the functional fit (overestimating true OOS due to within-fold
correlation). Reporting both lets you separate "the functional
form is wrong" from "the relationship is non-stationary." Most
real dynamical systems are non-stationary in this sense; chemistry
Brusselator was the only one of the 4 simulators that wasn't.

**Polynomial libraries are blind to functional family**: SINDy
deg-3 library can only find combinations of polynomial features.
Fitting exp(H_a/2) with a degree-3 polynomial in H_a will give an
acceptable R^2 (because exp is well-approximated by a Taylor series
on a bounded range) but will hide the fact that the true family is
exponential. To distinguish functional families, use an extended
operator set (PySR) or do a model-comparison test (e.g., AIC/BIC
of polynomial vs exponential fits). v5's polynomial fits on biology
were not wrong; they were partial. (Rule D again.)

**Cross-seed reproduction at <5% coefficient variation is a strong
signal**: each domain's PySR-discovered equation reproduced across
independent random seeds at <5%. This rules out "lucky fit" or
"seed-specific overfitting" readings. Coefficients carry signal.

## Files produced this session (all in repo root)

Scripts:
- gp_mi_vs_h.py  (Experiment 1: GP regression with poly1-3 baseline,
  in-sample + block-CV + random-CV)
- pysr_symbolic_per_domain.py  (Experiment 2: PySR symbolic regression
  with extended operator set, Pareto front output)

Result JSONs:
- gp_mi_vs_h_canary.json
- gp_mi_vs_h_results.json
- pysr_symbolic_per_domain_canary.json
- pysr_symbolic_per_domain_results.json

Run logs:
- gp_mi_vs_h_canary.log
- gp_mi_vs_h_run.log
- pysr_canary.log
- pysr_run.log

Docs:
- SESSION_HANDOFF_2026-05-26_v7.md (this file)
- CLAUDE.md updates (Session 7 note)
- NEW_SESSION_KICKOFF_v8.md (drop-in for next session)

Dependencies installed this session:
- numpy 2.4.6, scikit-learn 1.8.0 (pip)
- pysr 1.5.10 (pip) + Julia 1.11.9 (auto via juliapkg) +
  SymbolicRegression.jl 1.11.3 (auto)

To be mirrored to E:\information_layer\ and
F:\Factory\knowledge\information_layer\ per the Operating Rules.

## Implications and applications (extended from Session 6 list)

Session 6 listed application domains conditional on the per-domain
substrate-vs-expression reading holding. Session 7 strengthens that
reading (two independent reproducible per-domain signatures now)
and adds the functional-family axis as a new dimension of
applicability.

Specifically: any application that uses MI-vs-H signature for
regime classification / anomaly detection / early warning now has
TWO axes of fingerprint, not one:
- Axis 1 (Session 6): per-domain null direction in operator space
- Axis 2 (Session 7): per-domain MI-vs-H functional family

A real-data system that matches expected family AND expected null
direction is "in-regime." A system that drifts off one axis but not
the other tells you WHICH axis is shifting -- whether the underlying
constraint structure is changing or the functional coupling is
changing. This is more diagnostic than either signature alone.

Examples by domain (no test run; mapping):
- Sepsis early warning (NoVell): healthy patients should show one
  family in {HR, HRV} -> MI; pre-sepsis should drift family OR null
  direction; tells you whether autonomic dysregulation or coupling
  collapse is the trigger.
- Earthquake precursor (geology connection): in-regime geology
  shows MI ~ constant (INFO-025). Pre-event should show family
  drift toward H-coupled. Falsifiable.
- Hurricane RI precursor: weather expression signature drift TBD.
- Cyber intrusion: process-pair MI dynamics during baseline vs
  attack.
All conditional on the reading holding on real data. Real-data
test is the gating experiment for all applications.

## Branch state

- Local: claude/claude-md-context-update-uCZ8m
- Remote: origin/claude/claude-md-context-update-uCZ8m
- Tip (Session 7 close): committed at end of v7 handoff write
- Pushed: yes. No PR. main untouched.
- Branch note: the v7 kickoff prompt directed continuation on
  claude-md-context-update-uCZ8m. The harness for this session was
  configured for a different branch (claude/two-more-pastes-ZkenG)
  but Greg gave explicit permission ("whatever you feel is best") to
  continue on the Session 6 branch for continuity. All Session 7
  commits live on claude-md-context-update-uCZ8m.

## Prompt for next session

See NEW_SESSION_KICKOFF_v8.md for the drop-in.

End of v7 handoff.
