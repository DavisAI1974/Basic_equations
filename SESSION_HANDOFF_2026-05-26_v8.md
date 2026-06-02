# Information Layer Session Handoff -- 2026-05-26 Session 8

## Read me first

Read in this order:
1. CLAUDE.md (project root). Has Session 7 note at the top
   alongside Sessions 3-6 notes. Six Operating Rules in force from
   Sessions 4-7 (no pre-assigned meaning, probe-not-falsifier,
   speaking posture before/after, Rule D incomplete-not-wrong,
   They never stacked, Treat literature as conjecture by default).
   All six bear on Session 8. NO NEW RULE added Session 8.
2. This file. Session 8 ran items 1 (four-force probe), 4
   (robustness on INFO-025), and 3 (mapping campaign). Item 2
   (storm/waves real-data) NOT started -- queued for Session 9
   when real data is staged.
3. SESSION_HANDOFF_2026-05-26_v7.md (Session 7 record).
4. SESSION_HANDOFF_2026-05-26_v6.md (Session 6 record) and
   SESSION_HANDOFF_2026-05-26_v5.md (Session 5 record) for the
   per-domain stack methodology, literature scan, and probe
   discipline.

## What Session 8 did

Greg's directive at start of session: order 1, 4, 3, 2 from the v8
kickoff menu. Two of the top three items completed in full; storm/
waves real-data (item 2) deferred to Session 9 pending real data
download/staging.

All work pushed to branch `claude/two-more-tasks-O6ahb` (harness
branch; merged in Session 7 work from
`claude/claude-md-context-update-uCZ8m` at start of session for
continuity).

### Experiment 1 -- Four-force unification probe (four_force_probe.py + four_force_pysr.py)

Greg's "real test" from Session 7 close. Procedure: build four
2-channel toy simulators capturing the characteristic structural
feature of each force, then apply the per-domain stack from
Sessions 6/7 (KBK + AI Poincare + SINDy on ensemble-H operator
matrix) plus PySR for functional family.

Force-field caricatures (NOT QFT predictions, structural caricatures):
- **EM**: long-range bilinear linear coupling, weak radiation damping.
  F_i = -omega_i^2 * x_i + K * x_j; massless mediator analog.
- **Weak**: Yukawa-suppressed bilinear coupling via massive mediator.
  F_i = -omega^2 * x_i + K * x_j * exp(-M * (x_i - x_j)^2); strong
  damping (W/Z short-lived).
- **Strong**: confining nonlinear coupling. F_conf = K_lin*dx +
  K_conf*dx^3; force grows with separation (flux-tube caricature).
- **Gravity**: universal energy-density-mediated attractive coupling.
  F_grav = G * E_total * dx / (dx^2 + eps_soft); E_total = x1^2 +
  v1^2 + x2^2 + v2^2; softened 1/r^2.

#### KBK stack result (cross-seed 11, 22)

ALL four forces share the [2,3,4] null direction at cos > 0.997
cross-force, both seeds. The shared direction is approximately
(-0.40, -0.41, +0.81), which is Session 3's (-1, -1, +2)/sqrt(6)
attractor — the SIGN-FLIPPED twin of Session 5/6's (+1,+1,+2)/
sqrt(6) protocol artifact. (-1,-1,+2) corresponds to algebraic
identity -(H_a - H_b)^2 ~ 0, i.e., the two channels are correlated.

EM is slightly distinct in 6D: it has a nontrivial MI coefficient
(+0.27) in the null direction, giving cos ~0.96 in 6D to the
others (vs ~0.998 weak-strong-gravity to each other). In [2,3,4]
projection EM is at cos > 0.999 to the others.

Cross-seed coefficient reproduction <5% within force. Eigenvalue
floors: EM 9.2e-5, weak 1.8e-6, strong 2.7e-7, gravity 9.7e-8.

#### PySR result (cross-seed 11, 22, n=600)

Low-complexity Pareto front (c=5-6):
```
                seed 11                          seed 22
EM       (H_a - H_b)^2 + 0.20  c=6      (H_b - H_a)^2 + 0.21  c=6
weak     0.305 + 0.063*H_a     c=5      0.317 + 0.068*H_a     c=5
strong   exp(H_a^2) * 0.227    c=5      H_b^4 + 0.246          c=5
gravity  0.336 - 0.572*H_b^2   c=6      0.355 - (H_a-0.11)^2   c=6
```

EM and weak reproduce across seeds at coefficient level (<5%).
Strong and gravity reproduce only at functional-form level
(power-of-one-channel or constant-minus-quadratic) but PySR breaks
the H_a/H_b channel symmetry randomly between seeds. This is
because the strong and gravity caricatures are exactly symmetric
in H_a vs H_b at the simulator level (both channels enter the
force law symmetrically), while EM and weak are not (asymmetric
omegas, asymmetric coupling-amplitude factors).

**Cross-link to INFO-025**: EM's functional family is the SAME as
Session 7 physics Duffing (INFO-025: physics MI = (H_b - H_a)^2 +
0.28). Both are linear-restoring nonlinear-coupling oscillator
dynamics. Weak's functional family is the SAME as Session 7
chemistry Brusselator (INFO-025: chemistry MI = 0.73*H_a + 1.08).
Yukawa suppression collapses to effective weak linear coupling.

#### Joint reading (three levels, separated per Result Discipline)

A. *Data level*: Four force-field caricatures share substrate-
level signature (cos > 0.997 cross-force on [2,3,4] null) AND
differ at expression-level (four PySR functional families with the
caveat that strong/gravity are channel-symmetry-broken). The
shared substrate is the Session 3 attractor direction; expression
families differ.

B. *Methodological*: PySR functional-family reproducibility
across seeds requires the underlying dynamics to break the channel
symmetry. When the dynamics is exactly H_a/H_b symmetric (strong,
gravity), PySR breaks symmetry randomly per seed. EM and weak (and
all four Session 6/7 simulators) had asymmetric dynamics and
reproduced cleanly. Methodological consequence: when running PySR
on a new dataset, check the underlying symmetry before reading
coefficient drift across seeds as instability.

C. *Substantive on substrate-vs-expression frame*: this is now the
THIRD independent reproducible per-domain signature of the frame,
after INFO-023 (per-domain null direction) and INFO-025 (per-
domain MI-vs-H functional family). Three supporting data points;
still a frame, not a claim. The frame consistently survives
generalization from "biology/physics/chemistry/geology" to
"EM/weak/strong/gravity" caricatures.

D. *Substantive on four-force unification (Greg's "real test")*:
the toy four-force probe is consistent with the substrate-vs-
expression reading -- one shared substrate, four distinct
expressions -- but the caricatures are not QFT and do not
constitute evidence about real four-force unification. The
caricatures DO show that "share substrate but differ in
expression" is a coherent geometric possibility in the operator
basis. A real-data probe (PDG coupling-running data, LIGO strain
subsets, LHC kinematics) is needed to test on actual physics.
Queued for Session 9.

### Experiment 4 -- Robustness check on INFO-025 (robustness_info025.py)

Procedure: re-simulate biology Lotka-Volterra under a grid of
robustness conditions, fit PySR + matched exp baseline + matched
linear baseline. Conditions:
- T sweep: 10, 30 (baseline), 100
- N_ens sweep: 200, 600 (baseline), 1200
- Observation noise sweep: 0.0 (baseline), 0.10, 0.50

7 conditions x 2 seeds = 14 PySR fits.

#### Result table (exp baseline A * exp(H_a / B) coefficients)

```
                                    seed 11                seed 22
T=10  N=600 noise=0    A=0.26 B=3e7 R^2=0.00  A=0.24 B=17  R^2=0.11
T=30  N=600 noise=0    A=0.58 B=2.38 R^2=0.73 A=0.48 B=1.92 R^2=0.82  <-- baseline
T=100 N=600 noise=0    A=0.21 B=1.14 R^2=0.57 A=0.21 B=1.16 R^2=0.68
T=30  N=200 noise=0    A=0.73 B=2.95 R^2=0.67 A=0.82 B=3.81 R^2=0.64
T=30  N=1200 noise=0   A=0.50 B=2.12 R^2=0.75 A=0.56 B=2.37 R^2=0.71
T=30  N=600 noise=0.10 A=0.18 B=1.13 R^2=0.87 A=0.17 B=1.08 R^2=0.85
T=30  N=600 noise=0.50 A=0.04 B=0.96 R^2=0.78 A=0.22 B=-5e7 R^2=0
```

PySR low-complexity Pareto front consistently includes exp(x0)*const
at c=4 across conditions with signal.

#### Joint reading

A. *Data level*: The FUNCTIONAL FAMILY A * exp(H_a/B) survives all
conditions with sufficient signal (T >= 30). The COEFFICIENTS A
and B drift heavily with regime.

B. *INFO-025 reproduction*: baseline (T=30, N_ens=600, no noise)
gives cross-seed A = 0.53 +/- 0.05, B = 2.15 +/- 0.23, matching
INFO-025's A=0.5, B=2.0 within seed-to-seed variation. INFO-025
reproduces.

C. *Regime drift*:
   - T=10: insufficient phase exploration -> no signal
   - T=100: trajectory leaves small-deviation regime -> A and B
     both ~halve. Family persists but parameter regime shifts.
   - N_ens=200: noisier estimate, coefficients drift up (A~0.78)
   - N_ens=1200: TIGHTEST match to INFO-025 (A=0.50 exact)
   - noise=0.10: broadened H_a flattens apparent exp, A drops to
     ~0.18 but R^2 stays high
   - noise=0.50: marginal; one seed broken

D. *Methodological consequence*: when citing INFO-025 coefficients,
specify the regime. The structural family signature
(A * exp(H_a/B)) is robust; the specific A, B values require
N_ens >= 600 and T = 30 to land at the INFO-025 point.

### Experiment 3 -- Mapping campaign (mapping_campaign.py)

Procedure: vary one dominant knob per Session 6/7 domain across 3
values, run PySR, check whether the INFO-025 functional family
persists, mutates, or disappears.

Knobs:
- physics_duffing K (neighbor coupling): 0.05, 0.20 (baseline), 0.50
- biology_lotka beta (predation rate): 0.30, 0.50 (baseline), 0.80
- chemistry_brussel B (autocatalytic param): 2.0, 3.0 (baseline), 4.0
- geology_bk drift_rate: 0.02, 0.04 (baseline), 0.08

4 domains x 3 values x 1 seed = 12 PySR fits (100.7s total).
Cross-seed scaling deferred to follow-up.

#### Result (low-complexity Pareto front per condition)

```
physics_duffing K (baseline K=0.20):
  K=0.05  c=6 loss=1.3e-3  -0.17*H_a^2 + 0.204               [single-channel]
  K=0.20  c=6 loss=6.4e-3  (H_b - H_a)^2 + 0.275              [INFO-025 family!]
  K=0.50  c=6 loss=7.2e-3  -0.40*H_a^2 + 0.297                [single-channel]

biology_lotka beta (baseline beta=0.50):
  beta=0.30  c=5 loss=4.8e-2  sqrt(exp(H_a - 1.67))           [exp family]
  beta=0.50  c=5 loss=4.3e-2  sqrt(exp(H_a - 1.27))           [INFO-025 family!]
  beta=0.80  c=5 loss=3.2e-2  0.30*H_a + 0.50                 [LINEAR mutation]

chemistry_brussel B (baseline B=3.0):
  B=2.0  c=5 loss=4.0e-2  (1.29 / H_b) / H_a                  [RATIO mutation]
  B=3.0  c=5 loss=1.3e-2  0.71*H_a + 1.08                     [INFO-025 family!]
  B=4.0  c=5 loss=1.7e-1  1.32 + 0.25*H_a                     [linear, different coeff]

geology_bk drift_rate (baseline drift=0.04):
  drift=0.02  c=4 loss=4.3e-4  exp(H_b) + 0.132   (tiny loss diffs across complexity)
  drift=0.04  c=5 loss=4.4e-4  0.027*H_b + 0.273  (constant-ish, INFO-025 family)
  drift=0.08  c=5 loss=4.2e-4  0.196*H_a / H_b
```

#### Joint reading

A. *Data level*: the INFO-025 functional families ARE regime
signatures. At baseline values of the dominant knob, the family
PySR picks at low complexity matches the INFO-025 reading. At
off-baseline values, the family can mutate.

B. *Per-domain reading*:
   - Physics: (H_b - H_a)^2 + const is K=0.20-specific. At low K
     (weak coupling) the difference structure weakens; at high K
     (strong coupling) it also mutates toward single-channel
     quadratic. The family lives in a Goldilocks regime of
     coupling strength.
   - Biology: exp(H_a/B) persists across beta=0.30 -> 0.50; at
     beta=0.80 (strong predation), the family flips to linear in
     H_a. The exponential family is a low-to-moderate-predation
     signature; high predation linearizes it.
   - Chemistry: linear in H_a is B=3.0-specific. At B=2.0
     (sub-bifurcation), ratio content surfaces; at B=4.0
     (super-bifurcation), linear with shifted coefficients. The
     Brusselator's qualitative regime (B above/at/below the Hopf
     bifurcation at B = 1 + A^2 = 2) drives the family.
   - Geology: constant family is robust across drift rate; all
     three drift values give essentially flat MI (variance
     ~4e-4, no family stronger than constant + tiny noise).

C. *Methodological*: the dynamical knob selects the family, not
just the coefficients. INFO-025 is a fingerprint of the SPECIFIC
baseline configuration of each simulator. Off-baseline, families
can mutate qualitatively (exp -> linear; quadratic-difference ->
single-channel quadratic; linear -> ratio). This refines INFO-029
(functional family is structural; coefficients are regime-
dependent) to: the family ITSELF is regime-dependent for large
deviations from baseline.

D. *Implications for the substrate-vs-expression frame*: the
expression-level family is regime-conditional. The "per-domain
expression of substrate" frame (Session 6/7) should be stated as
"per-domain-AND-per-regime expression of substrate." Whether the
underlying SUBSTRATE (rank-3 null subspace; shared attractor
direction in operator space) survives across the same knob sweep
is a follow-up question. The mapping campaign here ran PySR but
not the KBK stack across knob values; that's the natural next
extension.

## Ledger updates

**INFO-027 -- LOCATED FINDING (Session 8, new)**: Four force-field
caricatures (EM, weak, strong, gravity) share the [2,3,4] null
direction at cos > 0.997 cross-force, sitting on (-1, -1, +2)/
sqrt(6) Session 3 attractor (the sign-flipped twin of Session 5/6
(+1, +1, +2)/sqrt(6) protocol artifact; corresponds to algebraic
identity -(H_a - H_b)^2 ~ 0). EM differs from weak/strong/gravity
in 6D space (cos ~0.96 to others) by carrying a nontrivial MI
coefficient (+0.27) in null direction. Cross-seed coefficient
reproduction <5%. PySR cross-seed low-complexity (c=5-6):
- EM:      MI = (H_a - H_b)^2 + (0.20-0.21)   [reproducible]
- Weak:    MI = (0.31-0.32) + (0.063-0.068)*H_a   [reproducible]
- Strong:  one-channel-power + const [form reproducible, channel
           not -- dynamics is H_a-H_b symmetric]
- Gravity: const - one-channel-quadratic [form reproducible,
           channel not -- dynamics is H_a-H_b symmetric]
EM family matches Session 7 physics Duffing (INFO-025); weak
family matches Session 7 chemistry Brusselator (INFO-025). Third
independent reproducible per-domain signature of substrate-vs-
expression frame after INFO-023 (per-domain null direction) and
INFO-025 (per-domain MI-vs-H functional family).

**INFO-028 -- METHODOLOGICAL (Session 8, new)**: For PySR
functional-family extraction, cross-seed reproducibility at
coefficient level requires the underlying dynamics to break
H_a-H_b channel symmetry. When the dynamics is exactly
channel-symmetric (Session 8 strong, gravity caricatures), PySR
breaks symmetry randomly per seed -- different seeds select
different channels, giving "non-reproducible" coefficients while
the FUNCTIONAL FORM is consistent. When the dynamics is
asymmetric (EM, weak, all four Session 6/7 simulators),
coefficients reproduce within <5%. Consequence: when reading
"seed-dependent" or "non-reproducible" output from PySR, check
the underlying symmetry of the dynamics before reading instability.

**INFO-029 -- METHODOLOGICAL (Session 8, new)**: INFO-025's
biology MI ~ 0.5 * exp(H_a/2) survives a robustness sweep across
T (10, 30, 100), N_ens (200, 600, 1200), and observation noise
(0.0, 0.10, 0.50) at the FUNCTIONAL FAMILY level (PySR
consistently surfaces exp(x0)*const at low complexity across all
conditions with signal). The COEFFICIENTS A and B drift heavily:
   T=10 -> no signal
   T=30, N_ens=600 (baseline) -> A=0.53+/-0.05, B=2.15+/-0.23
   T=30, N_ens=1200 -> A=0.50 exact (tightest match to INFO-025)
   T=100 -> A and B halve (parameter regime shifts as system
            explores more of phase space)
   noise=0.10 -> A drops to ~0.18 (broadened H_a flattens exp)
INFO-025 reproduces at baseline (T=30, N_ens=600, no noise).
Structural family signature is robust; specific A, B values
require the same regime. Operational consequence: when citing
INFO-025 coefficients, specify the regime.

**INFO-030 -- LOCATED FINDING (Session 8, new)**: PySR functional
families from INFO-025 are baseline-specific regime signatures
when the dominant simulator knob is varied. One-knob, three-value
sweep on each Session 6/7 domain (physics K, biology beta,
chemistry B, geology drift_rate), one seed (11). At baseline knob
values, the PySR low-complexity Pareto matches the INFO-025
reading. Off-baseline values mutate the family qualitatively:
- physics K=0.20: (H_b - H_a)^2 + const; K=0.05 and K=0.50
  mutate to single-channel quadratic.
- biology beta=0.30 to 0.50: exp(H_a/B) family; beta=0.80 mutates
  to linear in H_a.
- chemistry B=3.0: linear in H_a; B=2.0 mutates to ratio
  (1.29/H_b)/H_a; B=4.0 keeps linear but shifts coefficient.
- geology drift_rate: constant family holds across sweep (MI
  variance ~4e-4 throughout; sweep does not move the system
  out of the decoupled regime).
Refines INFO-029: the functional FAMILY itself is regime-
conditional for large deviations from baseline, not just the
coefficients. Operational consequence: the "per-domain expression
of substrate" frame (Session 6/7) should be stated as "per-
domain-AND-per-regime expression of substrate." Cross-seed scaling
is the next extension; this finding is one-seed.

**INFO-031 -- LOCATED FINDING (Session 8 post-close, new)**:
Symmetry-swap test of EM vs strong force (Greg's "quick glance at
em/sf" at session close). Four configs x 3 seeds:
- em_asym_baseline   (omega1=1.0, omega2=1.2)
- em_sym_swap        (omega1=omega2=1.0)
- strong_sym_baseline (omega1=omega2=0.5)
- strong_asym_swap   (omega1=0.5, omega2=0.7)

Substrate level (KBK v_null): em_asym, em_sym, strong_sym all hit
(-0.41, -0.41, +0.82); strong_asym shifts to (-0.26, -0.59, +0.77).
Cross-seed cos > 0.998 within each config (substrate stable per
config; symmetry change moves substrate slightly only for strong).

Expression level (PySR cross-seed):
- em_asym:    (H_a - H_b)^2 + 0.21    [reproducible 3/3 seeds]
- em_sym:     (const - channel^2)^2   [channel randomly assigned;
                                       2/3 chose H_b, 1/3 H_a]
- strong_sym: channel-power+const     [channel randomly assigned;
                                       2/3 H_a, 1/3 H_b]
- strong_asym: exp((H_b - H_a) - 1.4) [reproducible 3/3 seeds]

INFO-028 directly confirmed via swap: symmetric dynamics break
PySR cross-seed coefficient reproducibility through random channel
assignment; making dynamics asymmetric RESTORES reproducibility.

SUBSTANTIVE: EM and strong produce DIFFERENT functional families
on top of SHARED substrate even when channel symmetry is
controlled. EM = (H_a - H_b)^2 + const (POLYNOMIAL in difference).
Strong = exp((H_b - H_a) - const) (EXPONENTIAL in difference).
Both express the same substrate algebraic identity
-(H_a - H_b)^2 ~ 0 (the Session 3 attractor direction) but realize
it through different functional forms.

Reading: the dynamical structure is visible in the form of the
expression. EM's linear-restoring + linear-coupling produces
polynomial; strong's confining-cubic produces exponential. The
"same substrate, different expression" pattern is NOT arbitrary --
the form of the expression carries the linear vs nonlinear
character of the underlying dynamics.

Implication for the four-force question (Greg's "real test"): EM
and strong are NOT the same expression even at the toy-caricature
level, when dynamical structure is exposed. This is consistent
with the substrate-vs-expression frame and with the empirical
fact that electroweak unification is established (EM+weak share
dynamical structure) while strong sits separately (genuinely
different structure).

**INFO-032 -- LOCATED FINDING (Session 8 post-close, new)**:
Gravity caricature substrate signature is CONFIGURATION-DEPENDENT
in a way EM/strong (INFO-031) are not. Three-config probe with 3
seeds each: gravity_sym_baseline (omega1=omega2=0.8, universal
E_total coupling), gravity_asym_omegas (omega1=0.8, omega2=1.0,
universal), gravity_nonuniversal (asym omegas, coupling only to
E1, not E_total).

Cross-config v_null 6D mean cosines:
  sym <-> asym_omegas       = +0.11   (different substrate!)
  sym <-> nonuniversal      = +0.02   (orthogonal!)
  asym_omegas <-> nonuniv   = +0.99   (essentially same)

Cross-seed cos > 0.99 within each config -- the substrate shift
is BETWEEN configs, not noise.

Substrate identity per config:
  gravity_sym + universal:    v_null in (H_a^2, H_b^2, H_a*H_b);
                              hits Session 3 attractor
                              (-1, -1, +2)/sqrt(6); algebraic
                              identity -(H_a - H_b)^2 ~ 0.
                              MI std ~0.11.
  gravity_asym + universal:   v_null FLIPS so MI coefficient
                              is +0.99 (dominant). Algebraic
                              identity becomes MI ~ const.
                              MI std collapses to ~0.03 (4x
                              reduction).
  gravity_asym + nonuniversal: same MI-dominant substrate
                              (cos +0.99 to asym universal).
                              MI std ~0.04.

PySR per config (cross-seed):
  gravity_sym:        const - one_channel^2  (channel randomly
                      assigned, 2/3 H_b, 1/3 H_a -- matches
                      INFO-028 for symmetric dynamics)
  gravity_asym:       0.21 + weak power-of-H_a (nearly constant
                      MI; PySR finds weak linear/cubic content)
  gravity_nonuniv:    0.21 + small linear in H_a or in (H_a - H_b)

Contrast with INFO-031 (EM/SF glance): EM and strong both kept
their substrate direction across symmetry swap. Only their
EXPRESSION-level family changed. Gravity's substrate ITSELF is
configuration-fragile.

Frame-level reading (no claim, frame only): gravity's universal
energy-mediated coupling under mass asymmetry produces "MI is
approximately invariant" as the dominant structural constraint.
The information shared between channels saturates because
gravitational coupling preserves total energy without preferred
direction -- once the masses differ, the dynamics balance into
a configuration where MI doesn't grow or shrink over time.
For EM/strong/weak the coupling has a preferred mode (electric
field, color confinement, weak interaction) that gives MI a
specific functional shape regardless of channel mass equality;
gravity does not have such a preferred mode.

What gravity "is" in this framework (frame only, no claim):
a coupling whose substrate signature in operator space depends
on configuration -- "channel correlation" when channels are
mass-degenerate, "MI invariance" when they are not. Genuinely
different signature pattern from the other three forces.

Real-data probe (LIGO strain, orbital data, binary pulsar timing)
remains the gating test for whether this carries to actual
gravitational physics.

## Substantive notes carried forward

**Four-force probe is a frame-level test, not a claim-level test.**
The caricatures (EM bilinear linear, weak Yukawa, strong cubic,
gravity universal energy-mediated) are structural toy models. They
do NOT predict real QFT behavior. They DO show that "shared
substrate, distinct expression" is a coherent geometric pattern in
the operator basis for any set of coupled-pair dynamical systems
with sufficiently different coupling structures. Whether this
pattern holds on real four-force data is the next question --
queued for Session 9 with PDG coupling-running data as the
candidate first dataset.

**Real-data probe is the gating step for any application claim.**
The Session 8 work strengthens the substrate-vs-expression frame
within simulators (three independent signatures now). The
applications discussion from Session 6 (sepsis, earthquake,
hurricane RI, cyber intrusion) remains conditional on real-data
reproduction of the simulator findings.

## Experiments queued for Session 9 (priority order)

(1) **Storm/waves real-data probe**. Greg raised this in Session 7
    close and ranked it #4 in v8 kickoff (so #2 within the queue
    after four-force). Real-data per-domain stack on NOAA buoy /
    HRRR atmospheric reanalysis / hurricane track data. Tests
    substrate-to-expression frame on a real physical-expression
    system. Greg's framing: "how does pure physics get expressed
    as a storm or waves?"

(2) **Four-force real-data probe (PDG coupling-running data)**.
    The simulator-level four-force probe (Session 8 INFO-027)
    showed shared substrate + distinct expression on toy
    caricatures. The real-data version: pull PDG running of g1,
    g2, g3 with energy, construct ensemble via published
    measurement uncertainties or renormalization-scheme variations,
    apply the per-domain stack. Tests whether real EM/weak/strong
    behave like the caricatures.

(3) **Gravity probe (LIGO O3/O4 strain subsets)**. Within reach as
    a real-data extension of the four-force probe. Gravity strain
    signal vs detector noise; treat strain channels (H1, L1) as
    coupled pair; run per-domain stack.

(4) **PhysioNet ECG cross-domain probe (NoVell extension)**. Carl
    Saab angle. HR/HRV pair as coupled-pair data; tests whether
    biology functional family from INFO-025 survives on real ECG
    data.

(5) **Symmetry-asymmetry sweep on Session 8 caricatures**. From
    INFO-028: PySR cross-seed reproducibility depends on dynamics
    H_a-H_b symmetry. Add asymmetric perturbation to strong and
    gravity simulators (different omegas, different masses, etc.)
    and check whether functional family reproducibility recovers.

(6) **Mapping campaign extension** (item 3 follow-up): two-knob
    sweeps per domain, cross-seed scaling. The Session 8 mapping
    is one-seed and one-knob per domain; need second seed and
    second knob to nail down the dynamical property that selects
    family.

(7) **Five-physics-areas substrate hunt** (Greg's Session 5
    reorientation, still open). Pick 5 physics areas, find
    structural element that appears in all 5.

(8) **Methodology stacks from v5 menu** (still open):
    complexity-entropy plane (Rosso); Yamada 2023 SVD with
    entropy as field; Gomez-Herrero 2015 ensemble infrastructure;
    Stochastic Electrodynamics revisit; Jacobson 1995 cross-
    domain extension.

## Methodological notes from this session

(Not promoted to Rules. Domain-knowledge refinements.)

**Force-field caricatures as a probe template**. A 2-channel
toy simulator with the characteristic coupling structure of a
candidate domain is a useful first-pass probe before real-data.
The Session 8 four-force caricatures (EM linear, weak Yukawa,
strong confining, gravity universal) took ~5 minutes to write
and run end-to-end. Same template applies to other candidate
domains (e.g., turbulence with cascade coupling; financial
markets with leverage feedback; epidemic SIR with mass-action
coupling). Build the caricature first, get the family signature,
then go to real data.

**Symmetry check before reading PySR cross-seed variation**
(INFO-028 generalization). Before judging PySR coefficient
non-reproducibility as instability, check whether the underlying
dynamics has a relevant symmetry (channel-channel, time-reversal,
spatial). Symmetric dynamics -> PySR breaks symmetry randomly
between seeds; report at the functional-form level, not the
coefficient level.

**Regime-specification when citing coefficients** (INFO-029
generalization). Specify (T, N_ens, noise, dt) when citing PySR-
discovered coefficients. The structural family is what survives
across regimes; the specific numbers don't.

## Files produced this session (all in repo root)

Scripts:
- four_force_probe.py        (Experiment 1: four-force caricatures
                              + KBK stack)
- four_force_pysr.py         (Experiment 1 follow-up: PySR on the
                              four caricatures)
- robustness_info025.py      (Experiment 4: T / N_ens / noise
                              sweep on biology, PySR + exp + linear
                              baselines)
- mapping_campaign.py        (Experiment 3: per-domain knob sweeps,
                              one knob per domain, 3 values)

Result JSONs:
- four_force_probe_canary.json, four_force_probe_results_seed11.json,
  four_force_probe_results_seed22.json
- four_force_pysr_canary.json, four_force_pysr_results.json
- robustness_info025_canary.json, robustness_info025_results.json
- mapping_campaign_canary.json, mapping_campaign_results.json

Run logs (matching the above):
- four_force_probe_canary.log, four_force_probe_run_seed11.log,
  four_force_probe_run_seed22.log
- four_force_pysr_canary.log, four_force_pysr_run.log
- robustness_info025_canary.log, robustness_info025_run.log
- mapping_campaign_canary.log, mapping_campaign_run.log

Docs:
- SESSION_HANDOFF_2026-05-26_v8.md (this file)
- CLAUDE.md updates (Session 8 note added)
- NEW_SESSION_KICKOFF_v9.md (drop-in for next session)

Dependencies: same as Session 7 (numpy, scikit-learn, scipy, pysr
1.5.10 with Julia 1.11.9 backend). PySR Julia precompile took
107s on first invocation in this session; subsequent runs reuse
the compiled cache.

To be mirrored to E:\information_layer\ and
F:\Factory\knowledge\information_layer\ per the Operating Rules.

## Branch state

- Local: claude/two-more-tasks-O6ahb (harness branch)
- Remote: origin/claude/two-more-tasks-O6ahb
- Session 7 work pulled in from
  origin/claude/claude-md-context-update-uCZ8m at session start
  (fast-forward merge). All Session 8 commits live on the harness
  branch per the durable instructions.
- Pushed: yes. No PR. main untouched.

End of v8 handoff.
