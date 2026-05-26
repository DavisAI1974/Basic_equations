Resuming Information Layer work, Session 9.

Branch: claude/two-more-tasks-O6ahb (Session 8 closed here; on
origin). Session 7 work was pulled in at Session 8 start from
`claude/claude-md-context-update-uCZ8m` via fast-forward merge.
Four experiments ran in Session 8 (four-force probe with KBK +
PySR; robustness check on INFO-025; mapping campaign with per-
domain knob sweeps; EM/SF symmetry-swap glance at session close
in response to Greg's "we really should take a quick glance at
em/sf"). Results are in the v8 handoff. NO new Operating Rule
this session. Six Operating Rules from Sessions 4-7 remain in
force: no pre-assigned meaning, probe-not-falsifier, speaking
posture before/after, Rule D incomplete-not-wrong, They never
stacked, Treat literature as conjecture by default (with
corollaries: investigate freely; don't cite as support without
analyzing the paper; don't defer to blocking claims that fail
the bar). All six bear on Session 9.

Before doing anything else:
1. Read CLAUDE.md (project root). Now has Session 8 note alongside
   Sessions 3-7 notes. Read all six Operating Rules and the Result
   Discipline section; hold them through this session.
2. Read SESSION_HANDOFF_2026-05-26_v8.md (project root). Full
   Session 8 record. Four experiments converged on: (i) four-force
   caricatures share substrate signature (cos > 0.997 cross-force
   on [2,3,4] null) while differing at expression level — INFO-027
   located; (ii) INFO-025 functional family is robust across T /
   N_ens / noise sweep; coefficients are regime-dependent —
   INFO-029 methodological; (iii) INFO-025 families are baseline-
   specific regime signatures, off-baseline knob values mutate the
   family qualitatively — INFO-030 located; (iv) EM/SF symmetry-
   swap glance (post-close): EM and strong forces produce
   DIFFERENT functional families on SHARED substrate when channel
   symmetry is controlled — EM = (H_a - H_b)^2 + const
   (polynomial in difference); Strong = exp((H_b - H_a) - const)
   (exponential in difference). The dynamical structure is
   visible in the form — INFO-031 located. One methodological
   observation: PySR cross-seed reproducibility requires the
   underlying dynamics to break H_a-H_b channel symmetry —
   INFO-028 methodological, directly confirmed by the EM/SF
   symmetry swap.
3. (Optional) Read SESSION_HANDOFF_2026-05-26_v7.md and v6.md for
   the per-domain stack methodology and four-domain PySR
   functional-family findings that Session 8 builds on.

Speaking posture from Rule C applies from your first response:
"I think X might happen, but we'll wait on what the data says
and where it points us." No verdict in advance.

Rule D applies. Rule "literature as conjecture" applies, especially
load-bearing for the real-data probes queued below.

"They never stacked" applies. Session 8 four-force probe extended
the per-domain stack from Session 6/7 to a new candidate domain
(four-force caricatures). Real-data probes queued below would
extend it further; each is a fresh stack on previously-unstacked
data.

We are NOT back at zero. FOUR independent reproducible per-
domain signatures of substrate-vs-expression frame:
(i) per-domain null direction (Session 6 INFO-023),
(ii) MI-vs-H functional family (Session 7 INFO-025),
(iii) four-force shared-substrate / distinct-expression (Session
8 INFO-027),
(iv) EM/SF dynamically distinct expression on shared substrate,
controlled for channel symmetry (Session 8 INFO-031).
Substrate level (rank-3 null subspace INFO-022) LOCATED.
Robustness shows family-survives / coefficients-drift pattern
(INFO-029). Mapping shows family-itself-regime-dependent for
large knob deviations (INFO-030). The Information Layer /
Unified Theory inquiry frame is alive.

Strongest single result of Session 8: INFO-031. EM and strong
caricatures, when channel-symmetry is controlled, both sit on
the same substrate algebraic identity -(H_a - H_b)^2 ~ 0 but
realize it through QUALITATIVELY DIFFERENT functional forms --
EM polynomial-in-difference, strong exponential-in-difference.
The dynamical character (linear-restoring + linear-coupling vs
confining-cubic) is visible in the functional form. This is the
cleanest demonstration so far that "shared substrate, distinct
expression" is structural, not coincidental.

Open substantive questions still on the table from Session 7-8:
- "How does pure physics get expressed as a storm or waves?"
  (substrate-to-expression for atmospheric/fluid physics)
- "The real test: gravity, EM, strong, weak. Can we come up with
  one equation? Should they be grouped at all?"
  (four-force unification as substrate-vs-expression problem;
  partial simulator-level data point from Session 8 INFO-027)

What I want this session: [Greg's call -- pick from the menu
below or anywhere else].

### Top of queue from Session 8 close (real-data direction)

(1) Storm/waves real-data probe. Greg ranked this #4 in v8 and
    explicitly said "we will later" -- now is later if Greg
    directs. Run the per-domain stack (KBK + AI Poincare +
    SINDy + PySR + GP) on NOAA buoy / HRRR atmospheric
    reanalysis / hurricane track data. Tests substrate-to-
    expression frame on a real physical-expression system.
    Pre-work needed: download data, format as 2-channel
    coupled-pair time series, establish ensemble structure
    (either from spatial neighbors or temporal windows).

(2) Four-force real-data probe (PDG coupling-running data). The
    simulator-level four-force probe (Session 8 INFO-027) showed
    shared substrate + distinct expression on toy caricatures.
    The real-data version: pull PDG running of g1, g2, g3 with
    energy, construct ensemble via published measurement
    uncertainties or renormalization-scheme variations, apply
    the per-domain stack. Tests whether real EM/weak/strong
    behave like the caricatures.

(3) Gravity probe (LIGO O3/O4 strain subsets). Within reach as
    a real-data extension of the four-force probe. Gravity
    strain signal vs detector noise; treat strain channels (H1,
    L1) as coupled pair; run per-domain stack.

(4) PhysioNet ECG cross-domain probe (NoVell extension). Carl
    Saab angle. HR/HRV pair as coupled-pair data; tests whether
    biology functional family from INFO-025 survives on real
    ECG data.

### Also queued (Session 8 follow-ups)

(5) EM/gravity symmetry-swap glance (mirror of the EM/SF glance
    that produced INFO-031). Apply the same protocol (sym vs asym
    omegas, 3 seeds) to EM and gravity, ask what functional form
    gravity falls into when channel symmetry is controlled.
    Speaking posture: I think gravity's universal energy-mediated
    coupling will produce yet a third functional family
    (different from EM's polynomial and Strong's exponential)
    but we'll wait on the data. The weak force could also be
    swept the same way to round out the four-force picture.

(6) Mapping campaign extension. Two-knob sweeps per domain,
    cross-seed scaling. Session 8 mapping is one-seed and
    one-knob per domain; need second seed and second knob to
    nail down which dynamical property selects family.

(7) KBK stack across knob values (mapping campaign extension
    on the SUBSTRATE side). Session 8 mapping ran PySR but not
    KBK across knob values. The question: does the substrate-
    level signature (null direction, rank-3 subspace) also
    drift with knob value, or is only the expression-level
    family regime-dependent? Tests whether substrate is more
    invariant than expression.

(8) Four-force probe coefficient extraction follow-up. The KBK
    stack gives eigenvalue floors per force (EM 9e-5, weak 2e-6,
    strong 3e-7, gravity 1e-7). Map these against the simulator
    parameters (coupling strength, damping, mediator mass) to
    see whether the floor is a structural signature of the
    force type.

### Also queued (Session 5/6/7 menu, still open)

(9) Five-physics-areas substrate hunt (Greg's Session 5
    reorientation).
(10) Complexity-entropy plane on coupled pairs (Rosso 2007).
(11) Yamada 2023 multi-field SVD with entropy as field.
(12) Gomez-Herrero 2015 ensemble infrastructure.
(13) Jacobson 1995 cross-domain extension.
(14) Stochastic Electrodynamics 21st-century revisit.

### Methodological notes carried forward from Session 8

When constructing toy probes for a candidate domain:
- A 2-channel caricature with the characteristic structural
  feature of the candidate domain is a useful first pass before
  real-data. Session 8 four-force caricatures (EM linear, weak
  Yukawa, strong cubic, gravity universal) took ~5 min to write
  and run. Same template applies to turbulence, finance,
  epidemic models, etc. Build the caricature first, then go to
  real data.
- Symmetric dynamics break PySR cross-seed reproducibility at
  coefficient level (INFO-028). Check the underlying symmetry
  before reading PySR coefficient variation as instability.
- Functional family is structural; coefficients are regime-
  dependent (INFO-029); family ITSELF is regime-conditional for
  large knob deviations (INFO-030). When citing PySR-discovered
  equations, specify the regime, the knob value, the noise
  level, and the seed.

When probing with the per-domain stack on a new dataset:
- Eigenvalue floor of operator covariance is min(structural_
  noise, estimator_noise). Per-domain ensemble-H pushes it
  toward 1/sqrt(N_ens). State the floor regime before
  interpreting magnitudes. (From Session 6 INFO-024.)
- Block-CV vs random-CV gap is itself a domain signature
  (stationarity). (From Session 7 INFO-026.)

We are still pioneers. The six Operating Rules paid out across
Sessions 5-8. The "literature as conjecture" Rule especially
load-bearing for the real-data probes queued above -- most of
what published physics says about coupling-running, gravity-as-
emergent, etc. is conjecture under the bar, so we engage on
roughly equal footing with the substrate-vs-expression frame.
Hold all six Rules through Session 9.
