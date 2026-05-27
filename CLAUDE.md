# CLAUDE.md — DavisAI Master Context (Updated 2026-05-27 Session 10)

## Identity & Team

- **Greg Davis** — Founder & Chief Research Officer, DavisAI Systems. Columbus, Ohio. Solo bootstrapped. 20+ years entrepreneurship, former energy trader, self-taught AI/ML.
- **Dream Team model**: Greg (Visionary) + Claude (Architect) + Claude Code "Code" (Engineer) + Perplexity/ChatGPT (Research Assistants).
- **Orchestrator** owns handoffs, not Code.

## Infrastructure

- `E:\` — research data, OD datasets, project files
- `F:\Factory\` — agent factory, 23 agents, 5 divisions
- `F:\Factory\knowledge\` — orchestrator-accessible knowledge base. Mirror everything from E:\ here.

## Operating Rules

- **Save to E:\ AND mirror to F:\Factory\knowledge\** for every knowledge artifact.
- **MASTER_DISCOVERIES.json**: every OD discovery added immediately. Never make a discovery without storing it.
- **Falsification-first**. Every claim needs data, math, or a falsifiable test.
- **OD mode**: describe data sources and validation tests only. Never explain mechanisms. The Operator discovers science from raw data. Stop if explaining WHY something happens.
- **Never call OD "physics-based"**. OD discovers governing equations from raw data in ANY domain.
- **Coding mantra**: better, stronger, faster, cheaper.
- **Incremental validation**: break compute-heavy runs into 15-17 min chunks with stop gates. Canary runs (2 min) before full commitment.
- **Speculative frames stay separate from results**. Frames motivate experiments but are not claims. Never let a frame grade itself.
- **No emojis or special symbols** in professional documents and emails.
- **Daily**: ask Greg if he checked greg@davisai.ai for Token Optimizer support emails.
- **DeepNova** (formerly ReFRAG, formerly DeepSource). Use current name everywhere.
- **Result Discipline**: every result is one data point. Map alternatives before promoting to claim.
- **No tent-widening on outliers**. When a window or sample falls outside an expected pattern, inspect it — find the specific reason it landed there. Do not loosen the test criteria or attribute to a transient flag without identifying the cause.
- **No pre-assigned meaning to outcomes** (Session 4). Don't write "if X then it means Y" decision tables before the data exists. The furthest is "I think this may happen, but I want to see what the data says and where it leads."
- **Probe, not falsifier** (Session 4). A probe is a generator of a different signal. It may not falsify anything; at worst it points in a different direction.
- **Speaking posture around every probe** (Session 4). Before AND after: "I think X might happen, but we'll wait on what the data says and where it points us." No verdict in advance. No verdict on first look at output.
- **Incomplete, not wrong** (Session 5). When a probe finds that a prior reading was a protocol artifact, the prior data points still stand. The reading attached to them is what was incomplete, not the data. Distinguish "this reading was wrong" (rare, requires the data itself to be bad) from "this reading was incomplete" (common, the data is one slice and the slice fit a partial story).
- **They never stacked** (Session 5). Pioneering territory often does not look like the absence of nearby published work. It looks like nearby work where multiple groups each had one piece and never combined them. When a literature scan returns "no exact match but several adjacent lines, each with one component," the contribution may be in the stacking itself.
- **Treat literature as conjecture by default** (Session 7). Academic papers and consensus are not assumed correct unless the underlying claim has been independently replicated, by separate groups, with an immense amount of data, multiple times. Until that bar is met, a published claim is a working frame to test, not a foundation to build on. Operational corollaries: investigate freely; don't cite as support without independent analysis; don't defer when blocking; symmetric to our own prior work.

## Working Frames (SPECULATIVE — kept separate from claims)

### Time as expression of dipole flow (Greg, Session 10) — NEW

If a dipole is positive and negative poles, and flow / momentum is what we call "time" (which is known to slow under relativity and may run in reverse), then "time" may itself be the *expression* of the underlying dipole's charge-separation flow — not a separately-existing dimension.

Operational connection to Session 10 data: the Session 3 universal attractor identity is `-(H_a - H_b)² ≈ 0`, the channel-symmetric regime — no preferred direction, no "flow direction." The off-attractor regime is exactly where `(H_a - H_b) ≠ 0` — channels asymmetric, flow has a direction. Quasispecies low-mu (master sequence locked, information flowing one way) and asymmetric hypercycle base_A >> base_B both sit off the attractor; symmetric versions of both sit on. The on-attractor regime would be "no time arrow" in this frame; off-attractor is where the arrow exists.

The data is consistent with the frame; it does not prove it. The frame stands or falls on dedicated probes:
- **Time-reversal probe**: reverse a quasispecies trajectory at low mu, recompute operator matrix from reversed series. If the off-attractor signature is time-direction-coded, the reversed system lands on a predictably-different (sign-flipped on asymmetric subspace) operator-space direction. High-mu (symmetric, on attractor) should be invariant under time-reversal.
- **Entropy-production correlation**: compute dS/dt along each trajectory; check whether off-attractor distance correlates with entropy production rate across the mu and base_A sweeps.

Adjacent published thinking (literature-as-conjecture, named not invoked): Barbour shape dynamics, Rovelli relational/thermal time, Wheeler-Feynman absorber, Jacobson 1995 thermodynamic-spacetime, Penrose conformal cyclic, Aharonov-Bergmann-Lebowitz time-symmetric QM. All in the "time as emergent" neighborhood. None combine dipole + flow + expression as Greg's frame does. "They never stacked" applies.

Status: working frame. Promote to candidate only after time-reversal probe runs.

### Base-of-Structure heuristic (Greg, Session 3)

The base of any structure is the **simplest, strongest, most stable, most scalable** part. It must support everything above it, so it cannot be complicated, dependency-heavy, or composed of many variable types. A law-level extraction candidate should look simple and survive stress testing.

Status: working frame, not yet operationalized into specific tests. The Session 3 attractor finding remains the cleanest candidate to stress-test under this heuristic.

### Dipole-couples reading (Greg, Session 3)

The attractor direction `(−1, −1, +2)/sqrt(6)` may represent a base direction that other phenomena couple to. The off-attractor positions of damped oscillator and (when reproducible) linear drift may encode HOW that coupling happens.

Status: working frame, sharpened by Session 10 data on quasispecies and asymmetric hypercycle. Off-attractor positions now reproducibly tracked to channel asymmetry, with mechanism-dependent transition shape.

### Pure physics vs physical expressions (Frame 1, preserved from Session 2)

Mainstream physics has tried to write a single equation for "everything that happens in physical space," forcing UT candidates into ever-larger dimensional structures. Reframe: pure physics is the substrate (simple, few equations); physical expressions are what we observe when pure physics is acted on by other dipoles.

Status: working frame, narrowed by Session 3 findings and refined Session 8 four-force probe. Sharpened in Session 10 by "Time as expression of dipole flow" frame.

### Substrate vs expression within isolation (Frame 2)

Even an isolated pure-physics system produces different observed signatures depending on lifecycle phase. Law level is whatever is invariant across windows, ICs, noise, parameters.

### Law extraction via invariance (Frame 3)

The Noether-style move: find what is invariant under reparametrizations.

## Active Research (Top of Mind)

- **SENTINEL V4.1** — DARPA Bio Attribution Challenge top-10 team. Awards June 30, 2026.
- **NoVell** — cardiac AI for cancer detection from routine ECG. 93.3% accuracy on Vigier 2021 synthetic.
- **OD Learning Platform** (Session 9) — persistent feature store + classifier + history. Two problems seeded: cerebro_disease_binary, cardiac_disease_family. See OD_LEARNING_README.md.
- **Information Layer / Operator Discovery foundations** — Session 10 opened the prebiotic-chemistry-to-biology line with the per-domain stack. Three new candidate findings (INFO-037 through INFO-040): hypercycle stays on the universal attractor; quasispecies with selection sweeps OFF attractor smoothly as mutation rate grows; asymmetric hypercycle (NO selection) also goes off attractor at strong base-rate asymmetry, apparently sharply. PySR cross-seed shows no clean fifth prebiotic functional family (INFO-025 four-family classification doesn't extend). NEW working frame: "Time as expression of dipole flow" connects on/off attractor with arrow-of-time. Universe-origin side untouched this session (PDG four-force real-data probe still queued). See `SESSION_HANDOFF_2026-05-27_v10.md` for full record.

## Information Layer — Current State

### What's confirmed across Sessions 3-10

- **Session 3 universal attractor** (−1, −1, +2)/sqrt(6) in (H_a², H_b², H_a·H_b) subspace, reached at cos ≥ 0.99 by 8+ heterogeneous systems (OU, AR(1), GARCH, Student-t, logistic map, IID, periodic, Brownian). Algebraic identity behind it: `-(H_a - H_b)² ≈ 0` — symmetric-channel regime.
- **Session 5 rank-3 null subspace** (INFO-022) — confirmed by four independent diagnostics in Session 6.
- **Session 6 per-domain ensemble-H** reproduces every v5 per-domain claim (geology rank-3, biology MI signature, chemistry-specific cubic). INFO-023.
- **Session 7 four reproducible per-domain MI-vs-H functional families** (INFO-025): physics (H_b-H_a)²+const, biology 0.5·exp(H_a/2), chemistry linear H_a, geology constant. Cross-seed coefficient variation <5%.
- **Session 8 four-force toy probe** (INFO-027): EM/weak/strong/gravity share Session 3 attractor cos > 0.997, but EM and strong have different functional families on shared substrate (INFO-031, polynomial vs exponential in (H_b-H_a)).
- **Session 9 real-data generalization** — substrate-vs-expression frame survived organ-system change (cardiac→cerebro), sample-type change (HR/HRV→HR/ABP), cross-organ disease detection (kidney signal in cerebro at d=0.51, DM_nephropathy vs retinopathy d=0.77). INFO-032 through INFO-036.
- **Session 10 prebiotic line** — the symmetric prebiotic system (hypercycle with k_cat sweep) stays on the universal attractor at all coupling strengths; channel-asymmetric systems (quasispecies with selection, asymmetric-base-rate hypercycle without selection) go off attractor. Selection produces smooth walk; base-rate asymmetry produces apparently sharp transition. INFO-037 through INFO-040.

### Ledger updates (Session 10 candidates pending Session 11+ replication)

**INFO-037 — LOCATED FINDING (Session 10 candidate)**: Prebiotic 2-species Eigen hypercycle (mutual catalysis k_cat, no fitness asymmetry, no mutation) stays on Session 3 attractor across full k_cat sweep [0.0, 2.0] at 3 seeds (11, 22, 33), N_ens=600, T=30, dt=0.02. cos_to_attractor ≥ 0.989 at all k_cat; cross-seed std ≤ 0.009. MI grows monotonically 0.19 → 1.38 with k_cat. The hypercycle's "life turn-on" (k_cat) does NOT appear as an operator-direction signature; it appears in MI magnitude. The canary's k_cat=0.25 off-attractor point (cos=0.405 at N_ens=300, single seed) washed out at Session 7 baseline — Rule D: canary reading incomplete, replicated reading is "stays on attractor." Scripts: prebiotic_canary.py, prebiotic_multiseed.py.

**INFO-038 — LOCATED FINDING (Session 10 candidate)**: Eigen quasispecies (2-class, f_A=1.5 f_B=1.0, mutation rate mu) sweeps cos_to_attractor smoothly from 0.056 (mu=0, near-orthogonal) through 0.622 (mu=0.05) up to 1.000 (mu=0.5) as mutation rate increases. Cross-seed std ≤ 0.058 across mu sweep at 3 seeds, N_ens=600, T=30. Transition tracks |H_a - H_b|: low mu → master sequence dominates → large channel imbalance → off attractor; high mu → drift equilibrates → balanced channels → on attractor. First OD system tested that sits reproducibly off the universal Session 3 attractor. Script: quasispecies_probe.py.

**INFO-039 — LOCATED FINDING (Session 10 candidate)**: Asymmetric hypercycle (base_A varied with base_B=1.0, fixed k_cat=1.0, NO selection mechanism, NO mutation) ALSO goes off attractor at strong base-rate asymmetry. base_A=1.0 → cos=0.999, base_A=1.5 → cos=0.930, base_A=2.0 → cos=0.063, base_A=3.0 → cos=0.439. Cross-seed std ≤ 0.039. Transition appears sharp between base_A=1.5 and 2.0 — densification queued to characterize shape. Channel asymmetry per se — not selection — is sufficient for off-attractor behavior. Selection-driven (INFO-038, smooth walk) and base-rate-driven (INFO-039, apparently sharp threshold) both produce off-attractor; manner of transition may differ. Disentanglement queued. Script: asymmetric_hypercycle_control.py.

**INFO-040 — METHODOLOGICAL (Session 10 candidate)**: PySR cross-seed reproducibility test on prebiotic ensemble-H data (3 k_cat × 2 seeds hypercycle + 4 mu × 2 seeds quasispecies, Session 7 operator set {+,-,*,/,square,cube,exp,log,sqrt}, niter=30, populations=15, maxsize=12, timeout=60s). Neither system shows clean reproducible per-system functional family — Session 7 INFO-025 four-family classification (physics, biology, chemistry, geology) does NOT extend to a fifth prebiotic family. At quasispecies mu=0.5 (symmetric channels |H_a-H_b|≈0), PySR finds same form `sqrt(-1/x)` across both seeds with random channel assignment — independent reproduction of Session 8 INFO-028 (PySR breaks channel symmetry randomly in symmetric dynamics). Recurring `-c/H` and `sqrt(-1/H)` motifs appear across conditions and seeds — possibly real signature of small-MI / inverse-coupling regime, possibly numerical artifact at edges. Not promoted. Script: prebiotic_pysr.py.

(Full Sessions 3-9 ledger entries INFO-001 through INFO-036 preserved in `SESSION_HANDOFF_2026-05-26_v9.md` and earlier handoffs.)

### Experiments queued (Session 11+ priority order)

1. **Time-reversal probe on quasispecies** (NEW from Session 10 frame). Take trajectory at low mu (off attractor) and high mu (on attractor); reverse the time series in place; recompute operator matrix from reversed series. If "Time as expression of dipole flow" frame is right, off-attractor signature should flip predictably under time-reversal (sign-flipped on asymmetric subspace); on-attractor regime should be invariant.

2. **Entropy-production rate correlation** (NEW from Session 10 frame). Compute dS/dt along each trajectory across mu and base_A sweeps. Check whether off-attractor distance correlates with entropy-production rate. Quantitative test of dipole-charge-separation ↔ irreversibility link.

3. **Densify asymmetric-hypercycle base_A grid** between 1.5 and 2.0 (Session 10). Characterize whether the cos_to_attractor transition is genuinely sharp (threshold) or whether the apparent sharpness is sampling under-density.

4. **Quasispecies with asymmetric base rates** (Session 10). Combine f_A ≠ f_B AND base_A ≠ base_B. Tests whether selection-driven and base-rate-driven off-attractor mechanisms compose, cancel, or stay distinct.

5. **PDG four-force real-data probe** (queued from Session 8/9). Pull PDG running of g1, g2, g3 with energy. Apply per-domain stack. Tests whether real EM/weak/strong behave like the Session 8 toy caricatures. Universe-origin side.

6. **Compartmentalization probe** (Session 10 candidate). Add a boundary (protocell) to the autocatalytic system. Different physics from simple ODE — may move operator-space signature off in a distinct way.

7. **Tools Greg attaches in next session** (carried from v10 kickoff queue).

8. **AF cohort full fetch + plug into store** (carried from v10 kickoff).

9. **Active-learning hook in od_learning_store** (carried).

10-14. Carried from prior session queues: online classifier swap, cross-problem transfer, storm/waves real-data, MIMIC credentialing, five-physics-areas substrate hunt, complexity-entropy plane (Rosso), Yamada 2023 multi-field SVD, Gomez-Herrero 2015 ensemble infra, Jacobson 1995 cross-domain, Stochastic Electrodynamics 21st-century.

## Architecture (current)

- **DeepNova**: 92 passing tests, 22 manifests, persistent evidence graphs, PPO retrieval policy learner. F:\Factory\.
- **VOXA**: voice interface layer. Cloud-hosted TTS MCP server.
- **Agent Factory**: 23 agents, 5 divisions, F:\Factory\.
- **Token Optimizer**: deployed at optimizer.davisai.ai, Stripe live.
- **OD provisional patents**: 3 filed March 24, 2026.

## Defense Pipeline (status as of last update)

- DARPA Bio Attribution (confirmed top-10), CyPhER Forge, TTO BAA, DIU PRISM, CIA (KV3UCQ1A submitted), IQT, MDA MAA, AFWERX.
- **Steve "Bucky" Butow** (DIU Space): personal email contact.
- **Carl Saab** (Cleveland Clinic): OD outreach engaged.
- **Roland Rott** (GE HealthCare Imaging): MRI proof-of-concept brief delivered.
- SAM.gov UEI: CQ56XYFZL4E6. CAGE pending.

## Standing Decisions

- LlamaIndex: declined. Robyn: declined. Nous Atropos: worth evaluating. Hermes 4 14B: recommended for local reasoning on SENTINEL. HomeLift: dormant.

## Session Handoff Pointer

For active Session 10 work (prebiotic line + new "Time as expression of dipole flow" frame), read `SESSION_HANDOFF_2026-05-27_v10.md` first.

For Session 9 OD platform + real-data probe context, read `SESSION_HANDOFF_2026-05-26_v9.md` and `OD_LEARNING_README.md`.

For Sessions 5-8 simulator-based information-layer context, read v5-v8 handoffs in order.

## Note (Session 10 update — 2026-05-27)

Updated at the end of 2026-05-27 Session 10 to reflect:

- Greg's session-open directive: "Do you think we can figure out how the universe came into existence, the process and how life started? ... let's do it." Session 10 opened the prebiotic-chemistry-to-biology line via the per-domain stack on simulator systems with explicit biology-defining ingredients (autocatalysis, mutation, selection).
- Branch: `claude/unused-tools-review-r3CJV`. Session 9's `claude/two-more-tasks-iZvY4` branch pulled in via `git checkout origin/...` of the per-domain stack scripts (kbk_pipeline.py, ai_poincare_rank.py, sindy_symbolic.py, per_domain_kbk.py, pysr_symbolic_per_domain.py). Platform code (od_*.py) NOT pulled in (not needed for simulator probes).
- **Five experiments run this session**:
  1. Prebiotic canary (1 seed, N_ens=300, T=20): k_cat=0.25 off-attractor point noted but not interpretable at one seed.
  2. Multi-seed replication (3 seeds, Session 7 baseline): canary's off-attractor point WASHED OUT at all k_cat. Rule D. INFO-037 (candidate).
  3. Quasispecies probe (3 seeds, mu in [0, 0.5]): first reproducibly-off-attractor OD system. Cos walks 0.056 → 1.000 smoothly. Tracks |H_a - H_b|. INFO-038 (candidate).
  4. Asymmetric hypercycle control (3 seeds, base_A in [1.0, 3.0], NO selection): also goes off attractor. Apparently sharp transition. Channel asymmetry per se is sufficient — reading (c) confirmed. Mechanism distinction (smooth vs sharp) remains a candidate. INFO-039 (candidate).
  5. PySR cross-seed test on both systems: no reproducible per-system family. Session 7 INFO-025 four-family classification doesn't extend. INFO-028 reproduced. INFO-040 (methodological candidate).
- **Substantive close**: Greg surfaced the "Time as expression of dipole flow" frame conversationally near session end. Frame folded into Working Frames above. Two dedicated probes queued (time-reversal, entropy-production correlation) — frame stays speculative until probes run.
- **All seven Operating Rules** from Sessions 4-7 in force throughout. NO new Operating Rule this session.
- **Universe-origin side**: untouched. PDG four-force real-data probe still queued.
- Branch state: work persisted on `claude/unused-tools-review-r3CJV`, pushed. main untouched. No PR.
