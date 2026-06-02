# SESSION HANDOFF — 2026-06-02 (S16) — gravity / 4-force + CLAUDE workflow

Branch `claude/claude-md-strategy-pwfZr` (main untouched; no PR). All Operating
Rules in force. Scope this session: gravity (S15 thread) + 4 forces (S14 thread).
Greg's standing calls this session: (a) MI-in-null coupling discriminator is
DISPROVED — do not build on it; (b) ignore Markets for now; (c) no data discarded,
odd outputs are kept and diagnosed; (d) treat all literature as conjecture.

## What ran (all committed + pushed, with result JSONs)

1. **Time-blind RECHECK** (`recheck_time_blind.py`) — resolved the INFO-S15a vs
   INFO-046 contradiction on GW150914. SENSE 1 (within-window row shuffle): cos =
   1.00000000, covariance drift 9e-14 => row-order invariance is a MATHEMATICAL
   IDENTITY of the column-covariance extraction, NOT data time-blindness. SENSE 2
   (28-block null trajectory): the null DIRECTION moves strongly in time
   (consecutive |cos| 0.07-0.99), merger block distinct, merger quad-null on
   equal-entropy 0.958 vs noise 0.372. Conclusion: "level null is time-blind" must
   be restated as "the pooled extraction is row-permutation-invariant by
   construction." The TIME-RESOLVED null carries strong time structure.

2. **H-D trajectory PINNED** (`gravity_HD_trajectory.py`) — fixed the resolution
   wall with a stable pooled-noise reference + K=40-row window. Clean
   substrate->DEPART->return: |cos(window-null, noise-null)| 0.82 pre / 0.28 AT
   MERGER / 0.89 post. DEFLATIONARY (load-bearing): pre/post "substrate" IS the
   detector-noise reference, so "return" = return-to-NOISE (INFO-046 gate), not
   ringdown. The DEPARTURE at merger is the real signal. H-D ("structure zeroed to
   substrate") stays a FRAME, ungraded.
   (Superseded the K=16 `gravity_null_trajectory.py`, which was estimator-noise-
   dominated away from the merger — kept on record.)

3. **Gravity chirp Piece-1** (`gravity_chirp_law.py`) — NOT cleanly recovered, and
   recorded as such (no manufactured hit). Tried 5 robust extractions on GW150914
   (n wandered 0.18->0.58->2.0->5.0, R^2 <= 0.46) AND on GW170817 (fetched from
   GWOSC; long BNS inspiral tracked f 30->320 Hz over 6s/18k samples, yet f^-8/3
   linearity R^2 0.001). BOTH events fail => the limiter is the METHOD (per-sample
   whitened-Hilbert instantaneous freq too noisy to track the monotone chirp), NOT
   the event. RULE D correction: my earlier "wrong event, use GW170817" diagnosis
   was INCOMPLETE. Clean recovery needs a Q-transform/matched-filter (backlog #3).
   Qualitative chirp present on both events.

4. **WF Piece-1** (`s16_weak_zpropagator.py`) — CLEAN. Z Breit-Wigner recovered
   from 10227 real CMS dimuon events: M_Z = 90.75 GeV = 99.5% of PDG 91.1876
   (chi2/ndf 3.3-4.2); Gamma_Z over-wide (43-53%) consistent with detector mass-
   resolution broadening (honest instrumental effect). PySR independently recovered
   a BW-like rational lineshape. Same-charge control: no peak. Lit note: SymbolFit
   used this exact data but modeled background; we recovered the resonance.

5. **SF Piece-1** (`s16_strong_running.py`) — CLEAN. QCD running recovered from 13
   real measured alpha_s(Q) world-data points (1.78 GeV-1 TeV), no beta function
   assumed: 1/alpha_s LINEAR in ln(Q), slope +1.31 (POSITIVE => asymptotic freedom
   FORCED by data), b0 8.25, n_f eff ~4.1, Lambda_QCD ~150 MeV, chi2/ndf 0.81.
   PySR recovered the affine 1.24*lnQ+2.57 form. Femtoscopy-R path deferred (backlog #7).

6. **PROBE 1** (`probe1_gate_check.py`, `probe1_rate_dynamics_sim.py`) — gravity
   rate-dynamics gate: merger window ~8 independent samples, can't fit lagged
   dynamics (deferred). Sim half: dMI/dt carries lagged memory in chem (dR2 0.17) >
   phys (0.12) > bio (0.05), NULL in geology (random-CV dR2 negative; INFO-025
   constant-MI). Data level only.

## Literature scan (2 agents; "they never stacked"; all gravity-time refs CONJECTURE-level)
- Gravity-time: Rovelli-Smerlak "temperature as the speed of time" (Tolman-Ehrenfest);
  Castro Ruiz-Brukner "gravity entangles quantum clocks" (PNAS); Smith-Ahmadi
  (arXiv:2304.01263, derives time dilation from a global clock); de Freitas 2024
  (arXiv:2412.12211, Tsallis entropy of real GW150914 strain, single-detector, no MI).
- 4-force: Lemos-Cranmer 2022 (Newton from ephemerides via PySR); Moynihan 2026
  (arXiv:2602.15169, SR rediscovers gravity=gauge^2 double-copy, theory only); PhySO
  units-constrained SR; Hamber-Williams running-of-G; spin-2 universality (MEETS BAR).
- NEGATIVE FINDINGS (= our unstacked lanes): nobody re-derived all 4 force laws from
  raw data in one SR framework; nobody recovered the GW chirp law from real strain by
  SR; nobody recovered the Z Breit-Wigner from real dimuon by SR; "only gravity dilates
  clocks => gravity is the unique time-coupled force" is unwritten as a formal claim;
  inter-detector MI as a physics probe (vs glitch-vetoing) is untouched.
- CAUTION absorbed: MSSM two-loop shrinks the unification triangle — our "no single
  unification (triangle ~1e4)" rests on one-loop; map before leaning on it.

## Operator hunches this session (logged, NOT graded — frames)
- Algebra sharpening (Greg): only a like term cancels a like term, so if time zeroes
  at the horizon, the canceling term must itself be a TIME term inside the gravity
  equation => gravity must CONTAIN a time/flow term (H-A/H-B). Coherent; ungraded.
  Anchor (re-verify): horizon time-freeze is the EXTERNAL/observer-frame statement.

## Data-level analysis: "they all fall into that area" (S16)
Everything measured sits ON the equal-entropy region: 8 toy systems (>=0.99), 4 sim
flow axes (0.94-0.99), 4 real gauge forces (~1.0), LIGO noise (0.98), LIGO merger
no-MI view (0.96). The ONLY departures: LIGO merger in the full 6-op basis (0.23-0.28,
MI spike) and gravity's flow axis (0.61-0.78, loosely pinned). Both departures are
GRAVITY-associated. READING (direction-shift, frame): the universal clustering is the
signal — that region is the universal SUBSTRATE/base (Base-of-Structure; existence
pole), and the gravity-associated departures are where structure leaves the base.
DEFLATIONARY CRUX (= backlog #1): the clustering is partly forced by building
symmetric equal-marginal-entropy channels (INFO-038: asymmetry moves objects off).
Open question that decides nature-vs-bookkeeping: is equal-marginal-entropy a forced
property of physical 2-channel observables or our construction choice?

## State
Branch pushed. Repo CLAUDE.md still canonical only through S11 (master is S15, this
is S16) — the master-merge is backlog #12. See BACKLOG_tests_and_probes.md for the
full queue and the standing "clear backlog before new probes" rule.
