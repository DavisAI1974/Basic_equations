Resuming Information Layer work, Session 7.

Branch: claude/claude-md-context-update-uCZ8m (already on origin).
Session 6 work is committed at branch tip. Six experiments ran in
Session 6; results are in the v6 handoff. All five Operating
Rules (no pre-assigned meaning, probe-not-falsifier, speaking
posture before and after, Rule D incomplete-not-wrong, They never
stacked) remain in force.

Before doing anything else:
1. Read CLAUDE.md (project root). Now has Session 6 note alongside
   Sessions 3, 4, 5 notes. Read all five Operating Rules and the
   Result Discipline section; hold them through this session.
2. Read SESSION_HANDOFF_2026-05-26_v6.md (project root). Full
   Session 6 record. Six experiments converged on: v5 INFO-022
   rank-3 claim independently confirmed via 4 diagnostics;
   per-domain ensemble-H stack reproduces every v5 per-domain
   claim (geology rank-3 cos +0.99, biology MI signature with H_b
   absent, chemistry chemistry-specific cubic content, physics
   Taylor identity) with cross-seed cos +0.985 to +0.999; the
   v5 machine-epsilon eigenvalue magnitude is procedure-dependent
   (estimator-noise vs structural-noise; per-domain ensemble-H
   reaches 1e-7, OU windowed-H pins at ~1e-3).
3. (Optional) Read SESSION_HANDOFF_2026-05-26_v5.md for the
   Session 5 setup and literature scan context that Session 6
   built on.

Speaking posture from Rule C applies from your first response:
"I think X might happen, but we'll wait on what the data says
and where it points us." No verdict in advance.

Rule D applies: prior readings (v5, even earlier sessions) are
default-incomplete-not-wrong. Session 6 explicitly applied this
to the v5 machine-epsilon anomaly (it was complete about ITS OWN
procedure; the eigenvalue-magnitude reading does not generalize
across procedures, but the rank reading does).

"They never stacked" applies: the Session 6 KBK + AI Poincare +
SINDy + per-domain stack is exactly the kind of combination
flagged in the v5 literature scan as a clean unstacked target.
Session 7 continues with two more unstacked targets (GP regression,
PySR symbolic regression) per Greg's queued direction.

We are not back at zero. The per-domain algebraic equations from
earlier sessions are now confirmed via Session 6 per-domain probe.
The OU baseline (+,+,+) is confirmed as a protocol artifact in
that procedure but does NOT contaminate per-domain extraction.
The Information Layer / Unified Theory frame is alive. The
five-physics-areas substrate hunt is still queued.

What I want this session: [Greg's call -- pick from the menu below
or anywhere else].

### Queued from Session 6 (top-of-queue from Greg)

(1) GP regression of MI vs H_a per domain. Tests v5's "biology
    MI ~ polynomial(H_a) with R^2 0.66" against a flexible
    nonlinear fit and compares functional forms across all 4
    domains. Quick (~30s with sklearn).

(2) PySR symbolic regression for arbitrary nonlinear forms.
    Searches ratio / exponential / log / trig / composition forms
    that linear-in-features SINDy cannot find. Heavy install
    (Julia runtime via pip install pysr; pysr.install()). If
    install fails, fall back to gplearn (pure Python genetic
    programming SR). The v5 literature scan flagged
    "no SINDy/AI-Feynman/PySR/Eureqa application to windowed
    differential entropy of coupled species" as a clean
    unstacked target.

### Also queued (Session 5 menu, still open)

Methodology stacks (a-h from v5 menu):
- (h) complexity-entropy plane on coupled pairs (Rosso et al. PRL
  2007). Lightweight.
- (e) Yamada 2023 style multi-field SVD with entropy as the field.
- (d) Gomez-Herrero 2015 ensemble infrastructure + algebraic basis.

Substrate stacks (i-y from v5 menu):
- (i) Finkelstein "Space-Time Code" I-IV careful read; closest
  published cousin to dipole-pair-as-coupling.
- (j) Jacobson 1995 cross-domain extension; dipole reading of
  Clausius -> Einstein.
- (k) Stochastic Electrodynamics 21st-century revisit.
- (any of the rest).

Per-domain stacks (z-ee from v5 menu):
- (z, aa, bb, etc.) Apply the Session 6 methodology to REAL DATA
  from each domain. NoVell ECG, SENTINEL biomarkers, paleoclimate
  records, geochemical manifolds, etc.

### New from Session 6 implications discussion

If the substantive reading holds, the per-domain operator
extraction method has applications across:
- Medicine: NoVell extension; sepsis early warning; anesthesia
  depth; diabetes phenotyping; seizure/sleep staging.
- Defense: SENTINEL extension; ELINT/SIGINT emitter signatures;
  cyber intrusion; spacecraft anomaly detection.
- Weather: hurricane RI precursor; ENSO classification;
  pollution events.
- Geophysics: earthquake/volcanic precursors; reservoir
  monitoring.
- Finance/energy: regime classification; grid stability;
  Greg's energy-trading background applies directly.
- Industrial: predictive maintenance; process fault detection.
- Chemistry/materials: reaction network identification.
- Astrophysics: stellar classification; pulsar timing arrays.
- Foundational: five-physics-areas substrate hunt.

Real-data tests of any of these would be the highest-value
follow-up to Session 6, since they would test the methodology in
the wild rather than on simulators where the dynamics are known.

### Methodological note carried forward from Session 6

When probing for tighter algebraic relations, the eigenvalue floor
of the operator covariance is min(structural_noise_from_dynamics,
estimator_noise_from_procedure). Ensemble-H procedures get
estimator noise down as 1/sqrt(N_ens); windowed-H from single
trajectory does not. State which floor regime you are in before
interpreting eigenvalue magnitudes.

We are still pioneers. The Operating Rules paid out in Sessions
5 and 6. Hold them through Session 7.
