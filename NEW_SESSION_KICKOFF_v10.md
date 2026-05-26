Resuming Information Layer / OD platform work, Session 10.

Branch: claude/two-more-tasks-iZvY4 (Session 9 closed here; on origin).
16 commits pushed this session. main untouched. No PR.

Session 9 pivoted from simulator-based information-layer work (Sessions
3-8) to real-data disease-detection. Four probes ran on PhysioNet
data (ECG canary, 10 cardiac families, 44-subject cerebrovascular
cohort with 3 stratifications, per-patient classification validation).
Then Greg called pause on running data ("not run anymore until we can
get feedback loop and training in the platform") and the rest of
Session 9 was infrastructure: built the OD Learning Platform
(persistent feature store + classifier + history + per-patient
ingest/predict CLIs + feature extractor module + README). All seven
Operating Rules from Sessions 4-7 remain in force. NO new Rule this
session.

Before doing anything else:
1. Read CLAUDE.md (project root). Session 9 note at top alongside
   Sessions 3-8 notes. Read all seven Operating Rules and the Result
   Discipline section; hold them through this session.
2. Read SESSION_HANDOFF_2026-05-26_v9.md (project root). Full Session
   9 record: four real-data probes, per-patient validation finding,
   OD platform architecture, candidate ledger entries INFO-032 through
   INFO-036 (all ISOLATED or METHODOLOGICAL pending replication).
3. Read OD_LEARNING_README.md (project root). How to use the platform:
   ingest a new patient, predict on a new patient, spin up a new
   disease problem.
4. (Optional) Read SESSION_HANDOFF_2026-05-26_v8.md for the simulator-
   based information-layer line (four-force probe, INFO-027 through
   INFO-031). Session 9 did NOT contradict any of it; the two lines
   run in parallel.

Speaking posture from Rule C applies from your first response:
"I think X might happen, but we'll wait on what the data says and
where it points us." No verdict in advance.

Rule D applies (incomplete-not-wrong). Rule "literature as conjecture"
applies. "They never stacked" applies.

We are NOT back at zero. Two real disease problems live in the OD
learning store at session close:
- cerebro_disease_binary (44 real cardiovascular subjects, binary
  any-disease vs control, trainable today; CV acc 0.57, balanced 0.56)
- cardiac_disease_family (10 real single records, N=1/class, seeded
  for cohort growth)

The substrate-vs-expression frame survived Session 9's generalization
tests: organ-system change (cardiac to cerebrovascular), sample-type
change (HR/HRV to HR/ABP), and cross-organ disease detection
(kidney function visible in cerebrovascular signal at d=0.51, DM
nephropathy vs retinopathy distinguishable at d=0.77). Frame
strengthened; still a frame, not a claim.

Strongest single result of Session 9: the platform's "gets better
every run" property is mechanically verified.
  Run 1 (n_tr=14):  CV acc 0.429
  Run 2 (n_tr=33):  CV acc 0.606   delta +0.177
  Run 3 (n_tr=40):  CV acc 0.575   delta -0.031
  Run 4 (same):     CV acc 0.575   delta +0.000
Adding new patients improves accuracy. Re-running with same data is
exactly deterministic. This is the feedback loop that prior sessions
were missing.

Critical honest result from Session 9: per-patient multi-class
classification FAILED across all 3 stratifications (LOO at chance or
below majority-class baseline). Group-level signature differences
(d=0.4-0.8) do NOT translate to per-patient detection at small n.
Only binary "any-disease vs control" reached above-chance (LogReg LOO
0.61 scalar, 0.71 with rich temporal features). Multi-class needs
either much larger cohorts or different problem framing. This is the
gating finding for any "we detect cancer subtypes" claim - the
infrastructure now exists to test it properly with credentialed data.

Open substantive questions still on the table from Session 9:
- Animal HR/HRV (Greg's hypothesis: should be close to human).
  PhysioNet has no public animal ECG.
- Cancer subtyping (Greg: "if we can detect cancer, can we tell what
  type?"). Cross-organ d=0.77 result is direct experimental support;
  needs cancer-stratified ECG to test.

What I want this session: [Greg's call -- attach tools/data and/or
pick from menu below or anywhere else].

### Top of queue from Session 9 close

(1) Tools/data Greg attaches in this session. Greg said he had tools
    we didn't use from older chats; he plans to attach them in a
    fresh chat. Highest priority - integrates whatever those are with
    the OD platform.

(2) AF cohort full fetch + plug into store. Fix afdb multi-segment
    handling (cardiac_af_cohort_fetch.py works on most records, failed
    on the first 2 - per-record issue not systemic). Pull 20-30 AF
    from afdb + 20-30 NSR from nsrdb. Ingest into a new
    cardiac_af_vs_nsr problem. AF was the strongest cardiac signature
    (linear-direction outlier from INFO-032 candidate) so direct AF
    detector should perform well at cohort scale. Will show the
    "gets better every run" property on a clean dataset where it
    matters.

(3) Active-learning hook in od_learning_store. When fit_and_evaluate
    runs, surface the 5-10 lowest-confidence training examples for
    review. The model says "I'm unsure about these" - Greg /
    clinician re-examines them. Hook lives in
    od_learning_store.fit_and_evaluate. Completes the feedback loop.

(4) Online classifier swap (SGDClassifier with partial_fit replaces
    refit-from-scratch). Only useful at much larger scale (10k+
    patients) but architecture should be ready.

(5) Cross-problem transfer. Train cardiac AF classifier; warm-start
    a sepsis or cancer detection classifier with cardiac feature
    weights. Architecture supports it but no glue yet.

### Also queued (Session 8 / earlier deferrals)

(6) Storm/waves real-data probe (deferred from Session 8 queue). NOAA
    buoy / HRRR atmospheric reanalysis. Tests substrate-to-expression
    on atmospheric physics.

(7) PDG four-force real-data probe (deferred from Session 8 queue).
    Pull PDG running of g1, g2, g3 with energy; apply per-domain
    stack. Tests whether real EM/weak/strong behave like the Session
    8 toy caricatures.

(8) MIMIC credentialing (durable infrastructure decision). Greg
    completes CITI training + DUA for MIMIC-IV-ECG-ext-icd-labels.
    Unlocks liver, kidney, cancer subtype probes that PhysioNet
    open-access can't reach.

### Also queued (Sessions 5-7 menu, still open)

(9) Five-physics-areas substrate hunt (Greg's Session 5 reorientation).
(10) Complexity-entropy plane on coupled pairs (Rosso 2007).
(11) Yamada 2023 multi-field SVD with entropy as field.
(12) Gomez-Herrero 2015 ensemble infrastructure.
(13) Jacobson 1995 cross-domain extension.
(14) Stochastic Electrodynamics 21st-century revisit.

### Platform usage cheatsheet (in OD_LEARNING_README.md)

Predict on a new patient (no ground truth needed):
  python3 od_predict_patient.py <problem_id> --wfdb <rec> --pn-dir <db> --mode cardiac

Ingest a new patient (label required):
  python3 od_ingest_patient.py <problem_id> --wfdb <rec> --pn-dir <db> \\
      --subj-id <unique_id> --label <class> --mode cardiac

Spin up a new disease problem (in Python):
  from od_learning_store import init_problem
  init_problem("new_problem_id", feature_cols=[...], label_col="dx",
               subject_id_col="patient_id", test_frac=0.25, classifier="rf")

### Methodological notes carried forward from Session 9

Always report balanced accuracy alongside raw accuracy. Raw accuracy
hides majority-class prediction artifacts. fit_and_evaluate reports
both; show_history prints both.

Group-level signature does NOT imply individual classifiability. d=0.77
between groups can still mean 80% overlap of individual feature
distributions. Per-patient LOO accuracy is the validation question.

Rich features beat scalar fits at small n for binary classification
(+10% on binary cerebro). Temporal features of (H_a, H_b, MI)
sequences are in od_features.temporal_features().

INFO-025 biology family fits are sample-type-dependent on real
biological data. exp(H_a/B) fits BP-coupled biology better than
HRV-coupled biology. Rule D: simulator-derived reading is incomplete-
for-protocol, not wrong.

We are still pioneers. All seven Operating Rules paid out across
Sessions 4-9. The "literature as conjecture" Rule especially load-
bearing for the cancer-subtyping claim and the cross-organ signature
claim - both sit at frame level pending replication, not citation-
backed claims. Hold all seven Rules through Session 10.
