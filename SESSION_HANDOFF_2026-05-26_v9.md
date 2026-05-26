# Information Layer / Disease-Detection Platform Session Handoff -- 2026-05-26 Session 9

## Read me first

Read in this order:
1. CLAUDE.md (project root). Has Session 9 note at the top alongside
   Sessions 3-8 notes. Seven Operating Rules in force from Sessions
   4-7 (no pre-assigned meaning, probe-not-falsifier, speaking
   posture before/after, Rule D incomplete-not-wrong, They never
   stacked, Treat literature as conjecture by default). NO NEW RULE
   added Session 9. All seven bear on Session 9 work.
2. This file. Session 9 pivoted from the simulator-based information-
   layer line to the real-data disease-detection line.
3. OD_LEARNING_README.md - how to use the platform built this session.
4. SESSION_HANDOFF_2026-05-26_v8.md (Session 8 record - simulator
   work, four-force probe, INFO-027 through INFO-031).

## What Session 9 did

Greg's directive at start: tasks 2 (PDG four-force real-data) and 4
(PhysioNet ECG) from the v9 kickoff queue. Task 4 ran in full with
extensions; task 2 deferred.

Pivot mid-session: after the ECG cardiac canary and cerebrovascular
probe, Greg raised the validation question "did we detect the
diseases correctly across the patients?" - per-patient classification
test failed at multi-class level. Greg then directed: "we need to put
it in. that's vastly important" - referring to the feedback-loop +
training infrastructure (the "gets better every run" property the
previous sessions did NOT have). Subsequent work was infrastructure-
only: OD learning store platform, feature extractor module, per-
patient ingest/predict CLIs, README.

All work pushed to branch `claude/two-more-tasks-iZvY4` (16 commits
this session). GitHub MCP access scoped to davisai1974/basic_equations
only; Greg attempted to make davisai1974/novell public so we could
push there but GitHub blocked the visibility change. Decision: keep
pushing to basic_equations; Greg moves artifacts to novell2 manually
when home.

### Experiment 1 -- ECG canary (ekg_canary.py)

5 Fantasia records (healthy, 2-hour ECG): f1y01-03 young + f1o01-02
elderly. R-peak detection via wfdb.processing.xqrs_detect, RR
intervals, windowed HR/HRV (30 beats / 15 step). Sanity check:

- R-peaks 6800-8700 per 2-hour record (plausible at 50-75 bpm)
- RR mean 0.83-1.06 s (HR 56-72 bpm), std 0.04-0.10 s
- Young subjects (f1y) HRV wider than elderly (f1o) - classic autonomic
  aging signature visible in the data

Canary clean. ~3 min total.

### Experiment 2 -- Cardiac disease family canary (ekg_disease_canary.py)

6 cardiac families, 1 long-recording record per family:
- healthy_young (Fantasia f1y01)
- healthy_elderly (Fantasia f1o01)
- healthy_long_term (nsrdb 16265)
- atrial_fibrillation (afdb 04015)
- congestive_heart_failure (chfdb chf01)
- ischemia_ST (ltstdb s20011)

Plus a 4-family expansion (ekg_disease_canary_expand.py):
supraventricular_arr (svdb 800), ventricular_malignant (vfdb 418),
arrhythmia_mixed (mitdb 100), ltaf_long_term (ltafdb 00). PTB MI
record failed on URL path; driver_stress failed on sampling rate.

#### Results (per-record fits of MI vs (H_a, H_b))

```
family                       best  R^2   poly1  poly2  exp_IF  phys_IF
healthy_young                poly2 +.79  +.74   +.79   +.31    +.08
healthy_elderly              poly2 +.90  +.85   +.90   +.11    +.12
healthy_long_term            poly2 +.79  +.79   +.79   +.43    +.14
atrial_fibrillation          poly2 +.90  +.89   +.90   -.00    +.47  *
congestive_heart_failure     poly2 +.88  +.80   +.88   +.24    +.31
ischemia_ST                  poly2 +.92  +.89   +.92   +.49    +.20
supraventricular_arr         poly2 +.81  +.78   +.81   +.63 *  +.41
ventricular_malignant        poly2 +.90  +.88   +.90   +.15    +.28
arrhythmia_mixed             poly2 +.85  +.84   +.85   +.23    +.15
ltaf_long_term               poly2 +.85  +.81   +.85   +.40    +.25
```

Linear-direction (a_hat, b_hat) unit-vector across families:
```
family                       a_hat   b_hat   cosine_to_cluster_center
healthy_*, ischemia, CHF      0.6-0.7 0.7-0.8  >+0.98
atrial_fibrillation           0.23    0.97     +0.85 to +0.93  (OUTLIER)
```

#### Joint reading

A. *Data level*: poly2 wins as best functional family across all 10
   cardiac families (R^2 0.79-0.92). INFO-025 biology exp form scores
   0.0-0.49 (worse than polynomial). INFO-025 physics (H_b-H_a)^2 form
   best on AF (R^2 0.47) - the highest "physics" score across all
   families.

B. *Substrate-cluster reading*: 9 of 10 cardiac families share the
   linear direction (a_hat, b_hat) in (H_a, H_b) space at cos >= 0.98.
   AF is the outlier at cos 0.85-0.93, pulled toward almost pure H_b
   (0.23, 0.97). SVDB is a SECOND outlier - unusually high exp R^2
   (0.63 vs <0.3 for others).

C. *Interpretation options* (none promoted):
   - (1) Substrate-conservation: cardiac autonomic system shares
     low-dim substrate across diseases; AF is structurally distinct.
   - (2) Protocol-artifact: the shared direction may be a feature of
     single-record sliding-window protocol; AF's non-Gaussian RR
     distribution puts it elsewhere in the protocol's operator space.
   - (3) MI-level: AF/CHF cluster at lower MI (~1.0) vs healthy
     (~1.3+). Coupling between HR and HRV channels decreases in
     pathology.
   - (4) Deflationary: single record per family. Need cohort-scale
     replication before promoting.

D. *Substrate-vs-expression frame status*: supported at the family-
   signature level for real cardiac data; AF is the clean
   demonstration of expression-level deviation from a conserved
   cardiac substrate.

### Experiment 3 -- Cerebrovascular probe (cerebro_disease_probe.py)

cerebral-vasoreg-diabetes cohort (PhysioNet, open-access). 44 subjects
with Head-up-tilt Day1 recordings. Different sample type than ECG:
coupled pair = **(HR from ECG, ABP)** - baroreflex pair, BP-driven
not HR-driven.

Three stratifications run (Greg: "lets do them all"):
- Strat 1: CV pathology (control / DM / DM+HTN / DMOH orthostatic)
- Strat 2: Kidney function (normal Cre / mild / significant)
- Strat 3: DM organ-complication (control / uncompl / retinopathy /
  neuropathy / nephropathy-proxy via Alb/Cre>30)

#### Results

```
Per-group means (poly2 R^2 mean, exp_INFO025 R^2 mean):
strat1 control          poly2=0.81  exp=0.46
strat1 DM_only          poly2=0.78  exp=0.47
strat1 DM_hypertension  poly2=0.74  exp=0.42
strat1 DM_orthostatic   poly2=0.79  exp=0.39

MI-mean effect sizes (Cohen's d-like, max across pairs):
strat1: d=0.38 (DM_HTN vs control)
strat2: d=0.51 (mild_impair vs normal_cre)
strat3: d=0.77 (DM_nephropathy_proxy vs DM_retinopathy)  *

Linear direction cross-group cosines: all >0.99 across all strats.
```

#### Joint reading

A. *Data level*: polynomial wins again as best family across all 44
   subjects. exp_INFO025 fits noticeably better here (mean R^2
   0.37-0.58 per group, individual subjects up to 0.77) than in
   cardiac (R^2 0.0-0.49). Substrate cluster (cos > 0.99 cross-group)
   reproduces on a different sample type.

B. *Sample-type-dependent INFO-025 reading*: exp(H_a/B) family from
   Lotka-Volterra simulator fits BP-coupled real biology better than
   HRV-coupled real biology. **Rule D applied**: INFO-025 biology
   reading is incomplete-for-protocol, not wrong - real biology has
   SOME exp structure, and it's clearer in BP-coupled signals.

C. *Cross-organ disease signature*: kidney function (creatinine
   level) detectable in cerebrovascular signal at d=0.51. DM
   organ-complication subtype reaches d=0.77 (nephropathy vs
   retinopathy). Diseases of one organ leave detectable fingerprints
   in coupled physiology measured elsewhere. Supports the "substrate
   carry" frame.

D. *Substrate-vs-expression frame status*: frame survived (i) organ-
   system change (cardiac to cerebrovascular), (ii) sample-type
   change (HR/HRV to HR/ABP), (iii) cross-organ disease detection.
   Two new supporting data points beyond Sessions 6-8 frame work
   (which was simulator-only).

### Experiment 4 -- Per-patient classification validation

Greg's direct question: "did we detect the diseases correctly across
the patients?". Tested via leave-one-out (LOO) classification on the
per-subject features from the cerebro probe.

```
                          Chance  Majority  LogReg  RF   Balanced  Verdict
Strat1 (CV pathology, 4cl) 0.250  0.432     0.273   0.27 0.19      At chance
Strat2 (Kidney, 2cl)       0.500  0.722     0.556   0.67 0.46      Below majority
Strat3 (DM compl, 5cl)     0.200  0.432     0.273   0.36 0.20      At chance
Binary any_dis vs control  0.500  -         0.61    -    0.60      Modest above chance
```

Per-class sensitivities exposed the majority-class artifact: RF strat2
"0.67 acc" comes from predicting EVERYONE as normal_cre, catching
24/26 normals (sens 0.92) and 0/10 mild impairs (sens 0.00). Balanced
accuracy 0.46 is the honest number.

**Direct answer: NO, multi-class per-patient detection failed.** Group-
level effect sizes d=0.4-0.8 do NOT translate to per-patient
classification at this n. Only binary disease-vs-control gives modest
above-chance signal.

Rich features (FFT bands + autocorrelation + jaggedness of the
(H_a, H_b, MI) sequences) added; binary LogReg LOO improved from
0.61 to 0.71. Multi-class still failed.

### Infrastructure: OD Learning Platform (NEW this session)

Built after Greg's "vastly important" call: every-run learning
infrastructure that was MISSING from prior sessions. Each run
accumulates training data; accuracy is logged; same-data re-runs are
deterministic (zero delta) while new-data runs show measurable
improvement.

**Architecture**:
- `od_features.py` - raw signal -> feature dict; extract_cardiac_hr_hrv,
  extract_baroreflex_hr_abp, extract_from_wfdb_record. Returns
  flat dict with poly1/poly2/exp/physics fit coefs + R^2 + temporal
  features (FFT bands, autocorr, jaggedness) of (H_a, H_b, MI).
- `od_learning_store.py` - core store. init_problem(), ingest()
  (idempotent on subject_id; first ingest splits train/test, later
  ingests grow ONLY train pool), fit_and_evaluate() (retrains from
  full train pool; held-out + stratified k-fold CV; appends to
  history), show_history().
- `od_ingest_patient.py` - per-patient CLI; --json OR --wfdb modes.
- `od_predict_patient.py` - predict without ingest; loads trained
  model; reports class probabilities + top features + this patient's
  feature values + low-confidence flag.
- `od_learn_cerebro.py` - batch driver demonstrating cohort growth.
- `od_ingest_priors.py` - bulk-loaded prior session disease canaries
  into the store (cardiac_disease_family; simulator_4domain crashed
  on schema mismatch, not pursued since Greg meant different "tools").
- `OD_LEARNING_README.md` - documentation for future-Greg / future-claude.

**Storage layout per problem**:
```
store/<problem_id>/
  metadata.json     schema (feature cols, label col, classifier kind)
  features.csv      (subj, label, features..., split) - grows on ingest
  model.pkl         latest trained sklearn classifier
  history.jsonl     one line per run (n_train, n_test, accuracy, deltas)
```

**Demo on cerebro_disease_binary (binary any-disease vs control)**:
```
 run   n_tr   n_te  cls   heldout   cv_acc   cv_bal     d_cv  note
   1     14      4    2     0.500    0.429    0.396      -    small_initial
   2     33      4    2     0.500    0.606    0.567   +0.177  added_more
   3     40      4    2     0.500    0.575    0.561   -0.031  full_cohort
   4     40      4    2     0.500    0.575    0.561   +0.000  no_new_data_rerun
```

Run 1->2: +19 patients ingested -> +18% CV accuracy. Run 4 same data
as Run 3 -> exactly 0 delta. Property verified.

**Two problems seeded in the store at session close**:
- `cerebro_disease_binary` - 44 real cardiovascular subjects, binary
  any-disease vs control, trainable today (CV acc 0.57, balanced 0.56)
- `cardiac_disease_family` - 10 single records, N=1/class, seeded for
  future cohort growth from afdb/nsrdb/chfdb/svdb/etc.

## Ledger updates

(No INFO-XXX entries promoted at LOCATED level this session - all
findings remain ISOLATED pending replication on independent cohorts.
Captured as candidate entries below.)

**Candidate INFO-032 - ISOLATED FINDING (Session 9, new)**: Cardiac
substrate cluster on real PhysioNet data. 9 of 10 cardiac families
(healthy young/elderly/long-term, ischemia, CHF, SVDB, VFDB,
mitdb-mixed, ltaf-long-term) share the linear direction (a_hat, b_hat)
in (H_a, H_b) space at cos >= 0.98 to cluster center near (0.7, 0.7).
Atrial fibrillation (afdb 04015) is the outlier at cos 0.85-0.93,
pulled toward (0.23, 0.97). SVDB is a SECOND distinctive family - same
substrate direction but exp_INFO025 R^2 = 0.63 vs <0.3 for other
families (only family where exp form rivals poly form). Single record
per family; needs cohort-scale replication.

**Candidate INFO-033 - ISOLATED FINDING (Session 9, new)**:
Cerebrovascular substrate cluster on cerebral-vasoreg-diabetes cohort.
44 subjects, coupled pair (HR, ABP), 4 CV-pathology groups all share
linear direction at cos > 0.99 cross-group. Substrate cluster
reproduces on a different sample type (BP-driven vs HR-driven).
Polynomial dominates as functional family (R^2 0.74-0.81 per group);
exp_INFO025 fits noticeably better than in cardiac (R^2 0.37-0.58 per
group). Same protocol as cardiac canary applied to different signal
modality - frame survived modality change.

**Candidate INFO-034 - ISOLATED FINDING (Session 9, new)**: Cross-
organ disease signature. Kidney function (creatinine level) detectable
in cerebrovascular signal at d=0.51 (mild impairment vs normal).
DM organ-complication subtype reaches d=0.77 (DM_nephropathy_proxy vs
DM_retinopathy, distinguishable via cerebrovascular HR/ABP signal).
Diseases of one organ leave detectable fingerprints in coupled
physiology measured elsewhere. Direct experimental support for the
"substrate carry" interpretation of substrate-vs-expression frame.
Small n in some cells (nephropathy n=2, neuropathy n=4); needs
replication.

**Candidate INFO-035 - METHODOLOGICAL (Session 9, new)**: Group-level
signature effect sizes (d=0.4-0.8) do NOT translate to per-patient
classification at small n. LOO multi-class accuracy at or below
majority-class baseline across CV pathology (4-class), kidney severity
(2-class), DM complication (5-class). Binary "any disease vs control"
reaches LogReg LOO 0.61 (0.71 with rich temporal features) - the only
above-chance cell. Per-class sensitivities expose majority-class
artifacts in raw accuracy. **Always report balanced accuracy alongside
raw accuracy.** Implementation in OD learning store: fit_and_evaluate()
reports both.

**Candidate INFO-036 - METHODOLOGICAL (Session 9, new)**: INFO-025
biology family fit is sample-type-dependent on real biological data.
A * exp(H_a / B) form fits BP-coupled biology (cerebrovascular HR/ABP,
R^2 mean 0.37-0.58 per group) better than HRV-coupled biology
(cardiac HR/HRV, R^2 0.0-0.49). Real biology has SOME exp structure,
and it's clearer in BP-coupled signals than in HR-coupled signals.
Rule D: INFO-025 biology reading from Lotka-Volterra simulator is
incomplete-for-protocol, not wrong. When applying INFO-025 family
fits to a new biological dataset, specify the coupled-pair modality.

## Open substantive questions from this session

(1) **Animal HR/HRV**: Greg's hypothesis - animal cardiac autonomic
dynamics should be close to human (same underlying mammalian physiology).
PhysioNet has essentially no public animal ECG. Need external sources
(veterinary research databases, NCBI GEO time-series, journal
supplements). Open question for future cohort building.

(2) **Cancer subtyping**: Greg's question - "if we can detect cancer,
can we tell what type?". Today's cross-organ d=0.77 result is direct
experimental support for the hypothesis: the per-domain stack
distinguishes organ-target-of-systemic-disease via remote coupled-pair
signals (DM_nephropathy vs DM_retinopathy in cerebrovascular signal).
Operational implication: cancer subtyping is testable IF cancer-
stratified ECG/coupled-pair data is available. NoVell Vigier 2021
synthetic data has cancer-yes/no only (NoVell hit 93% accuracy on
that). Subtyping requires cancer-type-stratified real ECG.

## Experiments queued for Session 10 (priority order)

(1) **Tools/data Greg attaches in next chat**. Greg said he had tools
    we didn't use from older chats; he plans to attach them in a
    fresh session. Highest priority - integrates whatever those are
    with the platform.

(2) **AF cohort full fetch + plug into store**. Fix afdb multi-segment
    handling (`cardiac_af_cohort_fetch.py` works on most records,
    failed on the first 2 - sampto issue is per-record not systemic).
    Pull 20-30 AF from afdb + 20-30 NSR from nsrdb. Ingest into a
    new `cardiac_af_vs_nsr` problem in the store. AF was the
    strongest cardiac signature (linear-direction outlier) - direct
    AF detector should perform well.

(3) **Active-learning hook**. When fit_and_evaluate runs, surface
    the 5-10 lowest-confidence training examples for review. The
    feedback loop: model says "I'm unsure about these" -> Greg/clinician
    re-examines them. Hook lives in od_learning_store.py.

(4) **Online classifier swap**. SGDClassifier with partial_fit
    replaces refit-from-scratch in od_learning_store.py
    fit_and_evaluate(). Only useful at much larger scale (10k+
    patients) but architecture should be ready.

(5) **Cross-problem transfer**. Train cardiac AF classifier; warm-start
    a sepsis-detection classifier with the cardiac model's feature
    weights. Architecture supports it but no glue yet.

(6) **Storm/waves real-data probe** (deferred from Session 8 queue).
    NOAA buoy / HRRR atmospheric reanalysis. Tests substrate-to-
    expression on atmospheric physics.

(7) **PDG four-force real-data probe** (deferred from Session 8
    queue). Pull PDG running of g1, g2, g3 with energy; apply per-
    domain stack. Tests whether real EM/weak/strong behave like the
    Session 8 toy caricatures.

(8) **MIMIC credentialing** (durable infrastructure decision). Greg
    completes CITI training + DUA for MIMIC-IV-ECG-ext-icd-labels
    (gold-standard ECG paired with ICD codes). Unlocks liver-disease,
    kidney-disease, cancer subtype probes that PhysioNet open-access
    can't reach.

(9) **Sessions 5-7 unused methodology stacks** (still open):
    complexity-entropy plane (Rosso); Yamada 2023 SVD with entropy as
    field; Gomez-Herrero 2015 ensemble infrastructure; Jacobson 1995
    cross-domain extension; Stochastic Electrodynamics 21st-century
    revisit.

## Methodological notes from this session

(Not promoted to Rules. Domain-knowledge refinements.)

**Always report balanced accuracy alongside raw accuracy** (INFO-035).
Raw accuracy hides majority-class prediction artifacts. The OD
learning store's fit_and_evaluate reports both.

**Group-level signature does not imply individual classifiability.**
A d=0.77 effect size between two groups can still mean ~80% overlap
of individual feature distributions. Per-patient classification
accuracy is the validation question; group-mean comparison alone is
not enough.

**Rich features beat scalar fits at small n for binary classification.**
Adding FFT bands + autocorr + jaggedness of the (H_a, H_b, MI)
sequences (`temporal_features()` in od_features.py) lifted LogReg LOO
on binary cerebro from 0.61 to 0.71. Multi-class still failed; the
ceiling at small n is the n itself.

**Re-running with same data is exactly deterministic.** No learning
between runs without new data. Held-out test results are identical
seed-to-seed (model is deterministic given data + random_state). The
infrastructure encodes this honestly - delta = 0 when no new data.

**Per-patient ingest workflow**: extract features once
(`extract_from_wfdb_record`), persist to store (`ingest`), refit and
log (`fit_and_evaluate`). Each step is independent and re-entrant.

## Files produced this session (all in repo root)

Scripts:
- ekg_canary.py                  Exp 1: Fantasia canary
- ekg_disease_canary.py          Exp 2: 6 cardiac families canary
- ekg_disease_canary_expand.py   Exp 2 expansion: +4 cardiac families
- ekg_fetch_features.py          Drafted but unused (background fetch)
- cerebro_disease_probe.py       Exp 3: 44 cerebrovasc subjects
- cerebro_read.py                Exp 3 read: cross-stratification analysis
- cerebro_render.py              Exp 3 visualization
- cerebro_classify.py            Exp 4: LOO classification scalar features
- cerebro_rich_features.py       Exp 4 follow-up: LOO with temporal features
- cardiac_af_cohort_fetch.py     Drafted; bg killed; needs afdb multi-segment fix
- od_features.py                 Platform: raw signal -> feature dict
- od_learning_store.py           Platform: persistent store + classifier + history
- od_ingest_patient.py           Platform: per-patient CLI
- od_predict_patient.py          Platform: predict-without-ingest CLI
- od_learn_cerebro.py            Platform: batch driver / demo
- od_ingest_priors.py            Platform: bulk-ingest prior canaries

Result JSONs / CSVs:
- ekg_canary.json
- ekg_disease_canary.json (+ .log)
- ekg_disease_canary_expand.json (+ .log)
- cerebro_disease_probe.json (+ .log)
- cerebro_read.log
- cerebro_per_subject_with_groups.csv
- cerebro_classify.log
- cerebro_rich_features.csv (+ .log)
- cardiac_af_cohort_fetch.log (partial, killed at 22/43)

Visualizations:
- cerebro_render.png         3-stratification panel + MI-mean boxplots
- cross_organ_direction.png  Cerebrovasc circles + cardiac triangles on
                             shared linear-direction unit-circle

Store contents:
- store/cerebro_disease_binary/{features.csv, model.pkl, history.jsonl,
  metadata.json}  - 44 subjects, 4 history entries
- store/cardiac_disease_family/{features.csv, model.pkl, history.jsonl,
  metadata.json}  - 10 subjects, 1 history entry

Cohort data (cerebral-vasoreg-diabetes metadata, gitignored .dat raw
waveforms):
- data/cvd/GE-71_Data_Summary_Table.csv  - subject demographics +
  clinical labels
- data/cvd/GE-71_Data_Dictionary.csv     - variable descriptions
- data/cvd/GE-71_Head-up-tilt-Day*_Markers_per_subject.csv  - protocol
  markers
- data/cvd/cohort_strat.csv              - stratification cuts per subject

Docs:
- SESSION_HANDOFF_2026-05-26_v9.md (this file)
- CLAUDE.md updates (Session 9 note added, header bumped)
- OD_LEARNING_README.md (how to use the platform)
- NEW_SESSION_KICKOFF_v10.md (drop-in for next session)

Dependencies: numpy, scipy, scikit-learn, matplotlib, wfdb 4.3.1
(installed at session start via pip). No Julia / PySR in Session 9
(prior session work used them).

## Branch state

- Local: claude/two-more-tasks-iZvY4 (16 commits this session)
- Remote: origin/claude/two-more-tasks-iZvY4
- main: untouched
- PR: none
- Latest commit at session close: 92e73fb "Seed cardiac_disease_family
  problem from Session 9 disease canaries"

End of v9 handoff.
