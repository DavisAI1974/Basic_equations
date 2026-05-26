# OD Learning Platform

Persistent per-domain-stack disease-detection framework. Built Session 9.

## What it does

Given a coupled-pair physiological signal (e.g., ECG, ECG+ABP) per
patient, extract per-domain stack features, train a classifier, and
**accumulate training data across runs**. Every new patient grows the
training pool; accuracy is logged so improvement is tracked over time.

Re-running with the same data is exactly deterministic (no learning).
Adding a new patient is where the accuracy delta shows up.

## Architecture

```
od_features.py            Raw signal -> feature dict
                            extract_cardiac_hr_hrv(ecg, fs)
                            extract_baroreflex_hr_abp(ecg, abp, fs)
                            extract_from_wfdb_record(rec, ..., mode=)
od_learning_store.py      Persistent store + classifier
                            init_problem(...)
                            ingest(problem_id, df)
                            fit_and_evaluate(problem_id)
                            show_history(problem_id)
od_ingest_patient.py      CLI: ingest one patient (from JSON or WFDB)
od_predict_patient.py     CLI: predict on a new patient (no ingest)
od_learn_cerebro.py       Batch driver (demo on cerebro_disease_binary)
```

Storage layout per problem:

```
store/<problem_id>/
  metadata.json     schema (feature cols, label col, classifier type)
  features.csv      growing store of (subj, label, features, split)
  model.pkl         latest trained sklearn classifier
  history.jsonl     one line per run (n_train, accuracy, deltas)
```

## Quick start: existing problem (cerebro_disease_binary)

```bash
# Predict on a new patient (no label needed)
python3 od_predict_patient.py cerebro_disease_binary \
    --wfdb s0030DA --local-dir data/cvd/sample --mode baroreflex

# Ingest a new patient (label required)
python3 od_ingest_patient.py cerebro_disease_binary \
    --wfdb s0030DA --local-dir data/cvd/sample --mode baroreflex \
    --subj-id s0030_v2 --label control --note "second visit"

# Re-evaluate without ingesting anything new (same accuracy as last run)
python3 -c "from od_learning_store import fit_and_evaluate, show_history; \
    fit_and_evaluate('cerebro_disease_binary', note='manual_eval'); \
    show_history('cerebro_disease_binary')"
```

## Quick start: NEW problem

```python
from od_learning_store import init_problem, ingest, fit_and_evaluate
import pandas as pd

# 1. Decide the feature schema and label column
FEATURE_COLS = ["HR_mean", "HRV_mean", "MI_mean", ...]  # or use od_features
init_problem("cardiac_af_vs_nsr",
             feature_cols=FEATURE_COLS,
             label_col="diagnosis",
             subject_id_col="patient_id",
             test_frac=0.25,
             classifier="rf")  # or "logreg"

# 2. Build first batch and ingest
df = pd.DataFrame([...])  # subj, label, features
ingest("cardiac_af_vs_nsr", df)
fit_and_evaluate("cardiac_af_vs_nsr", note="first_batch")

# 3. Later, ingest individual patients via CLI
# python3 od_ingest_patient.py cardiac_af_vs_nsr --wfdb 04015 --pn-dir afdb \
#         --subj-id pn_04015 --label AF --mode cardiac
```

## How "gets better with new data" works

1. First ingest splits data into train + frozen held-out test (test_frac
   of total, stratified by class).
2. Every subsequent ingest adds ONLY to the train pool. Test set stays
   frozen so accuracy is comparable across runs.
3. `fit_and_evaluate()` retrains the classifier from scratch on the full
   train pool (no online learning yet) and evaluates on the test set +
   stratified k-fold CV on train.
4. History log captures (run_id, n_train, n_test, accuracy, balanced,
   per-class sensitivity).
5. `show_history()` prints the table with deltas — every run shows how
   much accuracy moved since the previous run.

Demo from Session 9:

```
 run   n_tr   n_te  cls   heldout   cv_acc   cv_bal     d_cv  note
   1     14      4    2     0.500    0.429    0.396      -    small_initial
   2     33      4    2     0.500    0.606    0.567   +0.177  added_more
   3     40      4    2     0.500    0.575    0.561   -0.031  full_cohort
   4     40      4    2     0.500    0.575    0.561   +0.000  no_new_data
```

Run 1->2: +14 more patients, +18% CV accuracy.
Run 4 same data as Run 3: exactly 0 delta (deterministic).

## What's NOT yet in

- **Online learning** (true partial_fit, no refit-from-scratch). The
  current refit-from-scratch is more accurate at small N but doesn't
  scale to millions of patients. Swap in SGDClassifier or
  MultinomialNB if you need true incremental.
- **Active learning** (uncertainty sampling -> ask the user to label
  the patient the model is most unsure about). Hook point would be
  in `fit_and_evaluate`, returning the lowest-confidence training
  examples for review.
- **Multi-coupled-pair features per patient** (e.g., HR/HRV *and*
  HR/ABP from the same recording). Today's extractor returns one
  feature dict per (record, mode); to fuse modes per patient you'd
  union the feature dicts before ingest.
- **Cross-problem transfer** (train on cardiac AF, use as warm-start
  for sepsis detection). Architecture supports it but no glue yet.

## Caveats from Session 9 work

- Effect sizes at the group level (d=0.4-0.8) do NOT guarantee per-patient
  classifiability. With small per-class n, classifiers default to
  majority-class. Watch balanced accuracy, not raw accuracy.
- The held-out test set is small (n=test_frac * total) unless the initial
  batch is large. k-fold CV on the train pool is the more stable metric
  at small N.
- Re-running with the SAME data shows zero delta. The system gets better
  only when NEW patients are ingested.

## Files in the platform

| File | Purpose |
|---|---|
| `od_features.py` | Raw signal -> feature dict |
| `od_learning_store.py` | Persistent store + classifier + history |
| `od_ingest_patient.py` | CLI: ingest one patient |
| `od_predict_patient.py` | CLI: predict on one patient |
| `od_learn_cerebro.py` | Batch driver for cerebro problem (demo) |
| `OD_LEARNING_README.md` | This file |
