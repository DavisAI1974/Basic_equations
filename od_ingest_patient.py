"""Ingest ONE patient into an OD learning problem and refit.

Primary use case: you have one new patient's per-domain stack features
(extracted from their ECG / cerebro / whatever signal pair) plus their
ground-truth label. Run this script to:
  1. Add them to the problem's feature store (deduplicated by subject_id)
  2. Retrain on the full updated train pool
  3. Re-evaluate on the frozen held-out test set + k-fold CV
  4. Log the run; print accuracy delta vs prior run

Usage:
  python3 od_ingest_patient.py <problem_id> <features.json>

features.json schema:
  {
    "subj":  "<unique subject id string>",
    "label": "<class label string>",
    "<feature_name_1>": <float>,
    ...
  }

The problem must already be initialized (an earlier batch ingest set
its schema). Subjects with IDs already present in the store are skipped
silently (idempotent).

Example:
  python3 od_ingest_patient.py cerebro_disease_binary patient_s0999.json

This is the feedback loop: every new patient -> bigger training pool ->
recorded accuracy delta. Re-running with no new patients shows zero
delta (model is deterministic given data + seed).
"""
import json
import sys
from pathlib import Path

import pandas as pd

from od_learning_store import (_load_meta, ingest, fit_and_evaluate,
                               show_history)

if len(sys.argv) < 3:
    print(__doc__)
    sys.exit(1)

problem_id = sys.argv[1]
features_path = Path(sys.argv[2])
note = sys.argv[3] if len(sys.argv) > 3 else f"ingest:{features_path.name}"

meta = _load_meta(problem_id)
print(f"Problem: {problem_id}")
print(f"Schema: {len(meta['feature_cols'])} features, label={meta['label_col']}")

# Load single patient
payload = json.loads(features_path.read_text())
print(f"\nIngesting: {payload.get(meta['subject_id_col'], '?')}")
print(f"  label: {payload.get(meta['label_col'], '?')}")

# Validate required fields present
missing = [c for c in meta["feature_cols"] if c not in payload]
if missing:
    print(f"ERROR: features.json missing required features: {missing[:5]}...")
    sys.exit(2)

# Convert to single-row DataFrame
df = pd.DataFrame([{**payload}])

# Ingest + refit + eval
n_added = ingest(problem_id, df)
entry = fit_and_evaluate(problem_id, note=note)

print(f"\n=== Run {entry['run_id']} ===")
print(f"  added: {n_added}  train pool: {entry['n_train']}  "
      f"held-out test: {entry['n_test']}")
print(f"  held-out: acc={entry.get('accuracy')}  bal={entry.get('balanced_accuracy')}")
print(f"  CV-on-train: acc={entry.get('cv_accuracy')}  bal={entry.get('cv_balanced_accuracy')}")
if "delta_acc" in entry and entry["delta_acc"] is not None:
    print(f"  delta held-out: {entry['delta_acc']:+.3f} (acc)  "
          f"{entry['delta_bal']:+.3f} (bal)")

show_history(problem_id)
