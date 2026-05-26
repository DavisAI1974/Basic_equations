"""Demonstrate the OD learning store on the cerebro disease-detection
problem.

Each run:
  1. Initializes the problem if not yet
  2. Ingests any new subjects from cerebro_rich_features.csv (idempotent;
     duplicates skipped)
  3. Retrains on the full train pool
  4. Evaluates on the frozen held-out test set
  5. Appends to history.jsonl
  6. Prints history table showing accuracy delta across runs

Re-running with the SAME data shows no improvement (we proved this is
fundamental: same data + same seed = same model). Running after a NEW
cohort is added should show accuracy delta.
"""
import json
import sys

import pandas as pd

from od_learning_store import (init_problem, ingest, fit_and_evaluate,
                               show_history)

PROBLEM = "cerebro_disease_binary"

# Load enriched cerebro features
df = pd.read_csv("cerebro_rich_features.csv")
print(f"Loaded {len(df)} cerebro subjects, {df.shape[1]} columns")

# Binary "any disease vs control" label (best per-patient signal we have)
df["label"] = df["strat1"].apply(
    lambda x: "control" if x == "control" else "any_disease")

# All numeric feature columns (everything except IDs and labels)
FEATURE_COLS = [c for c in df.columns
                if c not in ("subj", "label", "strat1", "strat2", "strat3")]
print(f"Using {len(FEATURE_COLS)} features per subject")

# Init problem (idempotent)
init_problem(PROBLEM, feature_cols=FEATURE_COLS, label_col="label",
             subject_id_col="subj", test_frac=0.25, random_state=0,
             classifier="rf")

# Allow caller to control ingest fraction (so a SECOND run with --extra N
# can simulate "new cohort arrived" by ingesting more than were ingested
# initially). Default behavior: ingest all available subjects.
fraction = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
if fraction < 1.0:
    n_keep = max(8, int(round(len(df) * fraction)))
    df = df.sample(n=n_keep, random_state=int(fraction * 1000))
    print(f"  Subsampled to {len(df)} subjects (fraction={fraction})")

# Ingest
n_added = ingest(PROBLEM, df[["subj", "label"] + FEATURE_COLS])

# Fit + evaluate + log
note = sys.argv[2] if len(sys.argv) > 2 else f"frac={fraction}"
entry = fit_and_evaluate(PROBLEM, note=note)
print(f"\n=== Run {entry['run_id']} ===")
print(f"  n_train={entry['n_train']}  n_test={entry['n_test']}")
print(f"  accuracy={entry['accuracy']}  balanced={entry['balanced_accuracy']}")
print(f"  per-class:")
for lbl, m in entry["per_class"].items():
    print(f"    {lbl:<14s}  n={m['n_test']}  correct={m['correct']}  "
          f"sens={m['sensitivity']}")
if "delta_acc" in entry:
    print(f"  delta_acc={entry['delta_acc']:+.3f}  delta_bal={entry['delta_bal']:+.3f}")

show_history(PROBLEM)
