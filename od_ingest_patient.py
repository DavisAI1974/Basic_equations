"""Ingest ONE patient into an OD learning problem and refit.

Primary use case: you have one new patient. Run this script to:
  1. Extract per-domain stack features (if --wfdb given) OR read them
     from a JSON file you already prepared
  2. Add them to the problem's feature store (deduplicated by subject id)
  3. Retrain on the full updated train pool
  4. Re-evaluate on the frozen held-out test set + k-fold CV
  5. Log the run; print accuracy delta vs prior run

Usage A: JSON dict you already prepared
  python3 od_ingest_patient.py <problem_id> --json patient_s0999.json

  JSON schema: {"subj": "<id>", "label": "<class>", <feature_name>: <float>, ...}

Usage B: WFDB record on PhysioNet
  python3 od_ingest_patient.py <problem_id> --wfdb 16265 --pn-dir nsrdb \\
      --subj-id pn_16265 --label NSR --mode cardiac

Usage C: local WFDB record
  python3 od_ingest_patient.py <problem_id> --wfdb s0030DA \\
      --local-dir data/cvd/sample --subj-id s0030 --label control \\
      --mode baroreflex

The problem must already be initialized (an earlier batch ingest set
its schema). Subjects with IDs already present in the store are
skipped silently (idempotent).
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from od_learning_store import (_load_meta, ingest, fit_and_evaluate,
                               show_history)
from od_features import extract_from_wfdb_record


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("problem_id")
    ap.add_argument("--json", help="JSON file with feature dict")
    ap.add_argument("--wfdb", help="WFDB record name (extract on the fly)")
    ap.add_argument("--pn-dir", help="PhysioNet db dir for wfdb fetch")
    ap.add_argument("--local-dir", help="Local dir containing .dat/.hea")
    ap.add_argument("--mode", default="cardiac",
                    choices=["cardiac", "baroreflex"])
    ap.add_argument("--subj-id", help="subject id (required when --wfdb used)")
    ap.add_argument("--label", help="ground-truth label (required when --wfdb used)")
    ap.add_argument("--note", default="", help="annotation for history log")
    args = ap.parse_args()

    meta = _load_meta(args.problem_id)
    print(f"Problem: {args.problem_id}")
    print(f"Schema: {len(meta['feature_cols'])} features, "
          f"label={meta['label_col']}, sid_col={meta['subject_id_col']}")

    if args.wfdb:
        if not args.subj_id or not args.label:
            raise SystemExit("--wfdb requires --subj-id and --label")
        print(f"\nExtracting features from {args.wfdb}  mode={args.mode}")
        payload = extract_from_wfdb_record(
            args.wfdb, pn_dir=args.pn_dir, local_dir=args.local_dir,
            mode=args.mode)
        payload[meta["subject_id_col"]] = args.subj_id
        payload[meta["label_col"]] = args.label
        print(f"  n_rr={payload['n_rr']}  n_sup={payload['n_sup']}  "
              f"HR_mean={payload['HR_mean']:.1f}  MI_mean={payload['MI_mean']:.3f}")
    elif args.json:
        payload = json.loads(Path(args.json).read_text())
    else:
        raise SystemExit("Provide either --json or --wfdb")

    print(f"\nIngesting: subj={payload.get(meta['subject_id_col'], '?')} "
          f"label={payload.get(meta['label_col'], '?')}")

    missing = [c for c in meta["feature_cols"] if c not in payload]
    if missing:
        print(f"WARN: {len(missing)} required features missing from payload: "
              f"{missing[:5]}...  filling with 0")
        for c in missing: payload[c] = 0.0

    df = pd.DataFrame([payload])
    n_added = ingest(args.problem_id, df)

    note = args.note or (f"wfdb:{args.wfdb}" if args.wfdb
                         else f"json:{Path(args.json).name}")
    entry = fit_and_evaluate(args.problem_id, note=note)

    print(f"\n=== Run {entry['run_id']} ===")
    print(f"  added: {n_added}  train pool: {entry['n_train']}  "
          f"held-out test: {entry['n_test']}")
    print(f"  held-out: acc={entry.get('accuracy')}  bal={entry.get('balanced_accuracy')}")
    print(f"  CV-on-train: acc={entry.get('cv_accuracy')}  bal={entry.get('cv_balanced_accuracy')}")

    show_history(args.problem_id)


if __name__ == "__main__":
    main()
