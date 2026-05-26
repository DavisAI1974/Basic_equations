"""Predict on one patient WITHOUT ingesting them.

Use this when you have a new patient and want to know what the trained
model says, but you don't have a ground-truth label yet (or you want
to verify the prediction before deciding to ingest).

Usage:
  python3 od_predict_patient.py <problem_id> <features.json>
  python3 od_predict_patient.py <problem_id> --wfdb <rec_name> --pn-dir <db> --mode cardiac
  python3 od_predict_patient.py <problem_id> --wfdb <rec_name> --local-dir <path> --mode baroreflex

The model is loaded from store/<problem_id>/model.pkl (must have been
trained by at least one prior fit_and_evaluate run).

Prints: predicted label + class probabilities + per-feature contribution
(for tree models: feature importances of the trained model).
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from od_learning_store import _load_meta, _problem_dir
from od_features import extract_from_wfdb_record


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("problem_id")
    ap.add_argument("features_json", nargs="?",
                    help="JSON file with feature dict; ignored if --wfdb given")
    ap.add_argument("--wfdb", help="WFDB record name; loads via wfdb")
    ap.add_argument("--pn-dir", help="PhysioNet db dir for wfdb fetch")
    ap.add_argument("--local-dir", help="Local dir containing .dat/.hea")
    ap.add_argument("--mode", default="cardiac",
                    choices=["cardiac", "baroreflex"],
                    help="extractor to use if --wfdb")
    args = ap.parse_args()

    meta = _load_meta(args.problem_id)
    d = _problem_dir(args.problem_id)
    model_path = d / "model.pkl"
    if not model_path.exists():
        raise RuntimeError(f"no trained model for {args.problem_id} - "
                           f"run a fit_and_evaluate first")
    with model_path.open("rb") as f:
        clf = pickle.load(f)

    if args.wfdb:
        print(f"Extracting features from WFDB record: {args.wfdb}  mode={args.mode}")
        feats = extract_from_wfdb_record(
            args.wfdb, pn_dir=args.pn_dir, local_dir=args.local_dir,
            mode=args.mode)
        print(f"  n_rr={feats['n_rr']}  n_sup={feats['n_sup']}  "
              f"HR_mean={feats['HR_mean']:.1f}  MI_mean={feats['MI_mean']:.3f}")
    elif args.features_json:
        feats = json.loads(Path(args.features_json).read_text())
    else:
        raise SystemExit("Provide either features.json or --wfdb")

    # Build feature row in schema order
    missing = [c for c in meta["feature_cols"] if c not in feats]
    if missing:
        # Fall back: use 0 for any feature the extractor didn't produce
        print(f"WARN: {len(missing)} features missing from extraction: {missing[:5]}...")
        for c in missing: feats[c] = 0.0
    X = np.array([[feats[c] for c in meta["feature_cols"]]])

    pred = clf.predict(X)[0]
    print(f"\nPredicted class: {pred}")
    if hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(X)[0]
        classes = clf.classes_
        print(f"Class probabilities:")
        for c, p in sorted(zip(classes, proba), key=lambda kv: -kv[1]):
            print(f"  {c:<24s}  {p:.3f}")
        # Calibration note
        conf = max(proba)
        if conf < 0.6:
            print(f"  -- LOW CONFIDENCE ({conf:.3f}) -- consider getting "
                  f"ground truth and ingesting this patient")
    # Feature importance (if RandomForest-like)
    if hasattr(clf, "feature_importances_"):
        imp = clf.feature_importances_
        top = np.argsort(imp)[-8:][::-1]
        print(f"\nTop 8 model features by importance:")
        for i in top:
            print(f"  {meta['feature_cols'][i]:<24s}  imp={imp[i]:.3f}  "
                  f"this_patient={feats[meta['feature_cols'][i]]:.4f}")


if __name__ == "__main__":
    main()
