"""Persistent learning store for OD disease-detection problems.

Each "problem" has:
  - A growing feature store (subjects x features + label)
  - A held-out test set (frozen at first ingest, never retrained on)
  - A persistent classifier that retrains from the full store each run
  - A history log of (run_id, n_train, n_test, accuracy, balanced_acc,
    per-class sensitivity)

Every run reads existing state, ingests new examples (if provided),
retrains, evaluates on the held-out test set, and writes:
  - Updated feature store
  - Updated trained model
  - New history entry showing accuracy delta vs prior runs

This is the "gets better every run" infrastructure the project was
missing. New cohorts → bigger training set → tracked improvement.

Storage layout (one folder per problem):
  store/<problem_id>/
    features.csv          (full feature store, includes train+test marker)
    model.pkl             (trained sklearn pipeline)
    history.jsonl         (one JSON line per run)
    metadata.json         (problem schema: feature cols, label col, etc.)
"""
import datetime as dt
import json
import os
import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             confusion_matrix)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

STORE_ROOT = Path("store")


def _problem_dir(problem_id: str) -> Path:
    p = STORE_ROOT / problem_id
    p.mkdir(parents=True, exist_ok=True)
    return p


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def init_problem(problem_id: str, feature_cols: list[str], label_col: str,
                 subject_id_col: str = "subj", test_frac: float = 0.25,
                 random_state: int = 0,
                 classifier: str = "rf") -> None:
    """Initialize a new problem. Idempotent: doesn't overwrite if exists."""
    d = _problem_dir(problem_id)
    meta_path = d / "metadata.json"
    if meta_path.exists():
        print(f"  problem={problem_id} already initialized")
        return
    meta = {
        "problem_id": problem_id,
        "created_at": _now(),
        "feature_cols": feature_cols,
        "label_col": label_col,
        "subject_id_col": subject_id_col,
        "test_frac": test_frac,
        "random_state": random_state,
        "classifier": classifier,
    }
    meta_path.write_text(json.dumps(meta, indent=2))
    print(f"  problem={problem_id} initialized: "
          f"{len(feature_cols)} features, label={label_col}, "
          f"test_frac={test_frac}, classifier={classifier}")


def _load_meta(problem_id: str) -> dict:
    p = _problem_dir(problem_id) / "metadata.json"
    if not p.exists():
        raise RuntimeError(f"problem {problem_id} not initialized (run init_problem first)")
    return json.loads(p.read_text())


def _make_classifier(kind: str):
    if kind == "rf":
        return RandomForestClassifier(
            n_estimators=200, max_depth=6, min_samples_leaf=2,
            random_state=0, class_weight="balanced", n_jobs=-1)
    elif kind == "logreg":
        return Pipeline([
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=4000, C=1.0,
                                       class_weight="balanced",
                                       solver="lbfgs",
                                       multi_class="auto"))])
    elif kind == "knn1":
        # k=1 nearest neighbor with standardized features. Good for tiny
        # cohorts where within-class tightness is much smaller than
        # between-class distance (e.g., simulator domain signatures).
        return Pipeline([
            ("scale", StandardScaler()),
            ("clf", KNeighborsClassifier(n_neighbors=1))])
    else:
        raise ValueError(f"unknown classifier kind {kind!r}")


def ingest(problem_id: str, df: pd.DataFrame) -> int:
    """Add new (subject_id, features, label) rows to the store.

    On the very first ingest, splits new rows into train + held-out
    test using meta['test_frac']. Subsequent ingests add ONLY to the
    train pool (never grow the test set), so accuracy comparisons stay
    apples-to-apples.

    Returns number of new rows added.
    """
    meta = _load_meta(problem_id)
    d = _problem_dir(problem_id)
    fp = d / "features.csv"
    sid_col = meta["subject_id_col"]

    cols = [sid_col, meta["label_col"]] + meta["feature_cols"] + ["split"]
    # Filter caller df to required columns; add 'split' as None for now
    if sid_col not in df.columns or meta["label_col"] not in df.columns:
        raise RuntimeError(f"df missing {sid_col} or {meta['label_col']}")
    missing = [c for c in meta["feature_cols"] if c not in df.columns]
    if missing:
        raise RuntimeError(f"df missing features: {missing[:5]}...")
    new = df[[sid_col, meta["label_col"]] + meta["feature_cols"]].copy()

    if fp.exists():
        existing = pd.read_csv(fp)
        existing_ids = set(existing[sid_col].astype(str))
        # Filter out duplicate subjects (idempotent ingest)
        new = new[~new[sid_col].astype(str).isin(existing_ids)]
        if len(new) == 0:
            print(f"  ingest: no new subjects (all {len(df)} already present)")
            return 0
        # All new rows go to TRAIN; test set is frozen
        new["split"] = "train"
        combined = pd.concat([existing, new[cols]], ignore_index=True)
    else:
        # First ingest: split into train / test by meta['test_frac'].
        # If test_frac == 0: everything goes to train (CV-only mode for
        # very small datasets).
        rng = np.random.RandomState(meta["random_state"])
        n = len(new)
        labels = new[meta["label_col"]].values
        unique = pd.Series(labels).value_counts().index.tolist()
        if meta["test_frac"] <= 0:
            new["split"] = "train"
        else:
            test_idx = []
            for lbl in unique:
                cls_idx = np.where(labels == lbl)[0]
                cls_test_n = max(1 if len(cls_idx) >= 2 else 0,
                                 int(round(len(cls_idx) * meta["test_frac"])))
                picked = rng.choice(cls_idx,
                                    size=min(cls_test_n, len(cls_idx)),
                                    replace=False)
                test_idx.extend(picked.tolist())
            test_mask = np.zeros(n, dtype=bool)
            test_mask[test_idx] = True
            new["split"] = np.where(test_mask, "test", "train")
        combined = new[cols]

    combined.to_csv(fp, index=False)
    print(f"  ingest: added {len(new)} new rows  total store size={len(combined)}")
    print(f"  split breakdown: {dict(combined['split'].value_counts())}")
    print(f"  label distribution (train): "
          f"{dict(combined[combined['split']=='train'][meta['label_col']].value_counts())}")
    return len(new)


def fit_and_evaluate(problem_id: str, note: str = "") -> dict:
    """Retrain on full train pool, evaluate on frozen test set, log to history."""
    meta = _load_meta(problem_id)
    d = _problem_dir(problem_id)
    fp = d / "features.csv"
    if not fp.exists():
        raise RuntimeError(f"no features yet for {problem_id}")
    df = pd.read_csv(fp)
    train = df[df["split"] == "train"].copy()
    test = df[df["split"] == "test"].copy()
    if len(train) < 5 or len(test) < 2:
        print(f"  WARN: train={len(train)} test={len(test)} too small")
    # Drop any classes that appear only in test (model can't learn them)
    train_cls = set(train[meta["label_col"]].unique())
    test = test[test[meta["label_col"]].isin(train_cls)]

    Xtr = train[meta["feature_cols"]].fillna(0.0).values
    ytr = train[meta["label_col"]].values
    Xte = test[meta["feature_cols"]].fillna(0.0).values
    yte = test[meta["label_col"]].values

    clf = _make_classifier(meta["classifier"])
    clf.fit(Xtr, ytr)
    yhat = clf.predict(Xte) if len(Xte) > 0 else np.array([])
    acc = float(accuracy_score(yte, yhat)) if len(yhat) > 0 else None
    bal = float(balanced_accuracy_score(yte, yhat)) if len(yhat) > 0 else None
    labels = sorted(train_cls)
    cm = (confusion_matrix(yte, yhat, labels=labels).tolist()
          if len(yhat) > 0 else None)

    # ALSO: stratified k-fold CV on train pool — gives a stable accuracy
    # estimate that scales smoothly with cohort size (small held-out set
    # has too much variance to track learning).
    cv_acc = cv_bal = None
    cv_cm = None
    cv_per_class = {}
    try:
        cls_counts = pd.Series(ytr).value_counts()
        if cls_counts.min() >= 2 and len(np.unique(ytr)) >= 2:
            n_splits = int(min(5, cls_counts.min(), len(ytr) // 2))
            n_splits = max(2, n_splits)
            cv = StratifiedKFold(n_splits=n_splits, shuffle=True,
                                 random_state=meta["random_state"])
            clf_cv = _make_classifier(meta["classifier"])
            ycv = cross_val_predict(clf_cv, Xtr, ytr, cv=cv)
            cv_acc = float(accuracy_score(ytr, ycv))
            cv_bal = float(balanced_accuracy_score(ytr, ycv))
            cv_cm = confusion_matrix(ytr, ycv, labels=labels).tolist()
            for i, l in enumerate(labels):
                n_t = int((ytr == l).sum())
                cv_per_class[l] = {
                    "n": n_t, "correct": int(cv_cm[i][i]) if n_t > 0 else 0,
                    "sensitivity": float(cv_cm[i][i] / n_t) if n_t > 0 else None,
                }
    except Exception as e:
        print(f"  CV failed: {e}")
    # Per-class sensitivity
    per_class = {}
    if cm is not None:
        for i, l in enumerate(labels):
            n_t = int((yte == l).sum())
            per_class[l] = {
                "n_test": n_t,
                "correct": int(cm[i][i]) if n_t > 0 else 0,
                "sensitivity": float(cm[i][i] / n_t) if n_t > 0 else None,
            }

    # Persist model
    with (d / "model.pkl").open("wb") as f:
        pickle.dump(clf, f)

    # Append to history
    hist_path = d / "history.jsonl"
    prev = []
    if hist_path.exists():
        prev = [json.loads(line) for line in hist_path.read_text().splitlines()
                if line.strip()]
    run_id = len(prev) + 1
    entry = {
        "run_id": run_id, "ts": _now(), "note": note,
        "n_train": int(len(train)), "n_test": int(len(test)),
        "n_classes": len(labels), "labels": labels,
        "accuracy": acc, "balanced_accuracy": bal,
        "confusion_matrix": cm, "per_class": per_class,
        "cv_accuracy": cv_acc, "cv_balanced_accuracy": cv_bal,
        "cv_confusion_matrix": cv_cm, "cv_per_class": cv_per_class,
    }
    with hist_path.open("a") as f:
        f.write(json.dumps(entry) + "\n")
    # Delta vs previous run
    prev_acc = prev[-1]["accuracy"] if prev and prev[-1].get("accuracy") is not None else None
    prev_bal = prev[-1]["balanced_accuracy"] if prev and prev[-1].get("balanced_accuracy") is not None else None
    if prev_acc is not None and acc is not None:
        entry["delta_acc"] = acc - prev_acc
        entry["delta_bal"] = bal - prev_bal
    return entry


def show_history(problem_id: str) -> None:
    d = _problem_dir(problem_id)
    hp = d / "history.jsonl"
    if not hp.exists():
        print(f"  no history for {problem_id} yet"); return
    runs = [json.loads(l) for l in hp.read_text().splitlines() if l.strip()]
    print(f"\n=== {problem_id} history ({len(runs)} runs) ===")
    print(f"{'run':>4s}  {'n_tr':>5s}  {'n_te':>5s}  {'cls':>3s}  "
          f"{'heldout':>8s}  {'cv_acc':>7s}  {'cv_bal':>7s}  "
          f"{'d_cv':>7s}  note")
    prev_cv = None
    for r in runs:
        a = r.get("accuracy"); cva = r.get("cv_accuracy"); cvb = r.get("cv_balanced_accuracy")
        a_s = f"{a:.3f}" if a is not None else "   -  "
        cva_s = f"{cva:.3f}" if cva is not None else "   -  "
        cvb_s = f"{cvb:.3f}" if cvb is not None else "   -  "
        d_cv = (cva - prev_cv) if prev_cv is not None and cva is not None else None
        d_cv_s = f"{d_cv:+.3f}" if d_cv is not None else "    -  "
        print(f"{r['run_id']:>4d}  {r['n_train']:>5d}  {r['n_test']:>5d}  "
              f"{r['n_classes']:>3d}  {a_s:>8s}  {cva_s:>7s}  {cvb_s:>7s}  "
              f"{d_cv_s:>7s}  {r.get('note', '')}")
        if cva is not None: prev_cv = cva
