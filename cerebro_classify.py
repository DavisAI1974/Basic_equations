"""Did we detect diseases correctly per patient?

Take the per-subject features from cerebro_disease_probe and run
leave-one-out classification per stratification. Compare to random
baseline and majority-class baseline.

This is the validation question: signature differences at the group
level don't necessarily translate to per-patient detection.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (accuracy_score, confusion_matrix,
                             classification_report, balanced_accuracy_score)

df = pd.read_csv("cerebro_per_subject_with_groups.csv")

# Feature columns: aggregate stats + fit coefficients (NOT raw arrays)
FEATURES = [
    "HR_mean", "HR_std", "ABP_mean", "ABP_std",
    "MI_mean", "MI_std", "Ha_mean", "Hb_mean",
    "r2_poly1", "r2_poly2", "r2_exp_INFO025", "r2_physics_INFO025",
    "poly1_a", "poly1_b", "poly1_c",
    "poly2_a", "poly2_b", "poly2_d", "poly2_e", "poly2_f", "poly2_c",
    "a_hat", "b_hat", "exp_A", "exp_B",
]
FEATURES = [f for f in FEATURES if f in df.columns]
print(f"Using {len(FEATURES)} features per subject")


def run_classification(strat_col, min_cell_size=2):
    print(f"\n{'=' * 80}")
    print(f"Stratification: {strat_col}")
    print('=' * 80)
    sub = df[df[strat_col].notna()].copy()
    # Drop cells smaller than min_cell_size
    counts = sub.groupby(strat_col).size()
    keep = counts[counts >= min_cell_size].index.tolist()
    sub = sub[sub[strat_col].isin(keep)].copy()
    print(f"Cells (n>={min_cell_size}): {dict(counts.loc[keep])}")
    print(f"Total subjects: {len(sub)}")

    # Replace any inf/nan features with column-median
    X = sub[FEATURES].replace([np.inf, -np.inf], np.nan).copy()
    for c in FEATURES:
        X[c] = X[c].fillna(X[c].median())
    y = sub[strat_col].values
    n_classes = len(np.unique(y))
    print(f"Classes: {sorted(np.unique(y))}")

    # Baselines
    majority = pd.Series(y).value_counts().iloc[0] / len(y)
    random_chance = 1.0 / n_classes
    print(f"Random chance: {random_chance:.3f}  Majority class: {majority:.3f}")

    # Classifiers
    classifiers = {
        "LogReg(L2)": Pipeline([
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=2000, class_weight="balanced",
                                       C=1.0)),
        ]),
        "RandomForest(100)": RandomForestClassifier(
            n_estimators=100, max_depth=4, random_state=0,
            class_weight="balanced"),
    }

    for name, clf in classifiers.items():
        loo = LeaveOneOut()
        y_pred = cross_val_predict(clf, X.values, y, cv=loo)
        acc = accuracy_score(y, y_pred)
        bal = balanced_accuracy_score(y, y_pred)
        print(f"\n  {name}  LOO acc={acc:.3f}  balanced={bal:.3f}  "
              f"(chance={random_chance:.3f}, majority={majority:.3f})")
        # Confusion matrix
        labels = sorted(np.unique(y))
        cm = confusion_matrix(y, y_pred, labels=labels)
        print(f"  Confusion matrix (rows=true, cols=pred):")
        print("  " + " " * 22 + "  ".join(f"{l[:14]:>14s}" for l in labels))
        for i, l in enumerate(labels):
            row = "  ".join(f"{cm[i, j]:14d}" for j in range(len(labels)))
            print(f"  {l[:20]:>22s}  {row}")
        # Per-class metrics
        print(f"  Per-class accuracy (sensitivity = recall):")
        for i, l in enumerate(labels):
            n_true = (y == l).sum()
            n_correct = cm[i, i]
            sens = n_correct / n_true if n_true > 0 else 0
            print(f"    {l:<24s} n={n_true:3d}  correct={n_correct:3d}  "
                  f"sensitivity={sens:.3f}")


for strat in ["strat1", "strat2", "strat3"]:
    run_classification(strat, min_cell_size=2)

# Also: try a binary "disease vs control" cut for max statistical power
print("\n" + "=" * 80)
print("Binary: any-disease vs control (largest cells)")
print("=" * 80)
sub = df.copy()
sub["binary"] = sub["strat1"].apply(lambda x: "control" if x == "control" else "any_disease")
X = sub[FEATURES].replace([np.inf, -np.inf], np.nan)
for c in FEATURES:
    X[c] = X[c].fillna(X[c].median())
y = sub["binary"].values
print(f"Cells: {dict(pd.Series(y).value_counts())}")
clf = Pipeline([("scale", StandardScaler()),
                ("clf", LogisticRegression(max_iter=2000,
                                           class_weight="balanced"))])
loo = LeaveOneOut()
y_pred = cross_val_predict(clf, X.values, y, cv=loo)
acc = accuracy_score(y, y_pred)
bal = balanced_accuracy_score(y, y_pred)
print(f"\n  LogReg LOO acc={acc:.3f}  balanced={bal:.3f}  (chance=0.500)")
labels = sorted(np.unique(y))
cm = confusion_matrix(y, y_pred, labels=labels)
print(f"  Confusion matrix:")
print("  " + " " * 18 + "  ".join(f"{l:>14s}" for l in labels))
for i, l in enumerate(labels):
    row = "  ".join(f"{cm[i, j]:14d}" for j in range(len(labels)))
    print(f"  {l:>18s}  {row}")
for i, l in enumerate(labels):
    n_true = (y == l).sum()
    n_correct = cm[i, i]
    sens = n_correct / n_true if n_true > 0 else 0
    print(f"    {l:<20s} sensitivity={sens:.3f}  (n={n_true})")
