"""Re-classify cerebro probe with RICH per-subject features.

Add temporal features from the (H_a, H_b, MI) sequences themselves:
  - FFT band powers (low/mid/hi)
  - Lag-1 / lag-3 autocorrelation
  - First-difference std (jaggedness)
Then re-run LOO with combined scalar + temporal features.

Tests whether per-patient classification improves beyond chance when
the temporal dynamics of the operator-space trajectory are included.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import (accuracy_score, confusion_matrix,
                             balanced_accuracy_score)


def temporal_features(seq, prefix):
    s = np.asarray(seq, float)
    s = s - s.mean()
    n = len(s)
    out = {}
    spec = np.abs(np.fft.rfft(s)) ** 2
    freqs = np.fft.rfftfreq(n)
    total_p = spec.sum() + 1e-12
    out[f"{prefix}_p_low"] = float(spec[(freqs < 0.05)].sum() / total_p)
    out[f"{prefix}_p_mid"] = float(spec[(freqs >= 0.05) & (freqs < 0.20)].sum() / total_p)
    out[f"{prefix}_p_hi"] = float(spec[(freqs >= 0.20)].sum() / total_p)
    if n > 4 and s.std() > 1e-9:
        out[f"{prefix}_ac1"] = float(np.corrcoef(s[:-1], s[1:])[0, 1])
        out[f"{prefix}_ac3"] = float(np.corrcoef(s[:-3], s[3:])[0, 1]) if n > 6 else 0.0
    else:
        out[f"{prefix}_ac1"] = 0.0; out[f"{prefix}_ac3"] = 0.0
    out[f"{prefix}_diff_std"] = float(np.diff(s).std()) if n > 1 else 0.0
    return out


with open("cerebro_disease_probe.json") as f:
    rec = json.load(f)
ok = [r for r in rec if r["ok"]]
print(f"Subjects: {len(ok)}")

# Build per-subject feature dict, including temporal features
rows = []
for r in ok:
    feats = {
        "subj": r["subj"],
        "HR_mean": r["HR_mean"], "HR_std": r["HR_std"],
        "ABP_mean": r["ABP_mean"], "ABP_std": r["ABP_std"],
        "Ha_mean": r["Ha_mean"], "Hb_mean": r["Hb_mean"],
        "MI_mean": r["MI_mean"], "MI_std": r["MI_std"],
        "n_sup": r["n_supertimes"],
    }
    f = r["fits"]
    feats.update({
        "r2_poly1": f["poly1"]["r2"], "r2_poly2": f["poly2"]["r2"],
        "r2_exp": f.get("exp_INFO025", {}).get("r2", 0.0),
        "r2_phys": f["physics_INFO025"]["r2"],
        "poly1_a": f["poly1"]["a"], "poly1_b": f["poly1"]["b"],
        "poly2_a": f["poly2"]["a"], "poly2_b": f["poly2"]["b"],
        "poly2_d": f["poly2"]["d"], "poly2_e": f["poly2"]["e"],
        "poly2_f": f["poly2"]["f"],
        "exp_A": f.get("exp_INFO025", {}).get("A", 0.0),
        "exp_B": f.get("exp_INFO025", {}).get("B", 0.0),
        "phys_alpha": f["physics_INFO025"]["alpha"],
    })
    feats.update(temporal_features(r["Ha"], "Ha"))
    feats.update(temporal_features(r["Hb"], "Hb"))
    feats.update(temporal_features(r["MI"], "MI"))
    rows.append(feats)
feat_df = pd.DataFrame(rows)

# Sanitize
for c in feat_df.columns:
    if c != "subj":
        feat_df[c] = pd.to_numeric(feat_df[c], errors="coerce")
        feat_df[c] = feat_df[c].replace([np.inf, -np.inf], np.nan)
        feat_df[c] = feat_df[c].fillna(feat_df[c].median())

cohort = pd.read_csv("data/cvd/cohort_strat.csv").rename(
    columns={"SubjectID_lc": "subj"})
df = feat_df.merge(cohort[["subj", "strat1", "strat2", "strat3"]],
                   on="subj", how="left")
print(f"Joined: {len(df)} subjects, {len(feat_df.columns)-1} features")

FEATURE_COLS = [c for c in feat_df.columns if c != "subj"]

CLASSIFIERS = {
    "LogReg": Pipeline([("scale", StandardScaler()),
                        ("clf", LogisticRegression(
                            max_iter=2000, C=0.5,
                            class_weight="balanced", solver="liblinear"))]),
    "RF(50)": RandomForestClassifier(n_estimators=50, max_depth=5,
                                     random_state=0, class_weight="balanced"),
}


def run_loo(strat_col, min_n=2):
    print(f"\n{'='*78}")
    print(f"Stratification: {strat_col}  (rich features)")
    print('='*78)
    sub = df[df[strat_col].notna()].copy()
    counts = sub.groupby(strat_col).size()
    keep = counts[counts >= min_n].index.tolist()
    sub = sub[sub[strat_col].isin(keep)].copy()
    print(f"Cells: {dict(counts.loc[keep])}  (n_total={len(sub)})")
    X = sub[FEATURE_COLS].values
    y = sub[strat_col].values
    n_cls = len(np.unique(y))
    chance = 1.0 / n_cls
    majority = pd.Series(y).value_counts().iloc[0] / len(y)
    print(f"Chance={chance:.3f}  Majority={majority:.3f}\n")

    for name, clf in CLASSIFIERS.items():
        try:
            y_pred = cross_val_predict(clf, X, y, cv=LeaveOneOut())
            acc = accuracy_score(y, y_pred)
            bal = balanced_accuracy_score(y, y_pred)
            labels = sorted(np.unique(y))
            cm = confusion_matrix(y, y_pred, labels=labels)
            print(f"  {name:8s}  LOO acc={acc:.3f}  balanced={bal:.3f}")
            print(f"    Per-class sensitivity:")
            for i, l in enumerate(labels):
                n_t = (y == l).sum()
                print(f"      {l:<24s} n={n_t:3d}  {cm[i,i]:3d}/{n_t:<3d}  "
                      f"sens={cm[i,i]/n_t:.3f}")
        except Exception as e:
            print(f"  {name:8s}  ERR: {e}")


for strat in ["strat1", "strat2", "strat3"]:
    run_loo(strat)

print("\n" + "="*78)
print("Binary any-disease vs control (rich features)")
print("="*78)
sub = df.copy()
sub["binary"] = sub["strat1"].apply(lambda x: "control" if x == "control" else "any_disease")
X = sub[FEATURE_COLS].values
y = sub["binary"].values
print(f"Cells: {dict(pd.Series(y).value_counts())}")
for name, clf in CLASSIFIERS.items():
    y_pred = cross_val_predict(clf, X, y, cv=LeaveOneOut())
    acc = accuracy_score(y, y_pred)
    bal = balanced_accuracy_score(y, y_pred)
    labels = sorted(np.unique(y))
    cm = confusion_matrix(y, y_pred, labels=labels)
    print(f"  {name:8s}  LOO acc={acc:.3f}  balanced={bal:.3f}")
    for i, l in enumerate(labels):
        n_t = (y == l).sum()
        print(f"    {l:<20s} sens={cm[i,i]/n_t:.3f}  ({cm[i,i]}/{n_t})")

# Save the enriched feature table for downstream use
df.to_csv("cerebro_rich_features.csv", index=False)
print("\nSaved cerebro_rich_features.csv")
