"""Cross-stratification read of cerebrovascular probe.

For each of the 3 stratifications:
- Per-group mean +/- inter-subject scatter on each fit's R^2 and key coefs
- Per-group linear-direction unit vector mean + cross-group cosines
- Effect-size markers (group separation in operator-space)
"""
import json
import numpy as np
import pandas as pd

with open("cerebro_disease_probe.json") as f:
    rec = json.load(f)
res = pd.DataFrame([r for r in rec if r["ok"]])

cohort = pd.read_csv("data/cvd/cohort_strat.csv")
cohort = cohort.rename(columns={"SubjectID_lc": "subj"})
df = res.merge(cohort, on="subj", how="left")
print(f"Joined: {len(df)} subjects")

# Extract scalar fit features for easy aggregation
def get_fit(row, family, key):
    f = row["fits"].get(family) if isinstance(row["fits"], dict) else None
    return None if f is None else f.get(key)

for fam in ["poly1", "poly2", "exp_INFO025", "physics_INFO025"]:
    df[f"r2_{fam}"] = df.apply(lambda r: get_fit(r, fam, "r2"), axis=1)
for c in ["a", "b", "c"]:
    df[f"poly1_{c}"] = df.apply(lambda r: get_fit(r, "poly1", c), axis=1)
for c in ["a", "b", "d", "e", "f", "c"]:
    df[f"poly2_{c}"] = df.apply(lambda r: get_fit(r, "poly2", c), axis=1)
df["exp_A"] = df.apply(lambda r: get_fit(r, "exp_INFO025", "A"), axis=1)
df["exp_B"] = df.apply(lambda r: get_fit(r, "exp_INFO025", "B"), axis=1)

# Normalize poly1 (a,b) into unit vector per subject
def unit_vec(row):
    v = np.array([row["poly1_a"], row["poly1_b"]])
    n = np.linalg.norm(v)
    return v / n if n > 0 else np.array([np.nan, np.nan])
uv = df.apply(unit_vec, axis=1, result_type="expand")
df["a_hat"], df["b_hat"] = uv[0], uv[1]


def summarize(strat_col):
    print("\n" + "=" * 88)
    print(f"  Stratification: {strat_col}")
    print("=" * 88)
    g = df.groupby(strat_col)
    sizes = g.size()
    rows = []
    for grp, sub in g:
        rows.append({
            "group": grp,
            "n": len(sub),
            "HR_mean": sub["HR_mean"].mean(),
            "HR_mean_std": sub["HR_mean"].std(),
            "ABP_mean": sub["ABP_mean"].mean(),
            "MI_mean": sub["MI_mean"].mean(),
            "MI_mean_std": sub["MI_mean"].std(),
            "r2_poly1_mean": sub["r2_poly1"].mean(),
            "r2_poly2_mean": sub["r2_poly2"].mean(),
            "r2_exp_mean": sub["r2_exp_INFO025"].mean(),
            "r2_phys_mean": sub["r2_physics_INFO025"].mean(),
            "a_hat_mean": sub["a_hat"].mean(),
            "b_hat_mean": sub["b_hat"].mean(),
            "a_hat_std": sub["a_hat"].std(),
            "b_hat_std": sub["b_hat"].std(),
        })
    summ = pd.DataFrame(rows)
    print(summ.to_string(index=False, float_format=lambda x: f"{x:+.3f}" if isinstance(x, float) else str(x)))

    # Pairwise cosines of mean linear direction
    grps = summ["group"].tolist()
    vecs = []
    for grp in grps:
        v = np.array([summ.loc[summ["group"] == grp, "a_hat_mean"].iloc[0],
                      summ.loc[summ["group"] == grp, "b_hat_mean"].iloc[0]])
        n = np.linalg.norm(v)
        vecs.append(v / n if n > 0 else v)
    vecs = np.array(vecs)
    print("\n  Pairwise cosines of group-mean linear direction (a_hat, b_hat):")
    print(" " * 20 + "  ".join(f"{str(g)[:14]:>14s}" for g in grps))
    for i, gi in enumerate(grps):
        row = "  ".join(f"{float(vecs[i] @ vecs[j]):+14.3f}" for j in range(len(grps)))
        print(f"  {str(gi)[:18]:>18s}  {row}")

    # MI-mean separation between groups (effect size = (mean1-mean2)/pooled_std)
    print("\n  Pairwise MI-mean separation (Cohen's d-like, |dmean|/pooled_std):")
    print(" " * 20 + "  ".join(f"{str(g)[:14]:>14s}" for g in grps))
    for i, gi in enumerate(grps):
        row = []
        for j, gj in enumerate(grps):
            si = df[df[strat_col] == gi]["MI_mean"]
            sj = df[df[strat_col] == gj]["MI_mean"]
            pooled = np.sqrt((si.var() * (len(si) - 1) + sj.var() * (len(sj) - 1)) /
                             max(1, len(si) + len(sj) - 2))
            d = abs(si.mean() - sj.mean()) / pooled if pooled > 1e-9 else 0.0
            row.append(f"{d:14.3f}")
        print(f"  {str(gi)[:18]:>18s}  {'  '.join(row)}")


for strat in ["strat1", "strat2", "strat3"]:
    summarize(strat)

# Save the joined table for downstream
df.to_csv("cerebro_per_subject_with_groups.csv", index=False)
print("\nSaved cerebro_per_subject_with_groups.csv")
