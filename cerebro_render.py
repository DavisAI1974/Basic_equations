"""Visualize the cerebro probe: per-subject linear direction + MI-mean by group, for each stratification."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("cerebro_per_subject_with_groups.csv")

# Cardiac canary (round 1+2) for cross-organ visual
with open("ekg_disease_canary.json") as f:
    cardiac1 = json.load(f)
with open("ekg_disease_canary_expand.json") as f:
    cardiac2 = json.load(f)
cardiac = [r for r in cardiac1 + cardiac2 if r.get("ok")]

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

colors_s1 = {"control": "#2c7fb8", "DM_only": "#41b6c4",
             "DM_hypertension": "#fdae61", "DM_orthostatic": "#d73027"}
colors_s3 = {"control": "#2c7fb8", "DM_uncomplicated": "#41b6c4",
             "DM_retinopathy": "#fec44f", "DM_neuropathy": "#fe9929",
             "DM_nephropathy_proxy": "#d95f0e"}

for col, strat, colors in [(0, "strat1", colors_s1),
                            (1, "strat2", None),
                            (2, "strat3", colors_s3)]:
    # Top row: linear-direction scatter
    ax = axes[0, col]
    grps = sorted(df[strat].dropna().unique())
    if colors is None:
        cmap = plt.cm.viridis
        colors = {g: cmap(i / max(1, len(grps) - 1)) for i, g in enumerate(grps)}
    for grp in grps:
        sub = df[df[strat] == grp]
        ax.scatter(sub["a_hat"], sub["b_hat"], s=80, alpha=0.7,
                   color=colors.get(grp, "gray"),
                   label=f"{grp} (n={len(sub)})", edgecolor="k", linewidth=0.5)
    # Draw unit circle
    t = np.linspace(0, np.pi / 2, 50)
    ax.plot(np.cos(t), np.sin(t), "k--", linewidth=0.5, alpha=0.3)
    ax.set_xlabel("a_hat (H_a direction)")
    ax.set_ylabel("b_hat (H_b direction)")
    ax.set_title(f"{strat}: per-subject linear direction in (H_a, H_b)")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_aspect("equal")
    ax.set_xlim(-0.1, 1.05)
    ax.set_ylim(-0.1, 1.05)
    ax.grid(alpha=0.3)

    # Bottom row: MI_mean by group (box+jitter)
    ax = axes[1, col]
    data = [df[df[strat] == g]["MI_mean"].values for g in grps]
    bp = ax.boxplot(data, labels=[str(g)[:14] for g in grps],
                    patch_artist=True, showmeans=True)
    for patch, g in zip(bp["boxes"], grps):
        patch.set_facecolor(colors.get(g, "gray"))
        patch.set_alpha(0.5)
    for g_i, g in enumerate(grps):
        y = df[df[strat] == g]["MI_mean"].values
        x = np.random.normal(g_i + 1, 0.05, size=len(y))
        ax.scatter(x, y, s=40, alpha=0.7, color="k", zorder=3)
    ax.set_ylabel("MI mean (HR ; ABP)")
    ax.set_title(f"{strat}: MI(HR;ABP) per subject, by group")
    ax.tick_params(axis="x", rotation=25)
    ax.grid(alpha=0.3, axis="y")

# Add cardiac canary at the bottom for visual comparison
# Replace bottom-right with cardiac comparison
plt.suptitle("Cerebrovascular probe on cerebral-vasoreg-diabetes cohort (44 subjects)\n"
             "Coupled pair: (HR from ECG, ABP)  -  Different sample type than cardiac HR/HRV",
             fontsize=14, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("cerebro_render.png", dpi=130, bbox_inches="tight")
print("Saved cerebro_render.png")

# Second plot: cross-organ comparison (cerebro + cardiac on same direction-space)
fig2, ax = plt.subplots(1, 1, figsize=(10, 9))
# Cerebro groups (use strat3 since it has most differentiation)
for grp in sorted(df["strat3"].dropna().unique()):
    sub = df[df["strat3"] == grp]
    ax.scatter(sub["a_hat"], sub["b_hat"], s=120, alpha=0.7,
               color=colors_s3.get(grp, "gray"),
               marker="o", label=f"cerebro: {grp} (n={len(sub)})",
               edgecolor="k", linewidth=0.5)
# Cardiac points (single-record per family)
cardiac_colors = {"healthy_young": "#2c7fb8", "healthy_elderly": "#3690c0",
                  "healthy_long_term": "#41b6c4",
                  "atrial_fibrillation": "#a50026",
                  "congestive_heart_failure": "#fc8d59",
                  "ischemia_ST": "#fdae61",
                  "supraventricular_arr": "#984ea3",
                  "ventricular_malignant": "#7a0177",
                  "arrhythmia_mixed": "#1b9e77",
                  "ltaf_long_term": "#d62728"}
for r in cardiac:
    p1 = r["fits"]["poly1"]["coef"]
    v = np.array([p1["a"], p1["b"]])
    n = np.linalg.norm(v)
    if n < 1e-9: continue
    ax.scatter(v[0]/n, v[1]/n, s=200, alpha=0.85,
               color=cardiac_colors.get(r["family"], "gray"),
               marker="^", edgecolor="k", linewidth=1.0,
               label=f"cardiac: {r['family']}")
t = np.linspace(0, np.pi / 2, 50)
ax.plot(np.cos(t), np.sin(t), "k--", linewidth=0.5, alpha=0.3)
ax.set_xlabel("a_hat (H_a direction)", fontsize=12)
ax.set_ylabel("b_hat (H_b direction)", fontsize=12)
ax.set_title("Cross-organ linear-direction scatter\n"
             "Cerebrovascular (HR;ABP, circles) + Cardiac (HR;HRV, triangles)\n"
             "Substrate cluster around (0.7, 0.7); AF and significant CKD pull toward pure-H_b",
             fontsize=12)
ax.legend(fontsize=8, loc="lower left", ncol=2)
ax.set_aspect("equal")
ax.set_xlim(-0.1, 1.05)
ax.set_ylim(-0.1, 1.05)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cross_organ_direction.png", dpi=130, bbox_inches="tight")
print("Saved cross_organ_direction.png")
