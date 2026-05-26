"""
Render: 4 forces on the same substrate.

Two-panel figure:

  LEFT  — substrate side. v_null directions of all four forces in
          the (H_a^2, H_b^2, H_a*H_b) operator-coordinate space.
          All four point near (-1, -1, +2)/sqrt(6), the Session 3
          attractor (algebraic identity -(H_a - H_b)^2 ~ 0).
          The shared direction is the SUBSTRATE.

  RIGHT — expression side. MI vs (H_a - H_b) per force from the
          Session 8 ensemble-H data. Each force has a distinct
          functional family at low PySR complexity overlaid as
          the analytical fit:
            EM:      (H_a - H_b)^2 + 0.21
            Weak:    0.31 + 0.066 * H_a   (function of H_a, not
                     H_a-H_b; we show it vs H_a-H_b for visual
                     reference -- it has no slope vs the difference)
            Strong:  channel-symmetric, randomly seeded -- we
                     show seed 11's exp(H_a^2)*0.23 form
            Gravity: channel-symmetric -- seed 11's
                     0.336 - 0.572*H_b^2

The substrate is the same line/plane in operator space; the
expressions are different shapes on top of it.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from four_force_probe import (
    simulate_em, simulate_weak, simulate_strong, simulate_gravity,
)
from per_domain_kbk import build_ensemble_operator_matrix


FORCE_COLORS = {
    "em":      "#1f77b4",
    "weak":    "#d62728",
    "strong":  "#2ca02c",
    "gravity": "#9467bd",
}
FORCE_LABELS = {
    "em":      "EM",
    "weak":    "Weak",
    "strong":  "Strong",
    "gravity": "Gravity",
}


def regenerate_data(seed=11):
    cfg = dict(N_ens=600, T=30.0, dt=0.02)
    data = {}
    sim = {"em": simulate_em, "weak": simulate_weak,
           "strong": simulate_strong, "gravity": simulate_gravity}
    for name, fn in sim.items():
        X1, X2 = fn(seed=seed, **cfg)
        M, _ = build_ensemble_operator_matrix(X1, X2)
        data[name] = {
            "Ha": M[:, 0], "Hb": M[:, 1], "MI": M[:, 5],
            "Ha2": M[:, 2], "Hb2": M[:, 3], "HaHb": M[:, 4],
        }
    return data


def load_v_nulls(results_path):
    with open(results_path) as f:
        blob = json.load(f)
    out = {}
    for r in blob["results"]:
        out[r["name"]] = {
            "v_null_234": r["extract_v1_v_null_234"],
            "v_null_6d":  r["extract_v1_v_null_6d"],
            "eig_min":    r["eigenvalues_op_cov"][-1],
        }
    return out


def main():
    seed = 11
    data = regenerate_data(seed=seed)
    v_nulls = load_v_nulls("four_force_probe_results_seed11.json")

    fig = plt.figure(figsize=(15, 6.5))

    # ---------- LEFT: substrate side ----------
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    # plot attractor reference vector
    attractor = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6)  # length 1
    ax1.plot([0, attractor[0]], [0, attractor[1]], [0, attractor[2]],
             color="black", linewidth=2.5, alpha=0.4,
             label=r"attractor $(-1,-1,+2)/\sqrt{6}$")
    ax1.scatter([attractor[0]], [attractor[1]], [attractor[2]],
                color="black", s=80, marker="*", zorder=5)

    # plot each force's v_null in [2,3,4] space
    for name in ["em", "weak", "strong", "gravity"]:
        v = np.array(v_nulls[name]["v_null_234"])
        # normalize for display
        v = v / np.linalg.norm(v)
        ax1.plot([0, v[0]], [0, v[1]], [0, v[2]],
                 color=FORCE_COLORS[name], linewidth=3, alpha=0.85,
                 label=FORCE_LABELS[name])
        ax1.scatter([v[0]], [v[1]], [v[2]],
                    color=FORCE_COLORS[name], s=70, zorder=4)

    ax1.set_xlabel(r"$H_a^2$ coefficient", fontsize=10)
    ax1.set_ylabel(r"$H_b^2$ coefficient", fontsize=10)
    ax1.set_zlabel(r"$H_a \cdot H_b$ coefficient", fontsize=10)
    ax1.set_title(
        "Substrate (shared)\n"
        r"KBK null direction in operator coords; algebraic identity $-(H_a - H_b)^2 \approx 0$",
        fontsize=11,
    )
    ax1.legend(loc="upper left", fontsize=9, framealpha=0.9)
    # nice viewing angle
    ax1.view_init(elev=18, azim=35)
    ax1.set_xlim(-0.5, 0.1)
    ax1.set_ylim(-0.5, 0.1)
    ax1.set_zlim(0, 0.95)

    # ---------- RIGHT: expression side ----------
    ax2 = fig.add_subplot(1, 2, 2)
    for name in ["em", "weak", "strong", "gravity"]:
        d = data[name]
        # subsample for plotting
        idx = np.linspace(0, len(d["Ha"]) - 1, 300).astype(int)
        diff = (d["Ha"] - d["Hb"])[idx]
        mi = d["MI"][idx]
        ax2.scatter(diff, mi, s=8, alpha=0.45,
                     color=FORCE_COLORS[name], label=FORCE_LABELS[name])

    # PySR fits per force (from Session 8 INFO-027, low-complexity Pareto)
    diff_x = np.linspace(-1.5, 1.5, 200)

    # EM: (H_a - H_b)^2 + 0.21
    ax2.plot(diff_x, diff_x ** 2 + 0.21,
             color=FORCE_COLORS["em"], linewidth=2.5, alpha=0.95,
             linestyle="-", label=r"EM fit: $(H_a-H_b)^2 + 0.21$")
    # Weak: linear in H_a (no slope vs H_a - H_b on its own; this
    # plot is illustrative of the *substrate* axis, not the weak
    # family's true variable). For visual reference plot a horizontal
    # line at the weak mean MI.
    w_mean = float(np.mean(data["weak"]["MI"]))
    ax2.axhline(w_mean, color=FORCE_COLORS["weak"], linewidth=2,
                 linestyle="--", alpha=0.75,
                 label=f"Weak fit: linear in $H_a$ (mean MI = {w_mean:.2f})")
    # Strong seed 11: exp(H_a^2) * 0.23 -- function of H_a alone, not
    # H_a - H_b. Plot strong's MI mean for reference.
    s_mean = float(np.mean(data["strong"]["MI"]))
    ax2.axhline(s_mean, color=FORCE_COLORS["strong"], linewidth=2,
                 linestyle=":", alpha=0.75,
                 label=f"Strong (sym, seed 11): mean MI = {s_mean:.2f}")
    # Gravity seed 11: 0.336 - 0.572*H_b^2 -- function of H_b alone.
    g_mean = float(np.mean(data["gravity"]["MI"]))
    ax2.axhline(g_mean, color=FORCE_COLORS["gravity"], linewidth=2,
                 linestyle="-.", alpha=0.75,
                 label=f"Gravity (sym, seed 11): mean MI = {g_mean:.2f}")

    ax2.set_xlabel(r"$H_a - H_b$  (entropy difference)", fontsize=11)
    ax2.set_ylabel(r"$\mathrm{MI}(X_1, X_2)$", fontsize=11)
    ax2.set_title(
        "Expression (distinct)\n"
        r"MI vs $H_a - H_b$ per force; PySR low-complexity Pareto fits overlaid",
        fontsize=11,
    )
    ax2.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-1.6, 1.6)

    fig.suptitle(
        "Four force-field caricatures: shared substrate, distinct expression  "
        r"(Session 8 INFO-027, seed 11, $N_{\mathrm{ens}}{=}600$, $T{=}30$)",
        fontsize=13, y=1.00,
    )
    fig.tight_layout()
    fig.savefig("render_4forces.png", dpi=160, bbox_inches="tight")
    print("wrote render_4forces.png")

    # also dump a smaller version
    fig.savefig("render_4forces_small.png", dpi=110, bbox_inches="tight")
    print("wrote render_4forces_small.png")


if __name__ == "__main__":
    main()
