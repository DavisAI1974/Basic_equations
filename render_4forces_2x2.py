"""
Render: 4 forces on the same substrate, 2x2 grid.

Each panel shows ONE force against its own preferred variable with
its PySR low-complexity Pareto family fit overlaid. The four panels
share the same substrate (same KBK null direction in operator
coords) but live in different expression-axes.

Panels (Session 8 INFO-027, seed 11, N_ens=600, T=30):

  EM:      MI vs (H_a - H_b)    fit = (H_a - H_b)^2 + 0.203
  Weak:    MI vs H_a            fit = 0.305 + 0.063 * H_a
  Strong:  MI vs H_a (sym dyn)  fit = 0.227 * exp(H_a^2)
  Gravity: MI vs H_b (sym dyn)  fit = 0.336 - 0.572 * H_b^2

Strong and gravity are channel-symmetric so seed 22 might pick the
other channel; we show seed 11 to keep one consistent slice. The
post-close em_strong_glance (INFO-031) showed that making the
dynamics asymmetric flips strong to exp((H_b - H_a) - 1.4),
revealing the underlying difference-variable expression.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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


def regenerate(seed=11):
    cfg = dict(N_ens=600, T=30.0, dt=0.02)
    out = {}
    sim = {"em": simulate_em, "weak": simulate_weak,
           "strong": simulate_strong, "gravity": simulate_gravity}
    for name, fn in sim.items():
        X1, X2 = fn(seed=seed, **cfg)
        M, _ = build_ensemble_operator_matrix(X1, X2)
        out[name] = {"Ha": M[:, 0], "Hb": M[:, 1], "MI": M[:, 5]}
    return out


def main():
    data = regenerate(seed=11)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # --- EM ---
    ax = axes[0, 0]
    Ha = data["em"]["Ha"]; Hb = data["em"]["Hb"]; MI = data["em"]["MI"]
    diff = Ha - Hb
    ax.scatter(diff, MI, s=10, alpha=0.4, color=FORCE_COLORS["em"])
    xline = np.linspace(diff.min(), diff.max(), 200)
    yfit = xline ** 2 + 0.203
    ax.plot(xline, yfit, color="black", linewidth=2.2, alpha=0.85,
             label=r"PySR: $(H_a - H_b)^2 + 0.203$  (c=6)")
    ax.set_title("EM  -- (asymmetric: omega_1=1.0, omega_2=1.2)",
                  fontsize=11)
    ax.set_xlabel(r"$H_a - H_b$  (entropy difference)")
    ax.set_ylabel(r"MI$(X_1, X_2)$")
    ax.legend(loc="upper center", fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Weak ---
    ax = axes[0, 1]
    Ha = data["weak"]["Ha"]; MI = data["weak"]["MI"]
    ax.scatter(Ha, MI, s=10, alpha=0.4, color=FORCE_COLORS["weak"])
    xline = np.linspace(Ha.min(), Ha.max(), 200)
    yfit = 0.305 + 0.063 * xline
    ax.plot(xline, yfit, color="black", linewidth=2.2, alpha=0.85,
             label=r"PySR: $0.305 + 0.063 \cdot H_a$  (c=5)")
    ax.set_title("Weak  -- (Yukawa suppression, strong damping)",
                  fontsize=11)
    ax.set_xlabel(r"$H_a$")
    ax.set_ylabel(r"MI$(X_1, X_2)$")
    ax.legend(loc="upper left", fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Strong (seed 11) ---
    ax = axes[1, 0]
    Ha = data["strong"]["Ha"]; MI = data["strong"]["MI"]
    ax.scatter(Ha, MI, s=10, alpha=0.4, color=FORCE_COLORS["strong"])
    xline = np.linspace(Ha.min(), Ha.max(), 200)
    yfit = 0.227 * np.exp(xline ** 2)
    ax.plot(xline, yfit, color="black", linewidth=2.2, alpha=0.85,
             label=r"PySR: $0.227 \cdot \exp(H_a^2)$  (c=5)")
    ax.set_title("Strong  -- (channel-symmetric; seed 11 picked $H_a$)",
                  fontsize=11)
    ax.set_xlabel(r"$H_a$")
    ax.set_ylabel(r"MI$(X_1, X_2)$")
    ax.legend(loc="upper center", fontsize=10)
    ax.grid(True, alpha=0.3)

    # --- Gravity (seed 11) ---
    ax = axes[1, 1]
    Hb = data["gravity"]["Hb"]; MI = data["gravity"]["MI"]
    ax.scatter(Hb, MI, s=10, alpha=0.4, color=FORCE_COLORS["gravity"])
    xline = np.linspace(Hb.min(), Hb.max(), 200)
    yfit = 0.336 - 0.572 * xline ** 2
    ax.plot(xline, yfit, color="black", linewidth=2.2, alpha=0.85,
             label=r"PySR: $0.336 - 0.572 \cdot H_b^2$  (c=6)")
    ax.set_title("Gravity  -- (channel-symmetric; seed 11 picked $H_b$)",
                  fontsize=11)
    ax.set_xlabel(r"$H_b$")
    ax.set_ylabel(r"MI$(X_1, X_2)$")
    ax.legend(loc="lower center", fontsize=10)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Same substrate, distinct expression  "
        "(Session 8 INFO-027, seed 11, $N_{\\mathrm{ens}}{=}600$, $T{=}30$)\n"
        r"All four sit on the (-1,-1,+2)/$\sqrt{6}$ null direction in $(H_a^2, H_b^2, H_a \cdot H_b)$; "
        "PySR low-complexity Pareto family per panel",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    fig.savefig("render_4forces_2x2.png", dpi=160, bbox_inches="tight")
    print("wrote render_4forces_2x2.png")


if __name__ == "__main__":
    main()
