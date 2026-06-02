"""Render the four INFO-025 per-domain functional families."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa


# H range covers the typical operating envelope across Sessions 6-10
Ha_grid = np.linspace(-0.5, 3.0, 80)
Hb_grid = np.linspace(-0.5, 3.0, 80)
Ha, Hb = np.meshgrid(Ha_grid, Hb_grid)


def mi_physics(Ha, Hb):
    """Physics Duffing: symmetric quadratic in difference (2D)."""
    return (Hb - Ha) ** 2 + 0.275


def mi_biology(Ha, Hb):
    """Biology Lotka-Volterra: exponential in H_a, H_b absent (1D in H_a)."""
    return 0.50 * np.exp(Ha / 2.0) + 0 * Hb


def mi_chemistry(Ha, Hb):
    """Chemistry Brusselator: linear in H_a, H_b absent (1D in H_a)."""
    return 0.73 * Ha + 1.075 + 0 * Hb


def mi_geology(Ha, Hb):
    """Geology Burridge-Knopoff: constant (0D, decoupled)."""
    return np.full_like(Ha, 0.199)


FAMILIES = [
    ("Physics Duffing",   "MI = (H_b - H_a)^2 + 0.275",
     "symmetric quadratic in difference", 2, mi_physics, "viridis"),
    ("Biology Lotka-Volterra", "MI = 0.50 * exp(H_a / 2)",
     "exponential ramp in H_a; H_b absent", 1, mi_biology, "plasma"),
    ("Chemistry Brusselator", "MI = 0.73 * H_a + 1.075",
     "tilted plane in H_a; H_b absent", 1, mi_chemistry, "cividis"),
    ("Geology Burridge-Knopoff", "MI = 0.199",
     "flat plane; decoupled (no H dependence)", 0, mi_geology, "magma"),
]


def render_2x2():
    fig = plt.figure(figsize=(14, 11))
    for i, (name, eq, shape, dims, fn, cmap) in enumerate(FAMILIES):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        MI = fn(Ha, Hb)
        ax.plot_surface(Ha, Hb, MI, cmap=cmap, alpha=0.85,
                        linewidth=0, antialiased=True)
        ax.set_xlabel("H_a")
        ax.set_ylabel("H_b")
        ax.set_zlabel("MI")
        ax.set_title(f"{name}\n{eq}\nshape: {shape}\n"
                     f"effective dims: {dims}",
                     fontsize=10)
        ax.view_init(elev=25, azim=-130)
    fig.suptitle("INFO-025: Four reproducible per-domain MI-vs-H families "
                 "(Session 7)", fontsize=13, y=0.995)
    plt.tight_layout()
    plt.savefig("info025_four_families.png", dpi=140, bbox_inches="tight")
    plt.close()
    print("wrote info025_four_families.png")


def render_individual():
    for name, eq, shape, dims, fn, cmap in FAMILIES:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
        MI = fn(Ha, Hb)
        ax.plot_surface(Ha, Hb, MI, cmap=cmap, alpha=0.9,
                        linewidth=0, antialiased=True)
        ax.set_xlabel("H_a")
        ax.set_ylabel("H_b")
        ax.set_zlabel("MI")
        ax.set_title(f"{name}\n{eq}\nshape: {shape}  |  effective dims: {dims}",
                     fontsize=11)
        ax.view_init(elev=25, azim=-130)
        slug = name.split()[0].lower()
        path = f"info025_{slug}.png"
        plt.tight_layout()
        plt.savefig(path, dpi=140, bbox_inches="tight")
        plt.close()
        print(f"wrote {path}")


if __name__ == "__main__":
    render_2x2()
    render_individual()
