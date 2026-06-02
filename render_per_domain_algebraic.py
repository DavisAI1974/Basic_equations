"""Render the four per-domain algebraic constraint surfaces.

These are the INFO-008 line of equations (chem / bio / geo / phy)
reproduced and located by Session 6 per-domain ensemble-H + KBK probe.
They are constraints of the form f(operators) approximately 0 that
the per-domain ensemble trajectories satisfy. Each constraint defines
a surface in the 6D operator basis [H_a, H_b, H_a^2, H_b^2,
H_a*H_b, MI].

We render each constraint residual r(H_a, H_b) = f(operators) over
the (H_a, H_b) plane. Where r is near zero, the constraint is
tightly satisfied. The shape of r over (H_a, H_b) reveals the
geometry of the per-domain law.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa


Ha_grid = np.linspace(-0.5, 3.0, 80)
Hb_grid = np.linspace(-0.5, 3.0, 80)
Ha, Hb = np.meshgrid(Ha_grid, Hb_grid)


def physics_residual(Ha, Hb):
    """Physics Duffing: Taylor quadratic identity.
    Constraint: H_a^2 + H_b^2 - 2*H_a*H_b approximately 0
    (the universal Session 3 attractor identity, equivalent to (H_a - H_b)^2 approximately 0).
    Residual surface is a parabolic bowl with zero-valley along H_a = H_b.
    """
    return Ha ** 2 + Hb ** 2 - 2 * Ha * Hb


def chemistry_residual(Ha, Hb):
    """Chemistry Brusselator: KBK null[0] from Session 6.
    Constraint: -0.30*H_a + 0.78*H_b - 0.43*H_b^2 + 0.33*H_a*H_b approximately 0
    Full mixed linear-plus-quadratic shape, breaks the symmetry of physics
    (the H_a-H_b non-equality of Brusselator dynamics).
    """
    return -0.30 * Ha + 0.78 * Hb - 0.43 * Hb ** 2 + 0.33 * Ha * Hb


def biology_residual(Ha, Hb):
    """Biology Lotka-Volterra: KBK smallest null direction from Session 6.
    Constraint (in operator basis): -0.27*H_a + 0.96*MI approximately 0; H_b absent.
    Equivalent reading: MI approximately 0.28 * H_a (1D structure in H_a).
    Since MI is not in (H_a, H_b) plane, we render the residual using
    the equivalent algebraic form: H_a*0.27 approximately MI, so
    residual is -0.27*H_a over the (H_a, H_b) plane (zero only at H_a=0
    when MI is set to 0); the H_b-absent structure shows as flatness
    along H_b.
    """
    # The residual against the actual MI is zero everywhere if the
    # constraint holds. Visualize the linear-in-H_a slope alone:
    return -0.27 * Ha + 0 * Hb


def geology_residual(Ha, Hb):
    """Geology Burridge-Knopoff: rank-3 null relation.
    Constraint: 0.724*(H_a*H_b) - 0.441*H_b^2 - 0.290*H_a^2 approximately 0
    Quadratic form in (H_a^2, H_b^2, H_a*H_b) only. No linear terms.
    Cross-seed cos +0.99 to v5 coefficients.
    """
    return 0.724 * (Ha * Hb) - 0.441 * Hb ** 2 - 0.290 * Ha ** 2


FAMILIES = [
    ("Physics Duffing",
     "H_a^2 + H_b^2 - 2*H_a*H_b ~ 0",
     "parabolic bowl; zero-valley along the H_a = H_b diagonal",
     "2 H-channels; symmetric quadratic",
     physics_residual, "viridis"),
    ("Biology Lotka-Volterra",
     "-0.27*H_a + 0.96*MI ~ 0  (H_b absent)",
     "linear ramp in H_a; flat in H_b (H_b-absent structure)",
     "1 H-channel (H_a only)",
     biology_residual, "plasma"),
    ("Chemistry Brusselator",
     "-0.30*H_a + 0.78*H_b - 0.43*H_b^2 + 0.33*H_a*H_b ~ 0",
     "mixed linear-plus-quadratic saddle; symmetry-broken",
     "2 H-channels; full mixed linear + quadratic + cubic (SINDy)",
     chemistry_residual, "cividis"),
    ("Geology Burridge-Knopoff",
     "0.724*H_a*H_b - 0.441*H_b^2 - 0.290*H_a^2 ~ 0",
     "quadratic-only saddle; zero-set is a conic in (H_a, H_b)",
     "2 H-channels; pure quadratic forms (no linear terms); rank-3 surface",
     geology_residual, "magma"),
]


def render_2x2():
    fig = plt.figure(figsize=(16, 12))
    for i, (name, eq, shape, dims, fn, cmap) in enumerate(FAMILIES):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        R = fn(Ha, Hb)
        ax.plot_surface(Ha, Hb, R, cmap=cmap, alpha=0.85,
                        linewidth=0, antialiased=True)
        # Draw the zero-plane for reference
        Z0 = np.zeros_like(R)
        ax.plot_surface(Ha, Hb, Z0, color="lightgray", alpha=0.15,
                        linewidth=0, antialiased=False)
        ax.set_xlabel("H_a")
        ax.set_ylabel("H_b")
        ax.set_zlabel("constraint residual")
        ax.set_title(f"{name}\n{eq}\nshape: {shape}\ndim signature: {dims}",
                     fontsize=10)
        ax.view_init(elev=25, azim=-130)
    fig.suptitle("Per-domain algebraic constraint surfaces "
                 "(INFO-008 line, located by Session 6 stack)",
                 fontsize=13, y=0.995)
    plt.tight_layout()
    plt.savefig("per_domain_algebraic_four.png", dpi=140,
                bbox_inches="tight")
    plt.close()
    print("wrote per_domain_algebraic_four.png")


def render_individual():
    for name, eq, shape, dims, fn, cmap in FAMILIES:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
        R = fn(Ha, Hb)
        ax.plot_surface(Ha, Hb, R, cmap=cmap, alpha=0.9,
                        linewidth=0, antialiased=True)
        Z0 = np.zeros_like(R)
        ax.plot_surface(Ha, Hb, Z0, color="lightgray", alpha=0.2,
                        linewidth=0, antialiased=False)
        ax.set_xlabel("H_a")
        ax.set_ylabel("H_b")
        ax.set_zlabel("constraint residual")
        ax.set_title(f"{name}\n{eq}\n{shape}",
                     fontsize=11)
        ax.view_init(elev=25, azim=-130)
        slug = name.split()[0].lower()
        path = f"per_domain_algebraic_{slug}.png"
        plt.tight_layout()
        plt.savefig(path, dpi=140, bbox_inches="tight")
        plt.close()
        print(f"wrote {path}")


if __name__ == "__main__":
    render_2x2()
    render_individual()
