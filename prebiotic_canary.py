"""
Prebiotic-chemistry canary on the per-domain stack.

Question: when chemistry crosses into autocatalytic-mutual-catalysis
("hypercycle" territory, Eigen-Schuster), do its operator-space
signatures move toward biology, stay on chemistry, or land on a
distinct direction?

Speaking posture (Rule C): I think we might see a transition between
the Session 7 chemistry signature (MI ~ linear in H_a) and the
biology signature (MI ~ 0.5 * exp(H_a/2)) as the catalytic coupling
k_cat increases. But we'll wait on the data and where it points.
No verdict in advance, no verdict on first look.

System: 2-species Eigen hypercycle in a CSTR-like outflow constraint.
A and B replicate autocatalytically with mutual catalysis:

  dA/dt = A * (f_A - phi)         f_A = base_A + k_cat * B
  dB/dt = B * (f_B - phi)         f_B = base_B + k_cat * A
  phi    = (A*f_A + B*f_B) / (A + B)   (keeps A + B near constant)

  k_cat = 0   -> pure competing autocatalysts (selection / chemistry)
  k_cat > 0   -> mutualistic cycle (closer to biology)

Stochastic forcing on each species as in Session 6-8 simulators.

Per-domain stack identical to Session 6 / per_domain_kbk.py:
  ensemble H_a(t), H_b(t), MI(t) -> 6-column operator matrix
  -> KBK null direction
  -> polynomial fit of MI vs (H_a, H_b) for INFO-025-family comparison

Canary: one seed, N_ens=300, T=20, dt=0.02, k_cat sweep across
[0.0, 0.25, 0.5, 1.0, 2.0]. ~5 conditions, expect ~1-2 min total.
Output: prebiotic_canary.json
"""
from __future__ import annotations
import json
import time
import numpy as np

from per_domain_kbk import (
    entropy_hist_1d,
    build_ensemble_operator_matrix,
)
from kbk_pipeline import mi_hist_2d, extract_v1, project_234


def simulate_hypercycle(
    seed: int,
    k_cat: float,
    N_ens: int = 300,
    T: float = 20.0,
    dt: float = 0.02,
    base_A: float = 1.0,
    base_B: float = 1.0,
    noise: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """Eigen 2-species hypercycle with mutual catalysis k_cat.

    Returns (X_A, X_B) each of shape (steps, N_ens).
    """
    rng = np.random.default_rng(seed)
    steps = int(T / dt)
    A = 0.5 + 0.05 * rng.standard_normal(N_ens)
    B = 0.5 + 0.05 * rng.standard_normal(N_ens)
    A = np.clip(A, 0.05, None)
    B = np.clip(B, 0.05, None)
    XA = np.zeros((steps, N_ens))
    XB = np.zeros((steps, N_ens))
    for t in range(steps):
        fA = base_A + k_cat * B
        fB = base_B + k_cat * A
        tot = A + B
        phi = (A * fA + B * fB) / tot
        dA = A * (fA - phi) * dt
        dB = B * (fB - phi) * dt
        A = A + dA + noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        B = B + dB + noise * rng.standard_normal(N_ens) * np.sqrt(dt)
        A = np.clip(A, 0.01, None)
        B = np.clip(B, 0.01, None)
        XA[t] = A
        XB[t] = B
    return XA, XB


def fit_poly_and_exp(M: np.ndarray) -> dict:
    """Polynomial + exp(H_a/B) family fits to MI vs (H_a, H_b)."""
    Ha = M[:, 0]
    Hb = M[:, 1]
    MI = M[:, 5]
    out = {}

    # poly1 in (H_a, H_b)
    A1 = np.column_stack([np.ones_like(Ha), Ha, Hb])
    c1, *_ = np.linalg.lstsq(A1, MI, rcond=None)
    pred1 = A1 @ c1
    r2_1 = 1.0 - np.var(MI - pred1) / (np.var(MI) + 1e-12)
    out["poly1"] = {"coef": c1.tolist(), "r2": float(r2_1)}

    # poly2 in (H_a, H_b)
    A2 = np.column_stack([
        np.ones_like(Ha), Ha, Hb, Ha * Ha, Hb * Hb, Ha * Hb,
    ])
    c2, *_ = np.linalg.lstsq(A2, MI, rcond=None)
    pred2 = A2 @ c2
    r2_2 = 1.0 - np.var(MI - pred2) / (np.var(MI) + 1e-12)
    out["poly2"] = {"coef": c2.tolist(), "r2": float(r2_2)}

    # Session 7 biology family: MI = A * exp(H_a / B). Fit log(MI) ~ a0 + a1*H_a.
    mi_pos = MI[MI > 1e-6]
    ha_pos = Ha[MI > 1e-6]
    if len(mi_pos) > 10:
        Aexp = np.column_stack([np.ones_like(ha_pos), ha_pos])
        cexp, *_ = np.linalg.lstsq(Aexp, np.log(mi_pos), rcond=None)
        a0, a1 = cexp
        A_amp = float(np.exp(a0))
        B_scale = float(1.0 / a1) if abs(a1) > 1e-9 else float("inf")
        pred_exp = A_amp * np.exp(Ha / B_scale) if np.isfinite(B_scale) else None
        if pred_exp is not None:
            r2_exp = 1.0 - np.var(MI - pred_exp) / (np.var(MI) + 1e-12)
        else:
            r2_exp = float("nan")
        out["exp_INFO025"] = {"A": A_amp, "B": B_scale, "r2": float(r2_exp)}
    else:
        out["exp_INFO025"] = {"A": float("nan"), "B": float("nan"),
                              "r2": float("nan")}

    # Session 7 physics family: MI = c + (H_b - H_a)^2.
    Aphys = np.column_stack([np.ones_like(Ha), (Hb - Ha) ** 2])
    cphys, *_ = np.linalg.lstsq(Aphys, MI, rcond=None)
    pred_phys = Aphys @ cphys
    r2_phys = 1.0 - np.var(MI - pred_phys) / (np.var(MI) + 1e-12)
    out["phys_diff_sq"] = {"coef": cphys.tolist(), "r2": float(r2_phys)}

    return out


SESSION3_ATTRACTOR_234 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)


def cos_to_attractor(v6: np.ndarray) -> tuple[float, list]:
    """Cosine of |sub| to Session 3 attractor in [2,3,4] subspace.

    project_234 returns (sub_list, cos_to_plus_plus_plus). We sign-flip
    the absolute value so symmetric +++ and --- both register as on
    the attractor (Session 8 noted the sign-flipped twin)."""
    sub_list, _cos_plus = project_234(v6)
    if _cos_plus is None:
        return float("nan"), sub_list
    sub = np.asarray(sub_list)
    return float(abs(np.dot(sub, SESSION3_ATTRACTOR_234))), sub_list


def run_condition(k_cat: float, seed: int = 11,
                  N_ens: int = 300, T: float = 20.0,
                  dt: float = 0.02) -> dict:
    t0 = time.time()
    XA, XB = simulate_hypercycle(seed, k_cat, N_ens=N_ens, T=T, dt=dt)
    M, t_idx = build_ensemble_operator_matrix(XA, XB)
    v_null, S, _, _ = extract_v1(M)
    fits = fit_poly_and_exp(M)
    eigs = (S ** 2) / max(1, len(M) - 1)
    cos_attr, sub234 = cos_to_attractor(v_null)
    return {
        "k_cat": k_cat,
        "seed": seed,
        "N_ens": N_ens, "T": T, "dt": dt,
        "M_shape": list(M.shape),
        "v_null6": v_null.tolist(),
        "v_null_234_normalized": sub234,
        "cos_to_session3_attractor": cos_attr,
        "eigs_top5": eigs[:5].tolist(),
        "eigs_bottom3": eigs[-3:].tolist(),
        "MI_mean": float(M[:, 5].mean()),
        "MI_std": float(M[:, 5].std()),
        "Ha_mean": float(M[:, 0].mean()),
        "Hb_mean": float(M[:, 1].mean()),
        "fits": fits,
        "elapsed_s": float(time.time() - t0),
    }


def main():
    t_start = time.time()
    results = {
        "speaking_posture_before": (
            "I think the prebiotic system might show a transition "
            "between Session 7 chemistry and biology signatures as "
            "k_cat increases. We'll wait on the data and look together."
        ),
        "reference_signatures_session7": {
            "chemistry_brussel": "MI = 0.73 * H_a + 1.075  (linear in H_a)",
            "biology_lotka":     "MI = 0.50 * exp(H_a / 2.0)",
            "physics_duffing":   "MI = (H_b - H_a)^2 + 0.275",
            "geology_bk":        "MI = 0.199 (constant)",
        },
        "session3_attractor_234_subspace": SESSION3_ATTRACTOR_234.tolist(),
        "conditions": [],
    }
    k_cat_grid = [0.0, 0.25, 0.5, 1.0, 2.0]
    for k in k_cat_grid:
        print(f"[prebiotic] k_cat={k}", flush=True)
        rec = run_condition(k_cat=k)
        results["conditions"].append(rec)
        print(
            f"  done in {rec['elapsed_s']:.1f}s.  "
            f"MI_mean={rec['MI_mean']:.3f}  "
            f"poly1_r2={rec['fits']['poly1']['r2']:.3f}  "
            f"poly2_r2={rec['fits']['poly2']['r2']:.3f}  "
            f"exp_INFO025_r2={rec['fits']['exp_INFO025']['r2']:.3f}  "
            f"phys_diff_sq_r2={rec['fits']['phys_diff_sq']['r2']:.3f}  "
            f"cos_attr={rec['cos_to_session3_attractor']:.3f}",
            flush=True,
        )
    results["total_elapsed_s"] = float(time.time() - t_start)
    with open("prebiotic_canary.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nTotal: {results['total_elapsed_s']:.1f}s")
    print("Wrote prebiotic_canary.json")


if __name__ == "__main__":
    main()
