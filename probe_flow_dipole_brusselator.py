"""
*** SUPERSEDED / UNTRUSTWORTHY (Greg, S19) -- DO NOT BUILD ON THIS ***
Toy Brusselator ODE simulator (per_domain_kbk.simulate_chemistry). Greg S19:
"bruss might be wrong tool ... consider it untrustworthy." Replaced by a REAL
chemistry dataset in probe_flow_dipole_chem.py. The opposition signature seen here
is a toy-simulator artifact -- it did not appear in the real-data forces (gravity,
weak). Kept only for the audit trail.
*** end deprecation notice ***

FLOW DIPOLE EQUATION -- system 1: Brusselator (native 2 species x,y).
Backlog STEP 1, "START HERE, cleanest" (Greg, S18).

A flow was, until now, a single characteristic vs an axis (INFO-062/063).
A DIPOLE needs TWO coupled channels in the info-dipole paper's form
(https://davisai.ai/dipole/):

    dMI/dt ~ sum_i c_self,i * H_i^2 + sum_{i<j} c_cross,ij * H_i*H_j + linear
             (opposition signature: self vs cross opposing sign)

Here the two channels are the Brusselator's native species x (X1) and y (X2).
Per-time-step ensemble entropies H_a(t)=H(x), H_b(t)=H(y) and MI(t)=I(x;y) are
the SAME operators used for the per-domain chemistry results (INFO-023/025/040),
so this dipole equation is directly comparable to that prior work.

DELIBERATE DEFLATIONARY CONTROLS (so the frame does not grade itself):
  (B) 1-D MI-relaxation control: dMI/dt ~ a*MI + b. If a simple 1-D relaxation
      already explains dMI/dt, the 2-channel H_a/H_b dipole terms are NOT needed
      -> "dipole" would be a redescription of relaxation, not a coupling.
  (C) self-only model (no cross term): does H_a*H_b carry signal beyond the
      self-terms? The cross term is the actual coupling channel.
  (S) time-shuffle null on the FEATURES (keep dMI/dt aligned to its own time,
      permute the operator rows): destroys any genuine instantaneous
      operator->rate relation; R^2 should collapse to ~0 if the fit is real
      (guards against a spurious high R^2 from sheer column count).
  Raw AND per-column standardized coefficients are both reported -- the
  marginal entropies are units/binning dependent (INFO-051: MI is the only
  scale-invariant operator), so standardized betas are the unit-robust read.

Result Discipline: >=3 seeds, inter-seed scatter reported. No claim from one seed.
Run:  python probe_flow_dipole_brusselator.py            (full, seeds 11/22/33)
      python probe_flow_dipole_brusselator.py --canary    (1 seed, short)
"""
import numpy as np
import json
import time
import argparse

from per_domain_kbk import simulate_chemistry, build_ensemble_operator_matrix

# operator column order from build_ensemble_operator_matrix:
# [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]
COL = dict(Ha=0, Hb=1, Ha2=2, Hb2=3, HaHb=4, MI=5)

FEAT_NAMES = ["H_a^2", "H_b^2", "H_a*H_b", "H_a", "H_b", "1"]


def moving_average(v, w):
    if w <= 1:
        return v.copy()
    k = np.ones(w) / w
    return np.convolve(v, k, mode="same")


def ols(X, y):
    """Return coefficients, R^2."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2


def standardized_betas(Xcols, y):
    """Per-column z-scored design (no const), z-scored target -> comparable betas."""
    Xs = np.column_stack([(c - c.mean()) / (c.std() + 1e-12) for c in Xcols])
    ys = (y - y.mean()) / (y.std() + 1e-12)
    beta, r2 = ols(Xs, ys)
    return beta, r2


def analyze_seed(seed, T, N_ens, smooth_w):
    X1, X2 = simulate_chemistry(seed, N_ens=N_ens, T=T, dt=0.02)
    M, t_idx = build_ensemble_operator_matrix(X1, X2)
    dt = 0.02
    t = t_idx * dt

    Ha = M[:, COL["Ha"]]
    Hb = M[:, COL["Hb"]]
    MI = M[:, COL["MI"]]

    # smooth MI (histogram estimator noise) then differentiate vs real time
    MIs = moving_average(MI, smooth_w)
    dMI = np.gradient(MIs, t)

    # trim convolution edges
    edge = smooth_w
    sl = slice(edge, len(t) - edge)
    Ha, Hb, MI, dMI = Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha * Ha, Hb * Hb, Ha * Hb

    # ---- Model A: full flow-dipole form ----
    XA = np.column_stack([Ha2, Hb2, HaHb, Ha, Hb, np.ones_like(Ha)])
    betaA, r2A = ols(XA, dMI)
    # standardized betas for the 5 non-constant terms
    sbeta, _ = standardized_betas([Ha2, Hb2, HaHb, Ha, Hb], dMI)

    # ---- Model B: 1-D MI-relaxation control ----
    XB = np.column_stack([MI, np.ones_like(MI)])
    betaB, r2B = ols(XB, dMI)

    # ---- Model C: self-only (no cross term) ----
    XC = np.column_stack([Ha2, Hb2, Ha, Hb, np.ones_like(Ha)])
    _, r2C = ols(XC, dMI)

    # ---- Model relaxation+full check: full over relaxation ----
    dR2_dipole_over_relax = r2A - r2B
    dR2_cross = r2A - r2C  # marginal value of the H_a*H_b cross term

    # ---- time-shuffle null on the features ----
    rng = np.random.default_rng(seed + 999)
    r2_null = []
    for _ in range(20):
        perm = rng.permutation(len(dMI))
        Xn = np.column_stack([Ha2[perm], Hb2[perm], HaHb[perm],
                              Ha[perm], Hb[perm], np.ones_like(Ha)])
        _, r2n = ols(Xn, dMI)
        r2_null.append(r2n)
    r2_null = np.array(r2_null)

    # opposition signature in the QUADRATIC subspace (self vs cross)
    c_self = np.array([betaA[0], betaA[1]])      # H_a^2, H_b^2
    c_cross = betaA[2]                            # H_a*H_b
    opposition = bool(np.sign(c_cross) != 0 and
                      np.all(np.sign(c_self) == -np.sign(c_cross)))

    return {
        "seed": seed,
        "n_points": int(len(dMI)),
        "MI_range": [float(MI.min()), float(MI.max())],
        "dMI_abs_mean": float(np.mean(np.abs(dMI))),
        "model_A_full": {
            "coeffs_raw": {FEAT_NAMES[i]: float(betaA[i]) for i in range(6)},
            "coeffs_std": {FEAT_NAMES[i]: float(sbeta[i]) for i in range(5)},
            "R2": float(r2A),
        },
        "model_B_MI_relax": {"a": float(betaB[0]), "b": float(betaB[1]),
                             "R2": float(r2B)},
        "model_C_self_only_R2": float(r2C),
        "dR2_dipole_over_relax": float(dR2_dipole_over_relax),
        "dR2_cross_term": float(dR2_cross),
        "shuffle_null_R2_mean": float(r2_null.mean()),
        "shuffle_null_R2_max": float(r2_null.max()),
        "opposition_signature_quad": opposition,
    }


def aggregate(seed_results):
    """mean +/- std of standardized betas and key R^2 across seeds."""
    keys_std = FEAT_NAMES[:5]
    agg = {}
    for k in keys_std:
        vals = [s["model_A_full"]["coeffs_std"][k] for s in seed_results]
        agg[k] = {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
    for metric in ["R2"]:
        vals = [s["model_A_full"][metric] for s in seed_results]
        agg["full_" + metric] = {"mean": float(np.mean(vals)),
                                 "std": float(np.std(vals))}
    for fld in ["model_C_self_only_R2", "dR2_dipole_over_relax",
                "dR2_cross_term", "shuffle_null_R2_mean"]:
        vals = [s[fld] for s in seed_results]
        agg[fld] = {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
    relax = [s["model_B_MI_relax"]["R2"] for s in seed_results]
    agg["relax_R2"] = {"mean": float(np.mean(relax)), "std": float(np.std(relax))}
    opp = [s["opposition_signature_quad"] for s in seed_results]
    agg["opposition_fraction"] = float(np.mean(opp))
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    if args.canary:
        seeds, T, N_ens, smooth_w = [11], 15.0, 300, 5
        out_path = "probe_flow_dipole_brusselator_canary.json"
    else:
        seeds, T, N_ens, smooth_w = [11, 22, 33], 30.0, 600, 5
        out_path = "probe_flow_dipole_brusselator_results.json"

    seed_results = [analyze_seed(s, T, N_ens, smooth_w) for s in seeds]
    agg = aggregate(seed_results)

    result = {
        "probe": "flow-dipole equation -- Brusselator (native 2 species)",
        "form": "dMI/dt ~ c_self*H_i^2 + c_cross*H_a*H_b + linear + const",
        "config": {"seeds": seeds, "T": T, "N_ens": N_ens, "dt": 0.02,
                   "smooth_w": smooth_w},
        "per_seed": seed_results,
        "aggregate": agg,
        "reading_guide": {
            "is_a_real_dipole_if": "dR2_cross_term > 0 (H_a*H_b adds signal) "
                                   "AND dR2_dipole_over_relax > 0 (beats 1-D MI "
                                   "relaxation) AND shuffle_null_R2 ~ 0",
            "deflationary_if": "relax_R2 ~= full R2 (just 1-D MI relaxation) OR "
                               "shuffle_null_R2 high (column-count artifact)",
            "scale_caveat": "raw coeffs are units-dependent (INFO-051); read the "
                            "standardized betas for the unit-robust shape",
        },
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"=== flow-dipole Brusselator ({'CANARY' if args.canary else 'FULL'}) ===")
    print(f"seeds={seeds}  runtime={result['runtime_s']}s")
    print(f"full-model R2      : {agg['full_R2']['mean']:.3f} "
          f"+/- {agg['full_R2']['std']:.3f}")
    print(f"1-D MI-relax R2    : {agg['relax_R2']['mean']:.3f} "
          f"+/- {agg['relax_R2']['std']:.3f}")
    print(f"self-only R2       : {agg['model_C_self_only_R2']['mean']:.3f}")
    print(f"dR2 dipole>relax   : {agg['dR2_dipole_over_relax']['mean']:+.3f}")
    print(f"dR2 cross term     : {agg['dR2_cross_term']['mean']:+.3f}")
    print(f"shuffle-null R2    : {agg['shuffle_null_R2_mean']['mean']:.3f}")
    print(f"opposition fraction: {agg['opposition_fraction']:.2f}")
    print("standardized betas (mean +/- std across seeds):")
    for k in FEAT_NAMES[:5]:
        print(f"   {k:>8} : {agg[k]['mean']:+.3f} +/- {agg[k]['std']:.3f}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
