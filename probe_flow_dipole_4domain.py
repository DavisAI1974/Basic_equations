"""
*** SUPERSEDED / UNTRUSTWORTHY (Greg, S19) -- DO NOT BUILD ON THIS ***
This uses the TOY ODE simulators in per_domain_kbk.py (Duffing / Lotka-Volterra /
Brusselator / Burridge-Knopoff). Greg's S19 call: "why are we using toy brusselator?
... that's an old file, consider it untrustworthy" + the S10 rule "no point
fine-tuning fake data." The apparent shared opposition signature here did NOT survive
contact with REAL data (gravity LIGO + weak CMS dimuon both collapsed: R2~=shuffle-
null, no opposition). The real-data replacements are probe_flow_dipole_{gravity,weak,
em,strong,chem}.py. Kept only for the audit trail (Result Discipline: misses cataloged).
*** end deprecation notice ***

FLOW DIPOLE EQUATION across all 4 simulated domains (native 2-channel systems).
Backlog STEP 1 (Greg, S18) -- generalizes probe_flow_dipole_brusselator.py.

Paper form (davisai.ai/dipole, Section 2.2, confirmed verbatim S19):
    dMI_total/dt ~ sum_i c_self,i * H_i^2 + sum_{i<j} c_cross,ij * H_i*H_j + linear
Algebraic ratio (Section 2.1):  C = H_self / H_cross ,  H_self = internal Shannon
entropy, H_cross = MI(A;B). Opposition signature = self-coeffs vs cross-coeffs
opposite sign (paper: cellular 57.1%, organ 43.3%).

Each domain's two native channels (X1, X2) come from per_domain_kbk.py
simulators (physics Duffing, biology Lotka-Volterra, chemistry Brusselator,
geology Burridge-Knopoff) -- the SAME systems behind INFO-023/025/040, so the
flow-dipole equation here is directly comparable to that per-domain work.

DECISIVE TEST (Greg): do the flows share the SAME equation FORM (substrate) with
only COEFFICIENTS differing (expression)? Measured as cosine between the per-domain
standardized-beta vectors. High pairwise cos = shared form; the coefficients differ
= per-domain expression.

DEFLATIONARY CONTROLS (frame must not grade itself):
  (B) 1-D MI-relaxation: dMI/dt ~ a*MI + b  -> if it matches the full model, the
      dipole is just relaxation, not coupling.
  (C) self-only (no cross term) -> does H_a*H_b carry signal?
  (S) feature time-shuffle null -> R^2 must collapse if the fit is real.
  Standardized betas reported (marginal entropies are units-dependent, INFO-051;
  MI is the only scale-invariant operator).
Result Discipline: >=3 seeds, inter-seed scatter reported.

Run:  python probe_flow_dipole_4domain.py            (full, seeds 11/22/33)
      python probe_flow_dipole_4domain.py --canary    (1 seed, short)
"""
import numpy as np
import json
import time
import argparse

from per_domain_kbk import DOMAINS, build_ensemble_operator_matrix

COL = dict(Ha=0, Hb=1, Ha2=2, Hb2=3, HaHb=4, MI=5)
FEAT_NAMES = ["H_a^2", "H_b^2", "H_a*H_b", "H_a", "H_b", "1"]


def moving_average(v, w):
    if w <= 1:
        return v.copy()
    return np.convolve(v, np.ones(w) / w, mode="same")


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2


def standardized_betas(cols, y):
    Xs = np.column_stack([(c - c.mean()) / (c.std() + 1e-12) for c in cols])
    ys = (y - y.mean()) / (y.std() + 1e-12)
    beta, r2 = ols(Xs, ys)
    return beta, r2


def analyze_seed(simulate, seed, T, N_ens, smooth_w):
    X1, X2 = simulate(seed, N_ens=N_ens, T=T, dt=0.02)
    M, t_idx = build_ensemble_operator_matrix(X1, X2)
    dt = 0.02
    t = t_idx * dt
    Ha, Hb, MI = M[:, COL["Ha"]], M[:, COL["Hb"]], M[:, COL["MI"]]
    MIs = moving_average(MI, smooth_w)
    dMI = np.gradient(MIs, t)
    edge = smooth_w
    sl = slice(edge, len(t) - edge)
    Ha, Hb, MI, dMI = Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha * Ha, Hb * Hb, Ha * Hb

    XA = np.column_stack([Ha2, Hb2, HaHb, Ha, Hb, np.ones_like(Ha)])
    betaA, r2A = ols(XA, dMI)
    sbeta, _ = standardized_betas([Ha2, Hb2, HaHb, Ha, Hb], dMI)

    XB = np.column_stack([MI, np.ones_like(MI)])
    _, r2B = ols(XB, dMI)
    XC = np.column_stack([Ha2, Hb2, Ha, Hb, np.ones_like(Ha)])
    _, r2C = ols(XC, dMI)

    rng = np.random.default_rng(seed + 999)
    r2_null = []
    for _ in range(20):
        p = rng.permutation(len(dMI))
        Xn = np.column_stack([Ha2[p], Hb2[p], HaHb[p], Ha[p], Hb[p],
                              np.ones_like(Ha)])
        _, r2n = ols(Xn, dMI)
        r2_null.append(r2n)
    r2_null = np.array(r2_null)

    # algebraic ratio C = H_self / H_cross  (H_self = mean marginal entropy)
    Hself = 0.5 * (Ha + Hb)
    C = float(np.mean(Hself / (MI + 1e-9)))

    c_self = np.array([betaA[0], betaA[1]])
    c_cross = betaA[2]
    opposition = bool(c_cross != 0 and
                      np.all(np.sign(c_self) == -np.sign(c_cross)))

    return {
        "seed": seed,
        "std_betas": {FEAT_NAMES[i]: float(sbeta[i]) for i in range(5)},
        "raw_betas": {FEAT_NAMES[i]: float(betaA[i]) for i in range(6)},
        "R2_full": float(r2A),
        "R2_relax": float(r2B),
        "R2_self_only": float(r2C),
        "dR2_over_relax": float(r2A - r2B),
        "dR2_cross": float(r2A - r2C),
        "shuffle_null_R2_mean": float(r2_null.mean()),
        "C_ratio": C,
        "opposition_quad": opposition,
        "MI_range": [float(MI.min()), float(MI.max())],
    }


def agg_domain(seed_res):
    out = {}
    for k in FEAT_NAMES[:5]:
        v = [s["std_betas"][k] for s in seed_res]
        out[k] = {"mean": float(np.mean(v)), "std": float(np.std(v))}
    for fld in ["R2_full", "R2_relax", "R2_self_only", "dR2_over_relax",
                "dR2_cross", "shuffle_null_R2_mean", "C_ratio"]:
        v = [s[fld] for s in seed_res]
        out[fld] = {"mean": float(np.mean(v)), "std": float(np.std(v))}
    out["opposition_fraction"] = float(np.mean([s["opposition_quad"]
                                                for s in seed_res]))
    out["std_beta_vec"] = [out[k]["mean"] for k in FEAT_NAMES[:5]]
    return out


def cos(a, b):
    a, b = np.array(a), np.array(b)
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    if args.canary:
        seeds, T, N_ens, smooth_w = [11], 15.0, 300, 5
        out_path = "probe_flow_dipole_4domain_canary.json"
    else:
        seeds, T, N_ens, smooth_w = [11, 22, 33], 30.0, 600, 5
        out_path = "probe_flow_dipole_4domain_results.json"

    per_domain = {}
    for name, sim in DOMAINS.items():
        sr = [analyze_seed(sim, s, T, N_ens, smooth_w) for s in seeds]
        per_domain[name] = {"per_seed": sr, "aggregate": agg_domain(sr)}

    # decisive test: pairwise cosine between standardized-beta vectors
    names = list(per_domain.keys())
    pairwise = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a = per_domain[names[i]]["aggregate"]["std_beta_vec"]
            b = per_domain[names[j]]["aggregate"]["std_beta_vec"]
            pairwise[f"{names[i]} vs {names[j]}"] = round(cos(a, b), 3)
    cosvals = list(pairwise.values())

    result = {
        "probe": "flow-dipole equation across 4 simulated domains",
        "form": "dMI/dt ~ c_self*H_i^2 + c_cross*H_a*H_b + linear + const "
                "(davisai.ai/dipole Sec 2.2, verbatim)",
        "config": {"seeds": seeds, "T": T, "N_ens": N_ens, "dt": 0.02,
                   "smooth_w": smooth_w},
        "per_domain": per_domain,
        "decisive_test_pairwise_cos_of_std_betas": pairwise,
        "decisive_test_summary": {
            "mean_pairwise_cos": round(float(np.mean(cosvals)), 3),
            "min_pairwise_cos": round(float(np.min(cosvals)), 3),
            "max_pairwise_cos": round(float(np.max(cosvals)), 3),
            "reading": "high pairwise cos => shared FORM (substrate); the differing "
                       "coefficients => per-domain EXPRESSION. low cos => not a "
                       "shared equation, 'flow' stays a description not a law.",
        },
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"=== flow-dipole 4-domain ({'CANARY' if args.canary else 'FULL'}) "
          f"seeds={seeds} {result['runtime_s']}s ===")
    hdr = f"{'domain':<18}{'R2full':>8}{'R2relax':>8}{'dR2crs':>8}{'shuf':>7}{'opp':>6}{'C':>7}"
    print(hdr)
    for name in names:
        a = per_domain[name]["aggregate"]
        print(f"{name:<18}{a['R2_full']['mean']:>8.3f}{a['R2_relax']['mean']:>8.3f}"
              f"{a['dR2_cross']['mean']:>+8.3f}{a['shuffle_null_R2_mean']['mean']:>7.3f}"
              f"{a['opposition_fraction']:>6.2f}{a['C_ratio']['mean']:>7.2f}")
    print("\nstandardized betas (mean) per domain:")
    print(f"{'domain':<18}" + "".join(f"{k:>10}" for k in FEAT_NAMES[:5]))
    for name in names:
        a = per_domain[name]["aggregate"]
        print(f"{name:<18}" + "".join(f"{a[k]['mean']:>+10.2f}"
                                      for k in FEAT_NAMES[:5]))
    print("\nDECISIVE TEST -- pairwise cos of std-beta vectors:")
    for k, v in pairwise.items():
        print(f"   {k:<42}: {v:+.3f}")
    print(f"   mean {result['decisive_test_summary']['mean_pairwise_cos']:+.3f} "
          f"min {result['decisive_test_summary']['min_pairwise_cos']:+.3f} "
          f"max {result['decisive_test_summary']['max_pairwise_cos']:+.3f}")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
