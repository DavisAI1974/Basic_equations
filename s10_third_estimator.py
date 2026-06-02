"""
S10 item (4): third-estimator triangulation on the INFO-033 flip.

Follow-up A (s9_std_vs_raw.py) showed the MI-dominant substrate flip is
PARTLY procedure-inflated (RAW |MI|~0.99), PARTLY structural (per-column
STANDARDIZED |MI|~0.65, still past the 0.5 flip line). Both extractions
used the SAME histogram MI estimator, so "procedure" there meant the
SVD scaling (raw-centered vs per-column-standardized), not the MI
estimator. Follow-up A flagged the open extension: a third, genuinely
different MI ESTIMATOR (kNN/KSG) to separate estimator-family artifact
from structure.

This script does exactly that. It rebuilds the SAME operator matrix as
the original flip extraction, with the entropy columns [H_a, H_b, H_a^2,
H_b^2, H_a*H_b] held byte-identical to the original (histogram entropy),
swapping ONLY the MI column between:
    hist : mi_hist_2d        (original; what INFO-033 / follow-up A used)
    ksg  : KSG estimator 1   (Kraskov-Stoegbauer-Grassberger 2004, kNN)
For each MI estimator we run BOTH extractions:
    RAW : extract_v1 (centered SVD; original constancy detector)
    STD : per-column standardized SVD (correlation detector; follow-up A)
and report |MI coef| in the smallest-null direction (the flip metric).

Reading rule (speaking posture, Rule C, BEFORE running): I think the
structural core (STD |MI|~0.65) MAY survive a kNN estimator if it is a
real linear dependency, and the spectacular RAW ~0.99 MAY shrink because
kNN MI does not collapse to a near-constant the way histogram MI does
under asymmetry. But we wait on the data and where it points. If the
kNN STD |MI| also lands ~0.5-0.7, the structural core is estimator-
robust (two independent estimators agree). If it drops to ~0, the ~0.65
was histogram-specific. No verdict in advance.

Same sims, same configs, same seeds as follow-up A so numbers line up.
Canary: 1 seed, N_ens=200, T=10. Full: 3 seeds, N_ens=600, T=30.
"""
import sys
import json
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.special import digamma

# entropy_hist_1d lives in per_domain_kbk, mi_hist_2d in kbk_pipeline
from per_domain_kbk import entropy_hist_1d
from kbk_pipeline import mi_hist_2d, extract_v1
from em_strong_glance import simulate_em_config, simulate_strong_config
from gravity_glance import simulate_gravity_config
from four_force_probe import simulate_weak

ATTR3 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22, 33)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)


def mi_ksg(x, y, k=4, rng=None):
    """Kraskov-Stoegbauer-Grassberger estimator 1 (nats), max-norm.

    MI = psi(k) + psi(N) - <psi(n_x+1) + psi(n_y+1)>
    n_x / n_y count points strictly within the k-th-neighbor joint
    distance along each marginal. Tiny jitter breaks ties (standard).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < k + 2:
        return 0.0
    if rng is None:
        rng = np.random.default_rng(0)
    # jitter at ~1e-10 of the data scale to break exact ties
    x = x + 1e-10 * np.std(x) * rng.standard_normal(n)
    y = y + 1e-10 * np.std(y) * rng.standard_normal(n)
    xy = np.column_stack([x, y])
    tree = cKDTree(xy)
    # Chebyshev (max-norm) distance to the k-th neighbor (exclude self).
    d, _ = tree.query(xy, k=k + 1, p=np.inf)
    eps = d[:, -1]  # distance to k-th neighbor
    tx = cKDTree(x[:, None])
    ty = cKDTree(y[:, None])
    nx = np.array(tx.query_ball_point(x[:, None], eps - 1e-15, p=np.inf,
                                      return_length=True)) - 1
    ny = np.array(ty.query_ball_point(y[:, None], eps - 1e-15, p=np.inf,
                                      return_length=True)) - 1
    nx = np.clip(nx, 0, None)
    ny = np.clip(ny, 0, None)
    mi = (digamma(k) + digamma(n)
          - np.mean(digamma(nx + 1) + digamma(ny + 1)))
    return float(max(0.0, mi))


def build_matrix(X1, X2, mi_mode="hist", bins_H=30, bins_MI=18,
                 burn_in_frac=0.2, ksg_k=4, seed=0):
    """Per time-step operator matrix. Entropy cols are histogram (as in
    the original build_ensemble_operator_matrix). Only the MI column
    swaps between 'hist' and 'ksg'."""
    rng = np.random.default_rng(seed + 9973)
    T_steps = X1.shape[0]
    t_start = int(burn_in_frac * T_steps)
    rows = []
    for t in range(t_start, T_steps):
        Ha = entropy_hist_1d(X1[t], bins=bins_H)
        Hb = entropy_hist_1d(X2[t], bins=bins_H)
        if mi_mode == "hist":
            MI = mi_hist_2d(X1[t], X2[t], bins=bins_MI)
        elif mi_mode == "ksg":
            MI = mi_ksg(X1[t], X2[t], k=ksg_k, rng=rng)
        else:
            raise ValueError(mi_mode)
        rows.append([Ha, Hb, Ha * Ha, Hb * Hb, Ha * Hb, MI])
    return np.array(rows)


def extract_std(M):
    mu = M.mean(0)
    sd = M.std(0) + 1e-12
    Z = (M - mu) / sd
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    v = Vt[-1]
    j = int(np.argmax(np.abs(v)))
    if v[j] < 0:
        v = -v
    return v, S


def metrics(M):
    vr, _, _, _ = extract_v1(M)        # raw
    vs, _ = extract_std(M)             # standardized
    subr = vr[2:5] / (np.linalg.norm(vr[2:5]) + 1e-30)
    subs = vs[2:5] / (np.linalg.norm(vs[2:5]) + 1e-30)
    return (abs(vr[5]), abs(vs[5]),
            float(abs(subr @ ATTR3)), float(abs(subs @ ATTR3)),
            float(M[:, 5].std()), float(M[:, 5].mean()))


def run(name, simfn, **skw):
    acc = {m: [] for m in
           ("h_raw", "h_std", "k_raw", "k_std",
            "h_mi_std", "k_mi_std", "h_mi_mean", "k_mi_mean")}
    for s in seeds:
        X1, X2 = simfn(seed=s, **cfg, **skw)
        Mh = build_matrix(X1, X2, mi_mode="hist", seed=s)
        Mk = build_matrix(X1, X2, mi_mode="ksg", seed=s)
        hr, hs, _, _, hstd, hmean = metrics(Mh)
        kr, ks, _, _, kstd, kmean = metrics(Mk)
        acc["h_raw"].append(hr); acc["h_std"].append(hs)
        acc["k_raw"].append(kr); acc["k_std"].append(ks)
        acc["h_mi_std"].append(hstd); acc["k_mi_std"].append(kstd)
        acc["h_mi_mean"].append(hmean); acc["k_mi_mean"].append(kmean)
    out = {k: (float(np.mean(v)), float(np.std(v))) for k, v in acc.items()}
    print(f"{name:20s}  "
          f"HIST raw|MI|={out['h_raw'][0]:.3f} std|MI|={out['h_std'][0]:.3f}   "
          f"KSG raw|MI|={out['k_raw'][0]:.3f}+/-{out['k_raw'][1]:.2f} "
          f"std|MI|={out['k_std'][0]:.3f}+/-{out['k_std'][1]:.2f}   "
          f"(MI sd hist={out['h_mi_std'][0]:.2e} ksg={out['k_mi_std'][0]:.2e})")
    return dict(name=name, **{k: list(v) for k, v in out.items()})


t0 = time.time()
R = {}
print("=" * 118)
print(f"S10 THIRD-ESTIMATOR (item 4) {'[CANARY]' if CANARY else '[FULL]'} "
      f"seeds={seeds} cfg={cfg}")
print("HIST = original histogram MI; KSG = Kraskov kNN MI. raw=centered SVD, "
      "std=per-column standardized SVD.")
print("Flip metric = |MI coef| in smallest null (>0.5 => MI-dominant substrate).")
print("=" * 118)

print("\n--- gravity sweep (universal coupling) ---")
for o2 in [0.8, 0.9, 1.0, 1.2]:
    R[f"grav_0.8/{o2}"] = run(f"gravity 0.8/{o2}", simulate_gravity_config,
                              omega1=0.8, omega2=o2, universal=True)
print("--- EM sweep ---")
for o2 in [1.0, 1.3, 1.5, 2.0]:
    R[f"em_1.0/{o2}"] = run(f"EM 1.0/{o2}", simulate_em_config,
                            omega1=1.0, omega2=o2)
print("--- strong sweep ---")
for o2 in [0.5, 1.0, 1.5, 2.0]:
    R[f"strong_0.5/{o2}"] = run(f"strong 0.5/{o2}", simulate_strong_config,
                                omega1=0.5, omega2=o2)
print("--- weak baseline ---")
R["weak_sym"] = run("weak 1.0/1.0", simulate_weak)

out = dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg,
                     ksg_k=4, elapsed_s=round(time.time() - t0, 1)),
           legend="tuple = (mean, std) across seeds; keys h_=hist k_=ksg, "
                  "_raw/_std = extraction, _mi_std/_mi_mean = MI column stats",
           results=R)
fn = "s10_third_estimator_canary.json" if CANARY else "s10_third_estimator_results.json"
json.dump(out, open(fn, "w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
