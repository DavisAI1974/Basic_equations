"""
Probe: how does the converged direction look with 30 seeds vs 3?

Runs the extract_v1 protocol with 30 master_seeds at N_REAL=500 on two cases:
  BASELINE:    v_a=v_b=1, sigma=0.1, rho=0   (Session 4 spec)
  GAUSS_V0:    v_a=v_b=0, sigma=0.1, rho=0   (zero-drift Gaussian)

Reports per case:
  - per-seed v3 (canonical sign on H_a*H_b)
  - mean v3 across 30 seeds
  - per-coordinate standard deviation across seeds
  - distribution of inter-seed cosines
  - cos of mean v3 to (+1,+1,+2)/sqrt(6) and (-1,-1,+2)/sqrt(6)

This is a different-signal generator. No verdict.
"""

import numpy as np
from scipy.stats import differential_entropy
import json
import time

# ---------------------------------------------------------------------------
# operators
# ---------------------------------------------------------------------------
def entropy_1d(x):
    return differential_entropy(x, method='vasicek')

def entropy_2d_kde(xy, n_grid=30):
    H, xedges, yedges = np.histogram2d(xy[:, 0], xy[:, 1], bins=n_grid, density=True)
    dx = xedges[1] - xedges[0]
    dy = yedges[1] - yedges[0]
    p = H.flatten()
    p = p[p > 0]
    return -np.sum(p * np.log(p)) * dx * dy

def mutual_info(a, b):
    H_a = entropy_1d(a)
    H_b = entropy_1d(b)
    H_ab = entropy_2d_kde(np.column_stack([a, b]))
    return H_a + H_b - H_ab

def operators(a, b):
    H_a = entropy_1d(a)
    H_b = entropy_1d(b)
    MI = mutual_info(a, b)
    return np.array([H_a, H_b, H_a**2, H_b**2, H_a*H_b, MI])

def linear_drift(T, dt, v, seed, sigma=0.1):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    a = v*t + rng.normal(0, sigma, n)
    b = v*t + rng.normal(0, sigma, n)
    return a, b

def windowed_operators(a, b, window, dt, overlap=0.5):
    win_samples = int(window / dt)
    step = int(win_samples * (1 - overlap))
    n = len(a)
    ops_list = []
    start = 0
    while start + win_samples <= n:
        ops_list.append(operators(a[start:start+win_samples],
                                  b[start:start+win_samples]))
        start += step
    return np.array(ops_list)

def extract_v1(v, N_REAL, master_seed, T=100, dt=0.1, window=40, sigma=0.1):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = linear_drift(T, dt, v, seed, sigma=sigma)
        ops = windowed_operators(a, b, window, dt)
        all_ops.append(ops)
    M = np.vstack(all_ops)
    M_c = M - M.mean(axis=0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return dict(null_full=Vt[-1], singular_values=S)

def project_3d(v6):
    sub = v6[[2, 3, 4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    if sub[2] < 0:
        sub = -sub
    return sub

def cos_sim(u, v):
    return float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

REF_PLUS  = np.array([+1.0, +1.0, +2.0]) / np.sqrt(6.0)
REF_MINUS = np.array([-1.0, -1.0, +2.0]) / np.sqrt(6.0)

def run_case(label, v, N_seeds=30, N_REAL=500):
    print()
    print("=" * 72)
    print(f"  CASE: {label}  (v={v}, N_REAL={N_REAL}, N_seeds={N_seeds})")
    print("=" * 72)
    t0 = time.time()
    v3_list = []
    for s in range(N_seeds):
        r = extract_v1(v=v, N_REAL=N_REAL, master_seed=s)
        v3 = project_3d(r['null_full'])
        v3_list.append(v3)
        print(f"  seed={s:2d}: v3=({v3[0]:+.4f}, {v3[1]:+.4f}, {v3[2]:+.4f})")
    v3_arr = np.array(v3_list)
    mean_v3 = v3_arr.mean(axis=0)
    mean_v3_norm = mean_v3 / np.linalg.norm(mean_v3)
    sd_per_coord = v3_arr.std(axis=0)

    inter_pairs = []
    for i in range(N_seeds):
        for j in range(i+1, N_seeds):
            inter_pairs.append(cos_sim(v3_list[i], v3_list[j]))
    inter_pairs = np.array(inter_pairs)

    print(f"  mean v3 (raw):                 ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
    print(f"  mean v3 (re-normalized):       ({mean_v3_norm[0]:+.4f}, {mean_v3_norm[1]:+.4f}, {mean_v3_norm[2]:+.4f})")
    print(f"  per-coord SD across seeds:     ({sd_per_coord[0]:+.4f}, {sd_per_coord[1]:+.4f}, {sd_per_coord[2]:+.4f})")
    print(f"  cos to (+1,+1,+2)/sqrt6:       {cos_sim(mean_v3_norm, REF_PLUS):+.4f}")
    print(f"  cos to (-1,-1,+2)/sqrt6:       {cos_sim(mean_v3_norm, REF_MINUS):+.4f}")
    print(f"  inter-seed cos: mean={inter_pairs.mean():+.4f}, "
          f"min={inter_pairs.min():+.4f}, max={inter_pairs.max():+.4f}, "
          f"sd={inter_pairs.std():+.4f}")
    print(f"  elapsed: {time.time()-t0:.1f}s")

    return dict(
        label=label, v=v, N_seeds=N_seeds, N_REAL=N_REAL,
        per_seed_v3=[v.tolist() for v in v3_list],
        mean_v3_raw=mean_v3.tolist(),
        mean_v3_normalized=mean_v3_norm.tolist(),
        per_coord_sd=sd_per_coord.tolist(),
        cos_plus_plus_plus=cos_sim(mean_v3_norm, REF_PLUS),
        cos_minus_minus_plus=cos_sim(mean_v3_norm, REF_MINUS),
        inter_seed_cos_mean=float(inter_pairs.mean()),
        inter_seed_cos_min=float(inter_pairs.min()),
        inter_seed_cos_max=float(inter_pairs.max()),
        inter_seed_cos_sd=float(inter_pairs.std()),
    )

def main():
    t_start = time.time()
    out = []
    out.append(run_case("BASELINE_v=1", v=1.0))
    out.append(run_case("GAUSS_V0",     v=0.0))

    out_path = '/home/user/Basic_equations/probe_highseed_convergence_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print()
    print(f"  Total runtime: {time.time()-t_start:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
