"""
Falsifier suite for Session 4 INFO-018 (linear-drift N_REAL sweep).

Uses the EXACT protocol from extract_v1_linear_drift.py:
  - Vasicek 1D differential entropy (scipy.stats.differential_entropy)
  - 6D operator basis [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI(a,b)]
  - MI from 30x30 2D histogram, joint entropy
  - Windowed: w=40s, dt=0.1, 50% overlap -> 4 windows per realization
  - Per realization: 4 rows in operator matrix (one per window)
  - Pool across N_REAL realizations -> 4*N_REAL x 6 matrix
  - Center columns, SVD, smallest right singular vector = null direction
  - Project to (H_a^2, H_b^2, H_a*H_b) subspace = indices [2,3,4]

Cases run (3 seeds each, N_REAL=500):

  BASELINE     v_a=v_b=1, sigma=0.1, rho=0      (Session 4 spec)
  FALSIFIER 1  v_a=v_b=0, sigma=0.1, rho=0      (no drift)
  FALSIFIER 2  v_a=1, v_b=2, sigma=0.1, rho=0   (asymmetric drift)
  FALSIFIER 3  v_a=v_b=1, sigma=0.1, rho in {0.3, 0.7, 0.95}

Predictions if INFO-018 is extractor artifact:
  - SV4,5,6 nearly equal -> non-unique null direction
  - FALSIFIER 1 direction differs from BASELINE (different distribution shape)
    but inter-seed cosines remain high -> the "convergence" is geometric, not
    physical
  - FALSIFIER 3 direction does NOT track rho (because the protocol's MI
    estimator IS sensitive to coupling, but the (H_a^2, H_b^2, H_a*H_b)
    subspace projection mostly washes that out)

Predictions if INFO-018 is real:
  - SV4 << SV5 (unique null direction)
  - FALSIFIER 3 direction tracks rho monotonically
"""

import numpy as np
from scipy.stats import differential_entropy
import json
import time

# ---------------------------------------------------------------------------
# operators -- matches extract_v1_linear_drift.py exactly
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

# ---------------------------------------------------------------------------
# system: linear drift with optional asymmetric drift and coupling
# ---------------------------------------------------------------------------
def linear_drift(T, dt, v_a, v_b, seed, noise_std=0.1, rho=0.0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    eps_a = rng.normal(0.0, noise_std, n)
    eps_b_indep = rng.normal(0.0, noise_std, n)
    eps_b = rho * eps_a + np.sqrt(max(0.0, 1.0 - rho**2)) * eps_b_indep
    a = v_a * t + eps_a
    b = v_b * t + eps_b
    return a, b

# ---------------------------------------------------------------------------
# windowed-null extraction
# ---------------------------------------------------------------------------
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

def extract_v1(N_REAL, master_seed, v_a, v_b, sigma=0.1, rho=0.0,
               T=100, dt=0.1, window=40):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = linear_drift(T, dt, v_a, v_b, seed, noise_std=sigma, rho=rho)
        ops = windowed_operators(a, b, window, dt)
        all_ops.append(ops)
    M = np.vstack(all_ops)
    M_c = M - M.mean(axis=0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return dict(null_full=Vt[-1], singular_values=S, M_shape=M.shape,
                M_mean=M.mean(axis=0))

def project_3d(v6):
    sub = v6[[2, 3, 4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    # canonical sign: positive H_a*H_b component
    if sub[2] < 0:
        sub = -sub
    return sub

def cos_sim(u, v):
    return float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

# ---------------------------------------------------------------------------
# reference directions
# ---------------------------------------------------------------------------
DIPOLE_PLUS  = np.array([1.0, 1.0, +2.0]) / np.sqrt(6.0)   # +(H_a+H_b)^2
DIPOLE_MINUS = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)   # -(H_a-H_b)^2, OU attractor
COMMON_MODE  = np.array([1.0, 1.0, +1.0]) / np.sqrt(3.0)

# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def run_case(name, v_a, v_b, sigma, rho, N_REAL=500, seeds=(0, 1, 2),
             verbose=True):
    if verbose:
        print()
        print("=" * 72)
        print(f"  {name}")
        print("=" * 72)
    t0 = time.time()
    per_seed = []
    for s in seeds:
        r = extract_v1(N_REAL=N_REAL, master_seed=s,
                       v_a=v_a, v_b=v_b, sigma=sigma, rho=rho)
        v3 = project_3d(r['null_full'])
        per_seed.append(dict(seed=s,
                             null_full=r['null_full'],
                             v3=v3,
                             singular_values=r['singular_values'],
                             M_mean=r['M_mean']))

    sv_stack = np.stack([p['singular_values'] for p in per_seed])
    sv_mean = sv_stack.mean(axis=0)
    v3_stack = np.stack([p['v3'] for p in per_seed])
    mean_v3 = v3_stack.mean(axis=0)
    mean_v3 = mean_v3 / np.linalg.norm(mean_v3)
    M_mean = np.mean([p['M_mean'] for p in per_seed], axis=0)

    inter_pairs = []
    for i in range(len(per_seed)):
        for j in range(i+1, len(per_seed)):
            inter_pairs.append(cos_sim(per_seed[i]['v3'], per_seed[j]['v3']))
    inter_mean = float(np.mean(inter_pairs))

    if verbose:
        print(f"  Per-realization op-vector means (over windows & realizations):")
        print(f"    <H_a>={M_mean[0]:+.4f}  <H_b>={M_mean[1]:+.4f}  "
              f"<H_a^2>={M_mean[2]:+.4f}  <H_b^2>={M_mean[3]:+.4f}  "
              f"<H_a*H_b>={M_mean[4]:+.4f}  <MI>={M_mean[5]:+.4f}")
        print(f"  Singular values (mean across seeds):")
        for i, s in enumerate(sv_mean):
            print(f"    SV_{i+1} = {s:.4e}")
        print(f"  Gap ratios SV_i / SV_(i+1):")
        for i in range(len(sv_mean)-1):
            print(f"    SV_{i+1}/SV_{i+2} = {sv_mean[i]/sv_mean[i+1]:8.2f}")
        print(f"  Per-seed v3 (after canonical sign on H_a*H_b coord):")
        for p in per_seed:
            print(f"    seed={p['seed']}: ({p['v3'][0]:+.4f}, {p['v3'][1]:+.4f}, {p['v3'][2]:+.4f})")
        print(f"  Mean v3:        ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
        print(f"  cos to (+1,+1,+2)/sqrt6 = {cos_sim(mean_v3, DIPOLE_PLUS):+.4f}")
        print(f"  cos to (+1,+1,-2)/sqrt6 = {cos_sim(mean_v3, DIPOLE_MINUS):+.4f}")
        print(f"  cos to (+1,+1,+1)/sqrt3 = {cos_sim(mean_v3, COMMON_MODE):+.4f}")
        print(f"  Inter-seed cos (3D): {[f'{c:+.4f}' for c in inter_pairs]}")
        print(f"  Mean inter-seed cos: {inter_mean:+.4f}")
        print(f"  Elapsed: {time.time()-t0:.1f}s")

    return dict(
        name=name, params=dict(v_a=v_a, v_b=v_b, sigma=sigma, rho=rho, N_REAL=N_REAL),
        singular_values=sv_mean.tolist(),
        per_seed_sv=[p['singular_values'].tolist() for p in per_seed],
        mean_v3=mean_v3.tolist(),
        per_seed_v3=[p['v3'].tolist() for p in per_seed],
        per_seed_null_full=[p['null_full'].tolist() for p in per_seed],
        op_means=M_mean.tolist(),
        cos_dipole_plus=cos_sim(mean_v3, DIPOLE_PLUS),
        cos_dipole_minus=cos_sim(mean_v3, DIPOLE_MINUS),
        cos_common=cos_sim(mean_v3, COMMON_MODE),
        inter_seed_cos_pairs=inter_pairs,
        inter_seed_cos_mean=inter_mean,
        elapsed_s=time.time()-t0,
    )

def main():
    t_start = time.time()
    cases = [
        ("BASELINE      v_a=v_b=1, sigma=0.1, rho=0    (Session 4 spec)",
            dict(v_a=1.0, v_b=1.0, sigma=0.1, rho=0.0)),
        ("FALSIFIER 1   v_a=v_b=0, sigma=0.1, rho=0    (no drift)",
            dict(v_a=0.0, v_b=0.0, sigma=0.1, rho=0.0)),
        ("FALSIFIER 2   v_a=1, v_b=2, sigma=0.1, rho=0 (asymmetric drift)",
            dict(v_a=1.0, v_b=2.0, sigma=0.1, rho=0.0)),
        ("FALSIFIER 3a  v_a=v_b=1, sigma=0.1, rho=0.30 (mild coupling)",
            dict(v_a=1.0, v_b=1.0, sigma=0.1, rho=0.30)),
        ("FALSIFIER 3b  v_a=v_b=1, sigma=0.1, rho=0.70 (strong coupling)",
            dict(v_a=1.0, v_b=1.0, sigma=0.1, rho=0.70)),
        ("FALSIFIER 3c  v_a=v_b=1, sigma=0.1, rho=0.95 (near-shared noise)",
            dict(v_a=1.0, v_b=1.0, sigma=0.1, rho=0.95)),
    ]

    summary = []
    for name, params in cases:
        rec = run_case(name, **params)
        summary.append(rec)

    out_path = '/home/user/Basic_equations/falsifier_suite_INFO018_results.json'
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print()
    print("=" * 78)
    print("  VERDICT TABLE")
    print("=" * 78)
    print(f"  {'case':<28s}{'mean v3':>26s}{'cos+2':>8s}{'cos-2':>8s}{'SV4/5':>8s}")
    for rec in summary:
        v = rec['mean_v3']
        short = rec['name'][:28]
        sv45 = rec['singular_values'][3] / rec['singular_values'][4]
        print(f"  {short:<28s}"
              f"({v[0]:+.3f},{v[1]:+.3f},{v[2]:+.3f})"
              f"{rec['cos_dipole_plus']:>+8.3f}"
              f"{rec['cos_dipole_minus']:>+8.3f}"
              f"{sv45:>8.2f}")
    print()
    print(f"  Total runtime: {time.time()-t_start:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
