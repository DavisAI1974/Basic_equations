"""
Probe: where do diverse systems land in the protocol's 3D subspace?

Runs the extract_v1 protocol (from extract_v1_linear_drift.py) on eight
two-stream systems, three seeds each, N_REAL=500. Reports per-system:
  - mean 3D direction (canonical sign on H_a*H_b coord)
  - full 6D null vector per seed
  - singular value spectrum
  - operator means
  - cosines against two reference directions from prior sessions:
       (-1,-1,+2)/sqrt(6)  (INFO-014 attractor)
       (+1,+1,+2)/sqrt(6)  (Session 4 linear-drift converged direction)
    Reported as reference distances, NOT as classification labels.

This is a different-signal generator. It may speak to INFO-018, INFO-014,
INFO-017, or somewhere else entirely. No verdict in script or output.

Systems:
  OU           dx = -gamma*x dt + sigma dW   (gamma=1, sigma=1)
  uniform      x ~ Uniform[-1, +1] iid per step
  gauss_small  x ~ N(0, 0.1^2) iid per step  (matches Falsifier 1)
  gauss_large  x ~ N(0, 1.0^2) iid per step  (scale control)
  student_t    x ~ t(df=3) iid per step
  laplace      x ~ Laplace(0, 1) iid per step
  sine         x_a = sin(omega t) + noise, x_b = sin(omega t + phi) + noise
  logistic     x_{n+1} = r x_n (1 - x_n), r=3.9, two independent orbits
"""

import numpy as np
from scipy.stats import differential_entropy
import json
import time

# ---------------------------------------------------------------------------
# operators (exact match to extract_v1_linear_drift.py)
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
# system generators
# ---------------------------------------------------------------------------
def ou_pair(T, dt, seed, gamma=1.0, sigma=1.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = np.zeros(n); b = np.zeros(n)
    for k in range(1, n):
        a[k] = a[k-1] - gamma*a[k-1]*dt + sigma*np.sqrt(dt)*rng.normal()
        b[k] = b[k-1] - gamma*b[k-1]*dt + sigma*np.sqrt(dt)*rng.normal()
    return a, b

def uniform_pair(T, dt, seed):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.uniform(-1.0, 1.0, n), rng.uniform(-1.0, 1.0, n)

def gauss_pair(T, dt, seed, sigma):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.normal(0.0, sigma, n), rng.normal(0.0, sigma, n)

def studentt_pair(T, dt, seed, df=3):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.standard_t(df, n), rng.standard_t(df, n)

def laplace_pair(T, dt, seed, scale=1.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.laplace(0.0, scale, n), rng.laplace(0.0, scale, n)

def sine_pair(T, dt, seed, omega=0.5, phi=0.7, noise=0.05):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    t = np.arange(n) * dt
    a = np.sin(omega*t) + noise*rng.normal(size=n)
    b = np.sin(omega*t + phi) + noise*rng.normal(size=n)
    return a, b

def logistic_pair(T, dt, seed, r=3.9):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = np.zeros(n); b = np.zeros(n)
    a[0] = rng.uniform(0.1, 0.9)
    b[0] = rng.uniform(0.1, 0.9)
    for k in range(1, n):
        a[k] = r*a[k-1]*(1 - a[k-1])
        b[k] = r*b[k-1]*(1 - b[k-1])
    return a, b

SYSTEMS = {
    'OU_gamma1_sigma1':  lambda T, dt, seed: ou_pair(T, dt, seed, 1.0, 1.0),
    'uniform_pm1':       lambda T, dt, seed: uniform_pair(T, dt, seed),
    'gauss_small_0p1':   lambda T, dt, seed: gauss_pair(T, dt, seed, 0.1),
    'gauss_large_1p0':   lambda T, dt, seed: gauss_pair(T, dt, seed, 1.0),
    'student_t_df3':     lambda T, dt, seed: studentt_pair(T, dt, seed, 3),
    'laplace_scale1':    lambda T, dt, seed: laplace_pair(T, dt, seed, 1.0),
    'sine_w0p5_phi0p7':  lambda T, dt, seed: sine_pair(T, dt, seed, 0.5, 0.7, 0.05),
    'logistic_r3p9':     lambda T, dt, seed: logistic_pair(T, dt, seed, 3.9),
}

# ---------------------------------------------------------------------------
# extraction (matches extract_v1)
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

def extract_v1(system_fn, N_REAL, master_seed, T=100, dt=0.1, window=40):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = system_fn(T, dt, seed)
        ops = windowed_operators(a, b, window, dt)
        all_ops.append(ops)
    M = np.vstack(all_ops)
    M_c = M - M.mean(axis=0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return dict(null_full=Vt[-1], singular_values=S, M_mean=M.mean(axis=0))

def project_3d(v6):
    sub = v6[[2, 3, 4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    if sub[2] < 0:
        sub = -sub
    return sub

def cos_sim(u, v):
    return float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

REF_MINUS_MINUS_PLUS = np.array([-1.0, -1.0, +2.0]) / np.sqrt(6.0)  # INFO-014
REF_PLUS_PLUS_PLUS   = np.array([+1.0, +1.0, +2.0]) / np.sqrt(6.0)  # Session 4

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    t_start = time.time()
    seeds = (0, 1, 2)
    N_REAL = 500
    out = {}

    for name, fn in SYSTEMS.items():
        print()
        print("=" * 72)
        print(f"  SYSTEM: {name}")
        print("=" * 72)
        t0 = time.time()
        per_seed = []
        for s in seeds:
            r = extract_v1(fn, N_REAL=N_REAL, master_seed=s)
            v3 = project_3d(r['null_full'])
            per_seed.append(dict(seed=s,
                                 null_full=r['null_full'].tolist(),
                                 v3=v3.tolist(),
                                 singular_values=r['singular_values'].tolist(),
                                 M_mean=r['M_mean'].tolist()))
            print(f"  seed={s}: v3=({v3[0]:+.4f}, {v3[1]:+.4f}, {v3[2]:+.4f})")
            print(f"          null_full={r['null_full'].round(4).tolist()}")
            print(f"          SVs={[f'{x:.3e}' for x in r['singular_values']]}")
            print(f"          op_means={r['M_mean'].round(4).tolist()}")

        v3_stack = np.array([p['v3'] for p in per_seed])
        mean_v3 = v3_stack.mean(axis=0)
        mean_v3 = mean_v3 / np.linalg.norm(mean_v3)
        pair_cos = []
        for i in range(len(per_seed)):
            for j in range(i+1, len(per_seed)):
                pair_cos.append(cos_sim(per_seed[i]['v3'], per_seed[j]['v3']))

        cos_mmp = cos_sim(mean_v3, REF_MINUS_MINUS_PLUS)
        cos_ppp = cos_sim(mean_v3, REF_PLUS_PLUS_PLUS)
        print(f"  mean_v3:                ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
        print(f"  cos to (-1,-1,+2)/sqrt6: {cos_mmp:+.4f}    (INFO-014 ref)")
        print(f"  cos to (+1,+1,+2)/sqrt6: {cos_ppp:+.4f}    (Session 4 ref)")
        print(f"  inter-seed cos pairs:    {[f'{c:+.4f}' for c in pair_cos]}")
        print(f"  mean inter-seed cos:     {np.mean(pair_cos):+.4f}")
        print(f"  elapsed: {time.time()-t0:.1f}s")

        out[name] = dict(
            per_seed=per_seed,
            mean_v3=mean_v3.tolist(),
            cos_to_minus_minus_plus=cos_mmp,
            cos_to_plus_plus_plus=cos_ppp,
            inter_seed_cos_pairs=pair_cos,
            inter_seed_cos_mean=float(np.mean(pair_cos)),
        )

    out_path = '/home/user/Basic_equations/probe_system_directions_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print()
    print(f"  Total runtime: {time.time()-t_start:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
