"""
Probe: what does extract_v1 return on inputs with no temporal / dynamical structure?

Different-signal generator (Rule B). Not a falsifier. Result reads as:
"on these inputs, the protocol returned these directions."

Protocol matches extract_v1_linear_drift.py exactly:
  basis = [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]
  T=100, dt=0.1, window=40, overlap=0.5
  N_REAL realizations, pool, center, smallest right SV = null direction.
  Project to indices [2, 3, 4] for the 3D subspace.

Cases (all 7 use the same seed conventions as probe_system_directions.py):

  scrambled_baseline_joint  linear-drift baseline, same random permutation
                            on a and b per realization (preserves joint pair)
  scrambled_baseline_indep  linear-drift baseline, independent permutations
                            on a and b (also kills a-b pairing)
  scrambled_OU_joint        OU(gamma=1, sigma=1), same joint scramble
  near_constant             x_a = c_a + 1e-3 * N(0,1) per step, c_a drawn
                            once per realization from N(0,1); same for b
  identical_streams         x_a = x_b = N(0,1) iid per step (extreme coupling)
  mixed_marginals           x_a Gaussian, x_b uniform (different marginals)
  random_walk_indep         x_a, x_b = cumsum of independent N(0, 0.1)

Reporting per case: per-seed v3, mean v3, cos to (+,+,+), cos to (-,-,+),
inter-seed cos, full singular spectrum, operator means.

Run modes:
  python probe_null_inputs.py canary   # N_REAL=50, 1 seed, cases [0, 3]
  python probe_null_inputs.py          # full: N_REAL=500, 3 seeds, all cases
"""
import sys
import json
import time
import numpy as np
from scipy.stats import differential_entropy

# ---------- operators (exact match to extract_v1_linear_drift.py) ----------

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

# ---------- system generators ----------

def linear_drift_baseline(T, dt, seed, v=1.0, noise_std=0.1):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    a = v * t + rng.normal(0, noise_std, n)
    b = v * t + rng.normal(0, noise_std, n)
    return a, b

def ou_pair(T, dt, seed, gamma=1.0, sigma=1.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = np.zeros(n); b = np.zeros(n)
    for k in range(1, n):
        a[k] = a[k-1] - gamma*a[k-1]*dt + sigma*np.sqrt(dt)*rng.normal()
        b[k] = b[k-1] - gamma*b[k-1]*dt + sigma*np.sqrt(dt)*rng.normal()
    return a, b

def scrambled_baseline_joint(T, dt, seed):
    rng = np.random.default_rng(seed ^ 0xABCDEF)
    a, b = linear_drift_baseline(T, dt, seed)
    perm = rng.permutation(len(a))
    return a[perm], b[perm]

def scrambled_baseline_indep(T, dt, seed):
    rng = np.random.default_rng(seed ^ 0xABCDEF)
    a, b = linear_drift_baseline(T, dt, seed)
    pa = rng.permutation(len(a))
    pb = rng.permutation(len(b))
    return a[pa], b[pb]

def scrambled_OU_joint(T, dt, seed):
    rng = np.random.default_rng(seed ^ 0xABCDEF)
    a, b = ou_pair(T, dt, seed)
    perm = rng.permutation(len(a))
    return a[perm], b[perm]

def near_constant(T, dt, seed, sigma_floor=1e-3):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    c_a = rng.normal()
    c_b = rng.normal()
    a = c_a + sigma_floor * rng.normal(size=n)
    b = c_b + sigma_floor * rng.normal(size=n)
    return a, b

def identical_streams(T, dt, seed):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    x = rng.normal(0.0, 1.0, n)
    return x, x.copy()

def mixed_marginals(T, dt, seed):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = rng.normal(0.0, 1.0, n)
    b = rng.uniform(-1.0, 1.0, n)
    return a, b

def random_walk_indep(T, dt, seed, step_sigma=0.1):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = np.cumsum(rng.normal(0, step_sigma, n))
    b = np.cumsum(rng.normal(0, step_sigma, n))
    return a, b

CASES = [
    ('scrambled_baseline_joint', scrambled_baseline_joint),
    ('scrambled_baseline_indep', scrambled_baseline_indep),
    ('scrambled_OU_joint',       scrambled_OU_joint),
    ('near_constant',            near_constant),
    ('identical_streams',        identical_streams),
    ('mixed_marginals',          mixed_marginals),
    ('random_walk_indep',        random_walk_indep),
]

# ---------- extraction (exact match to extract_v1) ----------

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
    return dict(null_full=Vt[-1], singular_values=S, M_mean=M.mean(axis=0),
                M_shape=M.shape)

def project_3d(v6):
    sub = v6[[2, 3, 4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    if sub[2] < 0:
        sub = -sub
    return sub

def cos_sim(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return float('nan')
    return float(np.dot(u, v) / (nu * nv))

REF_MINUS_MINUS_PLUS = np.array([-1.0, -1.0, +2.0]) / np.sqrt(6.0)
REF_PLUS_PLUS_PLUS   = np.array([+1.0, +1.0, +2.0]) / np.sqrt(6.0)

# ---------- main ----------

def run_case(name, fn, N_REAL, seeds, T, dt, window):
    print()
    print("=" * 72)
    print(f"  CASE: {name}")
    print("=" * 72)
    per_seed = []
    for s in seeds:
        r = extract_v1(fn, N_REAL=N_REAL, master_seed=s, T=T, dt=dt, window=window)
        v3 = project_3d(r['null_full'])
        per_seed.append(dict(
            seed=s,
            null_full=r['null_full'].tolist(),
            v3=v3.tolist(),
            singular_values=r['singular_values'].tolist(),
            M_mean=r['M_mean'].tolist(),
            M_shape=list(r['M_shape']),
        ))
        print(f"  seed={s}: v3=({v3[0]:+.4f}, {v3[1]:+.4f}, {v3[2]:+.4f})")
        print(f"          null_full={r['null_full'].round(4).tolist()}")
        print(f"          SVs={[f'{x:.3e}' for x in r['singular_values']]}")
        print(f"          op_means={r['M_mean'].round(4).tolist()}")

    v3_stack = np.array([p['v3'] for p in per_seed])
    mean_v3 = v3_stack.mean(axis=0)
    norm = np.linalg.norm(mean_v3)
    if norm > 0:
        mean_v3 = mean_v3 / norm

    pair_cos = []
    for i in range(len(per_seed)):
        for j in range(i+1, len(per_seed)):
            pair_cos.append(cos_sim(per_seed[i]['v3'], per_seed[j]['v3']))

    cos_mmp = cos_sim(mean_v3, REF_MINUS_MINUS_PLUS)
    cos_ppp = cos_sim(mean_v3, REF_PLUS_PLUS_PLUS)
    print(f"  mean_v3:                ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
    print(f"  cos to (+1,+1,+2)/sqrt6: {cos_ppp:+.4f}")
    print(f"  cos to (-1,-1,+2)/sqrt6: {cos_mmp:+.4f}")
    if pair_cos:
        print(f"  inter-seed cos pairs:    {[f'{c:+.4f}' for c in pair_cos]}")
        print(f"  mean inter-seed cos:     {np.mean(pair_cos):+.4f}")

    return dict(
        per_seed=per_seed,
        mean_v3=mean_v3.tolist(),
        cos_to_plus_plus_plus=cos_ppp,
        cos_to_minus_minus_plus=cos_mmp,
        inter_seed_cos_pairs=pair_cos,
        inter_seed_cos_mean=float(np.mean(pair_cos)) if pair_cos else float('nan'),
    )

def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'full'
    if mode == 'canary':
        N_REAL = 50
        seeds = (0,)
        cases = [CASES[0], CASES[3]]
        out_path = '/home/user/Basic_equations/probe_null_inputs_canary.json'
    else:
        N_REAL = 500
        seeds = (0, 1, 2)
        cases = CASES
        out_path = '/home/user/Basic_equations/probe_null_inputs_results.json'

    print(f"mode={mode}  N_REAL={N_REAL}  seeds={seeds}  cases={len(cases)}")
    print(f"out={out_path}")

    t_start = time.time()
    out = {}
    for name, fn in cases:
        t0 = time.time()
        out[name] = run_case(name, fn, N_REAL, seeds, T=100, dt=0.1, window=40)
        print(f"  elapsed: {time.time()-t0:.1f}s")

    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print()
    print(f"  Total runtime: {time.time()-t_start:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
