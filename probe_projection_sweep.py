"""
Probe (Thread 2): projection-rule sweep.

Same baseline data (linear drift, N_REAL=500, 3 seeds). Same protocol
through SVD. Same smallest right singular vector. But project the 6D
null vector to different 3D subspaces and report direction in each.

If (+,+,+) is a feature of subspace [2,3,4] alone, other subspaces
should give different shapes. If it's a feature of the data, multiple
subspaces should show something analogous.

For each subspace of size 3: sign-canonicalize so the last coordinate
is positive, report mean direction across seeds, inter-seed cos, and
cos to several reference patterns:
   (+1, +1, +2)/sqrt(6)   (canonical)
   (+1, +1, +1)/sqrt(3)   (symmetric simplex)
   (+1, +1, 0)/sqrt(2)    (two-axis sum)
   (0, 0, +1)             (last axis only)
"""
import json
import time
import numpy as np
from itertools import combinations
from scipy.stats import differential_entropy

OP_NAMES = ['H_a', 'H_b', 'H_a^2', 'H_b^2', 'H_a*H_b', 'MI']

# --- protocol pieces (exact match) ---
def entropy_1d(x): return differential_entropy(x, method='vasicek')
def entropy_2d_kde(xy, n_grid=30):
    H, xe, ye = np.histogram2d(xy[:,0], xy[:,1], bins=n_grid, density=True)
    dx, dy = xe[1]-xe[0], ye[1]-ye[0]
    p = H.flatten(); p = p[p>0]
    return -np.sum(p*np.log(p))*dx*dy
def mutual_info(a, b):
    return entropy_1d(a)+entropy_1d(b)-entropy_2d_kde(np.column_stack([a,b]))
def operators(a, b):
    Ha, Hb = entropy_1d(a), entropy_1d(b)
    return np.array([Ha, Hb, Ha**2, Hb**2, Ha*Hb, mutual_info(a,b)])
def linear_drift(T, dt, v, seed, ns=0.1):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt); n = len(t)
    return v*t+rng.normal(0,ns,n), v*t+rng.normal(0,ns,n)
def windowed_operators(a, b, window, dt, overlap=0.5):
    ws = int(window/dt); step = int(ws*(1-overlap))
    ops, start = [], 0
    while start+ws <= len(a):
        ops.append(operators(a[start:start+ws], b[start:start+ws]))
        start += step
    return np.array(ops)
def extract_v1_null(N_REAL, master_seed, T=100, dt=0.1, window=40, v=1.0):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = linear_drift(T, dt, v, seed)
        all_ops.append(windowed_operators(a, b, window, dt))
    M = np.vstack(all_ops)
    M_c = M - M.mean(0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return Vt[-1], S, M.mean(0)

# --- projection + analysis ---
def project_to_subspace(v6, idx_tuple):
    sub = v6[list(idx_tuple)]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    # sign-canonicalize: last component positive
    if sub[-1] < 0:
        sub = -sub
    return sub

REFS = {
    '(+1,+1,+2)/sqrt6': np.array([1,1,2])/np.sqrt(6),
    '(+1,+1,+1)/sqrt3': np.array([1,1,1])/np.sqrt(3),
    '(+1,+1, 0)/sqrt2': np.array([1,1,0])/np.sqrt(2),
    '( 0, 0,+1)':       np.array([0,0,1]),
}

def cos_to_refs(v3):
    return {k: float(np.dot(v3, r)) for k, r in REFS.items()}

# --- main ---
def main():
    t0 = time.time()
    seeds = (0, 1, 2)

    # collect null vectors per seed
    null_vecs = []
    for s in seeds:
        v6, S, op_means = extract_v1_null(500, s)
        null_vecs.append(v6)
        print(f"seed={s}: null_full={v6.round(4).tolist()}")
        print(f"        SVs={[f'{x:.3e}' for x in S]}")
        print(f"        op_means={op_means.round(4).tolist()}")

    subspaces = list(combinations(range(6), 3))  # 20 subspaces

    out = {}
    print()
    print(f"{'subspace':<25} {'mean v3 (canonical)':<35} {'cos(+,+,+2)':<10} {'cos(+++)':<10} {'cos(++0)':<10} {'cos(001)':<10} {'inter-seed':<10}")
    print("-" * 130)
    for idx in subspaces:
        label = '[' + ','.join(OP_NAMES[i] for i in idx) + ']'
        v3_per_seed = [project_to_subspace(v6, idx) for v6 in null_vecs]
        v3_stack = np.array(v3_per_seed)
        # global sign-canonicalize across seeds: flip seeds with neg inner-product to mean
        ref = v3_stack[0]
        for k in range(1, len(v3_stack)):
            if np.dot(v3_stack[k], ref) < 0:
                v3_stack[k] = -v3_stack[k]
        mean_v3 = v3_stack.mean(0)
        n = np.linalg.norm(mean_v3)
        mean_v3 = mean_v3 / n if n > 0 else mean_v3
        cs = cos_to_refs(mean_v3)
        # inter-seed agreement
        pair_cos = []
        for i in range(len(v3_stack)):
            for j in range(i+1, len(v3_stack)):
                pair_cos.append(abs(np.dot(v3_stack[i], v3_stack[j])))
        inter = float(np.mean(pair_cos))

        print(f"{label:<25} ({mean_v3[0]:+.3f},{mean_v3[1]:+.3f},{mean_v3[2]:+.3f})       "
              f"{cs['(+1,+1,+2)/sqrt6']:+.3f}     "
              f"{cs['(+1,+1,+1)/sqrt3']:+.3f}     "
              f"{cs['(+1,+1, 0)/sqrt2']:+.3f}     "
              f"{cs['( 0, 0,+1)']:+.3f}     "
              f"{inter:+.3f}")
        out[label] = dict(
            indices=list(idx),
            mean_v3=mean_v3.tolist(),
            per_seed_v3=v3_stack.tolist(),
            cos_to_refs=cs,
            inter_seed_cos_mean=inter,
        )

    with open('/home/user/Basic_equations/probe_projection_sweep_results.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nTotal: {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
