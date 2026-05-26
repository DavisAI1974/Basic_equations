"""
Probe (Thread 1): operator-noise bypass.

Skip the input/window/entropy path entirely. Generate the M matrix
(shape = N_pool x 6) directly from a chosen distribution. Center, SVD,
take smallest right singular vector, project to indices [2,3,4],
sign-canonicalize (coord 2 positive). Report direction and cos to
(+1,+1,+2)/sqrt(6) and (-1,-1,+2)/sqrt(6).

If (+,+,+) appears under random M, the direction is encoded in the
post-SVD geometry alone, not in any input structure.

Variants:
  iid_unit          M_ij ~ N(0, 1)            iid in row and column
  iid_col_scaled    columns scaled to match baseline operator SDs
                    (so SVD doesn't trivially pick the smallest-scale coord)
  cov_matched       M ~ N(0, Sigma) where Sigma = sample cov of baseline M
                    (preserves the inter-operator covariance from real data)

Baseline M source: rerun extract_v1 on linear-drift baseline N_REAL=500,
master_seed=0, exact protocol. Take the resulting pool, compute marginal
SDs and full covariance for the scaled / cov-matched variants.

10 seeds per variant, N_pool = 2000 rows (matches baseline 500 x 4 windows).
"""
import json
import time
import numpy as np
from scipy.stats import differential_entropy

# ---------- protocol pieces (exact match) ----------

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
    return entropy_1d(a) + entropy_1d(b) - entropy_2d_kde(np.column_stack([a, b]))

def operators(a, b):
    H_a = entropy_1d(a)
    H_b = entropy_1d(b)
    MI = mutual_info(a, b)
    return np.array([H_a, H_b, H_a**2, H_b**2, H_a*H_b, MI])

def linear_drift(T, dt, v, seed, noise_std=0.1):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    return v*t + rng.normal(0, noise_std, n), v*t + rng.normal(0, noise_std, n)

def windowed_operators(a, b, window, dt, overlap=0.5):
    win_samples = int(window / dt)
    step = int(win_samples * (1 - overlap))
    ops = []
    start = 0
    while start + win_samples <= len(a):
        ops.append(operators(a[start:start+win_samples], b[start:start+win_samples]))
        start += step
    return np.array(ops)

def build_baseline_M(N_REAL=500, master_seed=0, T=100, dt=0.1, window=40, v=1.0):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = linear_drift(T, dt, v, seed)
        all_ops.append(windowed_operators(a, b, window, dt))
    return np.vstack(all_ops)  # (N_REAL * n_windows, 6)

# ---------- projection ----------

REF_PPP = np.array([+1.0, +1.0, +2.0]) / np.sqrt(6.0)
REF_MMP = np.array([-1.0, -1.0, +2.0]) / np.sqrt(6.0)

def smallest_rsv(M):
    M_c = M - M.mean(axis=0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return Vt[-1], S

def project_3d(v6):
    sub = v6[[2, 3, 4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    if sub[2] < 0:
        sub = -sub
    return sub

def cos_sim(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(np.dot(u, v) / (nu*nv)) if nu*nv > 0 else float('nan')

# ---------- variants ----------

def gen_iid_unit(rng, n_rows):
    return rng.normal(0, 1, size=(n_rows, 6))

def gen_iid_col_scaled(rng, n_rows, col_sds):
    M = rng.normal(0, 1, size=(n_rows, 6))
    return M * col_sds  # broadcast (1,6)

def gen_cov_matched(rng, n_rows, cov):
    return rng.multivariate_normal(mean=np.zeros(6), cov=cov, size=n_rows)

# ---------- main ----------

def run_variant(name, gen_fn, n_seeds=10, n_rows=2000, **kwargs):
    print()
    print("=" * 72)
    print(f"  VARIANT: {name}")
    print("=" * 72)
    per_seed = []
    for s in range(n_seeds):
        rng = np.random.default_rng(1_000_000 * 13 + s)  # offset from data seeds
        M = gen_fn(rng, n_rows, **kwargs)
        v6, S = smallest_rsv(M)
        v3 = project_3d(v6)
        per_seed.append(dict(
            seed=s, null_full=v6.tolist(), v3=v3.tolist(),
            singular_values=S.tolist(),
        ))
        print(f"  seed={s}: v3=({v3[0]:+.4f}, {v3[1]:+.4f}, {v3[2]:+.4f})  "
              f"SV_gap={S[2]/S[3]:.2e}")

    v3_stack = np.array([p['v3'] for p in per_seed])
    mean_v3 = v3_stack.mean(axis=0)
    n = np.linalg.norm(mean_v3)
    mean_v3 = mean_v3 / n if n > 0 else mean_v3
    pair_cos = []
    for i in range(n_seeds):
        for j in range(i+1, n_seeds):
            pair_cos.append(cos_sim(per_seed[i]['v3'], per_seed[j]['v3']))
    cos_p = cos_sim(mean_v3, REF_PPP)
    cos_m = cos_sim(mean_v3, REF_MMP)
    print(f"  mean_v3:                ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
    print(f"  cos to (+1,+1,+2)/sqrt6: {cos_p:+.4f}")
    print(f"  cos to (-1,-1,+2)/sqrt6: {cos_m:+.4f}")
    print(f"  mean inter-seed cos:     {np.mean(pair_cos):+.4f}")
    print(f"  per-coord SD across seeds: "
          f"({v3_stack.std(0)[0]:.3f}, {v3_stack.std(0)[1]:.3f}, {v3_stack.std(0)[2]:.3f})")
    return dict(
        per_seed=per_seed,
        mean_v3=mean_v3.tolist(),
        cos_to_plus_plus_plus=cos_p,
        cos_to_minus_minus_plus=cos_m,
        inter_seed_cos_mean=float(np.mean(pair_cos)),
        per_coord_sd=v3_stack.std(0).tolist(),
    )

def main():
    t0 = time.time()
    print("Building baseline M (linear drift, N_REAL=500, master_seed=0)...")
    M_base = build_baseline_M()
    print(f"  M shape: {M_base.shape}")
    print(f"  op_means: {M_base.mean(0).round(4).tolist()}")
    col_sds = M_base.std(0, ddof=1)
    print(f"  col_sds:  {col_sds.round(4).tolist()}")
    cov = np.cov(M_base, rowvar=False)
    print(f"  cov diag: {np.diag(cov).round(4).tolist()}")

    out = {
        'baseline_M_shape': list(M_base.shape),
        'baseline_op_means': M_base.mean(0).tolist(),
        'baseline_col_sds': col_sds.tolist(),
        'baseline_cov': cov.tolist(),
    }
    out['iid_unit']       = run_variant('iid_unit',       gen_iid_unit)
    out['iid_col_scaled'] = run_variant('iid_col_scaled', gen_iid_col_scaled, col_sds=col_sds)
    out['cov_matched']    = run_variant('cov_matched',    gen_cov_matched, cov=cov)

    out_path = '/home/user/Basic_equations/probe_operator_noise_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n  Total: {time.time()-t0:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
