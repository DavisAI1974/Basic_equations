"""
Probe (Thread 3): high-seed convergence on scrambled_baseline_joint.

In the null-inputs probe (Session 5), this was the only "reasonable" case
(no degenerate input, no marginal mismatch) where the protocol drifted
off (+1,+1,+2)/sqrt(6). At 3 seeds, inter-seed cos was only +0.40 and
the mean direction landed closer to (-1,-1,+2)/sqrt(6) — the Session 3
INFO-014 reference. This probe asks: where does it actually want to go
when given 30 seeds, and does it converge?

Exact extract_v1 protocol. 30 seeds, N_REAL=500. Reports per-seed v3,
mean v3, inter-seed cos histogram, cos to (+,+,+) and (-,-,+) references.
"""
import json
import time
import numpy as np
from scipy.stats import differential_entropy

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

def linear_drift_baseline(T, dt, seed, v=1.0, noise_std=0.1):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt); n = len(t)
    return v*t+rng.normal(0,noise_std,n), v*t+rng.normal(0,noise_std,n)

def scrambled_baseline_joint(T, dt, seed):
    rng = np.random.default_rng(seed ^ 0xABCDEF)
    a, b = linear_drift_baseline(T, dt, seed)
    perm = rng.permutation(len(a))
    return a[perm], b[perm]

def windowed_operators(a, b, window, dt, overlap=0.5):
    ws = int(window/dt); step = int(ws*(1-overlap))
    ops, start = [], 0
    while start+ws <= len(a):
        ops.append(operators(a[start:start+ws], b[start:start+ws]))
        start += step
    return np.array(ops)

def extract_v1(fn, N_REAL, master_seed, T=100, dt=0.1, window=40):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = fn(T, dt, seed)
        all_ops.append(windowed_operators(a, b, window, dt))
    M = np.vstack(all_ops)
    M_c = M - M.mean(0, keepdims=True)
    _, S, Vt = np.linalg.svd(M_c, full_matrices=False)
    return Vt[-1], S

def project_3d(v6):
    sub = v6[[2,3,4]]
    n = np.linalg.norm(sub)
    sub = sub / (n if n > 0 else 1.0)
    if sub[2] < 0: sub = -sub
    return sub

REF_PPP = np.array([1,1,2])/np.sqrt(6)
REF_MMP = np.array([-1,-1,2])/np.sqrt(6)

def cos_sim(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    return float(np.dot(u, v)/(nu*nv)) if nu*nv > 0 else float('nan')

def main():
    n_seeds = 30
    N_REAL = 500
    t0 = time.time()
    per_seed = []
    print(f"Running scrambled_baseline_joint, N_REAL={N_REAL}, {n_seeds} seeds...")
    for s in range(n_seeds):
        v6, S = extract_v1(scrambled_baseline_joint, N_REAL, s)
        v3 = project_3d(v6)
        per_seed.append(dict(
            seed=s, null_full=v6.tolist(), v3=v3.tolist(),
            singular_values=S.tolist(),
            cos_ppp=cos_sim(v3, REF_PPP), cos_mmp=cos_sim(v3, REF_MMP),
        ))
        print(f"  s={s:2d}: v3=({v3[0]:+.3f},{v3[1]:+.3f},{v3[2]:+.3f})  "
              f"cos(+,+,+)={per_seed[-1]['cos_ppp']:+.3f}  "
              f"cos(-,-,+)={per_seed[-1]['cos_mmp']:+.3f}  "
              f"SV3/SV4={S[2]/S[3]:.2e}")

    v3_stack = np.array([p['v3'] for p in per_seed])
    mean_v3 = v3_stack.mean(0)
    n = np.linalg.norm(mean_v3)
    mean_v3 = mean_v3 / n if n > 0 else mean_v3
    pair_cos = []
    for i in range(n_seeds):
        for j in range(i+1, n_seeds):
            pair_cos.append(cos_sim(per_seed[i]['v3'], per_seed[j]['v3']))
    pair_cos = np.array(pair_cos)

    print()
    print(f"mean v3 across {n_seeds} seeds: ({mean_v3[0]:+.4f}, {mean_v3[1]:+.4f}, {mean_v3[2]:+.4f})")
    print(f"cos to (+,+,+)/sqrt6: {cos_sim(mean_v3, REF_PPP):+.4f}")
    print(f"cos to (-,-,+)/sqrt6: {cos_sim(mean_v3, REF_MMP):+.4f}")
    print(f"per-coord SD across seeds: "
          f"({v3_stack.std(0)[0]:.3f}, {v3_stack.std(0)[1]:.3f}, {v3_stack.std(0)[2]:.3f})")
    print(f"inter-seed cos: mean={pair_cos.mean():+.3f}  median={np.median(pair_cos):+.3f}  "
          f"min={pair_cos.min():+.3f}  max={pair_cos.max():+.3f}")
    print(f"inter-seed cos histogram bins [-1,-0.5,0,0.5,1]: "
          f"{np.histogram(pair_cos, bins=[-1, -0.5, 0, 0.5, 1.0])[0].tolist()}")

    # also check distribution of per-seed cos to refs
    cos_ppp_arr = np.array([p['cos_ppp'] for p in per_seed])
    cos_mmp_arr = np.array([p['cos_mmp'] for p in per_seed])
    print()
    print(f"per-seed cos to (+,+,+): mean={cos_ppp_arr.mean():+.3f}  "
          f"SD={cos_ppp_arr.std():.3f}  min={cos_ppp_arr.min():+.3f}  max={cos_ppp_arr.max():+.3f}")
    print(f"per-seed cos to (-,-,+): mean={cos_mmp_arr.mean():+.3f}  "
          f"SD={cos_mmp_arr.std():.3f}  min={cos_mmp_arr.min():+.3f}  max={cos_mmp_arr.max():+.3f}")

    out = dict(
        per_seed=per_seed,
        mean_v3=mean_v3.tolist(),
        cos_to_plus_plus_plus=cos_sim(mean_v3, REF_PPP),
        cos_to_minus_minus_plus=cos_sim(mean_v3, REF_MMP),
        per_coord_sd=v3_stack.std(0).tolist(),
        inter_seed_cos_pairs=pair_cos.tolist(),
        inter_seed_cos_mean=float(pair_cos.mean()),
        inter_seed_cos_median=float(np.median(pair_cos)),
        per_seed_cos_ppp_arr=cos_ppp_arr.tolist(),
        per_seed_cos_mmp_arr=cos_mmp_arr.tolist(),
    )
    with open('/home/user/Basic_equations/probe_scrambled_highseed_results.json', 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nTotal: {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
