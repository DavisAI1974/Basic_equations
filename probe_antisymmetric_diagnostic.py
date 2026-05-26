"""
Probe: INFO-017 antisymmetric energy fraction across systems and cases.

INFO-017 (Session 3) defines a portable diagnostic on the algebraic basis
  e1 = (H_a - H_b) / sqrt(2)
  e2 = (H_a^2 - H_b^2) / sqrt(2)
The antisymmetric energy fraction is the fraction of the centered operator
matrix's variance that lives in this 2D antisymmetric subspace, vs the
symmetric subspace.

Scales as 1/N_eff per INFO-017. Greg flagged it as a portable diagnostic
for finite-effective-sample noise.

Runs across:
  - Six cases from the earlier suite (linear drift baseline, no-drift,
    asymmetric drift, three coupling levels)
  - Eight systems from probe_system_directions (OU, uniform, Gauss-small,
    Gauss-large, Student-t, Laplace, sine, logistic)

Per case: anti_frac, sym_frac, and full per-window operator stats.

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

# ---------------------------------------------------------------------------
# system generators
# ---------------------------------------------------------------------------
def linear_drift(T, dt, v_a, v_b, seed, sigma=0.1, rho=0.0):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    eps_a = rng.normal(0, sigma, n)
    eps_b_i = rng.normal(0, sigma, n)
    eps_b = rho*eps_a + np.sqrt(max(0, 1-rho**2))*eps_b_i
    return v_a*t + eps_a, v_b*t + eps_b

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
    return rng.uniform(-1, 1, n), rng.uniform(-1, 1, n)

def gauss_pair(T, dt, seed, sigma):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.normal(0, sigma, n), rng.normal(0, sigma, n)

def studentt_pair(T, dt, seed, df=3):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.standard_t(df, n), rng.standard_t(df, n)

def laplace_pair(T, dt, seed, scale=1.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    return rng.laplace(0, scale, n), rng.laplace(0, scale, n)

def sine_pair(T, dt, seed, omega=0.5, phi=0.7, noise=0.05):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    t = np.arange(n) * dt
    return np.sin(omega*t) + noise*rng.normal(size=n), \
           np.sin(omega*t + phi) + noise*rng.normal(size=n)

def logistic_pair(T, dt, seed, r=3.9):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    a = np.zeros(n); b = np.zeros(n)
    a[0] = rng.uniform(0.1, 0.9); b[0] = rng.uniform(0.1, 0.9)
    for k in range(1, n):
        a[k] = r*a[k-1]*(1 - a[k-1])
        b[k] = r*b[k-1]*(1 - b[k-1])
    return a, b

# ---------------------------------------------------------------------------
# windowed operator collection
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

def collect_M(system_fn, N_REAL, master_seed, T=100, dt=0.1, window=40):
    all_ops = []
    for r in range(N_REAL):
        seed = master_seed * 1_000_000 + r
        a, b = system_fn(T, dt, seed)
        all_ops.append(windowed_operators(a, b, window, dt))
    return np.vstack(all_ops)

# ---------------------------------------------------------------------------
# antisymmetric energy fraction (INFO-017)
# ---------------------------------------------------------------------------
def antisymmetric_energy(M):
    """
    M has columns [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI].
    Compute energy in algebraic antisymmetric basis:
       e1 = (col_H_a - col_H_b)/sqrt(2)
       e2 = (col_H_a^2 - col_H_b^2)/sqrt(2)
    vs symmetric:
       s1 = (col_H_a + col_H_b)/sqrt(2)
       s2 = (col_H_a^2 + col_H_b^2)/sqrt(2)
    Also keep H_a*H_b (s-only) and MI separately for completeness.
    """
    Mc = M - M.mean(axis=0, keepdims=True)
    H_a, H_b, H_a2, H_b2, _, _ = [Mc[:, i] for i in range(6)]
    e1 = (H_a - H_b) / np.sqrt(2)
    e2 = (H_a2 - H_b2) / np.sqrt(2)
    s1 = (H_a + H_b) / np.sqrt(2)
    s2 = (H_a2 + H_b2) / np.sqrt(2)
    anti_energy = np.sum(e1**2) + np.sum(e2**2)
    sym_energy  = np.sum(s1**2) + np.sum(s2**2)
    total = anti_energy + sym_energy
    anti_frac = anti_energy / total if total > 0 else 0.0
    sym_frac  = sym_energy  / total if total > 0 else 0.0
    return dict(
        anti_energy=float(anti_energy),
        sym_energy=float(sym_energy),
        anti_frac=float(anti_frac),
        sym_frac=float(sym_frac),
        e1_var=float(np.var(e1)),
        e2_var=float(np.var(e2)),
        s1_var=float(np.var(s1)),
        s2_var=float(np.var(s2)),
    )

# ---------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------
def make_linear_drift_case(v_a, v_b, sigma, rho):
    return lambda T, dt, seed: linear_drift(T, dt, v_a, v_b, seed, sigma=sigma, rho=rho)

CASES = {
    'linear_drift_baseline_v1':       make_linear_drift_case(1.0, 1.0, 0.1, 0.0),
    'linear_drift_gauss_v0':          make_linear_drift_case(0.0, 0.0, 0.1, 0.0),
    'linear_drift_asym_v1v2':         make_linear_drift_case(1.0, 2.0, 0.1, 0.0),
    'linear_drift_rho_0p3':           make_linear_drift_case(1.0, 1.0, 0.1, 0.3),
    'linear_drift_rho_0p7':           make_linear_drift_case(1.0, 1.0, 0.1, 0.7),
    'linear_drift_rho_0p95':          make_linear_drift_case(1.0, 1.0, 0.1, 0.95),
    'OU_gamma1_sigma1':               lambda T, dt, seed: ou_pair(T, dt, seed, 1.0, 1.0),
    'uniform_pm1':                    lambda T, dt, seed: uniform_pair(T, dt, seed),
    'gauss_small_0p1':                lambda T, dt, seed: gauss_pair(T, dt, seed, 0.1),
    'gauss_large_1p0':                lambda T, dt, seed: gauss_pair(T, dt, seed, 1.0),
    'student_t_df3':                  lambda T, dt, seed: studentt_pair(T, dt, seed, 3),
    'laplace_scale1':                 lambda T, dt, seed: laplace_pair(T, dt, seed, 1.0),
    'sine_w0p5_phi0p7':               lambda T, dt, seed: sine_pair(T, dt, seed, 0.5, 0.7, 0.05),
    'logistic_r3p9':                  lambda T, dt, seed: logistic_pair(T, dt, seed, 3.9),
}

def main():
    t_start = time.time()
    seeds = (0, 1, 2)
    N_REAL = 500

    out = {}
    for name, fn in CASES.items():
        print()
        print("=" * 72)
        print(f"  CASE: {name}")
        print("=" * 72)
        t0 = time.time()
        per_seed = []
        for s in seeds:
            M = collect_M(fn, N_REAL, master_seed=s)
            diag = antisymmetric_energy(M)
            per_seed.append(dict(seed=s, **diag))
            print(f"  seed={s}: anti_frac={diag['anti_frac']:.5f}  "
                  f"sym_frac={diag['sym_frac']:.5f}  "
                  f"e1_var={diag['e1_var']:.3e}  e2_var={diag['e2_var']:.3e}  "
                  f"s1_var={diag['s1_var']:.3e}  s2_var={diag['s2_var']:.3e}")
        anti_mean = np.mean([p['anti_frac'] for p in per_seed])
        print(f"  mean anti_frac across seeds: {anti_mean:.5f}")
        print(f"  elapsed: {time.time()-t0:.1f}s")
        out[name] = dict(per_seed=per_seed, anti_frac_mean=float(anti_mean))

    out_path = '/home/user/Basic_equations/probe_antisymmetric_diagnostic_results.json'
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print()
    print(f"  Total runtime: {time.time()-t_start:.1f}s")
    print(f"  Results: {out_path}")

if __name__ == "__main__":
    main()
