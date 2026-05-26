"""
extract_v1 protocol — linear drift N_REAL sweep.

Disambiguates INFO-016: basis structural limit vs sample noise.

Protocol (per SESSION_HANDOFF_2026-05-25_v3.md):
- Operator basis: {H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI}
- Window: 40s, dt=0.1, T=100 per realization
- Pool operator matrices across N_REAL realizations
- Null direction = smallest right singular vector of centered pooled matrix
- Project to (H_a^2, H_b^2, H_a*H_b) subspace
- Compare to attractor direction (-1, -1, +2)/sqrt(6)
- N_REAL in {50, 100, 200, 500}, 3 seeds each
- Linear drift system: x_a(t) = v*t + noise_a, x_b(t) = v*t + noise_b, v=1
"""
import numpy as np
from scipy.stats import differential_entropy
import json
import time

# ---------- Operator computation ----------

def entropy_1d(x):
    """Vasicek-style 1D differential entropy (fast)."""
    return differential_entropy(x, method='vasicek')

def entropy_2d_kde(xy, n_grid=30):
    """2D differential entropy via histogram (fast, biased but consistent)."""
    H, xedges, yedges = np.histogram2d(xy[:, 0], xy[:, 1], bins=n_grid, density=True)
    dx = xedges[1] - xedges[0]
    dy = yedges[1] - yedges[0]
    p = H.flatten()
    p = p[p > 0]
    return -np.sum(p * np.log(p)) * dx * dy

def mutual_info(a, b):
    """MI = H(a) + H(b) - H(a,b)."""
    H_a = entropy_1d(a)
    H_b = entropy_1d(b)
    H_ab = entropy_2d_kde(np.column_stack([a, b]))
    return H_a + H_b - H_ab

def operators(a, b):
    """6-vector: [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]"""
    H_a = entropy_1d(a)
    H_b = entropy_1d(b)
    MI = mutual_info(a, b)
    return np.array([H_a, H_b, H_a**2, H_b**2, H_a*H_b, MI])

# ---------- System: linear drift ----------

def linear_drift(T, dt, v, seed, noise_std=0.1):
    """x_a = v*t + noise_a, x_b = v*t + noise_b (shared drift, independent noise)."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt)
    n = len(t)
    a = v * t + rng.normal(0, noise_std, n)
    b = v * t + rng.normal(0, noise_std, n)
    return a, b

# ---------- Windowed-null extraction ----------

def windowed_operators(a, b, window, dt, overlap=0.5):
    """Compute operator vector in sliding windows with given overlap."""
    win_samples = int(window / dt)
    step = int(win_samples * (1 - overlap))
    n = len(a)
    ops_list = []
    start = 0
    while start + win_samples <= n:
        ops_list.append(operators(a[start:start + win_samples],
                                   b[start:start + win_samples]))
        start += step
    return np.array(ops_list)  # shape (n_windows, 6)

def extract_v1(N_REAL, master_seed, T=100, dt=0.1, window=40, v=1.0):
    """Pool windowed operator matrices across N_REAL realizations, extract null."""
    all_ops = []
    for r in range(N_REAL):
        # Distinct seed per realization, deterministic in master_seed
        seed = master_seed * 1_000_000 + r
        a, b = linear_drift(T, dt, v, seed)
        ops = windowed_operators(a, b, window, dt)
        all_ops.append(ops)
    M = np.vstack(all_ops)  # (N_REAL * n_windows_per_real, 6)

    # Center
    M_centered = M - M.mean(axis=0, keepdims=True)

    # SVD: smallest singular vector = null direction
    U, S, Vt = np.linalg.svd(M_centered, full_matrices=False)
    null_vec = Vt[-1]  # 6-vector

    return {
        'null_full': null_vec,
        'singular_values': S,
        'M_shape': M.shape,
    }

# ---------- Analysis ----------

ATTRACTOR_SUB = np.array([-1, -1, 2]) / np.sqrt(6)  # in (H_a^2, H_b^2, H_a*H_b)

def project_and_compare(null_vec):
    """Pull (H_a^2, H_b^2, H_a*H_b) coefficients from full 6-vector, normalize, compare to attractor."""
    # Index mapping: [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]
    sub = null_vec[[2, 3, 4]]
    sub_norm = sub / (np.linalg.norm(sub) + 1e-12)
    # Sign-canonicalize: flip so largest |component| is positive
    idx_max = np.argmax(np.abs(sub_norm))
    if sub_norm[idx_max] < 0:
        sub_norm = -sub_norm
    cos_attractor = float(np.dot(sub_norm, ATTRACTOR_SUB))
    # Also try with attractor flipped — the sign of the attractor is itself a convention
    cos_attractor_flipped = float(np.dot(sub_norm, -ATTRACTOR_SUB))
    cos_attractor_signed = max(cos_attractor, cos_attractor_flipped)
    return sub_norm, cos_attractor_signed

def inter_seed_cosine(vecs):
    """Pairwise cosines, with sign canonicalization."""
    cosines = []
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            c = np.dot(vecs[i], vecs[j])
            cosines.append(abs(c))  # sign-blind for "direction" similarity
    return cosines

# ---------- Smoke test ----------

if __name__ == "__main__":
    print("=== SMOKE TEST: linear drift N_REAL=10, seed=0 ===")
    t0 = time.time()
    result = extract_v1(N_REAL=10, master_seed=0)
    print(f"M shape: {result['M_shape']}")
    print(f"Singular values: {result['singular_values']}")
    print(f"Null vector (full): {result['null_full']}")
    sub, cos_a = project_and_compare(result['null_full'])
    print(f"Subspace (H_a^2, H_b^2, H_a*H_b) normalized: {sub}")
    print(f"Cos to attractor (-1,-1,+2)/sqrt(6): {cos_a:.4f}")
    print(f"Wall time: {time.time() - t0:.1f}s")
