"""
LEVEL 2 PROBE - FOUR SCIENCES NETWORK TEST
============================================
Question: Is there structure BELOW the information dipole?

Previous Level 2 attempt (meta_information_layer.py) used
entropy-of-entropy observables. R² collapsed to 0.02-0.17.
Logged as OPEN, not closed: wrong tool, not empty floor.

This script applies the SAME five-step process that found Level 1,
but with a different observable class:

  STEP 1. New observable class: time-resolved dipole balance B_ij(t)
          computed across a NETWORK of subsystems, not nested entropy.
          B_ij(t) = |c_{H_i^2}| / |c_{H_i*H_j}| in sliding windows.
  STEP 2. Treat B_ij(t) as autonomous state variables.
  STEP 3. Build network operator library: B_k, B_k^2, B_net (mean of
          neighbor dipoles), B_net^2, B_k*B_net, <B^2>, const.
  STEP 4. Solve dB_k/dt = library @ coefficients per pair per science.
  STEP 5. Look for universal opposing structure across all 4 sciences.

Four sciences, each as a chain of N=6 coupled subsystems:
  PHYSICS:    Coupled Duffing oscillators with neighbor coupling
  BIOLOGY:    Coupled Lotka-Volterra patches with prey migration
  CHEMISTRY:  Coupled Brusselator cells with diffusion
  GEOLOGY:    Coupled Burridge-Knopoff fault segments

Falsification-first:
  - No universal opposing pair across 4 sciences -> Level 2 either
    truly absent in this observable class OR needs yet another lens.
  - Universal opposing pair in >=3/4 sciences -> candidate Level 2
    dipole, requires real-data validation (NHANES, PTB-XL).
"""
import numpy as np
from numpy.linalg import lstsq
from itertools import combinations
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
T_START = time.time()

# ============================================================
# INFO-THEORETIC PRIMITIVES
# ============================================================
def entropy_1d(samples, bins=30):
    h, e = np.histogram(samples, bins=bins, density=True)
    h = h[h > 0]
    dx = e[1] - e[0]
    return -np.sum(h * np.log(h) * dx)

def mi_2d(x, y, bins=18):
    hx = entropy_1d(x, bins)
    hy = entropy_1d(y, bins)
    h2, ex, ey = np.histogram2d(x, y, bins=bins, density=True)
    mask = h2 > 0
    if not mask.any():
        return 0.0
    dA = (ex[1] - ex[0]) * (ey[1] - ey[0])
    hxy = -np.sum(h2[mask] * np.log(h2[mask])) * dA
    return max(0.0, hx + hy - hxy)

# ============================================================
# FOUR NETWORK SIMULATIONS (N=6 coupled nodes each)
# ============================================================
def simulate_physics(N=6, N_ens=500, T=30.0, dt=0.02):
    steps = int(T / dt)
    x = 0.5 * np.random.randn(N_ens, N)
    v = 0.1 * np.random.randn(N_ens, N)
    K = 0.15
    out = np.zeros((steps, N_ens, N))
    for t in range(steps):
        nbr = np.zeros_like(x)
        nbr[:, 1:-1] = x[:, :-2] + x[:, 2:] - 2 * x[:, 1:-1]
        nbr[:, 0] = x[:, 1] - x[:, 0]
        nbr[:, -1] = x[:, -2] - x[:, -1]
        F = -x - 0.5 * x**3 + 0.1 * np.sin(0.5 * t * dt) + K * nbr
        v = v + F * dt - 0.05 * v * dt
        x = x + v * dt + 0.02 * np.random.randn(*x.shape) * np.sqrt(dt)
        out[t] = x
    return out, dt

def simulate_biology(N=6, N_ens=500, T=30.0, dt=0.02):
    steps = int(T / dt)
    prey = 1.0 + 0.1 * np.random.randn(N_ens, N)
    pred = 0.5 + 0.1 * np.random.randn(N_ens, N)
    K = 0.12
    out = np.zeros((steps, N_ens, N))
    for t in range(steps):
        nbr_prey = np.zeros_like(prey)
        nbr_prey[:, 1:-1] = prey[:, :-2] + prey[:, 2:] - 2 * prey[:, 1:-1]
        nbr_prey[:, 0] = prey[:, 1] - prey[:, 0]
        nbr_prey[:, -1] = prey[:, -2] - prey[:, -1]
        dprey = (1.0 * prey - 0.5 * prey * pred + K * nbr_prey) * dt
        dpred = (0.3 * prey * pred - 0.4 * pred) * dt
        prey = np.clip(prey + dprey + 0.01 * np.random.randn(*prey.shape) * np.sqrt(dt), 0.01, None)
        pred = np.clip(pred + dpred + 0.01 * np.random.randn(*pred.shape) * np.sqrt(dt), 0.01, None)
        out[t] = prey
    return out, dt

def simulate_chemistry(N=6, N_ens=500, T=30.0, dt=0.02):
    steps = int(T / dt)
    X = 1.0 + 0.1 * np.random.randn(N_ens, N)
    Y = 3.0 + 0.1 * np.random.randn(N_ens, N)
    A, B = 1.0, 3.0
    D = 0.18
    out = np.zeros((steps, N_ens, N))
    for t in range(steps):
        nbr_X = np.zeros_like(X)
        nbr_X[:, 1:-1] = X[:, :-2] + X[:, 2:] - 2 * X[:, 1:-1]
        nbr_X[:, 0] = X[:, 1] - X[:, 0]
        nbr_X[:, -1] = X[:, -2] - X[:, -1]
        dX = (A - (B + 1) * X + X**2 * Y + D * nbr_X) * dt
        dY = (B * X - X**2 * Y) * dt
        X = X + dX + 0.01 * np.random.randn(*X.shape) * np.sqrt(dt)
        Y = Y + dY + 0.01 * np.random.randn(*Y.shape) * np.sqrt(dt)
        out[t] = X
    return out, dt

def simulate_geology(N=6, N_ens=500, T=30.0, dt=0.02):
    steps = int(T / dt)
    # Burridge-Knopoff-like: chain of stick-slip blocks
    u = 0.1 * np.random.randn(N_ens, N)
    v = 0.05 * np.random.randn(N_ens, N)
    K = 0.4   # neighbor spring
    K0 = 0.6  # leaf spring to driver
    drive_rate = 0.04
    out = np.zeros((steps, N_ens, N))
    for t in range(steps):
        nbr = np.zeros_like(u)
        nbr[:, 1:-1] = u[:, :-2] + u[:, 2:] - 2 * u[:, 1:-1]
        nbr[:, 0] = u[:, 1] - u[:, 0]
        nbr[:, -1] = u[:, -2] - u[:, -1]
        # nonlinear friction
        fric = -0.6 * v / (1.0 + 8.0 * v**2)
        F = K * nbr - K0 * (u - drive_rate * t * dt) + fric
        v = v + F * dt
        u = u + v * dt + 0.015 * np.random.randn(*u.shape) * np.sqrt(dt)
        out[t] = u
    return out, dt

# ============================================================
# LEVEL 1 EXTRACTION (per-node H, pairwise MI over time)
# ============================================================
def compute_info_layer(states, dt, bins=18, subsample=2):
    T, N_ens, N = states.shape
    # subsample time to save compute (info layer varies slowly relative to dt)
    times = np.arange(0, T, subsample)
    pairs = list(combinations(range(N), 2))
    H = np.zeros((len(times), N))
    MI = np.zeros((len(times), len(pairs)))
    for idx, t in enumerate(times):
        for i in range(N):
            H[idx, i] = entropy_1d(states[t, :, i], bins)
        for k, (i, j) in enumerate(pairs):
            MI[idx, k] = mi_2d(states[t, :, i], states[t, :, j], bins)
    dt_info = dt * subsample
    return H, MI, pairs, dt_info

# ============================================================
# DIPOLE BALANCE TIME SERIES (per pair, sliding window OD)
# ============================================================
def extract_balance_timeseries(H, MI, pairs, dt_info, window=30, stride=4):
    T, N = H.shape
    n_pairs = len(pairs)
    n_windows = (T - window) // stride
    if n_windows < 20:
        return None, None
    B = np.zeros((n_windows, n_pairs))
    t_centers = np.zeros(n_windows)

    for w in range(n_windows):
        s = w * stride
        e = s + window
        for k, (i, j) in enumerate(pairs):
            Hi = H[s:e, i]
            Hj = H[s:e, j]
            Mij = MI[s:e, k]
            dHi = np.gradient(Hi, dt_info)
            # Level-1 dipole library: H_i, H_j, H_i^2, H_j^2, H_i*H_j, MI, const
            lib = np.column_stack([
                Hi, Hj, Hi**2, Hj**2, Hi * Hj, Mij, np.ones(len(Hi))
            ])
            c, *_ = lstsq(lib, dHi, rcond=None)
            num = abs(c[2])              # |c_{H_i^2}|
            den = abs(c[4]) + 1e-10      # |c_{H_i*H_j}|
            B[w, k] = num / den
        t_centers[w] = (s + window / 2.0) * dt_info
    return B, t_centers

# ============================================================
# LEVEL 2 OD ON THE DIPOLE NETWORK
# ============================================================
def od_on_dipole_network(B, t_centers):
    n_windows, n_pairs = B.shape
    if n_windows < 30:
        return None
    dt_b = t_centers[1] - t_centers[0]
    # log-transform to tame dynamic range of B
    Bl = np.log(B + 1e-6)

    labels = ['B_k', 'B_net', 'B_k^2', 'B_net^2', 'B_k*B_net', 'mean_B2', 'const']
    results = {}
    for k in range(n_pairs):
        dBk = np.gradient(Bl[:, k], dt_b)
        bk = Bl[:, k]
        mask = np.ones(n_pairs, bool)
        mask[k] = False
        b_net = Bl[:, mask].mean(axis=1)
        b_net2_mean = (Bl[:, mask] ** 2).mean(axis=1)
        lib = np.column_stack([
            bk, b_net, bk**2, b_net**2, bk * b_net, b_net2_mean, np.ones(len(bk))
        ])
        c, *_ = lstsq(lib, dBk, rcond=None)
        ss_res = np.sum((dBk - lib @ c) ** 2)
        ss_tot = np.sum((dBk - dBk.mean()) ** 2) + 1e-15
        r2 = 1.0 - ss_res / ss_tot
        results[k] = {'r2': r2, 'coeffs': dict(zip(labels, c))}
    return results, labels

# ============================================================
# DRIVER
# ============================================================
print("=" * 65)
print("  LEVEL 2 PROBE - FOUR SCIENCES NETWORK TEST")
print("=" * 65)
print("""
  Question:  Is there a Level 2 (sub-dipole) structure?
  Method:    Compute B_ij(t) per pair, run OD on network library.
  Falsifier: No universal opposing pair across 4 sciences = no
             Level 2 in THIS observable class.
""")

sciences = [
    ("PHYSICS",   simulate_physics),
    ("BIOLOGY",   simulate_biology),
    ("CHEMISTRY", simulate_chemistry),
    ("GEOLOGY",   simulate_geology),
]

all_results = {}
labels_ref = None

for name, sim_fn in sciences:
    print(f"\n{'-' * 50}\n  {name}\n{'-' * 50}")
    t0 = time.time()
    states, dt_sim = sim_fn()
    print(f"  Level 0 simulated:   {states.shape}        [{time.time()-t0:5.1f}s]")

    t0 = time.time()
    H, MI, pairs, dt_info = compute_info_layer(states, dt_sim)
    print(f"  Level 1 extracted:   H{H.shape} MI{MI.shape}  [{time.time()-t0:5.1f}s]")

    t0 = time.time()
    B, t_c = extract_balance_timeseries(H, MI, pairs, dt_info)
    if B is None:
        print(f"  Not enough windows for Level 2; skipping.")
        continue
    print(f"  Dipole B_ij(t):      {B.shape}             [{time.time()-t0:5.1f}s]")
    print(f"  B range: log10 [{np.log10(B.min()+1e-9):.2f}, {np.log10(B.max()+1e-9):.2f}]")

    t0 = time.time()
    res, labels_ref = od_on_dipole_network(B, t_c)
    if res is None:
        print(f"  OD failed; skipping.")
        continue
    print(f"  OD on network:       {len(res)} pairs           [{time.time()-t0:5.1f}s]")

    r2s = [r['r2'] for r in res.values()]
    print(f"  R^2 range: [{min(r2s):.3f}, {max(r2s):.3f}]  mean={np.mean(r2s):.3f}  median={np.median(r2s):.3f}")

    # tally signs within this science
    signs = {label: [] for label in labels_ref}
    for r in res.values():
        if r['r2'] < 0.3:
            continue  # skip low-R^2 fits
        for label, c in r['coeffs'].items():
            if abs(c) > 1e-6:
                signs[label].append(np.sign(c))

    print(f"\n  Operator sign tally (pairs with R^2 >= 0.3 only):")
    for label in labels_ref:
        sl = signs[label]
        if not sl:
            print(f"    {label:<12s} no significant fits")
            continue
        pos = sum(1 for s in sl if s > 0)
        neg = sum(1 for s in sl if s < 0)
        consistency = max(pos, neg) / len(sl) if sl else 0.0
        dominant = '+' if pos > neg else '-'
        print(f"    {label:<12s} +{pos:2d}/-{neg:2d}/n={len(sl):2d}  consistency={consistency:.2f}  dominant={dominant}")

    all_results[name] = {'r2s': r2s, 'signs': signs}

# ============================================================
# UNIVERSAL STRUCTURE SEARCH
# ============================================================
print(f"\n{'=' * 65}\n  UNIVERSAL STRUCTURE SEARCH\n{'=' * 65}")

# classify each operator's dominant sign per science (threshold 0.7)
ops_signs = {}
for label in labels_ref:
    sci_signs = {}
    for sci_name in all_results:
        sl = all_results[sci_name]['signs'][label]
        if not sl:
            sci_signs[sci_name] = 0
            continue
        pos = sum(1 for s in sl if s > 0)
        neg = sum(1 for s in sl if s < 0)
        if len(sl) >= 3 and pos / len(sl) >= 0.7:
            sci_signs[sci_name] = +1
        elif len(sl) >= 3 and neg / len(sl) >= 0.7:
            sci_signs[sci_name] = -1
        else:
            sci_signs[sci_name] = 0
    ops_signs[label] = sci_signs

print("\n  Operator dominant sign per science (threshold 70%):")
print(f"  {'operator':<14s} " + ' '.join(f"{s:>9s}" for s in all_results.keys()) + "  same?")
for label in labels_ref:
    row = ops_signs[label]
    vals = [row[s] for s in all_results.keys()]
    glyphs = ['+' if v > 0 else ('-' if v < 0 else '?') for v in vals]
    all_same = (vals.count(+1) == 4) or (vals.count(-1) == 4)
    print(f"  {label:<14s} " + ' '.join(f"{g:>9s}" for g in glyphs) + f"   {all_same}")

print("\n  Opposing pair search (Level 2 dipole candidate):")
found_any = False
for a, b in combinations(labels_ref, 2):
    sa = [ops_signs[a][s] for s in all_results.keys()]
    sb = [ops_signs[b][s] for s in all_results.keys()]
    opposing = sum(1 for x, y in zip(sa, sb) if x != 0 and y != 0 and x * y < 0)
    if opposing >= 3:
        print(f"    *** {a:<12s} opposes {b:<12s} in {opposing}/4 sciences  <- CANDIDATE")
        found_any = True
if not found_any:
    print("    none found at >=3/4 threshold")

print(f"\n{'=' * 65}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 65}")

print("""
INTERPRETATION GUIDE
====================
- Universal opposing pair (>=3/4)  -> candidate Level 2 dipole.
                                     Validate on NHANES (21 organ
                                     pairs) and PTB-XL next.
- No universal pair, R^2 strong   -> dynamics exist but no dipole
                                     structure; network observables
                                     are real but Level 2 is not
                                     dipole-shaped.
- No universal pair, R^2 weak     -> wrong observable class again;
                                     try transfer entropy, phase lag,
                                     or correlation length next.
""")
