"""
TRIPLE FOLLOW-UP
================
Three tests after the static algebraic dipole result.

  A. DNA codon-position-stratified.  H computed at codon
     position 1, 2, 3 separately.  These ARE independent
     observables, unlike H_R and H_S which shared substrate.
     Test the static algebraic dipole on properly independent
     partitions.

  B. Extract the geology equation.  R^2 = 0.997 from the
     prior run -> effectively a discovered law.  Write it
     in closed form.

  C. Level 2 static test.  Rerun the 6-node network, compute
     B_ij(t) per pair, then test whether log(B_k)^2 vs
     log(B_k) * log(B_net) holds algebraically.  This is
     the static reformulation of the Level 2 question we
     never actually ran.
"""
import numpy as np
from numpy.linalg import lstsq
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
T_START = time.time()
BASES = np.array(['A', 'C', 'G', 'T'])

# ============================================================
# SHARED PRIMITIVES
# ============================================================
def entropy_1d(samples, bins=30):
    h, e = np.histogram(samples, bins=bins, density=True)
    h = h[h > 0]; dx = e[1] - e[0]
    return -np.sum(h * np.log(h) * dx)

def binary_entropy(p):
    if p <= 0 or p >= 1: return 0.0
    return -p*np.log(p) - (1-p)*np.log(1-p)

def static_fit(x, y):
    """Fit y = polynomial(x) for degrees 1-3. Return R^2s and coefficients."""
    fits = {}
    coeffs = {}
    for deg in (1, 2, 3):
        lib = np.column_stack([x**k for k in range(deg+1)])
        c, *_ = lstsq(lib, y, rcond=None)
        pred = lib @ c
        ss_res = np.sum((y - pred)**2)
        ss_tot = np.sum((y - y.mean())**2) + 1e-15
        fits[deg] = 1 - ss_res/ss_tot
        coeffs[deg] = c
    pearson = np.corrcoef(x, y)[0, 1]
    return {'pearson': pearson, 'r2': fits, 'coeffs': coeffs}

# ============================================================
# A. DNA CODON-POSITION-STRATIFIED
# ============================================================
def gen_codon_biased_dna(n_codons):
    codons, weights = [], []
    for b1 in BASES:
        for b2 in BASES:
            for b3 in BASES:
                c = b1+b2+b3
                if c in ('TAA','TAG','TGA'):
                    w = 0.001
                else:
                    w = 1.0
                    if b3 in ('G','C'): w *= 1.6
                    if b1 in ('A','G'): w *= 1.15
                    if b2 == 'A': w *= 1.2
                codons.append(c); weights.append(w)
    weights = np.array(weights) / np.sum(weights)
    idx = np.random.choice(len(codons), size=n_codons, p=weights)
    seq = ''.join(codons[i] for i in idx)
    return np.array(list(seq))

def gen_random_dna(n_codons):
    return np.random.choice(BASES, size=n_codons*3)

def codon_position_streams(seq):
    """Return arrays of bases at codon position 1, 2, 3."""
    n = (len(seq) // 3) * 3
    s = seq[:n]
    p1 = s[0::3]
    p2 = s[1::3]
    p3 = s[2::3]
    return p1, p2, p3

def sliding_codon_entropies(p1, p2, p3, window_codons=80, stride=2):
    """
    At each window (in codon units), compute binary R/Y entropy
    at each codon position INDEPENDENTLY.
    Returns H1(t), H2(t), H3(t) - now genuinely independent observables
    because each is computed from disjoint nucleotides.
    """
    n = len(p1)
    starts = np.arange(0, n - window_codons + 1, stride)
    H1 = np.zeros(len(starts)); H2 = np.zeros(len(starts)); H3 = np.zeros(len(starts))
    for k, s in enumerate(starts):
        e = s + window_codons
        H1[k] = binary_entropy(np.isin(p1[s:e], ['A','G']).mean())
        H2[k] = binary_entropy(np.isin(p2[s:e], ['A','G']).mean())
        H3[k] = binary_entropy(np.isin(p3[s:e], ['A','G']).mean())
    return H1, H2, H3

print("=" * 70)
print("  A. DNA CODON-POSITION-STRATIFIED STATIC TEST")
print("=" * 70)
print("""
  H computed at codon positions 1, 2, 3 separately.  These
  partitions use DISJOINT nucleotides, so any algebraic
  relationship between them is not a substrate-sharing artifact.
""")

for name, seq in [
    ("DNA-RANDOM",  gen_random_dna(20000)),
    ("DNA-CODON",   gen_codon_biased_dna(20000)),
]:
    p1, p2, p3 = codon_position_streams(seq)
    H1, H2, H3 = sliding_codon_entropies(p1, p2, p3)
    print(f"\n  {name}")
    print(f"  {'-' * len(name)}")
    print(f"    H1: mean={H1.mean():.3f} std={H1.std():.4f}")
    print(f"    H2: mean={H2.mean():.3f} std={H2.std():.4f}")
    print(f"    H3: mean={H3.mean():.3f} std={H3.std():.4f}")

    # Test H_i^2 vs H_i*H_j for all 6 ordered pairs
    print(f"    Static dipole tests (H_a^2 vs H_a*H_b):")
    print(f"      {'pair':<12s} {'Pearson':>9s} {'R2 lin':>8s} {'R2 quad':>9s} {'R2 cubic':>10s}")
    for a_idx, b_idx in [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]:
        Ha = [H1,H2,H3][a_idx]; Hb = [H1,H2,H3][b_idx]
        r = static_fit(Ha*Hb, Ha**2)
        tag = f"H{a_idx+1}^2 ~ H{a_idx+1}*H{b_idx+1}"
        print(f"      {tag:<14s} {r['pearson']:>+9.3f} {r['r2'][1]:>8.3f} {r['r2'][2]:>9.3f} {r['r2'][3]:>10.3f}")

# ============================================================
# B. GEOLOGY EQUATION IN CLOSED FORM
# ============================================================
print(f"\n{'=' * 70}")
print(f"  B. GEOLOGY DIPOLE EQUATION (closed form)")
print(f"{'=' * 70}")

# Regenerate geology 2-system data
def geology_2sys(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    u1 = 0.1*np.random.randn(N_ens); v1 = 0.05*np.random.randn(N_ens)
    u2 = 0.1*np.random.randn(N_ens); v2 = 0.05*np.random.randn(N_ens)
    K = 0.5
    H_a = np.zeros(steps); H_b = np.zeros(steps)
    for t in range(steps):
        fric1 = -0.6*v1/(1.0 + 8.0*v1**2)
        fric2 = -0.6*v2/(1.0 + 8.0*v2**2)
        F1 = K*(u2-u1) - 0.6*(u1 - 0.04*t*dt) + fric1
        F2 = K*(u1-u2) - 0.6*(u2 - 0.04*t*dt) + fric2
        v1 += F1*dt; v2 += F2*dt
        u1 += v1*dt + 0.015*np.random.randn(N_ens)*np.sqrt(dt)
        u2 += v2*dt + 0.015*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(u1); H_b[t] = entropy_1d(u2)
    return H_a, H_b

H_a, H_b = geology_2sys()
X = H_a * H_b
Y = H_a ** 2

# Linear fit Y = m*X + b
lin = static_fit(X, Y)
m, b_lin = lin['coeffs'][1][1], lin['coeffs'][1][0]
print(f"\n  Linear fit:    H_a^2 = {m:+.6f} * (H_a*H_b) + {b_lin:+.6f}")
print(f"                 R^2 = {lin['r2'][1]:.6f}")

# Try alternative forms
# Power law: Y = a * X^b  -> log Y = log a + b * log X
lX = np.log(X + 1e-12); lY = np.log(Y + 1e-12)
lib = np.column_stack([lX, np.ones(len(lX))])
c, *_ = lstsq(lib, lY, rcond=None)
b_pow, log_a = c[0], c[1]
pred = lib @ c
r2_pow = 1 - np.sum((lY - pred)**2) / (np.sum((lY - lY.mean())**2) + 1e-15)
print(f"\n  Power law:     H_a^2 = {np.exp(log_a):.6f} * (H_a*H_b)^{b_pow:.6f}")
print(f"                 R^2 (in log space) = {r2_pow:.6f}")

# Identity ratio: is H_a^2 / (H_a*H_b) = H_a/H_b constant?
ratio = H_a / H_b
print(f"\n  Note: H_a^2 / (H_a*H_b) = H_a/H_b")
print(f"    H_a/H_b: mean={ratio.mean():.4f} std={ratio.std():.4f} CV={ratio.std()/abs(ratio.mean()):.4f}")
print(f"    If H_a/H_b is nearly constant, the dipole equation reduces to")
print(f"    H_a = k * H_b  with k = {ratio.mean():.4f}")

# Slope = H_a / H_b should match the linear-fit slope
print(f"\n  Cross-check: linear fit slope m = {m:.6f}, mean(H_a/H_b) = {ratio.mean():.6f}")

# ============================================================
# C. LEVEL 2 STATIC ALGEBRAIC TEST
# ============================================================
print(f"\n{'=' * 70}")
print(f"  C. LEVEL 2 STATIC ALGEBRAIC TEST")
print(f"{'=' * 70}")
print(f"""
  Same reframe applied one layer up.
  For each pair k in a network, scatter log(B_k)^2 vs
  log(B_k) * log(B_net).  If an algebraic constraint holds
  across the dipole network, R^2 should be high.
""")

def physics_network(N=6, N_ens=500, T=30.0, dt=0.02):
    steps = int(T/dt)
    x = 0.5*np.random.randn(N_ens, N); v = 0.1*np.random.randn(N_ens, N)
    K = 0.15
    out = np.zeros((steps, N_ens, N))
    for t in range(steps):
        nbr = np.zeros_like(x)
        nbr[:,1:-1] = x[:,:-2] + x[:,2:] - 2*x[:,1:-1]
        nbr[:,0] = x[:,1] - x[:,0]; nbr[:,-1] = x[:,-2] - x[:,-1]
        F = -x - 0.5*x**3 + 0.1*np.sin(0.5*t*dt) + K*nbr
        v = v + F*dt - 0.05*v*dt
        x = x + v*dt + 0.02*np.random.randn(*x.shape)*np.sqrt(dt)
        out[t] = x
    return out, dt

def mi_2d(x, y, bins=18):
    hx = entropy_1d(x, bins); hy = entropy_1d(y, bins)
    h2, ex, ey = np.histogram2d(x, y, bins=bins, density=True)
    m = h2 > 0
    if not m.any(): return 0.0
    dA = (ex[1]-ex[0])*(ey[1]-ey[0])
    hxy = -np.sum(h2[m] * np.log(h2[m])) * dA
    return max(0.0, hx + hy - hxy)

def info_layer(states, dt, bins=18, subsample=2):
    T, N_ens, N = states.shape
    times = np.arange(0, T, subsample)
    pairs = list(combinations(range(N), 2))
    H = np.zeros((len(times), N))
    MI = np.zeros((len(times), len(pairs)))
    for idx, t in enumerate(times):
        for i in range(N):
            H[idx, i] = entropy_1d(states[t, :, i], bins)
        for k, (i, j) in enumerate(pairs):
            MI[idx, k] = mi_2d(states[t, :, i], states[t, :, j], bins)
    return H, MI, pairs, dt*subsample

def balance_series(H, MI, pairs, dt_info, window=30, stride=4):
    T, N = H.shape
    n_pairs = len(pairs)
    n_w = (T - window) // stride
    B = np.zeros((n_w, n_pairs))
    for w in range(n_w):
        s = w*stride; e = s + window
        for k, (i, j) in enumerate(pairs):
            Hi, Hj, Mij = H[s:e,i], H[s:e,j], MI[s:e,k]
            dHi = np.gradient(Hi, dt_info)
            lib = np.column_stack([Hi, Hj, Hi**2, Hj**2, Hi*Hj, Mij, np.ones(len(Hi))])
            c, *_ = lstsq(lib, dHi, rcond=None)
            num = abs(c[2]); den = abs(c[4]) + 1e-10
            B[w, k] = num / den
    return B

states, dt_sim = physics_network()
H, MI, pairs, dt_info = info_layer(states, dt_sim)
B = balance_series(H, MI, pairs, dt_info)
print(f"  Network B shape: {B.shape}  (n_windows x n_pairs)")
print(f"  log10(B) range: [{np.log10(B.min()+1e-12):.2f}, {np.log10(B.max()+1e-12):.2f}]")

# log-transform for sanity
Bl = np.log(B + 1e-6)
n_w, n_pairs = Bl.shape

print(f"\n  Static algebraic test per pair (log(B_k)^2 vs log(B_k)*log(B_net)):")
print(f"  {'pair':<6s} {'Pearson':>9s} {'R2 lin':>8s} {'R2 quad':>9s} {'R2 cubic':>10s}")
results_l2 = []
for k in range(n_pairs):
    bk = Bl[:, k]
    mask = np.ones(n_pairs, bool); mask[k] = False
    b_net = Bl[:, mask].mean(axis=1)
    X = bk * b_net
    Y = bk ** 2
    r = static_fit(X, Y)
    results_l2.append(r)
    print(f"  {k:<6d} {r['pearson']:>+9.3f} {r['r2'][1]:>8.3f} {r['r2'][2]:>9.3f} {r['r2'][3]:>10.3f}")

# Aggregate
pearsons = [r['pearson'] for r in results_l2]
r2_lins = [r['r2'][1] for r in results_l2]
r2_cubs = [r['r2'][3] for r in results_l2]
print(f"\n  Aggregate over {n_pairs} pairs:")
print(f"    |Pearson|:  mean={np.mean(np.abs(pearsons)):.3f}  median={np.median(np.abs(pearsons)):.3f}")
print(f"    R^2 linear: mean={np.mean(r2_lins):.3f}  median={np.median(r2_lins):.3f}  max={np.max(r2_lins):.3f}")
print(f"    R^2 cubic:  mean={np.mean(r2_cubs):.3f}  median={np.median(r2_cubs):.3f}  max={np.max(r2_cubs):.3f}")

# Also test the AGGREGATE relationship across all pairs and times pooled
Xall = []; Yall = []
for k in range(n_pairs):
    mask = np.ones(n_pairs, bool); mask[k] = False
    b_net = Bl[:, mask].mean(axis=1)
    Xall.append(Bl[:, k] * b_net)
    Yall.append(Bl[:, k] ** 2)
Xall = np.concatenate(Xall); Yall = np.concatenate(Yall)
r_pool = static_fit(Xall, Yall)
print(f"\n  Pooled across all pairs (N={len(Xall)} points):")
print(f"    Pearson = {r_pool['pearson']:+.3f}")
print(f"    R^2 linear = {r_pool['r2'][1]:.3f}")
print(f"    R^2 quadratic = {r_pool['r2'][2]:.3f}")
print(f"    R^2 cubic = {r_pool['r2'][3]:.3f}")

# Figure: Level 2 scatter
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.scatter(Xall, Yall, s=3, alpha=0.25, color='#7c3aed')
xs = np.linspace(Xall.min(), Xall.max(), 100)
lib_data = np.column_stack([Xall**k for k in range(4)])
c_p, *_ = lstsq(lib_data, Yall, rcond=None)
lib_xs = np.column_stack([xs**k for k in range(4)])
ax.plot(xs, lib_xs @ c_p, 'r-', linewidth=1.5)
ax.set_xlabel("log(B_k) * log(B_net)")
ax.set_ylabel("log(B_k)^2")
ax.set_title(f"Level 2 static test (pooled)\nR^2 cubic = {r_pool['r2'][3]:.3f}  Pearson = {r_pool['pearson']:+.2f}")
ax.grid(alpha=0.3)

# Geology overlay for comparison
ax = axes[1]
ax.scatter(H_a*H_b, H_a**2, s=4, alpha=0.4, color='#059669')
xs = np.linspace((H_a*H_b).min(), (H_a*H_b).max(), 100)
lib_data = np.column_stack([(H_a*H_b)**k for k in range(2)])
c_p2, *_ = lstsq(lib_data, H_a**2, rcond=None)
ax.plot(xs, c_p2[0] + c_p2[1]*xs, 'r-', linewidth=1.5)
ax.set_xlabel("H_a * H_b")
ax.set_ylabel("H_a^2")
ax.set_title(f"Geology Level 1 reference\nR^2 linear = {lin['r2'][1]:.4f}  slope = {m:.3f}")
ax.grid(alpha=0.3)

plt.tight_layout()
out = '/mnt/user-data/outputs/triple_followup.png'
plt.savefig(out, dpi=120, bbox_inches='tight')
print(f"\n  Figure: {out}")

print(f"\n{'=' * 70}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 70}")
