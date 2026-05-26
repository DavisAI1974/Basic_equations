"""
STATIC ALGEBRAIC DIPOLE TEST
=============================
Reframe: a governing equation can be DIFFERENTIAL (dH/dt = f(...))
or ALGEBRAIC (f(H, ...) = const, holds at every instant).

The Level 1 finding was differential.  DNA gave R^2 = 0 on the
differential form BUT preserved the sign opposition of H_a^2 vs
H_a*H_b — suggesting the same dipole exists as an ALGEBRAIC
constraint, not a flow law.

Test: For each system, scatter H_a^2(t) vs H_a*H_b(t).
Fit best line and best low-order polynomial.  Measure R^2.

  R^2 >= 0.7  -> strong algebraic constraint (static dipole)
  R^2 0.3-0.7 -> partial constraint
  R^2 < 0.3   -> no algebraic constraint, dipole is flow-only

Run on:
  PHYSICS, BIOLOGY, CHEMISTRY, GEOLOGY (synthetic, 2-system)
  DNA-RANDOM, DNA-CODON-BIASED  (synthetic, sliding window)
"""
import numpy as np
from numpy.linalg import lstsq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
T_START = time.time()

# ============================================================
# INFO PRIMITIVES
# ============================================================
def entropy_1d(samples, bins=30):
    h, e = np.histogram(samples, bins=bins, density=True)
    h = h[h > 0]; dx = e[1] - e[0]
    return -np.sum(h * np.log(h) * dx)

def mi_2d(x, y, bins=20):
    hx = entropy_1d(x, bins)
    hy = entropy_1d(y, bins)
    h2, ex, ey = np.histogram2d(x, y, bins=bins, density=True)
    m = h2 > 0
    if not m.any():
        return 0.0
    dA = (ex[1]-ex[0])*(ey[1]-ey[0])
    hxy = -np.sum(h2[m] * np.log(h2[m])) * dA
    return max(0.0, hx + hy - hxy)

# ============================================================
# 4 SCIENCES — 2-SYSTEM TIME SERIES (gives H_a, H_b vs t)
# ============================================================
def physics_2sys(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    x1 = 0.5*np.random.randn(N_ens); v1 = 0.1*np.random.randn(N_ens)
    x2 = 0.5*np.random.randn(N_ens); v2 = 0.1*np.random.randn(N_ens)
    K = 0.2
    H_a = np.zeros(steps); H_b = np.zeros(steps)
    for t in range(steps):
        F1 = -x1 - 0.5*x1**3 + K*(x2-x1) + 0.1*np.sin(0.5*t*dt)
        F2 = -1.2*x2 - 0.4*x2**3 + K*(x1-x2) + 0.08*np.cos(0.4*t*dt)
        v1 += F1*dt - 0.05*v1*dt; v2 += F2*dt - 0.05*v2*dt
        x1 += v1*dt + 0.02*np.random.randn(N_ens)*np.sqrt(dt)
        x2 += v2*dt + 0.02*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(x1); H_b[t] = entropy_1d(x2)
    return H_a, H_b

def biology_2sys(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    prey = 1.0 + 0.1*np.random.randn(N_ens)
    pred = 0.5 + 0.1*np.random.randn(N_ens)
    H_a = np.zeros(steps); H_b = np.zeros(steps)
    for t in range(steps):
        dprey = (1.0*prey - 0.5*prey*pred)*dt
        dpred = (0.3*prey*pred - 0.4*pred)*dt
        prey = np.clip(prey + dprey + 0.01*np.random.randn(N_ens)*np.sqrt(dt), 0.01, None)
        pred = np.clip(pred + dpred + 0.01*np.random.randn(N_ens)*np.sqrt(dt), 0.01, None)
        H_a[t] = entropy_1d(prey); H_b[t] = entropy_1d(pred)
    return H_a, H_b

def chemistry_2sys(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    X = 1.0 + 0.1*np.random.randn(N_ens)
    Y = 3.0 + 0.1*np.random.randn(N_ens)
    A_, B_ = 1.0, 3.0
    H_a = np.zeros(steps); H_b = np.zeros(steps)
    for t in range(steps):
        dX = (A_ - (B_+1)*X + X**2*Y)*dt
        dY = (B_*X - X**2*Y)*dt
        X += dX + 0.01*np.random.randn(N_ens)*np.sqrt(dt)
        Y += dY + 0.01*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(X); H_b[t] = entropy_1d(Y)
    return H_a, H_b

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

# ============================================================
# DNA — SLIDING WINDOW H_R, H_S
# ============================================================
BASES = np.array(['A', 'C', 'G', 'T'])

def gen_random_dna(n):
    return np.random.choice(BASES, size=n)

def gen_codon_biased(n):
    n = (n // 3) * 3
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
    idx = np.random.choice(len(codons), size=n//3, p=weights)
    seq = ''.join(codons[i] for i in idx)
    return np.array(list(seq))

def binary_entropy(p):
    if p <= 0 or p >= 1: return 0.0
    return -p*np.log(p) - (1-p)*np.log(1-p)

def dna_sliding(seq, window=120, stride=4):
    is_R = np.isin(seq, ['A','G']).astype(np.float64)
    is_S = np.isin(seq, ['G','C']).astype(np.float64)
    starts = np.arange(0, len(seq)-window+1, stride)
    H_R = np.zeros(len(starts)); H_S = np.zeros(len(starts))
    for k, s in enumerate(starts):
        e = s + window
        H_R[k] = binary_entropy(is_R[s:e].mean())
        H_S[k] = binary_entropy(is_S[s:e].mean())
    return H_R, H_S

# ============================================================
# STATIC FIT — H_a^2 vs H_a*H_b
# ============================================================
def static_fit(H_a, H_b):
    """
    Test the algebraic dipole: is H_a^2 a function of H_a*H_b?
    Fit polynomial of degree 1, 2, 3.  Report R^2 each.
    """
    x = H_a * H_b
    y = H_a ** 2
    pearson = np.corrcoef(x, y)[0, 1]
    fits = {}
    for deg in (1, 2, 3):
        lib = np.column_stack([x**k for k in range(deg+1)])
        c, *_ = lstsq(lib, y, rcond=None)
        pred = lib @ c
        ss_res = np.sum((y - pred)**2)
        ss_tot = np.sum((y - y.mean())**2) + 1e-15
        fits[deg] = 1 - ss_res/ss_tot
    return {
        'pearson': pearson,
        'r2_linear': fits[1],
        'r2_quad':   fits[2],
        'r2_cubic':  fits[3],
        'x': x, 'y': y,
    }

# ============================================================
# RUN ALL SYSTEMS
# ============================================================
print("=" * 70)
print("  STATIC ALGEBRAIC DIPOLE TEST")
print("=" * 70)
print("""
  Hypothesis:  H_a^2 = f(H_a*H_b)  holds algebraically at every
               moment, independent of any flow.
  Falsifier:   Low R^2 across systems -> no algebraic constraint,
               dipole is flow-only.
""")

systems = []
for name, fn in [
    ("PHYSICS",   physics_2sys),
    ("BIOLOGY",   biology_2sys),
    ("CHEMISTRY", chemistry_2sys),
    ("GEOLOGY",   geology_2sys),
]:
    t0 = time.time()
    H_a, H_b = fn()
    print(f"  {name:<10s} simulated  [{time.time()-t0:.1f}s]")
    systems.append((name, H_a, H_b))

# DNA
print(f"  DNA generating...")
seq_rand  = gen_random_dna(60000)
seq_codon = gen_codon_biased(60000)
H_R_r, H_S_r = dna_sliding(seq_rand)
H_R_c, H_S_c = dna_sliding(seq_codon)
systems.append(("DNA-RANDOM",      H_R_r, H_S_r))
systems.append(("DNA-CODON",       H_R_c, H_S_c))

print(f"\n  {'system':<14s} {'Pearson':>9s} {'R^2 lin':>9s} {'R^2 quad':>9s} {'R^2 cubic':>10s} {'verdict':>12s}")
print(f"  {'-'*14} {'-'*9} {'-'*9} {'-'*9} {'-'*10} {'-'*12}")

results = {}
for name, H_a, H_b in systems:
    r = static_fit(H_a, H_b)
    results[name] = r
    best = max(r['r2_linear'], r['r2_quad'], r['r2_cubic'])
    if best >= 0.7:
        verdict = "ALGEBRAIC"
    elif best >= 0.3:
        verdict = "partial"
    else:
        verdict = "flow-only"
    print(f"  {name:<14s} {r['pearson']:>+9.3f} {r['r2_linear']:>9.3f} "
          f"{r['r2_quad']:>9.3f} {r['r2_cubic']:>10.3f} {verdict:>12s}")

# ============================================================
# SCATTER FIGURE
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()
plot_order = ["PHYSICS", "BIOLOGY", "CHEMISTRY", "GEOLOGY",
              "DNA-RANDOM", "DNA-CODON"]
for ax, name in zip(axes, plot_order):
    if name not in results:
        ax.axis('off'); continue
    r = results[name]
    ax.scatter(r['x'], r['y'], s=4, alpha=0.35, color='#2563eb')
    # best-fit cubic for visual
    xs = np.linspace(r['x'].min(), r['x'].max(), 100)
    lib = np.column_stack([xs**k for k in range(4)])
    lib_data = np.column_stack([r['x']**k for k in range(4)])
    c, *_ = lstsq(lib_data, r['y'], rcond=None)
    ax.plot(xs, lib @ c, 'r-', linewidth=1.5, alpha=0.8)
    best = max(r['r2_linear'], r['r2_quad'], r['r2_cubic'])
    ax.set_title(f"{name}\nR^2(cubic)={r['r2_cubic']:.3f}  Pearson={r['pearson']:+.2f}",
                 fontsize=10)
    ax.set_xlabel("H_a * H_b"); ax.set_ylabel("H_a^2")
    ax.grid(alpha=0.3)

plt.suptitle("Static Algebraic Dipole Test: H_a^2 vs H_a*H_b",
             fontsize=13, y=1.00)
plt.tight_layout()
out = '/mnt/user-data/outputs/static_dipole_scatter.png'
plt.savefig(out, dpi=120, bbox_inches='tight')
print(f"\n  Scatter figure: {out}")

print(f"\n{'=' * 70}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 70}")
print("""
INTERPRETATION GUIDE
====================
- All 4 sciences ALGEBRAIC -> dipole has BOTH forms; differential
  was what we found first because we asked for flow, but the
  underlying constraint is static.  This unifies dynamic and
  static substrates (DNA) under one architecture.
- DNA ALGEBRAIC + sciences flow-only -> dipole is purely static,
  the time-derivative version we found is a SHADOW of an algebraic
  invariant we hadn't isolated.
- Sciences algebraic + DNA flow-only -> reversed expectation,
  redesign needed.
- All weak -> the static reformulation is also the wrong framing;
  the dipole is something else entirely.
""")
