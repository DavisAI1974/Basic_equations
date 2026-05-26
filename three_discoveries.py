"""
THREE DISCOVERIES FOLLOW-UP
============================
After the nonlinear deep dive identified:
  - Chemistry quadratic in H_a*H_b      (R^2 = 0.943)
  - Geology rank-deficient (2 constraints)
  - Biology dipole lives in MI-space

This script:
  1. Extracts chemistry's exact quadratic equation
     and geology's SECOND (smaller) constraint direction.
  2. Multi-variable SVD on Level 2 dipole network — does the
     operator cloud across dipoles drop in rank?  This is the
     algebraic Level 2 search we have not actually run.
  3. Biology re-run with MI as the primary observable.  Test
     algebraic constraints in MI-space across the 4 sciences.
"""
import numpy as np
from numpy.linalg import lstsq, svd
from itertools import combinations
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
T_START = time.time()

# ============================================================
# PRIMITIVES
# ============================================================
def entropy_1d(samples, bins=30):
    h, e = np.histogram(samples, bins=bins, density=True)
    h = h[h > 0]; dx = e[1] - e[0]
    return -np.sum(h * np.log(h) * dx)

def mi_2d(x, y, bins=20):
    hx = entropy_1d(x, bins); hy = entropy_1d(y, bins)
    h2, ex, ey = np.histogram2d(x, y, bins=bins, density=True)
    m = h2 > 0
    if not m.any(): return 0.0
    dA = (ex[1]-ex[0])*(ey[1]-ey[0])
    hxy = -np.sum(h2[m] * np.log(h2[m])) * dA
    return max(0.0, hx + hy - hxy)

def r2(y, pred):
    ss_res = np.sum((y - pred)**2)
    ss_tot = np.sum((y - y.mean())**2) + 1e-15
    return 1 - ss_res/ss_tot

# ============================================================
# FOUR SCIENCES
# ============================================================
def physics(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    x1 = 0.5*np.random.randn(N_ens); v1 = 0.1*np.random.randn(N_ens)
    x2 = 0.5*np.random.randn(N_ens); v2 = 0.1*np.random.randn(N_ens)
    H_a = np.zeros(steps); H_b = np.zeros(steps); MI = np.zeros(steps)
    for t in range(steps):
        F1 = -x1 - 0.5*x1**3 + 0.2*(x2-x1) + 0.1*np.sin(0.5*t*dt)
        F2 = -1.2*x2 - 0.4*x2**3 + 0.2*(x1-x2) + 0.08*np.cos(0.4*t*dt)
        v1 += F1*dt - 0.05*v1*dt; v2 += F2*dt - 0.05*v2*dt
        x1 += v1*dt + 0.02*np.random.randn(N_ens)*np.sqrt(dt)
        x2 += v2*dt + 0.02*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(x1); H_b[t] = entropy_1d(x2)
        MI[t] = mi_2d(x1, x2)
    return H_a, H_b, MI

def biology(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    prey = 1.0 + 0.1*np.random.randn(N_ens)
    pred = 0.5 + 0.1*np.random.randn(N_ens)
    H_a = np.zeros(steps); H_b = np.zeros(steps); MI = np.zeros(steps)
    for t in range(steps):
        dprey = (1.0*prey - 0.5*prey*pred)*dt
        dpred = (0.3*prey*pred - 0.4*pred)*dt
        prey = np.clip(prey + dprey + 0.01*np.random.randn(N_ens)*np.sqrt(dt), 0.01, None)
        pred = np.clip(pred + dpred + 0.01*np.random.randn(N_ens)*np.sqrt(dt), 0.01, None)
        H_a[t] = entropy_1d(prey); H_b[t] = entropy_1d(pred)
        MI[t] = mi_2d(prey, pred)
    return H_a, H_b, MI

def chemistry(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    X = 1.0 + 0.1*np.random.randn(N_ens)
    Y = 3.0 + 0.1*np.random.randn(N_ens)
    A_, B_ = 1.0, 3.0
    H_a = np.zeros(steps); H_b = np.zeros(steps); MI = np.zeros(steps)
    for t in range(steps):
        dX = (A_ - (B_+1)*X + X**2*Y)*dt
        dY = (B_*X - X**2*Y)*dt
        X += dX + 0.01*np.random.randn(N_ens)*np.sqrt(dt)
        Y += dY + 0.01*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(X); H_b[t] = entropy_1d(Y)
        MI[t] = mi_2d(X, Y)
    return H_a, H_b, MI

def geology(N_ens=600, T=30.0, dt=0.02):
    steps = int(T/dt)
    u1 = 0.1*np.random.randn(N_ens); v1 = 0.05*np.random.randn(N_ens)
    u2 = 0.1*np.random.randn(N_ens); v2 = 0.05*np.random.randn(N_ens)
    H_a = np.zeros(steps); H_b = np.zeros(steps); MI = np.zeros(steps)
    for t in range(steps):
        fric1 = -0.6*v1/(1.0 + 8.0*v1**2)
        fric2 = -0.9*v2/(1.0 + 14.0*v2**2)
        F1 = 0.5*(u2-u1) - 0.6*(u1 - 0.04*t*dt) + fric1
        F2 = 0.5*(u1-u2) - 0.45*(u2 - 0.03*t*dt) + fric2
        v1 += F1*dt; v2 += F2*dt
        u1 += v1*dt + 0.015*np.random.randn(N_ens)*np.sqrt(dt)
        u2 += v2*dt + 0.020*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(u1); H_b[t] = entropy_1d(u2)
        MI[t] = mi_2d(u1, u2)
    return H_a, H_b, MI

print("=" * 70)
print("  THREE DISCOVERIES FOLLOW-UP")
print("=" * 70)

systems_data = {}
print("\n  Simulating 4 sciences...")
for name, fn in [('PHYSICS', physics), ('BIOLOGY', biology),
                 ('CHEMISTRY', chemistry), ('GEOLOGY', geology)]:
    t0 = time.time()
    H_a, H_b, MI = fn()
    print(f"    {name:<10s} done  [{time.time()-t0:.1f}s]")
    systems_data[name] = (H_a, H_b, MI)

# ============================================================
# DISCOVERY 1: EXACT EQUATIONS
# ============================================================
print(f"\n{'=' * 70}")
print("  DISCOVERY 1: EXACT ALGEBRAIC EQUATIONS")
print(f"{'=' * 70}")

# --- Chemistry quadratic ---
H_a, H_b, MI = systems_data['CHEMISTRY']
x = H_a * H_b
y = H_a ** 2
lib = np.column_stack([np.ones(len(x)), x, x**2])
c, *_ = lstsq(lib, y, rcond=None)
pred = lib @ c
r2_chem = r2(y, pred)
print(f"\n  CHEMISTRY:")
print(f"  H_a^2 = {c[0]:+.4f} {c[1]:+.4f}*(H_a*H_b) {c[2]:+.4f}*(H_a*H_b)^2")
print(f"  R^2 = {r2_chem:.4f}")
# range info
print(f"  Valid range of H_a*H_b: [{x.min():.3f}, {x.max():.3f}]")
print(f"  Range of H_a^2:        [{y.min():.3f}, {y.max():.3f}]")

# --- Geology: full SVD, all constraints below threshold ---
H_a, H_b, MI = systems_data['GEOLOGY']
ops = np.column_stack([H_a, H_b, H_a**2, H_b**2, H_a*H_b, MI])
op_names = ['H_a', 'H_b', 'H_a^2', 'H_b^2', 'H_a*H_b', 'MI']
ops_c = ops - ops.mean(axis=0)
norms = np.linalg.norm(ops_c, axis=0) + 1e-15
ops_n = ops_c / norms
U, S, Vt = svd(ops_n, full_matrices=False)
s_norm = S / S.max()

print(f"\n  GEOLOGY constraints (singular values):")
for i, s in enumerate(s_norm):
    flag = " <-- constraint" if s < 0.05 else ""
    print(f"    sigma_{i+1} = {s:.6f}{flag}")

# extract each constraint direction below threshold
print(f"\n  GEOLOGY constraint directions:")
for i in range(len(s_norm)):
    if s_norm[i] < 0.05:
        v = Vt[i]  # constraint direction
        # also rescale by 1/norms so the constraint is in original units
        v_orig = v / norms
        v_orig = v_orig / np.linalg.norm(v_orig)  # renormalize
        terms = sorted(zip(op_names, v_orig), key=lambda x: -abs(x[1]))
        eq = " + ".join(f"({w:+.3f})*{n}" for n, w in terms if abs(w) > 0.05)
        # also compute residual of constraint
        constraint_value = ops @ v_orig
        resid_std = constraint_value.std()
        mean_op = ops.mean()
        print(f"    sigma={s_norm[i]:.6f}: {eq}")
        print(f"      constraint std = {resid_std:.6f}  (relative to mean-op {mean_op:.3f})")

# ============================================================
# DISCOVERY 2: LEVEL 2 MULTI-VARIABLE SVD
# ============================================================
print(f"\n{'=' * 70}")
print("  DISCOVERY 2: LEVEL 2 MULTI-VARIABLE SVD")
print(f"{'=' * 70}")
print("""
  Build dipole network across 6 nodes, compute B_ij(t) per pair.
  Then run SVD on the operator cloud over the NETWORK of dipoles.
  If the cloud is rank-deficient, an algebraic Level 2 constraint
  exists.
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

def info_layer_network(states, dt, bins=18, subsample=2):
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

t0 = time.time()
states, dt_sim = physics_network()
H_net, MI_net, pairs, dt_info = info_layer_network(states, dt_sim)
B = balance_series(H_net, MI_net, pairs, dt_info)
Bl = np.log(B + 1e-6)
n_w, n_pairs = Bl.shape
print(f"  Network B: {Bl.shape}  ({time.time()-t0:.1f}s)")

# Build network operator cloud for each time window
# For each window, we have n_pairs values of B.
# Operators across the network: B_k, B_k^2, mean(B), var(B), mean(B^2), max(B), min(B), etc.
print(f"\n  Operator cloud across (windows x pairs) = ({n_w}, {n_pairs})")

# Build a per-pair operator matrix similar to single-system case
# Then SVD it across all windows pooled
all_rows = []
for w in range(n_w):
    bw = Bl[w]   # n_pairs values
    for k in range(n_pairs):
        bk = bw[k]
        mask = np.ones(n_pairs, bool); mask[k] = False
        b_net = bw[mask].mean()
        b_net2 = (bw[mask]**2).mean()
        # row: per-pair operators in their network context
        all_rows.append([bk, b_net, bk**2, b_net**2, bk*b_net, b_net2])
op_mat = np.array(all_rows)
op_names_l2 = ['B_k', 'B_net', 'B_k^2', 'B_net^2', 'B_k*B_net', 'mean(B^2)_net']
print(f"  Operator matrix shape: {op_mat.shape}")
print(f"  Range per operator:")
for nm, col in zip(op_names_l2, op_mat.T):
    print(f"    {nm:<14s}  [{col.min():+7.2f}, {col.max():+7.2f}]")

# normalize and SVD
op_c = op_mat - op_mat.mean(axis=0)
norms_l2 = np.linalg.norm(op_c, axis=0) + 1e-15
op_n = op_c / norms_l2
_, S_l2, Vt_l2 = svd(op_n, full_matrices=False)
s_l2 = S_l2 / S_l2.max()
rank_l2 = (s_l2 > 0.05).sum()
print(f"\n  Level 2 SVD singular values (normalized):")
for i, s in enumerate(s_l2):
    flag = " <-- constraint" if s < 0.05 else ""
    print(f"    sigma_{i+1} = {s:.6f}{flag}")
print(f"  Effective rank = {rank_l2} / {len(s_l2)}")

if rank_l2 < len(s_l2):
    print(f"\n  Level 2 constraint directions:")
    for i in range(len(s_l2)):
        if s_l2[i] < 0.05:
            v = Vt_l2[i]
            v_orig = v / norms_l2
            v_orig = v_orig / np.linalg.norm(v_orig)
            terms = sorted(zip(op_names_l2, v_orig), key=lambda x: -abs(x[1]))
            eq = " + ".join(f"({w:+.3f})*{n}" for n, w in terms if abs(w) > 0.05)
            print(f"    sigma={s_l2[i]:.6f}: {eq}")
else:
    print(f"\n  No Level 2 algebraic constraint detected in this operator basis.")
    print(f"  Smallest singular value = {s_l2[-1]:.4f} (above 0.05 threshold).")

# Also try a richer Level 2 operator set: include log(B_k)^2, log(B_k)*log(B_j) over neighbor pairs
print(f"\n  Trying richer operator basis (with neighbor cross-terms)...")
all_rows_rich = []
for w in range(n_w):
    bw = Bl[w]
    for k in range(n_pairs):
        bk = bw[k]
        mask = np.ones(n_pairs, bool); mask[k] = False
        b_neigh = bw[mask]
        b_net = b_neigh.mean()
        all_rows_rich.append([
            bk, b_net, bk**2, b_net**2, bk*b_net,
            bk**3, bk**2*b_net, bk*b_net**2,
            (b_neigh**2).mean(), b_neigh.std()
        ])
op_names_rich = ['B_k', 'B_net', 'B_k^2', 'B_net^2', 'B_k*B_net',
                 'B_k^3', 'B_k^2*B_net', 'B_k*B_net^2',
                 'mean(B^2)', 'std(B)']
op_mat_r = np.array(all_rows_rich)
op_c_r = op_mat_r - op_mat_r.mean(axis=0)
norms_r = np.linalg.norm(op_c_r, axis=0) + 1e-15
op_n_r = op_c_r / norms_r
_, S_r, Vt_r = svd(op_n_r, full_matrices=False)
s_r = S_r / S_r.max()
rank_r = (s_r > 0.05).sum()
print(f"  Singular values: " + " ".join(f"{s:.3f}" for s in s_r))
print(f"  Effective rank = {rank_r} / {len(s_r)}")
if rank_r < len(s_r):
    print(f"  Level 2 constraint directions (richer basis):")
    for i in range(len(s_r)):
        if s_r[i] < 0.05:
            v = Vt_r[i]
            v_orig = v / norms_r
            v_orig = v_orig / np.linalg.norm(v_orig)
            terms = sorted(zip(op_names_rich, v_orig), key=lambda x: -abs(x[1]))
            eq = " + ".join(f"({w:+.3f})*{n}" for n, w in terms if abs(w) > 0.05)
            print(f"    sigma={s_r[i]:.6f}: {eq}")

# ============================================================
# DISCOVERY 3: BIOLOGY IN MI-SPACE
# ============================================================
print(f"\n{'=' * 70}")
print("  DISCOVERY 3: BIOLOGY IN MI-SPACE")
print(f"{'=' * 70}")
print("""
  Prior result: biology multivar gain = +0.27, dominated by MI
  and H_a*MI operators.  Re-test using MI as the primary observable,
  not H_a^2.  Try MI^2 vs MI*H_a, MI^2 vs H_a*H_b, etc.
""")

for name, (H_a, H_b, MI) in systems_data.items():
    # Test several relationships in MI-space
    cands = [
        ('MI^2 ~ MI*H_a',     MI*H_a,        MI**2),
        ('MI^2 ~ MI*H_b',     MI*H_b,        MI**2),
        ('MI^2 ~ H_a*H_b',    H_a*H_b,       MI**2),
        ('MI ~ H_a*H_b',      H_a*H_b,       MI),
        ('H_a*MI ~ H_b*MI',   H_b*MI,        H_a*MI),
        ('MI^2 ~ H_a^2+H_b^2', H_a**2+H_b**2, MI**2),
    ]
    print(f"\n  {name}")
    for label, xv, yv in cands:
        # quadratic fit
        lib = np.column_stack([np.ones(len(xv)), xv, xv**2])
        c, *_ = lstsq(lib, yv, rcond=None)
        r2v = r2(yv, lib @ c)
        # also linear
        liblin = np.column_stack([np.ones(len(xv)), xv])
        clin, *_ = lstsq(liblin, yv, rcond=None)
        r2lin = r2(yv, liblin @ clin)
        pearson = np.corrcoef(xv, yv)[0,1]
        flag = " <--" if r2v > 0.8 else ""
        print(f"    {label:<25s}  Pearson={pearson:+.3f}  R^2(lin)={r2lin:.3f}  R^2(quad)={r2v:.3f}{flag}")

# Also test biology's multivar fit with MI*x as target instead of H_a^2
print(f"\n  Biology focus: regress MI itself on full operator library...")
H_a, H_b, MI = systems_data['BIOLOGY']
y = MI
lib = np.column_stack([
    H_a, H_b, H_a*H_b, H_a**2, H_b**2, H_a*H_b**2, H_a**2*H_b, np.ones(len(MI))
])
labels = ['H_a','H_b','H_a*H_b','H_a^2','H_b^2','H_a*H_b^2','H_a^2*H_b','const']
c, *_ = lstsq(lib, y, rcond=None)
r2_mi = r2(y, lib @ c)
print(f"  R^2 fitting MI(t) = polynomial in H_a, H_b: {r2_mi:.4f}")
items = [(k, v) for k, v in zip(labels, c) if k != 'const' and abs(v) > 1e-6]
items.sort(key=lambda x: -abs(x[1]))
print(f"  Top operators:")
for nm, cv in items[:5]:
    print(f"    {nm:<14s} {cv:+.4f}")

# ============================================================
# FIGURE
# ============================================================
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# top-left: chemistry quadratic
H_a, H_b, MI = systems_data['CHEMISTRY']
x = H_a * H_b; y = H_a ** 2
ax = axes[0, 0]
ax.scatter(x, y, s=4, alpha=0.4, color='#dc2626')
xs = np.linspace(x.min(), x.max(), 200)
lib = np.column_stack([np.ones(len(x)), x, x**2])
c, *_ = lstsq(lib, y, rcond=None)
ys = c[0] + c[1]*xs + c[2]*xs**2
ax.plot(xs, ys, 'k-', linewidth=1.5)
ax.set_title(f"CHEMISTRY: H_a^2 = quadratic(H_a*H_b)\nR^2 = {r2_chem:.4f}", fontsize=10)
ax.set_xlabel("H_a * H_b"); ax.set_ylabel("H_a^2"); ax.grid(alpha=0.3)

# top-middle: geology singular value spectrum
ax = axes[0, 1]
ax.semilogy(range(1, len(s_norm)+1), s_norm, 'o-', color='#059669', markersize=8)
ax.axhline(0.05, color='red', linestyle='--', alpha=0.5, label='constraint threshold')
ax.set_title(f"GEOLOGY: singular value spectrum\nrank = {(s_norm > 0.05).sum()}/6", fontsize=10)
ax.set_xlabel("singular value index"); ax.set_ylabel("sigma / sigma_max")
ax.grid(alpha=0.3); ax.legend()

# top-right: Level 2 singular value spectrum
ax = axes[0, 2]
ax.semilogy(range(1, len(s_l2)+1), s_l2, 'o-', color='#7c3aed', markersize=8, label='basic basis')
ax.semilogy(range(1, len(s_r)+1), s_r, 's-', color='#ea580c', markersize=6, label='rich basis')
ax.axhline(0.05, color='red', linestyle='--', alpha=0.5)
ax.set_title(f"LEVEL 2: dipole network SVD\nrank(basic)={rank_l2}, rank(rich)={rank_r}", fontsize=10)
ax.set_xlabel("singular value index"); ax.set_ylabel("sigma / sigma_max")
ax.grid(alpha=0.3); ax.legend()

# bottom-left: biology MI vs H_a*H_b (best MI relationship)
H_a, H_b, MI = systems_data['BIOLOGY']
ax = axes[1, 0]
ax.scatter(H_a*H_b, MI**2, s=4, alpha=0.4, color='#0891b2')
xs = np.linspace((H_a*H_b).min(), (H_a*H_b).max(), 200)
lib = np.column_stack([np.ones(len(H_a)), H_a*H_b, (H_a*H_b)**2])
c, *_ = lstsq(lib, MI**2, rcond=None)
ys = c[0] + c[1]*xs + c[2]*xs**2
ax.plot(xs, ys, 'k-', linewidth=1.5)
r2_bio = r2(MI**2, lib @ c)
ax.set_title(f"BIOLOGY in MI-space\nMI^2 = quadratic(H_a*H_b)\nR^2 = {r2_bio:.3f}", fontsize=10)
ax.set_xlabel("H_a * H_b"); ax.set_ylabel("MI^2"); ax.grid(alpha=0.3)

# bottom-middle: chemistry residuals
H_a, H_b, MI = systems_data['CHEMISTRY']
x = H_a * H_b; y = H_a**2
lib = np.column_stack([np.ones(len(x)), x, x**2])
c, *_ = lstsq(lib, y, rcond=None)
resid = y - lib @ c
ax = axes[1, 1]
ax.scatter(x, resid, s=3, alpha=0.4, color='#dc2626')
ax.axhline(0, color='black', linewidth=0.8)
ax.set_title(f"CHEMISTRY residuals\n(should be flat noise)", fontsize=10)
ax.set_xlabel("H_a * H_b"); ax.set_ylabel("residual"); ax.grid(alpha=0.3)

# bottom-right: biology MI from H_a, H_b polynomial
H_a, H_b, MI = systems_data['BIOLOGY']
lib = np.column_stack([H_a, H_b, H_a*H_b, H_a**2, H_b**2, H_a*H_b**2, H_a**2*H_b, np.ones(len(MI))])
c, *_ = lstsq(lib, MI, rcond=None)
pred_MI = lib @ c
ax = axes[1, 2]
ax.scatter(MI, pred_MI, s=4, alpha=0.4, color='#0891b2')
mn, mx = MI.min(), MI.max()
ax.plot([mn, mx], [mn, mx], 'k--', linewidth=1)
ax.set_title(f"BIOLOGY: MI from H_a,H_b polynomial\nR^2 = {r2_mi:.3f}", fontsize=10)
ax.set_xlabel("MI observed"); ax.set_ylabel("MI predicted"); ax.grid(alpha=0.3)

plt.tight_layout()
out = '/mnt/user-data/outputs/three_discoveries.png'
plt.savefig(out, dpi=120, bbox_inches='tight')
print(f"\n  Figure: {out}")

print(f"\n{'=' * 70}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 70}")
