"""
NON-LINEAR DEEP DIVE
=====================
The biggest signal so far: Chemistry R^2 linear = 0.05 ->
R^2 quadratic = 0.94.  A curve exists, polynomial captured it.
But polynomials are universal approximators - the real question
is what FUNCTIONAL FAMILY actually describes the relationship.

Plan:
  1. For each of 4 sciences, fit H_a^2 vs H_a*H_b under multiple
     functional forms.  Look for one form that wins universally.
  2. Test multi-variable static fits: H_a^2 = f(H_a*H_b, H_b, MI).
     The 2D projection may be too lossy - the real constraint
     surface may live in higher dimensions.
  3. SVD/manifold rank analysis: how many effective dimensions
     does the operator point cloud occupy?  If the rank drops,
     there is an algebraic constraint regardless of the form.

Functional forms tested per system:
  linear:        y = a*x + b
  quadratic:     y = a + b*x + c*x^2
  cubic:         y = a + b*x + c*x^2 + d*x^3
  power-law:     y = a*x^b   (log-log linear)
  log:           y = a*log(x) + b
  sqrt:          y = a*sqrt(x) + b
  reciprocal:    y = a/x + b
  exp:           y = a*exp(b*x)   (semi-log linear)

The functional form that wins universally is the dipole's true
algebraic equation.  If no form wins, the constraint lives in
multi-variable space, not 2D.
"""
import numpy as np
from numpy.linalg import lstsq, svd
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
# FOUR SCIENCES (asymmetric variants to avoid symmetry confound)
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
    """Asymmetric variant: different friction and drive between segments."""
    steps = int(T/dt)
    u1 = 0.1*np.random.randn(N_ens); v1 = 0.05*np.random.randn(N_ens)
    u2 = 0.1*np.random.randn(N_ens); v2 = 0.05*np.random.randn(N_ens)
    H_a = np.zeros(steps); H_b = np.zeros(steps); MI = np.zeros(steps)
    for t in range(steps):
        fric1 = -0.6*v1/(1.0 + 8.0*v1**2)
        fric2 = -0.9*v2/(1.0 + 14.0*v2**2)    # different friction
        F1 = 0.5*(u2-u1) - 0.6*(u1 - 0.04*t*dt) + fric1
        F2 = 0.5*(u1-u2) - 0.45*(u2 - 0.03*t*dt) + fric2  # different K0, drive
        v1 += F1*dt; v2 += F2*dt
        u1 += v1*dt + 0.015*np.random.randn(N_ens)*np.sqrt(dt)
        u2 += v2*dt + 0.020*np.random.randn(N_ens)*np.sqrt(dt)
        H_a[t] = entropy_1d(u1); H_b[t] = entropy_1d(u2)
        MI[t] = mi_2d(u1, u2)
    return H_a, H_b, MI

# ============================================================
# FUNCTIONAL FORM BATTERY
# ============================================================
def fit_forms(x, y):
    """Try multiple functional forms.  Return dict of R^2 per form."""
    # mask out any non-positive values for log/power
    pos = (x > 0) & (y > 0)
    xp, yp = x[pos], y[pos]
    results = {}

    # 1. Linear  y = a*x + b
    lib = np.column_stack([x, np.ones(len(x))])
    c, *_ = lstsq(lib, y, rcond=None)
    results['linear']    = (r2(y, lib @ c), c)

    # 2. Quadratic  y = a + b*x + c*x^2
    lib = np.column_stack([np.ones(len(x)), x, x**2])
    c, *_ = lstsq(lib, y, rcond=None)
    results['quadratic'] = (r2(y, lib @ c), c)

    # 3. Cubic
    lib = np.column_stack([np.ones(len(x)), x, x**2, x**3])
    c, *_ = lstsq(lib, y, rcond=None)
    results['cubic']     = (r2(y, lib @ c), c)

    # 4. Power-law  log y = log a + b * log x
    lib = np.column_stack([np.log(xp), np.ones(len(xp))])
    c, *_ = lstsq(lib, np.log(yp), rcond=None)
    pred_log = lib @ c
    # R^2 reported in original-space
    pred = np.exp(pred_log)
    results['power']     = (r2(yp, pred), c)

    # 5. Log  y = a*log(x) + b
    lib = np.column_stack([np.log(xp), np.ones(len(xp))])
    c, *_ = lstsq(lib, yp, rcond=None)
    results['log']       = (r2(yp, lib @ c), c)

    # 6. Sqrt  y = a*sqrt(x) + b
    lib = np.column_stack([np.sqrt(xp), np.ones(len(xp))])
    c, *_ = lstsq(lib, yp, rcond=None)
    results['sqrt']      = (r2(yp, lib @ c), c)

    # 7. Reciprocal  y = a/x + b
    lib = np.column_stack([1.0/xp, np.ones(len(xp))])
    c, *_ = lstsq(lib, yp, rcond=None)
    results['reciprocal']= (r2(yp, lib @ c), c)

    # 8. Exponential  log y = log a + b*x
    lib = np.column_stack([x[pos], np.ones(pos.sum())])
    c, *_ = lstsq(lib, np.log(yp), rcond=None)
    pred = np.exp(lib @ c)
    results['exp']       = (r2(yp, pred), c)

    return results

# ============================================================
# MULTI-VARIABLE STATIC FIT
# ============================================================
def multivar_fit(H_a, H_b, MI):
    """
    Fit H_a^2 against a richer operator library STATICALLY.
    No derivatives.  Test whether the algebraic constraint
    lives in higher-dimensional space.
    """
    y = H_a**2
    # full operator library (same shape as Level 1 differential search)
    lib = np.column_stack([
        H_a, H_b, H_a*H_b, H_b**2, MI, H_a*MI, H_b*MI, MI**2, np.ones(len(H_a))
    ])
    labels = ['H_a','H_b','H_a*H_b','H_b^2','MI','H_a*MI','H_b*MI','MI^2','const']
    c, *_ = lstsq(lib, y, rcond=None)
    pred = lib @ c
    return r2(y, pred), dict(zip(labels, c))

# ============================================================
# MANIFOLD RANK (SVD on operator space)
# ============================================================
def manifold_rank(H_a, H_b, MI, threshold=0.01):
    """
    Compute SVD of the operator point cloud.  Effective rank
    counts singular values above threshold * max.  If rank < #cols,
    there is an algebraic constraint binding the operators.
    """
    ops = np.column_stack([
        H_a, H_b, H_a**2, H_b**2, H_a*H_b, MI
    ])
    # center
    ops_c = ops - ops.mean(axis=0)
    # normalize each column
    norms = np.linalg.norm(ops_c, axis=0) + 1e-15
    ops_n = ops_c / norms
    U, S, Vt = svd(ops_n, full_matrices=False)
    s_norm = S / S.max()
    eff_rank = (s_norm > threshold).sum()
    return eff_rank, S, Vt

# ============================================================
# RUN
# ============================================================
print("=" * 70)
print("  NON-LINEAR DEEP DIVE")
print("=" * 70)

systems_data = {}
print("\n  Simulating 4 sciences (asymmetric variants where applicable)...")
for name, fn in [('PHYSICS', physics), ('BIOLOGY', biology),
                 ('CHEMISTRY', chemistry), ('GEOLOGY', geology)]:
    t0 = time.time()
    H_a, H_b, MI = fn()
    print(f"    {name:<10s} done  [{time.time()-t0:.1f}s]")
    systems_data[name] = (H_a, H_b, MI)

# ------------------------------------------------------------
# 1. FUNCTIONAL FORM TABLE
# ------------------------------------------------------------
print(f"\n{'-' * 70}")
print("  TABLE 1: R^2 by functional form for H_a^2 vs H_a*H_b")
print(f"{'-' * 70}")
forms_list = ['linear','quadratic','cubic','power','log','sqrt','reciprocal','exp']
header = f"  {'system':<10s} " + ' '.join(f"{f:>9s}" for f in forms_list) + "   winner"
print(header)

best_forms = {}
all_form_results = {}
for name, (H_a, H_b, MI) in systems_data.items():
    x = H_a * H_b
    y = H_a ** 2
    res = fit_forms(x, y)
    all_form_results[name] = res
    row = f"  {name:<10s} "
    best_form, best_r2 = None, -np.inf
    for f in forms_list:
        r = res[f][0]
        row += f"{r:>9.3f} "
        if r > best_r2:
            best_r2 = r; best_form = f
    best_forms[name] = (best_form, best_r2)
    row += f"  {best_form} ({best_r2:.3f})"
    print(row)

# Universal winner?
print(f"\n  Best functional form per system:")
for name, (f, r) in best_forms.items():
    print(f"    {name:<10s} -> {f} (R^2 = {r:.3f})")

# ------------------------------------------------------------
# 2. MULTI-VARIABLE STATIC FIT
# ------------------------------------------------------------
print(f"\n{'-' * 70}")
print("  TABLE 2: Multi-variable static fit H_a^2 = f(H_a,H_b,MI,...)")
print(f"{'-' * 70}")
print(f"  Compares 2D projection (best of forms) vs richer operator library.\n")
print(f"  {'system':<10s} {'best 2D':>10s} {'multivar':>10s} {'gain':>8s}")
for name, (H_a, H_b, MI) in systems_data.items():
    r2_mv, coeffs = multivar_fit(H_a, H_b, MI)
    best_2d = best_forms[name][1]
    gain = r2_mv - best_2d
    print(f"  {name:<10s} {best_2d:>10.3f} {r2_mv:>10.3f} {gain:>+8.3f}")
    # show top operators
    items = [(k, v) for k, v in coeffs.items() if k != 'const' and abs(v) > 1e-6]
    items.sort(key=lambda x: -abs(x[1]))
    print(f"    top operators: " + ", ".join(f"{k}({v:+.3f})" for k, v in items[:4]))

# ------------------------------------------------------------
# 3. MANIFOLD RANK
# ------------------------------------------------------------
print(f"\n{'-' * 70}")
print("  TABLE 3: Manifold rank (SVD of normalized operator cloud)")
print(f"{'-' * 70}")
print(f"  6 operators total: H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI")
print(f"  rank < 6  ->  algebraic constraint binds them\n")
print(f"  {'system':<10s} {'eff rank':>9s}  {'singular values (normalized)':<40s}")
for name, (H_a, H_b, MI) in systems_data.items():
    rank, S, Vt = manifold_rank(H_a, H_b, MI)
    s_norm = S / S.max()
    s_str = " ".join(f"{s:.3f}" for s in s_norm)
    print(f"  {name:<10s} {rank:>9d}  [{s_str}]")
    # the smallest singular direction is the constraint!
    if rank < 6:
        # null direction = last row of Vt
        constraint = Vt[-1]
        op_names = ['H_a', 'H_b', 'H_a^2', 'H_b^2', 'H_a*H_b', 'MI']
        # report the dominant terms in the constraint
        terms = sorted(zip(op_names, constraint), key=lambda x: -abs(x[1]))
        active = [f"{w:+.3f}*{n}" for n, w in terms if abs(w) > 0.15]
        print(f"    constraint direction: {' '.join(active)} ~ 0")

# ------------------------------------------------------------
# FIGURE
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 4, figsize=(18, 9))
for col, name in enumerate(['PHYSICS','BIOLOGY','CHEMISTRY','GEOLOGY']):
    H_a, H_b, MI = systems_data[name]
    x = H_a * H_b; y = H_a**2

    # top: scatter with best form overlay
    ax = axes[0, col]
    ax.scatter(x, y, s=4, alpha=0.35, color='#2563eb')
    best_f, best_r = best_forms[name]
    xs = np.linspace(x.min(), x.max(), 200)
    coeffs = all_form_results[name][best_f][1]
    if best_f == 'linear':
        ys = coeffs[0]*xs + coeffs[1]
    elif best_f == 'quadratic':
        ys = coeffs[0] + coeffs[1]*xs + coeffs[2]*xs**2
    elif best_f == 'cubic':
        ys = coeffs[0] + coeffs[1]*xs + coeffs[2]*xs**2 + coeffs[3]*xs**3
    elif best_f == 'power':
        ys = np.exp(coeffs[1]) * np.maximum(xs, 1e-12)**coeffs[0]
    elif best_f == 'log':
        ys = coeffs[0]*np.log(np.maximum(xs, 1e-12)) + coeffs[1]
    elif best_f == 'sqrt':
        ys = coeffs[0]*np.sqrt(np.maximum(xs, 0)) + coeffs[1]
    elif best_f == 'reciprocal':
        ys = coeffs[0]/np.maximum(xs, 1e-12) + coeffs[1]
    elif best_f == 'exp':
        ys = np.exp(coeffs[0]*xs + coeffs[1])
    else:
        ys = np.full_like(xs, y.mean())
    ax.plot(xs, ys, 'r-', linewidth=1.5)
    ax.set_title(f"{name}: best={best_f}\nR^2={best_r:.3f}", fontsize=10)
    ax.set_xlabel("H_a * H_b"); ax.set_ylabel("H_a^2")
    ax.grid(alpha=0.3)

    # bottom: residuals from best fit
    ax = axes[1, col]
    coeffs = all_form_results[name][best_f][1]
    if best_f == 'linear':
        pred = coeffs[0]*x + coeffs[1]
    elif best_f == 'quadratic':
        pred = coeffs[0] + coeffs[1]*x + coeffs[2]*x**2
    elif best_f == 'cubic':
        pred = coeffs[0] + coeffs[1]*x + coeffs[2]*x**2 + coeffs[3]*x**3
    elif best_f == 'power':
        pred = np.exp(coeffs[1]) * np.maximum(x, 1e-12)**coeffs[0]
    elif best_f == 'log':
        pred = coeffs[0]*np.log(np.maximum(x, 1e-12)) + coeffs[1]
    else:
        pred = np.full_like(x, y.mean())
    resid = y - pred
    ax.scatter(x, resid, s=3, alpha=0.35, color='#7c3aed')
    ax.axhline(0, color='red', linewidth=0.8, alpha=0.6)
    ax.set_title(f"residuals (best {best_f})", fontsize=10)
    ax.set_xlabel("H_a * H_b"); ax.set_ylabel("y - pred")
    ax.grid(alpha=0.3)

plt.suptitle("Non-linear deep dive: best functional form per science",
             fontsize=13, y=1.00)
plt.tight_layout()
out = '/mnt/user-data/outputs/nonlinear_deep_dive.png'
plt.savefig(out, dpi=120, bbox_inches='tight')
print(f"\n  Figure: {out}")

print(f"\n{'=' * 70}")
print(f"  TOTAL RUNTIME: {time.time() - T_START:.1f}s")
print(f"{'=' * 70}")
