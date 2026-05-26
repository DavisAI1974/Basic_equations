"""
Render: combined view of Session 5 threads 1-3.

Outputs:
  fig_thread1_bypass.png       per-seed v3 scatter for three M-noise
                                variants vs baseline reference (3D + bars)
  fig_thread1_eigspectrum.png  log eigenvalue spectrum of baseline op cov
                                showing rank-3 structure
  fig_thread2_projsweep.png    per-subspace projected direction, bar chart
  fig_thread3_bimodal.png      30-seed per-seed v3 for scrambled_baseline_joint
                                on the unit sphere (bimodal structure)
  fig_thread3_interseed_hist.png inter-seed cos histogram bimodality
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

REF_PPP = np.array([1, 1, 2]) / np.sqrt(6)
REF_MMP = np.array([-1, -1, 2]) / np.sqrt(6)

def draw_unit_sphere(ax, alpha=0.1):
    u = np.linspace(0, 2*np.pi, 30)
    v = np.linspace(0, np.pi, 15)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color='lightgray', alpha=alpha, linewidth=0.3)

def draw_refs(ax):
    for ref, color, label in [(REF_PPP, 'green', '(+1,+1,+2)/sqrt6'),
                              (REF_MMP, 'red',   '(-1,-1,+2)/sqrt6')]:
        ax.plot([0, ref[0]], [0, ref[1]], [0, ref[2]], color=color, lw=2.2, alpha=0.9)
        ax.scatter([ref[0]], [ref[1]], [ref[2]], color=color, s=160,
                   marker='*', label=label, edgecolor='black', zorder=10)

# ---------------------------------------------------------------
# Thread 1
# ---------------------------------------------------------------
def fig_thread1_bypass():
    d = json.load(open('probe_operator_noise_results.json'))
    variants = ['iid_unit', 'iid_col_scaled', 'cov_matched']
    colors = {'iid_unit': '#888888', 'iid_col_scaled': '#d97757', 'cov_matched': '#7c3aed'}

    fig = plt.figure(figsize=(15, 6))

    # Left: 3D scatter of per-seed v3
    ax3d = fig.add_subplot(1, 2, 1, projection='3d')
    draw_unit_sphere(ax3d)
    draw_refs(ax3d)
    for var in variants:
        per_seed = d[var]['per_seed']
        v3s = np.array([p['v3'] for p in per_seed])
        ax3d.scatter(v3s[:,0], v3s[:,1], v3s[:,2], color=colors[var],
                     s=45, alpha=0.85, edgecolor='black', linewidth=0.4,
                     label=var)
    ax3d.set_xlabel('H_a^2 (proj)')
    ax3d.set_ylabel('H_b^2 (proj)')
    ax3d.set_zlabel('H_a*H_b (proj)')
    ax3d.legend(loc='upper left', fontsize=9)
    ax3d.set_title('Thread 1: per-seed v3 for synthetic-M variants\n'
                   'cov_matched (purple) clusters at baseline direction; iid (gray) scattered',
                   fontsize=10)
    ax3d.view_init(elev=22, azim=42)

    # Right: bar chart of metrics
    ax = fig.add_subplot(1, 2, 2)
    metrics = ['cos to (+,+,+)', 'inter-seed cos']
    x = np.arange(len(variants))
    w = 0.35
    cos_vals = [d[v]['cos_to_plus_plus_plus'] for v in variants]
    inter_vals = [d[v]['inter_seed_cos_mean'] for v in variants]
    ax.bar(x - w/2, cos_vals, w, color=[colors[v] for v in variants],
           edgecolor='black', label='cos(mean v3, (+,+,+))')
    ax.bar(x + w/2, inter_vals, w, color=[colors[v] for v in variants],
           edgecolor='black', alpha=0.5, hatch='//', label='mean inter-seed cos')
    ax.axhline(0.95, ls='--', color='green', lw=1, alpha=0.6)
    ax.axhline(0, color='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(variants)
    ax.set_ylabel('cosine')
    ax.set_ylim(-0.1, 1.05)
    ax.legend(loc='upper left', fontsize=9)
    ax.set_title('Thread 1: alignment + inter-seed agreement\n'
                 '(green dashed = 0.95 reference)', fontsize=10)
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('/home/user/Basic_equations/fig_thread1_bypass.png', dpi=140)
    plt.close()

def fig_thread1_eigspectrum():
    d = json.load(open('probe_operator_noise_results.json'))
    cov = np.array(d['baseline_cov'])
    w, V = np.linalg.eigh(cov)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(range(1, 7), w, 'o-', color='#2b6cb0', markersize=12, lw=2)
    for i, val in enumerate(w):
        ax.annotate(f'{val:.2e}', (i+1, val), textcoords="offset points",
                    xytext=(15, 0), fontsize=9, color='black')
    ax.axhspan(1e-14, 1e-10, alpha=0.15, color='red', label='effectively zero')
    ax.set_xlabel('eigenvalue index (smallest to largest)', fontsize=11)
    ax.set_ylabel('eigenvalue (log scale)', fontsize=11)
    ax.set_xticks([1, 2, 3, 4, 5, 6])
    ax.set_title('Thread 1: eigenvalue spectrum of baseline operator covariance\n'
                 'three eigenvalues at machine epsilon  ->  effective rank 3',
                 fontsize=11)
    ax.legend(loc='center right')
    ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.savefig('/home/user/Basic_equations/fig_thread1_eigspectrum.png', dpi=140)
    plt.close()

# ---------------------------------------------------------------
# Thread 2
# ---------------------------------------------------------------
def fig_thread2_projsweep():
    d = json.load(open('probe_projection_sweep_results.json'))
    # Sort subspaces by how H_a / H_b appear
    labels = list(d.keys())
    cos_ppp_sqrt6 = [d[k]['cos_to_refs']['(+1,+1,+2)/sqrt6'] for k in labels]
    cos_ppz_sqrt2 = [d[k]['cos_to_refs']['(+1,+1, 0)/sqrt2'] for k in labels]
    cos_0_0_p1   = [d[k]['cos_to_refs']['( 0, 0,+1)'] for k in labels]

    # Mark which subspaces contain H_a and H_b
    def cat(k):
        has_a = 'H_a,' in k or k.startswith('[H_a,') or '[H_a]' in k or k.find('[H_a,') >= 0
        has_b = 'H_b' in k.replace('H_b^2', '').replace('H_b*', '') or 'H_b,' in k or '[H_b,' in k
        # cleaner: inspect indices
        idx = d[k]['indices']
        ha, hb = 0 in idx, 1 in idx
        if ha and hb: return 'both H_a, H_b'
        if ha or hb: return 'one of H_a, H_b'
        return 'neither H_a nor H_b'

    cats = [cat(k) for k in labels]
    color_map = {'both H_a, H_b': '#2b6cb0',
                 'one of H_a, H_b': '#888888',
                 'neither H_a nor H_b': '#c05621'}
    colors = [color_map[c] for c in cats]

    # Sort by category then by cos to (+,+,+2)
    order = sorted(range(len(labels)),
                   key=lambda i: (cats[i], -cos_ppp_sqrt6[i]))
    labels = [labels[i] for i in order]
    cos_ppp_sqrt6 = [cos_ppp_sqrt6[i] for i in order]
    cos_ppz_sqrt2 = [cos_ppz_sqrt2[i] for i in order]
    colors = [colors[i] for i in order]

    fig, ax = plt.subplots(figsize=(13, 9))
    y = np.arange(len(labels))
    w = 0.4
    ax.barh(y - w/2, cos_ppp_sqrt6, w, color=colors, edgecolor='black',
            linewidth=0.4, label='cos to (+1,+1,+2)/sqrt6')
    ax.barh(y + w/2, cos_ppz_sqrt2, w, color=colors, edgecolor='black',
            alpha=0.5, hatch='//', linewidth=0.4, label='cos to (+1,+1,0)/sqrt2')
    ax.axvline(0.95, ls='--', color='green', lw=1, alpha=0.6)
    ax.axvline(-0.95, ls='--', color='green', lw=1, alpha=0.6)
    ax.axvline(0, color='black', lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel('cosine', fontsize=11)

    # custom legend for color categories
    from matplotlib.patches import Patch
    legend_items = [Patch(facecolor=color_map[c], edgecolor='black', label=c)
                    for c in ['both H_a, H_b', 'one of H_a, H_b', 'neither H_a nor H_b']]
    ax.legend(handles=legend_items + [Patch(facecolor='white', edgecolor='black',
              hatch='//', label='hatched: cos to (+1,+1,0)/sqrt2')],
              loc='lower right', fontsize=9)
    ax.set_title('Thread 2: projection-rule sweep across all 20 3D subspaces\n'
                 'subspaces containing both H_a and H_b lock onto (+1,+1,0)/sqrt2\n'
                 'only [H_a^2, H_b^2, H_a*H_b] gives the (+,+,+2)/sqrt6 reading',
                 fontsize=11)
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig('/home/user/Basic_equations/fig_thread2_projsweep.png', dpi=140)
    plt.close()

# ---------------------------------------------------------------
# Thread 3
# ---------------------------------------------------------------
def fig_thread3_bimodal():
    d = json.load(open('probe_scrambled_highseed_results.json'))
    per_seed = d['per_seed']
    v3s = np.array([p['v3'] for p in per_seed])
    cos_mmp = np.array([p['cos_mmp'] for p in per_seed])

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')
    draw_unit_sphere(ax)
    draw_refs(ax)
    # color by cos to (-,-,+) to show the bimodal structure
    sc = ax.scatter(v3s[:,0], v3s[:,1], v3s[:,2], c=cos_mmp,
                    cmap='RdYlBu_r', s=80, alpha=0.92, edgecolor='black',
                    linewidth=0.5, vmin=-0.1, vmax=1.0)
    for i, v in enumerate(v3s):
        ax.text(v[0]*1.06, v[1]*1.06, v[2]*1.06, str(i),
                fontsize=7, color='black')
    cbar = plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label('cos to (-1,-1,+2)/sqrt(6)', fontsize=9)
    ax.set_xlabel('H_a^2 (proj)')
    ax.set_ylabel('H_b^2 (proj)')
    ax.set_zlabel('H_a*H_b (proj)')
    ax.set_title('Thread 3: scrambled_baseline_joint per-seed v3, 30 seeds\n'
                 'two clusters: H_a*H_b-dominant (red, near (-,-,+))\n'
                 'and H_a^2-H_b^2-dominant (blue, near (+1,-1,0) artifact)',
                 fontsize=10)
    ax.view_init(elev=22, azim=42)
    plt.tight_layout()
    plt.savefig('/home/user/Basic_equations/fig_thread3_bimodal.png', dpi=140)
    plt.close()

def fig_thread3_interseed_hist():
    d = json.load(open('probe_scrambled_highseed_results.json'))
    pair_cos = np.array(d['inter_seed_cos_pairs'])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(pair_cos, bins=40, color='#c05621', edgecolor='black', alpha=0.85)
    ax.axvline(0, color='black', lw=1)
    ax.axvline(pair_cos.mean(), ls='--', color='green',
               label=f'mean = {pair_cos.mean():+.3f}')
    ax.axvline(np.median(pair_cos), ls=':', color='blue',
               label=f'median = {np.median(pair_cos):+.3f}')
    ax.set_xlabel('inter-seed cosine (signed, 435 pairs)', fontsize=11)
    ax.set_ylabel('count')
    ax.set_title('Thread 3: inter-seed cos for scrambled_baseline_joint at 30 seeds\n'
                 'bimodal: 216 pairs > +0.5, 112 pairs < -0.5\n'
                 '(SVD picking arbitrarily within a degenerate 3D null subspace)',
                 fontsize=11)
    ax.legend(loc='upper center', fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig('/home/user/Basic_equations/fig_thread3_interseed_hist.png', dpi=140)
    plt.close()

def main():
    fig_thread1_bypass()
    fig_thread1_eigspectrum()
    fig_thread2_projsweep()
    fig_thread3_bimodal()
    fig_thread3_interseed_hist()
    print('Wrote 5 figures.')

if __name__ == "__main__":
    main()
