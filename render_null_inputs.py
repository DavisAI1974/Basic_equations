"""
Render: combine probe_null_inputs_results.json (Session 5) and
probe_system_directions_results.json (Session 4) into a single
visual landscape of where the protocol lands.

Output:
  fig_null_alignment.png        bar chart, cos to (+,+,+) for all cases
  fig_null_geometry.png         3D scatter of mean_v3 vs reference vectors
  fig_null_sv_spectrum.png      log-scale SV spectrum, small multiples
  fig_null_interseed.png        inter-seed agreement per case
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

NULL_PATH = '/home/user/Basic_equations/probe_null_inputs_results.json'
SYS_PATH  = '/home/user/Basic_equations/probe_system_directions_results.json'

REF_PPP = np.array([+1.0, +1.0, +2.0]) / np.sqrt(6.0)
REF_MMP = np.array([-1.0, -1.0, +2.0]) / np.sqrt(6.0)

def load_cases():
    with open(SYS_PATH) as f:
        sys_data = json.load(f)
    with open(NULL_PATH) as f:
        null_data = json.load(f)
    cases = []
    for name, d in sys_data.items():
        cases.append(dict(
            name=name, source='session4',
            mean_v3=np.array(d['mean_v3']),
            cos_ppp=d['cos_to_plus_plus_plus'],
            cos_mmp=d['cos_to_minus_minus_plus'],
            inter_seed=d['inter_seed_cos_mean'],
            svs=[p['singular_values'] for p in d['per_seed']],
        ))
    for name, d in null_data.items():
        cases.append(dict(
            name=name, source='session5_null',
            mean_v3=np.array(d['mean_v3']),
            cos_ppp=d['cos_to_plus_plus_plus'],
            cos_mmp=d['cos_to_minus_minus_plus'],
            inter_seed=d['inter_seed_cos_mean'],
            svs=[p['singular_values'] for p in d['per_seed']],
        ))
    return cases

def fig_alignment(cases, out_path):
    """Bar chart of cos to (+,+,+) per case, sorted descending. Color by source."""
    cases_sorted = sorted(cases, key=lambda c: -c['cos_ppp'])
    names = [c['name'] for c in cases_sorted]
    cos_p = [c['cos_ppp'] for c in cases_sorted]
    cos_m = [c['cos_mmp'] for c in cases_sorted]
    colors = ['#2b6cb0' if c['source'] == 'session4' else '#c05621' for c in cases_sorted]

    fig, ax = plt.subplots(figsize=(12, 7))
    x = np.arange(len(names))
    w = 0.4
    ax.bar(x - w/2, cos_p, w, color=colors, label='cos to (+1,+1,+2)/sqrt(6)')
    ax.bar(x + w/2, cos_m, w, color='lightgray',
           edgecolor='gray', label='cos to (-1,-1,+2)/sqrt(6)')
    ax.axhline(0.95, ls='--', color='green', lw=1, alpha=0.6, label='|cos|=0.95')
    ax.axhline(0.0, color='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=55, ha='right', fontsize=9)
    ax.set_ylabel('cosine similarity', fontsize=11)
    ax.set_ylim(-0.5, 1.05)
    ax.set_title('Where the protocol lands: cos to (+,+,+) vs (-,-,+)\n'
                 'blue = Session 4 system_directions probe   orange = Session 5 null_inputs probe',
                 fontsize=11)
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()

def fig_geometry(cases, out_path):
    """3D scatter of mean_v3 in [H_a^2, H_b^2, H_a*H_b] subspace, with refs."""
    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')

    # unit sphere outline for context (thin wireframe)
    u = np.linspace(0, 2*np.pi, 24)
    v = np.linspace(0, np.pi, 12)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color='lightgray', alpha=0.15, linewidth=0.4)

    # reference vectors
    for ref, label, color in [(REF_PPP, '(+1,+1,+2)/sqrt6', 'green'),
                              (REF_MMP, '(-1,-1,+2)/sqrt6', 'red')]:
        ax.plot([0, ref[0]], [0, ref[1]], [0, ref[2]],
                color=color, lw=2.5, alpha=0.9)
        ax.scatter([ref[0]], [ref[1]], [ref[2]],
                   color=color, s=120, marker='*', label=label, edgecolor='black')

    # case points
    for c in cases:
        v3 = c['mean_v3']
        col = '#2b6cb0' if c['source'] == 'session4' else '#c05621'
        ax.scatter([v3[0]], [v3[1]], [v3[2]], color=col, s=55,
                   edgecolor='black', linewidth=0.5, alpha=0.85)
        ax.text(v3[0]*1.06, v3[1]*1.06, v3[2]*1.06, c['name'],
                fontsize=7, color='black')

    ax.set_xlabel('H_a^2 (proj)')
    ax.set_ylabel('H_b^2 (proj)')
    ax.set_zlabel('H_a*H_b (proj)')
    ax.set_title('mean_v3 directions on the [H_a^2, H_b^2, H_a*H_b] unit sphere\n'
                 'green star = Session 4 attractor   red star = INFO-014 reference\n'
                 'blue = system probe   orange = null-input probe', fontsize=10)
    ax.legend(loc='upper left', fontsize=9)
    ax.view_init(elev=22, azim=42)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()

def fig_sv_spectrum(cases, out_path):
    """Log SVs per case, small multiples in a grid."""
    n = len(cases)
    cols = 4
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols*3.4, rows*2.6), sharey=True)
    axes = axes.flatten()
    for i, c in enumerate(cases):
        ax = axes[i]
        col = '#2b6cb0' if c['source'] == 'session4' else '#c05621'
        for sv in c['svs']:
            ax.semilogy(range(1, len(sv)+1), sv, '-o', color=col, alpha=0.6, ms=4)
        ax.set_title(c['name'], fontsize=8)
        ax.set_xticks([1, 2, 3, 4, 5, 6])
        ax.grid(True, which='both', alpha=0.3)
    for j in range(n, len(axes)):
        axes[j].axis('off')
    fig.suptitle('Singular value spectrum per case (log scale, 3 seeds overlaid)',
                 fontsize=12, y=1.00)
    fig.text(0.5, -0.01, 'singular value index', ha='center', fontsize=10)
    fig.text(-0.005, 0.5, 'singular value', va='center', rotation=90, fontsize=10)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140, bbox_inches='tight')
    plt.close()

def fig_interseed(cases, out_path):
    """Inter-seed cos per case, sorted."""
    cases_sorted = sorted(cases, key=lambda c: c['inter_seed'])
    names = [c['name'] for c in cases_sorted]
    vals  = [c['inter_seed'] for c in cases_sorted]
    colors = ['#2b6cb0' if c['source'] == 'session4' else '#c05621' for c in cases_sorted]

    fig, ax = plt.subplots(figsize=(9, 8))
    y = np.arange(len(names))
    ax.barh(y, vals, color=colors, edgecolor='black', linewidth=0.4)
    ax.axvline(0.95, ls='--', color='green', lw=1, alpha=0.6, label='cos=0.95')
    ax.axvline(0.0, color='black', lw=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('mean inter-seed cosine', fontsize=11)
    ax.set_xlim(-0.1, 1.05)
    ax.set_title('Inter-seed agreement per case (3 seeds, N_REAL=500)\n'
                 'blue = system probe   orange = null-input probe', fontsize=10)
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=140)
    plt.close()

def main():
    cases = load_cases()
    print(f"Loaded {len(cases)} cases.")
    fig_alignment(cases, '/home/user/Basic_equations/fig_null_alignment.png')
    fig_geometry (cases, '/home/user/Basic_equations/fig_null_geometry.png')
    fig_sv_spectrum(cases, '/home/user/Basic_equations/fig_null_sv_spectrum.png')
    fig_interseed(cases, '/home/user/Basic_equations/fig_null_interseed.png')
    print("Wrote 4 figures.")

if __name__ == "__main__":
    main()
