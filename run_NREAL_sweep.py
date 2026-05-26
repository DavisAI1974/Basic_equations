"""Run the N_REAL sweep on linear drift, three seeds per N_REAL."""
import numpy as np
import json
import time
import sys
sys.path.insert(0, '/home/claude')
from extract_v1_linear_drift import extract_v1, project_and_compare, inter_seed_cosine, ATTRACTOR_SUB

N_REAL_LIST = [50, 100, 200, 500]
SEEDS = [0, 1, 2]

results = {}
t_start = time.time()

for N_REAL in N_REAL_LIST:
    print(f"\n=== N_REAL = {N_REAL} ===")
    t0 = time.time()
    cell_results = []
    null_full_list = []
    sub_list = []
    cos_attr_list = []
    sing_vals_list = []
    for s in SEEDS:
        r = extract_v1(N_REAL=N_REAL, master_seed=s)
        null_full_list.append(r['null_full'])
        sub, cos_a = project_and_compare(r['null_full'])
        sub_list.append(sub)
        cos_attr_list.append(cos_a)
        sing_vals_list.append(r['singular_values'].tolist())
        print(f"  seed={s}: subspace={sub.round(3).tolist()}, cos_to_attractor={cos_a:+.4f}")
        cell_results.append({
            'seed': s,
            'null_full': r['null_full'].tolist(),
            'subspace': sub.tolist(),
            'cos_to_attractor': cos_a,
            'singular_values': r['singular_values'].tolist(),
        })

    # Inter-seed cosines (full 6-vec, sign-blind)
    full_inter = inter_seed_cosine(null_full_list)
    # Inter-seed cosines (subspace 3-vec, sign-blind)
    sub_inter = inter_seed_cosine(sub_list)

    print(f"  inter-seed cos FULL (6d, |.|): {[round(x, 4) for x in full_inter]}")
    print(f"  inter-seed cos SUB  (3d, |.|): {[round(x, 4) for x in sub_inter]}")
    print(f"  cos-to-attractor:              {[round(x, 4) for x in cos_attr_list]}")
    print(f"  time: {time.time() - t0:.1f}s")

    results[N_REAL] = {
        'cells': cell_results,
        'inter_seed_cos_full': full_inter,
        'inter_seed_cos_sub': sub_inter,
        'inter_seed_cos_full_mean': float(np.mean(full_inter)),
        'inter_seed_cos_sub_mean': float(np.mean(sub_inter)),
        'inter_seed_cos_full_min': float(min(full_inter)),
        'inter_seed_cos_sub_min': float(min(sub_inter)),
        'cos_to_attractor_values': cos_attr_list,
        'cos_to_attractor_mean': float(np.mean(cos_attr_list)),
    }

print(f"\nTotal wall time: {time.time() - t_start:.1f}s")

# Save
out_path = '/home/claude/L1_linear_drift_NREAL_sweep.json'
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Saved to {out_path}")

# Summary table
print("\n=== SUMMARY: monotonicity check ===")
print(f"{'N_REAL':>8} {'inter_full_mean':>16} {'inter_full_min':>16} {'inter_sub_mean':>16} {'inter_sub_min':>16} {'cos_attr_mean':>15}")
for N_REAL in N_REAL_LIST:
    r = results[N_REAL]
    print(f"{N_REAL:>8} "
          f"{r['inter_seed_cos_full_mean']:>16.4f} "
          f"{r['inter_seed_cos_full_min']:>16.4f} "
          f"{r['inter_seed_cos_sub_mean']:>16.4f} "
          f"{r['inter_seed_cos_sub_min']:>16.4f} "
          f"{r['cos_to_attractor_mean']:>15.4f}")

# Monotonicity decision
full_means = [results[N]['inter_seed_cos_full_mean'] for N in N_REAL_LIST]
sub_means = [results[N]['inter_seed_cos_sub_mean'] for N in N_REAL_LIST]
mono_full = all(full_means[i] <= full_means[i+1] for i in range(len(full_means)-1))
mono_sub = all(sub_means[i] <= sub_means[i+1] for i in range(len(sub_means)-1))
print(f"\nMonotonic improvement (full 6d): {mono_full}")
print(f"Monotonic improvement (sub 3d): {mono_sub}")
print(f"Range full means: {min(full_means):.4f} to {max(full_means):.4f}")
print(f"Range sub means:  {min(sub_means):.4f} to {max(sub_means):.4f}")
