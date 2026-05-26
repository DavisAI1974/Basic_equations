"""Post-hoc analysis: singular spectrum and converged direction structure."""
import json
import numpy as np

with open('/home/claude/L1_linear_drift_NREAL_sweep.json') as f:
    R = json.load(f)

ATTRACTOR_SUB = np.array([-1, -1, 2]) / np.sqrt(6)

print("=== Singular value spectrum per cell ===")
for N_str in ['50', '100', '200', '500']:
    print(f"\nN_REAL = {N_str}")
    for cell in R[N_str]['cells']:
        sv = np.array(cell['singular_values'])
        ratios = sv[:-1] / (sv[1:] + 1e-20)
        print(f"  seed={cell['seed']}: SV = {[f'{x:.2e}' for x in sv]}")
        print(f"           gaps SV[i]/SV[i+1]: {[f'{r:.1f}' for r in ratios]}")

print("\n=== Converged direction at N_REAL=500 (subspace H_a^2, H_b^2, H_a*H_b) ===")
sub_500 = np.array([c['subspace'] for c in R['500']['cells']])
print(f"Per seed:")
for i, s in enumerate(sub_500):
    print(f"  seed={i}: {s.round(4)}")
mean_dir = sub_500.mean(axis=0)
mean_dir = mean_dir / np.linalg.norm(mean_dir)
print(f"Mean (normalized): {mean_dir.round(4)}")
print(f"\nAttractor direction (-1,-1,+2)/sqrt(6): {ATTRACTOR_SUB.round(4)}")
print(f"Cos(converged_mean, attractor): {np.dot(mean_dir, ATTRACTOR_SUB):.4f}")
print(f"Cos(converged_mean, -attractor): {np.dot(mean_dir, -ATTRACTOR_SUB):.4f}")

print("\n=== T=30 Session 2 reference (asymmetric pre-correction) ===")
ref_t30 = np.array([0.45, 0.38, -0.80])
ref_t30 = ref_t30 / np.linalg.norm(ref_t30)
print(f"(+0.45, +0.38, -0.80) normalized: {ref_t30.round(4)}")
print(f"Cos(converged_mean, T=30_ref): {np.dot(mean_dir, ref_t30):.4f}")
print(f"Cos(converged_mean, -T=30_ref): {np.dot(mean_dir, -ref_t30):.4f}")

print("\n=== Specific candidate directions to test ===")
# Could it be (+1, +1, +2)/sqrt(6)?
d1 = np.array([1, 1, 2]) / np.sqrt(6)
print(f"(+1,+1,+2)/sqrt(6) = {d1.round(4)}")
print(f"  cos with mean direction: {np.dot(mean_dir, d1):.4f}")

# Could it be (+1, +1, sqrt(2)*something)?
# Mean was approx (0.47, 0.37, 0.80) -> ratios 0.47/0.80 ~ 0.59, 0.37/0.80 ~ 0.46
# Not a simple integer ratio
# Try (1, 1, sqrt(3))/sqrt(5)
d2 = np.array([1, 1, np.sqrt(3)]) / np.sqrt(5)
print(f"(+1,+1,sqrt(3))/sqrt(5) = {d2.round(4)}")
print(f"  cos with mean direction: {np.dot(mean_dir, d2):.4f}")

# What does this direction correspond to as an operator combination?
# 0.47*H_a^2 + 0.37*H_b^2 + 0.80*H_a*H_b
# This factorizes... if we ignored asymmetry: ~0.4*(H_a^2 + H_b^2) + 0.8*H_a*H_b
# That's 0.4*((H_a + H_b)^2 + (H_a - H_b)^2)/... let's see
# H_a^2 + H_b^2 + 2*H_a*H_b = (H_a + H_b)^2
# So 0.4*(H_a^2 + H_b^2) + 0.8*H_a*H_b = 0.4*((H_a+H_b)^2)
# Almost! The converged direction is approximately proportional to (H_a + H_b)^2
print(f"\nNote: the converged direction (0.47, 0.37, 0.80) is approximately")
print(f"proportional to coefficients of (H_a + H_b)^2 = H_a^2 + 2*H_a*H_b + H_b^2")
print(f"Coefficients of (H_a+H_b)^2: (1, 1, 2)/sqrt(6) = {np.array([1,1,2])/np.sqrt(6)}")
print(f"Linear drift seems to land at: positive (H_a+H_b)^2 direction")
print(f"OU attractor lands at:         negative (H_a-H_b)^2 direction = (-1,-1,+2)/sqrt(6)")
# Check: (H_a - H_b)^2 = H_a^2 - 2*H_a*H_b + H_b^2, coefficients (1,1,-2)
# So (-1,-1,+2)/sqrt(6) is the negative of (H_a-H_b)^2 coefficients

print("\nVerification:")
print(f"(H_a+H_b)^2 coefficients [H_a^2, H_b^2, H_a*H_b] = (1, 1, 2)")
print(f"(H_a-H_b)^2 coefficients [H_a^2, H_b^2, H_a*H_b] = (1, 1, -2)")
print(f"OU attractor:    -(H_a-H_b)^2 direction")
print(f"Linear drift:    +(H_a+H_b)^2 direction (within sign/numeric precision)")

# Test whether the linear drift direction is closer to (H_a+H_b)^2 = (1,1,2)/sqrt(6) than to anything else
d_plus = np.array([1, 1, 2]) / np.sqrt(6)
d_minus = np.array([1, 1, -2]) / np.sqrt(6)
print(f"\ncos(linear drift mean, +(H_a+H_b)^2 dir):  {np.dot(mean_dir, d_plus):.4f}")
print(f"cos(linear drift mean, +(H_a-H_b)^2 dir):  {np.dot(mean_dir, d_minus):.4f}")
print(f"cos(linear drift mean, OU attractor):      {np.dot(mean_dir, ATTRACTOR_SUB):.4f}")
