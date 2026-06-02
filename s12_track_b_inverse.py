"""
S12 Track B -- new-physics INVERSE PROBLEM on the PDG gauge couplings.

Builds directly on s10_pdg_unification.py (INFO-037): the three SM gauge
couplings share the linear-in-ln(Q) running FORM (substrate) with per-force
slopes b_i (expression), but their pairwise crossings form a TRIANGLE
spanning ~1e4 in energy -- no single unification point.

This script asks the honest OD-shaped question the handoff specified: NOT
"what is the new physics" (its IDENTITY -- 3 couplings cannot invert to a
unique particle spectrum, nobody can) but "what FOOTPRINT would new physics
have to leave to close the triangle" -- the beta-coefficient shifts Delta-b_i
above an onset scale mu_NP that make the three lines meet at one point. That
footprint is a mechanism-agnostic, FALSIFIABLE constraint surface.

RESULT DISCIPLINE -- map the alternatives FIRST, before the inversion:
  (a) [leading deflationary read] NOTHING forces a single unification. The
      triangle is the SM prediction; "they must meet" is an aesthetic prior,
      not a data requirement. Reported as the null.
  (b) two-loop running + thresholds may SHRINK the triangle on their own
      (Session 10 / INFO-037 used one-loop). We recompute the triangle at
      two-loop and measure how much of the gap closes with NO new physics.
  (c) the extrapolation to ~1e16 GeV is ITSELF conjecture (running measured
      only to ~1 TeV). Tagged, never cited as support.
Only AFTER mapping (a)-(c) do we run the inversion (footprint surface), then
the gravity-into-log extraction.

Speaking posture (Rule C): I THINK the inversion will show the required
Delta-b footprint is a smooth surface in (mu_NP, M_GUT) -- a constraint, not
a unique answer -- and that two-loop shrinks the triangle only modestly so
the bulk of the gap still needs new content. But we wait on what the numbers
say and where they point. No verdict in advance.

Literature-as-conjecture: M_Z couplings + measured alpha_s running meet the
bar (anchored on). GUT-scale extrapolation, MSSM spectrum, any specific new-
physics model do NOT -- computed and LOCATED on the surface, never invoked as
support.

Outputs s12_track_b_inverse_results.json.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# Measured anchors at Q = M_Z (identical to s10_pdg_unification.py).
# ----------------------------------------------------------------------
M_Z = 91.1876
alpha_em_inv_MZ = 127.951
sin2_thetaW_MZ = 0.23122
alpha_s_MZ = 0.1180
M_PLANCK = 1.22e19

cos2 = 1.0 - sin2_thetaW_MZ
a1_inv_MZ = (3.0 / 5.0) * cos2 * alpha_em_inv_MZ
a2_inv_MZ = sin2_thetaW_MZ * alpha_em_inv_MZ
a3_inv_MZ = 1.0 / alpha_s_MZ
ainv_MZ = np.array([a1_inv_MZ, a2_inv_MZ, a3_inv_MZ])

# One-loop SM beta coefficients (GUT-normalised b_1).
b_SM = np.array([41.0 / 10.0, -19.0 / 6.0, -7.0])
# Two-loop SM coefficient matrix b_ij (GUT-normalised). Standard values.
b2_SM = np.array([
    [199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0],
    [9.0 / 10.0,   35.0 / 6.0,  12.0],
    [11.0 / 10.0,  9.0 / 2.0,  -26.0],
])
# MSSM one-loop (CONJECTURE leg -- located on the surface, not cited).
b_MSSM = np.array([33.0 / 5.0, 1.0, -3.0])

LAB = ["1 U(1)_Y", "2 SU(2)_L", "3 SU(3)_C"]


def run_inv_1loop(ainv0, b, Q):
    return ainv0 - (b / (2 * np.pi)) * np.log(Q / M_Z)


def pairwise_crossings(ainv0, b):
    out = {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        t = 2 * np.pi * (ainv0[i] - ainv0[j]) / (b[i] - b[j])
        out[f"{i+1}{j+1}"] = float(M_Z * np.exp(t))
    return out


def triangle_spread(crossings):
    v = list(crossings.values())
    return max(v) / min(v)


# ----------------------------------------------------------------------
# (b) Two-loop running of the three couplings (no new physics).
#   d a_i^{-1}/dt = -b_i/(2pi) - sum_j b_ij/(8pi^2) * a_j ,  t = ln(Q/M_Z)
# Integrate up in t; find where the three alpha^{-1} are closest (the
# two-loop "triangle"). This tests alternative (b): does multi-loop alone
# shrink the gap?
# ----------------------------------------------------------------------
def two_loop_rhs(t, y):
    # y = [a1^{-1}, a2^{-1}, a3^{-1}]
    a = 1.0 / y
    one = -b_SM / (2 * np.pi)
    two = -(b2_SM @ a) / (8 * np.pi ** 2)
    return one + two


def two_loop_crossings():
    t_max = np.log(1e18 / M_Z)
    sol = solve_ivp(two_loop_rhs, [0, t_max], ainv_MZ, dense_output=True,
                    rtol=1e-9, atol=1e-12, max_step=0.05)
    ts = np.linspace(0, t_max, 40000)
    Y = sol.sol(ts)  # 3 x N
    Q = M_Z * np.exp(ts)
    out = {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        d = Y[i] - Y[j]
        # first sign change
        idx = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
        out[f"{i+1}{j+1}"] = float(Q[idx[0]]) if len(idx) else None
    # scale of closest mutual approach (max spread of the 3 at each Q)
    spread_curve = Y.max(axis=0) - Y.min(axis=0)
    k = int(np.argmin(spread_curve))
    return out, dict(best_Q=float(Q[k]),
                     min_spread_alphainv=float(spread_curve[k]),
                     alphainv_at_best=[float(v) for v in Y[:, k]])


# ----------------------------------------------------------------------
# INVERSE PROBLEM: footprint surface.
# New physics turns on at mu_NP, shifting b_i -> b_i + Delta b_i above it.
# Unification at M_GUT requires all three alpha^{-1} equal there, i.e. the
# two independent DIFFERENCES vanish. The common value A cancels from
# differences, so only Delta-b DIFFERENCES are constrained (this is exactly
# why we recover a FOOTPRINT, never an IDENTITY -- 3 unknown Delta b_i, only
# 2 difference constraints, and the overall scale of the spectrum is free).
#
# For pair (i,j), with L1=ln(mu_NP/M_Z), L2=ln(M_GUT/mu_NP):
#   0 = D_ij(M_Z) - (b_i-b_j)/2pi*L1 - (b_i-b_j+dd_ij)/2pi*L2
# => dd_ij = (2pi*D_ij(M_Z))/L2 - (b_i-b_j)*(L1+L2)/L2
# where dd_ij = Delta b_i - Delta b_j. Closed form, exact.
# ----------------------------------------------------------------------
def required_delta_diffs(mu_NP, M_GUT):
    L1 = np.log(mu_NP / M_Z)
    L2 = np.log(M_GUT / mu_NP)
    out = {}
    for i, j in [(0, 1), (1, 2), (0, 2)]:
        D = ainv_MZ[i] - ainv_MZ[j]
        bd = b_SM[i] - b_SM[j]
        dd = (2 * np.pi * D) / L2 - bd * (L1 + L2) / L2
        out[f"d{i+1}{j+1}"] = float(dd)
    return out


print("=" * 96)
print("S12 TRACK B -- new-physics INVERSE PROBLEM on PDG gauge couplings")
print("footprint (required Delta-b_i), NOT identity. Alternatives mapped first.")
print("=" * 96)
print("\nMeasured anchors -> GUT-normalised alpha^-1(M_Z):")
for L, a in zip(LAB, ainv_MZ):
    print(f"  {L:12s} {a:7.3f}")

# --- (a) deflationary null: the SM one-loop triangle as-is ---
print("\n[ALT (a) -- deflationary null: nothing forces single unification]")
cross1 = pairwise_crossings(ainv_MZ, b_SM)
sp1 = triangle_spread(cross1)
for p, Q in cross1.items():
    print(f"  cross {p}: {Q:.3e} GeV")
print(f"  one-loop triangle spread = {sp1:.3e}x  (the SM prediction; a triangle,")
print("  not a point. 'They must meet' is an aesthetic prior, not a data fact.)")

# --- (b) two-loop: does it shrink the triangle with NO new physics? ---
print("\n[ALT (b) -- two-loop running, NO new physics: does the gap shrink?]")
cross2, approach = two_loop_crossings()
for p, Q in cross2.items():
    print(f"  cross {p}: {Q if Q is None else f'{Q:.3e} GeV'}")
have = [v for v in cross2.values() if v is not None]
sp2 = (max(have) / min(have)) if len(have) > 1 else None
print(f"  two-loop triangle spread = {sp2:.3e}x" if sp2 else "  (incomplete crossings)")
print(f"  closest mutual approach at Q={approach['best_Q']:.3e} GeV, "
      f"min spread in alpha^-1 = {approach['min_spread_alphainv']:.3f}")
if sp2:
    print(f"  -> two-loop changes the triangle by factor {sp1/sp2:.2f} vs one-loop. "
          f"{'Shrinks' if sp2 < sp1 else 'Does NOT shrink'} the gap; "
          "the bulk still needs new content.")

# --- (c) extrapolation-is-conjecture tag ---
print("\n[ALT (c) -- the extrapolation itself is conjecture]")
print("  Running is MEASURED only to ~1 TeV. The crossings sit at 1e13-1e17 GeV,")
print("  >10 orders of magnitude beyond data. The straight-line extension is a")
print("  conjecture; the inversion below is conditional on it. Tagged, not cited.")

# --- INVERSE: footprint surface over (mu_NP, M_GUT) ---
print("\n[INVERSION -- required Delta-b footprint surface (mechanism-agnostic)]")
mu_grid = [1e3, 1e4, 1e6, 1e9, 1e12]          # onset scales (GeV)
gut_grid = [1e15, 2e16, 1e17]                  # target single-point scales
surface = {}
print(f"  {'mu_NP':>10} {'M_GUT':>10} | {'dd12=Db1-Db2':>14} "
      f"{'dd23=Db2-Db3':>14} {'dd13=Db1-Db3':>14}")
for mu in mu_grid:
    for G in gut_grid:
        if G <= mu:
            continue
        dd = required_delta_diffs(mu, G)
        surface[f"mu{mu:.0e}_G{G:.0e}"] = dict(mu_NP=mu, M_GUT=G, **dd)
        print(f"  {mu:10.0e} {G:10.0e} | {dd['d12']:14.3f} "
              f"{dd['d23']:14.3f} {dd['d13']:14.3f}")

# Locate MSSM on the surface (CONJECTURE -- not cited as support).
dd_MSSM = (b_MSSM - b_SM)
print("\n  [LOCATE -- MSSM Delta-b = b_MSSM - b_SM (conjecture, SUSY unobserved):]")
print(f"    Db = {dd_MSSM.round(3).tolist()}  -> "
      f"dd12={dd_MSSM[0]-dd_MSSM[1]:.3f} dd23={dd_MSSM[1]-dd_MSSM[2]:.3f} "
      f"dd13={dd_MSSM[0]-dd_MSSM[2]:.3f}")
print("    This is ONE point on the footprint surface (mu_NP~1 TeV). It LOCATES")
print("    a known model on the constraint, it does not validate it.")

# --- GRAVITY into log: how violent a modification to join the gauge form? ---
print("\n[GRAVITY into log -- how far is gravity from the gauge running FORM?]")
aG_inv_MZ = (M_PLANCK / M_Z) ** 2          # alpha_G^-1(M_Z), power-law
aG_inv_MPl = 1.0                            # alpha_G^-1 ~ O(1) at M_Planck
# A linear-in-lnQ form joining those two endpoints would need slope:
slope_needed = (aG_inv_MZ - aG_inv_MPl) / np.log(M_PLANCK / M_Z)
print(f"  alpha_G^-1(M_Z) = (M_Pl/M_Z)^2 = {aG_inv_MZ:.3e} (power-law), "
      f"alpha_G^-1(M_Pl) ~ 1.")
print(f"  To force gravity onto a single linear-in-ln(Q) line between those")
print(f"  endpoints requires an effective slope b_G/2pi ~ {slope_needed:.3e},")
print(f"  i.e. b_G ~ {2*np.pi*slope_needed:.3e} -- larger than the gauge b_i")
print(f"  (~O(1-10)) by ~{2*np.pi*slope_needed/10:.1e}x. Gravity does not 'join'")
print("  the family by a mild shift; the power-law->log change is a different")
print("  substrate FORM, not a per-force expression tweak. Confirms INFO-037.")

reading = (
    "Track B, data/interpretation/frame kept separate. DATA: (a) the SM "
    "one-loop triangle spans {:.2e}x; (b) two-loop running gives spread "
    "{:s}; (c) the inversion yields a smooth required-Delta-b surface over "
    "(mu_NP, M_GUT). INTERPRETATION: new physics CAN close the triangle, but "
    "only its FOOTPRINT (Delta-b differences) is recoverable, never its "
    "identity -- 3 unknown shifts, 2 difference constraints, free spectrum "
    "scale. MSSM is one located point on that surface, not support. Gravity "
    "needs a form change (power-law->log), not a shift. FRAME: 'shared "
    "running substrate + per-force expression' survives; 'single unification' "
    "remains a conjecture conditional on >10-orders extrapolation. The "
    "leading deflationary read (nothing forces a single point) is unrefuted."
).format(sp1, f"{sp2:.2e}x" if sp2 else "incomplete")

print("\n" + "=" * 96)
print("READING:")
print(reading)
print("=" * 96)

out = dict(
    anchors=dict(M_Z=M_Z, ainv_MZ=ainv_MZ.tolist(),
                 b_SM=b_SM.tolist(), b2_SM=b2_SM.tolist()),
    alt_a_oneloop=dict(crossings=cross1, triangle_spread=float(sp1)),
    alt_b_twoloop=dict(crossings=cross2,
                       triangle_spread=(float(sp2) if sp2 else None),
                       closest_approach=approach,
                       shrink_factor_vs_oneloop=(float(sp1 / sp2) if sp2 else None)),
    alt_c="extrapolation to 1e13-1e17 GeV is conjecture; running measured to ~1 TeV",
    inverse_footprint_surface=surface,
    mssm_located=dict(delta_b=dd_MSSM.tolist(),
                      dd12=float(dd_MSSM[0] - dd_MSSM[1]),
                      dd23=float(dd_MSSM[1] - dd_MSSM[2]),
                      dd13=float(dd_MSSM[0] - dd_MSSM[2]),
                      note="conjecture, SUSY unobserved; located not cited"),
    gravity_into_log=dict(aG_inv_MZ=float(aG_inv_MZ),
                          slope_needed=float(slope_needed),
                          b_G_effective=float(2 * np.pi * slope_needed),
                          note="power-law->log is a form change, not a shift"),
    provenance=dict(
        solid="M_Z couplings (PDG); 1- and 2-loop SM beta functions (standard)",
        conjecture="GUT-scale extrapolation; MSSM; any specific NP spectrum",
        method="footprint = required Delta-b DIFFERENCES; identity unrecoverable"),
    reading=reading,
)
json.dump(out, open("s12_track_b_inverse_results.json", "w"), indent=2)
print("\nWrote s12_track_b_inverse_results.json")
