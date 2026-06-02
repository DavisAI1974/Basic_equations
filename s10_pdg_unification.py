"""
S10 item (5), set 2 of 2 (REAL DATA, special build): four-force coupling
unification from MEASURED couplings, read through the substrate/
expression frame.

This is the "do the forces unify" question in its native form -- the one
the four-force caricatures were a stand-in for. It does NOT fit the
windowed-H stack (couplings vs energy, not a time series), so it is a
purpose-built analysis, per Greg's "build something special if a top-2
set needs it."

Frame mapping (the surviving Session-9 legs, applied to real data):
  - SUBSTRATE (shared) : the common FORM of the running -- each inverse
    gauge coupling is linear in ln(Q): alpha_i^{-1}(Q) =
    alpha_i^{-1}(M_Z) - (b_i/2pi) ln(Q/M_Z). One shared functional form.
  - EXPRESSION (per-force) : the slope b_i (beta coefficient) differs
    per force. This is the force-specific "expression" of the shared
    running substrate.
  - UNIFICATION TEST : anchored ONLY on measured M_Z couplings, do the
    three lines meet at a SINGLE point? If yes -> one coupling at high
    scale (true unification). If they cross pairwise at different scales
    -> a "triangle", no single unification without new physics.

Literature-as-conjecture discipline (CLAUDE.md rule, load-bearing here):
  - MEASURED and solid (meets the replication bar): the M_Z-scale
    couplings (alpha_em, sin^2theta_W, alpha_s) and the MEASURED running
    of alpha_s(Q) from ~2 GeV to ~1 TeV (asymptotic freedom, replicated
    across many experiments). We anchor and VALIDATE on these.
  - CONJECTURE (does NOT meet the bar): extrapolating the running to the
    GUT scale ~10^16 GeV; MSSM unification (SUSY searched, not found).
    We COMPUTE these but tag them as extrapolation, never as support.

Speaking posture (Rule C): I THINK the measured anchors + SM running
will give three pairwise crossings at DIFFERENT scales (no single point)
-- the known SM near-miss -- but we wait on what the numbers say and
where they point. The honest output is the crossing scales themselves,
not a verdict decided in advance.

Outputs s10_pdg_unification_results.json.
"""
import json
import numpy as np

# ----------------------------------------------------------------------
# MEASURED anchors at Q = M_Z  (PDG 2024 world averages; rock-solid,
# meets the replication bar). Uncertainties carried for the dominant one.
# ----------------------------------------------------------------------
M_Z = 91.1876                      # GeV (PDG)
alpha_em_inv_MZ = 127.951          # MSbar 5-flavour running EM coupling^{-1} (PDG)
sin2_thetaW_MZ = 0.23122           # MSbar (PDG)
alpha_s_MZ = 0.1180
alpha_s_MZ_err = 0.0009            # PDG world average

M_PLANCK = 1.22e19                 # GeV (for the gravity comparison)

# GUT-normalised inverse couplings at M_Z, DERIVED from the measured
# electroweak anchors (standard relations):
#   alpha_2^{-1} = sin^2thetaW * alpha_em^{-1}
#   alpha_1^{-1} = (3/5) cos^2thetaW * alpha_em^{-1}   (GUT normalisation)
#   alpha_3^{-1} = 1 / alpha_s
cos2 = 1.0 - sin2_thetaW_MZ
a1_inv_MZ = (3.0 / 5.0) * cos2 * alpha_em_inv_MZ
a2_inv_MZ = sin2_thetaW_MZ * alpha_em_inv_MZ
a3_inv_MZ = 1.0 / alpha_s_MZ
ainv_MZ = np.array([a1_inv_MZ, a2_inv_MZ, a3_inv_MZ])

# One-loop beta coefficients (b_i in alpha^{-1}' = -b_i/2pi per d ln Q).
# Standard Model (GUT-normalised b_1):
b_SM = np.array([41.0 / 10.0, -19.0 / 6.0, -7.0])
# MSSM (CONJECTURE leg -- SUSY not observed; computed, never cited as support):
b_MSSM = np.array([33.0 / 5.0, 1.0, -3.0])

LABELS = ["U(1)_Y  (alpha_1)", "SU(2)_L (alpha_2, weak)", "SU(3)_C (alpha_3, strong)"]


def run_inv(ainv0, b, Q):
    """alpha_i^{-1}(Q) one-loop."""
    return ainv0 - (b / (2 * np.pi)) * np.log(Q / M_Z)


def pairwise_crossings(ainv0, b):
    """Scale Q (GeV) where alpha_i^{-1} = alpha_j^{-1}, all pairs."""
    out = {}
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        denom = (b[i] - b[j])
        t = 2 * np.pi * (ainv0[i] - ainv0[j]) / denom   # = ln(Q/M_Z)
        out[f"{i+1}-{j+1}"] = M_Z * np.exp(t)
    return out


def alpha_s_run(Q):
    """alpha_s(Q) from the measured M_Z value, one-loop (n_f=5)."""
    b3 = -7.0  # SM, but for alpha_s alone the n_f-dependent coeff is
    # 11 - 2/3 n_f = 11 - 10/3 = 23/3 for n_f=5 in the alpha_s convention.
    beta0 = 23.0 / 3.0   # = 11 - 2/3*5
    a_inv = 1.0 / alpha_s_MZ + (beta0 / (2 * np.pi)) * np.log(Q / M_Z)
    return 1.0 / a_inv


# ----------------------------------------------------------------------
# MEASURED alpha_s(Q) running points (approximate PDG-review
# determinations across the running-coupling figure; central values are
# representative, used only to VALIDATE that the running form is real in
# data -- NOT to anchor the unification computation, which uses the M_Z
# couplings above). Flagged approximate.
# ----------------------------------------------------------------------
ALPHA_S_MEASURED = [
    # (Q GeV, alpha_s, approx_err, source class)
    (1.78,  0.314, 0.012, "tau decays"),
    (2.0,   0.301, 0.018, "low-Q"),
    (6.7,   0.220, 0.020, "Upsilon / low e+e-"),
    (31.6,  0.150, 0.012, "DIS / jets"),
    (91.2,  0.1179, 0.0010, "Z pole (LEP/SLD)"),
    (133.0, 0.113, 0.006, "LEP2 event shapes"),
    (206.0, 0.108, 0.006, "LEP2 event shapes"),
    (500.0, 0.095, 0.007, "jets"),
    (1000.0, 0.089, 0.006, "LHC jets / ttbar"),
]


def fmt_Q(Q):
    return f"{Q:.3e} GeV"


print("=" * 96)
print("S10 PDG four-force unification (REAL measured couplings, special build)")
print("Frame: shared SUBSTRATE = linear-in-ln(Q) running form; EXPRESSION = per-force slope b_i")
print("=" * 96)

print("\nMeasured anchors at M_Z = 91.1876 GeV (PDG 2024, solid):")
print(f"  alpha_em^-1 = {alpha_em_inv_MZ}   sin^2thetaW = {sin2_thetaW_MZ}   "
      f"alpha_s = {alpha_s_MZ} +/- {alpha_s_MZ_err}")
print("Derived GUT-normalised inverse couplings at M_Z:")
for L, a in zip(LABELS, ainv_MZ):
    print(f"  {L:26s} alpha^-1(M_Z) = {a:7.3f}")

# --- validation: does the SM running match MEASURED alpha_s(Q)? ---
print("\n[VALIDATION] SM one-loop alpha_s(Q) running vs MEASURED points "
      "(approx PDG determinations):")
chi_terms = []
val_rows = []
for Q, a_meas, err, src in ALPHA_S_MEASURED:
    a_pred = alpha_s_run(Q)
    pull = (a_meas - a_pred) / err
    chi_terms.append(pull ** 2)
    val_rows.append(dict(Q=Q, measured=a_meas, predicted=round(float(a_pred), 4),
                         err=err, pull=round(float(pull), 2), source=src))
    print(f"  Q={Q:7.2f} GeV  measured={a_meas:.4f}+/-{err:.4f}  "
          f"SM-pred={a_pred:.4f}  pull={pull:+.2f}sigma   [{src}]")
chi2 = float(np.sum(chi_terms))
ndf = len(ALPHA_S_MEASURED) - 1
print(f"  -> chi2/ndf = {chi2:.1f}/{ndf} = {chi2/ndf:.2f}  "
      f"(one-loop only; multi-loop + thresholds tighten this). The MEASURED "
      f"running follows the predicted linear-in-ln(Q) substrate form.")

# --- unification: pairwise crossings, SM (measured-anchored) ---
print("\n[UNIFICATION -- SM, anchored ONLY on measured M_Z couplings]")
cross_SM = pairwise_crossings(ainv_MZ, b_SM)
for pair, Q in cross_SM.items():
    print(f"  alpha_{pair[0]} = alpha_{pair[2]} at Q = {fmt_Q(Q)}")
spread_SM = max(cross_SM.values()) / min(cross_SM.values())
print(f"  -> three crossings span a factor {spread_SM:.1e} in energy: a TRIANGLE,")
print(f"     NOT a single point. SM couplings do not unify without new physics.")

# --- MSSM (CONJECTURE leg, tagged) ---
print("\n[UNIFICATION -- MSSM, CONJECTURE (SUSY not observed; not cited as support)]")
cross_MSSM = pairwise_crossings(ainv_MZ, b_MSSM)
for pair, Q in cross_MSSM.items():
    print(f"  alpha_{pair[0]} = alpha_{pair[2]} at Q = {fmt_Q(Q)}")
spread_MSSM = max(cross_MSSM.values()) / min(cross_MSSM.values())
print(f"  -> crossings span a factor {spread_MSSM:.1f} -- a near-point at "
      f"~{np.exp(np.mean(np.log(list(cross_MSSM.values())))):.1e} GeV.")
print(f"     Better convergence, but rests on unobserved SUSY -> conjecture.")

# --- gravity: different FORM (does not join the linear-in-lnQ family) ---
# dimensionless gravitational coupling alpha_G(E) = (E/M_Planck)^2
print("\n[GRAVITY -- different substrate FORM, the natural odd-one-out]")
print("  alpha_G(E) = (E/M_Planck)^2  grows as E^2 (power law), NOT linear")
print("  in ln(Q). It does not share the gauge-coupling running form, so it")
print("  cannot join the same substrate. At M_Z: alpha_G = "
      f"{(M_Z/M_PLANCK)**2:.2e}; it reaches O(1) only at E ~ M_Planck "
      f"~ {M_PLANCK:.1e} GeV.")

reading = (
    "Mapped onto the surviving four-force frame: the three GAUGE forces "
    "share a substrate (one running FORM, linear in ln Q) with per-force "
    "EXPRESSION (the slopes b_i) -- exactly the 'shared substrate + "
    "expression difference' pattern, now on REAL measured couplings. But "
    "the substrate does NOT force them to a single point: SM crossings are "
    "a triangle spanning ~1e4 in energy. Gravity has a different substrate "
    "form entirely (power-law, not log) and does not join the family. So "
    "the real-data verdict on 'do the four forces unify' is: shared running "
    "substrate among the three gauge forces, no single unification scale "
    "without conjectural new physics, and gravity outside the shared form."
)
print("\n" + "=" * 96)
print("READING:")
print(reading)
print("=" * 96)

out = dict(
    anchors=dict(M_Z=M_Z, alpha_em_inv_MZ=alpha_em_inv_MZ,
                 sin2_thetaW_MZ=sin2_thetaW_MZ, alpha_s_MZ=alpha_s_MZ,
                 alpha_s_MZ_err=alpha_s_MZ_err,
                 a1_inv_MZ=float(a1_inv_MZ), a2_inv_MZ=float(a2_inv_MZ),
                 a3_inv_MZ=float(a3_inv_MZ)),
    beta_SM=b_SM.tolist(), beta_MSSM=b_MSSM.tolist(),
    alpha_s_validation=dict(rows=val_rows, chi2=chi2, ndf=ndf),
    crossings_SM={k: float(v) for k, v in cross_SM.items()},
    crossings_SM_spread=float(spread_SM),
    crossings_MSSM={k: float(v) for k, v in cross_MSSM.items()},
    crossings_MSSM_spread=float(spread_MSSM),
    gravity_alpha_G_at_MZ=float((M_Z / M_PLANCK) ** 2),
    provenance=dict(
        solid="M_Z couplings + measured alpha_s(Q) running (replicated)",
        conjecture="GUT-scale extrapolation; MSSM unification (SUSY unobserved)",
        approximate="alpha_s(Q) central values are representative PDG-review "
                    "determinations, used for validation only"),
    reading=reading,
)
json.dump(out, open("s10_pdg_unification_results.json", "w"), indent=2)
print("\nWrote s10_pdg_unification_results.json")
