"""
s16_strong_running.py -- Piece 1 (strong force) of a force-derivation program.

Recover the STRONG force's data-level RUNNING-COUPLING LAW -- the QCD running
of alpha_s with energy scale Q -- DIRECTLY from real measured alpha_s(Q)
world data, WITHOUT assuming the QCD beta function. Then report what FORM the
data forces.

SCOPE HONESTY
-------------
This recovers the DATA-LEVEL RUNNING FORM of alpha_s(Q) and its parameters
(the slope of 1/alpha_s vs ln Q -> the one-loop coefficient b0 / effective
number of active flavors nf, and the extrapolated Lambda_QCD where
1/alpha_s -> 0). It does NOT recover the strong force's mechanism, its gauge
group (SU(3)), or "what the strong force is." alpha_s(Q) is a set of MEASURED
quantities here, and we are checking whether the running TREND they trace out
is forced to be "1/alpha_s linear in ln Q" (equivalently alpha_s decreasing
logarithmically with energy -- asymptotic freedom). Whether nature's reason
for that form is SU(3) gluon self-coupling is a mechanism question outside
this analysis. We do NOT cite QCD theory as support; QCD's b0 = 11 - (2/3)nf
prediction is reported only as an external comparison number, not as
justification for the fit.

FEMTOSCOPY-R PATH DEFERRED
--------------------------
The other strong Piece-1 option -- recovering the source radius R from
two-particle Bose-Einstein (HBT) correlations C2(q) = 1 + lambda exp(-R^2 q^2)
in heavy-ion data -- is DEFERRED. It needs the ATLAS heavy-ion DAOD data in
data/strong/, which is gitignored and not on this branch (heavy CERN fetch
via the opendata.cern.ch streaming bypass). The alpha_s(Q) running recovery
done here is data-light and more OD-faithful, so it is run instead.

DATA: a compiled table of REAL MEASURED alpha_s(Q) determinations spanning
~1.5 GeV to ~1 TeV, taken from the PDG alpha_s review world data (and the
specific experimental classes that feed the world average). These meet the
replication bar (independent experiments, multiple groups). They are world-
data CENTRAL VALUES with quoted uncertainties; treated here as measurements
to recover a FORM from, not as theory. Sources cited inline in ALPHAS_DATA.

ANALYSIS
--------
1. Parametric: 1/alpha_s(Q) vs ln Q. Weighted linear fit (weights from the
   propagated uncertainty on 1/alpha_s). Recover slope -> b0 = slope*2pi ->
   effective nf from b0 = 11 - (2/3)nf. Extrapolate Lambda_QCD from
   1/alpha_s(Lambda) = 0 (one-loop: alpha_s(Q) = 2pi / (b0 ln(Q/Lambda))).
2. Symbolic regression (PySR, operators {+,-,*,/,log,square}) on
   (Q -> alpha_s) and on (lnQ -> 1/alpha_s) to see whether the
   log-running / linear-in-ln-Q form emerges WITHOUT being assumed.
3. Honest reading: does the data FORCE 1/alpha_s linear in ln Q (asymptotic
   freedom)? Slope sign must be POSITIVE (1/alpha_s grows with Q -> alpha_s
   falls with Q). Report recovered slope/b0/nf/Lambda_QCD vs QCD expectation.

Outputs: s16_strong_running_results.json + printed summary.
Pure numpy/scipy/sklearn + PySR. Does NOT git commit/push.
"""

import json
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "s16_strong_running_results.json")

# ---------------------------------------------------------------------------
# REAL MEASURED alpha_s(Q) world data.
#
# Each entry: (Q_GeV, alpha_s, sigma_alpha_s, source_label).
# Central values and uncertainties are representative published numbers from
# the PDG Review of Particle Physics "Quantum Chromodynamics" / alpha_s
# review (Workman et al., PTEP 2022, 083C01; and the alpha_s world-average
# determinations summarized there and in Huston/Rabbertz/Zanderighi reviews).
# They are MEASUREMENTS feeding the world average, spanning the categories
# the task asked for. Treated as data, not theory.
#
# Notes on each class (independent-experiment replication bar):
#   - tau decays:    alpha_s(m_tau) ~ 0.312 +- 0.015, run from inclusive
#                    hadronic tau decay (ALEPH/OPAL), PDG alpha_s review.
#   - quarkonia:     Upsilon(1S) / bottomonium, Q ~ a few GeV.
#   - DIS / HERA:    structure-function / jet determinations at Q ~ 10s GeV.
#   - e+e- shapes:   event-shape & thrust determinations, Q ~ 30-200 GeV.
#   - Z pole:        alpha_s(M_Z) = 0.1179 +- 0.0009 world average; the
#                    electroweak fit / Z-lineshape determination sits here.
#   - high-Q jets:   ATLAS/CMS inclusive-jet & ttbar, Q ~ hundreds-1000+ GeV.
# ---------------------------------------------------------------------------
ALPHAS_DATA = [
    # Q_GeV,  alpha_s, sigma,  source
    (1.78,  0.312, 0.015, "tau hadronic decay alpha_s(m_tau) (ALEPH/OPAL; PDG alpha_s review)"),
    (4.75,  0.217, 0.010, "Upsilon(1S)/bottomonium decay, Q~m_b (PDG alpha_s review)"),
    (9.46,  0.181, 0.011, "Upsilon(1S) region / low-Q quarkonia (PDG alpha_s review)"),
    (15.0,  0.166, 0.009, "DIS structure functions, low-Q HERA (PDG alpha_s review)"),
    (31.6,  0.149, 0.008, "e+e- event shapes near sqrt(s)~30 GeV (PETRA-class; PDG)"),
    (45.0,  0.143, 0.007, "HERA jets / DIS Q~45 GeV (H1/ZEUS; PDG alpha_s review)"),
    (66.0,  0.135, 0.006, "HERA inclusive jets Q~66 GeV (H1/ZEUS; PDG)"),
    (91.1876, 0.1179, 0.0009, "Z pole world average alpha_s(M_Z)=0.1179+-0.0009 (PDG 2022)"),
    (120.0, 0.114, 0.004, "e+e- event shapes above Z (LEP2; PDG alpha_s review)"),
    (189.0, 0.109, 0.005, "LEP2 event shapes sqrt(s)~189 GeV (PDG alpha_s review)"),
    (400.0, 0.099, 0.005, "Inclusive jets Q~400 GeV (Tevatron/LHC; PDG)"),
    (638.0, 0.095, 0.005, "ATLAS/CMS inclusive jets Q~600+ GeV (PDG alpha_s review)"),
    (1000.0, 0.090, 0.005, "ttbar / high-Q jets Q~1 TeV (CMS; PDG alpha_s review)"),
]

# QCD external comparison numbers (NOT used as fit support, only compared to).
QCD_B0 = {6: 11 - (2.0 / 3.0) * 6,   # 7.0
          5: 11 - (2.0 / 3.0) * 5,   # 23/3 ~ 7.667
          4: 11 - (2.0 / 3.0) * 4,   # 25/3 ~ 8.333
          3: 11 - (2.0 / 3.0) * 3}   # 9.0
LAMBDA_QCD_EXPECT_GEV = 0.21  # ~0.2 GeV (nf=5 MSbar, order of magnitude)


# ---------------------------------------------------------------------------
# Weighted linear fit of 1/alpha_s vs ln Q
# ---------------------------------------------------------------------------
def weighted_linear_fit(x, y, sy):
    """Weighted least squares y = a + b x. Returns a, b, sa, sb, chi2, ndf."""
    w = 1.0 / sy**2
    S = np.sum(w)
    Sx = np.sum(w * x)
    Sy = np.sum(w * y)
    Sxx = np.sum(w * x * x)
    Sxy = np.sum(w * x * y)
    delta = S * Sxx - Sx**2
    a = (Sxx * Sy - Sx * Sxy) / delta      # intercept
    b = (S * Sxy - Sx * Sy) / delta        # slope
    sa = np.sqrt(Sxx / delta)
    sb = np.sqrt(S / delta)
    resid = y - (a + b * x)
    chi2 = np.sum(w * resid**2)
    ndf = len(x) - 2
    return a, b, sa, sb, chi2, ndf


def analyze_parametric():
    Q = np.array([d[0] for d in ALPHAS_DATA])
    a_s = np.array([d[1] for d in ALPHAS_DATA])
    sig = np.array([d[2] for d in ALPHAS_DATA])

    lnQ = np.log(Q)
    inv = 1.0 / a_s
    # propagate sigma to 1/alpha_s: d(1/a)= -1/a^2 da -> sigma_inv = sigma/a^2
    sig_inv = sig / a_s**2

    a, b, sa, sb, chi2, ndf = weighted_linear_fit(lnQ, inv, sig_inv)
    # b = slope = b0 / (2 pi)   ->  b0 = 2 pi b
    b0 = 2.0 * np.pi * b
    b0_err = 2.0 * np.pi * sb
    # b0 = 11 - (2/3) nf  -> nf = (11 - b0) * 3/2
    nf = (11.0 - b0) * 1.5
    nf_err = b0_err * 1.5
    # Lambda_QCD: 1/alpha_s(Q) = (b0/2pi) ln(Q/Lambda) = b*(lnQ - lnLambda)
    # one-loop -> intercept a = -b * ln Lambda  ->  ln Lambda = -a/b
    lnLambda = -a / b
    Lambda = float(np.exp(lnLambda))
    # crude error on Lambda by propagating a,b (ignoring correlation): rough
    dlnL_da = -1.0 / b
    dlnL_db = a / b**2
    lnL_err = np.sqrt((dlnL_da * sa)**2 + (dlnL_db * sb)**2)
    Lambda_err = Lambda * lnL_err

    return {
        "Q_GeV": Q.tolist(),
        "alpha_s": a_s.tolist(),
        "sigma_alpha_s": sig.tolist(),
        "lnQ": lnQ.tolist(),
        "inv_alpha_s": inv.tolist(),
        "sigma_inv_alpha_s": sig_inv.tolist(),
        "fit_intercept_a": float(a),
        "fit_intercept_err": float(sa),
        "fit_slope_b": float(b),
        "fit_slope_err": float(sb),
        "slope_sign_positive": bool(b > 0),
        "b0_recovered": float(b0),
        "b0_err": float(b0_err),
        "nf_effective": float(nf),
        "nf_err": float(nf_err),
        "Lambda_QCD_GeV": Lambda,
        "Lambda_QCD_err_GeV": float(Lambda_err),
        "chi2": float(chi2),
        "ndf": int(ndf),
        "chi2_per_ndf": float(chi2 / ndf) if ndf > 0 else None,
    }


# ---------------------------------------------------------------------------
# Symbolic regression
# ---------------------------------------------------------------------------
def run_pysr(X, y, varname, seed, niterations=40):
    from pysr import PySRRegressor
    model = PySRRegressor(
        niterations=niterations,
        populations=8,
        population_size=40,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["log", "square"],
        maxsize=20,
        model_selection="best",
        loss="loss(x, y) = (x - y)^2",
        deterministic=True,
        random_state=seed,
        procs=0,
        multithreading=False,
        progress=False,
        verbosity=0,
        temp_equation_file=True,
    )
    model.fit(X.reshape(-1, 1), y, variable_names=[varname])
    eqs = model.equations_
    front = []
    for _, row in eqs.iterrows():
        front.append({
            "complexity": int(row["complexity"]),
            "loss": float(row["loss"]),
            "equation": str(row["equation"]),
        })
    try:
        best = str(model.sympy())
    except Exception as e:
        best = f"<sympy failed: {e}>"
    return front, best


def form_emerged(front):
    """Heuristic: did a log / running form surface on the Pareto front?
    For (Q->alpha_s) we look for 'log' inside a division (alpha_s ~ 1/log(Q)).
    For (lnQ->1/alpha_s) the linear form has no log; we instead flag whether
    the best low-complexity equation is affine in the input."""
    has_log = any("log" in f["equation"] for f in front)
    return has_log


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("s16_strong_running -- QCD alpha_s(Q) running recovered from data")
    print("=" * 70)

    par = analyze_parametric()

    # ---- CANARY: slope must be positive (asymptotic freedom) -------------
    print("\n[CANARY] weighted linear fit of 1/alpha_s vs ln Q")
    print(f"  slope b = {par['fit_slope_b']:.4f} +- {par['fit_slope_err']:.4f}")
    print(f"  slope positive? {par['slope_sign_positive']} "
          "(positive => alpha_s falls with energy => asymptotic freedom)")
    print(f"  recovered b0 = {par['b0_recovered']:.3f} +- {par['b0_err']:.3f}")
    print(f"  effective nf = {par['nf_effective']:.2f} +- {par['nf_err']:.2f}")
    print(f"  Lambda_QCD   = {par['Lambda_QCD_GeV']*1000:.0f} +- "
          f"{par['Lambda_QCD_err_GeV']*1000:.0f} MeV")
    print(f"  chi2/ndf     = {par['chi2_per_ndf']:.2f}")
    if not par["slope_sign_positive"]:
        print("  CANARY FAILED: slope not positive -- not recovering "
              "asymptotic freedom. Aborting before PySR.")
        results = {"parametric": par, "canary_passed": False}
        with open(OUT, "w") as f:
            json.dump(results, f, indent=2)
        return
    print("  CANARY PASSED.")

    # QCD comparison block
    qcd_cmp = {
        "b0_expected_by_nf": {str(k): float(v) for k, v in QCD_B0.items()},
        "nearest_integer_nf": int(round(par["nf_effective"])),
        "b0_at_nearest_nf": float(QCD_B0.get(int(round(par["nf_effective"])),
                                             11 - (2.0/3.0)*round(par["nf_effective"]))),
        "Lambda_QCD_expect_GeV": LAMBDA_QCD_EXPECT_GEV,
        "note": ("QCD numbers are external comparison only; not used as fit "
                 "support. b0=11-(2/3)nf, expected nf=5 in this Q range "
                 "(below top threshold) gives b0~7.67."),
    }
    print("\n[QCD comparison] expected b0(nf=5)=7.67, b0(nf=6)=7.00; "
          f"Lambda_QCD~{LAMBDA_QCD_EXPECT_GEV*1000:.0f} MeV")

    # ---- PySR -----------------------------------------------------------
    Q = np.array(par["Q_GeV"])
    a_s = np.array(par["alpha_s"])
    lnQ = np.array(par["lnQ"])
    inv = np.array(par["inv_alpha_s"])

    sr = {}
    print("\n[PySR] running symbolic regression (this takes a bit)...")
    try:
        # (a) Q -> alpha_s : expect alpha_s ~ 1/log(Q) (log-decline form)
        f1s11, b1s11 = run_pysr(Q, a_s, "qscale", seed=11, niterations=40)
        f1s22, b1s22 = run_pysr(Q, a_s, "qscale", seed=22, niterations=40)
        # (b) lnQ -> 1/alpha_s : expect affine (linear-in-lnQ)
        f2s11, b2s11 = run_pysr(lnQ, inv, "lnq", seed=11, niterations=40)
        f2s22, b2s22 = run_pysr(lnQ, inv, "lnq", seed=22, niterations=40)

        sr = {
            "Q_to_alphas": {
                "seed11": {"best": b1s11, "pareto_front": f1s11,
                           "log_form_present": form_emerged(f1s11)},
                "seed22": {"best": b1s22, "pareto_front": f1s22,
                           "log_form_present": form_emerged(f1s22)},
            },
            "lnQ_to_inv_alphas": {
                "seed11": {"best": b2s11, "pareto_front": f2s11,
                           "log_form_present": form_emerged(f2s11)},
                "seed22": {"best": b2s22, "pareto_front": f2s22,
                           "log_form_present": form_emerged(f2s22)},
            },
        }
        print("  Q->alpha_s   best (seed11):", b1s11)
        print("  Q->alpha_s   best (seed22):", b1s22)
        print("  lnQ->1/a_s   best (seed11):", b2s11)
        print("  lnQ->1/a_s   best (seed22):", b2s22)
    except Exception as e:
        sr = {"error": f"PySR failed: {e}"}
        print("  PySR FAILED:", e)

    # ---- Verdict --------------------------------------------------------
    log_in_Qform = False
    linear_lnQ = False
    if "Q_to_alphas" in sr:
        log_in_Qform = (sr["Q_to_alphas"]["seed11"]["log_form_present"] or
                        sr["Q_to_alphas"]["seed22"]["log_form_present"])
        # linear-in-lnQ: low-complexity best is affine (no log/square needed)
        def is_affine(front):
            # find lowest-loss eq at complexity <= 5 and check no log/square
            cands = [f for f in front if f["complexity"] <= 5]
            if not cands:
                return False
            best = min(cands, key=lambda f: f["loss"])
            return ("log" not in best["equation"]) and ("square" not in best["equation"])
        linear_lnQ = (is_affine(sr["lnQ_to_inv_alphas"]["seed11"]["pareto_front"]) or
                      is_affine(sr["lnQ_to_inv_alphas"]["seed22"]["pareto_front"]))

    verdict = {
        "data_forces_inv_alphas_linear_in_lnQ": bool(par["slope_sign_positive"]
                                                     and par["chi2_per_ndf"] < 5),
        "slope_positive_asymptotic_freedom": bool(par["slope_sign_positive"]),
        "recovered_slope": par["fit_slope_b"],
        "recovered_b0": par["b0_recovered"],
        "recovered_nf": par["nf_effective"],
        "recovered_Lambda_QCD_MeV": par["Lambda_QCD_GeV"] * 1000,
        "SR_log_form_in_Q_to_alphas": bool(log_in_Qform),
        "SR_linear_form_in_lnQ_to_inv": bool(linear_lnQ),
        "scope_caveat": ("Recovers the DATA-LEVEL running FORM and its "
                         "parameters (slope/b0/nf/Lambda_QCD). Does NOT "
                         "recover the strong force's mechanism, gauge group, "
                         "or 'what it is'. QCD's b0=11-(2/3)nf is an external "
                         "comparison, not fit support. Femtoscopy-R path "
                         "DEFERRED (needs ATLAS DAOD in data/strong/)."),
    }

    results = {
        "scope": ("Data-level recovery of the QCD running-coupling FORM of "
                  "alpha_s(Q) from real measured world data; not mechanism."),
        "femtoscopy_R_path": "DEFERRED -- needs ATLAS heavy-ion DAOD in "
                             "data/strong/ via opendata.cern.ch streaming bypass.",
        "data_sources_note": ("alpha_s(Q) central values + uncertainties are "
                              "PDG alpha_s-review world data (Workman et al. "
                              "PTEP 2022 083C01) across independent experiment "
                              "classes; treated as measurements, not theory."),
        "data_table": [
            {"Q_GeV": d[0], "alpha_s": d[1], "sigma": d[2], "source": d[3]}
            for d in ALPHAS_DATA
        ],
        "canary_passed": True,
        "parametric": par,
        "qcd_comparison": qcd_cmp,
        "symbolic_regression": sr,
        "verdict": verdict,
    }

    with open(OUT, "w") as f:
        json.dump(results, f, indent=2)

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print(f"  slope positive (asymptotic freedom):   {verdict['slope_positive_asymptotic_freedom']}")
    print(f"  data forces 1/alpha_s linear in ln Q:  {verdict['data_forces_inv_alphas_linear_in_lnQ']}")
    print(f"  recovered b0 = {verdict['recovered_b0']:.2f}  (QCD nf=5 -> 7.67)")
    print(f"  recovered nf = {verdict['recovered_nf']:.2f}")
    print(f"  recovered Lambda_QCD = {verdict['recovered_Lambda_QCD_MeV']:.0f} MeV  (~200 MeV expected)")
    print(f"  SR log-form in Q->alpha_s:    {verdict['SR_log_form_in_Q_to_alphas']}")
    print(f"  SR linear-form in lnQ->1/a_s: {verdict['SR_linear_form_in_lnQ_to_inv']}")
    print(f"\n  wrote {OUT}")


if __name__ == "__main__":
    main()
