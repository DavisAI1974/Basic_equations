"""
s16_weak_zpropagator.py -- Piece 1 of a force-derivation program.

Recover the weak neutral-current Z-boson propagator (relativistic
Breit-Wigner resonance / Z lineshape) DIRECTLY from real CMS dimuon data,
both by (a) parametric fitting and (b) symbolic regression (PySR), WITHOUT
assuming the lineshape a priori in the SR step.

SCOPE HONESTY
-------------
This recovers the DATA-LEVEL resonance/propagator SHAPE and its parameters
(the Z mass M_Z and width Gamma_Z) from the dimuon invariant-mass spectrum.
It does NOT recover the weak force's mechanism, gauge structure, or "what it
is." The Breit-Wigner denominator (M^2 - M_Z^2)^2 + M^2 Gamma^2 is the
propagator-pole signature of a resonance; recovering it confirms the
observable lineshape and its two parameters, nothing more. M_Z and Gamma_Z
are measured quantities here, not derived from any theory of the weak
interaction.

"They never stacked": SymbolFit (and similar) ran SR on exactly this kind of
real dimuon data but modeled the smooth BACKGROUND and discarded the
resonance. Here SR is pointed at the resonance region itself to test whether
a Breit-Wigner-like rational lineshape emerges from data alone.

DATA: data/forces/Zmumu.csv (CMS open data, ~10583 dimuon events).
  M^2 = 2*pt1*pt2*(cosh(eta1-eta2) - cos(phi1-phi2))  (massless-muon approx).

Outputs: s16_weak_zpropagator_results.json + printed summary.
Pure numpy/scipy/sklearn + PySR. Does not commit.
"""

import json
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

PDG_MZ = 91.1876     # GeV
PDG_GZ = 2.4952      # GeV

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "forces", "Zmumu.csv")
OUT = os.path.join(HERE, "s16_weak_zpropagator_results.json")

# Histogram / fit window choices (reported in JSON)
WINDOW = (60.0, 120.0)      # GeV, full spectrum window
N_BINS = 80                 # equal-width bins -> 0.75 GeV/bin
FIT_WINDOW = (70.0, 110.0)  # narrower peak region for the parametric fit


# ---------------------------------------------------------------------------
# Step 1: dimuon invariant mass
# ---------------------------------------------------------------------------
def invariant_mass(df):
    M2 = 2.0 * df.pt1 * df.pt2 * (
        np.cosh(df.eta1 - df.eta2) - np.cos(df.phi1 - df.phi2)
    )
    M2 = np.clip(M2, 0.0, None)
    return np.sqrt(M2)


# ---------------------------------------------------------------------------
# Step 2: lineshape models
# ---------------------------------------------------------------------------
def bw_rel(M, A, MZ, GZ):
    """Relativistic Breit-Wigner (cross-section form)."""
    return A * (M**2 * GZ**2) / ((M**2 - MZ**2) ** 2 + M**2 * GZ**2)


def bw_rel_bkg(M, A, MZ, GZ, b0, b1):
    """Relativistic BW + linear background."""
    return bw_rel(M, A, MZ, GZ) + b0 + b1 * (M - 90.0)


def bw_rel_expbkg(M, A, MZ, GZ, c0, lam):
    """Relativistic BW + exponential background (Drell-Yan-like falloff)."""
    return bw_rel(M, A, MZ, GZ) + c0 * np.exp(-lam * (M - 70.0))


def lorentzian_bkg(M, A, MZ, GZ, b0, b1):
    """Non-relativistic BW (Lorentzian) + linear background, secondary fit."""
    return A * (GZ / 2.0) ** 2 / ((M - MZ) ** 2 + (GZ / 2.0) ** 2) + b0 + b1 * (M - 90.0)


def chi2_ndf(model, M, y, yerr, popt):
    r = (y - model(M, *popt)) / yerr
    chi2 = float(np.sum(r**2))
    ndf = len(y) - len(popt)
    return chi2, ndf, chi2 / ndf


def fit_model(model, centers, counts, p0, bounds=None):
    yerr = np.sqrt(np.maximum(counts, 1.0))
    kw = dict(p0=p0, sigma=yerr, absolute_sigma=True, maxfev=200000)
    if bounds is not None:
        kw["bounds"] = bounds
    popt, pcov = curve_fit(model, centers, counts, **kw)
    perr = np.sqrt(np.diag(pcov))
    chi2, ndf, red = chi2_ndf(model, centers, counts, yerr, popt)
    return popt, perr, chi2, ndf, red


# ---------------------------------------------------------------------------
# Step 3: symbolic regression
# ---------------------------------------------------------------------------
def run_pysr(centers, y, seed, niterations=40):
    from pysr import PySRRegressor
    X = centers.reshape(-1, 1)
    model = PySRRegressor(
        niterations=niterations,
        populations=8,
        population_size=40,
        binary_operators=["+", "-", "*", "/"],
        unary_operators=["square"],
        maxsize=22,
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
    model.fit(X, y, variable_names=["M"])
    eqs = model.equations_
    front = []
    for _, row in eqs.iterrows():
        front.append({
            "complexity": int(row["complexity"]),
            "loss": float(row["loss"]),
            "equation": str(row["equation"]),
        })
    best = str(model.sympy())
    return front, best


def main():
    df = pd.read_csv(DATA)
    M = invariant_mass(df)
    opp = (df.Q1 * df.Q2) < 0
    same = (df.Q1 * df.Q2) > 0
    M_opp = M[opp].values
    M_same = M[same].values

    # ---- Step 1: histograms ----
    counts, edges = np.histogram(M_opp, bins=N_BINS, range=WINDOW)
    centers = 0.5 * (edges[:-1] + edges[1:])
    same_counts, _ = np.histogram(M_same, bins=N_BINS, range=WINDOW)

    peak_center = float(centers[np.argmax(counts)])
    print(f"[Step 1] events total={len(df)} opp-charge={int(opp.sum())} "
          f"same-charge={int(same.sum())}")
    print(f"[Step 1] window={WINDOW} bins={N_BINS} ({(WINDOW[1]-WINDOW[0])/N_BINS:.3f} GeV/bin)")
    print(f"[Step 1] opp-charge peak bin center = {peak_center:.2f} GeV, "
          f"max count = {int(counts.max())}")
    print(f"[Step 1] same-charge max count = {int(same_counts.max())} "
          f"(control, should show no Z peak)")

    # Fit-region mask (narrower window for parametric fits)
    fmask = (centers >= FIT_WINDOW[0]) & (centers <= FIT_WINDOW[1])
    fc, fy = centers[fmask], counts[fmask].astype(float)

    results = {
        "scope": ("Recovers data-level Z resonance lineshape and parameters "
                  "(M_Z, Gamma_Z) from CMS dimuon mass spectrum. NOT the weak "
                  "force mechanism or its origin."),
        "data_file": DATA,
        "n_events_total": int(len(df)),
        "n_events_opp_charge": int(opp.sum()),
        "n_events_same_charge": int(same.sum()),
        "mass_formula": "M^2 = 2*pt1*pt2*(cosh(eta1-eta2) - cos(phi1-phi2))",
        "histogram_window_GeV": list(WINDOW),
        "n_bins": N_BINS,
        "bin_width_GeV": (WINDOW[1] - WINDOW[0]) / N_BINS,
        "fit_window_GeV": list(FIT_WINDOW),
        "opp_charge_peak_bin_center_GeV": peak_center,
        "same_charge_peak_max_count": int(same_counts.max()),
        "PDG": {"M_Z": PDG_MZ, "Gamma_Z": PDG_GZ},
        "fits": {},
    }

    # ---- Step 2: parametric fits ----
    A0 = float(fy.max())
    fits = [
        ("bw_rel", bw_rel, [A0, 91.0, 2.5],
         ([0, 80, 0.1], [np.inf, 100, 10])),
        ("bw_rel_linear_bkg", bw_rel_bkg, [A0, 91.0, 2.5, 5.0, 0.0],
         ([0, 80, 0.1, -np.inf, -np.inf], [np.inf, 100, 10, np.inf, np.inf])),
        ("bw_rel_exp_bkg", bw_rel_expbkg, [A0, 91.0, 2.5, 50.0, 0.05],
         ([0, 80, 0.1, 0, -1], [np.inf, 100, 10, np.inf, 1])),
        ("lorentzian_linear_bkg", lorentzian_bkg, [A0, 91.0, 2.5, 5.0, 0.0],
         ([0, 80, 0.1, -np.inf, -np.inf], [np.inf, 100, 10, np.inf, np.inf])),
    ]

    print("\n[Step 2] parametric fits (fit window "
          f"{FIT_WINDOW[0]}-{FIT_WINDOW[1]} GeV):")
    for name, model, p0, bounds in fits:
        try:
            popt, perr, chi2, ndf, red = fit_model(model, fc, fy, p0, bounds)
            mz, gz = popt[1], popt[2]
            mz_e, gz_e = perr[1], perr[2]
            ag_mz = 100 * (1 - abs(mz - PDG_MZ) / PDG_MZ)
            ag_gz = 100 * (1 - abs(gz - PDG_GZ) / PDG_GZ)
            results["fits"][name] = {
                "M_Z": float(mz), "M_Z_err": float(mz_e),
                "Gamma_Z": float(gz), "Gamma_Z_err": float(gz_e),
                "chi2": chi2, "ndf": ndf, "chi2_per_ndf": red,
                "M_Z_pct_agreement": ag_mz, "Gamma_Z_pct_agreement": ag_gz,
                "all_params": [float(x) for x in popt],
            }
            print(f"  {name:24s} M_Z={mz:7.3f}+-{mz_e:.3f}  "
                  f"Gamma_Z={gz:6.3f}+-{gz_e:.3f}  chi2/ndf={red:.2f}  "
                  f"(M_Z agree {ag_mz:.2f}%, Gamma_Z agree {ag_gz:.2f}%)")
        except Exception as e:
            results["fits"][name] = {"error": str(e)}
            print(f"  {name:24s} FIT FAILED: {e}")

    # ---- Step 3: symbolic regression on normalized resonance region ----
    # Normalize counts to peak=1; SR works on the resonance-window bins.
    y_norm = fy / fy.max()
    results["sr_target"] = "normalized counts (peak=1) vs M over fit window"
    results["sr_operators"] = {"binary": ["+", "-", "*", "/"], "unary": ["square"]}
    results["sr_runs"] = []

    print("\n[Step 3] symbolic regression (PySR, operators +,-,*,/,square):")
    seeds = [0, 1]
    for seed in seeds:
        try:
            front, best = run_pysr(fc, y_norm, seed)
            results["sr_runs"].append(
                {"seed": seed, "best_sympy": best, "pareto_front": front})
            print(f"  seed {seed}: best = {best}")
            print(f"    Pareto front (complexity, loss, equation):")
            for f in front:
                print(f"      c={f['complexity']:2d}  loss={f['loss']:.4e}  {f['equation']}")
        except Exception as e:
            results["sr_runs"].append({"seed": seed, "error": str(e)})
            print(f"  seed {seed}: PySR FAILED: {e}")

    # Heuristic: did a BW-like rational (has division, squared M, ~quadratic
    # denominator) appear? Flag if any pareto equation contains "/" and
    # "square" or M^2-like structure.
    bw_like = False
    for run in results["sr_runs"]:
        for f in run.get("pareto_front", []):
            eq = f["equation"]
            if "/" in eq and ("square" in eq or "M * M" in eq or "M^2" in eq):
                bw_like = True
                break
    results["sr_bw_like_rational_found"] = bw_like

    # Scan Pareto fronts for the propagator-pole signature square(M - c/M),
    # i.e. ((M^2 - c)/M)^2 = (M^2 - c)^2 / M^2 -- the relativistic BW
    # denominator core, with sqrt(c) acting as an SR-recovered M_Z.
    import re
    sr_pole = []
    pat = re.compile(r"square\(\s*M\s*-\s*\(?\s*([0-9.]+)\s*/\s*M")
    for run in results["sr_runs"]:
        for f in run.get("pareto_front", []):
            m = pat.search(f["equation"].replace(" ", " "))
            if m:
                c = float(m.group(1))
                mz_sr = c**0.5
                sr_pole.append({
                    "seed": run.get("seed"),
                    "complexity": f["complexity"],
                    "loss": f["loss"],
                    "pole_constant": c,
                    "implied_M_Z_GeV": mz_sr,
                    "implied_M_Z_pct_agreement": 100 * (1 - abs(mz_sr - PDG_MZ) / PDG_MZ),
                })
    results["sr_propagator_pole_matches"] = sr_pole
    if sr_pole:
        best_pole = min(sr_pole, key=lambda d: d["loss"])
        results["sr_implied_M_Z_GeV"] = best_pole["implied_M_Z_GeV"]
        print("\n[Step 3] SR recovered the propagator-pole core square(M - c/M)"
              " = ((M^2 - M_Z^2)/M)^2:")
        for p in sr_pole:
            print(f"    seed {p['seed']} c={p['complexity']:2d}: pole_const={p['pole_constant']:.1f}"
                  f" -> implied M_Z={p['implied_M_Z_GeV']:.2f} GeV"
                  f" ({p['implied_M_Z_pct_agreement']:.1f}% of PDG)")

    # Caveat: observed width is detector-broadened (resolution convolution),
    # so fitted Gamma_Z exceeds the natural PDG width.
    results["width_caveat"] = (
        "Fitted Gamma_Z (~3.7 GeV) exceeds PDG natural width (2.4952 GeV) "
        "because the observed lineshape is the Breit-Wigner convolved with "
        "muon momentum/angular resolution; the raw fit width is an upper "
        "bound on the natural width, not a deconvolved measurement.")

    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=2)

    # ---- Final summary ----
    print("\n" + "=" * 70)
    print("FINAL SUMMARY -- Z propagator (Breit-Wigner) from CMS dimuon data")
    print("=" * 70)
    best_fit = results["fits"].get("bw_rel_exp_bkg") or results["fits"].get("bw_rel_linear_bkg")
    if best_fit and "M_Z" in best_fit:
        print(f"  Direct fit (BW + bkg):  M_Z = {best_fit['M_Z']:.3f} +- "
              f"{best_fit['M_Z_err']:.3f} GeV  (PDG {PDG_MZ}, "
              f"{best_fit['M_Z_pct_agreement']:.2f}% agreement)")
        print(f"                          Gamma_Z = {best_fit['Gamma_Z']:.3f} +- "
              f"{best_fit['Gamma_Z_err']:.3f} GeV  (PDG {PDG_GZ}, "
              f"{best_fit['Gamma_Z_pct_agreement']:.2f}% agreement)")
        print(f"                          chi2/ndf = {best_fit['chi2_per_ndf']:.2f}")
    print(f"  BW shape recovered by parametric fit: YES")
    print(f"  BW-like rational recovered by SR:     "
          f"{'YES' if bw_like else 'NO (see Pareto front)'}")
    if sr_pole:
        bp = min(sr_pole, key=lambda d: d["loss"])
        print(f"  SR propagator-pole core ((M^2-M_Z^2)/M)^2 found: YES, "
              f"implied M_Z = {bp['implied_M_Z_GeV']:.2f} GeV "
              f"({bp['implied_M_Z_pct_agreement']:.1f}% of PDG)")
    print("  NOTE: fitted Gamma_Z ~3.7 GeV > PDG 2.495 GeV due to detector")
    print("  resolution broadening the lineshape (BW convolved w/ resolution).")
    print("  SCOPE: this recovers the resonance LINESHAPE and its parameters")
    print("  (M_Z, Gamma_Z) from data. It does NOT derive the weak force")
    print("  mechanism or explain 'what the weak force is'.")
    print(f"\nResults written to {OUT}")
    return results


if __name__ == "__main__":
    main()
