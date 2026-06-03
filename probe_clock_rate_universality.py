#!/usr/bin/env python3
"""
PROBE: Clock-rate law UNIVERSALITY  (Prong B #1 of the WHAT IS TIME backlog)
============================================================================

Question (backlog, S19):
  All of our POSITIVE "time" findings are CLOCK-RATE changes with gravitational
  potential / velocity:
    - LIGO  : gravitational-wave inspiral governed by GR proper time (chirp mass)
    - GPS   : the relativistic eccentricity term, coefficient -2/c^2
    - PULSAR: PSR B1913+16 Einstein-delay gamma (grav. redshift + 2nd-order Doppler)
              and orbital decay dP_b/dt
  Each was extracted INDEPENDENTLY, from a different raw dataset, by a different
  method, in a wildly different physical regime.  Are they ONE law -- the
  equivalence-principle / proper-time relation
        dtau/dt = 1 + Phi/c^2 - v^2/(2 c^2)         (weak field)
  ?

WHAT THIS SCRIPT DOES (the NEW content -- the cross-system synthesis):
  It does NOT re-extract the coefficients (those independent raw-data extractions
  are the prior probes: probe_6c_gps_positive / probe_time_dipole_gps_strengthen,
  probe_pulsar_time, probe_o3_gw170817_chirpmass / probe_gravity_chirp_ridge).
  It LOADS each independently-recovered coefficient, normalises it to the GR
  prediction for that system's PROJECTION of the single proper-time relation,
  forms the dimensionless ratio R_i = recovered / GR, places each on a common
  regime axis (the characteristic fractional clock-rate effect ~ (v/c)^2 / Phi/c^2
  the system probes), and tests whether {R_i} is ONE common value (= 1 if it is
  the GR proper-time law).

DEFLATIONARY BRAKE (kept on, per Operating Rules):
  The proper-time relation IS the equivalence principle -- well-established
  physics.  A positive result here is a POSITIVE CONTROL that the OD tool,
  applied independently to three raw datasets in three regimes, lands on ONE
  proper-time law.  It is NOT a new-physics claim.  The non-trivial content is
  (a) independence of the three extractions, (b) the ~9-decade span in field
  strength they cover, (c) that the SAME functional family fits all three.

Usage:
  python3 probe_clock_rate_universality.py            # full run, writes results
  python3 probe_clock_rate_universality.py --canary   # validate inputs only
"""
import json, sys, os, math

HERE = os.path.dirname(os.path.abspath(__file__))
CANARY = "--canary" in sys.argv

C = 299792458.0          # m/s
G = 6.67430e-11          # SI
M_SUN = 1.98892e30       # kg
M_EARTH = 5.9722e24      # kg
GM_EARTH = 3.986004418e14  # m^3/s^2 (WGS-84, more precise than G*M_earth)


def load(fname):
    p = os.path.join(HERE, fname)
    if not os.path.exists(p):
        raise FileNotFoundError(p)
    with open(p) as f:
        return json.load(f)


# ----------------------------------------------------------------------------
# Pull each system's INDEPENDENTLY-recovered coefficient + its GR prediction.
# Each entry records: the proper-time relation it projects, the recovered/GR
# ratio R (with an honest uncertainty), and the regime it probes ((v/c)^2-ish).
# ----------------------------------------------------------------------------
def build_systems():
    systems = []

    # ---- GPS: relativistic eccentricity term  dt_rel = -2(r.v)/c^2 = F e sqrt(a) sinE
    # The proper-time relation's PERIODIC part along an eccentric orbit. Coefficient -2/c^2.
    gps_pos = load("probe_6c_gps_positive_results.json")["summary"]
    gps_str = load("probe_time_dipole_gps_strengthen_results.json")["aggregate"]
    # cross-check the GR coefficient F = -2 sqrt(GM)/c^2 against the stored truth
    F_gr = -2.0 * math.sqrt(GM_EARTH) / C**2          # s / sqrt(m)
    F_stored = gps_pos["F_truth"]
    F_check = abs(F_gr - F_stored) / abs(F_stored)
    # robust aggregate ratio (multi-station/day eccentric Galileo, real scatter)
    ecc = gps_str["C1"]["eccentric_signal"]["k_over_truth"]
    R_gps, dR_gps, n_gps = ecc["mean"], ecc["std"], ecc["n"]
    # cleanest single-day eccentric pair (E18 z=380, E14) for the record
    gps_clean = list(gps_pos["gal_k_over_truth"])
    # regime: GNSS orbital speed ~3.874 km/s -> (v/c)^2
    v_gnss = 3874.0
    systems.append(dict(
        name="GPS (eccentric Galileo)",
        dataset="raw RINEX pseudorange, IGS stations BRUX/ALGO (2-3 days)",
        observable="relativistic eccentricity term  dt = -2(r.v)/c^2",
        projection="periodic part of dtau/dt along an eccentric orbit; coefficient -2/c^2",
        R=R_gps, dR=dR_gps, n_independent=n_gps,
        R_clean_singleday=gps_clean,
        regime_vc2=(v_gnss / C)**2,
        gr_coeff_check_rel_err=F_check,
        significance="z=380 vs scramble null (E18, detection); coeff scatter +/-%.2f over %d replicas" % (dR_gps, n_gps),
        use_in_fit=True,
    ))

    # ---- PULSAR: PSR B1913+16 Einstein-delay gamma (grav redshift + 2nd-order Doppler)
    pul = load("probe_pulsar_time_results.json")["experiments"]
    g = pul["exp3_einstein_delay_gamma"]
    p = pul["exp2_direct_pbdot"]
    R_pul_gamma = g["ratio_recovered_to_published"]
    R_pul_pbdot = p["ratio_recovered_to_GR"]
    # honest uncertainty: the probe flags formal errors as optimistic; use the
    # recovered formal fractional error but DO NOT treat its sigma as gospel.
    dR_pul_gamma = g["recovered_GAMMA_err_s"] / g["recovered_GAMMA_s"]
    # regime: companion potential at the pulsar ~ Gm_c/(c^2 a); pulsar v/c ~ 1e-3.
    # use the orbital 2nd-order-Doppler scale (v/c)^2 with v ~ 2 pi a / Pb.
    # Pb in days; a sin i ~ 2.34 ls; take v ~ 300 km/s typical for this binary.
    v_psr = 3.0e5
    systems.append(dict(
        name="Pulsar PSR B1913+16 (Einstein delay)",
        dataset="raw Arecibo TOAs 1981-2012 (Weisberg & Huang 2016)",
        observable="Einstein-delay gamma (grav. redshift + 2nd-order Doppler of pulsar clock)",
        projection="static/orbit-averaged part of dtau/dt in companion potential -> gamma",
        R=R_pul_gamma, dR=max(dR_pul_gamma, 1e-4), n_independent=1,
        R_secondary_pbdot=R_pul_pbdot,
        regime_vc2=(v_psr / C)**2,
        gr_coeff_check_rel_err=None,
        significance="gamma recovered to 0.014%% of published GR value; dP_b/dt ratio %.3f" % R_pul_pbdot,
        use_in_fit=True,
    ))

    # ---- LIGO: GW170817 BNS chirp mass (PN proper-time inspiral) -- matched filter
    lg = load("probe_o3_gw170817_chirpmass_results.json")
    rec = lg["recovery"]
    R_ligo = rec["vs_catalog_detframe_ratio"]
    # uncertainty: chirp-mass bank spacing 0.005 Msun around 1.2 -> ~0.4% ; plus
    # H1 vs L1 spread (1.200 vs 1.195) ~ 0.4%. Use the H1/L1 spread as honest dR.
    Mc_H1, Mc_L1 = rec["per_det_best_Mc"]["H1"], rec["per_det_best_Mc"]["L1"]
    dR_ligo = abs(Mc_H1 - Mc_L1) / lg["catalog"]["Mc_detframe_Msun"]
    # regime: BNS at merger v/c ~ (pi M f)^(1/3); at f~1 kHz, M~2.7 Msun -> v/c ~0.4
    systems.append(dict(
        name="LIGO GW170817 (BNS chirp)",
        dataset="GWOSC GW170817 H1+L1 strain, 4096 Hz",
        observable="chirp mass M_c via TaylorF2 matched filter",
        projection="integrated PN proper-time phase evolution of the inspiral",
        R=R_ligo, dR=max(dR_ligo, 0.004), n_independent=2,  # H1 + L1 agree independently
        regime_vc2=0.4**2,
        gr_coeff_check_rel_err=None,
        significance="net SNR 14.8 (L1 null 9.7 sigma); H1=L1 independently; 0.19%% from catalog",
        use_in_fit=True,
    ))

    # ---- LIGO GW150914 ridge (method-limited; reported, NOT used in strict fit) ----
    try:
        rg = load("probe_gravity_chirp_ridge_results.json")["events"]["GW150914"]["H1"]["fit"]
        Mc_ridge = rg["chirp_mass_Msun"]
        # catalog detector-frame chirp mass GW150914 ~ 30.8 Msun
        R_ridge = Mc_ridge / 30.8
        systems.append(dict(
            name="LIGO GW150914 (ridge, method-limited)",
            dataset="GWOSC GW150914 H1 strain (Morlet-CWT ridge)",
            observable="chirp mass via f^(-8/3)-linear-in-t ridge",
            projection="Newtonian late-inspiral proper-time chirp (known ~24% high bias)",
            R=R_ridge, dR=0.20, n_independent=1,
            regime_vc2=0.45**2,
            gr_coeff_check_rel_err=None,
            significance="ridge fit R^2=%.3f; absolute M_c biased high (Newtonian late-inspiral)" % rg["r2"],
            use_in_fit=False,
        ))
    except Exception as e:
        pass

    return systems


def universality_test(systems):
    fit = [s for s in systems if s["use_in_fit"]]
    Rs = [s["R"] for s in fit]
    dRs = [s["dR"] for s in fit]

    # (a) point-estimate agreement: each R vs 1.0
    dev_pct = [100.0 * (R - 1.0) for R in Rs]
    max_abs_dev = max(abs(d) for d in dev_pct)

    # (b) robust spread across systems (unweighted -- treats each system equally,
    #     immune to the heterogeneous/optimistic formal errors)
    mean_R = sum(Rs) / len(Rs)
    var_R = sum((R - mean_R)**2 for R in Rs) / (len(Rs) - 1)
    std_R = math.sqrt(var_R)

    # (c) single-common-law fit: weighted mean alpha (recovered = alpha * GR),
    #     and chi^2/dof of "all share one alpha".  Flag formal-error caveat.
    w = [1.0 / d**2 for d in dRs]
    alpha = sum(wi * R for wi, R in zip(w, Rs)) / sum(w)
    chi2_common = sum(wi * (R - alpha)**2 for wi, R in zip(w, Rs))
    dof_common = len(Rs) - 1
    chi2_vs_one = sum(wi * (R - 1.0)**2 for wi, R in zip(w, Rs))
    dof_one = len(Rs)

    # regime span
    regimes = [s["regime_vc2"] for s in fit]
    span_decades = math.log10(max(regimes) / min(regimes))

    return dict(
        systems_in_fit=[s["name"] for s in fit],
        ratios_R=Rs,
        deviation_from_unity_pct=dev_pct,
        max_abs_deviation_pct=max_abs_dev,
        robust_mean_R=mean_R,
        robust_std_R=std_R,
        weighted_alpha=alpha,
        chi2_common_law=chi2_common,
        dof_common_law=dof_common,
        chi2_vs_GR_unity=chi2_vs_one,
        dof_vs_GR_unity=dof_one,
        regime_span_decades_in_vc2=span_decades,
        formal_error_caveat=(
            "weighted chi^2 uses each probe's formal error; GPS uses real "
            "multi-replica scatter (honest), but pulsar gamma and LIGO bank "
            "errors are optimistic (per their own probes). The ROBUST headline "
            "is the point-estimate agreement (max_abs_deviation_pct) + the "
            "system-to-system spread (robust_std_R), NOT the formal chi^2."
        ),
    )


def main():
    systems = build_systems()

    if CANARY:
        out = {"canary": True,
               "n_systems": len(systems),
               "systems": [{"name": s["name"], "R": s["R"], "dR": s["dR"],
                            "regime_vc2": s["regime_vc2"], "use_in_fit": s["use_in_fit"]}
                           for s in systems]}
        print(json.dumps(out, indent=2))
        with open(os.path.join(HERE, "probe_clock_rate_universality_canary.json"), "w") as f:
            json.dump(out, f, indent=2)
        print("\n[canary] inputs load and ratios build OK.")
        return

    res = universality_test(systems)

    results = dict(
        probe="clock_rate_law_universality",
        question=("Are the independently-recovered LIGO / GPS / pulsar clock-rate "
                  "signatures ONE law -- the proper-time relation "
                  "dtau/dt = 1 + Phi/c^2 - v^2/(2c^2)?"),
        common_law="dtau/dt = 1 + Phi/c^2 - v^2/(2 c^2)  (weak-field proper time)",
        systems=systems,
        universality=res,
        reading=dict(
            DATA=("Three independent raw-data extractions, three methods, three "
                  "regimes spanning ~%.1f decades in (v/c)^2, each recover the "
                  "GR proper-time coefficient: ratios R = %s (each within %.1f%% "
                  "of unity); system-to-system spread %.3f."
                  % (res["regime_span_decades_in_vc2"],
                     ["%.4f" % r for r in res["ratios_R"]],
                     res["max_abs_deviation_pct"], res["robust_std_R"])),
            INTERPRETATION=("The three observables are different projections of ONE "
                  "relation: GPS = its periodic part (-2(r.v)/c^2), pulsar gamma = "
                  "its orbit-averaged grav-redshift+2nd-order-Doppler part, LIGO = "
                  "its integrated PN proper-time phase. They reduce to one law."),
            FRAME=("POSITIVE CONTROL, not new physics: this IS the equivalence "
                  "principle / GR proper time -- established. The non-trivial part "
                  "is that the OD tool, run independently on three raw datasets in "
                  "three regimes, lands on the same proper-time coefficient. "
                  "GW150914 ridge is reported but excluded from the strict fit "
                  "(known Newtonian-late-inspiral absolute-mass bias)."),
        ),
    )

    with open(os.path.join(HERE, "probe_clock_rate_universality_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    # human-readable summary
    print("=" * 74)
    print("CLOCK-RATE LAW UNIVERSALITY  --  is it ONE proper-time law?")
    print("=" * 74)
    print("common law:  dtau/dt = 1 + Phi/c^2 - v^2/(2 c^2)\n")
    print("%-38s %8s %8s %12s" % ("system", "R=rec/GR", "+/-", "regime (v/c)^2"))
    print("-" * 74)
    for s in systems:
        tag = "" if s["use_in_fit"] else "  [excl: method-limited]"
        print("%-38s %8.4f %8.3f %12.2e%s"
              % (s["name"][:38], s["R"], s["dR"], s["regime_vc2"], tag))
    print("-" * 74)
    u = res
    print("regime span in (v/c)^2          : %.1f orders of magnitude" % u["regime_span_decades_in_vc2"])
    print("each R within                   : %.2f%% of unity (GR)" % u["max_abs_deviation_pct"])
    print("system-to-system spread of R    : %.3f (robust, unweighted)" % u["robust_std_R"])
    print("weighted common-law scale alpha : %.4f" % u["weighted_alpha"])
    print("chi2/dof (all share one alpha)  : %.2f / %d" % (u["chi2_common_law"], u["dof_common_law"]))
    print("chi2/dof (all = GR, R=1)        : %.2f / %d  [formal errs optimistic -- see caveat]"
          % (u["chi2_vs_GR_unity"], u["dof_vs_GR_unity"]))
    print()
    print("READING (deflationary brake ON):")
    print("  DATA :", results["reading"]["DATA"])
    print("  FRAME:", results["reading"]["FRAME"])
    print("\nwrote probe_clock_rate_universality_results.json")


if __name__ == "__main__":
    main()
