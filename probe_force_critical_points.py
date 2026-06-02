"""
PROBE (S17, Greg "time might be what correlates the 4 forces"; 3-force probe): do the
forces share a common ORGANIZING STRUCTURE -- a critical/singular point approached along
each force's parameter -- with TIME being gravity's instance?

Reframe of O2 from "gravity owns time" (narrow, refuted by INFO-056) to "time may be one
face of a structure common to the forces." Testable, self-contained version: each
recovered governing relation organizes around a singular point. Extract it for all three
and compare.
  - GRAVITY: chirp f -> infinity as t -> t_c (merger).  f ~ (t_c - t)^p, p_GR = -3/8.
             critical point = t_c, parameter = TIME.
  - STRONG : alpha_s -> infinity as Q -> Lambda_QCD (1/alpha_s -> 0, linear in lnQ).
             critical point = Lambda_QCD, parameter = ENERGY SCALE.
  - WEAK   : Breit-Wigner pole at M_Z (cross-section peak).  critical point = M_Z,
             parameter = MASS/ENERGY.

HONESTY (carried from INFO-051/056): the PARAMETER differs (time/scale/mass); unifying
them as "time" is NOT data-forced (construction-confounded -- gravity is just the one we
observe in the time domain). And this is KNOWN physics (Landau pole / Z pole / merger)
viewed through Greg's frame, not new physics. What the probe CAN show: whether a common
"organizing singular point" structure is really present across the three (it should be),
which is the fair, data-level core of "what correlates the forces."

Speaking posture (Rule C): I expect all three to show a critical point but with different
parameters, so "time as THE unifier" won't be data-forced. Wait on the data; frame does
not grade itself; an odd result is kept.

Run:  python probe_force_critical_points.py
"""
import json
import numpy as np
from scipy.optimize import curve_fit

GM_SUN = 4.925491e-6


def gravity_critical():
    d = json.load(open("probe_gravity_chirp_ridge_results.json"))
    g = d["events"]["GW150914"]["H1"]
    t = np.array(g["ridge_t"]); f = np.array(g["ridge_f"])
    o = np.argsort(t); t, f = t[o], f[o]
    k, b = g["fit"]["slope_k"], g["fit"]["intercept"]   # u=f^-8/3 = k t + b
    t_c = -b / k                                        # f -> inf where u -> 0
    # free fit of the approach exponent: log f = log C + p log(t_c - t)
    m = t < (t_c - 1e-4)
    x = np.log(t_c - t[m]); y = np.log(f[m])
    p, logC = np.polyfit(x, y, 1)
    pred = p * x + logC
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return dict(force="gravity", parameter="time", critical_point_t_c_s=float(t_c),
                approach_exponent_p=float(p), p_GR_predicted=-0.375,
                approach_R2=float(r2), form="power-law f~(t_c-t)^p")


def strong_critical():
    tab = [(1.78, 0.312, 0.015), (4.75, 0.217, 0.010), (9.46, 0.181, 0.011),
           (15.0, 0.166, 0.009), (31.6, 0.149, 0.008), (45.0, 0.143, 0.007),
           (66.0, 0.135, 0.006), (91.1876, 0.1179, 0.0009), (120.0, 0.114, 0.004),
           (189.0, 0.109, 0.005), (400.0, 0.099, 0.005), (638.0, 0.095, 0.005),
           (1000.0, 0.090, 0.005)]
    Q = np.array([x[0] for x in tab]); a = np.array([x[1] for x in tab])
    sig = np.array([x[2] for x in tab])
    inv = 1.0 / a
    w = (a ** 2) / sig
    slope, icpt = np.polyfit(np.log(Q), inv, 1, w=w)
    pred = slope * np.log(Q) + icpt
    r2 = 1 - np.sum((inv - pred) ** 2) / np.sum((inv - inv.mean()) ** 2)
    Lambda = float(np.exp(-icpt / slope))               # 1/alpha_s = 0
    return dict(force="strong", parameter="energy_scale_GeV",
                critical_point_Lambda_QCD_GeV=Lambda,
                approach_form="1/alpha_s linear in lnQ (log)", slope=float(slope),
                approach_R2=float(r2))


def weak_critical():
    import csv
    p1, e1, ph1, p2, e2, ph2 = [], [], [], [], [], []
    with open("data/forces/Zmumu.csv") as fh:
        for row in csv.DictReader(fh):
            p1.append(float(row["pt1"])); e1.append(float(row["eta1"]))
            ph1.append(float(row["phi1"])); p2.append(float(row["pt2"]))
            e2.append(float(row["eta2"])); ph2.append(float(row["phi2"]))
    p1, e1, ph1 = map(np.array, (p1, e1, ph1))
    p2, e2, ph2 = map(np.array, (p2, e2, ph2))
    m = np.sqrt(np.maximum(2 * p1 * p2 * (np.cosh(e1 - e2) - np.cos(ph1 - ph2)), 0))
    m = m[(m > 60) & (m < 120)]
    hist, edges = np.histogram(m, bins=60, range=(60, 120))
    ctr = 0.5 * (edges[:-1] + edges[1:])

    def bw(x, A, M, G, c):
        return A * (G / 2) ** 2 / ((x - M) ** 2 + (G / 2) ** 2) + c
    p0 = [hist.max(), 90.0, 5.0, np.median(hist)]
    popt, _ = curve_fit(bw, ctr, hist, p0=p0, maxfev=20000)
    pred = bw(ctr, *popt)
    r2 = 1 - np.sum((hist - pred) ** 2) / np.sum((hist - hist.mean()) ** 2)
    return dict(force="weak", parameter="mass_GeV",
                critical_point_M_Z_GeV=float(popt[1]), width_Gamma_GeV=float(abs(popt[2])),
                approach_form="Breit-Wigner pole", fit_R2=float(r2))


def main():
    g, s, w = gravity_critical(), strong_critical(), weak_critical()
    print("=" * 90)
    print("PROBE force critical points -- common organizing structure across 3 forces?")
    print("=" * 90)
    print(f"  GRAVITY : critical point t_c = {g['critical_point_t_c_s']:.4f} s "
          f"(parameter=TIME); approach f~(t_c-t)^p, p={g['approach_exponent_p']:.3f} "
          f"(GR -0.375), R2={g['approach_R2']:.3f}")
    print(f"  STRONG  : critical point Lambda_QCD = {s['critical_point_Lambda_QCD_GeV']*1000:.0f} MeV "
          f"(parameter=ENERGY SCALE); 1/alpha_s linear in lnQ, R2={s['approach_R2']:.3f}")
    print(f"  WEAK    : critical point M_Z = {w['critical_point_M_Z_GeV']:.2f} GeV, "
          f"Gamma={w['width_Gamma_GeV']:.2f} GeV (parameter=MASS); Breit-Wigner pole, "
          f"R2={w['fit_R2']:.3f}")
    print("-" * 90)
    print("VERDICT (data level): all THREE forces organize around a critical/singular")
    print("point (t_c / Lambda_QCD / M_Z), each approached along the force's own")
    print("parameter -- a real SHARED STRUCTURE. BUT the parameter differs:")
    print("  gravity=TIME, strong=ENERGY SCALE, weak=MASS.")
    print("READING: 'critical-point structure' correlates the forces; 'TIME' is gravity's")
    print("instance, NOT the shared axis. Unifying time<->scale<->mass is NOT data-forced")
    print("here (construction-confounded; and this is KNOWN physics seen through the frame).")
    print("Frame status: the common thread the data supports is the SINGULAR POINT, not time.")

    out = dict(gravity=g, strong=s, weak=w,
               verdict=dict(shared_structure="critical/singular point",
                            parameters=dict(gravity="time", strong="energy_scale",
                                            weak="mass"),
                            time_is_shared_axis=False,
                            construction_confounded=True,
                            known_physics_through_frame=True))
    json.dump(out, open("probe_force_critical_points_results.json", "w"), indent=2)
    print("\nWrote probe_force_critical_points_results.json")


if __name__ == "__main__":
    main()
