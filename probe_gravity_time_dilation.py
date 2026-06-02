"""
PROBE (backlog 6c, reframed O2 / H-C) -- gravity-couples-to-TIME, the clean test.

Force-coupling objects are construction-confounded (INFO-056): gravity is our only
time-series observable, so "only gravity flows in time" cannot be decided by force
data. The clean test of gravity-touches-time is a TIME-DILATION / CLOCK measurement.

This recovers the relativistic clock law directly from raw IGS/CODE precise
products (GPS satellite ECEF position + onboard clock, SP3, 5-min cadence). An
eccentric satellite's clock carries a periodic offset

    dt_rel(t) = -(2/c^2) * (r_vec . v_vec) = -(2/c^2) * r * rdot

(combined gravitational + velocity time dilation over the orbit; the identity
r_vec.v_vec = r*rdot holds in any frame because |r| is frame-invariant, so we get
it from SP3 position magnitudes alone -- no inertial-velocity gymnastics).

OD shape: give the method raw clock+orbit data, recover the COEFFICIENT k and the
FUNCTIONAL FORM (linear in r*rdot). Ground truth k = -2/c^2 = -2.2256e-17 s^2/m^2,
checkable. This is intrinsically a TIME law (clock rate vs gravitational/kinematic
state) -- the thing no gauge force has. We do NOT predetermine the outcome; if CODE
pre-corrects the term the recovery is null and that is also a data point.

Joint least squares per satellite:  clock(t) = c0 + c1 t + c2 t^2 + c3 t^3 + k*X,
X = r*rdot. The cubic absorbs secular clock drift; X cannot be absorbed (period
~12 h, ~2 cycles/day). Scramble null: shuffle X vs clock -> |k|, incremental R^2
should collapse.
"""
import numpy as np
import json, sys

C_LIGHT = 299792458.0
K_TRUTH = -2.0 / C_LIGHT**2          # -2.2256e-17 s^2/m^2
BAD_CLK = 999999.0                    # SP3 bad-clock flag (microseconds)


def parse_sp3(path, systems=("G",)):
    """Return dict prn -> (t[s], xyz[m] (N,3), clk[s]) and epoch times."""
    epoch_times = []
    recs = {}           # prn -> list of (epoch_idx, x, y, z, clk)
    ei = -1
    with open(path) as f:
        for line in f:
            if line.startswith("*"):
                p = line.split()
                # * yyyy mm dd hh mm ss.ssss
                yr, mo, dy = int(p[1]), int(p[2]), int(p[3])
                hh, mm = int(p[4]), int(p[5]); ss = float(p[6])
                # seconds from file start (day-internal; absolute offset irrelevant)
                tsec = ((dy * 24 + hh) * 60 + mm) * 60.0 + ss
                epoch_times.append(tsec)
                ei += 1
            elif line.startswith("P") and len(line) > 4 and line[1] in systems:
                prn = line[1:4]
                try:
                    x = float(line[4:18]); y = float(line[18:32]); z = float(line[32:46])
                    clk = float(line[46:60])
                except ValueError:
                    continue
                recs.setdefault(prn, []).append((ei, x, y, z, clk))
    et = np.array(epoch_times)
    et -= et[0]
    out = {}
    for prn, rows in recs.items():
        idx = np.array([r[0] for r in rows])
        xyz = np.array([[r[1], r[2], r[3]] for r in rows]) * 1000.0   # km -> m
        clk = np.array([r[4] for r in rows])                          # microseconds
        good = clk < BAD_CLK
        if good.sum() < 50:
            continue
        t = et[idx][good]
        xyz = xyz[good]
        clk = clk[good] * 1e-6                                        # microsec -> s
        out[prn] = (t, xyz, clk)
    return out, et


def fit_sat(t, xyz, clk, poly_deg=3):
    r = np.linalg.norm(xyz, axis=1)
    rdot = np.gradient(r, t)
    X = r * rdot                                  # m^2/s  (= r_vec . v_vec)
    ecc = (r.max() - r.min()) / (r.max() + r.min())
    # normalise t for conditioning
    tn = (t - t.mean()) / (t.max() - t.min() + 1e-12)
    cols = [tn**d for d in range(poly_deg + 1)]
    A_poly = np.vstack(cols).T
    A = np.column_stack([A_poly, X])
    coef, *_ = np.linalg.lstsq(A, clk, rcond=None)
    k = coef[-1]
    pred = A @ coef
    # incremental R^2 of the X term: compare full model vs poly-only model
    coef_p, *_ = np.linalg.lstsq(A_poly, clk, rcond=None)
    resid_p = clk - A_poly @ coef_p
    resid_f = clk - pred
    ss_tot = np.sum((clk - clk.mean())**2)
    r2_full = 1 - np.sum(resid_f**2) / ss_tot
    # fraction of the poly-residual variance explained by X
    ss_p = np.sum(resid_p**2)
    incr = 1 - np.sum(resid_f**2) / ss_p if ss_p > 0 else 0.0
    amp = 2.0 * np.sqrt(3.986004418e14 * r.mean()) * ecc / C_LIGHT**2  # expected ns-scale
    pred_true = K_TRUTH * X                      # known-law term from raw orbit geometry
    pred_amp_true_ns = float(np.std(pred_true) * 1e9)
    resid_rms_ns = float(np.std(resid_p) * 1e9)
    removal_ratio = resid_rms_ns / (pred_amp_true_ns + 1e-12)  # ~1 if present, ~0 if removed
    return dict(prn_ecc=float(ecc), k=float(k), k_over_truth=float(k / K_TRUTH),
                r2_full=float(r2_full), incr_r2_X=float(incr),
                pred_amp_ns=float(amp * 1e9), pred_amp_true_ns=pred_amp_true_ns,
                resid_rms_ns=resid_rms_ns, removal_ratio=float(removal_ratio))


def scramble_null(t, xyz, clk, poly_deg=3, n=200, seed=0):
    r = np.linalg.norm(xyz, axis=1)
    rdot = np.gradient(r, t)
    X = r * rdot
    tn = (t - t.mean()) / (t.max() - t.min() + 1e-12)
    A_poly = np.vstack([tn**d for d in range(poly_deg + 1)]).T
    coef_p, *_ = np.linalg.lstsq(A_poly, clk, rcond=None)
    resid_p = clk - A_poly @ coef_p
    ss_p = np.sum(resid_p**2)
    rng = np.random.default_rng(seed)
    ks, incrs = [], []
    for _ in range(n):
        Xs = rng.permutation(X)
        A = np.column_stack([A_poly, Xs])
        coef, *_ = np.linalg.lstsq(A, clk, rcond=None)
        ks.append(coef[-1])
        resid_f = clk - A @ coef
        incrs.append(1 - np.sum(resid_f**2) / ss_p if ss_p > 0 else 0.0)
    return np.array(ks), np.array(incrs)


def main(canary=False):
    sp3 = "data/gps/COD_20230010000_05M_ORB.SP3"
    sats, et = parse_sp3(sp3, systems=("G", "E"))   # GPS + Galileo (incl. GREAT e~0.16)
    prns = sorted(sats)
    if canary:
        # pick the 5 highest-eccentricity GPS sats (largest signal)
        eccs = {p: (np.ptp(np.linalg.norm(sats[p][1], axis=1))
                    / (2 * np.linalg.norm(sats[p][1], axis=1).mean())) for p in prns}
        prns = sorted(prns, key=lambda p: -eccs[p])[:5]
    results = {}
    ks = []
    for p in prns:
        t, xyz, clk = sats[p]
        res = fit_sat(t, xyz, clk)
        results[p] = res
        ks.append(res["k"])
    ks = np.array(ks)
    # scramble null on the highest-ecc sat
    top = max(results, key=lambda p: results[p]["prn_ecc"])
    nk, nincr = scramble_null(*sats[top], n=(50 if canary else 300))
    # eccentricity-scaling check: predicted-term amplitude vs eccentricity (raw orbit)
    ecc_arr = np.array([results[p]["prn_ecc"] for p in prns])
    amp_arr = np.array([results[p]["pred_amp_true_ns"] for p in prns])
    slope_amp_ecc = float(np.polyfit(ecc_arr, amp_arr, 1)[0])
    rmean = float(np.mean([np.linalg.norm(sats[p][1], axis=1).mean() for p in prns]))
    slope_truth = 2.0 * np.sqrt(3.986004418e14 * rmean) / C_LIGHT**2 * 1e9  # ns per unit ecc
    removal = np.array([results[p]["removal_ratio"] for p in prns])
    summary = dict(
        k_truth=K_TRUTH, n_sats=len(prns),
        k_mean=float(ks.mean()), k_std=float(ks.std()),
        k_mean_over_truth=float(ks.mean() / K_TRUTH),
        median_incr_r2=float(np.median([results[p]["incr_r2_X"] for p in results])),
        median_r2_full=float(np.median([results[p]["r2_full"] for p in results])),
        median_removal_ratio=float(np.median(removal)),
        max_pred_amp_ns=float(amp_arr.max()), median_resid_rms_ns=float(
            np.median([results[p]["resid_rms_ns"] for p in results])),
        amp_vs_ecc_slope_recovered=slope_amp_ecc, amp_vs_ecc_slope_truth=slope_truth,
        amp_vs_ecc_slope_ratio=float(slope_amp_ecc / slope_truth),
        null_top_sat=top,
        null_incr_r2_mean=float(nincr.mean()), null_incr_r2_std=float(nincr.std()),
        real_incr_r2_top=float(results[top]["incr_r2_X"]),
    )
    out = dict(summary=summary, per_sat=results)
    tag = "_canary" if canary else "_results"
    fn = f"probe_gravity_time_dilation{tag}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {fn}  ({len(prns)} sats)")
    print("\nper-sat (prn ecc k/truth incrR2):")
    for p in sorted(results, key=lambda p: -results[p]["prn_ecc"])[:12]:
        r = results[p]
        print(f"  {p}  ecc={r['prn_ecc']:.4f}  k/truth={r['k_over_truth']:+.3f}  "
              f"incrR2={r['incr_r2_X']:.3f}  amp={r['pred_amp_ns']:.1f}ns")


if __name__ == "__main__":
    main(canary=("--canary" in sys.argv))
