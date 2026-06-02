"""
PROBE 6c (POSITIVE half) -- recover the GENERAL-RELATIVISTIC clock /
time-dilation governing law from RAW public GPS data, with an actual
recovered COEFFICIENT. Gravity analogue of recovering the Z Breit-Wigner
or the QCD running from data.

TARGET LAW (eccentric GNSS satellite onboard-clock periodic offset):
    dt_rel(t) = -(2/c^2) (r_vec . v_vec) = F * e * sqrt(a) * sin(E)
    F = -2 sqrt(mu)/c^2 = -4.442807633e-10 s/sqrt(m),  mu=3.986004418e14
    ground-truth coefficient -2/c^2 = -2.2256e-17 s^2/m^2.

WHY THE PRECISE-PRODUCT ROUTE WAS A NULL (already established):
The CODE/IGS PRECISE SP3 satellite clocks have the periodic relativistic
term MODELLED OUT in standard processing, so it is absent from the cleaned
product (recovered k ~0.4% of truth in probe_gravity_time_dilation.py).

WHY THE BROADCAST-DIFFERENCING ROUTE (Route 1) IS ALSO A NULL:
The broadcast nav clock polynomial (af0+af1 dt+af2 dt^2) EXCLUDES the
periodic term by convention -- and the PRECISE SP3 clock also has it removed.
So precise_clock - broadcast_poly contains NO relativistic term: the two
agree to ~3 ns over a day on E14/E18 while the term is ~760 ns ptp. This
script records that as a checked data point (route1_null in results), then
takes the working route:

ROUTE 2 (TERM-RETAINING -- raw observations). The pseudorange observation
PHYSICALLY contains the satellite's true apparent clock, which carries the
relativistic offset. Using the PRECISE SP3 satellite clock (term removed) to
correct a raw RINEX pseudorange leaves the relativistic term in the residual:

    resid_metres = P - rho_geom(SP3 orbit, station, Sagnac, travel-time)
                       + c * dt_sat_precise
                 = c*dt_recv(t) + tropo + iono + (-c*dt_rel) + noise

The receiver clock dt_recv is common to all sats at an epoch. We estimate it
per epoch from NEAR-CIRCULAR sats (whose dt_rel is negligible: e<0.03 ->
< ~15 ns), subtract it from the ECCENTRIC target sat, and the surviving
periodic content over the day IS the relativistic clock modulation. We fit

    cleaned_resid = poly2(t) + k * [ -c * e*sqrt(a)*sin(E) ]

so the recovered k matches F directly (the leading minus folds the
pseudorange sign convention P = rho + c(dt_recv - dt_sat) into the regressor).
The regressor e*sqrt(a)*sin(E) is computed PURELY from broadcast Keplerian
elements -- independent of the pseudorange data -- so a high-R^2 recovery is a
genuine cross-source agreement, not a tautology.

Targets: Galileo E14, E18 (e~0.162, ~275 ns / ~82 m signal) + high-ecc GPS.
Station: BRUX (Brussels), known IGS ECEF. One day (2023-001).
Null: scramble / time-slide the regressor vs the cleaned residual.

Pure numpy/scipy + hatanaka (RINEX Hatanaka decompression only).
"""
import numpy as np, json, sys, math

C_LIGHT = 299792458.0
MU = 3.986004418e14
F_TRUTH = -2.0 * math.sqrt(MU) / C_LIGHT**2        # -4.442807633e-10 s/sqrt(m)
K2_TRUTH = -2.0 / C_LIGHT**2                         # -2.2256e-17 s^2/m^2

from _route2_obs import (parse_sp3_full, parse_obs_c1, geom_residual,
                         dt_rel_reg, sp3_interp)

STATION = np.array([4027881.8468, 306998.2610, 4919498.6524])  # BRUX ECEF (m)
SP3 = "data/gps/COD_20230010000_05M_ORB.SP3"
NAV = "data/gps/BRDC00WRD_R_20230010000_01D_MN.rnx"
OBS = "data/gps/BRUX00BEL_R_20230010000_01D_30S_MO.rnx"


# --------------------------------------------------------------------------
# RINEX 3 nav parser (re-used by Route-1 null check + regressor elements)
# --------------------------------------------------------------------------
def _f(s):
    return float(s.replace("D", "E").replace("d", "e"))


def parse_rinex_nav(path, systems=("G", "E")):
    out = {}
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines) and "END OF HEADER" not in lines[i]:
        i += 1
    i += 1
    n = len(lines)
    while i < n:
        line = lines[i]
        if len(line) < 3 or line[0] not in systems or not line[1:3].strip().isdigit():
            i += 1; continue
        prn = line[0:3]
        try:
            dy = int(line[12:14]); hh = int(line[15:17]); mm = int(line[18:20]); ss = int(line[21:23])
            af0 = _f(line[23:42]); af1 = _f(line[42:61]); af2 = _f(line[61:80])
        except (ValueError, IndexError):
            i += 1; continue
        toc_sec = ((dy - 1) * 24 + hh) * 3600.0 + mm * 60.0 + ss
        if i + 7 >= n:
            break
        rows = lines[i + 1:i + 8]
        def vals(L):
            return [_f(L[4 + 19 * k:4 + 19 * (k + 1)]) for k in range(4)]
        try:
            r1 = vals(rows[0]); r2 = vals(rows[1]); r3 = vals(rows[2])
        except (ValueError, IndexError):
            i += 8; continue
        _, crs, dn, M0 = r1
        cuc, e, cus, sqrt_a = r2
        toe = r3[0]
        out.setdefault(prn, []).append(dict(
            toc_sec=toc_sec, af0=af0, af1=af1, af2=af2,
            toe_sec=toe % 86400.0, e=e, sqrt_a=sqrt_a, M0=M0, dn=dn))
        i += 8
    return out


def kepler_E(M, e, tol=1e-12, itmax=60):
    E = np.array(M, dtype=float)
    for _ in range(itmax):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E = E - dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


# --------------------------------------------------------------------------
# Route 1 null check (precise SP3 clock minus broadcast polynomial)
# --------------------------------------------------------------------------
def route1_null_check(nav):
    sats, _ = parse_sp3_full(SP3, systems=("E",))
    out = {}
    for prn in ("E14", "E18"):
        if prn not in sats or prn not in nav:
            continue
        d = sats[prn]; t = d["t"]; clk = d["clk"]; eph = nav[prn]
        tocs = np.array([e["toc_sec"] for e in eph])
        bc = np.empty_like(t)
        R = np.empty_like(t)
        for j, tj in enumerate(t):
            k = int(np.argmin(np.abs(tocs - tj)))
            e_ = eph[k]
            dtc = tj - e_["toc_sec"]
            if dtc > 43200: dtc -= 86400
            if dtc < -43200: dtc += 86400
            bc[j] = e_["af0"] + e_["af1"] * dtc + e_["af2"] * dtc**2
            n = math.sqrt(MU / e_["sqrt_a"]**6) + e_["dn"]
            tk = tj - e_["toe_sec"]
            if tk > 43200: tk -= 86400
            if tk < -43200: tk += 86400
            E = float(kepler_E(np.array([e_["M0"] + n * tk]), e_["e"])[0])
            R[j] = e_["e"] * e_["sqrt_a"] * math.sin(E)
        resid = clk - bc
        out[prn] = dict(
            precise_clk_ptp_ns=float(np.ptp(clk) * 1e9),
            broadcast_poly_ptp_ns=float(np.ptp(bc) * 1e9),
            precise_minus_broadcast_ptp_ns=float(np.ptp(resid) * 1e9),
            predicted_dt_rel_ptp_ns=float(np.ptp(F_TRUTH * R) * 1e9))
    return out


# --------------------------------------------------------------------------
# Route 2: per-epoch geometric residual, receiver-clock removal, recovery
# --------------------------------------------------------------------------
def build_residual(prn, obs, sats, dec):
    t, P = obs[prn]
    m = np.argsort(t); t, P = t[m], P[m]
    t = t[::dec]; P = P[::dec]
    tt, resid = geom_residual(prn, t, P, sats[prn], STATION)
    return tt, resid


def estimate_rx_clock(circ, obs, sats, dec):
    """Receiver clock (metres) per epoch from near-circular sats: remove each
    sat's own day-mean (folds ambiguity/hardware bias) then take epoch median."""
    rb = {}
    for p in circ:
        if p not in obs or p not in sats:
            continue
        tt, r = build_residual(p, obs, sats, dec)
        rb[p] = (tt, r)
    grid = np.unique(np.concatenate([rb[p][0] for p in rb]))
    M = np.full((len(rb), len(grid)), np.nan)
    for i, p in enumerate(rb):
        tmap = dict(zip(rb[p][0], rb[p][1]))
        M[i] = [tmap.get(g, np.nan) for g in grid]
    Mc = M - np.nanmean(M, axis=1, keepdims=True)
    rx = np.nanmedian(Mc, axis=0)
    g = np.isfinite(rx)
    return (grid[g], rx[g])   # (epoch grid, receiver-clock metres) for interpolation


def recover(prn, obs, sats, nav, rxmap, dec):
    tt, r = build_residual(prn, obs, sats, dec)
    r = r - np.nanmean(r)                      # remove sat day-mean (bias/ambiguity)
    rxgrid, rxval = rxmap
    # interpolate receiver clock onto this sat's epochs (within grid span only)
    rx = np.interp(tt, rxgrid, rxval, left=np.nan, right=np.nan)
    y = r - rx
    good = np.isfinite(y)
    te = tt[good]; ye = y[good]
    if len(te) < 60:
        return None
    R = dt_rel_reg(nav[prn], te)
    reg = -C_LIGHT * R                          # regressor folds pseudorange sign -> k matches F
    ecc = float(np.median([e["e"] for e in nav[prn]]))

    tn = (te - te.mean()) / (np.ptp(te) + 1e-9)
    Ap = np.vstack([tn**d for d in range(3)]).T
    A = np.column_stack([Ap, reg])
    coef, *_ = np.linalg.lstsq(A, ye, rcond=None)
    k = coef[-1]                                 # s/sqrt(m), compare to F_TRUTH
    pred = A @ coef
    ss_tot = np.sum((ye - ye.mean())**2)
    r2_full = 1 - np.sum((ye - pred)**2) / ss_tot
    coef_p, *_ = np.linalg.lstsq(Ap, ye, rcond=None)
    rp = ye - Ap @ coef_p
    ss_p = np.sum(rp**2)
    incr = 1 - np.sum((ye - pred)**2) / ss_p if ss_p > 0 else 0.0

    # secondary form: cleaned_resid = poly + k2 * (-c * r*rdot)  -> k2 ~ K2_TRUTH
    xs, _ = sp3_interp(sats[prn], te)
    rr = np.linalg.norm(xs, axis=1)
    rdot = np.gradient(rr, te)
    Xrr = -C_LIGHT * rr * rdot
    A2 = np.column_stack([Ap, Xrr])
    coef2, *_ = np.linalg.lstsq(A2, ye, rcond=None)
    k2 = coef2[-1]

    # scramble null -- track BOTH incremental-over-poly AND full-model R^2.
    # (incr-over-poly is unstable when the cubic already absorbs low-freq power;
    #  full-model R^2 vs scrambled full-model R^2 is the robust significance.)
    rng = np.random.default_rng(0)
    nk, nincr, nr2 = [], [], []
    for _ in range(300):
        rs = rng.permutation(reg)
        An = np.column_stack([Ap, rs])
        cf, *_ = np.linalg.lstsq(An, ye, rcond=None)
        nk.append(cf[-1])
        pn = An @ cf
        nincr.append(1 - np.sum((ye - pn)**2) / ss_p if ss_p > 0 else 0.0)
        nr2.append(1 - np.sum((ye - pn)**2) / ss_tot)
    nk = np.array(nk); nincr = np.array(nincr); nr2 = np.array(nr2)

    # time-slide null: roll the regressor (preserves its autocorrelation)
    slide = []
    for sh in range(10, len(reg) - 10, max(1, len(reg) // 40)):
        rs = np.roll(reg, sh)
        An = np.column_stack([Ap, rs])
        cf, *_ = np.linalg.lstsq(An, ye, rcond=None)
        pn = An @ cf
        slide.append(1 - np.sum((ye - pn)**2) / ss_p if ss_p > 0 else 0.0)
    slide = np.array(slide)

    return dict(
        prn=prn, n=int(len(te)), ecc=ecc,
        k=float(k), k_over_truth=float(k / F_TRUTH),
        k2_rrdot=float(k2), k2_over_truth=float(k2 / K2_TRUTH),
        r2_full=float(r2_full), incr_r2=float(incr),
        resid_ptp_m=float(np.ptp(ye)),
        pred_dt_rel_ptp_m=float(np.ptp(C_LIGHT * F_TRUTH * R)),
        corr_resid_pred=float(np.corrcoef(ye, C_LIGHT * F_TRUTH * R)[0, 1]),
        null_scramble_incr_mean=float(nincr.mean()),
        null_scramble_incr_std=float(nincr.std()),
        null_scramble_incr_max=float(nincr.max()),
        null_scramble_r2full_mean=float(nr2.mean()),
        null_scramble_r2full_std=float(nr2.std()),
        null_scramble_r2full_max=float(nr2.max()),
        null_scramble_k_std=float(nk.std()),
        null_timeslide_incr_mean=float(slide.mean()),
        null_timeslide_incr_max=float(slide.max()),
        # robust significance: full-model R^2 vs scrambled full-model R^2
        z_vs_scramble=float((r2_full - nr2.mean()) / (nr2.std() + 1e-12)),
        z_incr_vs_scramble=float((incr - nincr.mean()) / (nincr.std() + 1e-12)),
    )


def main(canary=False):
    nav = parse_rinex_nav(NAV, systems=("G", "E"))
    sats, _ = parse_sp3_full(SP3, systems=("G", "E"))

    r1 = route1_null_check(nav)

    # near-circular Galileo anchors for the receiver clock
    circ = ["E02", "E07", "E10", "E11", "E12", "E19", "E24", "E25", "E31", "E33"]
    if canary:
        targets = ["E18", "E14"]
        dec = 4
    else:
        # eccentric Galileo + highest-ecc GPS
        gps_ecc = sorted([p for p in nav if p.startswith("G")],
                         key=lambda p: -np.median([e["e"] for e in nav[p]]))
        targets = ["E18", "E14"] + gps_ecc[:4]   # G21,G02,G07,G15 ...
        dec = 2

    want = set(circ + targets)
    obs = parse_obs_c1(OBS, want_prns=want)
    circ_have = [p for p in circ if p in obs]
    rxmap = estimate_rx_clock(circ_have, obs, sats, dec)

    per_sat = {}
    for p in targets:
        if p not in obs or p not in sats or p not in nav:
            continue
        res = recover(p, obs, sats, nav, rxmap, dec)
        if res is not None:
            per_sat[p] = res

    kr = np.array([per_sat[p]["k_over_truth"] for p in per_sat])
    incr = np.array([per_sat[p]["incr_r2"] for p in per_sat])
    # headline = eccentric Galileo (largest, cleanest signal)
    gal = [p for p in per_sat if p in ("E14", "E18")]
    summary = dict(
        route="Route 2 -- raw RINEX C1 pseudorange (BRUX) + precise SP3 orbit/clock "
              "(term removed) + known station ECEF; receiver clock from near-circular "
              "sats; recover k in cleaned_resid = poly + k*(-c*e*sqrt(a)*sinE)",
        station="BRUX (Brussels) 4027881.8468 306998.2610 4919498.6524 m",
        day="2023-001",
        F_truth=F_TRUTH, K2_truth=K2_TRUTH,
        n_target_sats=len(per_sat), targets=list(per_sat),
        k_over_truth_mean=float(kr.mean()), k_over_truth_std=float(kr.std()),
        k_over_truth_median=float(np.median(kr)),
        eccentric_galileo=gal,
        gal_k_over_truth=[float(per_sat[p]["k_over_truth"]) for p in gal],
        gal_incr_r2=[float(per_sat[p]["incr_r2"]) for p in gal],
        gal_r2_full=[float(per_sat[p]["r2_full"]) for p in gal],
        gal_z_vs_scramble=[float(per_sat[p]["z_vs_scramble"]) for p in gal],
        gal_k2_over_truth=[float(per_sat[p]["k2_over_truth"]) for p in gal],
        median_incr_r2=float(np.median(incr)),
        median_z_vs_scramble=float(np.median([per_sat[p]["z_vs_scramble"] for p in per_sat])),
        grand_null_scramble_incr=float(np.mean([per_sat[p]["null_scramble_incr_mean"] for p in per_sat])),
        route1_null=dict(
            note="precise SP3 clock and broadcast polynomial agree to a few ns; "
                 "the ~hundreds-of-ns relativistic term is in NEITHER. Route 1 NULL.",
            per_sat=r1),
    )
    out = dict(summary=summary, per_sat=per_sat)
    tag = "_canary" if canary else "_results"
    fn = f"probe_6c_gps_positive{tag}.json"
    with open(fn, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {fn}")
    print("\nper-sat (prn ecc k/truth R2 incrR2 z_null corr resid_ptp pred_ptp):")
    for p in sorted(per_sat, key=lambda p: -per_sat[p]["ecc"]):
        r = per_sat[p]
        print(f"  {p}  e={r['ecc']:.4f}  k/truth={r['k_over_truth']:+.3f}  "
              f"R2={r['r2_full']:.3f}  incrR2={r['incr_r2']:.3f}  "
              f"z={r['z_vs_scramble']:.0f}  corr={r['corr_resid_pred']:+.3f}  "
              f"resid={r['resid_ptp_m']:.0f}m pred={r['pred_dt_rel_ptp_m']:.0f}m")


if __name__ == "__main__":
    main(canary=("--canary" in sys.argv))
