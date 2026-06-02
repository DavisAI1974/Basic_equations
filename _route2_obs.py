"""
Route 2 core: recover relativistic clock term from RAW RINEX pseudorange
observations + precise SP3 orbit/clock + known station ECEF.

The pseudorange model (single frequency C1C, in metres):
    P = rho_geom + c*(dt_recv - dt_sat_precise) + tropo + iono + c*dt_rel + eps
The precise SP3 satellite clock has the relativistic eccentricity term REMOVED
(standard IGS convention; confirmed null in the SP3-only probe). So if we use
dt_sat_precise to correct P, the leftover satellite-specific periodic term is
    c * dt_rel = c * F * e*sqrt(a)*sin(E).
Receiver clock dt_recv is common to all sats at an epoch -> remove it by
subtracting the per-epoch across-satellite mean of the residual. What survives
for the eccentric sat is (mostly) c*dt_rel. Fit resid = poly + k*(e*sqrt(a)*sinE).
"""
import numpy as np, math

C_LIGHT = 299792458.0
OMEGA_E = 7.2921151467e-5
MU = 3.986004418e14


# ---- SP3 full parse keeping ALL epochs incl. day-internal absolute seconds ----
def parse_sp3_full(path, systems=("G", "E")):
    epochs = []
    recs = {}
    ei = -1
    with open(path) as f:
        for line in f:
            if line.startswith("*"):
                p = line.split()
                yr, mo, dy = int(p[1]), int(p[2]), int(p[3])
                hh, mm = int(p[4]), int(p[5]); ss = float(p[6])
                tsec = ((dy - 1) * 24 + hh) * 3600.0 + mm * 60.0 + ss
                epochs.append(tsec); ei += 1
            elif line.startswith("P") and len(line) > 4 and line[1] in systems:
                prn = line[1:4]
                try:
                    x = float(line[4:18]); y = float(line[18:32]); z = float(line[32:46]); clk = float(line[46:60])
                except ValueError:
                    continue
                recs.setdefault(prn, []).append((ei, x, y, z, clk))
    et = np.array(epochs)
    sats = {}
    for prn, rows in recs.items():
        idx = np.array([r[0] for r in rows])
        xyz = np.array([[r[1], r[2], r[3]] for r in rows]) * 1000.0
        clk = np.array([r[4] for r in rows])
        good = clk < 999999.0
        if good.sum() < 50:
            continue
        sats[prn] = dict(t=et[idx][good], xyz=xyz[good], clk=clk[good] * 1e-6)
    return sats, et


# ---- RINEX 3 obs parser: extract C1C pseudorange per sat per epoch ----
def parse_obs_c1(path, want_prns):
    """Return dict prn -> (t_sec[N], P[N]) for the FIRST pseudorange (C1x) obs.
    Builds per-system obs-type maps from the header; picks the first C1* code."""
    obs_types = {}
    c1_idx = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    # header
    while "END OF HEADER" not in lines[i]:
        L = lines[i]
        if "SYS / # / OBS TYPES" in L:
            sysid = L[0]
            ntypes = int(L[3:6])
            types = L[7:60].split()
            j = i
            while len(types) < ntypes:
                j += 1
                types += lines[j][7:60].split()
            obs_types[sysid] = types
            i = j
        i += 1
    i += 1
    for s, ts in obs_types.items():
        idx = next((k for k, v in enumerate(ts) if v.startswith("C1")), None)
        c1_idx[s] = idx
    data = {p: ([], []) for p in want_prns}
    N = len(lines)
    while i < N:
        L = lines[i]
        if L.startswith(">"):
            p = L.split()
            yr, mo, dy = int(p[1]), int(p[2]), int(p[3])
            hh, mm = int(p[4]), int(p[5]); ss = float(p[6])
            nsat = int(L[32:35])   # cols 33-35 hold # of sats; p[7]=epoch flag, p[8]=nsat
            tsec = ((dy - 1) * 24 + hh) * 3600.0 + mm * 60.0 + ss
            for r in range(nsat):
                i += 1
                row = lines[i]
                prn = row[0:3]
                if prn not in want_prns:
                    continue
                s = prn[0]
                k = c1_idx.get(s)
                if k is None:
                    continue
                field = row[3 + 16 * k: 3 + 16 * k + 14]
                try:
                    val = float(field)
                except ValueError:
                    continue
                if val == 0.0:
                    continue
                data[prn][0].append(tsec); data[prn][1].append(val)
            i += 1
        else:
            i += 1
    return {p: (np.array(v[0]), np.array(v[1])) for p, v in data.items() if len(v[0]) > 50}


def sp3_interp(sats_prn, tq):
    """Lagrange/poly interp of SP3 xyz + clk to query times tq. Uses sliding
    10-point polynomial windows (standard for SP3 position interpolation)."""
    t = sats_prn["t"]; xyz = sats_prn["xyz"]; clk = sats_prn["clk"]
    xq = np.empty((len(tq), 3)); cq = np.empty(len(tq))
    deg = 10
    for j, tj in enumerate(tq):
        c = np.searchsorted(t, tj)
        lo = max(0, c - deg // 2); hi = min(len(t), lo + deg); lo = max(0, hi - deg)
        ts = t[lo:hi]
        # barycentric-ish via numpy polyfit on centered time (stable enough for 10 pts over 50 min)
        tc = ts - tj
        for d in range(3):
            xq[j, d] = np.polyval(np.polyfit(tc, xyz[lo:hi, d], min(deg - 1, len(ts) - 1)), 0.0)
        cq[j] = np.polyval(np.polyfit(tc, clk[lo:hi], min(deg - 1, len(ts) - 1)), 0.0)
    return xq, cq


def geom_residual(prn, obs_t, P, sats_prn, station_xyz):
    """Residual = P - rho_geom + c*dt_sat_precise, with travel-time + Sagnac.
    Returns (t, resid_metres, valid_mask)."""
    # initial guess travel time ~0.07 s; iterate transmit time
    tx = obs_t - 0.07
    for _ in range(3):
        xs, cs = sp3_interp(sats_prn, tx)
        # Sagnac: rotate sat ECEF by Earth rotation over travel time
        dt_travel = obs_t - tx
        theta = OMEGA_E * dt_travel
        xr = np.empty_like(xs)
        ct = np.cos(theta); st = np.sin(theta)
        xr[:, 0] = ct * xs[:, 0] + st * xs[:, 1]
        xr[:, 1] = -st * xs[:, 0] + ct * xs[:, 1]
        xr[:, 2] = xs[:, 2]
        rho = np.linalg.norm(xr - station_xyz, axis=1)
        tx = obs_t - (rho / C_LIGHT)
    xs, cs = sp3_interp(sats_prn, tx)
    # residual in metres: P - rho + c*dt_sat (dt_recv stays in)
    resid = P - rho + C_LIGHT * cs
    return obs_t, resid


def dt_rel_reg(eph_list, t):
    """e*sqrt(a)*sin(E) regressor (sqrt(m)) from broadcast Keplerian elems."""
    from probe_6c_gps_positive import kepler_E
    tocs = np.array([e["toc_sec"] for e in eph_list])
    R = np.zeros_like(t)
    for j, tj in enumerate(t):
        k = int(np.argmin(np.abs(tocs - tj)))
        e_ = eph_list[k]
        a = e_["sqrt_a"]**2
        n = math.sqrt(MU / a**3) + e_["dn"]
        tk = tj - e_["toe_sec"]
        if tk > 43200: tk -= 86400
        if tk < -43200: tk += 86400
        M = e_["M0"] + n * tk
        E = float(kepler_E(np.array([M]), e_["e"])[0])
        R[j] = e_["e"] * e_["sqrt_a"] * math.sin(E)
    return R
