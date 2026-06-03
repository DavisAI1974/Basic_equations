"""
STRENGTHEN the S18 GPS time-dilation dipole recovery -- multi-station, multi-day,
dual-frequency.

CONTEXT (S18, INFO-059): we recovered the GR relativistic clock coefficient -2/c^2
from RAW GPS data and showed two "other" dipole forms travel to it where the
windowed-entropy FLOW dipole was blind:
  (A) RAW cross-covariance dipole: eccentric Galileo E18 lag-0 |cc|~0.99, z~37.
  (B) STATIC ALGEBRAIC dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2 with a
      CIRCULAR-SHIFT tautology-killing null: E18 R2 0.989, excess +0.571, z=5.7.
But that rested on ONE station (BRUX), ONE day (2023 DOY 001), ~10 algebraic
windows -> thin stats.

THIS PROBE strengthens it to MULTIPLE (station, day) replicas and adds the
DUAL-FREQUENCY ionosphere-free combination P3 = (f1^2*C1 - f2^2*C2)/(f1^2-f2^2)
to try to rescue the near-circular GPS sats whose tiny relativistic signal
(17-33 m) was swamped by single-frequency C1 ionosphere in S18.

Channels (per station, day, sat), exactly as INFO-059:
  A = cleaned RINEX pseudorange relativistic residual
      (P - geometry(SP3 orbit + station ECEF + Sagnac + travel-time)
         + c*dt_sat_precise(SP3, term removed) - receiver clock(near-circular)).
      Single-freq (C1) OR dual-freq ionosphere-free (P3) per --combo.
  B = independent broadcast-element geometry driver c*F*e*sqrt(a)*sin(E(t))
      from broadcast Keplerian elements (different data source -> non-circular).

For each (station, day, sat) we run BOTH dipole forms:
  (A) raw cross-covariance over lag (raw_cov_dipole, reused from
      probe_time_dipole_gpspulsar.py): lag-0 |cc| vs scramble null.
  (B) static algebraic dipole with circular-shift tautology null
      (algebraic_dipole_with_null, reused): R2_real, shift-null mean, excess, z.

Aggregate: mean +/- std of the algebraic excess-over-null and z across the
(station, day) replicas, SEPARATELY for eccentric Galileo (signal) vs near-circular
(control). The question: does E18's algebraic excess (+0.571, z=5.7) hold up as a
MEAN across many replicas, and does dual-freq rescue the circular controls (do they
STAY low-excess = good control, or now ALSO show the coupling once iono is removed)?

HONESTY: NO novelty. This RECOVERS the known GR -2/c^2 coefficient -- a positive
control that the dipole tool travels, NOT new physics. Each replica is one data
point; the tautology-null excess is reported honestly per replica. If a station/day
won't fetch, it is skipped and reported -- nothing is fabricated. If dual-freq does
NOT rescue the controls, that is reported.

Run:  python probe_time_dipole_gps_strengthen.py --canary    (1 station/1 day)
      python probe_time_dipole_gps_strengthen.py             (full multi/multi)
      python probe_time_dipole_gps_strengthen.py --combo C1   (single-freq only)
      python probe_time_dipole_gps_strengthen.py --combo both (default: C1 + P3)
"""
import os
import sys
import gzip
import json
import math
import time
import argparse
import urllib.request
import numpy as np

from _route2_obs import (parse_sp3_full, geom_residual, dt_rel_reg, sp3_interp,
                         C_LIGHT, MU)
from probe_6c_gps_positive import parse_rinex_nav, F_TRUTH, K2_TRUTH
from probe_time_dipole_gpspulsar import raw_cov_dipole, algebraic_dipole_with_null

DATA = "data/gps"
os.makedirs(DATA, exist_ok=True)

# GNSS L-band frequencies (Hz) for the ionosphere-free combination
F_G1 = 1575.42e6      # GPS L1
F_G2 = 1227.60e6      # GPS L2
F_E1 = 1575.420e6     # Galileo E1
F_E5a = 1176.450e6    # Galileo E5a (the C5Q code in this RINEX)

# stations (IGS 9-char + ECEF read from RINEX header; we read header to be safe)
STATIONS = ["BRUX00BEL", "ONSA00SWE", "WTZR00DEU", "GRAZ00AUT", "MATE00ITA", "ALGO00CAN"]
DAYS = ["001", "002", "003"]
GPS_WEEK = "2243"     # 2023 DOY 001-003

BKG = "https://igs.bkg.bund.de/root_ftp/IGS"
BKG_NAV = "https://igs.bkg.bund.de/root_ftp/IGS/BRDC"
AIUB = "http://ftp.aiub.unibe.ch/CODE_MGEX/CODE/2023"

# near-circular Galileo anchors for the receiver-clock estimate (e < 0.03)
CIRC_E = ["E02", "E07", "E10", "E11", "E12", "E19", "E24", "E25", "E31", "E33"]
ECC_GAL = ["E18", "E14"]          # the eccentric Galileo SIGNAL sats (e~0.162)


# --------------------------------------------------------------------------
# fetch helpers (skip+report on failure; never fabricate)
# --------------------------------------------------------------------------
def _download(url, dst, timeout=120):
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        return True
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "od-gps/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        if len(data) < 1000:
            return False
        with open(dst, "wb") as f:
            f.write(data)
        return True
    except Exception as e:
        print(f"   [fetch fail] {url}\n      {type(e).__name__}: {e}")
        return False


def _gunzip(src_gz, dst):
    if os.path.exists(dst) and os.path.getsize(dst) > 1000:
        return True
    try:
        with gzip.open(src_gz, "rb") as fi, open(dst, "wb") as fo:
            fo.write(fi.read())
        return True
    except Exception as e:
        print(f"   [gunzip fail] {src_gz}: {e}")
        return False


def fetch_obs(station, day):
    """Hatanaka .crx.gz -> .crx -> RINEX .rnx (via hatanaka decoder).
    Returns local .rnx path or None."""
    rnx = f"{DATA}/{station}_2023{day}.rnx"
    if os.path.exists(rnx) and os.path.getsize(rnx) > 100000:
        return rnx
    crxgz = f"{DATA}/{station}_2023{day}.crx.gz"
    url = f"{BKG}/obs/2023/{day}/{station}_R_2023{day}0000_01D_30S_MO.crx.gz"
    if not _download(url, crxgz):
        return None
    # hatanaka decode (handles .crx.gz directly)
    try:
        import hatanaka
        with open(crxgz, "rb") as f:
            raw = f.read()
        txt = hatanaka.decompress(raw)
        if isinstance(txt, bytes):
            txt = txt.decode("ascii", "replace")
        with open(rnx, "w") as f:
            f.write(txt)
        return rnx
    except Exception as e:
        print(f"   [hatanaka fail] {station} {day}: {type(e).__name__}: {e}")
        return None


def fetch_nav(day):
    rnx = f"{DATA}/BRDC_2023{day}.rnx"
    if os.path.exists(rnx) and os.path.getsize(rnx) > 100000:
        return rnx
    gz = f"{DATA}/BRDC_2023{day}.rnx.gz"
    url = f"{BKG_NAV}/2023/{day}/BRDC00WRD_R_2023{day}0000_01D_MN.rnx.gz"
    if not _download(url, gz):
        return None
    return rnx if _gunzip(gz, rnx) else None


def fetch_sp3(day):
    sp3 = f"{DATA}/COD_2023{day}.SP3"
    if os.path.exists(sp3) and os.path.getsize(sp3) > 100000:
        return sp3
    gz = f"{DATA}/COD_2023{day}.SP3.gz"
    url = f"{AIUB}/COD0MGXFIN_2023{day}0000_01D_05M_ORB.SP3.gz"
    if not _download(url, gz):
        return None
    return sp3 if _gunzip(gz, sp3) else None


# --------------------------------------------------------------------------
# RINEX 3 obs parser -- single (C1) OR dual-freq ionosphere-free (P3)
# --------------------------------------------------------------------------
def _station_ecef(path):
    with open(path) as f:
        for line in f:
            if "APPROX POSITION XYZ" in line:
                return np.array([float(line[0:14]), float(line[14:28]), float(line[28:42])])
            if "END OF HEADER" in line:
                break
    return None


def parse_obs_dualfreq(path, want_prns, combo):
    """Parse RINEX 3 OBS. For each sat per epoch build the chosen pseudorange:
        combo='C1' : the first C1* code (matches S18 single-frequency).
        combo='P3' : ionosphere-free (f1^2*C_f1 - f2^2*C_f2)/(f1^2 - f2^2),
                     using the two best code observables for that system.
    Returns dict prn -> (t_sec[N], P[N]).
    GPS dual: C1W/C2W (P-code) preferred, else C1C/C2W or C1C/C2L.
    Galileo dual: C1C (E1) / C5Q (E5a)."""
    with open(path) as f:
        lines = f.readlines()
    obs_types = {}
    i = 0
    while i < len(lines) and "END OF HEADER" not in lines[i]:
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

    # choose code indices + frequencies per system
    def pick(sysid):
        ts = obs_types.get(sysid, [])
        def idx(code):
            return ts.index(code) if code in ts else None
        if combo == "C1":
            k = next((n for n, v in enumerate(ts) if v.startswith("C1")), None)
            return ("single", k, None, None, None)
        # P3 dual
        if sysid == "G":
            # prefer P-code C1W/C2W; fall back to C1C/C2W, C1C/C2L
            for c1, c2 in (("C1W", "C2W"), ("C1C", "C2W"), ("C1C", "C2L"),
                           ("C1W", "C2L"), ("C1C", "C2C")):
                a, b = idx(c1), idx(c2)
                if a is not None and b is not None:
                    return ("dual", a, b, F_G1, F_G2)
        elif sysid == "E":
            for c1, c2 in (("C1C", "C5Q"), ("C1C", "C8Q"), ("C1C", "C7Q")):
                a, b = idx(c1), idx(c2)
                if a is not None and b is not None:
                    return ("dual", a, b, F_E1, F_E5a)
        # no dual available -> fall back to single C1
        k = next((n for n, v in enumerate(ts) if v.startswith("C1")), None)
        return ("single", k, None, None, None)

    pickmap = {s: pick(s) for s in obs_types}
    data = {p: ([], []) for p in want_prns}
    N = len(lines)
    while i < N:
        L = lines[i]
        if L.startswith(">"):
            p = L.split()
            dy = int(p[3]); hh = int(p[4]); mm = int(p[5]); ss = float(p[6])
            nsat = int(L[32:35])
            tsec = ((dy - 1) * 24 + hh) * 3600.0 + mm * 60.0 + ss
            for _ in range(nsat):
                i += 1
                row = lines[i]
                prn = row[0:3]
                if prn not in want_prns:
                    continue
                s = prn[0]
                spec = pickmap.get(s)
                if spec is None:
                    continue
                mode, ka, kb, f1, f2 = spec
                if mode == "single":
                    if ka is None:
                        continue
                    try:
                        v = float(row[3 + 16 * ka: 3 + 16 * ka + 14])
                    except ValueError:
                        continue
                    if v == 0.0:
                        continue
                    P = v
                else:
                    try:
                        c1 = float(row[3 + 16 * ka: 3 + 16 * ka + 14])
                        c2 = float(row[3 + 16 * kb: 3 + 16 * kb + 14])
                    except (ValueError, IndexError):
                        continue
                    if c1 == 0.0 or c2 == 0.0:
                        continue
                    P = (f1 * f1 * c1 - f2 * f2 * c2) / (f1 * f1 - f2 * f2)
                data[prn][0].append(tsec)
                data[prn][1].append(P)
            i += 1
        else:
            i += 1
    return {p: (np.array(v[0]), np.array(v[1])) for p, v in data.items() if len(v[0]) > 50}


# --------------------------------------------------------------------------
# residual + receiver clock (mirror probe_6c_gps_positive, station-parametrized)
# --------------------------------------------------------------------------
def build_residual(prn, obs, sats_prn, station, dec):
    t, P = obs[prn]
    m = np.argsort(t); t, P = t[m], P[m]
    t = t[::dec]; P = P[::dec]
    return geom_residual(prn, t, P, sats_prn, station)


def estimate_rx_clock(circ, obs, sats, station, dec):
    rb = {}
    for p in circ:
        if p not in obs or p not in sats:
            continue
        tt, r = build_residual(p, obs, sats[p], station, dec)
        rb[p] = (tt, r)
    if not rb:
        return None
    grid = np.unique(np.concatenate([rb[p][0] for p in rb]))
    M = np.full((len(rb), len(grid)), np.nan)
    for i, p in enumerate(rb):
        tmap = dict(zip(rb[p][0], rb[p][1]))
        M[i] = [tmap.get(g, np.nan) for g in grid]
    Mc = M - np.nanmean(M, axis=1, keepdims=True)
    rx = np.nanmedian(Mc, axis=0)
    g = np.isfinite(rx)
    return (grid[g], rx[g])


def build_channels(prn, obs, sats, nav, rxmap, station, dec):
    """A = cleaned pseudorange residual [m]; B = c*F*e*sqrt(a)*sin(E) [m]."""
    tt, r = build_residual(prn, obs, sats[prn], station, dec)
    r = r - np.nanmean(r)
    rxgrid, rxval = rxmap
    rx = np.interp(tt, rxgrid, rxval, left=np.nan, right=np.nan)
    A = r - rx
    good = np.isfinite(A)
    te = tt[good]; A = A[good]
    if len(te) < 60:
        return None
    R = dt_rel_reg(nav[prn], te)
    B = C_LIGHT * F_TRUTH * R
    ecc = float(np.median([e["e"] for e in nav[prn]]))
    return te, A, B, ecc


def recover_k(te, A, B):
    """Single-lag linear recovery of k (matches the S18 coefficient recovery),
    reported alongside the dipoles so the coefficient travels too."""
    # regressor reg = -c*e*sqrt(a)*sin(E) = -B/F_TRUTH  (so k matches F_TRUTH)
    reg = -B / (C_LIGHT * F_TRUTH) * F_TRUTH    # == -B/C_LIGHT ; keep explicit
    reg = -B / C_LIGHT * 1.0                      # metres of dt_rel*c -> dt form
    # match probe_6c: reg = -C_LIGHT * R, R = B/(C_LIGHT*F_TRUTH); => reg = -B/F_TRUTH
    reg = -B / F_TRUTH
    tn = (te - te.mean()) / (np.ptp(te) + 1e-9)
    Ap = np.vstack([tn ** d for d in range(3)]).T
    Afull = np.column_stack([Ap, reg])
    coef, *_ = np.linalg.lstsq(Afull, A, rcond=None)
    k = coef[-1]
    pred = Afull @ coef
    ss_tot = np.sum((A - A.mean()) ** 2)
    r2 = 1 - np.sum((A - pred) ** 2) / ss_tot if ss_tot > 0 else float("nan")
    corr = float(np.corrcoef(A, B)[0, 1])
    return float(k), float(k / F_TRUTH), float(r2), corr


# --------------------------------------------------------------------------
def run_one(station9, day, combo, sats, nav, dec, win, stride, n_scram, n_null):
    """One (station, day) replica for the chosen pseudorange combo.
    Returns dict prn -> result for ECC_GAL + a near-circular GPS control set."""
    obs_path = fetch_obs(station9, day)
    if obs_path is None:
        return None, "obs_fetch_failed"
    station_ecef = _station_ecef(obs_path)
    if station_ecef is None:
        return None, "no_station_ecef"

    # control = near-circular GPS sats (tiny ~17-33 m signal; the dual-freq rescue target)
    gps_circ = sorted([p for p in nav if p.startswith("G")],
                      key=lambda p: np.median([e["e"] for e in nav[p]]))[:5]
    want = set(CIRC_E + ECC_GAL + gps_circ)
    obs = parse_obs_dualfreq(obs_path, want, combo)
    if not obs:
        return None, "obs_parse_empty"
    circ_have = [p for p in CIRC_E if p in obs]
    rxmap = estimate_rx_clock(circ_have, obs, sats, station_ecef, dec)
    if rxmap is None:
        return None, "no_rx_clock"

    out = {}
    for p in ECC_GAL + gps_circ:
        if p not in obs or p not in sats or p not in nav:
            continue
        ch = build_channels(p, obs, sats, nav, rxmap, station_ecef, dec)
        if ch is None:
            continue
        te, A, B, ecc = ch
        try:
            k, kot, kr2, corr = recover_k(te, A, B)
        except Exception:
            k = kot = kr2 = float("nan"); corr = float("nan")
        raw = raw_cov_dipole(te, A, B, n_scramble=n_scram, seed=11)
        # adaptive algebraic window: size to THIS sat's epoch count so stations
        # with shorter tracks still yield >=8 windows (the windowed entropy dipole
        # needs enough windows). target ~14 windows at 50% overlap; floor at the
        # S18-canary window so BRUX/E18 stays comparable to the anchor.
        n_ep = len(te)
        w_target = max(win // 3, min(win, int(n_ep / 7.5)))   # ~14 windows @ stride=w/2
        s_target = max(3, w_target // 2)
        alg = algebraic_dipole_with_null(A, B, w_target, s_target, n_null, seed=22)
        # also run the FIXED S18 window for direct anchor comparison
        alg_fixed = algebraic_dipole_with_null(A, B, win, stride, n_null, seed=22)
        out[p] = {
            "prn": p, "ecc": ecc, "n_epochs": int(len(te)),
            "signal_ptp_m": float(np.ptp(B)),
            "is_control": bool(p.startswith("G")),
            "k_over_truth": kot, "k_recov": k, "regression_r2": kr2,
            "corr_A_B": corr,
            "alg_window_used": int(w_target), "alg_stride_used": int(s_target),
            "raw_covariance_dipole": raw,
            "static_algebraic_dipole": alg,
            "static_algebraic_dipole_fixed_s18win": alg_fixed,
        }
    return out, "ok"


def aggregate(replicas, predicate, field_path):
    """Collect a nested field across all sats matching predicate, over replicas."""
    vals = []
    for rep in replicas:
        for p, d in rep["per_sat"].items():
            if not predicate(d):
                continue
            node = d
            ok = True
            for k in field_path:
                if node is None or (isinstance(node, dict) and k not in node):
                    ok = False; break
                node = node[k]
            if ok and node is not None and not (isinstance(node, float) and math.isnan(node)):
                vals.append(float(node))
    if not vals:
        return {"n": 0, "mean": None, "std": None, "min": None, "max": None}
    a = np.array(vals)
    return {"n": int(len(a)), "mean": float(a.mean()), "std": float(a.std()),
            "min": float(a.min()), "max": float(a.max())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true",
                    help="1 station (BRUX) / 1 day (001) only")
    ap.add_argument("--combo", choices=["C1", "P3", "both"], default="both",
                    help="pseudorange: single-freq C1, dual-freq P3, or both")
    args = ap.parse_args()
    t0 = time.time()

    stations = ["BRUX00BEL"] if args.canary else STATIONS
    days = ["001"] if args.canary else DAYS
    combos = ["C1"] if args.combo == "C1" else (["P3"] if args.combo == "P3"
                                                else ["C1", "P3"])
    dec = 4 if args.canary else 2
    win = 30 if args.canary else 90
    stride = win // 2
    n_scram = 100 if args.canary else 300
    n_null = 100 if args.canary else 400

    # pre-fetch NAV + SP3 per day (shared across stations)
    nav_by_day, sp3_by_day, skipped = {}, {}, []
    for day in days:
        navp = fetch_nav(day)
        sp3p = fetch_sp3(day)
        if navp is None or sp3p is None:
            skipped.append({"scope": f"day {day}", "reason":
                            f"nav={'ok' if navp else 'FAIL'} sp3={'ok' if sp3p else 'FAIL'}"})
            continue
        print(f"[parse] day {day} nav+sp3 ...")
        nav_by_day[day] = parse_rinex_nav(navp, systems=("G", "E"))
        sp3_by_day[day] = parse_sp3_full(sp3p, systems=("G", "E"))[0]

    replicas_by_combo = {c: [] for c in combos}
    for combo in combos:
        for day in days:
            if day not in nav_by_day:
                continue
            nav = nav_by_day[day]; sats = sp3_by_day[day]
            for st in stations:
                tag = f"{st[:4]}/{day}/{combo}"
                print(f"[run] {tag} ...", flush=True)
                per_sat, status = run_one(st, day, combo, sats, nav, dec,
                                          win, stride, n_scram, n_null)
                if per_sat is None or not per_sat:
                    skipped.append({"scope": tag, "reason": status})
                    print(f"   skipped: {status}")
                    continue
                replicas_by_combo[combo].append(
                    {"station": st[:4], "day": day, "combo": combo,
                     "per_sat": per_sat})
                # per-replica console line for the eccentric Galileo signal sats
                for p in ECC_GAL:
                    if p in per_sat:
                        d = per_sat[p]
                        alg = d["static_algebraic_dipole"] or {}
                        raw = d["raw_covariance_dipole"]
                        print(f"   {p} ECC e={d['ecc']:.3f} sig={d['signal_ptp_m']:.0f}m "
                              f"k/tru={d['k_over_truth']:+.2f} | ALG R2={alg.get('R2_real', float('nan')):.3f} "
                              f"excess={alg.get('excess_over_null', float('nan')):+.3f} "
                              f"z={alg.get('z_over_null', float('nan')):.1f} | "
                              f"RAW |cc0|={abs(raw['cc_at_lag0']):.3f} z={raw['z_over_scramble']:.0f}")

    # ---- aggregate stats per combo, eccentric (signal) vs circular (control) ----
    def is_ecc(d):  return not d["is_control"]
    def is_ctrl(d): return d["is_control"]

    agg = {}
    for combo, reps in replicas_by_combo.items():
        agg[combo] = {
            "n_replicas_station_day": len(reps),
            "n_sat_replicas_eccentric": sum(1 for r in reps for d in r["per_sat"].values() if is_ecc(d)),
            "n_sat_replicas_circular": sum(1 for r in reps for d in r["per_sat"].values() if is_ctrl(d)),
            "eccentric_signal": {
                "alg_excess_over_null": aggregate(reps, is_ecc, ["static_algebraic_dipole", "excess_over_null"]),
                "alg_z_over_null": aggregate(reps, is_ecc, ["static_algebraic_dipole", "z_over_null"]),
                "alg_R2_real": aggregate(reps, is_ecc, ["static_algebraic_dipole", "R2_real"]),
                "raw_cc_lag0_abs": aggregate(reps, is_ecc, ["raw_covariance_dipole", "cc_at_lag0"]),
                "raw_z_over_scramble": aggregate(reps, is_ecc, ["raw_covariance_dipole", "z_over_scramble"]),
                "k_over_truth": aggregate(reps, is_ecc, ["k_over_truth"]),
            },
            "circular_control": {
                "alg_excess_over_null": aggregate(reps, is_ctrl, ["static_algebraic_dipole", "excess_over_null"]),
                "alg_z_over_null": aggregate(reps, is_ctrl, ["static_algebraic_dipole", "z_over_null"]),
                "alg_R2_real": aggregate(reps, is_ctrl, ["static_algebraic_dipole", "R2_real"]),
                "raw_cc_lag0_abs": aggregate(reps, is_ctrl, ["raw_covariance_dipole", "cc_at_lag0"]),
                "raw_z_over_scramble": aggregate(reps, is_ctrl, ["raw_covariance_dipole", "z_over_scramble"]),
                "k_over_truth": aggregate(reps, is_ctrl, ["k_over_truth"]),
            },
        }

    result = {
        "probe": "STRENGTHEN the S18 GPS time-dilation dipole recovery -- "
                 "multi-station, multi-day, dual-frequency.",
        "honesty": "RECOVERY of the KNOWN GR -2/c^2 coefficient -- a positive "
                   "control that the dipole tool travels, NOT new physics. Each "
                   "(station,day,sat) is one data point; the tautology-null excess "
                   "is reported as-is per replica; failed fetches are skipped + "
                   "reported, never fabricated.",
        "s18_anchor": {"station": "BRUX", "day": "2023-001",
                       "E18_alg_excess": 0.571, "E18_alg_z": 5.7,
                       "E18_raw_cc_lag0": 0.99, "E18_raw_z": 37},
        "config": {
            "stations": [s[:4] for s in stations], "days": days, "combos": combos,
            "decimation": dec, "alg_window_epochs": win, "alg_stride_epochs": stride,
            "n_scramble": n_scram, "n_shift_null": n_null,
            "channel_A": "cleaned pseudorange relativistic residual [m] (C1 or "
                         "dual-freq P3 ionosphere-free)",
            "channel_B": "c*F*e*sqrt(a)*sin(E) broadcast-element driver [m] "
                         "(independent source -> non-circular)",
            "F_truth": F_TRUTH, "K2_truth": K2_TRUTH,
            "eccentric_galileo_signal": ECC_GAL,
            "circular_control": "5 lowest-e GPS sats per day (dual-freq rescue target)",
            "nav": "BRDC00WRD multi-GNSS broadcast (BKG)",
            "sp3": "CODE MGEX final COD0MGXFIN 05M ORB (AIUB)",
            "obs": "IGS station RINEX3 30S (BKG, Hatanaka)",
        },
        "aggregate": agg,
        "replicas": replicas_by_combo,
        "skipped": skipped,
        "runtime_s": round(time.time() - t0, 1),
    }

    tag = "_canary" if args.canary else "_results"
    out = f"probe_time_dipole_gps_strengthen{tag}.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)

    # ---- console summary ----
    print("\n" + "=" * 72)
    print(f"GPS STRENGTHEN  ({'CANARY' if args.canary else 'FULL'})  "
          f"{result['runtime_s']}s")
    print("=" * 72)
    for combo in combos:
        a = agg[combo]
        print(f"\n[combo {combo}]  replicas(station,day)={a['n_replicas_station_day']}  "
              f"sat-replicas: ecc={a['n_sat_replicas_eccentric']} circ={a['n_sat_replicas_circular']}")
        e = a["eccentric_signal"]; c = a["circular_control"]
        def fmt(s):
            if s["mean"] is None:
                return "n=0"
            return f"n={s['n']} mean={s['mean']:+.3f}+/-{s['std']:.3f} [{s['min']:+.2f},{s['max']:+.2f}]"
        print(f"  ECC  alg excess   : {fmt(e['alg_excess_over_null'])}")
        print(f"  ECC  alg z        : {fmt(e['alg_z_over_null'])}")
        print(f"  ECC  raw |cc0|    : {fmt(e['raw_cc_lag0_abs'])}")
        print(f"  ECC  raw z(scram) : {fmt(e['raw_z_over_scramble'])}")
        print(f"  ECC  k/truth      : {fmt(e['k_over_truth'])}")
        print(f"  CTRL alg excess   : {fmt(c['alg_excess_over_null'])}")
        print(f"  CTRL alg z        : {fmt(c['alg_z_over_null'])}")
        print(f"  CTRL raw |cc0|    : {fmt(c['raw_cc_lag0_abs'])}")
        print(f"  CTRL k/truth      : {fmt(c['k_over_truth'])}")
    if skipped:
        print(f"\nskipped ({len(skipped)}):")
        for s in skipped:
            print(f"  {s['scope']}: {s['reason']}")
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
