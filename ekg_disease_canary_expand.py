"""Expansion canary: 4 more disease families beyond Round 1.

Adjusted SUPERWINDOW/SUPERSTEP for shorter (~30 min) arrhythmia records.
Tests whether the 5-of-6 substrate cluster holds when more pathologies
are added; whether AF stays as the outlier or other arrhythmias also
break the cluster.
"""
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import wfdb
import wfdb.processing as wp
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

WINDOW_BEATS = 25      # smaller window for shorter records
STEP_BEATS = 10
SUPERWINDOW = 30
SUPERSTEP = 6
RR_MIN_S = 0.3
RR_MAX_S = 2.0
N_BINS_H = 16

TARGETS = [
    ("supraventricular_arr", "svdb", "800"),
    ("ventricular_malignant", "vfdb", "418"),
    ("arrhythmia_mixed", "mitdb", "100"),
    ("ptb_myocardial_inf", "ptbdb", "patient001/s0010_re"),
    ("ltaf_long_term", "ltafdb", "00"),
    ("driver_stress_v2", "drivedb", "drive01"),
]


def shannon_entropy(x, bins=N_BINS_H):
    x = np.asarray(x, float)
    if x.size < 4 or np.allclose(x.std(), 0.0):
        return 0.0
    lo, hi = np.percentile(x, [1.0, 99.0])
    if hi - lo < 1e-12:
        return 0.0
    hist, _ = np.histogram(x, bins=bins, range=(lo, hi))
    p = hist.astype(float) / hist.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def mutual_info(x, y, bins=N_BINS_H):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 4 or x.std() < 1e-12 or y.std() < 1e-12:
        return 0.0
    lo_x, hi_x = np.percentile(x, [1.0, 99.0])
    lo_y, hi_y = np.percentile(y, [1.0, 99.0])
    if hi_x - lo_x < 1e-12 or hi_y - lo_y < 1e-12:
        return 0.0
    Hxy, _, _ = np.histogram2d(x, y, bins=bins,
                               range=[[lo_x, hi_x], [lo_y, hi_y]])
    Pxy = Hxy / Hxy.sum()
    Px = Pxy.sum(axis=1, keepdims=True)
    Py = Pxy.sum(axis=0, keepdims=True)
    mask = (Pxy > 0) & (Px > 0) & (Py > 0)
    return float(np.sum(Pxy[mask] * np.log(Pxy[mask] / (Px * Py)[mask])))


def extract_ecg(family, db, rec_name):
    t0 = time.time()
    try:
        rec = wfdb.rdrecord(rec_name, pn_dir=db)
        fs = float(rec.fs)
        ecg_idx = None
        for i, n in enumerate(rec.sig_name):
            nm = n.upper()
            if "ECG" in nm or "MLII" in nm or nm.startswith("II") or nm == "V5":
                ecg_idx = i
                break
        if ecg_idx is None:
            ecg_idx = 0
        ecg = rec.p_signal[:, ecg_idx]
        nans = np.isnan(ecg)
        if nans.any():
            ecg = np.where(nans, 0.0, ecg)
        # Handle drivedb's high fs and resample to 250 Hz for xqrs
        if fs > 1000:
            from scipy.signal import resample_poly
            up, down = 250, int(round(fs))
            from math import gcd
            g = gcd(up, down)
            ecg = resample_poly(ecg, up // g, down // g)
            fs = 250.0

        cap_samples = int(2 * 3600 * fs)
        if len(ecg) > cap_samples:
            ecg = ecg[:cap_samples]

        qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
        rr_all = np.diff(qrs) / fs
        rr = rr_all[(rr_all > RR_MIN_S) & (rr_all < RR_MAX_S)]
        return {
            "family": family, "db": db, "rec_name": rec_name, "ok": True,
            "fs": fs, "duration_s": float(len(ecg) / fs),
            "ecg_idx": ecg_idx, "ecg_name": rec.sig_name[ecg_idx],
            "n_rpeaks": int(len(qrs)), "n_rr_valid": int(len(rr)),
            "rr": rr.tolist(), "elapsed_s": float(time.time() - t0),
        }
    except Exception as e:
        return {
            "family": family, "db": db, "rec_name": rec_name, "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "elapsed_s": float(time.time() - t0),
        }


def analyze_record(payload):
    rr = np.array(payload["rr"], float)
    n_rr = len(rr)
    if n_rr < WINDOW_BEATS + STEP_BEATS * 5:
        return {"ok": False, "reason": f"too few RR ({n_rr})"}
    n_win = (n_rr - WINDOW_BEATS) // STEP_BEATS + 1
    HR = np.empty(n_win)
    HRV = np.empty(n_win)
    for i in range(n_win):
        s = i * STEP_BEATS
        win = rr[s:s + WINDOW_BEATS]
        HR[i] = 60.0 / win.mean()
        HRV[i] = win.std() * 1000.0
    if n_win < SUPERWINDOW + SUPERSTEP * 5:
        return {"ok": False, "reason": f"too few windows ({n_win})"}
    n_sup = (n_win - SUPERWINDOW) // SUPERSTEP + 1
    Ha = np.empty(n_sup); Hb = np.empty(n_sup); MI = np.empty(n_sup)
    for k in range(n_sup):
        s = k * SUPERSTEP
        sub_hr = HR[s:s + SUPERWINDOW]
        sub_hrv = HRV[s:s + SUPERWINDOW]
        Ha[k] = shannon_entropy(sub_hr)
        Hb[k] = shannon_entropy(sub_hrv)
        MI[k] = mutual_info(sub_hr, sub_hrv)

    fits = {}
    X1 = np.column_stack([Ha, Hb])
    m1 = LinearRegression().fit(X1, MI); pred1 = m1.predict(X1)
    fits["poly1"] = {"form": "a*H_a + b*H_b + c",
                     "coef": {"a": float(m1.coef_[0]), "b": float(m1.coef_[1]),
                              "c": float(m1.intercept_)},
                     "r2": float(r2_score(MI, pred1))}
    X2 = np.column_stack([Ha, Hb, Ha**2, Hb**2, Ha * Hb])
    m2 = LinearRegression().fit(X2, MI); pred2 = m2.predict(X2)
    fits["poly2"] = {"form": "a*H_a + b*H_b + d*H_a^2 + e*H_b^2 + f*H_a*H_b + c",
                     "coef": {"a": float(m2.coef_[0]), "b": float(m2.coef_[1]),
                              "d": float(m2.coef_[2]), "e": float(m2.coef_[3]),
                              "f": float(m2.coef_[4]), "c": float(m2.intercept_)},
                     "r2": float(r2_score(MI, pred2))}
    pos = MI > 1e-6
    if pos.sum() >= 5:
        try:
            log_mi = np.log(MI[pos])
            mexp = LinearRegression().fit(Ha[pos].reshape(-1, 1), log_mi)
            B = 1.0 / mexp.coef_[0] if abs(mexp.coef_[0]) > 1e-9 else np.inf
            A = float(np.exp(mexp.intercept_))
            pred_exp = A * np.exp(Ha / B) if np.isfinite(B) else np.full_like(MI, A)
            fits["exp_INFO025"] = {"form": "A * exp(H_a / B)",
                                   "coef": {"A": float(A), "B": float(B)},
                                   "r2": float(r2_score(MI, pred_exp))}
        except Exception as e:
            fits["exp_INFO025"] = {"form": "A*exp(H_a/B)", "error": str(e)}
    diff_sq = (Hb - Ha) ** 2
    m_phys = LinearRegression().fit(diff_sq.reshape(-1, 1), MI)
    pred_phys = m_phys.predict(diff_sq.reshape(-1, 1))
    fits["physics_INFO025"] = {"form": "alpha*(H_b - H_a)^2 + c",
                               "coef": {"alpha": float(m_phys.coef_[0]),
                                        "c": float(m_phys.intercept_)},
                               "r2": float(r2_score(MI, pred_phys))}
    valid = {k: v for k, v in fits.items() if isinstance(v.get("r2"), float)}
    best = max(valid, key=lambda k: valid[k]["r2"]) if valid else None
    return {"ok": True, "n_windows": int(n_win), "n_supertimes": int(n_sup),
            "HR": {"mean": float(HR.mean()), "std": float(HR.std()),
                   "min": float(HR.min()), "max": float(HR.max())},
            "HRV": {"mean": float(HRV.mean()), "std": float(HRV.std()),
                    "min": float(HRV.min()), "max": float(HRV.max())},
            "Ha": {"mean": float(Ha.mean()), "std": float(Ha.std())},
            "Hb": {"mean": float(Hb.mean()), "std": float(Hb.std())},
            "MI": {"mean": float(MI.mean()), "std": float(MI.std())},
            "fits": fits, "best_fit": best,
            "best_r2": float(valid[best]["r2"]) if best else None,
            "Ha_arr": Ha.tolist(), "Hb_arr": Hb.tolist(), "MI_arr": MI.tolist()}


def main():
    START = time.time()
    print(f"Expansion canary: {len(TARGETS)} more families")
    print(f"window={WINDOW_BEATS} step={STEP_BEATS} super={SUPERWINDOW} super_step={SUPERSTEP}\n")
    with Pool(processes=min(6, len(TARGETS))) as pool:
        payloads = pool.starmap(extract_ecg, TARGETS)
    out = []
    for p in payloads:
        if p["ok"]:
            print(f"  {p['family']:25s} {p['db']:8s}/{p['rec_name']:24s}  "
                  f"fs={p['fs']:.0f} dur={p['duration_s']/60:5.1f}min  "
                  f"RR={p['n_rr_valid']:5d}  {p['elapsed_s']:5.1f}s")
            ana = analyze_record(p)
            if not ana["ok"]:
                print(f"      SKIP: {ana['reason']}")
                out.append({"family": p["family"], "ok": False, "reason": ana["reason"]})
                continue
            r = {"family": p["family"], "db": p["db"], "rec_name": p["rec_name"],
                 "duration_s": p["duration_s"], "fs": p["fs"], "ok": True, **ana}
            out.append(r)
            f = ana["fits"]
            print(f"      n_sup={ana['n_supertimes']:3d}  HR={ana['HR']['mean']:5.1f}+/-{ana['HR']['std']:4.1f}  "
                  f"HRV={ana['HRV']['mean']:5.1f}+/-{ana['HRV']['std']:4.1f}ms  "
                  f"MI={ana['MI']['mean']:.3f}+/-{ana['MI']['std']:.3f}  best={ana['best_fit']} R^2={ana['best_r2']:.3f}")
            for k, v in f.items():
                if isinstance(v.get("r2"), float):
                    print(f"      {k:20s} R^2={v['r2']:+6.3f}")
        else:
            print(f"  {p['family']:25s} FAIL: {p['error']}")
            out.append({"family": p["family"], "ok": False, "error": p["error"]})
    print(f"\nDone in {time.time()-START:.1f}s")
    Path("ekg_disease_canary_expand.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
