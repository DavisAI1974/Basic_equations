"""Cerebrovascular probe on cerebral-vasoreg-diabetes cohort.

Different sample type than ECG: coupled pair = (HR from ECG, ABP).
Tests substrate-vs-expression frame at the cardiovascular-disease
stratum, on a different signal modality, with 3 stratifications:

  Strat 1: CV pathology (control / DM / DM+HTN / DMOH orthostatic)
  Strat 2: Kidney function (normal / mild / significant impairment via Cre)
  Strat 3: DM complication subtype (control / uncompl / retin / neuro / nephr)

Per-subject pipeline mirrors ekg_disease_canary.py:
  ECG -> R-peaks -> RR -> (HR_window, ABP_window) per beat-window
  -> sliding window of windows -> (H_a, H_b, MI) per supertime
  -> fit poly1/poly2/exp/physics baselines.
"""
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import wfdb
import wfdb.processing as wp
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

DATA_DIR = "Data/Labview/Converted/Head-up-tilt_Day1"
PN_DIR = f"cerebral-vasoreg-diabetes/1.0.0/{DATA_DIR}"
WINDOW_BEATS = 20
STEP_BEATS = 8
SUPERWINDOW = 20
SUPERSTEP = 5
RR_MIN_S = 0.3
RR_MAX_S = 2.0
N_BINS_H = 12


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


def process_subject(subj_id):
    t0 = time.time()
    rec_name = f"{subj_id}DA"
    try:
        rec = wfdb.rdrecord(rec_name, pn_dir=PN_DIR)
        fs = float(rec.fs)
        ecg_idx = next(i for i, n in enumerate(rec.sig_name)
                       if n.lower() == "ecg")
        abp_idx = next(i for i, n in enumerate(rec.sig_name)
                       if n.lower() == "abp")
        ecg = rec.p_signal[:, ecg_idx]
        abp = rec.p_signal[:, abp_idx]
        nans = np.isnan(ecg) | np.isnan(abp)
        if nans.any():
            ecg = np.where(np.isnan(ecg), 0.0, ecg)
            abp = np.where(np.isnan(abp), 0.0, abp)

        qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
        if len(qrs) < WINDOW_BEATS + STEP_BEATS * 3:
            return {"subj": subj_id, "ok": False,
                    "reason": f"too few R-peaks ({len(qrs)})",
                    "elapsed_s": time.time() - t0}
        rr = np.diff(qrs) / fs
        rr = rr[(rr > RR_MIN_S) & (rr < RR_MAX_S)]
        # Per-beat ABP at the R-peak instant
        abp_per_beat = abp[qrs[:-1]]
        # Trim alignment to valid rr
        n = min(len(rr), len(abp_per_beat))
        rr = rr[:n]; abp_per_beat = abp_per_beat[:n]

        n_win = (n - WINDOW_BEATS) // STEP_BEATS + 1
        if n_win < SUPERWINDOW + SUPERSTEP * 2:
            return {"subj": subj_id, "ok": False,
                    "reason": f"too few windows ({n_win})",
                    "elapsed_s": time.time() - t0}
        HR = np.empty(n_win)
        ABP = np.empty(n_win)
        for i in range(n_win):
            s = i * STEP_BEATS
            HR[i] = 60.0 / rr[s:s + WINDOW_BEATS].mean()
            ABP[i] = abp_per_beat[s:s + WINDOW_BEATS].mean()

        n_sup = (n_win - SUPERWINDOW) // SUPERSTEP + 1
        Ha = np.empty(n_sup); Hb = np.empty(n_sup); MI = np.empty(n_sup)
        for k in range(n_sup):
            s = k * SUPERSTEP
            sub_hr = HR[s:s + SUPERWINDOW]
            sub_abp = ABP[s:s + SUPERWINDOW]
            Ha[k] = shannon_entropy(sub_hr)
            Hb[k] = shannon_entropy(sub_abp)
            MI[k] = mutual_info(sub_hr, sub_abp)

        # Fits
        fits = {}
        X1 = np.column_stack([Ha, Hb])
        m1 = LinearRegression().fit(X1, MI); pred = m1.predict(X1)
        fits["poly1"] = {"a": float(m1.coef_[0]), "b": float(m1.coef_[1]),
                         "c": float(m1.intercept_),
                         "r2": float(r2_score(MI, pred))}
        X2 = np.column_stack([Ha, Hb, Ha**2, Hb**2, Ha * Hb])
        m2 = LinearRegression().fit(X2, MI); pred = m2.predict(X2)
        fits["poly2"] = {"a": float(m2.coef_[0]), "b": float(m2.coef_[1]),
                         "d": float(m2.coef_[2]), "e": float(m2.coef_[3]),
                         "f": float(m2.coef_[4]), "c": float(m2.intercept_),
                         "r2": float(r2_score(MI, pred))}
        pos = MI > 1e-6
        if pos.sum() >= 5:
            log_mi = np.log(MI[pos])
            mexp = LinearRegression().fit(Ha[pos].reshape(-1, 1), log_mi)
            B = 1.0 / mexp.coef_[0] if abs(mexp.coef_[0]) > 1e-9 else np.inf
            A = float(np.exp(mexp.intercept_))
            pred = A * np.exp(Ha / B) if np.isfinite(B) else np.full_like(MI, A)
            fits["exp_INFO025"] = {"A": A, "B": float(B),
                                   "r2": float(r2_score(MI, pred))}
        diff_sq = (Hb - Ha) ** 2
        m_phys = LinearRegression().fit(diff_sq.reshape(-1, 1), MI)
        pred = m_phys.predict(diff_sq.reshape(-1, 1))
        fits["physics_INFO025"] = {"alpha": float(m_phys.coef_[0]),
                                   "c": float(m_phys.intercept_),
                                   "r2": float(r2_score(MI, pred))}
        return {"subj": subj_id, "ok": True,
                "fs": fs, "duration_s": float(len(ecg) / fs),
                "n_rpeaks": len(qrs), "n_windows": int(n_win),
                "n_supertimes": int(n_sup),
                "HR_mean": float(HR.mean()), "HR_std": float(HR.std()),
                "ABP_mean": float(ABP.mean()), "ABP_std": float(ABP.std()),
                "MI_mean": float(MI.mean()), "MI_std": float(MI.std()),
                "Ha_mean": float(Ha.mean()), "Hb_mean": float(Hb.mean()),
                "fits": fits,
                "Ha": Ha.tolist(), "Hb": Hb.tolist(), "MI": MI.tolist(),
                "elapsed_s": float(time.time() - t0)}
    except Exception as e:
        return {"subj": subj_id, "ok": False,
                "error": f"{type(e).__name__}: {str(e)[:120]}",
                "elapsed_s": float(time.time() - t0)}


def main():
    START = time.time()
    df = pd.read_csv("data/cvd/cohort_strat.csv")
    subjects = df[df["has_file"]]["SubjectID_lc"].tolist()
    print(f"Processing {len(subjects)} subjects (parallel)")
    print(f"Window={WINDOW_BEATS} step={STEP_BEATS}  super={SUPERWINDOW} super_step={SUPERSTEP}  bins={N_BINS_H}")

    with Pool(processes=8) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(process_subject, subjects)):
            status = "OK" if r["ok"] else f"FAIL {r.get('reason') or r.get('error', '?')}"
            extra = ""
            if r["ok"]:
                f = r["fits"]
                extra = (f"n_sup={r['n_supertimes']:3d}  "
                         f"HR={r['HR_mean']:5.1f}+/-{r['HR_std']:4.1f}  "
                         f"ABP={r['ABP_mean']:5.1f}  MI={r['MI_mean']:.3f}  "
                         f"poly2={f['poly2']['r2']:+.3f} "
                         f"exp={f.get('exp_INFO025', {}).get('r2', 0):+.3f}")
            print(f"  [{i+1:2d}/{len(subjects)}] {r['subj']:6s} "
                  f"{r['elapsed_s']:4.1f}s  {status}  {extra}")
            out.append(r)

    elapsed = time.time() - START
    ok = [r for r in out if r["ok"]]
    print(f"\nDone in {elapsed:.1f}s ({elapsed/60:.1f} min)  OK: {len(ok)}/{len(out)}")

    Path("cerebro_disease_probe.json").write_text(json.dumps(out, indent=2))
    print("Saved cerebro_disease_probe.json")


if __name__ == "__main__":
    main()
