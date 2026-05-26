"""Fetch 25 AF records (afdb) + 18 NSR records (nsrdb) for proper
AF-vs-NSR classification test. Compute per-subject features matching
the cerebro probe (HR;HRV pair, sliding window of windows, fit poly1/
poly2/exp/physics + temporal features of the (H_a, H_b, MI) sequences).
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

WINDOW_BEATS = 30
STEP_BEATS = 15
SUPERWINDOW = 40
SUPERSTEP = 10
RR_MIN_S = 0.3
RR_MAX_S = 2.0
N_BINS_H = 16
CAP_HOURS = 2.0


def shannon_entropy(x, bins=N_BINS_H):
    x = np.asarray(x, float)
    if x.size < 4 or np.allclose(x.std(), 0.0):
        return 0.0
    lo, hi = np.percentile(x, [1.0, 99.0])
    if hi - lo < 1e-12: return 0.0
    hist, _ = np.histogram(x, bins=bins, range=(lo, hi))
    p = hist.astype(float) / hist.sum()
    p = p[p > 0]
    return float(-np.sum(p * np.log(p)))


def mutual_info(x, y, bins=N_BINS_H):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if x.size < 4 or x.std() < 1e-12 or y.std() < 1e-12: return 0.0
    lo_x, hi_x = np.percentile(x, [1.0, 99.0])
    lo_y, hi_y = np.percentile(y, [1.0, 99.0])
    if hi_x - lo_x < 1e-12 or hi_y - lo_y < 1e-12: return 0.0
    Hxy, _, _ = np.histogram2d(x, y, bins=bins,
                               range=[[lo_x, hi_x], [lo_y, hi_y]])
    Pxy = Hxy / Hxy.sum()
    Px = Pxy.sum(axis=1, keepdims=True)
    Py = Pxy.sum(axis=0, keepdims=True)
    mask = (Pxy > 0) & (Px > 0) & (Py > 0)
    return float(np.sum(Pxy[mask] * np.log(Pxy[mask] / (Px * Py)[mask])))


def temporal_features(seq, prefix):
    """FFT bands, autocorrelation, sample entropy of a 1D sequence."""
    s = np.asarray(seq, float)
    s = s - s.mean()
    n = len(s)
    out = {}
    # FFT power in 3 bands
    spec = np.abs(np.fft.rfft(s)) ** 2
    freqs = np.fft.rfftfreq(n)
    total_p = spec.sum() + 1e-12
    out[f"{prefix}_p_low"] = float(spec[(freqs < 0.05)].sum() / total_p)
    out[f"{prefix}_p_mid"] = float(spec[(freqs >= 0.05) & (freqs < 0.20)].sum() / total_p)
    out[f"{prefix}_p_hi"] = float(spec[(freqs >= 0.20)].sum() / total_p)
    # Lag-1 and lag-3 autocorrelation
    if n > 4 and s.std() > 1e-9:
        out[f"{prefix}_ac1"] = float(np.corrcoef(s[:-1], s[1:])[0, 1])
        out[f"{prefix}_ac3"] = float(np.corrcoef(s[:-3], s[3:])[0, 1]) if n > 6 else 0.0
    else:
        out[f"{prefix}_ac1"] = 0.0; out[f"{prefix}_ac3"] = 0.0
    # Variance of finite differences (proxy for jaggedness)
    if n > 1:
        out[f"{prefix}_diff_std"] = float(np.diff(s).std())
    else:
        out[f"{prefix}_diff_std"] = 0.0
    return out


def process_record(args):
    label, db, rec_name = args
    t0 = time.time()
    try:
        rec = wfdb.rdrecord(rec_name, pn_dir=db)
        fs = float(rec.fs)
        ecg_idx = None
        for i, n in enumerate(rec.sig_name):
            nm = n.upper()
            if "ECG" in nm or "MLII" in nm or nm.startswith("II"):
                ecg_idx = i; break
        if ecg_idx is None: ecg_idx = 0
        ecg = rec.p_signal[:, ecg_idx]
        ecg = np.where(np.isnan(ecg), 0.0, ecg)
        cap = int(CAP_HOURS * 3600 * fs)
        if len(ecg) > cap: ecg = ecg[:cap]
        qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
        rr = np.diff(qrs) / fs
        rr = rr[(rr > RR_MIN_S) & (rr < RR_MAX_S)]
        if len(rr) < WINDOW_BEATS + STEP_BEATS * 3:
            raise RuntimeError(f"too few rr ({len(rr)})")
        n_win = (len(rr) - WINDOW_BEATS) // STEP_BEATS + 1
        HR = np.empty(n_win); HRV = np.empty(n_win)
        for i in range(n_win):
            s = i * STEP_BEATS
            win = rr[s:s + WINDOW_BEATS]
            HR[i] = 60.0 / win.mean()
            HRV[i] = win.std() * 1000.0
        if n_win < SUPERWINDOW + SUPERSTEP * 2:
            raise RuntimeError(f"too few windows ({n_win})")
        n_sup = (n_win - SUPERWINDOW) // SUPERSTEP + 1
        Ha = np.empty(n_sup); Hb = np.empty(n_sup); MI = np.empty(n_sup)
        for k in range(n_sup):
            s = k * SUPERSTEP
            Ha[k] = shannon_entropy(HR[s:s + SUPERWINDOW])
            Hb[k] = shannon_entropy(HRV[s:s + SUPERWINDOW])
            MI[k] = mutual_info(HR[s:s + SUPERWINDOW], HRV[s:s + SUPERWINDOW])
        # Fits
        X1 = np.column_stack([Ha, Hb])
        m1 = LinearRegression().fit(X1, MI); r2_p1 = r2_score(MI, m1.predict(X1))
        X2 = np.column_stack([Ha, Hb, Ha**2, Hb**2, Ha * Hb])
        m2 = LinearRegression().fit(X2, MI); r2_p2 = r2_score(MI, m2.predict(X2))
        pos = MI > 1e-6
        if pos.sum() >= 5:
            log_mi = np.log(MI[pos])
            mexp = LinearRegression().fit(Ha[pos].reshape(-1, 1), log_mi)
            B = 1.0 / mexp.coef_[0] if abs(mexp.coef_[0]) > 1e-9 else np.inf
            A = float(np.exp(mexp.intercept_))
            pred = A * np.exp(Ha / B) if np.isfinite(B) else np.full_like(MI, A)
            r2_exp = float(r2_score(MI, pred))
        else:
            A, B, r2_exp = 0.0, 0.0, 0.0
        diff_sq = (Hb - Ha) ** 2
        m_phys = LinearRegression().fit(diff_sq.reshape(-1, 1), MI)
        r2_phys = r2_score(MI, m_phys.predict(diff_sq.reshape(-1, 1)))

        feats = {
            "label": label, "rec": rec_name, "ok": True,
            "n_rr": int(len(rr)), "n_windows": int(n_win),
            "n_sup": int(n_sup), "elapsed_s": float(time.time() - t0),
            "HR_mean": float(HR.mean()), "HR_std": float(HR.std()),
            "HRV_mean": float(HRV.mean()), "HRV_std": float(HRV.std()),
            "Ha_mean": float(Ha.mean()), "Ha_std": float(Ha.std()),
            "Hb_mean": float(Hb.mean()), "Hb_std": float(Hb.std()),
            "MI_mean": float(MI.mean()), "MI_std": float(MI.std()),
            "poly1_a": float(m1.coef_[0]), "poly1_b": float(m1.coef_[1]),
            "poly1_c": float(m1.intercept_), "r2_poly1": float(r2_p1),
            "poly2_a": float(m2.coef_[0]), "poly2_b": float(m2.coef_[1]),
            "poly2_d": float(m2.coef_[2]), "poly2_e": float(m2.coef_[3]),
            "poly2_f": float(m2.coef_[4]), "poly2_c": float(m2.intercept_),
            "r2_poly2": float(r2_p2),
            "exp_A": float(A), "exp_B": float(B), "r2_exp": float(r2_exp),
            "phys_alpha": float(m_phys.coef_[0]),
            "phys_c": float(m_phys.intercept_), "r2_phys": float(r2_phys),
        }
        # Temporal features on each sequence
        feats.update(temporal_features(Ha, "Ha"))
        feats.update(temporal_features(Hb, "Hb"))
        feats.update(temporal_features(MI, "MI"))
        return feats
    except Exception as e:
        return {"label": label, "rec": rec_name, "ok": False,
                "error": f"{type(e).__name__}: {str(e)[:80]}",
                "elapsed_s": float(time.time() - t0)}


def main():
    START = time.time()
    # Build target list
    afdb_recs = wfdb.get_record_list("afdb")
    nsrdb_recs = wfdb.get_record_list("nsrdb")
    print(f"afdb: {len(afdb_recs)} records  nsrdb: {len(nsrdb_recs)} records")

    targets = []
    for r in afdb_recs[:25]:
        targets.append(("AF", "afdb", r))
    for r in nsrdb_recs[:18]:
        targets.append(("NSR", "nsrdb", r))
    print(f"Total targets: {len(targets)}\n")

    with Pool(processes=8) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(process_record, targets)):
            status = "OK" if r["ok"] else f"FAIL {r.get('error', '?')}"
            extra = ""
            if r["ok"]:
                extra = (f"HR={r['HR_mean']:5.1f}+/-{r['HR_std']:4.1f}  "
                         f"HRV={r['HRV_mean']:5.1f}  MI={r['MI_mean']:.3f}  "
                         f"poly2={r['r2_poly2']:+.3f}")
            print(f"  [{i+1:2d}/{len(targets)}] {r['label']:3s} {r['rec']:10s} "
                  f"{r['elapsed_s']:5.1f}s  {status}  {extra}")
            out.append(r)

    elapsed = time.time() - START
    ok = [r for r in out if r["ok"]]
    print(f"\nDone in {elapsed:.1f}s ({elapsed/60:.1f} min)  OK: {len(ok)}/{len(out)}")
    Path("cardiac_af_cohort.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
