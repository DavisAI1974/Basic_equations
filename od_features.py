"""Feature extractors: raw signal -> per-domain stack feature dict.

Bridges the gap between "I have an ECG record" and "I have a feature
dict ready for the OD learning store". Reusable across the cerebro,
cardiac, and any future probe.

Two extractors so far:
  - extract_cardiac_hr_hrv(ecg, fs)          coupled pair: (HR, HRV)
  - extract_baroreflex_hr_abp(ecg, abp, fs)  coupled pair: (HR, ABP)

Both produce the same feature schema (so the same trained classifier
can predict on either):
  HR_mean, HR_std, HRV_mean (or ABP_mean), MI_mean, ...,
  poly1_a, poly1_b, ..., poly2_*, r2_poly1, r2_poly2,
  exp_A, exp_B, r2_exp, phys_alpha, r2_phys,
  Ha_p_low, Ha_p_mid, Ha_p_hi, Ha_ac1, Ha_ac3, Ha_diff_std,
  Hb_*, MI_*  (temporal features on each sequence)

Sequences (H_a, H_b, MI) are NOT returned (they would bloat the
feature dict); only their scalar summaries.
"""
from typing import Optional

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
    x = np.asarray(x, float); y = np.asarray(y, float)
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


def temporal_features(seq, prefix):
    """FFT band powers + autocorrelation + jaggedness of a 1D sequence."""
    s = np.asarray(seq, float)
    s = s - s.mean()
    n = len(s)
    out = {}
    spec = np.abs(np.fft.rfft(s)) ** 2
    freqs = np.fft.rfftfreq(n)
    total_p = spec.sum() + 1e-12
    out[f"{prefix}_p_low"] = float(spec[(freqs < 0.05)].sum() / total_p)
    out[f"{prefix}_p_mid"] = float(spec[(freqs >= 0.05) & (freqs < 0.20)].sum() / total_p)
    out[f"{prefix}_p_hi"] = float(spec[(freqs >= 0.20)].sum() / total_p)
    if n > 4 and s.std() > 1e-9:
        out[f"{prefix}_ac1"] = float(np.corrcoef(s[:-1], s[1:])[0, 1])
        out[f"{prefix}_ac3"] = float(np.corrcoef(s[:-3], s[3:])[0, 1]) if n > 6 else 0.0
    else:
        out[f"{prefix}_ac1"] = 0.0
        out[f"{prefix}_ac3"] = 0.0
    out[f"{prefix}_diff_std"] = float(np.diff(s).std()) if n > 1 else 0.0
    return out


def _fit_baselines(Ha, Hb, MI):
    """Fit poly1, poly2, exp_INFO025, physics_INFO025; return flat dict."""
    fits = {}
    X1 = np.column_stack([Ha, Hb])
    m1 = LinearRegression().fit(X1, MI)
    p1 = m1.predict(X1)
    fits.update({
        "poly1_a": float(m1.coef_[0]), "poly1_b": float(m1.coef_[1]),
        "poly1_c": float(m1.intercept_),
        "r2_poly1": float(r2_score(MI, p1)),
    })
    X2 = np.column_stack([Ha, Hb, Ha**2, Hb**2, Ha * Hb])
    m2 = LinearRegression().fit(X2, MI)
    p2 = m2.predict(X2)
    fits.update({
        "poly2_a": float(m2.coef_[0]), "poly2_b": float(m2.coef_[1]),
        "poly2_d": float(m2.coef_[2]), "poly2_e": float(m2.coef_[3]),
        "poly2_f": float(m2.coef_[4]), "poly2_c": float(m2.intercept_),
        "r2_poly2": float(r2_score(MI, p2)),
    })
    pos = MI > 1e-6
    if pos.sum() >= 5:
        log_mi = np.log(MI[pos])
        mexp = LinearRegression().fit(Ha[pos].reshape(-1, 1), log_mi)
        B = 1.0 / mexp.coef_[0] if abs(mexp.coef_[0]) > 1e-9 else float("inf")
        A = float(np.exp(mexp.intercept_))
        pred = A * np.exp(Ha / B) if np.isfinite(B) else np.full_like(MI, A)
        fits.update({"exp_A": float(A), "exp_B": float(B),
                     "r2_exp": float(r2_score(MI, pred))})
    else:
        fits.update({"exp_A": 0.0, "exp_B": 0.0, "r2_exp": 0.0})
    diff_sq = (Hb - Ha) ** 2
    m_phys = LinearRegression().fit(diff_sq.reshape(-1, 1), MI)
    fits.update({
        "phys_alpha": float(m_phys.coef_[0]),
        "phys_c": float(m_phys.intercept_),
        "r2_phys": float(r2_score(MI, m_phys.predict(diff_sq.reshape(-1, 1)))),
    })
    return fits


def _ab_a_hat(fits):
    """Convenience: unit-vector form of (poly1_a, poly1_b)."""
    v = np.array([fits["poly1_a"], fits["poly1_b"]])
    n = np.linalg.norm(v)
    if n < 1e-9: return 0.0, 0.0
    return float(v[0] / n), float(v[1] / n)


def extract_cardiac_hr_hrv(ecg, fs, *, cap_seconds: Optional[float] = 7200,
                           window_beats=WINDOW_BEATS, step_beats=STEP_BEATS,
                           superwindow=SUPERWINDOW, superstep=SUPERSTEP):
    """Build per-record feature dict from a single ECG channel.

    Pipeline: ECG -> R-peaks -> RR -> (HR, HRV) per beat-window ->
    (H_a, H_b, MI) per supertime via sliding window -> baseline fits +
    temporal features.

    Returns dict on success, raises on insufficient data.
    """
    ecg = np.asarray(ecg, float)
    if np.isnan(ecg).any():
        ecg = np.where(np.isnan(ecg), 0.0, ecg)
    if cap_seconds is not None and len(ecg) > int(cap_seconds * fs):
        ecg = ecg[:int(cap_seconds * fs)]
    qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
    rr = np.diff(qrs) / fs
    rr = rr[(rr > RR_MIN_S) & (rr < RR_MAX_S)]
    if len(rr) < window_beats + step_beats * 3:
        raise RuntimeError(f"too few RR intervals ({len(rr)})")
    n_win = (len(rr) - window_beats) // step_beats + 1
    HR = np.empty(n_win); HRV = np.empty(n_win)
    for i in range(n_win):
        s = i * step_beats
        win = rr[s:s + window_beats]
        HR[i] = 60.0 / win.mean()
        HRV[i] = win.std() * 1000.0
    if n_win < superwindow + superstep * 2:
        raise RuntimeError(f"too few windows ({n_win})")
    n_sup = (n_win - superwindow) // superstep + 1
    Ha = np.empty(n_sup); Hb = np.empty(n_sup); MI = np.empty(n_sup)
    for k in range(n_sup):
        s = k * superstep
        Ha[k] = shannon_entropy(HR[s:s + superwindow])
        Hb[k] = shannon_entropy(HRV[s:s + superwindow])
        MI[k] = mutual_info(HR[s:s + superwindow], HRV[s:s + superwindow])

    feats = {
        "n_rr": int(len(rr)), "n_windows": int(n_win), "n_sup": int(n_sup),
        "HR_mean": float(HR.mean()), "HR_std": float(HR.std()),
        "ABP_mean": float(HRV.mean()), "ABP_std": float(HRV.std()),  # alias for schema compatibility
        "MI_mean": float(MI.mean()), "MI_std": float(MI.std()),
        "Ha_mean": float(Ha.mean()), "Ha_std": float(Ha.std()),
        "Hb_mean": float(Hb.mean()), "Hb_std": float(Hb.std()),
    }
    feats.update(_fit_baselines(Ha, Hb, MI))
    a_hat, b_hat = _ab_a_hat(feats)
    feats["a_hat"] = a_hat; feats["b_hat"] = b_hat
    feats.update(temporal_features(Ha, "Ha"))
    feats.update(temporal_features(Hb, "Hb"))
    feats.update(temporal_features(MI, "MI"))
    return feats


def extract_baroreflex_hr_abp(ecg, abp, fs, *, cap_seconds=None,
                              window_beats=20, step_beats=8,
                              superwindow=20, superstep=5):
    """Build per-record feature dict from coupled (ECG, ABP) signals.

    Matches the cerebro probe protocol: per-beat ABP sampled at the
    R-peak instant; HR and ABP windowed per N beats; sliding window of
    windows produces (H_a, H_b, MI) per supertime.
    """
    ecg = np.asarray(ecg, float); abp = np.asarray(abp, float)
    ecg = np.where(np.isnan(ecg), 0.0, ecg)
    abp = np.where(np.isnan(abp), 0.0, abp)
    if cap_seconds is not None:
        cap = int(cap_seconds * fs)
        if len(ecg) > cap: ecg = ecg[:cap]
        if len(abp) > cap: abp = abp[:cap]
    qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
    if len(qrs) < window_beats + step_beats * 3:
        raise RuntimeError(f"too few R-peaks ({len(qrs)})")
    rr = np.diff(qrs) / fs
    valid = (rr > RR_MIN_S) & (rr < RR_MAX_S)
    rr = rr[valid]
    abp_per_beat = abp[qrs[:-1]]
    n = min(len(rr), len(abp_per_beat))
    rr = rr[:n]; abp_per_beat = abp_per_beat[:n]
    n_win = (n - window_beats) // step_beats + 1
    if n_win < superwindow + superstep * 2:
        raise RuntimeError(f"too few windows ({n_win})")
    HR = np.empty(n_win); ABP = np.empty(n_win)
    for i in range(n_win):
        s = i * step_beats
        HR[i] = 60.0 / rr[s:s + window_beats].mean()
        ABP[i] = abp_per_beat[s:s + window_beats].mean()
    n_sup = (n_win - superwindow) // superstep + 1
    Ha = np.empty(n_sup); Hb = np.empty(n_sup); MI = np.empty(n_sup)
    for k in range(n_sup):
        s = k * superstep
        Ha[k] = shannon_entropy(HR[s:s + superwindow])
        Hb[k] = shannon_entropy(ABP[s:s + superwindow])
        MI[k] = mutual_info(HR[s:s + superwindow], ABP[s:s + superwindow])
    feats = {
        "n_rr": int(n), "n_windows": int(n_win), "n_sup": int(n_sup),
        "HR_mean": float(HR.mean()), "HR_std": float(HR.std()),
        "ABP_mean": float(ABP.mean()), "ABP_std": float(ABP.std()),
        "MI_mean": float(MI.mean()), "MI_std": float(MI.std()),
        "Ha_mean": float(Ha.mean()), "Ha_std": float(Ha.std()),
        "Hb_mean": float(Hb.mean()), "Hb_std": float(Hb.std()),
    }
    feats.update(_fit_baselines(Ha, Hb, MI))
    a_hat, b_hat = _ab_a_hat(feats)
    feats["a_hat"] = a_hat; feats["b_hat"] = b_hat
    feats.update(temporal_features(Ha, "Ha"))
    feats.update(temporal_features(Hb, "Hb"))
    feats.update(temporal_features(MI, "MI"))
    return feats


def extract_from_wfdb_record(rec_name, *, pn_dir=None, local_dir=None,
                             mode="cardiac"):
    """Convenience: load a WFDB record (local or PhysioNet) and run the
    appropriate extractor.

    mode = 'cardiac'  -> HR/HRV pair on first ECG-named channel
    mode = 'baroreflex' -> HR/ABP pair (ecg + abp channels required)
    """
    if local_dir is not None:
        import os
        rec = wfdb.rdrecord(os.path.join(local_dir, rec_name))
    else:
        rec = wfdb.rdrecord(rec_name, pn_dir=pn_dir)
    fs = float(rec.fs)
    ecg_idx = None
    for i, n in enumerate(rec.sig_name):
        nm = n.upper()
        if "ECG" in nm or "MLII" in nm or nm.startswith("II"):
            ecg_idx = i; break
    if ecg_idx is None: ecg_idx = 0
    ecg = rec.p_signal[:, ecg_idx]
    if mode == "cardiac":
        feats = extract_cardiac_hr_hrv(ecg, fs)
    elif mode == "baroreflex":
        abp_idx = None
        for i, n in enumerate(rec.sig_name):
            if "ABP" in n.upper() or "BP" in n.upper():
                abp_idx = i; break
        if abp_idx is None:
            raise RuntimeError(f"no ABP/BP channel found in {rec.sig_name}")
        abp = rec.p_signal[:, abp_idx]
        feats = extract_baroreflex_hr_abp(ecg, abp, fs)
    else:
        raise ValueError(f"unknown mode {mode!r}")
    feats["_rec_name"] = rec_name
    feats["_fs"] = fs
    feats["_duration_s"] = float(len(ecg) / fs)
    return feats
