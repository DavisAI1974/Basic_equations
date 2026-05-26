"""ECG canary: pull 5 Fantasia records, detect R-peaks, build HR/HRV windows.

Sanity check before full per-domain stack. Patient = ensemble member,
window-index = time. At each window t across N patients we will later
compute H_a(t)=entropy(HR distribution), H_b(t)=entropy(HRV), MI(t).

Canary only verifies: data fetches, R-peak detection works, RR intervals
are physiologic, HR/HRV windows have variation across patients.
"""
import json
import time
from pathlib import Path

import numpy as np
import wfdb
import wfdb.processing as wp

START = time.time()
print("=" * 60)
print("ECG canary: 5 Fantasia records, R-peak detection, HR/HRV windows")
print("=" * 60)

RECORDS = ["f1y01", "f1y02", "f1y03", "f1o01", "f1o02"]
WINDOW_BEATS = 30
STEP_BEATS = 15

results = {}

for rec_name in RECORDS:
    t0 = time.time()
    print(f"\n[{rec_name}] fetching...")
    rec = wfdb.rdrecord(rec_name, pn_dir="fantasia")
    fs = rec.fs
    print(f"  signals={rec.sig_name} fs={fs} duration={rec.sig_len/fs:.0f}s")

    # ECG is typically channel 1 in Fantasia (channel 0 is respiration).
    # Detect on whichever signal name matches.
    ecg_idx = None
    for i, name in enumerate(rec.sig_name):
        if "ECG" in name.upper():
            ecg_idx = i
            break
    if ecg_idx is None:
        ecg_idx = 1 if rec.p_signal.shape[1] > 1 else 0
    ecg = rec.p_signal[:, ecg_idx]
    print(f"  ECG channel idx={ecg_idx} name={rec.sig_name[ecg_idx]}")

    qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
    rr = np.diff(qrs) / fs
    rr_valid = rr[(rr > 0.3) & (rr < 2.0)]
    print(f"  R-peaks={len(qrs)} valid_RR={len(rr_valid)} "
          f"RR_mean={rr_valid.mean():.3f}s RR_std={rr_valid.std():.3f}s")

    n_win = (len(rr_valid) - WINDOW_BEATS) // STEP_BEATS + 1
    hr_w = np.zeros(n_win)
    hrv_w = np.zeros(n_win)
    for i in range(n_win):
        s = i * STEP_BEATS
        e = s + WINDOW_BEATS
        win = rr_valid[s:e]
        hr_w[i] = 60.0 / win.mean()
        hrv_w[i] = win.std() * 1000.0

    print(f"  windows={n_win} HR=[{hr_w.min():.1f}, {hr_w.max():.1f}] bpm "
          f"HRV=[{hrv_w.min():.1f}, {hrv_w.max():.1f}] ms (SDNN-30beat)")

    results[rec_name] = {
        "fs": float(fs),
        "duration_s": float(rec.sig_len / fs),
        "ecg_idx": ecg_idx,
        "ecg_name": rec.sig_name[ecg_idx],
        "n_rpeaks": int(len(qrs)),
        "n_rr_valid": int(len(rr_valid)),
        "rr_mean_s": float(rr_valid.mean()),
        "rr_std_s": float(rr_valid.std()),
        "n_windows": int(n_win),
        "hr_min_bpm": float(hr_w.min()),
        "hr_max_bpm": float(hr_w.max()),
        "hrv_min_ms": float(hrv_w.min()),
        "hrv_max_ms": float(hrv_w.max()),
        "elapsed_s": float(time.time() - t0),
    }

elapsed = time.time() - START
print(f"\nCanary done in {elapsed:.1f}s ({elapsed/60:.1f} min)")

out = Path("ekg_canary.json")
out.write_text(json.dumps(results, indent=2))
print(f"Wrote {out}")
