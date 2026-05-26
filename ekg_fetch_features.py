"""Fetch all 40 Fantasia records, extract RR intervals + windowed HR/HRV.

Patient = ensemble member, window-index = time. Saves per-record RR
arrays + windowed (HR, HRV) features to ekg_features.npz so the
per-domain stack can run without re-fetching.

Parallel via multiprocessing.Pool to keep wall-clock under ~10 min.
"""
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import wfdb
import wfdb.processing as wp

WINDOW_BEATS = 30
STEP_BEATS = 15
RR_MIN_S = 0.3
RR_MAX_S = 2.0


def extract_one(rec_name):
    t0 = time.time()
    try:
        rec = wfdb.rdrecord(rec_name, pn_dir="fantasia")
        fs = rec.fs
        ecg_idx = next(
            (i for i, n in enumerate(rec.sig_name) if "ECG" in n.upper()),
            1 if rec.p_signal.shape[1] > 1 else 0,
        )
        ecg = rec.p_signal[:, ecg_idx]
        qrs = wp.xqrs_detect(ecg, fs=fs, verbose=False)
        rr_all = np.diff(qrs) / fs
        valid = (rr_all > RR_MIN_S) & (rr_all < RR_MAX_S)
        rr = rr_all[valid]

        n_win = (len(rr) - WINDOW_BEATS) // STEP_BEATS + 1
        if n_win < 10:
            raise RuntimeError(f"too few windows ({n_win})")
        hr_w = np.empty(n_win)
        hrv_w = np.empty(n_win)
        for i in range(n_win):
            s = i * STEP_BEATS
            win = rr[s:s + WINDOW_BEATS]
            hr_w[i] = 60.0 / win.mean()
            hrv_w[i] = win.std() * 1000.0

        return {
            "rec_name": rec_name,
            "ok": True,
            "fs": float(fs),
            "rr": rr,
            "hr_w": hr_w,
            "hrv_w": hrv_w,
            "n_rpeaks": int(len(qrs)),
            "n_rr_valid": int(len(rr)),
            "n_windows": int(n_win),
            "elapsed_s": float(time.time() - t0),
        }
    except Exception as e:
        return {
            "rec_name": rec_name,
            "ok": False,
            "error": f"{type(e).__name__}: {e}",
            "elapsed_s": float(time.time() - t0),
        }


def main():
    START = time.time()
    records = wfdb.get_record_list("fantasia")
    print(f"Fantasia: {len(records)} records")
    print(f"Window: {WINDOW_BEATS} beats, step {STEP_BEATS} beats")

    N_WORKERS = 6
    print(f"Pool: {N_WORKERS} workers\n")

    with Pool(N_WORKERS) as pool:
        out = []
        for i, r in enumerate(pool.imap_unordered(extract_one, records)):
            status = "OK" if r["ok"] else f"FAIL {r.get('error', '?')}"
            extra = (
                f"n_rr={r.get('n_rr_valid', '-')} "
                f"n_win={r.get('n_windows', '-')}"
                if r["ok"]
                else ""
            )
            print(
                f"  [{i+1:2d}/{len(records)}] {r['rec_name']:6s} "
                f"{r['elapsed_s']:5.1f}s {status} {extra}"
            )
            out.append(r)

    elapsed = time.time() - START
    print(f"\nFetch done in {elapsed:.1f}s ({elapsed/60:.1f} min)")

    ok = [r for r in out if r["ok"]]
    print(f"OK: {len(ok)}/{len(out)}")

    # Save per-record arrays. Pack RR + windowed features as object arrays
    # so each record keeps its own length.
    rec_names = [r["rec_name"] for r in ok]
    cohort = ["young" if n.startswith("f1y") else "elderly" for n in rec_names]
    rr_list = np.array([r["rr"] for r in ok], dtype=object)
    hr_list = np.array([r["hr_w"] for r in ok], dtype=object)
    hrv_list = np.array([r["hrv_w"] for r in ok], dtype=object)

    np.savez(
        "ekg_features.npz",
        rec_names=np.array(rec_names),
        cohort=np.array(cohort),
        rr=rr_list,
        hr_w=hr_list,
        hrv_w=hrv_list,
        window_beats=WINDOW_BEATS,
        step_beats=STEP_BEATS,
        allow_pickle=True,
    )
    print(f"Saved ekg_features.npz")

    summary = {
        "n_records_total": len(records),
        "n_records_ok": len(ok),
        "elapsed_s": float(elapsed),
        "window_beats": WINDOW_BEATS,
        "step_beats": STEP_BEATS,
        "per_record": [
            {k: v for k, v in r.items() if k not in ("rr", "hr_w", "hrv_w")}
            for r in out
        ],
    }
    Path("ekg_features_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"Saved ekg_features_summary.json")


if __name__ == "__main__":
    main()
