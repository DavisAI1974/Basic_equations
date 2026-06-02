"""
PROBE (S17, "MI next" -- Greg): characterize the inter-detector MI merger signal --
the ONE genuinely physical (scale-invariant) footing left after S17 PROBE 1 showed
the equal-entropy substrate is bookkeeping.

THE question: is the MI-at-merger (INFO-036/045) genuine GRAVITATIONAL information
that carries structure BEYOND loudness, or just "two detectors saw the same loud
thing" (detection bookkeeping)? Four tests on GW150914 H1/L1 (data in hand), all
using the SAME hist-MI estimator (kbk_pipeline.mi_hist_2d) we use everywhere.

NOTE: MI is invariant under negation of one channel (negation is invertible), so the
GW150914 L1 inversion is irrelevant to MI; only the inter-detector LAG matters.

T1. MI(t) through the merger -- reproduce the peak (INFO-036).
T2. LAG SCAN (decisive geometry test): MI(H1_win, L1_win shifted by lag) vs lag in
    the merger window. GENUINE common gravitational signal peaks at the PHYSICAL
    H1-L1 light-travel lag (~7 ms for GW150914); a loudness coincidence is lag-
    independent. Off-merger (noise) window = flat control.
T3. TIME-SLIDE NULL (significance): MI in the merger window at large, non-physical
    lags (0.2-0.8 s >> light travel) decorrelates any real signal -> null
    distribution; report the merger-MI z-score vs that null.
T4. PHASE-SCRAMBLE control (beyond loudness/PSD): phase-randomize L1 in the merger
    window (preserves amplitude spectrum, destroys waveform phase/time structure),
    recompute MI over realizations. If MI collapses toward the null, MI carries
    waveform PHASE/TIME structure, not just amplitude.

Speaking posture (Rule C): no verdict in advance; an "it's just loudness" outcome is
real data, kept, not predetermined bad.

Run:  python probe_mi_merger_axis.py --canary
      python probe_mi_merger_axis.py
"""
import sys
import json
import time
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch
from scipy.interpolate import interp1d

from kbk_pipeline import mi_hist_2d

FS = 4096
CANARY = "--canary" in sys.argv
MERGER_T = 16.42          # s into the 32 s GW150914 segment (our t_c)


def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35.0, hi=350.0, order=4):
    b, a = butter(order, [lo / (FS / 2), hi / (FS / 2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x):
    n = len(x)
    fr, psd = welch(x, fs=FS, nperseg=min(4 * FS, n))
    ip = interp1d(fr, psd, bounds_error=False, fill_value=(psd[0], psd[-1]))
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / FS)
    asd = np.sqrt(np.maximum(ip(f), 1e-50))
    xw = np.fft.irfft(X / asd, n=n)
    return xw / np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def win_slice(x, center_t, win_s):
    c = int(center_t * FS)
    half = int(win_s * FS / 2)
    return x[c - half:c + half]


def mi_at_lag(h, l, center_t, win_s, lag_samples, bins=16):
    """MI between H1 window and L1 window shifted by lag_samples."""
    c = int(center_t * FS)
    half = int(win_s * FS / 2)
    hw = h[c - half:c + half]
    lw = l[c - half - lag_samples:c + half - lag_samples]
    n = min(len(hw), len(lw))
    return mi_hist_2d(hw[:n], lw[:n], bins=bins)


def main():
    t0 = time.time()
    h = preprocess(load_strain("data/ligo/H-H1_GW150914_32s.hdf5"))
    l = preprocess(load_strain("data/ligo/L-L1_GW150914_32s.hdf5"))
    win_s = 0.15
    bins = 16
    n_scramble = 50 if not CANARY else 10

    print("=" * 92)
    print(f"PROBE MI merger axis {'[CANARY]' if CANARY else '[FULL]'}  "
          "(GW150914 H1/L1; is merger MI physics or loudness?)")
    print("=" * 92, flush=True)

    out = {}

    # ---- T1: MI(t) through the merger (L1 pre-aligned at physical lag) ----
    phys_lag = int(round(0.0069 * FS))   # ~7 ms; sign set by L1-leads-H1
    stride = 0.0156
    t_lo, t_hi = MERGER_T - 0.6, MERGER_T + 0.4
    ts, mis = [], []
    t = t_lo
    while t <= t_hi:
        mis.append(mi_at_lag(h, l, t, win_s, phys_lag, bins))
        ts.append(t); t += stride
    ts, mis = np.array(ts), np.array(mis)
    i_pk = int(np.argmax(mis))
    out["T1_MI_time"] = dict(t=ts.tolist(), mi=mis.tolist(),
                             peak_t=float(ts[i_pk]), peak_mi=float(mis[i_pk]),
                             baseline_mi=float(np.median(mis)))
    print(f"\nT1  MI(t): peak MI={mis[i_pk]:.4f} at t={ts[i_pk]:.3f}s "
          f"(merger {MERGER_T:.2f}s); baseline(median)={np.median(mis):.4f}",
          flush=True)

    # ---- T2: lag scan at merger vs off-merger noise window ----
    lags_ms = np.arange(-20, 20.1, 0.5)
    lags = (lags_ms / 1000.0 * FS).round().astype(int)
    mi_merger = np.array([mi_at_lag(h, l, MERGER_T, win_s, lg, bins) for lg in lags])
    mi_noise = np.array([mi_at_lag(h, l, MERGER_T - 4.0, win_s, lg, bins) for lg in lags])
    j = int(np.argmax(mi_merger))
    out["T2_lag_scan"] = dict(lags_ms=lags_ms.tolist(),
                              mi_merger=mi_merger.tolist(),
                              mi_noise=mi_noise.tolist(),
                              peak_lag_ms=float(lags_ms[j]),
                              peak_mi=float(mi_merger[j]))
    print(f"T2  lag scan: merger MI peaks at lag={lags_ms[j]:+.1f} ms "
          f"(physical |lag|~7 ms), MI={mi_merger[j]:.4f}; "
          f"noise-window MI max={mi_noise.max():.4f} (flat control)", flush=True)

    # ---- T3: time-slide null (large non-physical lags) ----
    slide_lags = np.array([int(s * FS) for s in
                           ([0.2, 0.35, 0.5, 0.65, 0.8] if not CANARY
                            else [0.3, 0.5, 0.7])])
    slide_lags = np.concatenate([slide_lags, -slide_lags])
    mi_slides = np.array([mi_at_lag(h, l, MERGER_T, win_s, lg, bins)
                          for lg in slide_lags])
    null_mu, null_sd = float(mi_slides.mean()), float(mi_slides.std() + 1e-12)
    mi_phys = float(mi_merger[j])
    z = (mi_phys - null_mu) / null_sd
    out["T3_time_slide_null"] = dict(slide_mi=mi_slides.tolist(),
                                     null_mean=null_mu, null_std=null_sd,
                                     merger_mi=mi_phys, z_score=z)
    print(f"T3  time-slide null: mean={null_mu:.4f} sd={null_sd:.4f}; "
          f"merger MI={mi_phys:.4f}  => z={z:.1f}", flush=True)

    # ---- T4: phase-scramble control (preserve amplitude spectrum) ----
    c = int(MERGER_T * FS); half = int(win_s * FS / 2)
    hw = h[c - half:c + half]
    lw = l[c - half - phys_lag:c + half - phys_lag]
    n = min(len(hw), len(lw)); hw, lw = hw[:n], lw[:n]
    mi_true = mi_hist_2d(hw, lw, bins=bins)
    rng = np.random.default_rng(11)
    L = np.fft.rfft(lw)
    amp = np.abs(L)
    mi_scr = []
    for _ in range(n_scramble):
        ph = np.exp(1j * rng.uniform(0, 2 * np.pi, len(L)))
        ph[0] = 1.0
        lscr = np.fft.irfft(amp * ph, n=n)
        mi_scr.append(mi_hist_2d(hw, lscr, bins=bins))
    mi_scr = np.array(mi_scr)
    out["T4_phase_scramble"] = dict(mi_true=float(mi_true),
                                    mi_scramble_mean=float(mi_scr.mean()),
                                    mi_scramble_std=float(mi_scr.std()),
                                    mi_scramble_max=float(mi_scr.max()),
                                    n_scramble=int(n_scramble))
    print(f"T4  phase-scramble: true MI={mi_true:.4f}; "
          f"phase-scrambled MI={mi_scr.mean():.4f}+/-{mi_scr.std():.4f} "
          f"(max {mi_scr.max():.4f})  => drop factor "
          f"{mi_true / (mi_scr.mean() + 1e-12):.1f}x", flush=True)

    result = dict(meta=dict(canary=CANARY, win_s=win_s, bins=bins,
                            merger_t=MERGER_T, elapsed_s=round(time.time() - t0, 1)),
                  tests=out)
    fn = "probe_mi_merger_axis_canary.json" if CANARY \
        else "probe_mi_merger_axis_results.json"
    json.dump(result, open(fn, "w"), indent=2)
    print(f"\nWrote {fn}  ({result['meta']['elapsed_s']}s)", flush=True)


if __name__ == "__main__":
    main()
