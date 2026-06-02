"""
S10 item (5), set 1 of 2 (REAL DATA): two-detector LIGO strain through
the SAME windowed-H / MI / operator-extraction stack we built on the
toy systems. First time the method touches real physics data.

Why this is the four-force "shared substrate" test on real gravity:
our four-force leg (i) said all caricatures share the channel-
correlation identity -(H_a - H_b)^2 ~ 0 (the (-1,-1,+2)/sqrt(6)
attractor) WHEN the two channels are correlated. Real H1/L1 strain is
mostly INDEPENDENT detector noise (two sites), EXCEPT for the ~0.2 s
around a real gravitational-wave merger, when both detectors see the
SAME astrophysical signal. So the falsifiable real-data prediction is:
the operator-space null should sit on the channel-correlation substrate
only in windows that contain the common GW signal, not in noise-only
windows. If it sits there everywhere, or nowhere, that is informative
against the toy-model reading.

Speaking posture (Rule C, BEFORE running): I THINK the merger windows
MAY show the channel-correlation substrate that noise windows do not,
but we wait on what the data says. A null result (no difference) is a
real data point against the toy reading, not a failure to explain away.
Literature-as-conjecture also load-bearing: GW150914 detection is
replicated/solid (meets the bar); we use the public strain directly and
draw our own extraction, not lean on published interpretation.

Reuses kbk_pipeline.compute_operator_matrix (windowed Vasicek H + hist
MI) verbatim -- the SAME extraction used on OU/caricatures.

Data: GWOSC public 32 s / 4096 Hz strain, GW150914 (GPSstart
1126259446; merger ~1126259462.4 -> ~16.4 s into the segment).
Canary: coarse stride, event region only. Full: whole 32 s.
"""
import sys
import json
import time
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch
from scipy.interpolate import interp1d

from kbk_pipeline import compute_operator_matrix, extract_v1, project_234, OP_NAMES

FS = 4096
GPS_START = 1126259446
MERGER_GPS = 1126259462.4
MERGER_T = MERGER_GPS - GPS_START          # ~16.4 s into segment
ATTR3 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
CANARY = "--canary" in sys.argv


def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35.0, hi=350.0, fs=FS, order=4):
    b, a = butter(order, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x, fs=FS):
    """Whiten by the Welch ASD (standard GWOSC approach)."""
    n = len(x)
    freqs, psd = welch(x, fs=fs, nperseg=min(4 * fs, n))
    interp_psd = interp1d(freqs, psd, bounds_error=False,
                          fill_value=(psd[0], psd[-1]))
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    asd = np.sqrt(np.maximum(interp_psd(f), 1e-50))
    Xw = X / asd
    xw = np.fft.irfft(Xw, n=n)
    return xw / np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


# Whitening + filtfilt corrupt the segment edges (large transients).
# Crop EDGE_S seconds from each end before windowing so a filter edge
# artifact cannot masquerade as signal (the canary's peak MI landed at
# t=0.06s, a segment edge -- exactly this artifact).
EDGE_S = 2.0


def windowed_matrix(a, b, win_s, stride_s):
    win = int(win_s * FS)
    stride = int(stride_s * FS)
    M = compute_operator_matrix(a, b, win, stride, mi_bins=16)
    # window center times (s into segment)
    n_w = M.shape[0]
    centers = (np.arange(n_w) * stride + win / 2) / FS
    return M, centers


def null_metrics(M):
    v, S, _, _ = extract_v1(M)
    sub, cos_a = project_234(v)
    return v, float(abs(cos_a)), float(abs(v[5]))   # |cos to attractor|, |MI coef|


def analyze_block(a, b, centers_mask, label):
    """Extract the null over a subset of windows (mask over rows)."""
    return None  # placeholder, not used; kept for clarity


t0 = time.time()
print("=" * 100)
print(f"S10 LIGO real-data extraction {'[CANARY]' if CANARY else '[FULL]'}  "
      f"(GW150914, merger ~{MERGER_T:.1f}s into 32s segment)")
print("Test: does the operator null sit on the channel-correlation substrate")
print("(-1,-1,+2)/sqrt6 in WINDOWS WITH the common GW signal vs NOISE-ONLY windows?")
print("=" * 100)

h1 = load_strain("data/ligo/H-H1_GW150914_32s.hdf5")
l1 = load_strain("data/ligo/L-L1_GW150914_32s.hdf5")
print(f"loaded H1 {len(h1)} samples, L1 {len(l1)} samples ({len(h1)/FS:.0f}s @ {FS}Hz)")

print("preprocessing (bandpass 35-350 Hz + whiten)...")
h1p = preprocess(h1)
l1p = preprocess(l1)

# GW150914: L1 leads H1 by ~7 ms and is inverted. Align so a common
# signal shows as positive correlation (helps MI/channel-corr read).
shift = int(round(0.0069 * FS))
l1a = -np.roll(l1p, shift)

# crop filter-corrupted edges; keep a common time origin for centers
crop = int(EDGE_S * FS)
h1c = h1p[crop:-crop]
l1c = l1a[crop:-crop]
print(f"cropped {EDGE_S}s from each edge -> {len(h1c)/FS:.1f}s analysed")

win_s = 0.125
stride_s = 0.5 if CANARY else 0.03125  # 500ms canary / 31ms full
M, centers = windowed_matrix(h1c, l1c, win_s, stride_s)
centers = centers + EDGE_S   # shift back to segment-relative time
print(f"operator matrix: {M.shape[0]} windows of {win_s*1000:.0f}ms "
      f"(stride {stride_s*1000:.0f}ms)")

# Split windows into "event" (within +/-0.5 s of merger) and "noise"
event_mask = np.abs(centers - MERGER_T) <= 0.5
noise_mask = ~event_mask

results = {}
for label, mask in [("ALL", np.ones_like(event_mask)),
                    ("EVENT (+/-0.5s of merger)", event_mask),
                    ("NOISE-ONLY", noise_mask)]:
    Msub = M[mask]
    if Msub.shape[0] < 6:
        print(f"\n{label}: too few windows ({Msub.shape[0]})")
        continue
    v, cos_a, mi_coef = null_metrics(Msub)
    # also per-window MI mean as a sanity time-series stat
    mi_mean = float(Msub[:, 5].mean())
    mi_in_event = float(M[event_mask, 5].mean()) if event_mask.sum() else float("nan")
    print(f"\n{label}: {Msub.shape[0]} windows")
    print(f"  null 6D = {[f'{x:+.3f}' for x in v]}")
    print(f"  |cos to channel-corr attractor (-1,-1,+2)| = {cos_a:.3f}")
    print(f"  |MI coef in null| = {mi_coef:.3f}   mean MI in block = {mi_mean:.4f}")
    results[label] = dict(n_windows=int(Msub.shape[0]), null_6d=v.tolist(),
                          cos_to_attractor=cos_a, mi_coef=mi_coef,
                          mean_MI=mi_mean)

# peak-MI window (does MI peak at the merger?)
ipk = int(np.argmax(M[:, 5]))
print(f"\npeak-MI window at t={centers[ipk]:.2f}s (merger at {MERGER_T:.1f}s), "
      f"MI={M[ipk,5]:.4f}")
results["peak_MI_time_s"] = float(centers[ipk])
results["merger_t_s"] = float(MERGER_T)

out = dict(meta=dict(canary=CANARY, event="GW150914", fs=FS, win_s=win_s,
                     stride_s=stride_s, elapsed_s=round(time.time() - t0, 1)),
           results=results)
fn = "s10_ligo_extract_canary.json" if CANARY else "s10_ligo_extract_results.json"
json.dump(out, open(fn, "w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
