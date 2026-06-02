"""
S11 canary: does the S10 LIGO pipeline GENERALIZE to a second real event?
Reuses s10_ligo_extract preprocessing + kbk_pipeline extraction VERBATIM.
Only change: per-event metadata + GENERIC cross-correlation alignment
(GW150914's 7 ms lead / inversion was sky-location specific; we must not
hardcode it for other events). This validates mechanics on new data
before the full 12-event batch + off-source null distribution.

Speaking posture (Rule C): I THINK the method will RUN and produce sane
output on GW151226, but GW151226 is a quiet, long-inspiral event (much
softer than GW150914's loud merger), so a weak/absent MI-at-merger peak
here is informative, not a failure. Significance needs the null
distribution + louder events; this canary only asks "does it run + is
the output sane."
"""
import sys, json, time
import numpy as np
from scipy.signal import correlate
from s10_ligo_extract import load_strain, preprocess, windowed_matrix, null_metrics, FS, EDGE_S

# --- GW151226 metadata (GWOSC) ---
EVENT = "GW151226"
SEG_START = 1135136334          # 32s clip start GPS
MERGER_GPS = 1135136350.6
MERGER_T = MERGER_GPS - SEG_START   # ~16.6 s into segment
H1_PATH = "data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5"
L1_PATH = "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5"


def generic_align(h1p, l1p, max_lag_s=0.015):
    """Find integer lag + sign that maximise |xcorr| in a window around the
    merger. Returns aligned L1 so a common signal reads as +correlation."""
    c = int(MERGER_T * FS)
    half = int(0.25 * FS)
    a = h1p[c - half:c + half]
    b = l1p[c - half:c + half]
    xc = correlate(a, b, mode="full")
    lags = np.arange(-len(b) + 1, len(a))
    m = np.abs(lags) <= int(max_lag_s * FS)
    best = lags[m][np.argmax(np.abs(xc[m]))]
    sign = np.sign(xc[m][np.argmax(np.abs(xc[m]))])
    return np.roll(l1p, best) * sign, int(best), float(sign)


t0 = time.time()
print("=" * 90)
print(f"S11 LIGO CANARY: generalization test on {EVENT} (merger ~{MERGER_T:.1f}s)")
print("=" * 90)
h1, l1 = load_strain(H1_PATH), load_strain(L1_PATH)
print(f"loaded H1 {len(h1)}, L1 {len(l1)} samples ({len(h1)/FS:.0f}s @ {FS}Hz)")
h1p, l1p = preprocess(h1), preprocess(l1)
l1a, lag, sign = generic_align(h1p, l1p)
print(f"generic alignment: lag={lag} samp ({lag/FS*1e3:+.1f} ms), sign={sign:+.0f}")

crop = int(EDGE_S * FS)
h1c, l1c = h1p[crop:-crop], l1a[crop:-crop]
M, centers = windowed_matrix(h1c, l1c, 0.125, 0.03125)
centers = centers + EDGE_S

event_mask = np.abs(centers - MERGER_T) <= 0.5
noise_mask = ~event_mask
res = {}
for label, mask in [("ALL", np.ones_like(event_mask)),
                    ("EVENT (+/-0.5s)", event_mask), ("NOISE-ONLY", noise_mask)]:
    Msub = M[mask]
    if Msub.shape[0] < 6:
        continue
    v, cos_a, mi = null_metrics(Msub)
    res[label] = dict(n=int(Msub.shape[0]), cos_to_attractor=cos_a,
                      mi_coef=mi, mean_MI=float(Msub[:, 5].mean()))
    print(f"\n{label}: {Msub.shape[0]} windows")
    print(f"  |cos to (-1,-1,+2) attractor| = {cos_a:.3f}   mean MI = {Msub[:,5].mean():.4f}")

ipk = int(np.argmax(M[:, 5]))
mi_noise = float(M[noise_mask, 5].mean())
mi_evt = float(M[event_mask, 5].mean())
print(f"\npeak-MI window at t={centers[ipk]:.2f}s (merger {MERGER_T:.1f}s); "
      f"MI peak={M[ipk,5]:.4f}")
print(f"mean MI: event={mi_evt:.4f} vs noise={mi_noise:.4f} "
      f"({mi_evt/mi_noise:.2f}x)")
res["peak_MI_time_s"] = float(centers[ipk])
res["merger_t_s"] = MERGER_T
res["mi_event_over_noise"] = mi_evt / mi_noise
json.dump(dict(meta=dict(event=EVENT, align_lag=lag, align_sign=sign,
                         elapsed_s=round(time.time() - t0, 1)), results=res),
          open("s11_ligo_canary.json", "w"), indent=2)
print(f"\nWrote s11_ligo_canary.json ({time.time()-t0:.1f}s)")
