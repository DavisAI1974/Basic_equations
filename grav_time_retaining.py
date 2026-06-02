"""
Time-retaining extraction on REAL black-hole merger data (LIGO).

Tests the within-gravity half of Greg's S15-eve prediction: if gravity's
distinction lives in TIME (its flow/info component = its time-effect), then
(1) the LEVEL null we normally extract is time-BLIND and cannot see it, and
(2) when we RETAIN time, gravity's information/flow reveals itself as
time-LOCALIZED, in a way a stationary object would not.

We have real strain for 3 mergers (GW150914 + two GPS-tagged segments). No
gauge-force objects on this branch, so this is the gravity SIDE only --
within-construction, NO synthesis with the gauge forces (that comparison is
S15 real-data work). The honest deflationary reading is stated in the readout:
a merger is a transient, so time-localized flow may be a construction-type fact
(transient time-series vs stationary ensemble), not gravity ontology.

Four questions per event:
  A. ORDER-BLINDNESS: shuffle the time-rows of the operator matrix, re-extract.
     cos should be exactly 1 -> the level-null cannot see time-order, on real
     data too. (Confirms why the 4/4 self-pole result can't speak to gravity's
     time-specialness.)
  B. FLOW TIME-LOCALIZATION: is MI(t) a time-localized spike at the merger?
     peak-MI time vs merger; peak/median ratio; excess-MI fraction in the
     merger window; kurtosis of MI(t).
  C. MERGER-MOMENT DEPENDENCE: does removing the time-localized merger window
     change MI's variance and its participation in the null? (INFO-041: MI is
     "self-pole" because the merger spike makes it high-variance.)
  D. INFORMATION RATE: is dMI/dt itself a time-localized event (sharp spike at
     the merger)?

Reuses the S10 LIGO preprocessing + the kbk operator stack verbatim.
"""
import json
import time

import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch
from scipy.interpolate import interp1d

from kbk_pipeline import compute_operator_matrix, extract_v1, project_234, OP_NAMES

FS = 4096
EDGE_S = 2.0
WIN_S = 0.125
STRIDE_S = 0.03125

# (H1 file, L1 file, label, merger_t_into_segment_s, l1_shift_s, l1_invert)
EVENTS = [
    ("data/ligo/H-H1_GW150914_32s.hdf5", "data/ligo/L-L1_GW150914_32s.hdf5",
     "GW150914", 16.4, 0.0069, True),
    ("data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5",
     "GW170104(approx-align)", 16.6, 0.0, True),
    ("data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5",
     "GW151226(approx-align)", 16.6, 0.0, True),
]


def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35.0, hi=350.0, fs=FS, order=4):
    b, a = butter(order, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x, fs=FS):
    n = len(x)
    freqs, psd = welch(x, fs=fs, nperseg=min(4 * fs, n))
    interp_psd = interp1d(freqs, psd, bounds_error=False,
                          fill_value=(psd[0], psd[-1]))
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    asd = np.sqrt(np.maximum(interp_psd(f), 1e-50))
    xw = np.fft.irfft(X / asd, n=n)
    return xw / np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def cos(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else float("nan")


def kurtosis(x):
    x = np.asarray(x, float); m = x.mean(); s = x.std() + 1e-12
    return float(np.mean(((x - m) / s) ** 4) - 3.0)


def run_event(h1f, l1f, label, merger_t, shift_s, invert, seed=11):
    h1 = load_strain(h1f); l1 = load_strain(l1f)
    h1p = preprocess(h1); l1p = preprocess(l1)
    shift = int(round(shift_s * FS))
    l1a = (-1.0 if invert else 1.0) * np.roll(l1p, shift)
    crop = int(EDGE_S * FS)
    h1c, l1c = h1p[crop:-crop], l1a[crop:-crop]

    win = int(WIN_S * FS); stride = int(STRIDE_S * FS)
    M = compute_operator_matrix(h1c, l1c, win, stride, mi_bins=16)
    n_w = M.shape[0]
    centers = (np.arange(n_w) * stride + win / 2) / FS + EDGE_S
    mi = M[:, 5]

    # ---- A: order-blindness on real data ----
    v_full, _, _, _ = extract_v1(M)
    rng = np.random.default_rng(seed)
    v_shuf, _, _, _ = extract_v1(M[rng.permutation(n_w)])
    cos_order = abs(cos(v_full, v_shuf))

    # ---- B: flow time-localization ----
    ipk = int(np.argmax(mi))
    peak_t = float(centers[ipk])
    peak_over_median = float(mi[ipk] / (np.median(mi) + 1e-12))
    ev_mask = np.abs(centers - merger_t) <= 0.5
    noise_mask = ~ev_mask
    mi_floor = float(np.median(mi[noise_mask]))
    excess = np.maximum(mi - mi_floor, 0.0)
    excess_frac_in_merger = float(excess[ev_mask].sum() / (excess.sum() + 1e-12))
    mi_kurt = kurtosis(mi)

    # ---- C: merger-moment dependence of MI's null participation ----
    def mi_coef_and_var(mask):
        Msub = M[mask]
        v, _, _, _ = extract_v1(Msub)
        return float(abs(v[5])), float(Msub[:, 5].var())
    mi_coef_full, mi_var_full = float(abs(v_full[5])), float(mi.var())
    mi_coef_noise, mi_var_noise = mi_coef_and_var(noise_mask)
    _, cos_attr_full = project_234(v_full)

    # ---- D: information RATE dMI/dt ----
    dmi = np.diff(mi) / STRIDE_S
    idmax = int(np.argmax(np.abs(dmi)))
    drate_peak_t = float(centers[idmax])
    drate_kurt = kurtosis(dmi)

    return {
        "label": label, "n_windows": n_w, "merger_t_s": merger_t,
        "A_cos_order_invariance": round(cos_order, 6),
        "A_null_6d": [round(float(x), 3) for x in v_full],
        "B_peak_MI_t_s": round(peak_t, 3),
        "B_peak_over_median_MI": round(peak_over_median, 3),
        "B_excess_MI_frac_in_merger_window": round(excess_frac_in_merger, 3),
        "B_MI_kurtosis": round(mi_kurt, 2),
        "C_mi_coef_full": round(mi_coef_full, 3),
        "C_mi_coef_noise_only": round(mi_coef_noise, 3),
        "C_mi_var_full": float(f"{mi_var_full:.3e}"),
        "C_mi_var_noise_only": float(f"{mi_var_noise:.3e}"),
        "C_mi_var_ratio_full_over_noise": round(mi_var_full / (mi_var_noise + 1e-18), 2),
        "C_abs_cos_to_attractor_full": round(float(abs(cos_attr_full)), 3),
        "D_dMI_dt_peak_t_s": round(drate_peak_t, 3),
        "D_dMI_dt_kurtosis": round(drate_kurt, 2),
    }


def main():
    t0 = time.time()
    print(f"[grav_time_retaining] win={WIN_S*1000:.0f}ms stride={STRIDE_S*1000:.0f}ms "
          f"band 35-350Hz | OP_NAMES={OP_NAMES}", flush=True)
    out = []
    for h1f, l1f, label, mt, sh, inv in EVENTS:
        r = run_event(h1f, l1f, label, mt, sh, inv)
        out.append(r)
        print(f"\n=== {label} ({r['n_windows']} windows) ===", flush=True)
        print(f"  A order-shuffle cos = {r['A_cos_order_invariance']:.6f}  "
              f"(null time-blind on real data)", flush=True)
        print(f"  B peak-MI t={r['B_peak_MI_t_s']}s (merger {mt}s) | "
              f"peak/median={r['B_peak_over_median_MI']} | "
              f"excess-MI in merger window={r['B_excess_MI_frac_in_merger_window']} | "
              f"MI kurtosis={r['B_MI_kurtosis']}", flush=True)
        print(f"  C |MI coef| full={r['C_mi_coef_full']} noise-only={r['C_mi_coef_noise_only']} | "
              f"MI var full/noise ratio={r['C_mi_var_ratio_full_over_noise']} | "
              f"|cos attractor| full={r['C_abs_cos_to_attractor_full']}", flush=True)
        print(f"  D dMI/dt peak t={r['D_dMI_dt_peak_t_s']}s | "
              f"kurtosis={r['D_dMI_dt_kurtosis']}", flush=True)

    json.dump({"config": dict(win_s=WIN_S, stride_s=STRIDE_S, fs=FS, edge_s=EDGE_S),
               "events": out, "elapsed_s": round(time.time() - t0, 1)},
              open("grav_time_retaining_results.json", "w"), indent=2)
    print(f"\n[grav_time_retaining] wrote grav_time_retaining_results.json "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
