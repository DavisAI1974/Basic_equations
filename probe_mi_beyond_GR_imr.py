"""
PROBE O1-real (S17): does inter-detector MI carry common structure BEYOND the FULL GR
waveform? The rigorous beyond-GR test O1 (INFO-054) said was needed -- now with IMR
matched-filter templates (pycbc IMRPhenomD), not just the Newtonian inspiral.

Logic: remove the maximum-likelihood GR waveform from each GW150914 detector, then
recompute inter-detector MI on the residuals vs a time-slide null.
  residual MI -> null            => MI is FULLY the GR waveform (no new structure;
                                    closes the new-law road on this axis).
  residual MI significant at the
  physical 7 ms lag after FULL
  IMR subtraction                => candidate structure BEYOND GR -> then must control
                                    for template mismatch / calibration / glitches
                                    before any claim. NOT a new law on first look.

Method (parallels probe_mi_beyond_chirp, full IMR template instead of Newtonian):
  - ASD per detector from the raw 32 s (welch); whiten data AND the IMRPhenomD template
    with the SAME detector ASD; bandpass both (identical pipeline).
  - Template (m1=36,m2=29 Msun, IMRPhenomD) placed with its merger at the data merger
    sample; quadratures = template and its Hilbert 90-deg copy (covers phase, and L1's
    inversion is just a sign of the fit coefficient).
  - Per detector: fine time-shift search; least-squares fit data ~ a*cos + b*sin over
    the merger window (maximum-likelihood amplitude+phase); subtract -> residual.
  - MI(full) vs MI(residual): lag scan (+/-20 ms) + time-slide null. MI is scale-
    invariant (INFO-051) so amplitude normalization is irrelevant.

Speaking posture (Rule C): I expect the full IMR template to remove much more than the
Newtonian inspiral did; if residual MI drops to the null that is the clean "MI = GR
waveform" outcome (no new law). A surviving residual is the interesting case but is NOT
a new-law claim on first look -- template mismatch is the leading deflationary reading.
No verdict in advance; a messy result is data.

Run:  python probe_mi_beyond_GR_imr.py
"""
import json
import time
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch, hilbert
from scipy.interpolate import interp1d
from pycbc.waveform import get_td_waveform

from kbk_pipeline import mi_hist_2d

FS = 4096
MERGER = {"H1": 16.4236, "L1": 16.4187}   # INFO-052 recovered t_c per detector
BINS = 16


def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35.0, hi=350.0, order=4):
    b, a = butter(order, [lo / (FS / 2), hi / (FS / 2)], btype="band")
    return filtfilt(b, a, x)


def asd_fn_from(x):
    fr, psd = welch(x, fs=FS, nperseg=min(4 * FS, len(x)))
    return interp1d(fr, psd, bounds_error=False, fill_value=(psd[0], psd[-1]))


def whiten_with(x, asd_fn):
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), 1.0 / FS)
    asd = np.sqrt(np.maximum(asd_fn(f), 1e-50))
    return np.fft.irfft(X / asd, n=len(x))


def imr_template():
    hp, _ = get_td_waveform(approximant="IMRPhenomD", mass1=36, mass2=29,
                            spin1z=0.0, spin2z=0.0, delta_t=1.0 / FS,
                            f_lower=30.0)
    t = np.asarray(hp.numpy(), dtype=float)
    return t, int(np.argmax(np.abs(t)))


def place_template(tmpl, peak_idx, total_len, merger_sample):
    """Embed template in a zero array of total_len with its peak at merger_sample."""
    out = np.zeros(total_len)
    start = merger_sample - peak_idx
    s0 = max(0, start); s1 = min(total_len, start + len(tmpl))
    t0 = s0 - start; t1 = t0 + (s1 - s0)
    out[s0:s1] = tmpl[t0:t1]
    return out


def matched_subtract(d_white, cos_full, sin_full, merger_sample, fit_win_s=0.30,
                     search_ms=30.0):
    """Fine time-shift search + 2-quadrature lstsq fit over the merger window;
    return full-length residual + diagnostics."""
    half = int(fit_win_s * FS / 2)
    w = slice(merger_sample - half, merger_sample + half)
    best = None
    for sh in range(-int(search_ms / 1000 * FS), int(search_ms / 1000 * FS) + 1, 4):
        c = np.roll(cos_full, sh); s = np.roll(sin_full, sh)
        A = np.vstack([c[w], s[w]]).T
        coef, *_ = np.linalg.lstsq(A, d_white[w], rcond=None)
        resid = d_white[w] - A @ coef
        ss = float(np.sum(resid ** 2))
        if best is None or ss < best[0]:
            best = (ss, sh, coef)
    ss, sh, coef = best
    c = np.roll(cos_full, sh); s = np.roll(sin_full, sh)
    model = coef[0] * c + coef[1] * s
    resid_full = d_white - model
    var_removed = 1 - np.var(resid_full[w]) / np.var(d_white[w])
    return resid_full, float(var_removed), int(sh)


def mi_lag(a, b, lag):
    if lag >= 0:
        aa, bb = a[lag:], b[:len(b) - lag]
    else:
        aa, bb = a[:len(a) + lag], b[-lag:]
    n = min(len(aa), len(bb))
    return mi_hist_2d(aa[:n], bb[:n], bins=BINS)


def mi_window_lags(hser, lser, c_t, win_s, lags):
    c = int(c_t * FS); half = int(win_s * FS / 2)
    hw = hser[c - half:c + half]
    out = []
    for lg in lags:
        lw = lser[c - half - lg:c + half - lg]
        n = min(len(hw), len(lw))
        out.append(mi_hist_2d(hw[:n], lw[:n], bins=BINS))
    return np.array(out)


def time_slide_null(hser, lser, c_t, win_s):
    c = int(c_t * FS); half = int(win_s * FS / 2)
    hw = hser[c - half:c + half]
    null = []
    for off_s in [0.2, 0.35, 0.5, 0.65, 0.8, -0.2, -0.35, -0.5, -0.65, -0.8]:
        off = int(off_s * FS)
        lw = lser[c - half + off:c + half + off]
        n = min(len(hw), len(lw))
        null.append(mi_hist_2d(hw[:n], lw[:n], bins=BINS))
    return np.array(null)


def main():
    t0 = time.time()
    raw = {d: load_strain(p) for d, p in
           [("H1", "data/ligo/H-H1_GW150914_32s.hdf5"),
            ("L1", "data/ligo/L-L1_GW150914_32s.hdf5")]}
    N = len(raw["H1"])
    tmpl, peak_idx = imr_template()

    white, resid = {}, {}
    diag = {}
    for d in ("H1", "L1"):
        asd_fn = asd_fn_from(raw[d])
        wd = bandpass(whiten_with(bandpass(raw[d]), asd_fn))
        # whiten the merger-aligned template with the SAME detector ASD
        tpl_full = place_template(tmpl, peak_idx, N, int(MERGER[d] * FS))
        tcos = bandpass(whiten_with(tpl_full, asd_fn))
        tsin = np.imag(hilbert(tcos))
        rfull, vr, sh = matched_subtract(wd, tcos, tsin, int(MERGER[d] * FS))
        white[d] = wd; resid[d] = rfull
        diag[d] = dict(var_removed=vr, shift_samples=sh)
        print(f"{d}: IMR template removed {vr*100:.1f}% of merger-window variance "
              f"(time shift {sh} samp = {sh/FS*1000:+.1f} ms)", flush=True)

    c_t, win_s = MERGER["H1"], 0.15
    lags_ms = np.arange(-20, 20.1, 0.5)
    lags = (lags_ms / 1000 * FS).round().astype(int)

    mi_full = mi_window_lags(white["H1"], white["L1"], c_t, win_s, lags)
    mi_res = mi_window_lags(resid["H1"], resid["L1"], c_t, win_s, lags)
    jf, jr = int(np.argmax(mi_full)), int(np.argmax(mi_res))

    null_full = time_slide_null(white["H1"], white["L1"], c_t, win_s)
    null_res = time_slide_null(resid["H1"], resid["L1"], c_t, win_s)
    zf = (mi_full[jf] - null_full.mean()) / (null_full.std() + 1e-12)
    zr = (mi_res[jr] - null_res.mean()) / (null_res.std() + 1e-12)
    i7 = int(np.argmin(np.abs(lags_ms - 7)))

    print("=" * 92)
    print("PROBE beyond-GR (IMR matched-filter subtraction) -- GW150914")
    print("=" * 92)
    print(f"FULL     : peak MI {mi_full[jf]:.4f} at {lags_ms[jf]:+.1f} ms  "
          f"(null {null_full.mean():.4f}+/-{null_full.std():.4f}, z={zf:.1f})")
    print(f"RESIDUAL : peak MI {mi_res[jr]:.4f} at {lags_ms[jr]:+.1f} ms  "
          f"(null {null_res.mean():.4f}+/-{null_res.std():.4f}, z={zr:.1f})")
    print(f"MI at physical +7 ms: full {mi_full[i7]:.4f} -> residual {mi_res[i7]:.4f}")
    print(f"residual excess over its null: {mi_res[jr]-null_res.mean():+.4f}  "
          f"(=> {'STILL ABOVE null' if zr > 3 else 'consistent with null'})",
          flush=True)

    out = dict(meta=dict(approximant="IMRPhenomD", masses=[36, 29], bins=BINS,
                         win_s=win_s, elapsed_s=round(time.time() - t0, 1)),
               template_fit=diag, lags_ms=lags_ms.tolist(),
               mi_full=mi_full.tolist(), mi_residual=mi_res.tolist(),
               full=dict(peak_mi=float(mi_full[jf]), peak_lag_ms=float(lags_ms[jf]),
                         null_mean=float(null_full.mean()),
                         null_std=float(null_full.std()), z=float(zf)),
               residual=dict(peak_mi=float(mi_res[jr]), peak_lag_ms=float(lags_ms[jr]),
                             null_mean=float(null_res.mean()),
                             null_std=float(null_res.std()), z=float(zr),
                             mi_at_7ms=float(mi_res[i7])))
    json.dump(out, open("probe_mi_beyond_GR_imr_results.json", "w"), indent=2)
    print(f"\nWrote probe_mi_beyond_GR_imr_results.json  ({out['meta']['elapsed_s']}s)",
          flush=True)


if __name__ == "__main__":
    main()
