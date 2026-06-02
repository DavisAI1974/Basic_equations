"""
PROBE (backlog #3, S17): recover gravity's governing law -- the inspiral chirp --
from REAL strain by a time-frequency RIDGE, not the per-sample Hilbert that failed
in S16.

WHY this is gravity-FORWARD (not edge): S17 established the equal-entropy substrate
is bookkeeping; the only genuinely physical content is scale-invariant. A governing
LAW recovered from raw data is exactly that -- the gravity analogue of the S16 WF
(Z propagator M_Z 99.5%) and SF (QCD asymptotic freedom) wins. The Newtonian
inspiral predicts:
    df/dt = (96/5) pi^(8/3) (G M_c / c^3)^(5/3) f^(11/3)
so u(t) = f(t)^(-8/3) is LINEAR in t with slope
    k = -(256/5) pi^(8/3) (G M_c / c^3)^(5/3)
=> chirp mass  G M_c/c^3 = [ |k| * 5 / (256 * pi^(8/3)) ]^(3/5)  [seconds];
   M_c [Msun] = (G M_c/c^3) / 4.925491e-6.
Falsifiable targets (independent measured truth, treated as conjecture but solid):
  GW150914 detector-frame chirp mass ~ 30.4 Msun;  GW170817 ~ 1.197 Msun (BNS).
A clean linear u-vs-t (high R^2) with M_c in the right ballpark = the law recovered.

METHOD: own FFT-based Morlet CWT (scipy.signal.cwt removed in 1.17) -> scalogram
|W|(f,t); ridge = per-time peak frequency (parabolic-interpolated) above an SNR
floor; fit u=f^(-8/3) vs t over the rising-frequency inspiral up to merger
(t_c = ridge-frequency peak). Reports R^2 and recovered M_c per event + per detector.

Speaking posture (Rule C): the ridge should track the chirp far better than the
S16 Hilbert; I expect roughly linear u-vs-t with M_c in the right ballpark, but real
whitened strain is noisy and the ridge may wander -- wait on the data, no verdict on
first look, and a messy ridge is still data (don't predetermine good/bad).

Run:  python probe_gravity_chirp_ridge.py --canary
      python probe_gravity_chirp_ridge.py
"""
import sys
import json
import time
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch
from scipy.interpolate import interp1d

FS = 4096
CANARY = "--canary" in sys.argv
GM_SUN_OVER_C3 = 4.925491e-6   # seconds


def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=30.0, hi=400.0, order=4):
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


def morlet_cwt(x, freqs, n_cyc=6.0):
    """FFT-based Morlet CWT. Returns |W| of shape (len(freqs), len(x)).

    For center frequency f the wavelet is
        psi_f(t) = exp(2 pi i f t) exp(-t^2 / (2 sigma_t^2)),  sigma_t = n_cyc/(2 pi f)
    normalized to unit L2. Convolution done in the frequency domain.
    """
    n = len(x)
    Xf = np.fft.fft(x)
    out = np.empty((len(freqs), n))
    t = (np.arange(n) - n // 2) / FS
    for i, f in enumerate(freqs):
        sigma_t = n_cyc / (2 * np.pi * f)
        psi = np.exp(2j * np.pi * f * t) * np.exp(-t * t / (2 * sigma_t * sigma_t))
        psi /= np.linalg.norm(psi)
        Psi = np.fft.fft(np.fft.ifftshift(psi))
        W = np.fft.ifft(Xf * np.conj(Psi))
        out[i] = np.abs(W)
    return out


def _peak_logf(col, logf, i0=None, band=None):
    """Peak (parabolic-interp) log-frequency of a column, optionally restricted
    to indices within `band` of the previous index i0."""
    if i0 is not None and band is not None:
        lo = max(0, i0 - band)
        hi = min(len(col), i0 + band + 1)
        seg = col[lo:hi]
        ipk = lo + int(np.argmax(seg))
    else:
        ipk = int(np.argmax(col))
    if 0 < ipk < len(col) - 1:
        y0, y1, y2 = col[ipk - 1], col[ipk], col[ipk + 1]
        denom = (y0 - 2 * y1 + y2)
        delta = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        delta = float(np.clip(delta, -1, 1))
        lf = logf[ipk] + delta * (logf[1] - logf[0])
    else:
        lf = logf[ipk]
    return ipk, lf


def ridge_track_backward(scal, freqs, times, i_merger, snr_floor=3.0,
                         band_frac=0.25):
    """Continuity-constrained ridge, tracked BACKWARD from the merger column.

    Start at the merger column's global peak; walk to earlier times, each step
    restricting the search to frequencies within band_frac (in log-index) of the
    previous ridge point, so the track follows the descending inspiral arc and
    does not jump to unrelated noise peaks. Stop when column SNR < floor.
    Returns (t_ridge, f_ridge) ordered in increasing time.
    """
    logf = np.log(freqs)
    band = max(2, int(band_frac * len(freqs)))
    t_r, f_r = [], []
    # seed at merger
    col = scal[:, i_merger]
    i_prev, lf = _peak_logf(col, logf)
    t_r.append(times[i_merger]); f_r.append(np.exp(lf))
    gaps = 0
    for j in range(i_merger - 1, -1, -1):
        col = scal[:, j]
        med = np.median(col) + 1e-30
        ipk, lf = _peak_logf(col, logf, i0=i_prev, band=band)
        if col[ipk] / med < snr_floor:
            gaps += 1
            if gaps > 12:          # tolerate brief dropouts, then stop
                break
            continue
        gaps = 0
        i_prev = ipk
        t_r.append(times[j]); f_r.append(np.exp(lf))
    order = np.argsort(t_r)
    return np.array(t_r)[order], np.array(f_r)[order]


def fit_chirp(t_r, f_r):
    """Fit u = f^(-8/3) linear in t over the inspiral arc (tracked to merger).

    The backward tracker already returns the rising inspiral arc; t_c = last
    (merger) time. Drop the final ~merger/ringdown point and keep the
    monotone-rising body, then linear-fit u vs t. Returns R^2, slope, M_c.
    """
    if len(t_r) < 8:
        return {"ok": False, "reason": "too few ridge points"}
    t_c = float(t_r[-1])
    # drop the merger column itself (ringdown contamination), fit the inspiral
    ti, fi = t_r[:-1], f_r[:-1]
    if len(ti) < 8:
        return {"ok": False, "reason": "too few inspiral points", "t_c": t_c}
    # keep the monotone-rising body (drop tracker dropouts)
    keep = fi >= (np.maximum.accumulate(fi) * 0.6)
    ti, fi = ti[keep], fi[keep]
    u = fi ** (-8.0 / 3.0)
    A = np.vstack([ti, np.ones_like(ti)]).T
    coef, *_ = np.linalg.lstsq(A, u, rcond=None)
    k, b = float(coef[0]), float(coef[1])
    pred = A @ coef
    ss_res = float(np.sum((u - pred) ** 2))
    ss_tot = float(np.sum((u - u.mean()) ** 2)) + 1e-300
    r2 = 1.0 - ss_res / ss_tot
    # chirp mass from |slope| (u decreases as t->t_c, so k<0 for a real chirp)
    Mc = None
    if k < 0:
        gmc = (abs(k) * 5.0 / (256.0 * np.pi ** (8.0 / 3.0))) ** (3.0 / 5.0)
        Mc = float(gmc / GM_SUN_OVER_C3)
    return {"ok": True, "t_c": t_c, "n_inspiral": int(len(ti)),
            "slope_k": k, "intercept": b, "r2": r2,
            "f_lo": float(fi.min()), "f_hi": float(fi.max()),
            "chirp_mass_Msun": Mc}


def analyze_event(name, h1_path, l1_path, freqs, n_cyc, t_lo, t_hi, snr,
                  merger_t=None):
    res = {}
    for det, path in [("H1", h1_path), ("L1", l1_path)]:
        if path is None:
            continue
        x = preprocess(load_strain(path))
        crop = int(2.0 * FS)
        x = x[crop:-crop]
        seg_t0 = 2.0
        n = len(x)
        full_t = seg_t0 + np.arange(n) / FS
        # restrict to the analysis window before CWT (speed + focus)
        sel = (full_t >= t_lo) & (full_t <= t_hi)
        xa = x[sel]
        ta = full_t[sel]
        scal = morlet_cwt(xa, freqs, n_cyc=n_cyc)
        if merger_t is not None:
            # seed at the known merger column (robust for faint signals)
            i_merger = int(np.argmin(np.abs(ta - merger_t)))
        else:
            # merger column = max broadband power (sum over freqs), smoothed
            bb = scal.sum(axis=0)
            kk = max(1, len(bb) // 200)
            bb_s = np.convolve(bb, np.ones(kk) / kk, mode="same")
            i_merger = int(np.argmax(bb_s))
        t_r, f_r = ridge_track_backward(scal, freqs, ta, i_merger,
                                        snr_floor=snr)
        fit = fit_chirp(t_r, f_r)
        res[det] = {"n_ridge": int(len(t_r)), "fit": fit,
                    "ridge_t": [float(v) for v in t_r],
                    "ridge_f": [float(v) for v in f_r]}
        mc = fit.get("chirp_mass_Msun")
        print(f"  {name} {det}: ridge pts={len(t_r)}  "
              f"R^2={fit.get('r2', float('nan')):.3f}  "
              f"M_c={('%.2f' % mc) if mc else 'n/a'} Msun  "
              f"f {fit.get('f_lo', 0):.0f}->{fit.get('f_hi', 0):.0f} Hz  "
              f"t_c={fit.get('t_c', float('nan')):.3f}s", flush=True)
    return res


def find_merger_t(path):
    """Coarse merger time = peak of |whitened strain| envelope (smoothed)."""
    x = preprocess(load_strain(path))
    crop = int(2.0 * FS)
    x = x[crop:-crop]
    env = np.abs(x)
    k = FS // 50
    sm = np.convolve(env, np.ones(k) / k, mode="same")
    i = int(np.argmax(sm))
    return 2.0 + i / FS


def main():
    t0 = time.time()
    freqs = np.geomspace(30, 400, 60 if CANARY else 120)

    print("=" * 96)
    print(f"PROBE gravity chirp ridge {'[CANARY]' if CANARY else '[FULL]'}  "
          "(recover the inspiral law f^(-8/3) ~ linear in t from real strain)")
    print("=" * 96, flush=True)

    out = {}

    # GW150914: fast BBH chirp; merger ~16.4s; analyze a tight window, few cycles
    print("\nGW150914 (BBH, fast chirp):", flush=True)
    out["GW150914"] = analyze_event(
        "GW150914",
        "data/ligo/H-H1_GW150914_32s.hdf5",
        "data/ligo/L-L1_GW150914_32s.hdf5",
        freqs, n_cyc=4.0, t_lo=16.0, t_hi=16.50, snr=2.5)

    if not CANARY:
        # GW170817: long BNS inspiral. Merger from file metadata + known GPS
        # (envelope-peak detection is unreliable for the faint H1 BNS signal).
        print("\nGW170817 (BNS, long inspiral):", flush=True)
        with h5py.File("data/ligo/H-H1_GW170817_32s.hdf5", "r") as f:
            gps0 = float(f["meta/GPSstart"][()])
        mt = 1187008882.43 - gps0   # known GW170817 merger GPS -> segment time
        print(f"  merger from metadata ~{mt:.2f}s into segment", flush=True)
        out["GW170817"] = analyze_event(
            "GW170817",
            "data/ligo/H-H1_GW170817_32s.hdf5",
            None,  # GWOSC 32s file is H1 only here
            freqs, n_cyc=8.0, t_lo=max(2.1, mt - 12.0), t_hi=min(29.9, mt + 0.3),
            snr=1.6, merger_t=mt)

    result = dict(meta=dict(canary=CANARY, freqs_hz=[float(freqs[0]), float(freqs[-1])],
                            n_freqs=len(freqs), elapsed_s=round(time.time() - t0, 1)),
                  events=out)
    fn = "probe_gravity_chirp_ridge_canary.json" if CANARY \
        else "probe_gravity_chirp_ridge_results.json"
    json.dump(result, open(fn, "w"), indent=2)
    print(f"\nWrote {fn}  ({result['meta']['elapsed_s']}s)", flush=True)


if __name__ == "__main__":
    main()
