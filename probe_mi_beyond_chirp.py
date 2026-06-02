"""
PROBE O1 (S17): does the inter-detector MI carry common structure BEYOND the recovered
inspiral chirp? The only road to a NEW gravity law -- everything else this session
either validated known physics or cleared an artifact.

HONEST SCOPE: a rigorous "beyond GR" test needs full IMR matched-filter templates
(pycbc/lalsuite, not in our self-contained stack). This tests "beyond the recovered
INSPIRAL LAW (INFO-052)": fit + subtract the Newtonian chirp model we recovered from
each detector over the late inspiral (where the MI peak lives, T1 t~16.366 s < t_c),
then recompute the inter-detector MI on the RESIDUALS. If residual MI collapses to the
time-slide null => the recovered inspiral chirp accounts for ALL the common
information (no room for new structure in the inspiral). If residual MI persists at
the physical 7 ms lag => common structure remains -- but most plausibly merger/
ringdown (still GR), which only IMR templates can attribute. So a persisting residual
is NOT a new-law claim; a vanishing residual is a clean "MI = the inspiral chirp."

Model per detector (from INFO-052 fit u=f^(-8/3)=k t + b, f=u^(-3/8), t<t_c):
  phi(t) = 2 pi cumulative-integral of f;  Newtonian envelope e(t)=f^(2/3).
  Basis = e*cos(phi)*{1,t,t^2} and e*sin(phi)*{1,t,t^2}; least-squares fit to whitened
  strain over the inspiral window; subtract.

Speaking posture (Rule C): I expect subtracting the inspiral model to reduce MI but
likely not to zero (merger/ringdown remain in the wider window) -- and any residual is
GR, not new physics, absent templates. A residual that drops to the null is the clean
informative outcome. No verdict on first look; a messy result is data.

Run:  python probe_mi_beyond_chirp.py [--canary]
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
# INFO-052 recovered inspiral fit (u = k t + b ; f = u^(-3/8) ; t < t_c)
FIT = {"H1": dict(k=-6.7539e-04, b=1.1091e-02, t_c=16.4236),
       "L1": dict(k=-6.6908e-04, b=1.0983e-02, t_c=16.4187)}


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
    return np.fft.irfft(X / asd, n=n) / np.std(np.fft.irfft(X / asd, n=n))


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def chirp_basis(t, det):
    """Newtonian-envelope quadrature basis for the recovered inspiral chirp."""
    p = FIT[det]
    u = p["k"] * t + p["b"]
    valid = u > 1e-6
    f = np.where(valid, np.maximum(u, 1e-6) ** (-3.0 / 8.0), 0.0)
    dt = np.gradient(t)
    phi = 2 * np.pi * np.cumsum(f * dt)
    e = np.where(valid, f ** (2.0 / 3.0), 0.0)
    e = e / (np.max(e) + 1e-30)
    tc = t - t.mean()
    cols = []
    for k in range(3):
        cols.append(e * np.cos(phi) * tc ** k)
        cols.append(e * np.sin(phi) * tc ** k)
    return np.vstack(cols).T, valid


def remove_chirp(x, t, det):
    B, valid = chirp_basis(t, det)
    coef, *_ = np.linalg.lstsq(B[valid], x[valid], rcond=None)
    model = B @ coef
    return x - model, model, valid


def mi_lag(a, b, lag, bins=16):
    if lag >= 0:
        aa, bb = a[lag:], b[:len(b) - lag]
    else:
        aa, bb = a[:len(a) + lag], b[-lag:]
    n = min(len(aa), len(bb))
    return mi_hist_2d(aa[:n], bb[:n], bins=bins)


def main():
    t0 = time.time()
    h = preprocess(load_strain("data/ligo/H-H1_GW150914_32s.hdf5"))
    l = preprocess(load_strain("data/ligo/L-L1_GW150914_32s.hdf5"))

    # inspiral-centered window (MI peak T1 ~16.366 s, before t_c ~16.42 s)
    c_t = 16.37
    win_s = 0.12
    bins = 16
    c = int(c_t * FS); half = int(win_s * FS / 2)
    sl = slice(c - half, c + half)
    t_win = np.arange(c - half, c + half) / FS

    hw, lw = h[sl], l[sl]
    hr, hmod, hv = remove_chirp(hw, t_win, "H1")
    lr, lmod, lv = remove_chirp(lw, t_win, "L1")

    # how much of each detector's window variance the chirp model removed
    frac_h = 1 - np.var(hr) / np.var(hw)
    frac_l = 1 - np.var(lr) / np.var(lw)

    phys_lag = int(round(0.0069 * FS))   # ~7 ms
    lags_ms = np.arange(-20, 20.1, 0.5)
    lags = (lags_ms / 1000.0 * FS).round().astype(int)

    mi_full = np.array([mi_lag(hw, lw, lg, bins) for lg in lags])
    mi_res = np.array([mi_lag(hr, lr, lg, bins) for lg in lags])
    jf, jr = int(np.argmax(mi_full)), int(np.argmax(mi_res))

    # time-slide null: H1 merger window (full / residual) vs L1 windows pulled
    # from far-away times in the FULL array (noise) -- big offsets >> light travel.
    slide = np.array([int(s * FS) for s in [0.2, 0.35, 0.5, 0.65, 0.8]])
    slide = np.concatenate([slide, -slide])
    null_full, null_res = [], []
    for off in slide:
        lwin = l[c - half + off:c + half + off]
        n = min(len(hw), len(lwin))
        null_full.append(mi_hist_2d(hw[:n], lwin[:n], bins=bins))
        null_res.append(mi_hist_2d(hr[:n], lwin[:n], bins=bins))
    null_full = np.array(null_full); null_res = np.array(null_res)
    zf = (mi_full[jf] - null_full.mean()) / (null_full.std() + 1e-12)
    zr = (mi_res[jr] - null_res.mean()) / (null_res.std() + 1e-12)

    print("=" * 92)
    print(f"PROBE O1 MI beyond chirp {'[CANARY]' if CANARY else '[FULL]'} "
          "(GW150914; subtract recovered inspiral, MI of residual)")
    print("=" * 92)
    print(f"chirp model removed variance: H1 {frac_h*100:.1f}%  L1 {frac_l*100:.1f}%")
    print(f"FULL     : peak MI {mi_full[jf]:.4f} at lag {lags_ms[jf]:+.1f} ms  "
          f"(null {null_full.mean():.4f}+/-{null_full.std():.4f}, z={zf:.1f})")
    print(f"RESIDUAL : peak MI {mi_res[jr]:.4f} at lag {lags_ms[jr]:+.1f} ms  "
          f"(null {null_res.mean():.4f}+/-{null_res.std():.4f}, z={zr:.1f})")
    print(f"MI at physical +7ms lag: full {mi_full[np.argmin(np.abs(lags_ms-7))]:.4f}"
          f" -> residual {mi_res[np.argmin(np.abs(lags_ms-7))]:.4f}")
    drop = mi_full[jf] / (mi_res[jr] + 1e-12)
    print(f"peak-MI reduction after chirp removal: {drop:.2f}x  "
          f"(residual excess over null: {mi_res[jr]-null_res.mean():+.4f})", flush=True)

    out = dict(meta=dict(canary=CANARY, center_t=c_t, win_s=win_s, bins=bins,
                         elapsed_s=round(time.time()-t0, 1)),
               chirp_removed_var=dict(H1=float(frac_h), L1=float(frac_l)),
               lags_ms=lags_ms.tolist(),
               mi_full=mi_full.tolist(), mi_residual=mi_res.tolist(),
               full=dict(peak_mi=float(mi_full[jf]), peak_lag_ms=float(lags_ms[jf]),
                         null_mean=float(null_full.mean()),
                         null_std=float(null_full.std()), z=float(zf)),
               residual=dict(peak_mi=float(mi_res[jr]), peak_lag_ms=float(lags_ms[jr]),
                             null_mean=float(null_res.mean()),
                             null_std=float(null_res.std()), z=float(zr)),
               peak_reduction_x=float(drop))
    fn = "probe_mi_beyond_chirp_canary.json" if CANARY \
        else "probe_mi_beyond_chirp_results.json"
    json.dump(out, open(fn, "w"), indent=2)
    print(f"\nWrote {fn}  ({out['meta']['elapsed_s']}s)", flush=True)


if __name__ == "__main__":
    main()
