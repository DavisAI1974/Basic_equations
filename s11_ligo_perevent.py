"""
S11 per-event LIGO readout: run each event SEPARATELY, no merge / no average.
Greg's question: do the different per-event ENTROPY readings mean something?

So for each event we report, side by side but never pooled:
  - mean windowed entropy H_a (H1) and H_b (L1), in EVENT vs NOISE windows
  - entropy ASYMMETRY |H_a - H_b| (the (-1,-1,+2) attractor is the equal-
    entropy identity H_a~=H_b; the event is where it should break)
  - mean MI, peak-MI time vs merger, event/noise MI ratio
  - |cos to (-1,-1,+2) attractor| in EVENT vs NOISE windows

Self-contained: uses only kbk_pipeline pure functions (no s10/s11 import
side-effects). Runs the 3 events currently reachable via 32s clips; the
same per_event() generalizes to the full-12 bulk-fetch batch.

Speaking posture (Rule C): I THINK louder mergers (GW150914) will show a
bigger event/noise entropy-asymmetry break than quiet ones (GW151226),
but we wait on what the per-event numbers say -- and we do NOT average
them, because the spread across events IS the signal Greg is asking about.
"""
import json, numpy as np, h5py
from scipy.signal import butter, filtfilt, welch, correlate
from scipy.interpolate import interp1d
from kbk_pipeline import compute_operator_matrix, extract_v1, project_234

FS = 4096
ATTR3 = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)
EDGE_S = 2.0

EVENTS = [
    # name, h1, l1, merger_t (s into 32s segment)
    ("GW150914", "data/ligo/H-H1_GW150914_32s.hdf5",
                 "data/ligo/L-L1_GW150914_32s.hdf5", 16.4),
    ("GW151226", "data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
                 "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5", 16.6),
    ("GW170104", "data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
                 "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5", 16.6),
]


def load(p):
    with h5py.File(p, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35., hi=350., fs=FS, order=4):
    b, a = butter(order, [lo/(fs/2), hi/(fs/2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x, fs=FS):
    n = len(x); freqs, psd = welch(x, fs=fs, nperseg=min(4*fs, n))
    ip = interp1d(freqs, psd, bounds_error=False, fill_value=(psd[0], psd[-1]))
    X = np.fft.rfft(x); f = np.fft.rfftfreq(n, 1./fs)
    xw = np.fft.irfft(X/np.sqrt(np.maximum(ip(f), 1e-50)), n=n)
    return xw/np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def align(h, l, mt, max_lag_s=0.015):
    c = int(mt*FS); half = int(0.25*FS)
    xc = correlate(h[c-half:c+half], l[c-half:c+half], mode="full")
    lags = np.arange(-half*2+1, half*2)
    m = np.abs(lags) <= int(max_lag_s*FS)
    j = np.argmax(np.abs(xc[m]))
    return np.roll(l, lags[m][j])*np.sign(xc[m][j]), int(lags[m][j]), float(np.sign(xc[m][j]))


def per_event(name, h1p, l1p, mt):
    crop = int(EDGE_S*FS)
    a, b = h1p[crop:-crop], l1p[crop:-crop]
    win, stride = int(0.125*FS), int(0.03125*FS)
    M = compute_operator_matrix(a, b, win, stride, mi_bins=16)
    centers = (np.arange(M.shape[0])*stride + win/2)/FS + EDGE_S
    ev = np.abs(centers - mt) <= 0.5
    no = ~ev
    out = {"merger_t": mt}
    for lbl, mask in [("event", ev), ("noise", no)]:
        Ms = M[mask]
        v, *_ = extract_v1(Ms); _, cosa = project_234(v)
        out[lbl] = dict(
            n=int(mask.sum()),
            H_a=float(Ms[:, 0].mean()), H_b=float(Ms[:, 1].mean()),
            asym_absHaHb=float(np.abs(Ms[:, 0]-Ms[:, 1]).mean()),
            MI=float(Ms[:, 5].mean()),
            cos_attr=float(abs(cosa)))
    ipk = int(np.argmax(M[:, 5]))
    out["peak_MI_t"] = float(centers[ipk]); out["peak_MI"] = float(M[ipk, 5])
    out["MI_event_over_noise"] = out["event"]["MI"]/out["noise"]["MI"]
    return out


results = {}
print("="*108)
print("PER-EVENT LIGO entropy readout (each event SEPARATE -- no pooling, no average)")
print("="*108)
hdr = f"{'event':10s} {'align':>9s} | {'win':>5s} {'H_a':>7s} {'H_b':>7s} {'|Ha-Hb|':>8s} {'MI':>7s} {'cosAttr':>8s}"
for name, hp, lp, mt in EVENTS:
    h1p = preprocess(load(hp)); l1p = preprocess(load(lp))
    l1a, lag, sign = align(h1p, l1p, mt)
    r = per_event(name, h1p, l1a, mt)
    r["align"] = dict(lag=lag, sign=sign)
    results[name] = r
    print("\n"+"-"*108)
    print(f"{name}  (merger {mt:.1f}s, align lag {lag:+d} samp {lag/FS*1e3:+.1f}ms sign {sign:+.0f})")
    print(hdr)
    for w in ("event", "noise"):
        d = r[w]
        print(f"{'':10s} {'':>9s} | {w:>5s} {d['H_a']:7.3f} {d['H_b']:7.3f} "
              f"{d['asym_absHaHb']:8.4f} {d['MI']:7.4f} {d['cos_attr']:8.3f}")
    print(f"   peak-MI at t={r['peak_MI_t']:.2f}s (merger {mt:.1f}s)  "
          f"peak={r['peak_MI']:.4f}  event/noise MI={r['MI_event_over_noise']:.2f}x")

json.dump(results, open("s11_ligo_perevent_results.json", "w"), indent=2)
print("\n"+"="*108)
print("Wrote s11_ligo_perevent_results.json  (3 events; full 12 via bulk-fetch batch next session)")
