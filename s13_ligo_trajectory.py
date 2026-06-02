"""
S13 #3b -- "does anything right BEFORE or AFTER the event tell us something?"
(Greg's question). The binary event/noise split (|t-mt|<=0.5 vs rest) throws
away the approach (inspiral) and the aftermath (ringdown). Here we resolve
TIME relative to the merger into concentric zones and read MI / entropy-
asymmetry / attractor-cos in each, per event, no pooling.

Zones (tau = window_center - merger_t):
  merger    |tau| <= 0.125
  pre_near  -0.50 <= tau < -0.125     (late inspiral, just before)
  pre_far   -2.00 <= tau < -0.50      (early inspiral / approach)
  post_near  0.125 < tau <= 0.50      (ringdown / just after)
  post_far   0.50  < tau <= 2.00      (settling / just after, later)
  far        |tau| > 3.0              (clean noise reference)

Readout per zone: mean MI, mean |H_a-H_b|, and (where >=8 windows) the
attractor cos in the full 6-op and no-MI 5-op bases. The before/after question:
does MI/structure RAMP in pre_near vs pre_far/far (a precursor riding the
inspiral), and does it persist in post_near (aftermath), and is the profile
TIME-ASYMMETRIC around the merger (gradual chirp before, sharp ringdown after)?

Speaking posture (Rule C): I THINK loud events will show MI rising through
pre_far -> pre_near -> merger (the chirp building) and dropping quickly after
(ringdown brief), so BEFORE carries more than AFTER; quiet O1 events may be
indistinguishable from far noise. We wait on the per-event spread; no verdict
in advance, loud and quiet reported alike, no pooling.

Caches the per-event operator matrix M + window centers to data/ligo_M/{ev}.npz
so future per-window analyses need no re-download. Bulk strain deleted after.

Run:  python3 s13_ligo_trajectory.py [--canary]
"""
import sys, os, json, time, gc, urllib.request
import numpy as np, h5py
from scipy.signal import butter, filtfilt, welch, correlate
from scipy.interpolate import interp1d
from kbk_pipeline import compute_operator_matrix, extract_v1, project_234

CANARY = "--canary" in sys.argv
FS = 4096; EDGE_S = 2.0; SEG = 32
ATTR3 = np.array([-1., -1., 2.]) / np.sqrt(6.)
EVENTS = ["GW150914", "GW151012", "GW151226", "GW170104", "GW170608",
          "GW170729", "GW170809", "GW170814", "GW170817", "GW170818",
          "GW170823", "GW190521"]
ZONES = [("pre_far", -2.0, -0.5), ("pre_near", -0.5, -0.125),
         ("merger", -0.125, 0.125), ("post_near", 0.125, 0.5),
         ("post_far", 0.5, 2.0)]


def event_meta(ev):
    d = json.load(urllib.request.urlopen(
        f"https://gwosc.org/eventapi/json/event/{ev}/", timeout=60))
    k = list(d["events"])[0]; e = d["events"][k]; gps = e["GPS"]; url = {}
    for s in e["strain"]:
        if (s["sampling_rate"] == 4096 and s["duration"] == 4096
                and s["format"] == "hdf5" and s["detector"] in ("H1", "L1")):
            url.setdefault(s["detector"], s)
    return gps, {det: (s["url"], s["GPSstart"]) for det, s in url.items()}


def fetch_bulk(ev, det, url):
    os.makedirs("data/ligo_bulk", exist_ok=True)
    fn = f"data/ligo_bulk/{ev}_{det}.hdf5"
    if not os.path.exists(fn) or os.path.getsize(fn) < 1e6:
        urllib.request.urlretrieve(url, fn)
    return fn


def seg_strain(path, gps_start_file, gps_center):
    with h5py.File(path, "r") as f:
        full = f["strain/Strain"]
        i0 = int((gps_center - SEG/2 - gps_start_file) * FS)
        return full[i0:i0 + SEG*FS].astype(float)


def bandpass(x, lo=35., hi=350., order=4):
    b, a = butter(order, [lo/(FS/2), hi/(FS/2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x):
    n = len(x); fr, ps = welch(x, fs=FS, nperseg=min(4*FS, n))
    ip = interp1d(fr, ps, bounds_error=False, fill_value=(ps[0], ps[-1]))
    X = np.fft.rfft(x); f = np.fft.rfftfreq(n, 1./FS)
    return (np.fft.irfft(X/np.sqrt(np.maximum(ip(f), 1e-50)), n=n)) / np.std(
        np.fft.irfft(X/np.sqrt(np.maximum(ip(f), 1e-50)), n=n))


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def align(h, l, mt, max_lag_s=0.015):
    c = int(mt*FS); half = int(0.25*FS)
    xc = correlate(h[c-half:c+half], l[c-half:c+half], mode="full")
    lags = np.arange(-half*2+1, half*2); m = np.abs(lags) <= int(max_lag_s*FS)
    j = np.argmax(np.abs(xc[m]))
    return np.roll(l, lags[m][j]) * np.sign(xc[m][j]), int(lags[m][j])


def cos_full(M):
    if M.shape[0] < 8:
        return None
    v, *_ = extract_v1(M); _, c = project_234(v); return abs(float(c))


def cos_noMI(M):
    if M.shape[0] < 8:
        return None
    M5c = M[:, :5] - M[:, :5].mean(0)
    _, _, Vt = np.linalg.svd(M5c, full_matrices=False)
    v5 = Vt[-1]; sub = v5[2:5]; n = np.linalg.norm(sub)
    return abs(float(sub @ ATTR3 / n)) if n > 0 else 0.0


def run_event(ev):
    gps, urls = event_meta(ev)
    if "H1" not in urls or "L1" not in urls:
        return {"skipped": f"missing detector(s): {list(urls)}"}
    paths = {d: fetch_bulk(ev, d, urls[d][0]) for d in ("H1", "L1")}
    starts = {d: urls[d][1] for d in ("H1", "L1")}
    mt = SEG/2
    h = preprocess(seg_strain(paths["H1"], starts["H1"], gps))
    l = preprocess(seg_strain(paths["L1"], starts["L1"], gps))
    la, lag = align(h, l, mt)
    crop = int(EDGE_S*FS); a, b = h[crop:-crop], la[crop:-crop]
    win, stride = int(0.125*FS), int(0.03125*FS)
    M = compute_operator_matrix(a, b, win, stride, mi_bins=16)
    centers = (np.arange(M.shape[0])*stride + win/2)/FS + EDGE_S
    tau = centers - mt
    os.makedirs("data/ligo_M", exist_ok=True)
    np.savez_compressed(f"data/ligo_M/{ev}.npz", M=M, tau=tau, mt=mt, lag=lag)

    far = M[np.abs(tau) > 3.0]
    far_MI = float(far[:, 5].mean()); far_asym = float(np.abs(far[:, 0]-far[:, 1]).mean())
    out = {"GPS": gps, "align_lag": lag, "merger_t": mt,
           "far_MI": far_MI, "far_asym": far_asym,
           "far_cos6": cos_full(far), "far_cosNoMI": cos_noMI(far), "zones": {}}
    for zname, lo, hi in ZONES:
        msk = (tau >= lo) & (tau < hi)
        Ms = M[msk]
        if Ms.shape[0] == 0:
            out["zones"][zname] = {"n": 0}; continue
        out["zones"][zname] = dict(
            n=int(msk.sum()),
            MI=float(Ms[:, 5].mean()),
            MI_over_far=float(Ms[:, 5].mean()/far_MI),
            asym=float(np.abs(Ms[:, 0]-Ms[:, 1]).mean()),
            cos6=cos_full(Ms), cosNoMI=cos_noMI(Ms))
    # peak-MI window time (precursor vs merger-centered?)
    ipk = int(np.argmax(M[:, 5]))
    out["peak_MI_tau"] = float(tau[ipk]); out["peak_MI"] = float(M[ipk, 5])
    for p in paths.values():
        try: os.remove(p)
        except OSError: pass
    return out


t0 = time.time()
evs = EVENTS[:2] if CANARY else EVENTS
fn = "s13_ligo_trajectory_canary.json" if CANARY else "s13_ligo_trajectory_results.json"
print("=" * 110)
print(f"S13 #3b LIGO peri-event trajectory {'[CANARY]' if CANARY else ''}: "
      f"MI / |Ha-Hb| / attractor-cos in zones around the merger (no pooling)")
print("=" * 110)


def save(results):
    json.dump(dict(meta=dict(canary=CANARY, zones=ZONES,
                             elapsed_s=round(time.time()-t0, 1),
                             n_done=len(results)), events=results),
              open(fn, "w"), indent=2)


results = {}
if os.path.exists(fn):
    try:
        prev = json.load(open(fn)).get("events", {})
        results = {k: v for k, v in prev.items()
                   if isinstance(v, dict) and "error" not in v}
        if results:
            print(f"[resume] keeping {len(results)} scored: {list(results)}")
    except Exception:
        results = {}

for ev in evs:
    if ev in results:
        continue
    try:
        r = run_event(ev)
        results[ev] = r
        if "skipped" in r:
            print(f"\n{ev}: SKIPPED ({r['skipped']})"); save(results); continue
        print(f"\n{ev} (peak-MI at tau={r['peak_MI_tau']:+.3f}s, "
              f"far MI={r['far_MI']:.3f} asym={r['far_asym']:.3f}):")
        print(f"  {'zone':9s} {'n':>3s} {'MI':>7s} {'MI/far':>7s} "
              f"{'|Ha-Hb|':>8s} {'cos6':>6s} {'cosNoMI':>8s}")
        for zname, _, _ in ZONES:
            z = r["zones"][zname]
            if z["n"] == 0:
                continue
            c6 = f"{z['cos6']:.3f}" if z['cos6'] is not None else "  -"
            cn = f"{z['cosNoMI']:.3f}" if z['cosNoMI'] is not None else "  -"
            print(f"  {zname:9s} {z['n']:3d} {z['MI']:7.4f} {z['MI_over_far']:7.3f} "
                  f"{z['asym']:8.4f} {c6:>6s} {cn:>8s}")
    except Exception as ex:
        results[ev] = {"error": str(ex)}; print(f"\n{ev}: ERROR {ex}")
    save(results); gc.collect()

save(results)
print(f"\nWrote {fn} ({time.time()-t0:.0f}s)")
