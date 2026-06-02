"""
S13 #3 -- LIGO no-MI-basis decomposition across the FULL 12 events.

S11 (s11_first_run_entropy.py, GW150914 only) found a REVERSAL: in the
full 6-op basis NOISE windows sit ON the (-1,-1,+2) equal-entropy attractor
(cos 0.980) and the EVENT (merger) windows LEAVE it (0.231) -- the merger
departure is MI-driven. But in the no-MI 5-op basis [H_a,H_b,H_a^2,H_b^2,
H_a*H_b] the picture FLIPPED: NOISE sat OFF the attractor (0.301) and EVENT
sat ON it (0.959). That reversal was only checked on one event.

Kickoff #3: re-run the no-MI decomposition per event across all 12. Does the
noise-OFF / event-ON reversal hold or flip across events?

Speaking posture (Rule C): I THINK the 6-op picture (noise-ON / event-OFF,
MI-driven merger departure) will hold broadly because it is the INFO-036
equal-entropy geometry; the no-MI REVERSAL I am much less sure about -- on
one event it could be a quirk of how MI removal reshuffles the smallest
singular vector. We wait on the 12-event spread; no verdict in advance, and
we report it per event (no pooling), loud and quiet alike.

Reuses kbk_pipeline pure funcs only (self-contained helpers copied to avoid
the s11_* modules' import-time side effects). Each event uses ONE 32s segment
around the merger (event windows |t-mt|<=0.5, noise = the rest of the clip);
no off-source null needed for this decomposition. Bulk files are deleted
after each event to keep peak disk low. Incremental save + resume-safe.

Run:  python3 s13_ligo_nomi_batch.py            # all 12
      python3 s13_ligo_nomi_batch.py --canary   # first 2
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


# ---- fetch helpers (copied from s11_ligo_batch to avoid import side effects) --
def event_meta(ev):
    d = json.load(urllib.request.urlopen(
        f"https://gwosc.org/eventapi/json/event/{ev}/", timeout=60))
    k = list(d["events"])[0]; e = d["events"][k]
    gps = e["GPS"]; url = {}
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


# ---- preprocess + align (copied from s11_ligo_perevent) ----------------------
def bandpass(x, lo=35., hi=350., order=4):
    b, a = butter(order, [lo/(FS/2), hi/(FS/2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x):
    n = len(x); fr, ps = welch(x, fs=FS, nperseg=min(4*FS, n))
    ip = interp1d(fr, ps, bounds_error=False, fill_value=(ps[0], ps[-1]))
    X = np.fft.rfft(x); f = np.fft.rfftfreq(n, 1./FS)
    xw = np.fft.irfft(X/np.sqrt(np.maximum(ip(f), 1e-50)), n=n)
    return xw/np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


def align(h, l, mt, max_lag_s=0.015):
    c = int(mt*FS); half = int(0.25*FS)
    xc = correlate(h[c-half:c+half], l[c-half:c+half], mode="full")
    lags = np.arange(-half*2+1, half*2); m = np.abs(lags) <= int(max_lag_s*FS)
    j = np.argmax(np.abs(xc[m]))
    return np.roll(l, lags[m][j]) * np.sign(xc[m][j]), int(lags[m][j])


# ---- the two cos readouts (copied from s11_first_run_entropy) -----------------
def cos_full(M):
    v, *_ = extract_v1(M); _, c = project_234(v); return abs(float(c))


def cos_noMI(M):
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
    crop = int(EDGE_S*FS)
    a, b = h[crop:-crop], la[crop:-crop]
    win, stride = int(0.125*FS), int(0.03125*FS)
    M = compute_operator_matrix(a, b, win, stride, mi_bins=16)
    centers = (np.arange(M.shape[0])*stride + win/2)/FS + EDGE_S
    evm = np.abs(centers - mt) <= 0.5; nom = ~evm
    out = {"GPS": gps, "align_lag": lag, "merger_t": mt}
    for lbl, mask in [("event", evm), ("noise", nom)]:
        Ms = M[mask]
        out[lbl] = dict(
            n=int(mask.sum()),
            H_a=float(Ms[:, 0].mean()), H_b=float(Ms[:, 1].mean()),
            asym=float(np.abs(Ms[:, 0]-Ms[:, 1]).mean()),
            MI=float(Ms[:, 5].mean()),
            cos6=cos_full(Ms), cosNoMI=cos_noMI(Ms))
    # the reversal test, per event
    out["reversal_6op_noiseON_eventOFF"] = bool(
        out["noise"]["cos6"] - out["event"]["cos6"] > 0.2)
    out["reversal_noMI_noiseOFF_eventON"] = bool(
        out["event"]["cosNoMI"] - out["noise"]["cosNoMI"] > 0.2)
    # cleanup bulk to keep disk low
    for p in paths.values():
        try: os.remove(p)
        except OSError: pass
    return out


t0 = time.time()
evs = EVENTS[:2] if CANARY else EVENTS
fn = "s13_ligo_nomi_canary.json" if CANARY else "s13_ligo_nomi_results.json"
print("=" * 104)
print(f"S13 #3 LIGO no-MI decomposition {'[CANARY]' if CANARY else ''}: "
      f"per-event 6-op vs 5-op(no MI) attractor cos (no pooling)")
print("=" * 104)
print(f"{'event':10s} | {'block':6s} {'|Ha-Hb|':>8s} {'MI':>7s} "
      f"{'cos6':>7s} {'cosNoMI':>8s}")


def save(results):
    json.dump(dict(meta=dict(canary=CANARY, elapsed_s=round(time.time()-t0, 1),
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
        for w in ("noise", "event"):
            d = r[w]
            print(f"{ev if w=='noise' else '':10s} | {w:6s} {d['asym']:8.4f} "
                  f"{d['MI']:7.4f} {d['cos6']:7.3f} {d['cosNoMI']:8.3f}")
        print(f"{'':10s} | 6op noise-ON/event-OFF: "
              f"{r['reversal_6op_noiseON_eventOFF']}   "
              f"noMI noise-OFF/event-ON: {r['reversal_noMI_noiseOFF_eventON']}")
    except Exception as ex:
        results[ev] = {"error": str(ex)}; print(f"\n{ev}: ERROR {ex}")
    save(results); gc.collect()

# summary across events
ok = {k: v for k, v in results.items() if "error" not in v and "skipped" not in v}
n6 = sum(v["reversal_6op_noiseON_eventOFF"] for v in ok.values())
nm = sum(v["reversal_noMI_noiseOFF_eventON"] for v in ok.values())
print("\n" + "=" * 104)
print(f"SUMMARY ({len(ok)} events scored):")
print(f"  6-op noise-ON / event-OFF (INFO-036 geometry): {n6}/{len(ok)}")
print(f"  no-MI noise-OFF / event-ON REVERSAL (GW150914 thread): {nm}/{len(ok)}")
save(results)
print(f"\nWrote {fn} ({time.time()-t0:.0f}s)")
