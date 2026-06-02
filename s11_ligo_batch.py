"""
S11 FULL-12 per-event LIGO batch (run in a fresh session; needs ~130MB
bulk download per detector per event = ~1.5GB for 12 events, 29G free).

Design (per Greg: run each event SEPARATELY, no merge / no average):
  For each event we fetch the 4096s/4096Hz bulk strain (H1,L1) ONCE. That
  single file gives BOTH:
    (a) the 32s EVENT segment around the merger, and
    (b) ~100 off-source 32s NOISE segments (away from merger) -> a per-event
        NULL DISTRIBUTION of peak-MI and event/noise-MI ratio, so the
        merger peak gets a p-value instead of a bare ratio. (This is the
        load-bearing piece the 3-event down payment lacked.)
  Per-event JSON only; nothing is pooled or averaged across events.

Reuses the SAME extraction as the toy systems + the 3-event readout
(kbk_pipeline.compute_operator_matrix / extract_v1 / project_234, and the
identical preprocess + generic cross-correlation align from
s11_ligo_perevent). Speaking posture (Rule C) carried from s11_ligo_perevent.

Run:  python3 s11_ligo_batch.py            # all 12, full null
      python3 s11_ligo_batch.py --canary   # first 2 events, small null
Cached downloads in data/ligo_bulk/ (gitignored-large; keep JSON only).
"""
import sys, os, json, time, urllib.request
import numpy as np, h5py
from kbk_pipeline import compute_operator_matrix, extract_v1, project_234
from s11_ligo_perevent import preprocess, align, FS, EDGE_S, per_event

CANARY = "--canary" in sys.argv
N_NULL = 8 if CANARY else 100          # off-source 32s segments per event
SEG = 32                                # segment length (s)

# 12 confident events spanning O1/O2/O3, loud + quiet + BNS + 3-detector.
EVENTS = ["GW150914", "GW151012", "GW151226", "GW170104", "GW170608",
          "GW170729", "GW170809", "GW170814", "GW170817", "GW170818",
          "GW170823", "GW190521"]


def event_meta(ev):
    d = json.load(urllib.request.urlopen(
        f"https://gwosc.org/eventapi/json/event/{ev}/", timeout=30))
    k = list(d["events"])[0]; e = d["events"][k]
    gps = e["GPS"]
    url = {}
    for s in e["strain"]:
        if (s["sampling_rate"] == 4096 and s["duration"] == 4096
                and s["format"] == "hdf5" and s["detector"] in ("H1", "L1")):
            url.setdefault(s["detector"], s)  # first match
    return gps, {det: (s["url"], s["GPSstart"]) for det, s in url.items()}


def fetch_bulk(ev, det, url):
    os.makedirs("data/ligo_bulk", exist_ok=True)
    fn = f"data/ligo_bulk/{ev}_{det}.hdf5"
    if not os.path.exists(fn) or os.path.getsize(fn) < 1e6:
        urllib.request.urlretrieve(url, fn)
    return fn


def seg_strain(path, gps_start_file, gps_center):
    """Pull a SEG-second slice centered on gps_center from a bulk file."""
    with h5py.File(path, "r") as f:
        full = f["strain/Strain"]
        i0 = int((gps_center - SEG/2 - gps_start_file) * FS)
        return full[i0:i0 + SEG*FS].astype(float)


def run_event(ev):
    gps, urls = event_meta(ev)
    if "H1" not in urls or "L1" not in urls:
        return {"skipped": f"missing detector(s): {list(urls)}"}
    paths = {d: fetch_bulk(ev, d, urls[d][0]) for d in ("H1", "L1")}
    starts = {d: urls[d][1] for d in ("H1", "L1")}
    mt = SEG/2  # merger centered in the 32s segment

    # --- event segment ---
    h = preprocess(seg_strain(paths["H1"], starts["H1"], gps))
    l = preprocess(seg_strain(paths["L1"], starts["L1"], gps))
    la, lag, sign = align(h, l, mt)
    ev_read = per_event(ev, h, la, mt)
    ev_read["align"] = dict(lag=lag, sign=sign)

    # --- per-event NULL: N off-source segments (>= 100s from merger) ---
    rng = np.random.default_rng(0)
    span = 4096
    null_peakMI, null_ratio = [], []
    for _ in range(N_NULL):
        off = gps + float(rng.choice([-1, 1])) * rng.uniform(150, span/2 - SEG)
        try:
            hh = preprocess(seg_strain(paths["H1"], starts["H1"], off))
            ll = preprocess(seg_strain(paths["L1"], starts["L1"], off))
            if len(hh) < SEG*FS - 10 or len(ll) < SEG*FS - 10:
                continue
            lla, _, _ = align(hh, ll, mt)
            r = per_event(ev, hh, lla, mt)  # "merger_t" is just the segment center
            null_peakMI.append(r["peak_MI"]); null_ratio.append(r["MI_event_over_noise"])
        except Exception:
            continue
    null_peakMI = np.array(null_peakMI); null_ratio = np.array(null_ratio)
    p_peak = float((null_peakMI >= ev_read["peak_MI"]).mean()) if len(null_peakMI) else None
    ev_read["null"] = dict(
        n=len(null_peakMI),
        peakMI_mean=float(null_peakMI.mean()) if len(null_peakMI) else None,
        peakMI_std=float(null_peakMI.std()) if len(null_peakMI) else None,
        p_peakMI=p_peak,
        ratio_mean=float(null_ratio.mean()) if len(null_ratio) else None)
    ev_read["GPS"] = gps
    return ev_read


t0 = time.time()
evs = EVENTS[:2] if CANARY else EVENTS
print("="*100)
print(f"S11 FULL-12 LIGO batch {'[CANARY]' if CANARY else ''}: per-event, no pooling, "
      f"with off-source null (N={N_NULL})")
print("="*100)
results = {}
for ev in evs:
    try:
        r = run_event(ev)
        results[ev] = r
        if "skipped" in r:
            print(f"\n{ev}: SKIPPED ({r['skipped']})"); continue
        e, n = r["event"], r["noise"]
        print(f"\n{ev} (GPS {r['GPS']:.1f}, align {r['align']['lag']:+d}samp):")
        print(f"   H_a={e['H_a']:.3f} H_b={e['H_b']:.3f} |Ha-Hb|={e['asym_absHaHb']:.4f} "
              f"cosAttr ev/noise={e['cos_attr']:.3f}/{n['cos_attr']:.3f}")
        print(f"   peak-MI t={r['peak_MI_t']:.2f}s (merger {r['merger_t']:.1f}s) "
              f"peak={r['peak_MI']:.4f}  null peakMI={r['null']['peakMI_mean']}"
              f"+/-{r['null']['peakMI_std']}  p={r['null']['p_peakMI']}")
    except Exception as ex:
        results[ev] = {"error": str(ex)}; print(f"\n{ev}: ERROR {ex}")

fn = "s11_ligo_batch_canary.json" if CANARY else "s11_ligo_batch_results.json"
json.dump(dict(meta=dict(canary=CANARY, n_null=N_NULL,
                         elapsed_s=round(time.time()-t0, 1)), events=results),
          open(fn, "w"), indent=2)
print(f"\nWrote {fn} ({time.time()-t0:.0f}s)")
