"""
Cross-event test of gravity's off-substrate flow axis, fixed alignment.

Earlier (grav_flow_reproduce) only GW150914 aligned; GW170104/GW151226 failed
because cross-correlation hunted over a wide noisy window and the event window
was set by peak-MI (circular). Fixes:
  - anchor alignment on the KNOWN merger GPS time (tight +/-0.12 s window),
    searching lag in +/-12 ms AND the L1 sign (antenna pattern), pick max |corr|.
  - define the EVENT window by the KNOWN merger time (+/-0.5 s), not peak-MI.
Then measure the flow axis (dMI/dt regressed on [H_a,H_b,H_a^2,H_b^2,H_a*H_b]) in
event vs noise vs full, per event, and ask: does the off-substrate flow axis
(cos ~0.67 on GW150914) REPRODUCE across events in the event window with real R^2?

Held-fixed prediction (Greg, locked before the data): structured/existing domains
cluster ON the flow-substrate (high cos); gravity sits FAR OFF (low cos). So
event-window cos(flow->substrate) should be LOW and similar across events, with
event R^2 >> noise R^2, and the event flow axes should point a CONSISTENT way.
Deflationary alternatives kept live: (1) detector-state fingerprint (off-axis also
in noise, INFO-038); (2) construction-type (real strain vs sim ensembles) -- not
decidable here, needs a non-gravity real detector-pair object (HBT, S15).
"""
import json
import time

import numpy as np

from kbk_pipeline import compute_operator_matrix, extract_v1
from grav_time_retaining import load_strain, preprocess, FS, EDGE_S, WIN_S, STRIDE_S

SUBSTRATE_Q = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)

# (H1, L1, label, known merger time into 32s segment)
EVENTS = [
    ("data/ligo/H-H1_GW150914_32s.hdf5", "data/ligo/L-L1_GW150914_32s.hdf5", "GW150914", 16.4),
    ("data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5", "GW170104", 16.6),
    ("data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5", "GW151226", 16.6),
]


def cos(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else float("nan")


def align(h1p, l1p, merger_t, fs=FS, search_s=0.012, half=0.12):
    """Find (lag, sign) maximizing |corr| in a tight window around the known
    merger. Returns aligned L1 = sign*roll(l1p, lag)."""
    lo, hi = int((merger_t - half) * fs), int((merger_t + half) * fs)
    a = h1p[lo:hi]
    best, bl, bs = -np.inf, 0, 1.0
    for sign in (1.0, -1.0):
        ls = sign * l1p
        for lag in range(-int(search_s * fs), int(search_s * fs) + 1):
            c = abs(np.corrcoef(a, np.roll(ls, lag)[lo:hi])[0, 1])
            if c > best:
                best, bl, bs = c, lag, sign
    return bs * np.roll(l1p, bl), int(bl), float(bs), float(best)


def flow(MI, X, dt_row, mask):
    dMI = np.gradient(MI) / dt_row
    Xb, yb = X[mask], dMI[mask]
    if Xb.shape[0] < Xb.shape[1] + 2:
        return None
    beta, *_ = np.linalg.lstsq(Xb, yb, rcond=None)
    pred = Xb @ beta
    r2 = 1.0 - np.sum((yb - pred) ** 2) / (np.sum((yb - yb.mean()) ** 2) + 1e-18)
    quad = np.array([beta[2], beta[3], beta[4]])
    return {"n": int(Xb.shape[0]), "R2": round(float(r2), 4),
            "cos_to_substrate": round(abs(cos(quad, SUBSTRATE_Q)), 3),
            "quad_axis": [round(float(x), 4) for x in quad / (np.linalg.norm(quad) + 1e-18)],
            "_q": quad}


def run(h1f, l1f, label, mt):
    h1 = preprocess(load_strain(h1f)); l1 = preprocess(load_strain(l1f))
    l1a, lag, sign, corr = align(h1, l1, mt)
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1a[crop:-crop],
                                int(WIN_S * FS), int(STRIDE_S * FS), mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * int(STRIDE_S * FS) + int(WIN_S * FS) / 2) / FS + EDGE_S
    MI = M[:, 5]
    X = np.column_stack([M[:, 0], M[:, 1], M[:, 2], M[:, 3], M[:, 4], np.ones(n)])
    ev = np.abs(centers - mt) <= 0.5
    noise = np.abs(centers - mt) > 2.0
    blocks = {"full": flow(MI, X, STRIDE_S, np.ones(n, bool)),
              "event": flow(MI, X, STRIDE_S, ev),
              "noise": flow(MI, X, STRIDE_S, noise)}
    # MI excess at merger as a detection sanity check
    mi_excess = float(MI[ev].mean() / (np.median(MI[noise]) + 1e-12))
    return {"label": label, "merger_t": mt, "align_lag_samp": lag,
            "align_sign": sign, "align_corr": round(corr, 3),
            "MI_event_over_noise": round(mi_excess, 3), "blocks": blocks}


def main():
    t0 = time.time()
    print("[grav_flow_crossevent] fixed-alignment cross-event flow-axis test\n", flush=True)
    out = [run(*e) for e in EVENTS]
    for r in out:
        print(f"=== {r['label']}  align(lag={r['align_lag_samp']},sign={r['align_sign']:+.0f},"
              f"corr={r['align_corr']})  MI_event/noise={r['MI_event_over_noise']} ===", flush=True)
        for bn in ("full", "event", "noise"):
            b = r["blocks"][bn]
            if b:
                print(f"   {bn:6s}: n={b['n']:4d} R2={b['R2']:+.4f} "
                      f"cos(flow->substrate)={b['cos_to_substrate']:.3f} axis={b['quad_axis']}",
                      flush=True)
        print(flush=True)
    print("cross-event EVENT-window |cos| of flow axes (does direction reproduce?):", flush=True)
    ev = {r["label"]: np.array(r["blocks"]["event"]["_q"]) for r in out if r["blocks"]["event"]}
    L = list(ev); cross = {}
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            c = abs(cos(ev[L[i]], ev[L[j]])); cross[f"{L[i]}__{L[j]}"] = round(c, 3)
            print(f"   |cos|({L[i]},{L[j]}) = {c:.3f}", flush=True)
    print("\nevent-window cos(flow->substrate) summary:", flush=True)
    for r in out:
        b = r["blocks"]["event"]
        if b:
            print(f"   {r['label']}: cos->substrate={b['cos_to_substrate']}  "
                  f"R2={b['R2']}  (MI_event/noise={r['MI_event_over_noise']})", flush=True)
    for r in out:
        for bn in r["blocks"]:
            if r["blocks"][bn]:
                r["blocks"][bn].pop("_q", None)
    json.dump({"events": out, "cross_event_cos": cross,
               "elapsed_s": round(time.time() - t0, 1)},
              open("grav_flow_crossevent_results.json", "w"), indent=2)
    print(f"\n[grav_flow_crossevent] wrote results ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
