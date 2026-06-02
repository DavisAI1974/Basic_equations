"""
Test whether gravity's OFF-substrate flow axis (cos 0.608, the flow_dipole_axis
whole-segment number) REPRODUCES. Greg's rule: never discard data -- the weak
whole-segment R^2 means "measure it where it has signal", not "throw it away".

The whole-segment regression of dMI/dt was ~99% noise (R^2 0.024). Here we:
  1. Align each event CORRECTLY (cross-correlate whitened H1 vs inverted L1 in the
     central region to find the real inter-detector lag, not assume it).
  2. Measure the flow axis (dMI/dt regressed on [H_a,H_b,H_a^2,H_b^2,H_a*H_b]) in
     THREE blocks: FULL segment, EVENT window (+/-0.5 s of peak-MI, where dMI/dt
     has real signal), NOISE-only (outside +/-2 s).
  3. Report per block: R^2, cos(flow-quad -> equal-entropy substrate), the axis.
  4. Cross-event: does the EVENT-window flow axis point the SAME way across events?
     And a bootstrap on GW150914's event-window axis: is 0.608-ish stable or does
     it swing (i.e. is the number determined by the data or not)?

Readings to distinguish (no pre-assigned meaning):
  - off-substrate axis reproduces across events IN THE EVENT WINDOW with real R^2
    -> a genuine gravity flow signature worth chasing.
  - appears in NOISE too / event-specific direction -> detector-state fingerprint
    (INFO-038/046), a real WHY, not random, but not source physics.
  - swings under bootstrap -> the data does not determine the axis; 0.608 was
    undetermined (still logged, not discarded).
"""
import json
import time

import numpy as np

from kbk_pipeline import compute_operator_matrix, extract_v1
from grav_time_retaining import load_strain, preprocess, FS, EDGE_S, WIN_S, STRIDE_S

SUBSTRATE_Q = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)

EVENTS = [
    ("data/ligo/H-H1_GW150914_32s.hdf5", "data/ligo/L-L1_GW150914_32s.hdf5", "GW150914"),
    ("data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5", "GW170104"),
    ("data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
     "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5", "GW151226"),
]


def cos(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(np.dot(a, b) / (na * nb)) if na > 0 and nb > 0 else float("nan")


def best_lag(h1p, l1_inv, fs=FS, search_s=0.02, center_lo=14.0, center_hi=18.0):
    """Lag (samples) maximizing |corr(h1, shifted l1_inv)| in the central region
    where a loud common merger signal dominates."""
    lo, hi = int(center_lo * fs), int(center_hi * fs)
    a = h1p[lo:hi]
    best, bl = -np.inf, 0
    for lag in range(-int(search_s * fs), int(search_s * fs) + 1):
        b = np.roll(l1_inv, lag)[lo:hi]
        c = abs(np.corrcoef(a, b)[0, 1])
        if c > best:
            best, bl = c, lag
    return bl, float(best)


def flow_in_block(MI, X, dt_row, mask):
    """Regress dMI/dt (local derivative on the full series) on X within mask."""
    dMI = np.gradient(MI) / dt_row
    Xb, yb = X[mask], dMI[mask]
    if Xb.shape[0] < Xb.shape[1] + 2:
        return None
    beta, *_ = np.linalg.lstsq(Xb, yb, rcond=None)
    pred = Xb @ beta
    ss_res = float(np.sum((yb - pred) ** 2))
    ss_tot = float(np.sum((yb - yb.mean()) ** 2)) + 1e-18
    quad = np.array([beta[2], beta[3], beta[4]])
    return {
        "n": int(Xb.shape[0]),
        "R2": round(1.0 - ss_res / ss_tot, 4),
        "cos_to_substrate": round(abs(cos(quad, SUBSTRATE_Q)), 3),
        "quad_axis": [round(float(x), 4) for x in quad / (np.linalg.norm(quad) + 1e-18)],
        "_quad_raw": quad,
    }


def run_event(h1f, l1f, label):
    h1 = preprocess(load_strain(h1f))
    l1 = preprocess(load_strain(l1f))
    l1_inv = -l1
    lag, corr = best_lag(h1, l1_inv)
    l1a = np.roll(l1_inv, lag)
    crop = int(EDGE_S * FS)
    h1c, l1c = h1[crop:-crop], l1a[crop:-crop]
    M = compute_operator_matrix(h1c, l1c, int(WIN_S * FS), int(STRIDE_S * FS), mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * int(STRIDE_S * FS) + int(WIN_S * FS) / 2) / FS + EDGE_S
    MI = M[:, 5]
    X = np.column_stack([M[:, 0], M[:, 1], M[:, 2], M[:, 3], M[:, 4], np.ones(n)])

    ipk = int(np.argmax(MI))
    tpk = float(centers[ipk])
    ev = np.abs(centers - tpk) <= 0.5
    noise = np.abs(centers - tpk) > 2.0
    full = np.ones(n, bool)

    blocks = {"full": flow_in_block(MI, X, STRIDE_S, full),
              "event": flow_in_block(MI, X, STRIDE_S, ev),
              "noise": flow_in_block(MI, X, STRIDE_S, noise)}

    # bootstrap the EVENT-window axis (is it determined by the data?)
    dMI = np.gradient(MI) / STRIDE_S
    ev_idx = np.where(ev)[0]
    rng = np.random.default_rng(0)
    boot_cos = []
    for _ in range(200):
        s = rng.choice(ev_idx, size=len(ev_idx), replace=True)
        try:
            beta, *_ = np.linalg.lstsq(X[s], dMI[s], rcond=None)
            boot_cos.append(abs(cos(beta[2:5], SUBSTRATE_Q)))
        except Exception:
            pass
    boot = (float(np.median(boot_cos)), float(np.percentile(boot_cos, 16)),
            float(np.percentile(boot_cos, 84))) if boot_cos else (None, None, None)

    return {"label": label, "align_lag_samples": int(lag),
            "align_corr": round(corr, 3), "peak_MI_t_s": round(tpk, 3),
            "n_windows": n, "blocks": blocks,
            "event_axis_bootstrap_cos_to_substrate": {
                "median": round(boot[0], 3) if boot[0] is not None else None,
                "p16": round(boot[1], 3) if boot[1] is not None else None,
                "p84": round(boot[2], 3) if boot[2] is not None else None}}


def main():
    t0 = time.time()
    print("[grav_flow_reproduce] testing the 0.608 off-substrate flow axis\n", flush=True)
    out = []
    for h1f, l1f, label in EVENTS:
        r = run_event(h1f, l1f, label)
        out.append(r)
        print(f"=== {label}  (align lag={r['align_lag_samples']} samp, "
              f"corr={r['align_corr']}, peak-MI {r['peak_MI_t_s']}s) ===", flush=True)
        for bn in ("full", "event", "noise"):
            b = r["blocks"][bn]
            if b:
                print(f"   {bn:6s}: n={b['n']:4d}  R2={b['R2']:+.4f}  "
                      f"cos(flow->substrate)={b['cos_to_substrate']:.3f}  "
                      f"axis={b['quad_axis']}", flush=True)
        bs = r["event_axis_bootstrap_cos_to_substrate"]
        print(f"   event-axis bootstrap cos->substrate: "
              f"median={bs['median']} [p16 {bs['p16']}, p84 {bs['p84']}]\n", flush=True)

    # cross-event: do the EVENT-window flow axes point the same way?
    print("cross-event EVENT-window flow-axis cosines (does the direction reproduce?):",
          flush=True)
    evaxes = {r["label"]: np.array(r["blocks"]["event"]["_quad_raw"])
              for r in out if r["blocks"]["event"]}
    labels = list(evaxes)
    cross = {}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            c = abs(cos(evaxes[labels[i]], evaxes[labels[j]]))
            cross[f"{labels[i]}__{labels[j]}"] = round(c, 3)
            print(f"   |cos|({labels[i]}, {labels[j]}) = {c:.3f}", flush=True)

    # strip non-serializable raw vectors
    for r in out:
        for bn in r["blocks"]:
            if r["blocks"][bn] and "_quad_raw" in r["blocks"][bn]:
                r["blocks"][bn].pop("_quad_raw")
    json.dump({"events": out, "cross_event_event_axis_cos": cross,
               "elapsed_s": round(time.time() - t0, 1)},
              open("grav_flow_reproduce_results.json", "w"), indent=2)
    print(f"\n[grav_flow_reproduce] wrote results ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
