"""
S14 robustness follow-up to s14_force_strong (NO tent-widening on the q-range).

Two checks the main run could not make:
 (1) IS the Bose-Einstein (strong-interaction femtoscopy) correlation actually
     PRESENT in this data? -- the INFO-049 lesson: confirm the signal is real
     before reading the null. Build C(q) = N_same-event(q) / N_mixed-event(q) for
     same-charge identical-hadron pairs (mixed-event = tracks paired across
     different events, destroying the BE correlation). BE -> C(q) > 1 at low q.
 (2) The main run's equal-count q-bins bottomed out at q ~ 0.26 GeV; the BE peak
     lives at q -> 0. Re-run the 6-op extraction RESTRICTED to the low-q
     femtoscopy region (q < Q_LOW) with fine equal-count bins, so the operator
     object actually samples where the correlation is. Does MI enter the null
     THERE? (self-pole reading is only solid if it survives at low q.)

Speaking posture: I expect (1) C(q) rises at low q (BE is well established in
Pb-Pb) and (2) MI still stays out of the null even at low q (self-pole). But if
MI enters the null at low q where the genuine correlation peaks, that is the
interesting outcome. No verdict in advance.

Outputs s14_strong_lowq_becheck_results.json.
"""
import json, time
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d
from s14_force_strong import (read_tracks, decompose, OPS, N_BINS,
                              d_lin, d_quad, d_mi)

Q_LOW = 0.4          # GeV; low-q femtoscopy window for the focused extraction
N_BINS_LOW = 20
MIN_PER_BIN = 200


def same_charge_q(events, cs, observable, rng, qmax):
    """same-event same-charge pairs -> (obsA,obsB,q) with random A/B labels."""
    obsA, obsB, qall = [], [], []
    for ev in events:
        m = ev["charge"] == cs
        if m.sum() < 2:
            continue
        E, px, py, pz = ev["E"][m], ev["px"][m], ev["py"][m], ev["pz"][m]
        ob = ev[observable][m]
        ii, jj = np.triu_indices(len(E), k=1)
        q2 = ((px[ii]-px[jj])**2 + (py[ii]-py[jj])**2
              + (pz[ii]-pz[jj])**2 - (E[ii]-E[jj])**2)
        keep = (q2 > 0) & (q2 < qmax*qmax)
        q = np.sqrt(q2[keep]); oi, oj = ob[ii][keep], ob[jj][keep]
        sw = rng.random(len(q)) < 0.5
        obsA.append(np.where(sw, oj, oi)); obsB.append(np.where(sw, oi, oj)); qall.append(q)
    if not obsA:
        return (np.array([]),)*3
    return np.concatenate(obsA), np.concatenate(obsB), np.concatenate(qall)


def mixed_event_q(events, cs, rng, qmax, ntarget):
    """mixed-event same-charge pairs: track from event e vs track from e'!=e.
    Destroys BE correlation -> femtoscopy denominator."""
    pool = []  # (E,px,py,pz, evt_id)
    for k, ev in enumerate(events):
        m = ev["charge"] == cs
        for E, px, py, pz in zip(ev["E"][m], ev["px"][m], ev["py"][m], ev["pz"][m]):
            pool.append((E, px, py, pz, k))
    pool = np.array(pool)
    if len(pool) < 2:
        return np.array([])
    qs = []
    tries = 0
    while len(qs) < ntarget and tries < ntarget*4:
        i, j = rng.integers(0, len(pool), 2)
        tries += 1
        if pool[i, 4] == pool[j, 4]:
            continue  # same event -> skip (we want mixed)
        E1, px1, py1, pz1, _ = pool[i]; E2, px2, py2, pz2, _ = pool[j]
        q2 = (px1-px2)**2 + (py1-py2)**2 + (pz1-pz2)**2 - (E1-E2)**2
        if 0 < q2 < qmax*qmax:
            qs.append(np.sqrt(q2))
    return np.array(qs)


def operator_lowq(obsA, obsB, q, qlow, n_bins):
    sel = q < qlow
    a, b, qq = obsA[sel], obsB[sel], q[sel]
    order = np.argsort(qq); a, b, qq = a[order], b[order], qq[order]
    edges = np.linspace(0, len(qq), n_bins+1).astype(int)
    rows, centers, ns = [], [], []
    for i in range(n_bins):
        s, e = edges[i], edges[i+1]
        if e - s < MIN_PER_BIN:
            continue
        Ha = entropy_hist_1d(a[s:e], bins=20); Hb = entropy_hist_1d(b[s:e], bins=20)
        MIv = mi_hist_2d(a[s:e], b[s:e], bins=14)
        rows.append([Ha, Hb, Ha*Ha, Hb*Hb, Ha*Hb, MIv])
        centers.append(float(qq[s:e].mean())); ns.append(int(e-s))
    return np.array(rows), np.array(centers), ns


def main():
    t0 = time.time()
    import glob
    files = sorted(glob.glob("data/strong/*.root"))
    rng0 = np.random.default_rng(0)
    events = []
    for f in files:
        events.extend(read_tracks(f, 300, rng0))
    print(f"loaded {len(events)} events from {len(files)} files ({time.time()-t0:.1f}s)")

    out = dict(meta=dict(n_events=len(events), q_low=Q_LOW))

    # ---- (1) BE presence: C(q) = same-event / mixed-event, same-charge ----
    rng = np.random.default_rng(7)
    a, b, q_same = same_charge_q(events, +1, "pt", rng, qmax=0.5)
    q_mix = mixed_event_q(events, +1, rng, qmax=0.5, ntarget=len(q_same))
    qbins = np.linspace(0, 0.5, 26)
    hs, _ = np.histogram(q_same, bins=qbins)
    hm, _ = np.histogram(q_mix, bins=qbins)
    cen = 0.5*(qbins[1:]+qbins[:-1])
    # normalize C(q) to 1 in the high-q plateau (0.35-0.5)
    norm = (hs[cen > 0.35].sum() / max(hm[cen > 0.35].sum(), 1))
    Cq = (hs / np.maximum(hm, 1)) / norm
    lowq = Cq[cen < 0.1]
    print("\n(1) BE presence  C(q)=same/mixed (same-charge, normalized at q>0.35):")
    for c, v in zip(cen[:8], Cq[:8]):
        print(f"    q={c:.3f}  C={v:.3f}")
    print(f"    --> C(q<0.1) mean = {lowq.mean():.3f}  "
          f"({'BE enhancement PRESENT' if lowq.mean() > 1.1 else 'no clear enhancement'})")
    out["BE_presence"] = dict(q_centers=cen.tolist(), Cq=Cq.tolist(),
                              Cq_lowq_mean=float(lowq.mean()),
                              n_same=int(len(q_same)), n_mixed=int(len(q_mix)))

    # ---- (2) low-q focused operator extraction (does MI enter null there?) ----
    print(f"\n(2) low-q operator extraction (q < {Q_LOW} GeV, fine bins):")
    cfgs = []
    for cs in (+1, -1):
        for obs in ("pt", "eta"):
            mis, coss, eqs, asyms = [], [], [], []
            for sd in range(3):
                r = np.random.default_rng(200+sd)
                A, B, Q = same_charge_q(events, cs, obs, r, qmax=Q_LOW+0.05)
                Mop, cen2, ns = operator_lowq(A, B, Q, Q_LOW, N_BINS_LOW)
                if Mop.shape[0] < 6:
                    continue
                v, S, _, _ = extract_v1(Mop)
                _, cosa = project_234(v)
                eq, mi, res = decompose(v)
                mis.append(mi); coss.append(abs(cosa)); eqs.append(eq)
                asyms.append(float(np.abs(Mop[:, 0]-Mop[:, 1]).mean()))
            if not mis:
                print(f"    charge{cs:+d} {obs}: insufficient low-q pairs"); continue
            mis, coss, eqs, asyms = map(np.array, (mis, coss, eqs, asyms))
            verdict = "MI IN null (COUPLING!)" if mis.mean() > 0.5 else "MI NOT in null (self-pole)"
            print(f"    charge{cs:+d} {obs} (q in [{cen2.min():.3f},{cen2.max():.3f}], "
                  f"{int(np.mean(ns))}/bin): |Ha-Hb| {asyms.mean():.3f}  "
                  f"eqEnt {eqs.mean():.3f}  MI-coupling {mis.mean():.3f}+/-{mis.std():.3f}  "
                  f"--> {verdict}")
            cfgs.append(dict(config=f"charge{cs:+d}_{obs}",
                             MI_coupling_mean=float(mis.mean()),
                             MI_coupling_std=float(mis.std()),
                             equal_entropy_mean=float(eqs.mean()),
                             cos_attractor_mean=float(coss.mean()),
                             mean_abs_HaHb=float(asyms.mean()),
                             q_lo=float(cen2.min()), q_hi=float(cen2.max())))
    out["lowq_extraction"] = cfgs
    out["meta"]["elapsed_s"] = round(time.time()-t0, 2)
    json.dump(out, open("s14_strong_lowq_becheck_results.json", "w"), indent=2)
    print(f"\nWrote s14_strong_lowq_becheck_results.json ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
