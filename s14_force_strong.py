"""
S14 -- STRONG force real-data operator object (the 4th and last within-type force).

Decision gate PASSED (this session): ATLAS DAOD_HION14 2015 Pb-Pb open data
(record 80036, /eos/opendata/atlas/rucio/data15_hi/, CC0) reads with pure-Python
uproot -- NO CMSSW/AliPhysics VM. (CMS HI RECO was rejected: uproot cannot
reconstruct the vector<reco::Track> member-wise counts; ATLAS DAOD aux-store
branches are flat and readable, per ATLAS's own usage note "can be used ...
using uproot".) The inter-hadron Bose-Einstein (two-pion femtoscopy) correlation
of identical same-charge hadrons is the intrinsic strong-interaction signal --
not invented.

Construction = the SAME family as WEAK (s13_force_dimuon), so the strong object
is commensurable with the weak one (within-type comparison; both are
particle-pair event ensembles binned by an energy axis). Per Greg's S13 rule:
NO synthesis across construction types -- strong is compared to WEAK only.

  - CHANNELS = the two identical same-charge hadrons of a pair, assigned to
    A / B AT RANDOM (symmetric labeling). Identical bosons carry no charge label
    to distinguish them, so a random A/B assignment avoids manufacturing an
    entropy asymmetry -- the strong analogue of weak's physical mu+/mu- charge
    symmetry. Observable per channel = hadron pT (eta as a robustness observable).
  - ENERGY AXIS = pair relative momentum q_inv = sqrt(-(p1-p2)^2), equal-count
    bins -> a trajectory. This plays the role dimuon mass M plays in the weak
    object and TIME plays in LIGO / per_domain_kbk ensemble-H. q is the
    femtoscopy scale: the Bose-Einstein enhancement lives at low q.
  - Per q-bin, over the pair sub-ensemble, compute the 6-op vector
    [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI] on (obsA,obsB). Stack -> operator matrix.
    extract_v1 (smallest null) -> project_234 -> cos to (-1,-1,+2)/sqrt6
    equal-entropy attractor + INFO-040 decomposition (equal-entropy / MI-coupling
    / residual) + MI-vs-H family probe.

Speaking posture (Rule C): I EXPECT strong to land on the equal-entropy SELF-POLE
(MI NOT in the null, like weak/EM/gravity-noise and physics/geology) -- which
would complete a clean 4/4 real-data statement that the gauge forces are
bookkeeping-side. BUT if MI ENTERS the null, that is the most interesting outcome
of the whole thread (the one real coupler so far is simulated biology). No
verdict in advance; wait on where the null sits.

Robustness (real-data analogue of the >=3-seed rule): agreement across
observables (pT, eta), charge sign (++ vs --), and random A/B-assignment seeds.

Usage:
  python s14_force_strong.py --canary          # 1 file, few events, quick
  python s14_force_strong.py --files f1 f2 ...  # full run on given local files
Outputs s14_force_strong_results.json (or _canary.json).
"""
import sys, json, time, glob, argparse
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d

M_PI = 0.13957039          # charged-pion mass (GeV); identical-pion femtoscopy
PT_LO, PT_HI = 0.2, 2.0    # GeV; standard soft-hadron femtoscopy window
ETA_MAX = 2.5
Q_MAX = 2.0                # GeV; keep the low-q femtoscopy region (bounds pairs)
N_BINS = 24
MIN_PER_BIN = 200
MAX_TRK = 150              # cap tracks/event (random subsample) to bound O(N^2)

d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]
PRE = "InDetTrackParticlesAuxDyn."
TREE = "CollectionTree"


def read_tracks(path, n_events, rng):
    """Yield per-event (pt, eta, phi, charge, E, px, py, pz) for selected tracks."""
    import uproot
    ct = uproot.open(path)[TREE]
    stop = ct.num_entries if n_events is None else min(n_events, ct.num_entries)
    a = ct.arrays([PRE + "phi", PRE + "theta", PRE + "qOverP", PRE + "HITight"],
                  entry_start=0, entry_stop=stop)
    out = []
    for i in range(stop):
        phi = np.asarray(a[PRE + "phi"][i], dtype=np.float64)
        th = np.asarray(a[PRE + "theta"][i], dtype=np.float64)
        q = np.asarray(a[PRE + "qOverP"][i], dtype=np.float64)
        hit = np.asarray(a[PRE + "HITight"][i])
        if len(q) == 0:
            continue
        good = (hit != 0) & (np.abs(q) > 0) & (th > 0) & (th < np.pi)
        phi, th, q = phi[good], th[good], q[good]
        if len(q) == 0:
            continue
        p = 1.0 / np.abs(q) / 1000.0          # GeV
        pt = np.sin(th) * p
        eta = -np.log(np.tan(th / 2))
        sel = (pt > PT_LO) & (pt < PT_HI) & (np.abs(eta) < ETA_MAX)
        phi, th, p, pt, eta = phi[sel], th[sel], p[sel], pt[sel], eta[sel]
        charge = np.sign(q[sel])
        if len(pt) > MAX_TRK:
            idx = rng.choice(len(pt), MAX_TRK, replace=False)
            phi, th, p, pt, eta, charge = (x[idx] for x in (phi, th, p, pt, eta, charge))
        px, py, pz = pt * np.cos(phi), pt * np.sin(phi), p * np.cos(th)
        E = np.sqrt(p * p + M_PI * M_PI)
        out.append(dict(pt=pt, eta=eta, charge=charge, E=E, px=px, py=py, pz=pz))
    return out


def build_pairs(events, charge_sign, observable, rng):
    """Same-charge pairs within each event; return (obsA, obsB, q_inv).
    A/B assignment randomized per pair (symmetric labeling)."""
    obsA, obsB, qall = [], [], []
    for ev in events:
        m = ev["charge"] == charge_sign
        if m.sum() < 2:
            continue
        E, px, py, pz = ev["E"][m], ev["px"][m], ev["py"][m], ev["pz"][m]
        ob = ev[observable][m]
        n = len(E)
        ii, jj = np.triu_indices(n, k=1)
        dE = E[ii] - E[jj]
        dpx, dpy, dpz = px[ii] - px[jj], py[ii] - py[jj], pz[ii] - pz[jj]
        q2 = dpx * dpx + dpy * dpy + dpz * dpz - dE * dE  # = -(p1-p2)^2
        keep = (q2 > 0) & (q2 < Q_MAX * Q_MAX)
        ii, jj = ii[keep], jj[keep]
        q = np.sqrt(q2[keep])
        oi, oj = ob[ii], ob[jj]
        # randomize channel assignment per pair
        swap = rng.random(len(q)) < 0.5
        a = np.where(swap, oj, oi)
        b = np.where(swap, oi, oj)
        obsA.append(a); obsB.append(b); qall.append(q)
    if not obsA:
        return np.array([]), np.array([]), np.array([])
    return np.concatenate(obsA), np.concatenate(obsB), np.concatenate(qall)


def build_operator_matrix(obsA, obsB, q, n_bins, bins_H=20, bins_MI=14):
    order = np.argsort(q)
    obsA, obsB, q = obsA[order], obsB[order], q[order]
    edges = np.linspace(0, len(q), n_bins + 1).astype(int)
    rows, centers, ns = [], [], []
    for i in range(n_bins):
        s, e = edges[i], edges[i + 1]
        if e - s < MIN_PER_BIN:
            continue
        a, b = obsA[s:e], obsB[s:e]
        Ha = entropy_hist_1d(a, bins=bins_H)
        Hb = entropy_hist_1d(b, bins=bins_H)
        MIv = mi_hist_2d(a, b, bins=bins_MI)
        rows.append([Ha, Hb, Ha * Ha, Hb * Hb, Ha * Hb, MIv])
        centers.append(float(q[s:e].mean())); ns.append(int(e - s))
    return np.array(rows), np.array(centers), ns


def decompose(v):
    v = np.array(v) / np.linalg.norm(v)
    eq = (v @ d_lin) ** 2 + (v @ d_quad) ** 2
    mi = (v @ d_mi) ** 2
    return float(eq), float(mi), float(max(0.0, 1 - eq - mi))


def run_config(events, charge_sign, observable, n_bins, rng):
    obsA, obsB, q = build_pairs(events, charge_sign, observable, rng)
    if len(q) < MIN_PER_BIN * 4:
        return None
    Mop, centers, ns = build_operator_matrix(obsA, obsB, q, n_bins)
    if Mop.shape[0] < 6:
        return None
    v_null, S, _, _ = extract_v1(Mop)
    vn = v_null / np.linalg.norm(v_null)
    sub234, cos_a = project_234(v_null)
    eq, mi, res = decompose(v_null)
    asym = float(np.abs(Mop[:, 0] - Mop[:, 1]).mean())
    Ha, MIv = Mop[:, 0], Mop[:, 5]
    A1 = np.polyfit(Ha, MIv, 1)
    ss = np.sum((MIv - np.polyval(A1, Ha)) ** 2); st = np.sum((MIv - MIv.mean()) ** 2)
    lin_r2 = float(1 - ss / st) if st > 0 else 0.0
    return dict(charge=int(charge_sign), observable=observable,
                n_pairs=int(len(q)), n_bins=int(Mop.shape[0]),
                pairs_per_bin=int(np.mean(ns)),
                q_range=[float(centers.min()), float(centers.max())],
                mean_abs_HaHb=asym, null_6d=vn.tolist(),
                sub234=sub234, cos_to_attractor=float(abs(cos_a)),
                decomposition=dict(equal_entropy=eq, MI_coupling=mi, residual=res),
                MI_lin_Ha_R2=lin_r2, MI_lin_Ha_slope=float(A1[0]),
                MI_mean=float(MIv.mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    ap.add_argument("--files", nargs="*", default=None)
    ap.add_argument("--events", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()
    t0 = time.time()

    if args.canary:
        files = ["data/strong/f000001.root"]
        n_events = args.events or 80
        seeds = 1
        configs = [(+1, "pt")]
        outname = "s14_force_strong_canary.json"
    else:
        files = args.files or sorted(glob.glob("data/strong/*.root"))
        n_events = args.events
        seeds = args.seeds
        configs = [(+1, "pt"), (-1, "pt"), (+1, "eta"), (-1, "eta")]
        outname = "s14_force_strong_results.json"

    print("=" * 92)
    print(f"S14 STRONG force-operator object: ATLAS DAOD_HION14 Pb-Pb femtoscopy "
          f"({len(files)} file(s), events={n_events or 'all'})")
    print("=" * 92)

    rng0 = np.random.default_rng(0)
    events = []
    for f in files:
        ev = read_tracks(f, n_events, rng0)
        events.extend(ev)
        print(f"  {f}: +{len(ev)} events with selected tracks "
              f"(total {len(events)}); {time.time()-t0:.1f}s")
    mult = [len(e["pt"]) for e in events]
    print(f"selected-track multiplicity: mean {np.mean(mult):.0f} "
          f"median {np.median(mult):.0f} max {np.max(mult)}")

    results = []
    for (cs, obs) in configs:
        per_seed = []
        for sd in range(seeds):
            rng = np.random.default_rng(100 + sd)
            r = run_config(events, cs, obs, N_BINS, rng)
            if r:
                r["seed"] = sd
                per_seed.append(r)
        if not per_seed:
            print(f"  config charge={cs:+d} {obs}: insufficient pairs"); continue
        # cross-seed summary
        cos = np.array([r["cos_to_attractor"] for r in per_seed])
        mi = np.array([r["decomposition"]["MI_coupling"] for r in per_seed])
        eqv = np.array([r["decomposition"]["equal_entropy"] for r in per_seed])
        asym = np.array([r["mean_abs_HaHb"] for r in per_seed])
        r0 = per_seed[0]
        print(f"\n config charge={cs:+d} obs={obs}: "
              f"{r0['n_pairs']} pairs, {r0['n_bins']} q-bins "
              f"(~{r0['pairs_per_bin']}/bin), q in {np.round(r0['q_range'],3)}")
        print(f"   |H_a-H_b| {asym.mean():.3f}+/-{asym.std():.3f}   "
              f"cos->attractor {cos.mean():.3f}+/-{cos.std():.3f}")
        print(f"   decomposition: equal-entropy {eqv.mean():.3f}  "
              f"MI-coupling {mi.mean():.3f}+/-{mi.std():.3f}  "
              f"residual {1-eqv.mean()-mi.mean():.3f}")
        print(f"   --> MI {'IN' if mi.mean() > 0.5 else 'NOT in'} the null "
              f"(self-pole)" if mi.mean() <= 0.5 else
              f"   --> MI IN the null (COUPLING pole!)")
        results.append(dict(config=f"charge{cs:+d}_{obs}",
                            cos_mean=float(cos.mean()), cos_std=float(cos.std()),
                            MI_coupling_mean=float(mi.mean()),
                            MI_coupling_std=float(mi.std()),
                            equal_entropy_mean=float(eqv.mean()),
                            mean_abs_HaHb=float(asym.mean()),
                            per_seed=per_seed))

    out = dict(meta=dict(force="strong",
                         dataset="ATLAS DAOD_HION14 record 80036 data15_hi",
                         construction="two-identical-hadron femtoscopy, q-binned",
                         n_files=len(files), n_events=len(events),
                         pt_window=[PT_LO, PT_HI], eta_max=ETA_MAX, q_max=Q_MAX,
                         max_trk=MAX_TRK, seeds=seeds,
                         elapsed_s=round(time.time() - t0, 2)),
               configs=results)
    json.dump(out, open(outname, "w"), indent=2)
    print(f"\nWrote {outname} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
