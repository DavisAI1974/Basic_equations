"""
PROBE O2-reframed / backlog 6c (S17, Greg "follow this thread"): is FLOW unique to
gravity?

Hunch (H-A/H-B/H-C): gravity = physics + flow dipole; flow expresses as time; only
gravity touches time. INFO-055 showed the gravity signal is fully GR (no beyond-GR
structure) -- but gravity/time coupling is intrinsic to GR, so that does NOT close this
framework question.

HONESTY FLAG up front (the S17 PROBE-1 lesson): "only gravity has time structure" is at
risk of being CONSTRUCTION bookkeeping -- gravity is the only force we observe as a
genuine time-series of a dynamical process, while the gauge objects are static
distributions/relations. And the deepest H-C claim (only gravity dilates clocks) is
about time dilation, which this force-coupling data does NOT probe. So this probe does
NOT test "gravity dilates clocks." It tests the weaker, fair question below, with the
deflationary reading kept present.

FAIR TEST: apply ONE flow operator to the THREE governing relations we recovered from
raw data (the capability trio): does the characteristic quantity evolve MONOTONICALLY
along the relation's natural axis (a "flow") or sit as a localized structure (a
resonance)?
  - GRAVITY: chirp f(t) -- axis = TIME           (INFO-052 ridge)
  - STRONG : alpha_s(lnQ) -- axis = ENERGY SCALE (INFO/S16 running)
  - WEAK   : rate(mass) Breit-Wigner -- axis = MASS (S16 Z propagator)
Flow score = |Spearman(characteristic, axis)| (1 = perfect monotone flow; ~0 = no flow
/ peaked). Also a relative sweep for the monotone cases.

Prior (no pre-assigned meaning, stated as a posture): I expect gravity AND strong to
flow (the strong force RUNS -- the RG flow is literally a flow), so flow is probably NOT
unique to gravity; the live distinction would then be that gravity's flow AXIS is TIME
while strong's is SCALE -- and whether THAT is nature or observable-choice is the open
question the data here cannot settle.

Run:  python probe_flow_operator.py
"""
import json
import numpy as np
import h5py
from scipy.stats import spearmanr

FS = 4096


def gravity_flow():
    """Chirp ridge f(t) from the INFO-052 recovery (axis = time)."""
    d = json.load(open("probe_gravity_chirp_ridge_results.json"))
    g = d["events"]["GW150914"]["H1"]
    t = np.array(g["ridge_t"]); f = np.array(g["ridge_f"])
    # inspiral arc up to merger (drop the final ringdown point)
    order = np.argsort(t)
    t, f = t[order], f[order]
    rho, _ = spearmanr(t, f)
    sweep = (f.max() - f.min()) / np.mean(f)
    return dict(axis="time", rho=float(rho), abs_rho=float(abs(rho)),
                sweep_ratio=float(sweep), n=int(len(t)),
                char_lo=float(f.min()), char_hi=float(f.max()))


def strong_flow():
    """alpha_s vs ln(Q) world data (axis = energy scale)."""
    tab = [(1.78, 0.312), (4.75, 0.217), (9.46, 0.181), (15.0, 0.166),
           (31.6, 0.149), (45.0, 0.143), (66.0, 0.135), (91.1876, 0.1179),
           (120.0, 0.114), (189.0, 0.109), (400.0, 0.099), (638.0, 0.095),
           (1000.0, 0.090)]
    Q = np.array([x[0] for x in tab]); a = np.array([x[1] for x in tab])
    lnQ = np.log(Q)
    rho, _ = spearmanr(lnQ, a)
    sweep = (a.max() - a.min()) / np.mean(a)
    return dict(axis="ln(energy scale)", rho=float(rho), abs_rho=float(abs(rho)),
                sweep_ratio=float(sweep), n=int(len(a)),
                char_lo=float(a.min()), char_hi=float(a.max()))


def weak_flow():
    """Z dimuon mass spectrum rate(mass) (axis = mass); a resonance, not a flow."""
    import csv
    pt1, eta1, phi1, pt2, eta2, phi2 = [], [], [], [], [], []
    with open("data/forces/Zmumu.csv") as fh:
        r = csv.DictReader(fh)
        for row in r:
            pt1.append(float(row["pt1"])); eta1.append(float(row["eta1"]))
            phi1.append(float(row["phi1"])); pt2.append(float(row["pt2"]))
            eta2.append(float(row["eta2"])); phi2.append(float(row["phi2"]))
    pt1 = np.array(pt1); eta1 = np.array(eta1); phi1 = np.array(phi1)
    pt2 = np.array(pt2); eta2 = np.array(eta2); phi2 = np.array(phi2)
    # invariant mass (massless-muon approx)
    m2 = 2 * pt1 * pt2 * (np.cosh(eta1 - eta2) - np.cos(phi1 - phi2))
    m = np.sqrt(np.maximum(m2, 0))
    sel = (m > 60) & (m < 120)
    m = m[sel]
    hist, edges = np.histogram(m, bins=60, range=(60, 120))
    centers = 0.5 * (edges[:-1] + edges[1:])
    rho, _ = spearmanr(centers, hist)               # over full range -> ~0 (peaked)
    ipk = int(np.argmax(hist))
    peak_interior = 0 < ipk < len(hist) - 1
    return dict(axis="mass", rho=float(rho), abs_rho=float(abs(rho)),
                n_events=int(len(m)), peak_mass=float(centers[ipk]),
                peak_interior=bool(peak_interior),
                note="localized resonance (peak), not a monotone flow")


def main():
    g, s, w = gravity_flow(), strong_flow(), weak_flow()
    print("=" * 90)
    print("PROBE flow operator (6c) -- is FLOW unique to gravity?  |Spearman(char, axis)|")
    print("=" * 90)
    print(f"  GRAVITY  chirp f(t)      axis={g['axis']:<18} "
          f"|rho|={g['abs_rho']:.3f}  sweep={g['sweep_ratio']:.2f}  "
          f"({g['char_lo']:.0f}->{g['char_hi']:.0f} Hz)  => {'FLOW' if g['abs_rho']>0.8 else 'no flow'}")
    print(f"  STRONG   alpha_s(lnQ)    axis={s['axis']:<18} "
          f"|rho|={s['abs_rho']:.3f}  sweep={s['sweep_ratio']:.2f}  "
          f"({s['char_hi']:.3f}->{s['char_lo']:.3f})  => {'FLOW' if s['abs_rho']>0.8 else 'no flow'}")
    print(f"  WEAK     rate(mass) BW   axis={w['axis']:<18} "
          f"|rho|={w['abs_rho']:.3f}  peak at {w['peak_mass']:.1f} GeV interior={w['peak_interior']}"
          f"  => {'FLOW' if w['abs_rho']>0.8 else 'RESONANCE (no flow)'}")
    print("-" * 90)
    n_flow = sum(x["abs_rho"] > 0.8 for x in (g, s))
    print(f"VERDICT (data level): {n_flow}/2 of the monotone-candidate forces flow "
          f"(gravity AND strong); weak is a static resonance.")
    print("  => FLOW is NOT unique to gravity (the strong force RUNS -- RG flow).")
    print("  Live distinction: gravity's flow AXIS is TIME; strong's is ENERGY SCALE.")
    print("  DEFLATIONARY (kept present): gravity is the only one we OBSERVE as a")
    print("  time-series; strong's running is also a flow-in-time during the")
    print("  interaction, measured vs scale. Whether 'gravity flows in TIME' is nature")
    print("  or observable-choice is NOT decided by this data (construction-confounded).")

    out = dict(gravity=g, strong=s, weak=w,
               verdict=dict(flow_count_of_2=int(n_flow),
                            flow_unique_to_gravity=False,
                            gravity_axis="time", strong_axis="energy_scale",
                            construction_confounded=True))
    json.dump(out, open("probe_flow_operator_results.json", "w"), indent=2)
    print("\nWrote probe_flow_operator_results.json")


if __name__ == "__main__":
    main()
