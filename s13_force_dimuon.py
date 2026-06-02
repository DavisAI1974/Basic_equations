"""
S13 HEADLINE build #1 -- first REAL gauge-force operator-space object (WEAK).

Decision gate PASSED: CMS Open Data dimuon (record 545, Zmumu_Run2011A) is a
real per-event 2-channel object -- two muons from one Z (weak neutral-current)
decay; their kinematic correlation is intrinsic (Z propagator + energy-momentum
conservation), not invented. This builds the force-operator-space object the
SAME way per_domain_kbk builds the per-domain objects, so it is commensurable
with INFO-023/025/040 and with the LIGO gravity object (INFO-036).

Construction (stated choices):
  - CHANNELS BY CHARGE: channel A = mu+ (Q=+1), channel B = mu- (Q=-1). Physical
    and symmetric, NOT the arbitrary leading/subleading labeling (which would
    manufacture an entropy asymmetry). Observable per channel = muon pT.
  - ENERGY AXIS = dimuon invariant mass M (the Z-exchange scale), equal-count
    bins -> a trajectory. This plays the role that TIME plays in the LIGO object
    and per_domain_kbk's ensemble-H (which bins by time step).
  - Per M-bin, over the event sub-ensemble, compute the 6-op vector
    [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI] on (pT_+, pT_-). Stack -> operator
    matrix. Run extract_v1 (smallest null), project to (Ha^2,Hb^2,HaHb), cos to
    the (-1,-1,+2)/sqrt6 equal-entropy attractor, and the INFO-040 decomposition
    (equal-entropy / MI-coupling / residual). Also fit MI-vs-H family (INFO-025).

Speaking posture (Rule C): I EXPECT the weak object to land on the equal-entropy
attractor (mu+ / mu- channels are symmetric by charge, so H_a ~= H_b), which
would confirm INFO-036's "attractor = equal-marginal-entropy geometry" on a
THIRD kind of real data (after LIGO). But MI may carry real coupling content
(the Z mass constraint links the two muons) -- we wait on where the null sits
and what the MI-vs-H family is. No verdict in advance.

Outputs s13_force_dimuon_results.json.
"""
import sys, json, time
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d

CSV = "data/forces/Zmumu.csv"
MU = 0.1056583745  # muon mass (GeV)
N_BINS = 24
MIN_PER_BIN = 200
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
ATTR3 = np.array([-1., -1., 2.]) / np.sqrt(6.)
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def load_events():
    raw = np.genfromtxt(CSV, delimiter=",", names=True)
    pt1, eta1, phi1, Q1 = raw["pt1"], raw["eta1"], raw["phi1"], raw["Q1"]
    pt2, eta2, phi2, Q2 = raw["pt2"], raw["eta2"], raw["phi2"], raw["Q2"]
    # dimuon invariant mass (massless-muon approx is fine at ~90 GeV; keep mass)
    px1, py1, pz1 = pt1*np.cos(phi1), pt1*np.sin(phi1), pt1*np.sinh(eta1)
    px2, py2, pz2 = pt2*np.cos(phi2), pt2*np.sin(phi2), pt2*np.sinh(eta2)
    E1 = np.sqrt(px1**2+py1**2+pz1**2+MU**2)
    E2 = np.sqrt(px2**2+py2**2+pz2**2+MU**2)
    M = np.sqrt(np.maximum((E1+E2)**2 - (px1+px2)**2 - (py1+py2)**2 - (pz1+pz2)**2, 0.0))
    # assign channels by charge: A=mu+, B=mu-
    plus_is_1 = Q1 > 0
    ptA = np.where(plus_is_1, pt1, pt2)   # mu+
    ptB = np.where(plus_is_1, pt2, pt1)   # mu-
    # drop any event with ambiguous charge (Q1*Q2 not -1)
    ok = (Q1 * Q2) < 0
    return ptA[ok], ptB[ok], M[ok]


def build_operator_matrix(ptA, ptB, M, n_bins, bins_H=20, bins_MI=14):
    """Equal-count M-bins; per bin compute the 6-op vector on (ptA,ptB)."""
    order = np.argsort(M)
    ptA, ptB, M = ptA[order], ptB[order], M[order]
    edges = np.linspace(0, len(M), n_bins + 1).astype(int)
    rows, mbin_centers, ns = [], [], []
    for i in range(n_bins):
        s, e = edges[i], edges[i+1]
        if e - s < MIN_PER_BIN:
            continue
        a, b = ptA[s:e], ptB[s:e]
        Ha = entropy_hist_1d(a, bins=bins_H)
        Hb = entropy_hist_1d(b, bins=bins_H)
        MIv = mi_hist_2d(a, b, bins=bins_MI)
        rows.append([Ha, Hb, Ha*Ha, Hb*Hb, Ha*Hb, MIv])
        mbin_centers.append(float(M[s:e].mean())); ns.append(int(e-s))
    return np.array(rows), np.array(mbin_centers), ns


def decompose(v):
    v = np.array(v) / np.linalg.norm(v)
    eq = (v @ d_lin)**2 + (v @ d_quad)**2
    mi = (v @ d_mi)**2
    res = max(0.0, 1 - eq - mi)
    return float(eq), float(mi), float(res)


def main():
    t0 = time.time()
    ptA, ptB, M = load_events()
    print("=" * 92)
    print(f"S13 WEAK force-operator object: CMS Z->mu+mu- ({len(M)} events, "
          f"M in [{M.min():.1f},{M.max():.1f}] GeV)")
    print("=" * 92)
    Mop, centers, ns = build_operator_matrix(ptA, ptB, M, N_BINS)
    print(f"operator matrix: {Mop.shape[0]} M-bins x 6 ops "
          f"(~{int(np.mean(ns))} events/bin)")
    means = Mop.mean(0)
    print(f"op means: " + "  ".join(f"{OPS[i]}={means[i]:+.3f}" for i in range(6)))
    # marginal entropy symmetry (the equal-entropy premise)
    asym = float(np.abs(Mop[:, 0] - Mop[:, 1]).mean())
    print(f"mean |H_a - H_b| across bins (charge-symmetric premise): {asym:.4f}")

    v_null, S, _, _ = extract_v1(Mop)
    sub234, cos_a = project_234(v_null)
    eq, mi, res = decompose(v_null)
    relation = " ".join(f"{v_null[i]/np.linalg.norm(v_null):+.2f}*{OPS[i]}"
                        for i in range(6)
                        if abs(v_null[i]/np.linalg.norm(v_null)) > 0.15) + " ~ 0"
    print(f"\nnull[0] 6D = {[f'{x:+.3f}' for x in v_null/np.linalg.norm(v_null)]}")
    print(f"null relation: {relation}")
    print(f"[2,3,4] projection cos to (-1,-1,+2)/sqrt6 attractor = {abs(cos_a):.4f}")
    print(f"decomposition: equal-entropy={eq:.3f}  MI-coupling={mi:.3f}  residual={res:.3f}")

    # MI-vs-H family (INFO-025 comparison): MI vs H_a across bins
    Ha, Hb, MIv = Mop[:, 0], Mop[:, 1], Mop[:, 5]
    # simple family probes
    def r2(pred):
        ss = np.sum((MIv - pred)**2); st = np.sum((MIv - MIv.mean())**2)
        return 1 - ss/st if st > 0 else 0.0
    A1 = np.polyfit(Ha, MIv, 1); fit_linHa = r2(np.polyval(A1, Ha))
    diff = (Hb - Ha)
    A2 = np.polyfit(diff, MIv, 2); fit_quaddiff = r2(np.polyval(A2, diff))
    fit_const = r2(np.full_like(MIv, MIv.mean()))
    print(f"\nMI-vs-H family probes (R^2):")
    print(f"  linear in H_a            : {fit_linHa:+.3f}   (slope {A1[0]:+.3f})")
    print(f"  quadratic in (H_b-H_a)   : {fit_quaddiff:+.3f}")
    print(f"  constant                 : {fit_const:+.3f}   (MI mean {MIv.mean():.3f})")

    out = dict(
        meta=dict(force="weak", dataset="CMS Open Data record 545 Zmumu_Run2011A",
                  n_events=int(len(M)), n_bins=int(Mop.shape[0]),
                  events_per_bin=int(np.mean(ns)), elapsed_s=round(time.time()-t0, 2)),
        op_means=means.tolist(), mean_abs_HaHb=asym,
        null_6d=(v_null/np.linalg.norm(v_null)).tolist(),
        null_relation=relation, sub234=sub234, cos_to_attractor=float(abs(cos_a)),
        decomposition=dict(equal_entropy=eq, MI_coupling=mi, residual=res),
        MI_vs_H=dict(linear_Ha_R2=float(fit_linHa), linear_Ha_slope=float(A1[0]),
                     quad_diff_R2=float(fit_quaddiff), const_R2=float(fit_const),
                     MI_mean=float(MIv.mean())),
        M_bin_centers=centers.tolist())
    json.dump(out, open("s13_force_dimuon_results.json", "w"), indent=2)
    print(f"\nWrote s13_force_dimuon_results.json ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
