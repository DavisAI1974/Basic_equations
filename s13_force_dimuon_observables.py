"""
S13 headline build #1A -- weak object across OBSERVABLES (resolve the pt caveat).

Build #1 (s13_force_dimuon.py) used muon pT and found the weak Z->mumu object
sits on the equal-MARGINAL-entropy self-pole (H_a~=H_b, MI not in the null) --
like physics/geology (INFO-040), NOT the coupling pole, and NOT the caricature
"weak~chemistry" hit. Caveat: pT is parity-blind; the weak-specific signature
(parity violation / forward-backward asymmetry) lives in the ANGULAR / signed-
longitudinal observables, where mu+ and mu- distributions could differ.

Mechanism being tested: the self-pole appears BECAUSE mu+/mu- are charge-
symmetric so H(obs_+) ~= H(obs_-). The ONLY way to break it is an observable
whose mu+ vs mu- MARGINAL distributions differ -- exactly the A_FB / parity
signature. So we rebuild the same 6-op M-binned extraction on several
observables and check, per observable: (a) is H_a ~= H_b still, (b) does MI
enter the null, (c) where does the null sit (equal-entropy / MI / residual).

Observables per charge channel (A=mu+, B=mu-):
  pt        transverse momentum (build #1, parity-blind)
  E         energy = pt*cosh(eta)
  pz        SIGNED longitudinal momentum (parity-sensitive: A_FB lives here)
  eta       SIGNED pseudorapidity (parity-sensitive)
  costh_cs  Collins-Soper cos(theta*) per muon (mu+ gets +, mu- gets the value
            for its own direction) -- the canonical weak decay angle

Speaking posture (Rule C): I THINK the self-pole will hold across observables
because A_FB is a small correlation-with-boost effect, not a big marginal
difference, so H_a~=H_b likely survives even for signed pz/eta. But if any
parity-sensitive observable breaks H_a~=H_b or pulls MI into the null, that is
the weak coupling showing up. We wait on the per-observable numbers.

Outputs s13_force_dimuon_observables_results.json.
"""
import json, time
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d

CSV = "data/forces/Zmumu.csv"
MU = 0.1056583745
N_BINS = 24
MIN_PER_BIN = 200
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def load():
    raw = np.genfromtxt(CSV, delimiter=",", names=True)
    pt1, eta1, phi1, Q1 = raw["pt1"], raw["eta1"], raw["phi1"], raw["Q1"]
    pt2, eta2, phi2, Q2 = raw["pt2"], raw["eta2"], raw["phi2"], raw["Q2"]
    px1, py1, pz1 = pt1*np.cos(phi1), pt1*np.sin(phi1), pt1*np.sinh(eta1)
    px2, py2, pz2 = pt2*np.cos(phi2), pt2*np.sin(phi2), pt2*np.sinh(eta2)
    E1 = np.sqrt(px1**2+py1**2+pz1**2+MU**2)
    E2 = np.sqrt(px2**2+py2**2+pz2**2+MU**2)
    M = np.sqrt(np.maximum((E1+E2)**2-(px1+px2)**2-(py1+py2)**2-(pz1+pz2)**2, 0.0))
    QT2 = (px1+px2)**2 + (py1+py2)**2
    # Collins-Soper cos(theta*) (standard dilepton definition), one per event,
    # defined for the negative lepton; mu+ takes -that.
    p1p, p1m = (E1+pz1)/np.sqrt(2), (E1-pz1)/np.sqrt(2)
    p2p, p2m = (E2+pz2)/np.sqrt(2), (E2-pz2)/np.sqrt(2)
    Pz = pz1 + pz2
    costh_cs_neg = (2.0/(M*np.sqrt(M**2+QT2)+1e-12)) * (p1p*p2m - p1m*p2p)
    costh_cs_neg *= np.sign(Pz + 1e-12)  # orient by boost direction
    ok = (Q1*Q2) < 0
    plus1 = Q1 > 0
    def chan(v1, v2):
        return np.where(plus1, v1, v2)[ok], np.where(plus1, v2, v1)[ok]
    obs = {}
    obs["pt"] = chan(pt1, pt2)
    obs["E"] = chan(E1, E2)
    obs["pz"] = chan(pz1, pz2)
    obs["eta"] = chan(eta1, eta2)
    # costh_cs: mu+ channel = -costh_cs_neg, mu- channel = +costh_cs_neg
    cthA = np.where(plus1, -costh_cs_neg, costh_cs_neg)[ok]
    cthB = np.where(plus1, costh_cs_neg, -costh_cs_neg)[ok]
    obs["costh_cs"] = (cthA, cthB)
    return obs, M[ok]


def build_M(a, b, M, n_bins, bins_H=20, bins_MI=14):
    order = np.argsort(M)
    a, b, Ms = a[order], b[order], M[order]
    edges = np.linspace(0, len(Ms), n_bins+1).astype(int)
    rows = []
    for i in range(n_bins):
        s, e = edges[i], edges[i+1]
        if e-s < MIN_PER_BIN:
            continue
        aa, bb = a[s:e], b[s:e]
        Ha = entropy_hist_1d(aa, bins=bins_H); Hb = entropy_hist_1d(bb, bins=bins_H)
        MIv = mi_hist_2d(aa, bb, bins=bins_MI)
        rows.append([Ha, Hb, Ha*Ha, Hb*Hb, Ha*Hb, MIv])
    return np.array(rows)


def decompose(v):
    v = np.array(v)/np.linalg.norm(v)
    eq = (v@d_lin)**2 + (v@d_quad)**2
    mi = (v@d_mi)**2
    return float(eq), float(mi), float(max(0.0, 1-eq-mi))


t0 = time.time()
obs, M = load()
print("="*96)
print(f"S13 WEAK object across observables (CMS Z->mumu, {len(M)} events, "
      f"channels mu+/mu-)")
print("="*96)
print(f"{'observable':10s} {'|Ha-Hb|':>8s} {'MI_mean':>8s} {'cos->attr':>9s} "
      f"{'eqEnt':>7s} {'MI':>6s} {'resid':>7s}  null relation")
results = {}
for name in ["pt", "E", "pz", "eta", "costh_cs"]:
    a, b = obs[name]
    Mop = build_M(a, b, M, N_BINS)
    asym = float(np.abs(Mop[:, 0]-Mop[:, 1]).mean())
    v, S, _, _ = extract_v1(Mop)
    vn = v/np.linalg.norm(v)
    _, cos_a = project_234(v)
    eq, mi, res = decompose(v)
    rel = " ".join(f"{vn[i]:+.2f}*{OPS[i]}" for i in range(6) if abs(vn[i]) > 0.2) + " ~ 0"
    results[name] = dict(n_bins=int(Mop.shape[0]), mean_abs_HaHb=asym,
                         MI_mean=float(Mop[:, 5].mean()),
                         cos_to_attractor=float(abs(cos_a)),
                         eqEnt=eq, MI_coupling=mi, residual=res,
                         null_6d=vn.tolist(), null_relation=rel)
    print(f"{name:10s} {asym:8.4f} {Mop[:,5].mean():8.4f} {abs(cos_a):9.3f} "
          f"{eq:7.3f} {mi:6.3f} {res:7.3f}  {rel}")

reading = (
    "Weak Z->mumu 2-channel object across observables. Self-pole = equal-entropy "
    "+ MI-not-in-null. Parity-sensitive observables (pz, eta, costh_cs) test "
    "whether the weak coupling breaks H_a~=H_b or pulls MI into the null vs the "
    "parity-blind pt/E. See per-observable eqEnt / MI / residual."
)
json.dump(dict(meta=dict(force="weak", n_events=int(len(M)), n_bins=N_BINS,
                         elapsed_s=round(time.time()-t0, 2)),
               results=results, reading=reading),
          open("s13_force_dimuon_observables_results.json", "w"), indent=2)
print(f"\nWrote s13_force_dimuon_observables_results.json ({time.time()-t0:.1f}s)")
