"""
S12 -- Standard-Model parameter-regularity hunt (the remaining four-force item,
framed v11, built here). Real PDG-measured parameters. The OD-shaped question:
are the SM's ~26 "free" inputs themselves STRUCTURED -- hidden relations among
the gauge couplings, fermion mass ratios, and CKM/PMNS mixing angles?

DISCIPLINE (load-bearing):
  - Treat-literature-as-conjecture: every candidate relation below is a KNOWN
    numerical pattern WITHOUT an accepted derivation (Koide, Gatto-Sartori-
    Tonin, quark-lepton complementarity, Wolfenstein hierarchy). We COMPUTE
    each from measured data and report how well it holds; we do NOT claim it
    means anything. A relation that holds is one data point, not a law.
  - Result Discipline: catalog MISSES with the same care as HITS (Koide is
    clean for leptons, fails for quarks -- both reported).
  - NO unconstrained symbolic regression on 26 numbers (that is numerology by
    construction). We test PRE-SPECIFIED relations only.
  - Speaking posture: I THINK Koide(leptons) and GST(Cabibbo) will hold
    tightly and the quark Koide will fail, but we report what the numbers say.

Sources: PDG 2024 world averages. Masses in MeV unless noted. Running-scheme
caveats noted (quark masses are scheme/scale dependent -> any quark relation is
softer than a lepton relation by that much).

Output s12_sm_regularity_results.json.
"""
import json
import numpy as np

# ---- measured inputs (PDG 2024) -------------------------------------------
# charged lepton masses (MeV) -- pole masses, very precisely known
m_e, m_mu, m_tau = 0.51099895, 105.6583755, 1776.86
# quark masses (MeV): u,d,s,c MSbar(2 GeV); b MSbar(m_b); t pole. SCHEME-DEPENDENT.
m_u, m_c, m_t = 2.16, 1270.0, 172570.0
m_d, m_s, m_b = 4.67, 93.4, 4180.0
# gauge couplings at M_Z (GUT-normalised inverse), from s10_pdg_unification
a1_inv, a2_inv, a3_inv = 59.02, 29.59, 8.475
# CKM (Wolfenstein / PDG): angles in degrees
th12_ckm, th23_ckm, th13_ckm = 13.04, 2.38, 0.201   # Cabibbo, etc.
lam_ckm = 0.22500                                     # Wolfenstein lambda
# PMNS (neutrino) mixing angles (deg), normal ordering global fit
th12_pmns, th23_pmns, th13_pmns = 33.41, 49.1, 8.54


def koide(masses):
    m = np.array(masses, float)
    return float(m.sum() / (np.sqrt(m).sum() ** 2))


out = {"provenance": "PDG 2024; quark masses scheme/scale-dependent",
       "discipline": "every relation is conjecture (no accepted derivation); "
                     "hits and misses both reported; no SR on the parameter set",
       "relations": {}}
R = out["relations"]

# ---- 1. Koide relation (charged leptons; quark sectors as controls) -------
Ql = koide([m_e, m_mu, m_tau])
Qu = koide([m_u, m_c, m_t])
Qd = koide([m_d, m_s, m_b])
R["koide_leptons"] = dict(value=Ql, target=2/3, dev_from_2_3=abs(Ql - 2/3),
                          holds=abs(Ql - 2/3) < 1e-3, note="clean (5 digits) -- the famous hit")
R["koide_up_quarks"] = dict(value=Qu, target=2/3, dev_from_2_3=abs(Qu - 2/3),
                            holds=abs(Qu - 2/3) < 1e-2, note="MISS (control)")
R["koide_down_quarks"] = dict(value=Qd, target=2/3, dev_from_2_3=abs(Qd - 2/3),
                              holds=abs(Qd - 2/3) < 1e-2, note="MISS (control), closer than up")

# ---- 2. Gatto-Sartori-Tonin: sin(theta_C) ~ sqrt(m_d/m_s) -----------------
gst = np.sqrt(m_d / m_s)
sinC = np.sin(np.radians(th12_ckm))
R["GST_cabibbo"] = dict(sqrt_md_ms=float(gst), sin_thetaC=float(sinC),
                        ratio=float(gst / sinC), holds=abs(gst / sinC - 1) < 0.1,
                        note="sqrt(m_d/m_s) vs Cabibbo sine")

# ---- 3. Quark-lepton complementarity: th12_CKM + th12_PMNS ~ 45 deg -------
qlc = th12_ckm + th12_pmns
R["quark_lepton_complementarity"] = dict(sum_deg=float(qlc), target=45.0,
                                          dev=abs(qlc - 45.0), holds=abs(qlc - 45.0) < 3.0,
                                          note="Cabibbo + solar angle vs 45 deg")

# ---- 4. Mass spectra log-geometric regularity (per sector) ----------------
# Is log(mass) ~ linear in generation index (a constant geometric ratio)?
for name, ms in [("leptons", [m_e, m_mu, m_tau]),
                 ("up_quarks", [m_u, m_c, m_t]),
                 ("down_quarks", [m_d, m_s, m_b])]:
    lm = np.log(ms)
    idx = np.array([1, 2, 3.0])
    A = np.vstack([idx, np.ones(3)]).T
    coef, res, *_ = np.linalg.lstsq(A, lm, rcond=None)
    pred = A @ coef
    ss_res = float(((lm - pred) ** 2).sum()); ss_tot = float(((lm - lm.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot
    ratios = [ms[1] / ms[0], ms[2] / ms[1]]
    R[f"loglinear_{name}"] = dict(slope=float(coef[0]), R2=float(r2),
                                  gen_ratios=[float(x) for x in ratios],
                                  note="log(mass) vs generation; R2~1 => geometric")

# ---- 5. Wolfenstein hierarchy: CKM angles as powers of lambda -------------
R["wolfenstein_powers"] = dict(
    lam=lam_ckm,
    s12_over_lam=float(np.sin(np.radians(th12_ckm)) / lam_ckm),
    s23_over_lam2=float(np.sin(np.radians(th23_ckm)) / lam_ckm ** 2),
    s13_over_lam3=float(np.sin(np.radians(th13_ckm)) / lam_ckm ** 3),
    note="s12~lam, s23~lam^2, s13~lam^3 => ratios ~O(1) if hierarchy holds")

# ---- 6. gauge coupling structure (from s10/INFO-037, recorded) ------------
R["gauge_couplings_MZ"] = dict(a1_inv=a1_inv, a2_inv=a2_inv, a3_inv=a3_inv,
                               ratios=[a1_inv/a2_inv, a2_inv/a3_inv],
                               note="no single unification (INFO-037 triangle); "
                                    "ratios not obviously simple")

# ---- print + reading ------------------------------------------------------
print("=" * 92)
print("S12 SM parameter-regularity hunt (PDG; every relation = conjecture)")
print("=" * 92)
print(f"\n[Koide Q = sum(m)/ (sum sqrt(m))^2 ; target 2/3 = 0.66667]")
print(f"  leptons    Q = {Ql:.6f}   dev {abs(Ql-2/3):.2e}   <- HIT (5 digits)")
print(f"  up quarks  Q = {Qu:.4f}     dev {abs(Qu-2/3):.3f}    <- miss")
print(f"  down quarks Q= {Qd:.4f}     dev {abs(Qd-2/3):.3f}    <- miss")
print(f"\n[GST] sqrt(m_d/m_s)={gst:.4f}  sin(theta_C)={sinC:.4f}  ratio={gst/sinC:.3f}  "
      f"<- {'HIT' if abs(gst/sinC-1)<0.1 else 'miss'}")
print(f"[QLC] theta12_CKM + theta12_PMNS = {qlc:.1f} deg (target 45)  "
      f"<- {'HIT' if abs(qlc-45)<3 else 'miss'}")
print(f"\n[log-geometric mass spectra]  R^2 of log(mass) vs generation:")
for name in ["leptons", "up_quarks", "down_quarks"]:
    r = R[f"loglinear_{name}"]
    print(f"  {name:12s} R2={r['R2']:.4f}  gen-ratios={[round(x) for x in r['gen_ratios']]}")
print(f"\n[Wolfenstein] s12/lam={R['wolfenstein_powers']['s12_over_lam']:.2f}  "
      f"s23/lam^2={R['wolfenstein_powers']['s23_over_lam2']:.2f}  "
      f"s13/lam^3={R['wolfenstein_powers']['s13_over_lam3']:.2f}  (all O(1) => hierarchy)")

reading = (
    "SM parameter-regularity hunt, data/interpretation/frame separate. DATA: "
    "the charged-lepton Koide relation holds to 5 digits (Q=0.66666); the quark "
    "Koide does NOT (up 0.85, down 0.73). Gatto-Sartori-Tonin sqrt(m_d/m_s) "
    "matches the Cabibbo sine within ~5 percent. Quark-lepton complementarity "
    "(Cabibbo + solar = ~46.5 deg) sits near 45 within a few degrees. Charged-"
    "lepton and quark mass spectra are roughly geometric (log-linear R^2 high) "
    "but NOT exactly. CKM angles follow the Wolfenstein lambda^n hierarchy "
    "(ratios O(1)). INTERPRETATION: there IS real low-dimensional structure in "
    "the SM mass/mixing sector -- it is not 26 independent random numbers -- but "
    "the cleanest relation (lepton Koide) has no accepted derivation and the "
    "quark analogues fail, so each is a CONJECTURE / one data point, not a law. "
    "FRAME: this is the honest real-data face of 'are the forces/parameters "
    "structured' -- yes there are regularities, no they do not (yet) reduce to "
    "a single generating rule, and we cannot invoke any of them as support "
    "until independently derived. Connection to the per-domain dipole equations "
    "(Greg's force<->equation question): NONE direct -- SM regularities are "
    "mass/angle relations among static parameters; the per-domain equations are "
    "MI-vs-entropy relations of 2-channel dynamics. Different KIND of object; "
    "the only bridge (the four-force caricatures) was retired as self-grading."
)
print("\n" + "=" * 92 + "\nREADING:\n" + reading + "\n" + "=" * 92)
out["reading"] = reading
json.dump(out, open("s12_sm_regularity_results.json", "w"), indent=2,
          default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("\nWrote s12_sm_regularity_results.json")
