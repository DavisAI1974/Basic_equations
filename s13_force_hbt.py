"""
S13 headline build #2 -- first REAL EM force-operator object via Hanbury
Brown-Twiss (HBT) intensity interferometry. SAME construction type as LIGO
(two physical detectors sampling one common field), so this is the clean
WITHIN-TYPE test of INFO-036: is the (-1,-1,+2) equal-entropy attractor a
property of the two-detector construction geometry, or of gravity?

Data: Zenodo 5113016 (g2data), raw two-channel photon timetags (channel,
clock-cycle; res=156.25ps). Two sources:
  SPLIT thermal  -- one thermal field split 50/50 to detectors A(ch11),B(ch15):
                    a COMMON field -> HBT bunching -> the two count series are
                    correlated (measured corr ~0.33 at ~20-100us). The "signal".
  UNCORRELATED thermal -- detectors A(ch15),B(ch16) on independent fields: the
                    CONTROL (no common field). The EM analogue of LIGO off-source
                    noise.

Construction (matches LIGO/per_domain): bin each detector's photon arrivals into
fine bins of width dt -> count series n_A[k], n_B[k] (the intensity). Window the
series; per window compute H(n_A), H(n_B) and MI(n_A,n_B); stack -> 6-op operator
matrix [H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI]; extract_v1 null, project to
(Ha^2,Hb^2,HaHb), cos to (-1,-1,+2)/sqrt6, and the INFO-040 decomposition.

Speaking posture (Rule C): I THINK the COMMON-field (split) object will carry
real MI (shared intensity) the way LIGO's merger does, and -- because the two
detectors have DIFFERENT rates (A ~1.8x B) -- H_a may NOT equal H_b, so it may
sit OFF the pure equal-entropy attractor, unlike the charge-symmetric dimuon
object. The uncorrelated control should show much weaker inter-channel MI. But
we wait on the numbers; no verdict in advance, signal and control reported alike.

Outputs s13_force_hbt_results.json.
"""
import json, time
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d

RES = 156.25e-12
DT_GRID = [20e-6, 50e-6, 100e-6]
N_WIN = 200
SOURCES = {
    "split_correlated":   ("data/forces/em/ext/splitThermal_rawtags.txt", 11, 15),
    "uncorrelated_control": ("data/forces/em/ext/uncorrelatedThermal_rawtags.txt", 15, 16),
}
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def load_counts(path, chA, chB, dt):
    d = np.loadtxt(path)
    ch = d[:, 0].astype(int); t = d[:, 1]
    t0 = t.min(); T = (t.max() - t0) * RES
    tA = (t[ch == chA] - t0) * RES
    tB = (t[ch == chB] - t0) * RES
    nb = int(T / dt)
    a = np.histogram(tA, bins=nb, range=(0, T))[0].astype(float)
    b = np.histogram(tB, bins=nb, range=(0, T))[0].astype(float)
    return a, b, T, len(tA), len(tB)


def build_M(a, b, n_win, bins_H=24, bins_MI=14):
    W = len(a) // n_win
    rows = []
    for i in range(n_win):
        sa = a[i*W:(i+1)*W]; sb = b[i*W:(i+1)*W]
        if len(sa) < 10:
            continue
        Ha = entropy_hist_1d(sa, bins=bins_H); Hb = entropy_hist_1d(sb, bins=bins_H)
        MIv = mi_hist_2d(sa, sb, bins=bins_MI)
        rows.append([Ha, Hb, Ha*Ha, Hb*Hb, Ha*Hb, MIv])
    return np.array(rows)


def decompose(v):
    v = np.array(v)/np.linalg.norm(v)
    eq = (v@d_lin)**2 + (v@d_quad)**2
    mi = (v@d_mi)**2
    return float(eq), float(mi), float(max(0.0, 1-eq-mi))


t0 = time.time()
print("="*100)
print("S13 EM force object via HBT (Zenodo 5113016): split(common field) vs uncorrelated(control)")
print("="*100)
out = {"meta": dict(force="EM", dataset="Zenodo 5113016 HBT timetags",
                    dt_grid=DT_GRID, n_win=N_WIN), "results": {}}
for src, (path, chA, chB) in SOURCES.items():
    out["results"][src] = {}
    print(f"\n--- {src} (ch{chA}/ch{chB}) ---")
    print(f"{'dt(us)':>7s} {'corr(nA,nB)':>11s} {'|Ha-Hb|':>8s} {'MI_mean':>8s} "
          f"{'cos->attr':>9s} {'eqEnt':>7s} {'MI':>6s} {'resid':>7s}  null relation")
    for dt in DT_GRID:
        a, b, T, NA, NB = load_counts(path, chA, chB, dt)
        cc = float(np.corrcoef(a, b)[0, 1])
        Mop = build_M(a, b, N_WIN)
        asym = float(np.abs(Mop[:, 0]-Mop[:, 1]).mean())
        v, S, _, _ = extract_v1(Mop)
        vn = v/np.linalg.norm(v)
        _, cos_a = project_234(v)
        eq, mi, res = decompose(v)
        rel = " ".join(f"{vn[i]:+.2f}*{OPS[i]}" for i in range(6) if abs(vn[i]) > 0.2) + " ~ 0"
        out["results"][src][f"dt_{dt*1e6:.0f}us"] = dict(
            corr_nA_nB=cc, rate_A=NA/T, rate_B=NB/T, n_bins=int(len(a)),
            n_win=int(Mop.shape[0]), mean_abs_HaHb=asym, MI_mean=float(Mop[:, 5].mean()),
            cos_to_attractor=float(abs(cos_a)), eqEnt=eq, MI_coupling=mi, residual=res,
            null_6d=vn.tolist(), null_relation=rel)
        print(f"{dt*1e6:7.0f} {cc:11.4f} {asym:8.4f} {Mop[:,5].mean():8.4f} "
              f"{abs(cos_a):9.3f} {eq:7.3f} {mi:6.3f} {res:7.3f}  {rel}")

out["reading"] = (
    "EM HBT 2-channel object, same construction type as LIGO. split=common field "
    "(HBT bunching), uncorrelated=control. Compare cos-to-attractor and the "
    "eq-entropy/MI/residual decomposition between split and control, and against "
    "LIGO gravity (INFO-036): does the equal-entropy attractor reflect detector-"
    "pair geometry, and does the common-field MI enter the null or stay active?"
)
json.dump(out, open("s13_force_hbt_results.json", "w"), indent=2)
print(f"\nWrote s13_force_hbt_results.json ({time.time()-t0:.1f}s)")
