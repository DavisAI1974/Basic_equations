"""
S13 headline build #2B -- EM g2(0) coincidence object (closes the INFO-048 gap).

INFO-048 built the EM HBT object from windowed count-rate entropy at 20-100us
bins and flagged a caveat: that scale might capture slow intensity drift rather
than the HBT g2(0) quantum bunching. This build measures g2(tau) DIRECTLY and
rebuilds the operator object ACROSS bin scales spanning the coherence time, to
answer: (1) is the bunching real (g2(0)>1)? (2) at the coherence scale where the
bunching is strongest, does the inter-detector MI ENTER the null (genuine
coupling, like simulated biology) or stay an ACTIVE variable (like LIGO merger,
weak Z, and the 20-100us EM build)?

Rule D note (correcting INFO-048): recon shows g2(0)=1.83-1.92 with coherence
half-width ~2us (decays to 1.00 by +-400us). So the 20-100us build sat ABOVE
the coherence scale -- it captured the bunching only partially (the slow tail).
The caveat "20-100us = drift not g2" was INCOMPLETE; the bunching coherence is
~2us, and we test the object AT that scale here.

g2(tau) via FFT-free binned cross-correlation:
  g2(tau) = <n_A(t) n_B(t+tau)> / (<n_A><n_B>),  n binned at dt.
Operator object: window the coherence-scale count series; per window
[H_a,H_b,H_a^2,H_b^2,H_a*H_b,MI] on (n_A,n_B); extract_v1 null; project_234;
INFO-040 decomposition. Same machinery as INFO-047/048.

Speaking posture (Rule C): I THINK the bunching MI will STILL stay out of the
null even at the coherence scale (matching every real object so far -- MI active,
not a constraint), with biology remaining the lone MI-in-null case. But the
coherence scale is where it has the best chance to enter; we wait on the numbers.

Outputs s13_force_hbt_g2_results.json.
"""
import json, time
import numpy as np
from kbk_pipeline import mi_hist_2d, extract_v1, project_234
from per_domain_kbk import entropy_hist_1d

RES = 156.25e-12
DT_GRID = [1e-6, 2e-6, 4e-6, 50e-6]   # span the ~2us coherence scale
N_WIN = 200
SOURCES = {
    "split":   ("data/forces/em/ext/splitThermal_rawtags.txt", 11, 15),
    "uncorr":  ("data/forces/em/ext/uncorrelatedThermal_rawtags.txt", 15, 16),
}
d_lin = np.array([1, -1, 0, 0, 0, 0]) / np.sqrt(2)
d_quad = np.array([0, 0, -1, -1, 2, 0]) / np.sqrt(6)
d_mi = np.array([0, 0, 0, 0, 0, 1.0])
OPS = ["H_a", "H_b", "H_a^2", "H_b^2", "H_a*H_b", "MI"]


def counts(path, chA, chB, dt):
    d = np.loadtxt(path); ch = d[:, 0].astype(int); t = d[:, 1]
    t0 = t.min(); T = (t.max() - t0) * RES; nb = int(T / dt)
    a = np.histogram((t[ch == chA]-t0)*RES, bins=nb, range=(0, T))[0].astype(float)
    b = np.histogram((t[ch == chB]-t0)*RES, bins=nb, range=(0, T))[0].astype(float)
    return a, b, T


def g2_curve(a, b, L):
    ma, mb, nb = a.mean(), b.mean(), len(a)
    lags = np.arange(-L, L+1)
    g = np.array([np.mean(a[max(0, -k):nb-max(0, k)] * b[max(0, k):nb-max(0, -k)])
                  for k in lags]) / (ma*mb)
    return lags, g


def build_M(a, b, n_win, bins_H=20, bins_MI=12):
    W = len(a)//n_win; rows = []
    for i in range(n_win):
        sa, sb = a[i*W:(i+1)*W], b[i*W:(i+1)*W]
        if len(sa) < 10:
            continue
        Ha = entropy_hist_1d(sa, bins=bins_H); Hb = entropy_hist_1d(sb, bins=bins_H)
        rows.append([Ha, Hb, Ha*Ha, Hb*Hb, Ha*Hb, mi_hist_2d(sa, sb, bins=bins_MI)])
    return np.array(rows)


def decompose(v):
    v = np.array(v)/np.linalg.norm(v)
    eq = (v@d_lin)**2 + (v@d_quad)**2; mi = (v@d_mi)**2
    return float(eq), float(mi), float(max(0.0, 1-eq-mi))


t0 = time.time()
print("="*100)
print("S13 EM g2(0) coincidence object -- HBT bunching + operator object across the coherence scale")
print("="*100)
out = {"meta": dict(force="EM", dataset="Zenodo 5113016 HBT", dt_grid=DT_GRID,
                    n_win=N_WIN), "results": {}}
for src, (path, chA, chB) in SOURCES.items():
    out["results"][src] = {"g2": {}, "object": {}}
    # --- g2(0) + coherence at 1us bins, +-400us ---
    a1, b1, T = counts(path, chA, chB, 1e-6)
    lags, g = g2_curve(a1, b1, 400)
    g0 = float(g[len(g)//2]); gfar = float(np.mean(np.r_[g[:20], g[-20:]]))
    half = 1 + (g0-1)/2
    ic = int(np.argmin(np.abs(g[len(g)//2:] - half)))
    out["results"][src]["g2"] = dict(g2_0=g0, g2_far=gfar,
                                     coherence_halfwidth_us=float(abs(lags[len(g)//2+ic])))
    print(f"\n--- {src} (ch{chA}/ch{chB}) : g2(0)={g0:.3f}  g2(far)={gfar:.3f}  "
          f"coherence~{abs(lags[len(g)//2+ic])}us ---")
    print(f"{'dt(us)':>7s} {'meanA':>7s} {'|Ha-Hb|':>8s} {'MI_mean':>8s} "
          f"{'cos->attr':>9s} {'eqEnt':>7s} {'MI':>6s} {'resid':>7s}  null relation")
    for dt in DT_GRID:
        a, b, _ = counts(path, chA, chB, dt)
        Mop = build_M(a, b, N_WIN)
        asym = float(np.abs(Mop[:, 0]-Mop[:, 1]).mean())
        v, S, _, _ = extract_v1(Mop); vn = v/np.linalg.norm(v)
        _, cos_a = project_234(v); eq, mi, res = decompose(v)
        rel = " ".join(f"{vn[i]:+.2f}*{OPS[i]}" for i in range(6) if abs(vn[i]) > 0.2) + " ~ 0"
        out["results"][src]["object"][f"dt_{dt*1e6:.0f}us"] = dict(
            mean_count_A=float(a.mean()), mean_abs_HaHb=asym, MI_mean=float(Mop[:, 5].mean()),
            cos_to_attractor=float(abs(cos_a)), eqEnt=eq, MI_coupling=mi, residual=res,
            null_relation=rel)
        print(f"{dt*1e6:7.0f} {a.mean():7.3f} {asym:8.4f} {Mop[:,5].mean():8.4f} "
              f"{abs(cos_a):9.3f} {eq:7.3f} {mi:6.3f} {res:7.3f}  {rel}")

out["reading"] = (
    "EM HBT g2(0) confirmed >1 (Bose bunching, real). Operator object tested "
    "across bin scales spanning the ~2us coherence. Key: does the bunching MI "
    "enter the null at the coherence scale, or stay active like every other real "
    "object (LIGO/weak/20-100us-EM)? Compare MI_coupling across dt. Corrects "
    "INFO-048's coherence-scale caveat (Rule D)."
)
json.dump(out, open("s13_force_hbt_g2_results.json", "w"), indent=2)
print(f"\nWrote s13_force_hbt_g2_results.json ({time.time()-t0:.1f}s)")
