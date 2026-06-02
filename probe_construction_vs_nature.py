"""
PROBE (backlog #1, S17): construction-vs-nature of the equal-entropy substrate.

THE question (Greg, S16; flagged "clean, decisive"): everything we measure sits ON
the equal-entropy region (the (1,1,2)/sqrt6 protocol direction / -(H_a-H_b)^2
channel-correlation identity). Is equal-marginal-entropy a FORCED property of
physical 2-channel observables (nature), or our construction choice (bookkeeping)?

THE LEVER (clean + decisive): scale ONE channel by a constant factor s, b -> s*b.
This transform is:
  * MI-INVARIANT.  mi_hist_2d uses adaptive (data-range) bin edges, so scaling b
    rescales the edges proportionally and the 2D histogram is IDENTICAL =>
    MI(a, s*b) = MI(a, b) exactly. The coupling / physics is untouched.
  * a PURE UNITS / REPRESENTATION choice. Differential entropy shifts by a known
    constant: H(s*b) = H(b) + ln|s|. So scaling injects pure construction asymmetry
    into the marginal entropies while leaving the underlying system invariant.
If a physically-meaningless units choice moves objects OFF the substrate, then
"sitting on the substrate" depends on representation => bookkeeping (at least in
part). If the null direction is robust to per-channel scale, that points toward
something representation-independent => more nature-like.

NOTE on centering (honest caveat): extract_v1 CENTERS columns before SVD, so the
constant ln(s) added to the H_b column vanishes; the effect on the null reaches the
quadratic/cross terms (H_b^2, H_a*H_b) through their centered fluctuations. So the
test is whether the centered-covariance null direction is sensitive to a per-channel
scale, which is exactly the representation-dependence we care about.

We report the CANONICAL prior metric cos(v_null_234, (1,1,2)/sqrt6) [= project_234,
what s10/s11 reported as cos-to-attractor] AND cos to (-1,-1,2)/sqrt6 (the
-(H_a-H_b)^2 identity), the realized mean |H_a - H_b|, |MI coef in null|, and mean
MI (to confirm scale-invariance empirically). Sim: 4 heterogeneous generators x
>=3 seeds. Real: LIGO H1/L1 noise-only segments, same scaling lever post-whitening.

Speaking posture (Rule C): no verdict in advance; no verdict on first look.

Run:  python probe_construction_vs_nature.py --canary
      python probe_construction_vs_nature.py
"""
import sys
import json
import time
import numpy as np

from kbk_pipeline import (
    vasicek_entropy, mi_hist_2d, compute_operator_matrix, extract_v1,
)

CANARY = "--canary" in sys.argv

T_TARGET = np.array([1.0, 1.0, 2.0]) / np.sqrt(6.0)     # (+1,+1,+2) protocol dir
T_IDENT = np.array([-1.0, -1.0, 2.0]) / np.sqrt(6.0)    # -(H_a-H_b)^2 identity dir


# ----------------------------------------------------------------------------
# Heterogeneous two-channel generators (Session-3 attractor members).
# Each returns (a, b), two correlated 1-D streams, BEFORE the scaling lever.
# ----------------------------------------------------------------------------
def gen_ou(rng, n, rho=0.6, gamma=0.5, sigma=1.0, dt=0.05):
    a = np.zeros(n); b = np.zeros(n)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    for t in range(1, n):
        z = L @ rng.standard_normal(2)
        a[t] = a[t-1] - gamma*a[t-1]*dt + sigma*np.sqrt(dt)*z[0]
        b[t] = b[t-1] - gamma*b[t-1]*dt + sigma*np.sqrt(dt)*z[1]
    return a, b


def gen_ar1_laplace(rng, n, phi=0.7, rho=0.6):
    # AR(1) with correlated Laplace innovations (non-Gaussian, heavy-ish tails)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    a = np.zeros(n); b = np.zeros(n)
    for t in range(1, n):
        u = rng.uniform(-0.5, 0.5, 2)
        lap = -np.sign(u) * np.log1p(-2*np.abs(u))   # Laplace via inverse-CDF
        z = L @ lap
        a[t] = phi*a[t-1] + z[0]
        b[t] = phi*b[t-1] + z[1]
    return a, b


def gen_logistic(rng, n, r=3.9, eps=0.15):
    # Two coupled logistic maps (deterministic chaos + diffusive coupling)
    a = np.zeros(n); b = np.zeros(n)
    a[0], b[0] = rng.uniform(0.1, 0.9, 2)
    for t in range(1, n):
        fa = r*a[t-1]*(1-a[t-1]); fb = r*b[t-1]*(1-b[t-1])
        a[t] = (1-eps)*fa + eps*fb
        b[t] = (1-eps)*fb + eps*fa
    return a - a.mean(), b - b.mean()


def gen_sine_noise(rng, n, f=0.02, amp=1.0, noise=0.3):
    t = np.arange(n)
    common = amp*np.sin(2*np.pi*f*t)
    a = common + noise*rng.standard_normal(n)
    b = common + noise*rng.standard_normal(n)
    return a, b


GENERATORS = {
    "OU_corr": gen_ou,
    "AR1_laplace": gen_ar1_laplace,
    "logistic_chaos": gen_logistic,
    "sine_plus_noise": gen_sine_noise,
}


# ----------------------------------------------------------------------------
def null_report(M):
    """Extract null; return dict of metrics including cos to both targets."""
    v, S, mu, Mc = extract_v1(M)
    sub = v[[2, 3, 4]].copy()
    nrm = np.linalg.norm(sub)
    if nrm > 0:
        sub = sub / nrm
    cos_target = abs(float(np.dot(sub, T_TARGET)))
    cos_ident = abs(float(np.dot(sub, T_IDENT)))
    return dict(
        v6=[float(x) for x in v],
        sub234=[float(x) for x in sub],
        cos_to_112=cos_target,
        cos_to_m1m12=cos_ident,
        mi_coef=abs(float(v[5])),
        mean_Ha=float(mu[0]),
        mean_Hb=float(mu[1]),
        mean_MI=float(M[:, 5].mean()),
    )


def per_window_asym(M):
    """Mean |H_a - H_b| across windows (the INFO-038 asymmetry knob)."""
    return float(np.mean(np.abs(M[:, 0] - M[:, 1])))


def run_sim(scales, seeds, n, win, stride, mi_bins):
    out = {}
    for gname, gfn in GENERATORS.items():
        out[gname] = {}
        for s in scales:
            cells = []
            for seed in seeds:
                rng = np.random.default_rng(seed)
                a, b = gfn(rng, n)
                b_scaled = s * b
                M = compute_operator_matrix(a, b_scaled, win, stride,
                                            mi_bins=mi_bins)
                if M.shape[0] < 6:
                    continue
                rep = null_report(M)
                rep["asym_meanabs"] = per_window_asym(M)
                rep["n_windows"] = int(M.shape[0])
                rep["seed"] = seed
                cells.append(rep)
            # cross-seed aggregate
            if cells:
                ct = [c["cos_to_112"] for c in cells]
                ci = [c["cos_to_m1m12"] for c in cells]
                asy = [c["asym_meanabs"] for c in cells]
                mi = [c["mean_MI"] for c in cells]
                out[gname][f"{s:g}"] = dict(
                    scale=s, ln_s=float(np.log(s)),
                    cos112_mean=float(np.mean(ct)), cos112_std=float(np.std(ct)),
                    cosident_mean=float(np.mean(ci)), cosident_std=float(np.std(ci)),
                    asym_mean=float(np.mean(asy)),
                    meanMI_mean=float(np.mean(mi)),
                    per_seed=cells,
                )
    return out


# ----------------------------------------------------------------------------
# Real LIGO: scale one channel post-whitening, same lever.
# ----------------------------------------------------------------------------
def run_ligo(scales, win_s, stride_s):
    import h5py
    from scipy.signal import butter, filtfilt, welch
    from scipy.interpolate import interp1d
    FS = 4096

    def load(path):
        with h5py.File(path, "r") as f:
            return f["strain/Strain"][:].astype(float)

    def bandpass(x, lo=35.0, hi=350.0, order=4):
        bb, aa = butter(order, [lo/(FS/2), hi/(FS/2)], btype="band")
        return filtfilt(bb, aa, x)

    def whiten(x):
        nfft = min(4*FS, len(x))
        fr, psd = welch(x, fs=FS, nperseg=nfft)
        ip = interp1d(fr, psd, bounds_error=False, fill_value=(psd[0], psd[-1]))
        X = np.fft.rfft(x); f = np.fft.rfftfreq(len(x), 1.0/FS)
        asd = np.sqrt(np.maximum(ip(f), 1e-50))
        xw = np.fft.irfft(X/asd, n=len(x))
        return xw/np.std(xw)

    def prep(x):
        return bandpass(whiten(bandpass(x)))

    # (name, H1 path, L1 path, align?) -- GW150914 gets the s10 alignment
    # (L1 leads H1 ~7 ms and inverted) so its s=1 baseline engages the exact
    # s10 noise-on-substrate result; the V1/V2 are independent noise segments.
    segs = [
        ("GW150914_seg", "data/ligo/H-H1_GW150914_32s.hdf5",
         "data/ligo/L-L1_GW150914_32s.hdf5", True),
        ("noise_V1_1167559920", "data/ligo/H-H1_LOSC_4_V1-1167559920-32.hdf5",
         "data/ligo/L-L1_LOSC_4_V1-1167559920-32.hdf5", False),
        ("noise_V2_1135136334", "data/ligo/H-H1_LOSC_4_V2-1135136334-32.hdf5",
         "data/ligo/L-L1_LOSC_4_V2-1135136334-32.hdf5", False),
    ]
    win = int(win_s*FS); stride = int(stride_s*FS); crop = int(2.0*FS)
    out = {}
    for sname, hp, lp, align in segs:
        try:
            h = prep(load(hp))[crop:-crop]
            l = prep(load(lp))
            if align:
                l = -np.roll(l, int(round(0.0069*FS)))
            l = l[crop:-crop]
        except Exception as e:
            out[sname] = {"error": str(e)}
            continue
        out[sname] = {}
        for s in scales:
            M = compute_operator_matrix(h, s*l, win, stride, mi_bins=16)
            if M.shape[0] < 6:
                continue
            rep = null_report(M)
            rep["asym_meanabs"] = per_window_asym(M)
            rep["n_windows"] = int(M.shape[0])
            out[sname][f"{s:g}"] = rep
    return out


# ----------------------------------------------------------------------------
def main():
    t0 = time.time()
    if CANARY:
        scales = [1.0, 3.0, 10.0]
        seeds = [11, 22, 33]
        n = 4000; win = 200; stride = 100; mi_bins = 12
        ligo_scales = [1.0, 3.0, 10.0]
        out_fn = "probe_construction_vs_nature_canary.json"
    else:
        scales = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        seeds = [11, 22, 33, 44, 55]
        n = 20000; win = 400; stride = 100; mi_bins = 16
        ligo_scales = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        out_fn = "probe_construction_vs_nature_results.json"

    print("=" * 96)
    print(f"PROBE construction-vs-nature {'[CANARY]' if CANARY else '[FULL]'}  "
          "(scale one channel; MI-invariant units lever)")
    print("lever: b -> s*b  (MI exactly invariant; H_b -> H_b + ln s)")
    print("=" * 96, flush=True)

    print("\n[SIM] 4 heterogeneous generators x seeds, scale sweep...", flush=True)
    sim = run_sim(scales, seeds, n, win, stride, mi_bins)
    for gname, sd in sim.items():
        print(f"\n  {gname}:")
        print(f"    {'scale':>6} {'ln s':>6} {'asym':>7} {'cos112':>8} "
              f"{'cosIdent':>9} {'meanMI':>8}")
        for sk, c in sd.items():
            print(f"    {c['scale']:>6g} {c['ln_s']:>6.2f} {c['asym_mean']:>7.3f} "
                  f"{c['cos112_mean']:>8.3f} {c['cosident_mean']:>9.3f} "
                  f"{c['meanMI_mean']:>8.4f}", flush=True)

    print("\n[REAL LIGO] noise-only segments, scale L1 post-whitening...",
          flush=True)
    ligo = run_ligo(ligo_scales, win_s=0.125, stride_s=0.03125)
    for sname, sd in ligo.items():
        print(f"\n  {sname}:")
        if "error" in sd:
            print(f"    ERROR: {sd['error']}")
            continue
        print(f"    {'scale':>6} {'asym':>7} {'cos112':>8} {'cosIdent':>9} "
              f"{'meanMI':>8} {'nwin':>5}")
        for sk, c in sd.items():
            print(f"    {float(sk):>6g} {c['asym_meanabs']:>7.3f} "
                  f"{c['cos_to_112']:>8.3f} {c['cos_to_m1m12']:>9.3f} "
                  f"{c['mean_MI']:>8.4f} {c['n_windows']:>5d}", flush=True)

    # ---- scale-(in)variance summary: what moves with a pure units choice? ----
    # MI std ~ 0 across scales = scale-INVARIANT (candidate nature).
    # cos / asym large std across scales = representation-dependent (bookkeeping).
    print("\n" + "=" * 96)
    print("SCALE-(IN)VARIANCE SUMMARY  (across the scale sweep, per system)")
    print("  asym_range : max-min of mean|H_a-H_b| over scales  (units knob)")
    print("  cos112_std : std of cos-to-(1,1,2) over scales      (substrate metric)")
    print("  MI_cv      : coeff of variation of mean MI over scales (~0 => invariant)")
    print("=" * 96)
    summary = {}
    print(f"  {'system':>22} {'asym_range':>11} {'cos112_std':>11} "
          f"{'cosId_std':>10} {'MI_cv':>9}")
    for gname, sd in sim.items():
        asy = [c["asym_mean"] for c in sd.values()]
        c1 = [c["cos112_mean"] for c in sd.values()]
        ci = [c["cosident_mean"] for c in sd.values()]
        mi = [c["meanMI_mean"] for c in sd.values()]
        mi_cv = float(np.std(mi)/ (abs(np.mean(mi))+1e-12))
        summary[gname] = dict(asym_range=float(np.ptp(asy)),
                              cos112_std=float(np.std(c1)),
                              cosident_std=float(np.std(ci)), MI_cv=mi_cv)
        print(f"  {gname:>22} {np.ptp(asy):>11.3f} {np.std(c1):>11.3f} "
              f"{np.std(ci):>10.3f} {mi_cv:>9.2e}", flush=True)
    for sname, sd in ligo.items():
        if "error" in sd:
            continue
        asy = [c["asym_meanabs"] for c in sd.values()]
        c1 = [c["cos_to_112"] for c in sd.values()]
        ci = [c["cos_to_m1m12"] for c in sd.values()]
        mi = [c["mean_MI"] for c in sd.values()]
        mi_cv = float(np.std(mi)/(abs(np.mean(mi))+1e-12))
        summary[sname] = dict(asym_range=float(np.ptp(asy)),
                              cos112_std=float(np.std(c1)),
                              cosident_std=float(np.std(ci)), MI_cv=mi_cv)
        print(f"  {sname:>22} {np.ptp(asy):>11.3f} {np.std(c1):>11.3f} "
              f"{np.std(ci):>10.3f} {mi_cv:>9.2e}", flush=True)

    result = dict(
        meta=dict(canary=CANARY, scales=scales, seeds=seeds,
                  elapsed_s=round(time.time()-t0, 1)),
        sim=sim, ligo=ligo, scale_invariance_summary=summary,
    )
    json.dump(result, open(out_fn, "w"), indent=2)
    print(f"\nWrote {out_fn}  ({result['meta']['elapsed_s']}s)", flush=True)


if __name__ == "__main__":
    main()
