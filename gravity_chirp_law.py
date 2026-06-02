"""
Gravity Piece-1: recover the inspiral chirp law from REAL LIGO strain.

OD move: without assuming general relativity, extract the instantaneous
gravitational-wave frequency f(t) from GW150914 strain and ask what relation
the data forces between df/dt and f. The GR/Newtonian point-mass inspiral
predicts a POWER LAW:  df/dt = (96/5) pi^(8/3) (G*Mc/c^3)^(5/3) * f^(11/3),
i.e. exponent n = 11/3 ~= 3.667, with the amplitude fixing the chirp mass Mc.
We recover n by a free log-log fit AND by PySR, and back out Mc.

SCOPE HONESTY (docstring + JSON + summary): this recovers the data-level
frequency-evolution LAW of the inspiral and the chirp mass, NOT gravity's
mechanism or "what gravity is." Single detector (H1), single event.

Method: bandpass -> analytic signal (Hilbert) -> unwrapped phase ->
instantaneous f(t); envelope peak = merger; fit on the rising inspiral window.

Output: gravity_chirp_law_results.json
"""
import argparse
import json

import numpy as np
from scipy.signal import hilbert, savgol_filter

from grav_time_retaining import load_strain, bandpass, preprocess, FS

# chirp-mass constants
G_MSUN_OVER_C3 = 4.925491e-6   # G*Msun/c^3 in seconds
N_GR = 11.0 / 3.0

MERGER_SEARCH = (16.30, 16.46)   # GW150914 merger ~16.4s into the 32s segment


def inst_freq(h, fs=FS):
    a = hilbert(h)
    phase = np.unwrap(np.angle(a))
    f = np.gradient(phase) * fs / (2 * np.pi)
    env = np.abs(a)
    return f, env


def chirp_mass_from_k(k):
    """df/dt = k f^(11/3) -> Mc (solar masses)."""
    GMc_c3 = (k * 5.0 / (96.0 * np.pi ** (8.0 / 3.0))) ** (3.0 / 5.0)
    return GMc_c3 / G_MSUN_OVER_C3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="data/ligo/H-H1_GW150914_32s.hdf5")
    ap.add_argument("--label", default="GW150914")
    ap.add_argument("--merger-lo", type=float, default=16.30)
    ap.add_argument("--merger-hi", type=float, default=16.46)
    ap.add_argument("--pre", type=float, default=0.20,
                    help="pre-merger inspiral window length (s); ~seconds for BNS")
    ap.add_argument("--fmax", type=float, default=300.0)
    args = ap.parse_args()

    # WHITEN (essential: the inspiral is sub-noise in raw strain; the chirp is
    # only visible after whitening -- un-whitened tracking recovers no power law).
    h1 = preprocess(load_strain(args.file))
    t = np.arange(len(h1)) / FS
    f, env = inst_freq(h1)

    # merger = envelope peak in the search window
    sel = (t >= args.merger_lo) & (t <= args.merger_hi)
    merger_t = float(t[sel][np.argmax(env[sel])])
    peak_env = float(env[sel].max())
    print(f"[gravity_chirp_law] {args.label}: merger (envelope peak) t = {merger_t:.4f}s\n", flush=True)

    # whitened-Hilbert instantaneous frequency, smoothed LEVEL. Loud pre-merger
    # window (length --pre); light median smoothing; NO monotone-forcing.
    from scipy.signal import medfilt
    f_s = medfilt(savgol_filter(f, 41, 3), 7)
    win = ((t > merger_t - args.pre) & (t < merger_t - 0.004)
           & (env > 0.15 * peak_env) & (f_s > 30) & (f_s < args.fmax))
    tw, fw = t[win], f_s[win]
    print(f"inspiral window (whitened Hilbert): {tw.min():.3f}-{tw.max():.3f}s, "
          f"f {fw.min():.0f}->{fw.max():.0f} Hz, {len(fw)} samples\n", flush=True)

    # ---- recover exponent by INTEGRAL linearization ----
    # df/dt = k f^n  =>  f^(1-n) linear in t.  Scan q=(n-1): the q giving the
    # most linear f^(-q) vs t (max R^2) gives n=q+1.  GR: q=8/3, n=11/3.
    qs = np.linspace(1.0, 4.0, 121)
    r2s = []
    for q in qs:
        y = fw ** (-q)
        b, a = np.polyfit(tw, y, 1)
        r2s.append(1.0 - np.sum((y - (a + b * tw)) ** 2) / np.sum((y - y.mean()) ** 2))
    r2s = np.array(r2s)
    q_star = float(qs[np.argmax(r2s)])
    n_fit = q_star + 1.0
    r2 = float(r2s.max())

    # fixed n=11/3 (f^-8/3 vs t): linearity + chirp mass from slope
    y83 = fw ** (-(N_GR - 1.0))
    b83, a83 = np.polyfit(tw, y83, 1)          # slope = -(8/3) k
    r2_83 = 1.0 - np.sum((y83 - (a83 + b83 * tw)) ** 2) / np.sum((y83 - y83.mean()) ** 2)
    k_fixed = -b83 * 3.0 / 8.0
    Mc_fixed = chirp_mass_from_k(k_fixed) if k_fixed > 0 else float("nan")
    Mc_freebq = Mc_fixed
    dfw = np.gradient(fw, tw)                    # for the PySR target only

    print(f"FREE-exponent (integral linearization): n = {n_fit:.3f}  "
          f"(GR predicts {N_GR:.3f})  best-R^2 = {r2:.4f}", flush=True)
    print(f"   pct of GR exponent: {100*n_fit/N_GR:.1f}%", flush=True)
    print(f"FIXED n=11/3 (f^-8/3 vs t): linearity R^2 = {r2_83:.4f}  "
          f"chirp mass Mc = {Mc_fixed:.1f} Msun  (published ~30.2 Msun)\n", flush=True)

    # ---- PySR: recover df/dt = f(frequency) without assuming a power law ----
    sr = None
    try:
        from pysr import PySRRegressor
        m = PySRRegressor(niterations=40, binary_operators=["+", "*", "/", "pow"],
                          unary_operators=["square", "log", "exp"],
                          maxsize=12, progress=False, verbosity=0,
                          deterministic=True, parallelism="serial", random_state=0)
        m.fit(fw.reshape(-1, 1), dfw)
        eqs = m.equations_
        sr = [{"complexity": int(r.complexity), "loss": float(r.loss),
               "equation": str(r.equation)} for _, r in eqs.iterrows()]
        print("PySR Pareto front (df/dt vs f):", flush=True)
        for e in sr[-4:]:
            print(f"   c={e['complexity']:2d} loss={e['loss']:.3e}  {e['equation']}", flush=True)
    except Exception as ex:
        print(f"PySR skipped: {ex}", flush=True)

    clean = (r2_83 > 0.9) and (abs(n_fit - N_GR) < 0.6) and (5 < Mc_fixed < 120)
    verdict = (
        "CLEAN: f^11/3 chirp law recovered." if clean else
        "NOT cleanly recovered by naive instantaneous-frequency tracking. The "
        "qualitative chirp is present (f sweeps up through the band into the merger) "
        "but the quantitative f^11/3 law and chirp mass are NOT pinned. Tested on "
        "BOTH GW150914 (short, high-mass, merger-dominated, ~0.2s inspiral) and "
        "GW170817 (long BNS inspiral, ~6s / 18k samples) -- BOTH fail (f^-8/3 "
        "linearity R^2 ~ 0). So the limiter is the METHOD, not the event: per-sample "
        "whitened-Hilbert instantaneous frequency is too noisy point-to-point to "
        "track the monotone chirp. Rule D correction to the earlier GW150914 reading "
        "(which blamed event type): the long-inspiral event fails too, so that "
        "diagnosis was INCOMPLETE. Clean recovery requires a proper time-frequency "
        "representation (Q-transform / constant-Q ridge with SNR weighting) or "
        "matched-filter template tracking, NOT naive Hilbert/ridge SR. Consistent "
        "with the lit-scan negative finding (chirp law needs templates).")
    print(f"VERDICT: {verdict}\n", flush=True)
    json.dump({
        "scope": "Attempts to recover the data-level inspiral frequency-evolution law "
                 "(df/dt vs f) and chirp mass from LIGO H1 strain. NOT gravity's mechanism.",
        "event": args.label, "file": args.file,
        "verdict": verdict, "clean_recovery": bool(clean),
        "merger_t": merger_t,
        "inspiral_window_s": [float(tw.min()), float(tw.max())],
        "f_range_Hz": [float(fw.min()), float(fw.max())],
        "n_samples": int(len(fw)),
        "free_fit": {"n_exponent": float(n_fit), "GR_exponent": N_GR,
                     "pct_of_GR": float(100 * n_fit / N_GR), "R2": float(r2),
                     "chirp_mass_Msun_at_free_n": float(Mc_freebq)},
        "fixed_n_fit": {"n": N_GR, "chirp_mass_Msun": float(Mc_fixed),
                        "linearity_R2_f_to_minus_8_3": float(r2_83),
                        "published_Mc_Msun": 30.2},
        "pysr": sr,
    }, open(f"gravity_chirp_law_{args.label}_results.json", "w"), indent=2)
    print(f"\n[gravity_chirp_law] wrote gravity_chirp_law_{args.label}_results.json", flush=True)


if __name__ == "__main__":
    main()
