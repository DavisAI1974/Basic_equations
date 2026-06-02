"""
FLOW DIPOLE EQUATION -- GRAVITY via the two LIGO detectors (H1, L1).

Two-channel system: channel A = H1 whitened strain, channel B = L1 whitened
strain. Per sliding window (~125 ms) we compute the SAME operators used
everywhere in this repo (kbk_pipeline.compute_operator_matrix: windowed
Vasicek differential entropy H_a, H_b, and histogram MI I(H1;L1)). The
paper's flow-dipole form (https://davisai.ai/dipole/ Sec 2.2):

    dMI/dt ~ sum_i c_self,i*H_i^2 + sum_{i<j} c_cross,ij*H_i*H_j
             + linear terms + const
    (opposition signature: self-terms (H_a^2,H_b^2) opposite sign to cross H_a*H_b)
    algebraic ratio C = H_self/H_cross = mean(H_a,H_b)/MI

Mirrors probe_flow_dipole_brusselator.py for the regression, the DEFLATIONARY
controls, and the JSON schema. Reuses s10_ligo_extract.py loader / whitening
verbatim (so operators are consistent with INFO-036/053).

DEFLATIONARY CONTROLS (the frame must not grade itself):
  (B) 1-D MI-relaxation: dMI/dt ~ a*MI + b. If it matches the full model the
      "dipole" is just relaxation.
  (C) self-only (drop H_a*H_b): dR2_cross = R2_full - R2_self_only is the
      marginal value of the coupling term.
  (S) feature time-shuffle null: permute operator rows, keep dMI/dt fixed,
      refit 20x. R^2 should collapse to ~0 if the fit is real.

SCATTER (Result Discipline, one real event): we generate scatter by a
window-START-JITTER sweep (5 offsets) on GW150914 H1/L1. Each offset is one
"realization"; we report inter-realization scatter on standardized betas / R^2.
There is no second 2-detector event available locally (only H1 exists for
GW170817), so we DO NOT fabricate a second event -- noted as a caveat.

GUARD (merger transient may dominate, cf INFO-036): we fit the dipole on
(i) the FULL 32 s segment AND (ii) a NOISE-ONLY pre-merger window separately
and report both.

Run:  python probe_flow_dipole_gravity.py --canary   (1 offset, coarse stride)
      python probe_flow_dipole_gravity.py            (full: 5 offsets)
"""
import sys
import json
import time
import argparse
import numpy as np
import h5py
from scipy.signal import butter, filtfilt, welch
from scipy.interpolate import interp1d

from kbk_pipeline import compute_operator_matrix

# operator column order from compute_operator_matrix:
# [H_a, H_b, H_a^2, H_b^2, H_a*H_b, MI]
COL = dict(Ha=0, Hb=1, Ha2=2, Hb2=3, HaHb=4, MI=5)
FEAT_NAMES = ["H_a^2", "H_b^2", "H_a*H_b", "H_a", "H_b", "1"]

FS = 4096
GPS_START = 1126259446
MERGER_GPS = 1126259462.4
MERGER_T = MERGER_GPS - GPS_START          # ~16.4 s into segment
EDGE_S = 2.0                               # crop filter-corrupted edges

H1_PATH = "data/ligo/H-H1_GW150914_32s.hdf5"
L1_PATH = "data/ligo/L-L1_GW150914_32s.hdf5"


# ---------------- LIGO loader / preprocessing (verbatim from s10) ----------
def load_strain(path):
    with h5py.File(path, "r") as f:
        return f["strain/Strain"][:].astype(float)


def bandpass(x, lo=35.0, hi=350.0, fs=FS, order=4):
    b, a = butter(order, [lo / (fs / 2), hi / (fs / 2)], btype="band")
    return filtfilt(b, a, x)


def whiten(x, fs=FS):
    n = len(x)
    freqs, psd = welch(x, fs=fs, nperseg=min(4 * fs, n))
    interp_psd = interp1d(freqs, psd, bounds_error=False,
                          fill_value=(psd[0], psd[-1]))
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    asd = np.sqrt(np.maximum(interp_psd(f), 1e-50))
    xw = np.fft.irfft(X / asd, n=n)
    return xw / np.std(xw)


def preprocess(x):
    return bandpass(whiten(bandpass(x)))


# ---------------- regression helpers (mirror brusselator template) ----------
def moving_average(v, w):
    if w <= 1:
        return v.copy()
    k = np.ones(w) / w
    return np.convolve(v, k, mode="same")


def ols(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2


def standardized_betas(Xcols, y):
    Xs = np.column_stack([(c - c.mean()) / (c.std() + 1e-12) for c in Xcols])
    ys = (y - y.mean()) / (y.std() + 1e-12)
    beta, r2 = ols(Xs, ys)
    return beta, r2


def fit_block(Ha, Hb, MI, t, smooth_w, seed):
    """Run the full flow-dipole fit + controls on one block of windows.

    Ha,Hb,MI,t are 1-D arrays over windows (already trimmed to the block).
    Returns the per-block result dict, mirroring the brusselator schema.
    """
    MIs = moving_average(MI, smooth_w)
    dMI = np.gradient(MIs, t)

    # trim convolution edges
    edge = smooth_w
    if len(t) <= 2 * edge + 6:
        return None
    sl = slice(edge, len(t) - edge)
    Ha, Hb, MI, dMI = Ha[sl], Hb[sl], MI[sl], dMI[sl]
    Ha2, Hb2, HaHb = Ha * Ha, Hb * Hb, Ha * Hb

    # ---- Model A: full flow-dipole ----
    XA = np.column_stack([Ha2, Hb2, HaHb, Ha, Hb, np.ones_like(Ha)])
    betaA, r2A = ols(XA, dMI)
    sbeta, _ = standardized_betas([Ha2, Hb2, HaHb, Ha, Hb], dMI)

    # ---- Model B: 1-D MI-relaxation ----
    XB = np.column_stack([MI, np.ones_like(MI)])
    betaB, r2B = ols(XB, dMI)

    # ---- Model C: self-only (no cross term) ----
    XC = np.column_stack([Ha2, Hb2, Ha, Hb, np.ones_like(Ha)])
    _, r2C = ols(XC, dMI)

    dR2_dipole_over_relax = r2A - r2B
    dR2_cross = r2A - r2C

    # ---- feature time-shuffle null ----
    rng = np.random.default_rng(seed + 999)
    r2_null = []
    for _ in range(20):
        perm = rng.permutation(len(dMI))
        Xn = np.column_stack([Ha2[perm], Hb2[perm], HaHb[perm],
                              Ha[perm], Hb[perm], np.ones_like(Ha)])
        _, r2n = ols(Xn, dMI)
        r2_null.append(r2n)
    r2_null = np.array(r2_null)

    # opposition signature in the quadratic subspace
    c_self = np.array([betaA[0], betaA[1]])
    c_cross = betaA[2]
    opposition = bool(np.sign(c_cross) != 0 and
                      np.all(np.sign(c_self) == -np.sign(c_cross)))

    # algebraic ratio C = mean(H_a,H_b)/MI (guard small MI)
    meanH = 0.5 * (Ha + Hb)
    C_ratio = float(np.mean(meanH) / (np.mean(MI) + 1e-12))

    return {
        "n_points": int(len(dMI)),
        "MI_range": [float(MI.min()), float(MI.max())],
        "dMI_abs_mean": float(np.mean(np.abs(dMI))),
        "model_A_full": {
            "coeffs_raw": {FEAT_NAMES[i]: float(betaA[i]) for i in range(6)},
            "coeffs_std": {FEAT_NAMES[i]: float(sbeta[i]) for i in range(5)},
            "R2": float(r2A),
        },
        "model_B_MI_relax": {"a": float(betaB[0]), "b": float(betaB[1]),
                             "R2": float(r2B)},
        "model_C_self_only_R2": float(r2C),
        "dR2_dipole_over_relax": float(dR2_dipole_over_relax),
        "dR2_cross_term": float(dR2_cross),
        "shuffle_null_R2_mean": float(r2_null.mean()),
        "shuffle_null_R2_max": float(r2_null.max()),
        "opposition_signature_quad": opposition,
        "C_ratio_mean_over_MI": C_ratio,
    }


def analyze_offset(h1c, l1c, centers, win_s, stride_s, offset_windows,
                   smooth_w, seed):
    """One 'realization' = window-start jitter by offset_windows windows.

    We jitter by dropping the first `offset_windows` rows of the operator
    matrix (equivalent to starting the window grid offset_windows*stride
    later). Returns FULL-segment and NOISE-ONLY block fits.
    """
    n = len(centers)
    sl = slice(offset_windows, n)
    Ha = h1c[0][sl]
    Hb = h1c[1][sl]
    MI = h1c[2][sl]
    cen = centers[sl]
    t = cen  # window center time in seconds

    # FULL segment
    full = fit_block(Ha.copy(), Hb.copy(), MI.copy(), t.copy(),
                     smooth_w, seed)

    # NOISE-ONLY: pre-merger, ending >=0.5 s before the merger
    noise_mask = cen < (MERGER_T - 0.5)
    if noise_mask.sum() > 2 * smooth_w + 6:
        noise = fit_block(Ha[noise_mask].copy(), Hb[noise_mask].copy(),
                          MI[noise_mask].copy(), t[noise_mask].copy(),
                          smooth_w, seed)
    else:
        noise = None

    # EVENT-FOCUSED: +/-1.5 s around the merger (the regime INFO-036 found
    # structure; small-N so read with care)
    event_mask = np.abs(cen - MERGER_T) <= 1.5
    if event_mask.sum() > 2 * smooth_w + 6:
        event = fit_block(Ha[event_mask].copy(), Hb[event_mask].copy(),
                          MI[event_mask].copy(), t[event_mask].copy(),
                          smooth_w, seed)
    else:
        event = None

    return {"offset_windows": int(offset_windows),
            "n_windows_full": int(len(cen)),
            "full_segment": full,
            "noise_only": noise,
            "event_window": event}


def build_operators(win_s, stride_s):
    """Load + preprocess GW150914 H1/L1, build windowed operator series.

    Returns (Ha, Hb, MI) arrays over windows and window-center times (s).
    """
    h1 = load_strain(H1_PATH)
    l1 = load_strain(L1_PATH)
    h1p = preprocess(h1)
    l1p = preprocess(l1)

    # GW150914: L1 leads H1 by ~7 ms and is inverted. Align (same as s10).
    shift = int(round(0.0069 * FS))
    l1a = -np.roll(l1p, shift)

    crop = int(EDGE_S * FS)
    h1c = h1p[crop:-crop]
    l1c = l1a[crop:-crop]

    win = int(win_s * FS)
    stride = int(stride_s * FS)
    M = compute_operator_matrix(h1c, l1c, win, stride, mi_bins=16)
    n_w = M.shape[0]
    centers = (np.arange(n_w) * stride + win / 2) / FS + EDGE_S
    return M[:, COL["Ha"]], M[:, COL["Hb"]], M[:, COL["MI"]], centers


def aggregate_betas(block_results):
    """mean +/- std of standardized betas and key R^2 across realizations."""
    blocks = [b for b in block_results if b is not None]
    if not blocks:
        return {}
    agg = {}
    for k in FEAT_NAMES[:5]:
        vals = [b["model_A_full"]["coeffs_std"][k] for b in blocks]
        agg[k] = {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
    vals = [b["model_A_full"]["R2"] for b in blocks]
    agg["full_R2"] = {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
    for fld in ["model_C_self_only_R2", "dR2_dipole_over_relax",
                "dR2_cross_term", "shuffle_null_R2_mean", "C_ratio_mean_over_MI"]:
        vals = [b[fld] for b in blocks]
        agg[fld] = {"mean": float(np.mean(vals)), "std": float(np.std(vals))}
    relax = [b["model_B_MI_relax"]["R2"] for b in blocks]
    agg["relax_R2"] = {"mean": float(np.mean(relax)), "std": float(np.std(relax))}
    opp = [b["opposition_signature_quad"] for b in blocks]
    agg["opposition_fraction"] = float(np.mean(opp))
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    win_s = 0.125
    smooth_w = 5
    if args.canary:
        stride_s = 0.0625          # 62.5 ms -> ~448 windows over 28 s
        offsets = [0]
        out_path = "probe_flow_dipole_gravity_canary.json"
    else:
        stride_s = 0.03125         # 31 ms -> ~896 windows
        offsets = [0, 1, 2, 3, 4]  # window-start jitter sweep (5 realizations)
        out_path = "probe_flow_dipole_gravity_results.json"

    print("=" * 92)
    print(f"FLOW-DIPOLE GRAVITY (GW150914 H1/L1) {'[CANARY]' if args.canary else '[FULL]'}")
    print("=" * 92, flush=True)

    Ha, Hb, MI, centers = build_operators(win_s, stride_s)
    print(f"operators: {len(centers)} windows of {win_s*1000:.0f}ms "
          f"(stride {stride_s*1000:.0f}ms), merger at {MERGER_T:.1f}s", flush=True)

    realizations = []
    for off in offsets:
        seed = 11 + off
        r = analyze_offset((Ha, Hb, MI), None, centers, win_s, stride_s,
                           off, smooth_w, seed)
        realizations.append(r)

    full_blocks = [r["full_segment"] for r in realizations]
    noise_blocks = [r["noise_only"] for r in realizations]
    event_blocks = [r["event_window"] for r in realizations]
    agg_full = aggregate_betas(full_blocks)
    agg_noise = aggregate_betas(noise_blocks)
    agg_event = aggregate_betas(event_blocks)

    result = {
        "probe": "flow-dipole equation -- GRAVITY (LIGO H1/L1, GW150914)",
        "form": "dMI/dt ~ c_self*H_i^2 + c_cross*H_a*H_b + linear + const",
        "channels": {"A": "H1 whitened strain", "B": "L1 whitened strain (lag/inv aligned)"},
        "config": {"event": "GW150914", "fs": FS, "win_s": win_s,
                   "stride_s": stride_s, "edge_crop_s": EDGE_S,
                   "smooth_w": smooth_w, "offsets": offsets,
                   "merger_t_s": MERGER_T},
        "scatter_method": "window-start jitter sweep (5 offsets); one 2-detector "
                          "event only -- GW170817 has H1 but no L1 locally",
        "per_realization": realizations,
        "aggregate_full_segment": agg_full,
        "aggregate_noise_only": agg_noise,
        "aggregate_event_window": agg_event,
        "reading_guide": {
            "is_a_real_dipole_if": "dR2_cross_term > 0 (H_a*H_b adds signal) AND "
                                   "dR2_dipole_over_relax > 0 (beats 1-D MI relaxation) "
                                   "AND shuffle_null_R2 ~ 0 AND opposition holds",
            "deflationary_if": "relax_R2 ~= full R2 (just 1-D MI relaxation) OR "
                               "shuffle_null_R2 high (column-count artifact) OR "
                               "dR2_cross ~ 0 (cross term carries nothing)",
            "guard_full_vs_noise": "compare aggregate_full_segment vs "
                                   "aggregate_noise_only: the merger transient may "
                                   "drive the full-segment fit (INFO-036: noise sits "
                                   "on the equal-entropy attractor, the chirp leaves it)",
            "scale_caveat": "raw coeffs units-dependent (INFO-051: MI is the only "
                            "scale-invariant operator); read standardized betas",
            "data_caveat": "one real event, one segment; whitening sets the entropy "
                           "scale; short lever arm / finite SNR; no novelty claim "
                           "(literature/novelty scan pending)",
        },
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    def show(tag, agg):
        if not agg:
            print(f"  [{tag}] no valid block")
            return
        print(f"  [{tag}] full-R2={agg['full_R2']['mean']:.3f}+/-{agg['full_R2']['std']:.3f}  "
              f"relax-R2={agg['relax_R2']['mean']:.3f}  "
              f"self-only-R2={agg['model_C_self_only_R2']['mean']:.3f}")
        print(f"        dR2 dipole>relax={agg['dR2_dipole_over_relax']['mean']:+.3f}  "
              f"dR2 cross={agg['dR2_cross_term']['mean']:+.3f}  "
              f"shuffle-null={agg['shuffle_null_R2_mean']['mean']:.3f}  "
              f"opp={agg['opposition_fraction']:.2f}  "
              f"C={agg['C_ratio_mean_over_MI']['mean']:.2f}")
        print("        std betas:", "  ".join(
            f"{k}={agg[k]['mean']:+.3f}+/-{agg[k]['std']:.3f}" for k in FEAT_NAMES[:5]))

    print(f"seeds/offsets={offsets}  runtime={result['runtime_s']}s")
    print("FULL SEGMENT:")
    show("full", agg_full)
    print("NOISE-ONLY (pre-merger):")
    show("noise", agg_noise)
    print("EVENT-WINDOW (+/-1.5s of merger):")
    show("event", agg_event)
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
