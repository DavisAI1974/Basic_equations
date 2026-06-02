"""
THE OTHER DIPOLE, TO SOLVE FOR TIME -- gravity (real LIGO H1/L1, GW150914).
Greg (S19): the flow/entropy dipole came back blind on the forces, but our TIME
findings were POSITIVE -- don't throw them away; "time is something else"; try one
of the OTHER dipole forms to solve for time.

WHY: the EM tool-battery diagnosed exactly why the windowed-entropy->MI->dMI/dt
(flow) dipole is blind -- the per-channel marginal entropies H_a,H_b are
(here) independent of the inter-channel TIME LAG, so the coupling/time information
never enters the entropy operators. It lives in the RAW covariance. The raw-count
covariance dipole HIT (R2 0.49) where the entropy dipole was blind. And the repo's
static_dipole_test.py already noted the same pattern (DNA: R2=0 differential,
preserved as an ALGEBRAIC constraint). So we try the OTHER dipole forms:

  (1) RAW cross-covariance dipole vs lag  -> SOLVE FOR TIME: the inter-detector
      light-travel delay. GW150914 reached L1 ~6.9 ms before H1 (a real, known
      gravity TIME quantity, INFO-053, positive). Does the raw dipole recover it,
      and is the ENTROPY-MI version blind to it (lag-flat)?  This is the head-to-
      head: same data, two dipole forms, the time answer we already know.
  (2) STATIC ALGEBRAIC dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2 on the
      windowed operators -- does the algebraic constraint hold (R2) where the
      DIFFERENTIAL flow dipole was blind (flow R2 was 0.022)?  (DNA precedent.)

Deflationary guards (frame must not grade itself):
  - lag null: time-slide / off-source windows -> the recovered lag must beat a
    scramble, not just exist.
  - the raw cross-correlation recovering the light-travel time is STANDARD matched-
    interferometry, NOT new physics -- it is the positive-control that proves the
    raw dipole sees what the entropy dipole cannot. Reported as such.
  - algebraic R2 compared against a shuffle/AR null so a high R2 from smoothness is
    not mistaken for a constraint.

Reuses the verbatim s10 LIGO loader/whitening from probe_flow_dipole_gravity.py.
Run:  python probe_time_dipole_gravity.py [--canary]
"""
import json
import time
import argparse
import numpy as np

from probe_flow_dipole_gravity import (
    load_strain, preprocess, FS, MERGER_T, H1_PATH, L1_PATH,
)
from kbk_pipeline import compute_operator_matrix


# ----------------------------------------------------------------------------
# (1) RAW cross-covariance dipole vs lag  ->  solve for the inter-detector TIME
# ----------------------------------------------------------------------------
def solve_lag_raw(h1, l1, t_center, half_win_s, max_lag_ms):
    """Normalized cross-correlation of the two RAW whitened strains in a window
    around t_center, scanned over lag. Returns the lag (ms) that maximizes
    coupling = the solved inter-detector light-travel time."""
    i0 = int((t_center - half_win_s) * FS)
    i1 = int((t_center + half_win_s) * FS)
    a = h1[i0:i1].copy()
    b = l1[i0:i1].copy()
    a = (a - a.mean()) / (a.std() + 1e-30)
    b = (b - b.mean()) / (b.std() + 1e-30)
    max_lag = int(max_lag_ms * 1e-3 * FS)
    lags = np.arange(-max_lag, max_lag + 1)
    n = len(a)
    cc = np.empty(len(lags))
    for k, lag in enumerate(lags):
        if lag >= 0:
            x, y = a[lag:], b[: n - lag]
        else:
            x, y = a[: n + lag], b[-lag:]
        cc[k] = np.dot(x, y) / len(x)
    # L1 is sign-inverted vs H1 for GW150914; the physical coupling is |cc| peak
    kbest = int(np.argmax(np.abs(cc)))
    lag_ms = lags[kbest] / FS * 1e3
    return lag_ms, float(cc[kbest]), lags / FS * 1e3, cc


def solve_lag_entropy(h1, l1, t_center, half_win_s, max_lag_ms, mwin_s=0.01):
    """Same window, but does the ENTROPY-MI dipole carry the lag? Slide L1 vs H1,
    at each lag compute windowed (H_a,H_b,MI) operator matrix and use MI-mean as
    the coupling score. If entropy/MI is lag-flat -> the entropy dipole is blind
    to the TIME structure (the EM diagnosis)."""
    i0 = int((t_center - half_win_s) * FS)
    i1 = int((t_center + half_win_s) * FS)
    a = h1[i0:i1].copy()
    b = l1[i0:i1].copy()
    max_lag = int(max_lag_ms * 1e-3 * FS)
    lags = np.arange(-max_lag, max_lag + 1, max(1, int(0.001 * FS)))  # ~1ms grid
    n = len(a)
    mwin = max(8, int(mwin_s * FS))
    stride = mwin // 2
    mi_score = []
    for lag in lags:
        if lag >= 0:
            x, y = a[lag:], b[: n - lag]
        else:
            x, y = a[: n + lag], b[-lag:]
        M = compute_operator_matrix(x, y, mwin, stride, mi_bins=12)
        mi_score.append(float(np.mean(M[:, 5])))  # MI column
    mi_score = np.array(mi_score)
    lag_ms = lags / FS * 1e3
    kbest = int(np.argmax(mi_score))
    # flatness: how much does MI vary across lag vs its level?
    flatness = float(np.std(mi_score) / (np.mean(mi_score) + 1e-12))
    return lag_ms[kbest], lag_ms, mi_score, flatness


def lag_null(h1, l1, half_win_s, max_lag_ms, n_off, rng):
    """Off-source windows: recover a lag where there is no event -> null
    distribution of |best cc| and of the recovered lag scatter."""
    dur = len(h1) / FS
    best_abscc, best_lag = [], []
    for _ in range(n_off):
        tc = rng.uniform(half_win_s + 1.0, MERGER_T - 2.0)
        lag_ms, cc, _, _ = solve_lag_raw(h1, l1, tc, half_win_s, max_lag_ms)
        best_abscc.append(abs(cc))
        best_lag.append(lag_ms)
    return (float(np.mean(best_abscc)), float(np.std(best_abscc)),
            float(np.std(best_lag)))


# ----------------------------------------------------------------------------
# (2) STATIC ALGEBRAIC dipole on the windowed operators
# ----------------------------------------------------------------------------
def algebraic_dipole(h1, l1, win_s, stride_s, mask_lo, mask_hi):
    win = int(win_s * FS)
    stride = int(stride_s * FS)
    M = compute_operator_matrix(h1, l1, win, stride, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * stride + win / 2) / FS
    sel = (centers >= mask_lo) & (centers < mask_hi)
    Msel = M[sel]
    Ha, Hb, MI = Msel[:, 0], Msel[:, 1], Msel[:, 5]
    HaHb = Ha * Hb
    y = Ha * Ha
    X = np.column_stack([np.ones_like(HaHb), HaHb, HaHb * HaHb])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    ssr = float(np.sum((y - yhat) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ssr / sst if sst > 0 else float("nan")
    # shuffle null
    rng = np.random.default_rng(0)
    null = []
    for _ in range(50):
        p = rng.permutation(len(y))
        Xn = np.column_stack([np.ones_like(HaHb), HaHb[p], (HaHb * HaHb)[p]])
        bn, *_ = np.linalg.lstsq(Xn, y, rcond=None)
        rn = 1 - np.sum((y - Xn @ bn) ** 2) / sst
        null.append(rn)
    return float(r2), float(np.mean(null)), float(np.percentile(null, 95)), \
        beta.tolist(), int(sel.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    h1 = preprocess(load_strain(H1_PATH))
    l1 = preprocess(load_strain(L1_PATH))

    half_win_s = 0.10 if args.canary else 0.15
    max_lag_ms = 15.0
    n_off = 8 if args.canary else 30

    # (1a) RAW dipole solves for the lag at the merger
    lag_raw, cc_raw, lag_axis, cc_curve = solve_lag_raw(
        h1, l1, MERGER_T, half_win_s, max_lag_ms)
    # (1b) ENTROPY-MI dipole: does it carry the lag?
    lag_ent, ent_axis, ent_curve, ent_flatness = solve_lag_entropy(
        h1, l1, MERGER_T, half_win_s, max_lag_ms)
    # (1c) null
    rng = np.random.default_rng(7)
    cc_null_mean, cc_null_std, lag_null_scatter = lag_null(
        h1, l1, half_win_s, max_lag_ms, n_off, rng)

    # (2) static algebraic dipole on event window
    alg_r2_event, alg_null_mean, alg_null_p95, alg_beta, n_ev = algebraic_dipole(
        h1, l1, 0.125, 0.03125, MERGER_T - 1.5, MERGER_T + 1.5)
    alg_r2_noise, _, _, _, n_no = algebraic_dipole(
        h1, l1, 0.125, 0.03125, 2.0, MERGER_T - 0.5)

    PHYS_LAG_MS = 6.9   # GW150914 known inter-detector light-travel delay (L before H)
    raw_solves_time = bool(abs(abs(lag_raw) - PHYS_LAG_MS) <= 3.0 and
                           abs(cc_raw) > cc_null_mean + 3 * cc_null_std)

    result = {
        "probe": "the OTHER dipole, to solve for TIME -- gravity LIGO H1/L1 GW150914",
        "idea": "flow/entropy dipole blind (entropy is lag-independent); try raw & "
                "algebraic dipole; SOLVE FOR TIME = recover the inter-detector "
                "light-travel delay; do not discard the positive time findings.",
        "config": {"half_win_s": half_win_s, "max_lag_ms": max_lag_ms,
                   "n_off_null": n_off, "known_phys_lag_ms": PHYS_LAG_MS},
        "raw_covariance_dipole": {
            "solved_lag_ms": round(lag_raw, 2),
            "peak_abs_cc": round(abs(cc_raw), 4),
            "null_abs_cc_mean": round(cc_null_mean, 4),
            "null_abs_cc_std": round(cc_null_std, 4),
            "off_source_lag_scatter_ms": round(lag_null_scatter, 2),
            "z_over_null": round((abs(cc_raw) - cc_null_mean) /
                                 (cc_null_std + 1e-12), 1),
            "SOLVES_FOR_TIME": raw_solves_time,
        },
        "entropy_MI_dipole_same_data": {
            "best_lag_ms": round(lag_ent, 2),
            "MI_vs_lag_flatness": round(ent_flatness, 4),
            "blind_to_lag": bool(ent_flatness < 0.05),
            "note": "if MI is lag-flat the entropy dipole carries no TIME info -- "
                    "the EM-battery diagnosis, here on gravity.",
        },
        "static_algebraic_dipole": {
            "form": "H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2",
            "R2_event": round(alg_r2_event, 3),
            "R2_noise": round(alg_r2_noise, 3),
            "shuffle_null_mean": round(alg_null_mean, 3),
            "shuffle_null_p95": round(alg_null_p95, 3),
            "beta_event": [round(b, 4) for b in alg_beta],
            "flow_dipole_R2_for_reference": 0.022,
            "holds_where_flow_blind": bool(alg_r2_event > alg_null_p95 and
                                           alg_r2_event > 0.3),
        },
        "curves": {
            "raw_cc_lag_axis_ms": [round(x, 2) for x in lag_axis.tolist()[::4]],
            "raw_cc_values": [round(x, 4) for x in cc_curve.tolist()[::4]],
            "entropy_MI_lag_axis_ms": [round(x, 2) for x in ent_axis.tolist()],
            "entropy_MI_values": [round(x, 4) for x in ent_curve.tolist()],
        },
        "runtime_s": round(time.time() - t0, 1),
    }
    out = ("probe_time_dipole_gravity_canary.json" if args.canary
           else "probe_time_dipole_gravity_results.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)

    r = result
    print(f"=== OTHER dipole / solve-for-time: gravity ({'CANARY' if args.canary else 'FULL'}) "
          f"{r['runtime_s']}s ===")
    rd = r["raw_covariance_dipole"]
    print(f"\n[1] RAW covariance dipole -> SOLVE FOR TIME (inter-detector lag):")
    print(f"    solved lag = {rd['solved_lag_ms']:+.2f} ms   (known physical ~ -6.9 ms)")
    print(f"    peak|cc| = {rd['peak_abs_cc']}  vs off-source null "
          f"{rd['null_abs_cc_mean']}+/-{rd['null_abs_cc_std']}  (z={rd['z_over_null']})")
    print(f"    off-source lag scatter = {rd['off_source_lag_scatter_ms']} ms")
    print(f"    SOLVES FOR TIME: {rd['SOLVES_FOR_TIME']}")
    ed = r["entropy_MI_dipole_same_data"]
    print(f"\n[1b] ENTROPY-MI dipole on the SAME data:")
    print(f"    MI-vs-lag flatness = {ed['MI_vs_lag_flatness']}  "
          f"blind_to_lag={ed['blind_to_lag']}  (best lag {ed['best_lag_ms']} ms, spurious)")
    ad = r["static_algebraic_dipole"]
    print(f"\n[2] STATIC ALGEBRAIC dipole  H_a^2 = a + b*H_aH_b + c*(H_aH_b)^2:")
    print(f"    R2 event = {ad['R2_event']}  noise = {ad['R2_noise']}  "
          f"shuffle-null p95 = {ad['shuffle_null_p95']}")
    print(f"    (flow dipole R2 was {ad['flow_dipole_R2_for_reference']})  "
          f"holds_where_flow_blind = {ad['holds_where_flow_blind']}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
