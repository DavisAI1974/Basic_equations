"""
TAUTOLOGY-KILLING NULL for the static algebraic dipole on gravity (LIGO H1/L1).
Follow-up to probe_time_dipole_gravity.py: the algebraic dipole
  H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2
read R2 0.93 (event) but ALSO 0.71 (noise). H_a^2 and H_a*H_b SHARE the factor
H_a, so a high R2 can be a shared-factor tautology, not a genuine joint H_a<->H_b
relation. The shuffle-null (0.045) is too weak -- it destroys autocorrelation.

DECISIVE NULL -- circular-shift of H_b:
  Roll H_b(t) -> H_b(t+s) by a random lag s, recompute H_a*H_b_shifted, refit.
  This PRESERVES H_b's marginal + autocorrelation AND preserves the shared H_a(t)
  factor (so the tautological part survives in BOTH real and null). The only thing
  it destroys is the SPECIFIC instantaneous H_a<->H_b pairing. So:
    real R2 >> shift-null R2  -> genuine H_b-specific joint structure (real dipole)
    real R2 ~= shift-null R2  -> the whole 'constraint' is the shared-H_a tautology
Report event AND noise, each vs its own shift-null. Also report whether the
event EXCESS over its null exceeds the noise excess over its null.
Run:  python probe_time_dipole_gravity_null.py [--canary]
"""
import json
import time
import argparse
import numpy as np

from probe_flow_dipole_gravity import load_strain, preprocess, FS, MERGER_T, \
    H1_PATH, L1_PATH
from kbk_pipeline import compute_operator_matrix


def r2_fit(Ha2, HaHb):
    X = np.column_stack([np.ones_like(HaHb), HaHb, HaHb * HaHb])
    beta, *_ = np.linalg.lstsq(X, Ha2, rcond=None)
    sst = np.sum((Ha2 - Ha2.mean()) ** 2)
    ssr = np.sum((Ha2 - X @ beta) ** 2)
    return (1 - ssr / sst) if sst > 0 else float("nan"), beta


def window_block(h1, l1, lo, hi, win_s=0.125, stride_s=0.03125):
    win = int(win_s * FS); stride = int(stride_s * FS)
    M = compute_operator_matrix(h1, l1, win, stride, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * stride + win / 2) / FS
    sel = (centers >= lo) & (centers < hi)
    Ms = M[sel]
    return Ms[:, 0], Ms[:, 1], Ms[:, 5]   # Ha, Hb, MI


def analyze(Ha, Hb, n_null, seed):
    Ha2 = Ha * Ha
    HaHb = Ha * Hb
    r2_real, beta = r2_fit(Ha2, HaHb)
    rng = np.random.default_rng(seed)
    npts = len(Hb)
    min_shift = max(3, npts // 10)
    nulls = []
    for _ in range(n_null):
        s = rng.integers(min_shift, npts - min_shift)
        Hb_s = np.roll(Hb, s)
        r2n, _ = r2_fit(Ha2, Ha * Hb_s)
        nulls.append(r2n)
    nulls = np.array(nulls)
    # control: H_a^2 explained by H_a alone (deterministic upper bound ~1) and the
    # shared-factor baseline = how well ANY smooth regressor (shifted) does
    return {
        "n_points": int(npts),
        "Ha_std": float(np.std(Ha)), "Hb_std": float(np.std(Hb)),
        "corr_Ha_Hb": float(np.corrcoef(Ha, Hb)[0, 1]),
        "R2_real": float(r2_real),
        "shift_null_mean": float(nulls.mean()),
        "shift_null_p95": float(np.percentile(nulls, 95)),
        "shift_null_std": float(nulls.std()),
        "excess_over_null": float(r2_real - nulls.mean()),
        "z_over_null": float((r2_real - nulls.mean()) / (nulls.std() + 1e-9)),
        "p_null_ge_real": float(np.mean(nulls >= r2_real)),
        "beta": [float(b) for b in beta],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    n_null = 100 if args.canary else 400

    h1 = preprocess(load_strain(H1_PATH))
    l1 = preprocess(load_strain(L1_PATH))

    Ha_e, Hb_e, _ = window_block(h1, l1, MERGER_T - 1.5, MERGER_T + 1.5)
    Ha_n, Hb_n, _ = window_block(h1, l1, 2.0, MERGER_T - 0.5)

    event = analyze(Ha_e, Hb_e, n_null, 11)
    noise = analyze(Ha_n, Hb_n, n_null, 22)

    # decisive comparison: does the EVENT carry structure beyond its tautology null
    # MORE than the noise does beyond its own?
    event_genuine = bool(event["p_null_ge_real"] < 0.05 and
                         event["excess_over_null"] > 0.10)
    noise_genuine = bool(noise["p_null_ge_real"] < 0.05 and
                         noise["excess_over_null"] > 0.10)

    result = {
        "probe": "tautology-killing null (circular-shift of H_b) for the gravity "
                 "static algebraic dipole H_a^2 = a + b*H_aH_b + c*(H_aH_b)^2",
        "logic": "shift preserves H_b smoothness + shared H_a factor; kills only the "
                 "instantaneous H_a<->H_b pairing. real >> shift-null => genuine.",
        "event": event,
        "noise": noise,
        "verdict": {
            "event_beats_tautology_null": event_genuine,
            "noise_beats_tautology_null": noise_genuine,
            "event_excess_minus_noise_excess": round(
                event["excess_over_null"] - noise["excess_over_null"], 3),
            "reading": "if event_genuine and excess>>noise excess -> the algebraic "
                       "dipole carries real event-specific joint structure (not just "
                       "the shared-H_a tautology). if both ~0 excess -> tautology.",
        },
        "runtime_s": round(time.time() - t0, 1),
    }
    out = ("probe_time_dipole_gravity_null_canary.json" if args.canary
           else "probe_time_dipole_gravity_null_results.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)

    print(f"=== tautology-killing null (gravity algebraic dipole) "
          f"{'CANARY' if args.canary else 'FULL'} {result['runtime_s']}s ===")
    for nm, blk in [("EVENT", event), ("NOISE", noise)]:
        print(f"\n[{nm}]  n={blk['n_points']}  Ha_std={blk['Ha_std']:.3f} "
              f"Hb_std={blk['Hb_std']:.3f}  corr(Ha,Hb)={blk['corr_Ha_Hb']:+.3f}")
        print(f"   R2_real        = {blk['R2_real']:.3f}")
        print(f"   shift-null     = {blk['shift_null_mean']:.3f} "
              f"(p95 {blk['shift_null_p95']:.3f}, std {blk['shift_null_std']:.3f})")
        print(f"   excess         = {blk['excess_over_null']:+.3f}  "
              f"z={blk['z_over_null']:.1f}  p(null>=real)={blk['p_null_ge_real']:.3f}")
    v = result["verdict"]
    print(f"\nVERDICT: event_genuine={v['event_beats_tautology_null']}  "
          f"noise_genuine={v['noise_beats_tautology_null']}  "
          f"event_excess - noise_excess = {v['event_excess_minus_noise_excess']:+.3f}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
