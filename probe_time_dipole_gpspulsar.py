"""
THE OTHER DIPOLE, TO SOLVE FOR TIME -- carried to the S18 GPS time-dilation
recovery (and the decision-gate call on the pulsar recovery).

CONTEXT (S19). The windowed-entropy "flow dipole" (dMI/dt ~ ...) came back BLIND
on real force data, but two OTHER dipole forms travelled to real gravity (LIGO
H1/L1) and solved for TIME:
  (A) RAW cross-covariance dipole -> recovered the inter-detector light-travel
      TIME (~7 ms, z~14).
  (B) STATIC ALGEBRAIC dipole  H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2  -> held
      R2 0.93 at the event AND survived a circular-shift tautology-killing null
      (event excess > 0; noise excess ~0 = pure tautology).
See probe_time_dipole_gravity.py / probe_time_dipole_gravity_null.py.

GOAL: do the OTHER dipole forms also capture the GR time-coupling at our OTHER
positive TIME findings -- the S18 gravity-couples-to-time recoveries?
  - GPS time dilation (INFO-059): recovered -2/c^2 from raw RINEX pseudorange
    (BRUX) vs an INDEPENDENT broadcast-element regressor e*sqrt(a)*sin(E).
  - Pulsar (INFO-060): recovered dP_b/dt + Einstein-delay gamma from raw Arecibo
    TOAs of PSR B1913+16 via PINT (a global WLS PARAMETER fit, not a 2-channel
    time series).

DECISION GATE (Result Discipline -- stated, then acted on):
  GPS = GENUINE 2-channel object. The S18 recovery already builds two
  SEPARATELY-SOURCED, co-sampled time series over the orbit:
    channel A = cleaned pseudorange residual  ye(t)  = the measured
                relativistic clock-residual (RINEX pseudorange - geometry -
                precise SP3 clock - receiver clock), in metres.
    channel B = the orbital-geometry driver  c*F*e*sqrt(a)*sin(E(t))  computed
                from BROADCAST Keplerian elements (a different data source),
                in metres (so A,B are commensurate time-equivalent signals).
  Their coupling over the orbit IS the time-dilation law. The S18 probe did a
  single-lag linear REGRESSION of A on B (recovered coefficient k). The dipole
  question is genuinely different: (a) the RAW cross-covariance over LAG is a
  different operator (it asks WHEN, not just how much), and (b) the STATIC
  ALGEBRAIC dipole operates on the windowed ENTROPIES of A and B, which the S18
  probe never computed. So GPS is a legitimate 2-channel dipole object, not the
  regression relabeled -> PROCEED.

  PULSAR = CATEGORY STRETCH (written up, NOT fabricated -- like the chemistry
  call). dP_b/dt and gamma are SCALAR parameters of a global WLS timing-model
  fit; there is no co-sampled SECOND PHYSICAL channel. dP_b/dt is a SECULAR
  cumulative-parabola effect, not an instantaneous channel coupling; any
  "orbital-phase driver" channel would be synthesized from the SAME model the
  residuals were reduced against. Forcing a 2-channel entropy/covariance dipole
  onto a 1-parameter secular fit is the chemistry-into-the-4-forces mismatch.
  Recorded as a clean negative in pulsar_decision_gate below; no pulsar dipole
  is fabricated.

HONESTY: NO novelty / new-law claims. The S18 GPS recovery is a RECOVERY of a
KNOWN GR law (-2/c^2); if the OTHER dipole forms capture that coupling it is a
POSITIVE-CONTROL demonstration that the dipole tool travels, NOT a discovery.
Each result is one data point; the tautology-null excess is reported honestly.

Run:  python probe_time_dipole_gpspulsar.py [--canary]
Mirrors the JSON style of probe_time_dipole_gravity_null_results.json.
"""
import json
import time
import argparse
import math
import numpy as np

from _route2_obs import (parse_sp3_full, parse_obs_c1, geom_residual,
                         dt_rel_reg, sp3_interp)
from probe_6c_gps_positive import (parse_rinex_nav, estimate_rx_clock,
                                   build_residual, STATION, SP3, NAV, OBS,
                                   F_TRUTH, C_LIGHT)
from kbk_pipeline import compute_operator_matrix


# ----------------------------------------------------------------------------
# Build the two co-sampled channels for a satellite (the S18 GPS pipeline,
# stopping at the two time series instead of running the single-lag regression)
# ----------------------------------------------------------------------------
def build_channels_v2(prn, obs, sats, nav, rxmap, dec):
    """Returns (t, A, B, ecc) sampled on the satellite's epochs:
        A = cleaned pseudorange residual ye(t)            [metres]
        B = c*F*e*sqrt(a)*sin(E(t))  (= predicted dt_rel) [metres]
    Both are time-equivalent clock signals in metres (c*dt). B is computed
    purely from BROADCAST Keplerian elements -> independent of the pseudorange
    that builds A (non-circular)."""
    tt, r = build_residual(prn, obs, sats, dec)
    r = r - np.nanmean(r)
    rxgrid, rxval = rxmap
    rx = np.interp(tt, rxgrid, rxval, left=np.nan, right=np.nan)
    A = r - rx                                  # cleaned residual (metres)
    good = np.isfinite(A)
    te = tt[good]; A = A[good]
    if len(te) < 60:
        return None
    R = dt_rel_reg(nav[prn], te)                # e*sqrt(a)*sin(E)  (sqrt(m))
    B = C_LIGHT * F_TRUTH * R                    # predicted c*dt_rel (metres)
    ecc = float(np.median([e["e"] for e in nav[prn]]))
    return te, A, B, ecc


# ----------------------------------------------------------------------------
# (A) RAW cross-covariance dipole vs lag  ->  does the relativistic coupling
#     sit at lag 0 (in-phase A<->B), beating an off-phase / scramble null?
# ----------------------------------------------------------------------------
def raw_cov_dipole(t, A, B, max_lag_frac=0.25, n_scramble=300, seed=0):
    """Normalized cross-correlation of the two channels over integer-sample lag.
    The relativistic clock signal makes A and B coupled IN PHASE (the S18 fit
    found corr ~ -0.99 at lag 0). The dipole 'solve for time' here = recover that
    the coupling is at lag 0 (no inter-channel delay), and that |cc(0)| beats a
    scramble null. (Unlike LIGO's two detectors there is no light-travel offset
    between A and B -- they are the same epoch grid -- so the time answer is
    'coupling is instantaneous at lag 0', and the test is whether the RAW
    covariance sees the relativistic coupling at all where the entropy dipole is
    blind.)"""
    a = (A - A.mean()) / (A.std() + 1e-30)
    b = (B - B.mean()) / (B.std() + 1e-30)
    n = len(a)
    max_lag = max(1, int(max_lag_frac * n))
    lags = np.arange(-max_lag, max_lag + 1)
    cc = np.empty(len(lags))
    for k, lag in enumerate(lags):
        if lag >= 0:
            x, y = a[lag:], b[:n - lag]
        else:
            x, y = a[:n + lag], b[-lag:]
        cc[k] = np.dot(x, y) / len(x)
    kbest = int(np.argmax(np.abs(cc)))
    best_lag = int(lags[kbest])
    cc0 = float(cc[max_lag])              # lag 0
    # scramble null: permute B, recompute |cc| at lag 0
    rng = np.random.default_rng(seed)
    null = np.empty(n_scramble)
    for i in range(n_scramble):
        bp = rng.permutation(b)
        null[i] = abs(np.dot(a, bp) / n)
    return {
        "peak_abs_cc": float(abs(cc[kbest])),
        "cc_at_lag0": cc0,
        "best_lag_samples": best_lag,
        "coupling_is_at_lag0": bool(best_lag == 0),
        "scramble_null_mean": float(null.mean()),
        "scramble_null_std": float(null.std()),
        "z_over_scramble": float((abs(cc0) - null.mean()) / (null.std() + 1e-12)),
        "p_null_ge": float(np.mean(null >= abs(cc0))),
        "lag_axis_samples": lags[::max(1, len(lags) // 60)].tolist(),
        "cc_values": [round(float(v), 4) for v in cc[::max(1, len(lags) // 60)]],
    }


# ----------------------------------------------------------------------------
# (B) STATIC ALGEBRAIC dipole on windowed entropies of A and B
#     H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2  + circular-shift tautology null
# ----------------------------------------------------------------------------
def r2_fit(Ha2, HaHb):
    X = np.column_stack([np.ones_like(HaHb), HaHb, HaHb * HaHb])
    beta, *_ = np.linalg.lstsq(X, Ha2, rcond=None)
    sst = np.sum((Ha2 - Ha2.mean()) ** 2)
    ssr = np.sum((Ha2 - X @ beta) ** 2)
    return ((1 - ssr / sst) if sst > 0 else float("nan")), beta


def algebraic_dipole_with_null(A, B, win, stride, n_null, seed):
    """Windowed (H_a, H_b) of the two channels, fit the algebraic dipole, and
    apply the DECISIVE circular-shift tautology null (roll H_b, refit). Mirrors
    probe_time_dipole_gravity_null.py analyze()."""
    M = compute_operator_matrix(A, B, win, stride, mi_bins=16)
    if M.shape[0] < 8:
        return None
    Ha, Hb = M[:, 0], M[:, 1]
    Ha2 = Ha * Ha
    HaHb = Ha * Hb
    r2_real, beta = r2_fit(Ha2, HaHb)
    rng = np.random.default_rng(seed)
    npts = len(Hb)
    min_shift = max(3, npts // 10)
    nulls = []
    for _ in range(n_null):
        s = rng.integers(min_shift, npts - min_shift)
        nulls.append(r2_fit(Ha2, Ha * np.roll(Hb, s))[0])
    nulls = np.array(nulls)
    return {
        "n_windows": int(npts),
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


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    n_scram = 100 if args.canary else 300
    n_null = 100 if args.canary else 400
    dec = 4 if args.canary else 2

    nav = parse_rinex_nav(NAV, systems=("G", "E"))
    sats, _ = parse_sp3_full(SP3, systems=("G", "E"))

    circ = ["E02", "E07", "E10", "E11", "E12", "E19", "E24", "E25", "E31", "E33"]
    # eccentric Galileo (signal-rich) + a near-circular GPS CONTROL (tiny signal)
    targets = ["E18", "E14"]
    control = ["G15"]          # e~0.015, ~17 m relativistic signal -> no-coupling control
    want = set(circ + targets + control)
    obs = parse_obs_c1(OBS, want_prns=want)
    circ_have = [p for p in circ if p in obs]
    rxmap = estimate_rx_clock(circ_have, obs, sats, dec)

    # window for the algebraic dipole: ~ a few % of the day, in epoch-samples.
    # 30 s sampling, decimated by `dec`. Use ~80-window with ~50% overlap.
    win = 30 if args.canary else 90
    stride = win // 2

    per_sat = {}
    for p in targets + control:
        if p not in obs or p not in sats or p not in nav:
            continue
        ch = build_channels_v2(p, obs, sats, nav, rxmap, dec)
        if ch is None:
            continue
        te, A, B, ecc = ch
        raw = raw_cov_dipole(te, A, B, n_scramble=n_scram, seed=11)
        alg = algebraic_dipole_with_null(A, B, win, stride, n_null, seed=22)
        per_sat[p] = {
            "prn": p, "ecc": ecc, "n_epochs": int(len(te)),
            "signal_ptp_m": float(np.ptp(B)),
            "is_control": bool(p in control),
            "raw_covariance_dipole": raw,
            "static_algebraic_dipole": alg,
        }

    # headline comparison: eccentric Galileo (signal) vs near-circular control
    ecc_sats = [p for p in per_sat if p in targets]
    ctrl_sats = [p for p in per_sat if p in control]

    def avg(ps, path):
        vals = []
        for p in ps:
            d = per_sat[p]
            for k in path:
                if d is None:
                    break
                d = d[k]
            if d is not None and not (isinstance(d, float) and math.isnan(d)):
                vals.append(d)
        return float(np.mean(vals)) if vals else None

    result = {
        "probe": "the OTHER dipole, to solve for TIME -- carried to the S18 GPS "
                 "time-dilation recovery (+ pulsar decision-gate call)",
        "idea": "flow/entropy dMI/dt dipole was blind; the RAW cross-covariance "
                "and the STATIC ALGEBRAIC H_a^2=a+b*H_aH_b+c*(H_aH_b)^2 dipoles "
                "travelled to LIGO gravity. Do they capture the GR time-coupling "
                "of the S18 GPS recovery? (positive control, NOT a new-law claim)",
        "decision_gate": {
            "GPS": "GENUINE 2-channel object. A=cleaned pseudorange residual "
                   "(measured relativistic clock-residual, metres); B=independent "
                   "broadcast-element driver c*F*e*sqrt(a)*sin(E) (metres). "
                   "Co-sampled, separately-sourced. Dipole operators (cross-cov "
                   "over lag; entropy-algebraic) differ from the S18 single-lag "
                   "regression -> proceeded.",
            "PULSAR": "CATEGORY STRETCH -- NOT fabricated. dP_b/dt and gamma are "
                      "scalar params of a global WLS timing-model fit; no "
                      "co-sampled second physical channel; dP_b/dt is a secular "
                      "cumulative-parabola effect, not an instantaneous channel "
                      "coupling; any 'orbital-phase driver' would be synthesized "
                      "from the same model. Forcing a 2-channel entropy/cov dipole "
                      "onto a 1-parameter secular fit is the chemistry-into-the-4-"
                      "forces mismatch. See pulsar_decision_gate.",
        },
        "config": {"station": "BRUX 2023-001", "decimation": dec,
                   "alg_window_epochs": win, "alg_stride_epochs": stride,
                   "n_scramble": n_scram, "n_shift_null": n_null,
                   "channel_A": "cleaned pseudorange residual ye(t) [m]",
                   "channel_B": "c*F*e*sqrt(a)*sin(E(t)) broadcast driver [m]"},
        "gps_eccentric_galileo_targets": ecc_sats,
        "gps_near_circular_control": ctrl_sats,
        "per_sat": per_sat,
        "headline": {
            "ecc_raw_z_over_scramble_mean": avg(ecc_sats, ["raw_covariance_dipole", "z_over_scramble"]),
            "ctrl_raw_z_over_scramble_mean": avg(ctrl_sats, ["raw_covariance_dipole", "z_over_scramble"]),
            "ecc_alg_R2_mean": avg(ecc_sats, ["static_algebraic_dipole", "R2_real"]),
            "ecc_alg_excess_over_null_mean": avg(ecc_sats, ["static_algebraic_dipole", "excess_over_null"]),
            "ctrl_alg_R2_mean": avg(ctrl_sats, ["static_algebraic_dipole", "R2_real"]),
            "ctrl_alg_excess_over_null_mean": avg(ctrl_sats, ["static_algebraic_dipole", "excess_over_null"]),
        },
        "pulsar_decision_gate": {
            "verdict": "category stretch -- no defensible 2-channel object",
            "reasoning": "PSR B1913+16 recovery (probe_pulsar_time.py) gets "
                         "dP_b/dt and gamma as SCALAR parameters of a global PINT "
                         "WLS fit. There is no second co-sampled PHYSICAL channel: "
                         "the post-reduction residual series is white measurement "
                         "noise + the secular parameter signal, and the only "
                         "'second channel' available (orbital phase) is computed "
                         "from the SAME timing model. dP_b/dt is a cumulative "
                         "PARABOLA (secular), not an instantaneous A<->B coupling, "
                         "so neither the raw cross-covariance-over-lag dipole nor "
                         "the windowed-entropy algebraic dipole has a legitimate "
                         "two-channel input. A constructed channel-B sinusoid "
                         "would make any dipole 'fit' a tautology of the model, "
                         "not a recovery. Clean negative; no pulsar dipole run.",
            "what_would_be_needed": "a genuinely independent co-sampled second "
                         "observable of the same system (e.g. a contemporaneous, "
                         "independently-derived orbital-geometry time series) -- "
                         "not available from the single-pulsar TOA stream.",
        },
        "honesty": "RECOVERY of a KNOWN GR law (-2/c^2), not new physics. A dipole "
                   "hit here is a positive-control that the tool travels. Each "
                   "result is one data point; tautology-null excess reported as-is.",
        "runtime_s": round(time.time() - t0, 1),
    }

    out = ("probe_time_dipole_gpspulsar_canary.json" if args.canary
           else "probe_time_dipole_gpspulsar_results.json")
    with open(out, "w") as f:
        json.dump(result, f, indent=2)

    # ---- console summary ----
    print(f"=== OTHER dipole / solve-for-time: GPS+pulsar "
          f"({'CANARY' if args.canary else 'FULL'}) {result['runtime_s']}s ===")
    print("\nDECISION GATE: GPS = genuine 2-channel (proceeded); "
          "PULSAR = category stretch (written up, not fabricated).")
    for p in targets + control:
        if p not in per_sat:
            continue
        d = per_sat[p]
        raw = d["raw_covariance_dipole"]; alg = d["static_algebraic_dipole"]
        tag = "CONTROL" if d["is_control"] else "ECC-GAL"
        print(f"\n[{p}] {tag}  e={d['ecc']:.4f}  signal={d['signal_ptp_m']:.0f} m  "
              f"n={d['n_epochs']}")
        print(f"   RAW cov dipole : |cc(0)|={abs(raw['cc_at_lag0']):.3f} "
              f"peak|cc|={raw['peak_abs_cc']:.3f} @lag {raw['best_lag_samples']}  "
              f"z_scramble={raw['z_over_scramble']:.0f}  lag0={raw['coupling_is_at_lag0']}")
        if alg:
            print(f"   ALG dipole     : R2={alg['R2_real']:.3f}  "
                  f"shift-null={alg['shift_null_mean']:.3f}  "
                  f"excess={alg['excess_over_null']:+.3f}  z={alg['z_over_null']:.1f}  "
                  f"p={alg['p_null_ge_real']:.3f}")
    h = result["headline"]
    print(f"\nHEADLINE  ecc raw z(scramble) mean = {h['ecc_raw_z_over_scramble_mean']}  "
          f"vs control = {h['ctrl_raw_z_over_scramble_mean']}")
    print(f"          ecc alg R2 mean = {h['ecc_alg_R2_mean']}  "
          f"excess mean = {h['ecc_alg_excess_over_null_mean']}  "
          f"(control R2 {h['ctrl_alg_R2_mean']}, excess {h['ctrl_alg_excess_over_null_mean']})")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
