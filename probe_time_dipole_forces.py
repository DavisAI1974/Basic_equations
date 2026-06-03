"""
DOES THE *OTHER* DIPOLE TRAVEL TO THE GAUGE FORCES?  (gravity-time -> EM / strong / weak)

Context (Greg, S19+): the windowed-entropy "flow dipole" (dMI/d-axis) came back BLIND
on the real forces, but TWO OTHER dipole forms were POSITIVE on gravity AND on GPS
time-dilation, and survived a tautology-killing null:
  (A) RAW cross-covariance dipole over LAG  (for time-series 2-channel data).
  (B) STATIC ALGEBRAIC dipole  H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2  with a
      CIRCULAR-SHIFT tautology-killing null (shift channel B, refit, compare real R2
      to the shift-null; report excess / z / p).
References mirrored EXACTLY: probe_time_dipole_gravity.py + probe_time_dipole_gravity_null.py
(the circular-shift tautology null is mandatory -- it preserves H_b's marginal +
autocorrelation AND the shared H_a factor, killing ONLY the instantaneous H_a<->H_b
pairing, so real >> shift-null => genuine joint structure; real ~= shift-null => the
"constraint" is the shared-H_a tautology).

THE QUESTION: is the OTHER dipole a UNIVERSAL force tool, or gravity/time-specific?
Apply to THREE real systems whose KNOWN 2-channel correlation demonstrably FIRES:
  1. EM    -- HBT photon timetags (split thermal = g2(0)~1.85 correlated; uncorrelated
              = control). TIME SERIES -> BOTH dipole forms apply.
  2. STRONG-- CMS dimuon (100k events: J/psi + Upsilon + Z), channels = leading vs
              subleading muon pT, AXIS = pair invariant mass. EVENT ENSEMBLE binned
              along mass -> ONLY the static algebraic dipole applies (no time lag).
  3. WEAK  -- CMS Z->mumu (10k events), same construction, axis = invariant mass,
              control = same-charge (no resonance).

HONESTY (strict): NO novelty claims. Each system = one data point, no verdicts. The
KEY metric is the tautology-null EXCESS: does the algebraic dipole carry GENUINE
structure beyond the shared-H_a tautology (excess clearly > 0, z > ~3, beating the
system's control), or collapse to tautology (excess ~0) like the gravity NOISE window?
On gravity the genuine excess was the COMMON-SIGNAL coupling at the event (corr Ha,Hb
0.73) -- so we also watch whether corr(Ha,Hb) rises where the known correlation is
strong (at the resonance mass / at the g2 peak).

Run:  python probe_time_dipole_forces.py [--canary]
"""
import json
import time
import argparse
import numpy as np
import pandas as pd

# --- reuse the verbatim EM HBT loader (same as probe_flow_dipole_em_battery.py) ---
from probe_flow_dipole_em_battery import load_tags, bin_counts, RES
# --- reuse the verbatim entropy/MI estimators from the dimuon loader ---
from probe_flow_dipole_weak import entropy_1d, mi_2d


# ============================================================================
# STATIC ALGEBRAIC DIPOLE  +  circular-shift tautology-killing null
# (verbatim logic from probe_time_dipole_gravity_null.py: r2_fit + analyze)
# ============================================================================
def r2_fit(Ha2, HaHb):
    X = np.column_stack([np.ones_like(HaHb), HaHb, HaHb * HaHb])
    beta, *_ = np.linalg.lstsq(X, Ha2, rcond=None)
    sst = np.sum((Ha2 - Ha2.mean()) ** 2)
    ssr = np.sum((Ha2 - X @ beta) ** 2)
    return ((1 - ssr / sst) if sst > 0 else float("nan")), beta


def algebraic_dipole_with_null(Ha, Hb, n_null, seed):
    """Static algebraic dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2, with the
    circular-shift tautology-killing null on H_b along the axis. Returns the same
    fields as probe_time_dipole_gravity_null.analyze."""
    Ha = np.asarray(Ha, dtype=float)
    Hb = np.asarray(Hb, dtype=float)
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
    return {
        "n_points": int(npts),
        "Ha_std": float(np.std(Ha)),
        "Hb_std": float(np.std(Hb)),
        "corr_Ha_Hb": float(np.corrcoef(Ha, Hb)[0, 1]) if npts > 1 else float("nan"),
        "R2_real": float(r2_real),
        "shift_null_mean": float(nulls.mean()),
        "shift_null_p95": float(np.percentile(nulls, 95)),
        "shift_null_std": float(nulls.std()),
        "excess_over_null": float(r2_real - nulls.mean()),
        "z_over_null": float((r2_real - nulls.mean()) / (nulls.std() + 1e-9)),
        "p_null_ge_real": float(np.mean(nulls >= r2_real)),
        "beta": [float(b) for b in beta],
    }


def genuine(blk, excess_thresh=0.10):
    """The verdict criterion from the gravity null script."""
    return bool(blk["p_null_ge_real"] < 0.05 and
                blk["excess_over_null"] > excess_thresh)


# ============================================================================
# (A) EM: windowed entropy series along TIME, + raw cross-cov over LAG
# ============================================================================
def windowed_entropy_series(a, b, win_bins, max_windows, h_bins=12, mi_bins=8):
    """Slide a window of win_bins across the two binned count series; per window
    estimate H_a, H_b, MI. Returns (Ha, Hb, MI) along the TIME axis."""
    n = min(len(a), len(b))
    k = min(max_windows, n // win_bins)
    Ha, Hb, MI = [], [], []
    for i in range(k):
        sa = a[i * win_bins:(i + 1) * win_bins]
        sb = b[i * win_bins:(i + 1) * win_bins]
        Ha.append(entropy_1d(sa, bins=h_bins))
        Hb.append(entropy_1d(sb, bins=h_bins))
        MI.append(mi_2d(sa, sb, bins=mi_bins))
    return np.array(Ha), np.array(Hb), np.array(MI)


def raw_crosscov_lag(a, b, max_lag, n_shift_null, rng):
    """RAW cross-covariance dipole over lag: normalized cross-correlation of the
    two raw binned count series scanned over lag. peak|cc| = the zero-lag coupling
    (HBT same-time bunching). time-shift null = large circular shift destroys the
    same-time alignment while preserving marginals (mirrors the gravity lag_null +
    the EM battery's within-source shift null)."""
    n = min(len(a), len(b))
    a = a[:n].astype(float)
    b = b[:n].astype(float)
    a = (a - a.mean()) / (a.std() + 1e-30)
    b = (b - b.mean()) / (b.std() + 1e-30)
    lags = np.arange(-max_lag, max_lag + 1)
    cc = np.empty(len(lags))
    for k, lag in enumerate(lags):
        if lag >= 0:
            x, y = a[lag:], b[: n - lag]
        else:
            x, y = a[: n + lag], b[-lag:]
        cc[k] = np.dot(x, y) / len(x)
    kbest = int(np.argmax(np.abs(cc)))
    lag_best = int(lags[kbest])
    peak_abscc = float(abs(cc[kbest]))
    # time-shift null: large circular shift of b, recompute zero-lag |cc|
    nsh = max(max_lag * 4, n // 20)
    null = []
    for _ in range(n_shift_null):
        s = int(rng.integers(nsh, n - nsh))
        b_sh = np.roll(b, s)
        null.append(abs(float(np.dot(a, b_sh) / n)))
    null = np.array(null)
    return {
        "best_lag_bins": lag_best,
        "peak_abs_cc": peak_abscc,
        "zero_lag_cc": float(cc[lags == 0][0]),
        "shift_null_mean": float(null.mean()),
        "shift_null_std": float(null.std()),
        "z_over_null": float((peak_abscc - null.mean()) / (null.std() + 1e-12)),
        "cc_curve_lags": lags.tolist()[::max(1, len(lags) // 60)],
        "cc_curve": [round(float(v), 4) for v in cc.tolist()][::max(1, len(lags) // 60)],
    }


def run_em(cfg):
    out = {}
    for src_name, (path, chA, chB) in cfg["em_sources"].items():
        tA, tB, T = load_tags(path, chA, chB)
        dt = cfg["em_dt"]
        a, b = bin_counts(tA, tB, T, dt)
        meta = {"n_tags_A": int(len(tA)), "n_tags_B": int(len(tB)),
                "rate_A_hz": float(len(tA) / T), "rate_B_hz": float(len(tB) / T),
                "T_s": float(T), "bin_dt_s": dt, "n_bins": int(len(a))}
        # (A1) static algebraic dipole on windowed entropies along TIME
        Ha, Hb, MI = windowed_entropy_series(a, b, cfg["em_win_bins"],
                                             cfg["em_max_windows"])
        alg = algebraic_dipole_with_null(Ha, Hb, cfg["n_null"], seed=11)
        alg["mi_mean"] = float(np.mean(MI)) if len(MI) else float("nan")
        # (A2) raw cross-cov over lag (the form the EM battery already found HITS)
        rng = np.random.default_rng(7)
        raw = raw_crosscov_lag(a, b, cfg["em_max_lag_bins"], cfg["n_null"], rng)
        out[src_name] = {"meta": meta, "static_algebraic_dipole": alg,
                         "raw_crosscov_lag_dipole": raw,
                         "algebraic_genuine": genuine(alg)}
    return out


# ============================================================================
# (B) dimuon: bin events along MASS axis -> windowed H_a/H_b/MI -> algebraic dipole
# ============================================================================
def build_mass_axis(M, lead, sub, m_lo, m_hi, ev_per_bin, min_per_bin,
                    h_bins=12, mi_bins=8):
    """Equal-EVENT (quantile) bins along the mass axis; per bin estimate
    H_a=H(pT_lead), H_b=H(pT_sub), MI(lead;sub). Mirrors build_axis_series
    from probe_flow_dipole_weak.py (quantile binning so n_bins >> params)."""
    sel0 = (M >= m_lo) & (M < m_hi)
    Ms, Ls, Ss = M[sel0], lead[sel0], sub[sel0]
    order = np.argsort(Ms)
    Ms, Ls, Ss = Ms[order], Ls[order], Ss[order]
    n = len(Ms)
    k = max(8, n // ev_per_bin)
    edges = np.linspace(0, n, k + 1).astype(int)
    Mc, Ha, Hb, MI = [], [], [], []
    for i in range(k):
        a, b = edges[i], edges[i + 1]
        if b - a < max(40, min_per_bin):
            continue
        Mc.append(float(np.median(Ms[a:b])))
        Ha.append(entropy_1d(Ls[a:b], bins=h_bins))
        Hb.append(entropy_1d(Ss[a:b], bins=h_bins))
        MI.append(mi_2d(Ls[a:b], Ss[a:b], bins=mi_bins))
    return np.array(Mc), np.array(Ha), np.array(Hb), np.array(MI)


def load_dimuon_strong(path="data/strong/dimuon.csv"):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    pt1, pt2 = df.pt1.values, df.pt2.values
    Q1, Q2 = df.Q1.values, df.Q2.values
    M = df.M.values.astype(float)
    lead = np.maximum(pt1, pt2)
    sub = np.minimum(pt1, pt2)
    opp = (Q1 * Q2) < 0
    return dict(M=M, lead=lead, sub=sub, opp=opp)


def load_dimuon_weak(path="data/forces/Zmumu.csv"):
    df = pd.read_csv(path)
    pt1, eta1, phi1 = df.pt1.values, df.eta1.values, df.phi1.values
    pt2, eta2, phi2 = df.pt2.values, df.eta2.values, df.phi2.values
    Q1, Q2 = df.Q1.values, df.Q2.values
    M2 = 2 * pt1 * pt2 * (np.cosh(eta1 - eta2) - np.cos(phi1 - phi2))
    M = np.sqrt(np.clip(M2, 0, None))
    lead = np.maximum(pt1, pt2)
    sub = np.minimum(pt1, pt2)
    opp = (Q1 * Q2) < 0
    return dict(M=M, lead=lead, sub=sub, opp=opp)


def run_dimuon(label, ev, charge, m_lo, m_hi, cfg, resonance=None, seed=11):
    mask = ev["opp"] if charge == "opposite" else ~ev["opp"]
    M, lead, sub = ev["M"][mask], ev["lead"][mask], ev["sub"][mask]
    win = (M >= m_lo) & (M < m_hi)
    if win.sum() < cfg["min_per_bin"] * 8:
        return {"n_events": int(win.sum()), "charge": charge,
                "m_window": [m_lo, m_hi],
                "note": f"too few {charge}-charge events in [{m_lo},{m_hi}] "
                        f"to build a stable mass axis ({int(win.sum())} events).",
                "static_algebraic_dipole": None, "algebraic_genuine": None}
    Mc, Ha, Hb, MI = build_mass_axis(M, lead, sub, m_lo, m_hi,
                                     cfg["ev_per_bin"], cfg["min_per_bin"])
    if len(Mc) < 8:
        return {"n_events": int(win.sum()), "charge": charge,
                "m_window": [m_lo, m_hi],
                "note": f"only {len(Mc)} stable mass bins.",
                "static_algebraic_dipole": None, "algebraic_genuine": None}
    alg = algebraic_dipole_with_null(Ha, Hb, cfg["n_null"], seed=seed)
    # corr(Ha,Hb) at the resonance bins vs off-resonance bins (the "common-signal
    # coupling rises where the known correlation is strong" analog of gravity's
    # corr 0.73 at the event)
    corr_res = None
    if resonance is not None and len(Mc) >= 8:
        near = np.abs(Mc - resonance) <= cfg.get("res_halfwidth", 6.0)
        far = ~near
        if near.sum() >= 3 and far.sum() >= 3:
            cr = (float(np.corrcoef(Ha[near], Hb[near])[0, 1])
                  if near.sum() > 1 else float("nan"))
            cf = (float(np.corrcoef(Ha[far], Hb[far])[0, 1])
                  if far.sum() > 1 else float("nan"))
            corr_res = {"resonance_GeV": resonance, "n_near": int(near.sum()),
                        "n_far": int(far.sum()), "corr_HaHb_near_resonance": cr,
                        "corr_HaHb_off_resonance": cf}
    return {"n_events": int(win.sum()), "charge": charge,
            "n_axis_points": int(len(Mc)), "m_window": [m_lo, m_hi],
            "mass_axis_GeV": [round(float(x), 2) for x in Mc.tolist()],
            "MI_along_axis": [round(float(x), 4) for x in MI.tolist()],
            "corr_HaHb_resonance_vs_offres": corr_res,
            "static_algebraic_dipole": alg,
            "algebraic_genuine": genuine(alg)}


# ============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canary", action="store_true")
    args = ap.parse_args()
    t0 = time.time()

    if args.canary:
        cfg = {
            "em_sources": {  # canary: signal only
                "split_correlated": ("data/em/splitThermal_rawtags.txt", 11, 15),
            },
            "em_dt": 2e-6, "em_win_bins": 4000, "em_max_windows": 80,
            "em_max_lag_bins": 30,
            "ev_per_bin": 200, "min_per_bin": 60, "res_halfwidth": 6.0,
            "n_null": 100,
            "strong_path": "data/strong/dimuon.csv",
            "weak_path": "data/forces/Zmumu.csv",
            "strong_max_events": 30000,
        }
        out_path = "probe_time_dipole_forces_canary.json"
    else:
        cfg = {
            "em_sources": {
                "split_correlated": ("data/em/splitThermal_rawtags.txt", 11, 15),
                "uncorrelated_control": ("data/em/uncorrelatedThermal_rawtags.txt", 15, 16),
            },
            "em_dt": 2e-6, "em_win_bins": 4000, "em_max_windows": 300,
            "em_max_lag_bins": 60,
            "ev_per_bin": 150, "min_per_bin": 60, "res_halfwidth": 6.0,
            "n_null": 400,
            "strong_path": "data/strong/dimuon.csv",
            "weak_path": "data/forces/Zmumu.csv",
            "strong_max_events": None,
        }
        out_path = "probe_time_dipole_forces_results.json"

    results = {}

    # ----------------------------- EM (time series: both dipole forms) ----
    print("[EM] HBT photon timetags ...", flush=True)
    results["EM_HBT"] = run_em(cfg)

    # ----------------------------- STRONG (dimuon event ensemble) ---------
    print("[STRONG] CMS dimuon (J/psi+Upsilon+Z) ...", flush=True)
    ev_s = load_dimuon_strong(cfg["strong_path"])
    if cfg["strong_max_events"] and len(ev_s["M"]) > cfg["strong_max_events"]:
        n = cfg["strong_max_events"]
        for kk in ("M", "lead", "sub", "opp"):
            ev_s[kk] = ev_s[kk][:n]
    # full spectrum (J/psi 3.1, Upsilon 9.5, Z 91), opposite-charge = real dimuons.
    # control = same-charge (combinatorial, no resonance). watch corr(Ha,Hb) near
    # the J/psi resonance (the dominant peak in this 2-110 GeV sample).
    strong = {
        "opposite_full_spectrum": run_dimuon(
            "strong", ev_s, "opposite", 2.0, 110.0, cfg,
            resonance=3.1, seed=11),
        "same_charge_control": run_dimuon(
            "strong", ev_s, "same", 2.0, 110.0, cfg, resonance=None, seed=12),
    }
    results["STRONG_dimuon"] = strong

    # ----------------------------- WEAK (Z->mumu event ensemble) ----------
    print("[WEAK] CMS Z->mumu ...", flush=True)
    ev_w = load_dimuon_weak(cfg["weak_path"])
    weak = {
        "opposite_Z_window": run_dimuon(
            "weak", ev_w, "opposite", 60.0, 120.0, cfg,
            resonance=91.1876, seed=11),
        "same_charge_control": run_dimuon(
            "weak", ev_w, "same", 60.0, 120.0, cfg, resonance=None, seed=12),
    }
    results["WEAK_Zmumu"] = weak

    # ----------------------------- summary table --------------------------
    def row(name, alg, ctl_alg, raw=None):
        if alg is None:
            return {"system": name, "note": "no stable axis / too few events"}
        r = {
            "system": name,
            "R2_real": round(alg["R2_real"], 3),
            "tautology_null_excess": round(alg["excess_over_null"], 3),
            "z": round(alg["z_over_null"], 1),
            "p": round(alg["p_null_ge_real"], 4),
            "corr_Ha_Hb": round(alg["corr_Ha_Hb"], 3),
            "control_excess": (round(ctl_alg["excess_over_null"], 3)
                               if ctl_alg else None),
            "algebraic_genuine": genuine(alg),
            "travels": bool(genuine(alg) and
                            (ctl_alg is None or
                             alg["excess_over_null"] >
                             ctl_alg["excess_over_null"] + 0.05)),
        }
        if raw is not None:
            r["raw_crosscov_peak_abscc"] = round(raw["peak_abs_cc"], 3)
            r["raw_crosscov_z"] = round(raw["z_over_null"], 1)
            r["raw_crosscov_best_lag_bins"] = raw["best_lag_bins"]
        return r

    em = results["EM_HBT"]
    em_sig = em.get("split_correlated", {})
    em_ctl = em.get("uncorrelated_control", {})
    summary = []
    summary.append(row(
        "EM_HBT_split(algebraic)",
        em_sig.get("static_algebraic_dipole"),
        em_ctl.get("static_algebraic_dipole") if em_ctl else None,
        raw=em_sig.get("raw_crosscov_lag_dipole")))
    summary.append(row(
        "STRONG_dimuon_oppcharge",
        strong["opposite_full_spectrum"].get("static_algebraic_dipole"),
        strong["same_charge_control"].get("static_algebraic_dipole")))
    summary.append(row(
        "WEAK_Zmumu_oppcharge",
        weak["opposite_Z_window"].get("static_algebraic_dipole"),
        weak["same_charge_control"].get("static_algebraic_dipole")))

    out = {
        "probe": "does the OTHER dipole (raw-cross-cov + static-algebraic, "
                 "tautology-null-tested) TRAVEL from gravity-time to the GAUGE "
                 "FORCES (EM HBT / strong+weak dimuon)?",
        "honesty": "NO novelty claims. each system = one data point, no verdicts. "
                   "the KEY metric is the tautology-null EXCESS: genuine joint "
                   "structure (excess>0, z>~3, beats control) vs collapse-to-"
                   "tautology (excess~0, like the gravity NOISE window). report "
                   "excess honestly even if ~0.",
        "method": "static algebraic dipole H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2 "
                  "with CIRCULAR-SHIFT tautology null on H_b (verbatim from "
                  "probe_time_dipole_gravity_null.py). EM is a TIME SERIES so the "
                  "RAW cross-cov-over-lag dipole + time-shift null also applies; "
                  "dimuon are EVENT ENSEMBLES binned along the mass axis so only "
                  "the static algebraic form applies (no time lag).",
        "gravity_reference": {"event_genuine": True,
                              "event_corr_Ha_Hb": 0.73,
                              "noise_collapsed_to_tautology": True},
        "config": {k: v for k, v in cfg.items() if k != "em_sources"},
        "results": results,
        "summary_table": summary,
        "runtime_s": round(time.time() - t0, 1),
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    # console
    print(f"\n=== OTHER dipole -> gauge forces ({'CANARY' if args.canary else 'FULL'}) "
          f"{out['runtime_s']}s ===")
    hdr = (f"{'system':30s} {'R2':>6s} {'excess':>8s} {'z':>7s} {'p':>7s} "
           f"{'corr':>6s} {'ctl_exc':>8s} {'travels':>8s}")
    print(hdr)
    for r in summary:
        if "note" in r:
            print(f"{r['system']:30s}  {r['note']}")
            continue
        ce = r["control_excess"]
        ce = f"{ce:+.3f}" if ce is not None else "  n/a "
        print(f"{r['system']:30s} {r['R2_real']:6.3f} "
              f"{r['tautology_null_excess']:+8.3f} {r['z']:7.1f} {r['p']:7.4f} "
              f"{r['corr_Ha_Hb']:+6.3f} {ce:>8s} {str(r['travels']):>8s}")
    # EM raw-cov line
    rsig = em_sig.get("raw_crosscov_lag_dipole")
    if rsig:
        print(f"\n[EM raw cross-cov over lag] peak|cc|={rsig['peak_abs_cc']:.3f} "
              f"z={rsig['z_over_null']:.1f} best_lag={rsig['best_lag_bins']} bins "
              f"(KNOWN HBT same-time bunching; the time-series form)")
    print(f"-> {out_path}")


if __name__ == "__main__":
    main()
