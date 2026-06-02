"""
PROBE (backlog O3, S18): NAIL the GW170817 binary-neutron-star chirp mass M_c
from RAW LIGO strain by matched filtering against an inspiral template bank.

WHY: S17's CWT-ridge recovered the inspiral LAW FORM (u = f^(-8/3) linear in t,
R^2 0.88) for GW170817 but the ABSOLUTE M_c was biased (15.7 vs catalog 1.20)
because the visible H1 arc was a narrow low-frequency band (51-76 Hz, short
lever arm). O3 = get the BNS chirp mass right. A BNS inspiral is in-band for
~60-100 s sweeping 30 -> ~2000 Hz; the repo's 32 s file is too short/narrow. So
we (a) FETCH the GWOSC 4096 s / 4096 Hz H1+L1 files (merger GPS 1187008882.43),
(b) matched-filter a TaylorF2 inspiral bank spanning chirp mass over PSD-whitened
data, glitch-GATING the famous L1 transient ~1.1 s before merger, and (c) report
the M_c that maximizes network SNR.

CATALOG (conjecture-to-check, NOT cited as support): detector-frame M_c ~ 1.1977
Msun; source-frame ~ 1.186-1.188 Msun; network SNR ~ 32.4 (H1 ~18.8, L1 ~25.9 in
the discovery paper, V1 ~2).

METHOD (matched filtering, the standard CBC detection statistic):
  - Load 4096 s strain; crop a SEG_DUR-second analysis segment ending ~2 s after
    the merger so the full low-freq inspiral is present.
  - Estimate PSD by Welch (pycbc) on the segment; truncate inverse spectrum.
  - L1: GATE out the loud glitch window before whitening (taper to zero).
  - Build a TaylorF2 template bank: scan chirp mass M_c on a grid with the
    component masses set equal-ish (q=1 for BNS; M_c dominates inspiral phasing).
    For each template: matched_filter -> complex SNR time series; record peak |SNR|
    within a tight window around the known merger time (and report the unconstrained
    peak too as a cross-check).
  - The M_c maximizing per-detector SNR is the single-detector recovery; combine
    H1+L1 in quadrature at matched peak time (within light-travel + tolerance) for
    the network statistic. Refine around the coarse-bank peak.
  - Background/null: distribution of off-source peak |SNR| (windows far from the
    event) gives the noise floor the on-source peak must stand above.

Speaking posture (Rule C): I expect the SNR-vs-M_c curve to peak near ~1.20 Msun
with high network SNR if the method is sound; but real BNS data is long and the L1
glitch + PSD choices matter -- wait on the data, no verdict on first look, keep odd
outputs and diagnose. Catalog values are targets to compare against, not support.

Run:  python probe_o3_gw170817_chirpmass.py --canary   (coarse bank, short seg)
      python probe_o3_gw170817_chirpmass.py             (full bank + both dets + null)
"""
import sys
import json
import time
import numpy as np

import h5py
from pycbc.types import TimeSeries
from pycbc.waveform import get_fd_waveform
from pycbc.filter import matched_filter, sigma
from pycbc.psd import welch, interpolate, inverse_spectrum_truncation

CANARY = "--canary" in sys.argv

# ---- constants / catalog -------------------------------------------------
MERGER_GPS = 1187008882.43
CAT_MC_DET = 1.1977     # detector-frame chirp mass, Msun (target, conjecture)
CAT_MC_SRC = 1.187      # source-frame chirp mass, Msun (z~0.0099)
LIGHT_TRAVEL_HL = 0.010 # s, max Hanford-Livingston light travel time

BULK = "/home/user/Basic_equations/data/ligo_bulk"
FILES = {
    "H1": f"{BULK}/H-H1_GWOSC_4KHZ_R1-1187006835-4096.hdf5",
    "L1": f"{BULK}/L-L1_GWOSC_4KHZ_R1-1187006835-4096.hdf5",
}
# L1 loud glitch: ~1.1 s before merger (LVC GW170817 discovery paper)
L1_GLITCH_CENTER = MERGER_GPS - 1.05  # glitch peak measured at t-merger = -1.051 s
L1_GLITCH_HALFWIDTH = 0.50   # s, gate +/- this around the glitch (cosine-tapered)
L1_GLITCH_TAPER = 0.25       # s, taper width for pycbc .gate()


def load_strain_ts(path):
    """Load GWOSC hdf5 strain into a pycbc TimeSeries with correct epoch."""
    with h5py.File(path, "r") as f:
        d = f["strain/Strain"]
        x = d[:].astype(np.float64)
        dt = float(d.attrs["Xspacing"])
        t0 = float(d.attrs["Xstart"])
    return TimeSeries(x, delta_t=dt, epoch=t0)


def chirp_to_components(mc, q=1.0):
    """Chirp mass + mass ratio -> (m1, m2). q = m2/m1 <= 1."""
    eta = q / (1.0 + q) ** 2
    mtot = mc / eta ** 0.6
    m1 = mtot / (1.0 + q)
    m2 = mtot - m1
    return float(m1), float(m2)


def gate_glitch(ts, center, halfwidth, taper):
    """Gate out [center-halfwidth, center+halfwidth] using pycbc's tapered
    gate (cosine taper of width `taper` s on each side). A hand-rolled
    zero-with-narrow-taper rings the matched filter; pycbc's gate is clean."""
    return ts.gate(center, window=halfwidth, taper_width=taper,
                   method="taper", copy=True)


def prep_segment(det, seg_dur, f_low, do_gate):
    """Load, crop a seg_dur-s segment ending ~2 s after merger, estimate PSD."""
    full = load_strain_ts(FILES[det])
    fs = int(round(1.0 / full.delta_t))
    seg_end = MERGER_GPS + 2.0
    seg_start = seg_end - seg_dur
    i0 = int(round((seg_start - float(full.start_time)) * fs))
    i1 = i0 + seg_dur * fs
    seg = full.time_slice(seg_start, seg_end)
    # remove non-finite just in case
    seg.data[~np.isfinite(seg.data)] = 0.0
    if det == "L1" and do_gate:
        seg = gate_glitch(seg, L1_GLITCH_CENTER, L1_GLITCH_HALFWIDTH, L1_GLITCH_TAPER)
    # PSD via Welch (4 s segments), interpolate to data df, inverse-truncate
    seg_psd_len = 4 * fs
    psd = seg.psd(4)
    psd = interpolate(psd, seg.delta_f)
    psd = inverse_spectrum_truncation(psd, int(4 * fs), low_frequency_cutoff=f_low)
    return seg, psd, fs


def snr_for_mc(seg, psd, mc, fs, f_low, f_final, q=1.0):
    """Generate TaylorF2 template at chirp mass mc, matched-filter, return
    the complex SNR TimeSeries (in detector epoch)."""
    m1, m2 = chirp_to_components(mc, q)
    delta_f = seg.delta_f
    flen = len(seg) // 2 + 1
    hp, _ = get_fd_waveform(
        approximant="TaylorF2",
        mass1=m1, mass2=m2,
        delta_f=delta_f, f_lower=f_low, f_final=f_final,
    )
    hp.resize(flen)
    snr = matched_filter(hp, seg, psd=psd, low_frequency_cutoff=f_low)
    return snr, (m1, m2)


def peak_in_window(snr, fs, center, halfwidth):
    """Return (peak_abs_snr, peak_gps) of |snr| within [center-hw, center+hw]."""
    t0 = float(snr.start_time)
    idx = np.arange(len(snr)) / fs + t0
    m = (idx >= center - halfwidth) & (idx <= center + halfwidth)
    a = np.abs(snr.numpy())
    a_in = np.where(m, a, -1.0)
    ip = int(np.argmax(a_in))
    return float(a[ip]), float(idx[ip])


def global_peak(snr, fs, crop_edge=4.0):
    """Unconstrained peak away from segment edges."""
    a = np.abs(snr.numpy())
    n = len(a)
    e = int(crop_edge * fs)
    seg = a[e:n - e]
    ip = e + int(np.argmax(seg))
    t0 = float(snr.start_time)
    return float(a[ip]), float(ip / fs + t0)


def run():
    t_start = time.time()
    if CANARY:
        seg_dur = 64
        f_low = 30.0
        f_final = 512.0
        mc_grid = np.round(np.arange(1.00, 2.01, 0.10), 4)
        dets = ["H1", "L1"]
        refine = False
    else:
        seg_dur = 256
        f_low = 25.0
        f_final = 1024.0
        mc_grid = np.round(np.arange(1.05, 1.55, 0.02), 4)
        dets = ["H1", "L1"]
        refine = True

    print(f"[{'CANARY' if CANARY else 'FULL'}] seg_dur={seg_dur}s f_low={f_low} "
          f"f_final={f_final} grid={mc_grid[0]}..{mc_grid[-1]} ({len(mc_grid)} pts) "
          f"dets={dets}", flush=True)

    out = {
        "probe": "O3_gw170817_chirpmass",
        "mode": "canary" if CANARY else "full",
        "method": "matched filtering, TaylorF2 inspiral chirp-mass bank, "
                  "PSD-whitened, L1 glitch gated, peak SNR within +/-0.1s of "
                  "merger GPS",
        "data": {
            "source": "GWOSC GWTC-1-confident GW170817 v3",
            "files": {d: FILES[d].split('/')[-1] for d in dets},
            "merger_gps": MERGER_GPS,
            "seg_dur_s": seg_dur, "f_low": f_low, "f_final": f_final,
        },
        "catalog": {"Mc_detframe_Msun": CAT_MC_DET, "Mc_srcframe_Msun": CAT_MC_SRC,
                    "note": "catalog treated as conjecture-to-check"},
        "per_detector": {},
        "mc_grid": mc_grid.tolist(),
    }

    prepped = {}
    for det in dets:
        print(f"  prep {det} ...", flush=True)
        seg, psd, fs = prep_segment(det, seg_dur, f_low, do_gate=True)
        prepped[det] = (seg, psd, fs)

    # ---- coarse bank per detector ----
    per_det_curves = {}
    for det in dets:
        seg, psd, fs = prepped[det]
        snrs = []
        for mc in mc_grid:
            snr, (m1, m2) = snr_for_mc(seg, psd, mc, fs, f_low, f_final)
            pk, pkt = peak_in_window(snr, fs, MERGER_GPS, 0.10)
            gpk, gpkt = global_peak(snr, fs)
            snrs.append((float(mc), pk, pkt, gpk, gpkt, m1, m2))
            print(f"    {det} Mc={mc:.3f} peak|SNR|@merger={pk:.2f} "
                  f"(t={pkt-MERGER_GPS:+.4f}s)  global={gpk:.2f}", flush=True)
        per_det_curves[det] = snrs
        best = max(snrs, key=lambda r: r[1])
        out["per_detector"][det] = {
            "coarse_best_Mc": best[0],
            "coarse_best_snr": best[1],
            "coarse_best_peak_dt": best[2] - MERGER_GPS,
            "curve": [{"Mc": r[0], "snr_at_merger": r[1],
                       "peak_dt": r[2] - MERGER_GPS,
                       "global_snr": r[3], "global_dt": r[4] - MERGER_GPS}
                      for r in snrs],
        }

    # ---- refine around the combined coarse peak (full mode) ----
    if refine:
        # network coarse-best: sum of H1,L1 SNR^2 at each Mc (peaks aligned by window)
        net = []
        for k, mc in enumerate(mc_grid):
            ss = 0.0
            for det in dets:
                ss += per_det_curves[det][k][1] ** 2
            net.append((float(mc), float(np.sqrt(ss))))
        net_best_mc = max(net, key=lambda r: r[1])[0]
        out["network_coarse"] = {"curve": [{"Mc": m, "net_snr": s} for m, s in net],
                                 "best_Mc": net_best_mc}
        lo = max(mc_grid[0], net_best_mc - 0.05)
        hi = min(mc_grid[-1], net_best_mc + 0.05)
        fine_grid = np.round(np.arange(lo, hi + 1e-9, 0.005), 4)
        print(f"  REFINE around Mc={net_best_mc}: {fine_grid[0]}..{fine_grid[-1]} "
              f"({len(fine_grid)} pts)", flush=True)
        fine = {det: [] for det in dets}
        for det in dets:
            seg, psd, fs = prepped[det]
            for mc in fine_grid:
                snr, _ = snr_for_mc(seg, psd, mc, fs, f_low, f_final)
                pk, pkt = peak_in_window(snr, fs, MERGER_GPS, 0.10)
                fine[det].append((float(mc), pk, pkt))
                print(f"    fine {det} Mc={mc:.3f} |SNR|={pk:.2f} "
                      f"(t={pkt-MERGER_GPS:+.4f}s)", flush=True)
        net_fine = []
        for k, mc in enumerate(fine_grid):
            ss = sum(fine[det][k][1] ** 2 for det in dets)
            net_fine.append((float(mc), float(np.sqrt(ss))))
        best_net = max(net_fine, key=lambda r: r[1])
        out["network_fine"] = {
            "grid": fine_grid.tolist(),
            "curve": [{"Mc": m, "net_snr": s} for m, s in net_fine],
            "best_Mc": best_net[0], "best_net_snr": best_net[1],
        }
        for det in dets:
            b = max(fine[det], key=lambda r: r[1])
            out["per_detector"][det]["fine_best_Mc"] = b[0]
            out["per_detector"][det]["fine_best_snr"] = b[1]
            out["per_detector"][det]["fine_best_peak_dt"] = b[2] - MERGER_GPS
            out["per_detector"][det]["fine_curve"] = [
                {"Mc": r[0], "snr": r[1], "peak_dt": r[2] - MERGER_GPS}
                for r in fine[det]]

        rec_mc = best_net[0]
        out["recovery"] = {
            "Mc_recovered_detframe": rec_mc,
            "network_snr": best_net[1],
            "vs_catalog_detframe_ratio": rec_mc / CAT_MC_DET,
            "vs_catalog_detframe_pct_err": 100 * (rec_mc - CAT_MC_DET) / CAT_MC_DET,
            "per_det_best_Mc": {d: out["per_detector"][d]["fine_best_Mc"] for d in dets},
            "per_det_best_snr": {d: out["per_detector"][d]["fine_best_snr"] for d in dets},
        }

        # ---- null / background: off-source peak SNR at the recovered Mc ----
        # slide the matched-filter peak search to windows far from the event,
        # to characterize the noise floor the on-source peak stands above.
        null_peaks = {}
        for det in dets:
            seg, psd, fs = prepped[det]
            snr, _ = snr_for_mc(seg, psd, rec_mc, fs, f_low, f_final)
            a = np.abs(snr.numpy())
            t0 = float(snr.start_time)
            idx = np.arange(len(snr)) / fs + t0
            # off-source: > 0.5 s from merger, away from the inverse-spectrum-
            # truncation-corrupted edges (crop EDGE_CROP s each end), and for L1
            # also exclude the gated-glitch region (its residual is not clean).
            EDGE_CROP = 16.0  # PSD inverse-truncation corruption spans ~filter length
            offmask = (np.abs(idx - MERGER_GPS) > 0.5) & \
                      (idx > t0 + EDGE_CROP) & (idx < t0 + (seg_dur - EDGE_CROP))
            if det == "L1":
                offmask &= (np.abs(idx - L1_GLITCH_CENTER) >
                            (L1_GLITCH_HALFWIDTH + L1_GLITCH_TAPER + 0.5))
            off = a[offmask]
            on_pk, on_t = peak_in_window(snr, fs, MERGER_GPS, 0.10)
            null_peaks[det] = {
                "on_source_peak_snr": on_pk,
                "edge_crop_s": EDGE_CROP,
                "off_source_mean": float(np.mean(off)),
                "off_source_std": float(np.std(off)),
                "off_source_p99": float(np.percentile(off, 99)),
                "off_source_max": float(np.max(off)),
                "on_over_off_p99": float(on_pk / np.percentile(off, 99)),
                "on_minus_offmean_in_sigma": float((on_pk - np.mean(off)) / np.std(off)),
            }
            print(f"  NULL {det}: on={on_pk:.2f} off_p99={null_peaks[det]['off_source_p99']:.2f} "
                  f"off_max={null_peaks[det]['off_source_max']:.2f} "
                  f"(on is {null_peaks[det]['on_minus_offmean_in_sigma']:.1f} sigma over off-mean)",
                  flush=True)
        out["null"] = null_peaks
    else:
        # canary: just report per-det coarse bests + a crude network
        net = []
        for k, mc in enumerate(mc_grid):
            ss = sum(per_det_curves[det][k][1] ** 2 for det in dets)
            net.append((float(mc), float(np.sqrt(ss))))
        best_net = max(net, key=lambda r: r[1])
        out["network_coarse"] = {"curve": [{"Mc": m, "net_snr": s} for m, s in net],
                                 "best_Mc": best_net[0], "best_net_snr": best_net[1]}
        out["recovery"] = {
            "Mc_recovered_detframe_COARSE": best_net[0],
            "network_snr_coarse": best_net[1],
            "vs_catalog_detframe_ratio": best_net[0] / CAT_MC_DET,
        }

    out["elapsed_s"] = round(time.time() - t_start, 1)
    suffix = "_canary" if CANARY else "_results"
    path = f"/home/user/Basic_equations/probe_o3_gw170817_chirpmass{suffix}.json"
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWROTE {path}  (elapsed {out['elapsed_s']}s)", flush=True)
    if "recovery" in out:
        print("RECOVERY:", json.dumps(out["recovery"], indent=2), flush=True)


if __name__ == "__main__":
    run()
