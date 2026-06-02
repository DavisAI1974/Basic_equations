"""
PROBE 1 decision gate (Result Discipline, mapped before building):
Is the event-window sampling enough to fit a rate-DYNAMICS model for gravity
(real LIGO), or is finer windowing / a different phase / louder events needed?

A rate-DYNAMICS fit (unlike P3's rate-AXIS regression) needs LAGGED structure:
dMI/dt(t) ~ f(operators(t-1..t-p)). That needs enough *independent* samples
relative to the number of lag parameters. This script measures, on the one
cleanly-aligned event (GW150914), exactly how many windows the event window
gives, the overlap-driven effective-independent count, and the dMI/dt
autocorrelation noise floor -- the numbers that decide ordering.

NO interpretation of gravity; this is a feasibility/sampling measurement only.
"""
import numpy as np

from kbk_pipeline import compute_operator_matrix
from grav_time_retaining import (load_strain, preprocess, FS, EDGE_S,
                                  WIN_S, STRIDE_S, EVENTS)

WIN = int(WIN_S * FS)        # 512 samples
STR = int(STRIDE_S * FS)     # 128 samples
MERGER_T = 16.4              # GW150914 into the 32 s segment
HALF_EVENT = 0.5             # +/- 0.5 s event window (S15 convention)


def autocorr(x, maxlag):
    x = x - x.mean()
    v = np.dot(x, x)
    return [1.0] + [float(np.dot(x[:-k], x[k:]) / v) for k in range(1, maxlag + 1)]


def main():
    print("[probe1_gate] sampling feasibility for a rate-DYNAMICS fit\n", flush=True)
    print(f"windowing: WIN={WIN} samp ({WIN_S*1000:.0f} ms), STRIDE={STR} samp "
          f"({STRIDE_S*1000:.2f} ms), overlap={100*(1-STR/WIN):.0f}%", flush=True)
    print(f"independent-window stride factor = WIN/STRIDE = {WIN/STR:.0f} "
          f"(only ~1 of every {WIN//STR} windows is non-overlapping)\n", flush=True)

    # GW150914 (cleanly aligned per S15)
    h1f, l1f, label, mt, shift, inv = EVENTS[0]
    h1 = preprocess(load_strain(h1f))
    l1 = preprocess(load_strain(l1f))
    l1 = np.roll(l1, int(shift * FS)) * (-1.0 if inv else 1.0)
    crop = int(EDGE_S * FS)
    M = compute_operator_matrix(h1[crop:-crop], l1[crop:-crop], WIN, STR, mi_bins=16)
    n = M.shape[0]
    centers = (np.arange(n) * STR + WIN / 2) / FS + EDGE_S
    ev = np.abs(centers - mt) <= HALF_EVENT

    n_event = int(ev.sum())
    n_indep_event = n_event * STR / WIN  # heavy-overlap -> effective independent
    print(f"FULL segment : {n} windows over {28.0:.0f}s usable", flush=True)
    print(f"EVENT (+/-{HALF_EVENT}s): {n_event} windows  ~= "
          f"{n_indep_event:.0f} INDEPENDENT windows (after 75% overlap)\n", flush=True)

    # dMI/dt in the event window + its autocorrelation noise floor
    MI = M[:, 5]
    dMI = np.gradient(MI) / STRIDE_S
    dMI_ev = dMI[ev]
    maxlag = min(8, n_event - 2)
    ac = autocorr(dMI_ev, maxlag)
    se = 1.0 / np.sqrt(n_indep_event)  # autocorr SE on INDEPENDENT count
    print(f"dMI/dt event-window autocorrelation (lag: value):", flush=True)
    for k, v in enumerate(ac):
        flag = "  <- > 2*SE" if k > 0 and abs(v) > 2 * se else ""
        print(f"   lag {k}: {v:+.3f}{flag}", flush=True)
    print(f"\n2*SE detectability floor on ~{n_indep_event:.0f} indep samples = "
          f"+/-{2*se:.3f}", flush=True)
    print("(any |autocorr| below this floor is indistinguishable from zero)\n", flush=True)

    # parameter budget for a lagged model
    print("rate-DYNAMICS parameter budget (6 operators x p lags):", flush=True)
    for p in (1, 2, 3):
        params = 6 * p
        print(f"   p={p} lags -> {params} params; "
              f"samples/param on event-indep = {n_indep_event/params:.1f}", flush=True)
    print("\n[gate verdict printed in chat]", flush=True)


if __name__ == "__main__":
    main()
