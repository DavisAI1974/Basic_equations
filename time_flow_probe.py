"""
Quick probe of Greg's S15-eve hunch: "the flow dipole expresses itself AS time."

Two sharp, cheap questions, run on the simulated per-domain harness
(no external data). Tests the hunch in two distinct senses:

  Q1 (time-ORDER): is the per-domain null direction sensitive to the
      time-ORDERING of the operator rows?  We permute the rows of the
      operator matrix M and re-extract.  The null comes from the row
      covariance M^T M, which is permutation-invariant over rows, so
      this is expected to be an EXACT no-op -- demonstrated numerically.
      If true: the flow/coupling, as we extract it, does NOT encode
      time-as-sequence.  "flow = time-order" is falsified by construction.

  Q2 (time-DERIVATIVE / rate): does the flow/coupling prefer to live in
      the RATE operators [dH_a/dt, dH_b/dt, dMI/dt] rather than the LEVEL
      operators?  We build an extended, per-column-standardized basis
      [levels(6) | rates(3)] and extract the smallest null, then report
      what FRACTION of the null's squared loading sits on the rate block.
      If "flow = time(-derivative)", biology's coupling constraint should
      load heavily on the rate block.  Foil: physics (self-pole) should
      stay on the level block.

Deflationary readings kept present:
  - Q1 being an exact no-op is partly DEFINITIONAL (extraction is a row
    covariance). It still answers the as-stated hunch: the null is not
    the time-axis.
  - Per-column standardization in Q2 is a procedure choice (INFO-024:
    null DIRECTION is procedure-dependent when column variances differ).
    We therefore report the rate-fraction under standardization AND the
    raw operator means/stds so the scaling is auditable.

Speaking posture: no verdict in advance, no verdict on first look.
This is a canary to sharpen the hunch, not to settle it.

Usage:
  python time_flow_probe.py            # full-ish (N_ens=400, T=20)
  python time_flow_probe.py --canary   # fast (N_ens=200, T=10)
"""
import argparse
import json
import time

import numpy as np

from per_domain_kbk import (
    simulate_physics,
    simulate_biology,
    build_ensemble_operator_matrix,
    entropy_hist_1d,
)
from kbk_pipeline import extract_v1, mi_hist_2d, OP_NAMES

RATE_NAMES = ["dH_a", "dH_b", "dMI"]

DOMAINS = {
    "physics_duffing": simulate_physics,   # self-pole foil
    "biology_lotka": simulate_biology,     # the real coupler
}


def smallest_null(M):
    """Smallest-singular-vector null direction of (centered) M."""
    out = extract_v1(M)
    v = out[0] if isinstance(out, (tuple, list)) else out
    return np.asarray(v, dtype=float)


def cos(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na <= 0 or nb <= 0:
        return float("nan")
    return float(np.dot(a, b) / (na * nb))


def build_rates(M, dt):
    """Finite-difference time derivatives of the 3 marginal/MI levels
    (columns H_a=0, H_b=1, MI=5) along the operator trajectory.

    Returns rates aligned to M[:-1] (forward difference)."""
    dHa = np.diff(M[:, 0]) / dt
    dHb = np.diff(M[:, 1]) / dt
    dMI = np.diff(M[:, 5]) / dt
    return np.column_stack([dHa, dHb, dMI])


def standardize(A):
    mu = A.mean(axis=0)
    sd = A.std(axis=0) + 1e-12
    return (A - mu) / sd, mu, sd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--canary", action="store_true")
    p.add_argument("--seeds", type=int, nargs="+", default=[11, 22, 33])
    p.add_argument("--out", type=str, default="time_flow_probe_results.json")
    args = p.parse_args()

    if args.canary:
        cfg = dict(N_ens=200, T=10.0, dt=0.02)
        args.out = "time_flow_probe_canary.json"
    else:
        cfg = dict(N_ens=400, T=20.0, dt=0.02)
    dt = cfg["dt"]

    print(f"[time_flow_probe] cfg={cfg} seeds={args.seeds}", flush=True)
    t0 = time.time()
    results = {}

    for name, simfn in DOMAINS.items():
        per_seed = []
        for seed in args.seeds:
            X1, X2 = simfn(seed=seed, **cfg)
            M, _ = build_ensemble_operator_matrix(X1, X2)

            # ---- Q1: time-ORDER sensitivity ----------------------------
            v_lvl = smallest_null(M)
            rng = np.random.default_rng(seed + 999)
            Mp = M[rng.permutation(M.shape[0])]
            v_lvl_shuf = smallest_null(Mp)
            cos_order = abs(cos(v_lvl, v_lvl_shuf))

            # dominant level operator in the null (|loading|)
            dom_idx = int(np.argmax(np.abs(v_lvl)))
            dom_op = OP_NAMES[dom_idx] if dom_idx < len(OP_NAMES) else str(dom_idx)

            # ---- Q2: rate (time-derivative) vs level -------------------
            rates = build_rates(M, dt)
            M_ext_raw = np.column_stack([M[:-1], rates])     # 9 columns
            M_ext, mu_ext, sd_ext = standardize(M_ext_raw)
            v_ext = smallest_null(M_ext)
            v_ext = v_ext / (np.linalg.norm(v_ext) + 1e-12)
            level_frac = float(np.sum(v_ext[:6] ** 2))
            rate_frac = float(np.sum(v_ext[6:] ** 2))
            # which rate op (if any) carries the most among rates
            rate_load = v_ext[6:]
            top_rate_idx = int(np.argmax(np.abs(rate_load)))

            per_seed.append({
                "seed": seed,
                "n_rows": int(M.shape[0]),
                "q1_cos_order_invariance": cos_order,
                "level_null_6d": [round(float(x), 4) for x in v_lvl],
                "level_null_dominant_op": dom_op,
                "q2_level_frac": round(level_frac, 4),
                "q2_rate_frac": round(rate_frac, 4),
                "q2_top_rate_op": RATE_NAMES[top_rate_idx],
                "q2_ext_null_9d": [round(float(x), 4) for x in v_ext],
                "op_means": [round(float(x), 4) for x in M.mean(axis=0)],
            })
            print(
                f"  {name} seed={seed}: "
                f"Q1 cos(order-shuffle)={cos_order:.6f}  "
                f"null~{dom_op}  "
                f"Q2 level_frac={level_frac:.3f} rate_frac={rate_frac:.3f} "
                f"(top rate {RATE_NAMES[top_rate_idx]})",
                flush=True,
            )

        # cross-seed summary
        q1 = [s["q1_cos_order_invariance"] for s in per_seed]
        rf = [s["q2_rate_frac"] for s in per_seed]
        results[name] = {
            "per_seed": per_seed,
            "q1_cos_order_mean": float(np.mean(q1)),
            "q1_cos_order_min": float(np.min(q1)),
            "q2_rate_frac_mean": float(np.mean(rf)),
            "q2_rate_frac_std": float(np.std(rf)),
        }
        print(
            f"  >> {name}: Q1 order-invariance cos "
            f"mean={np.mean(q1):.6f} min={np.min(q1):.6f} | "
            f"Q2 rate_frac {np.mean(rf):.3f} +/- {np.std(rf):.3f}\n",
            flush=True,
        )

    blob = {"config": cfg, "seeds": args.seeds,
            "op_names": OP_NAMES, "rate_names": RATE_NAMES,
            "results": results, "elapsed_s": time.time() - t0}
    with open(args.out, "w") as f:
        json.dump(blob, f, indent=2)
    print(f"[time_flow_probe] wrote {args.out} "
          f"elapsed={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
