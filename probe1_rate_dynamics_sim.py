"""
PROBE 1 (sim half) -- dMI/dt rate-DYNAMICS baseline across the 4 sim domains.

The gate (probe1_gate_check.py) showed the gravity real-LIGO event window has
only ~8 independent samples -> a lagged rate-dynamics fit is not viable there;
the gravity comparison is deferred to 2nd. The sim domains have ~1200 ordered
ensemble-H time steps, so the rate-dynamics METHOD is fully fittable here.

Question (data-level only, no mechanism): does dMI/dt carry LAGGED / temporal
("memory") structure that the time-blind LEVEL null cannot see? For each domain:
  - INST  : dMI/dt(t) ~ contemporaneous level ops(t)            (level can see this)
  - PAST  : dMI/dt(t) ~ PAST ops(t-1..t-p) only                 (pure memory)
  - LAG   : dMI/dt(t) ~ ops(t) + ops(t-1..t-p)                  (level + memory)
 delta_R2 = R2(LAG) - R2(INST) isolates the memory contribution the level
relation discards. Also report dMI/dt autocorrelation (temporal structure).

Reporting: in-sample R2 (the question is STRUCTURAL -- is lagged structure
present -- not OOS prediction) PLUS a random 5-fold R2 per INFO-026 (random-CV
isolates functional-form fit on non-stationary series; block-CV is catastrophic
here and is NOT the right read for "is the structure present"). 3 seeds, scatter
reported (>=3-seed rule). Speaking posture: I think lagged structure MAY differ
by domain, but the reading waits on the numbers.
"""
import json
import time

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

from per_domain_kbk import (DOMAINS, build_ensemble_operator_matrix)

SEEDS = [11, 22, 33]
P_LAGS = 2          # number of past lags for PAST / LAG models
LEVEL_COLS = [0, 1, 2, 3, 4]   # H_a, H_b, H_a^2, H_b^2, H_a*H_b (MI excluded: it's the target's parent)
MI_COL = 5


def standardize(M):
    return (M - M.mean(0)) / (M.std(0) + 1e-12)


def r2_insample(X, y):
    m = LinearRegression().fit(X, y)
    return float(m.score(X, y))


def r2_randomcv(X, y, k=5, seed=0):
    kf = KFold(n_splits=k, shuffle=True, random_state=seed)
    sc = []
    for tr, te in kf.split(X):
        m = LinearRegression().fit(X[tr], y[tr])
        sc.append(m.score(X[te], y[te]))
    return float(np.mean(sc))


def lag_design(Z, p):
    """Return (X_inst, X_past, X_lag, y_index_slice). Rows aligned to t in [p:]."""
    n = Z.shape[0]
    inst = Z[p:]                                   # ops(t)
    past = np.column_stack([Z[p - k: n - k] for k in range(1, p + 1)])  # ops(t-1..t-p)
    lag = np.column_stack([inst, past])
    return inst, past, lag


def autocorr(x, maxlag=10):
    x = x - x.mean()
    v = np.dot(x, x) + 1e-18
    return [round(float(np.dot(x[:-k], x[k:]) / v), 3) for k in range(1, maxlag + 1)]


def run_domain(name, simfn):
    per_seed = []
    for sd in SEEDS:
        X1, X2 = simfn(sd)
        M, _ = build_ensemble_operator_matrix(X1, X2)
        Z = standardize(M)
        lvl = Z[:, LEVEL_COLS]
        mi = M[:, MI_COL]
        dMI = np.gradient(mi)        # dt constant -> absorbed into coefficients
        dMI = (dMI - dMI.mean()) / (dMI.std() + 1e-12)
        inst, past, lag = lag_design(lvl, P_LAGS)
        y = dMI[P_LAGS:]
        rec = {
            "seed": sd,
            "R2_inst": round(r2_insample(inst, y), 4),
            "R2_past": round(r2_insample(past, y), 4),
            "R2_lag": round(r2_insample(lag, y), 4),
            "R2_inst_cv": round(r2_randomcv(inst, y), 4),
            "R2_lag_cv": round(r2_randomcv(lag, y), 4),
            "acf_dMIdt": autocorr(y, maxlag=8),
        }
        rec["delta_R2_memory"] = round(rec["R2_lag"] - rec["R2_inst"], 4)
        rec["delta_R2_memory_cv"] = round(rec["R2_lag_cv"] - rec["R2_inst_cv"], 4)
        per_seed.append(rec)

    def ms(key):
        v = np.array([r[key] for r in per_seed], float)
        return [round(float(v.mean()), 4), round(float(v.std()), 4)]

    return {
        "domain": name,
        "n_points": int(M.shape[0]),
        "per_seed": per_seed,
        "R2_inst_mean_std": ms("R2_inst"),
        "R2_past_mean_std": ms("R2_past"),
        "R2_lag_mean_std": ms("R2_lag"),
        "delta_R2_memory_mean_std": ms("delta_R2_memory"),
        "delta_R2_memory_cv_mean_std": ms("delta_R2_memory_cv"),
    }


def main():
    t0 = time.time()
    print(f"[probe1_rate_dynamics_sim] p_lags={P_LAGS}, seeds={SEEDS}\n", flush=True)
    results = []
    for name, simfn in DOMAINS.items():
        print(f"--- {name} ---", flush=True)
        r = run_domain(name, simfn)
        results.append(r)
        print(f"   n_pts={r['n_points']}  "
              f"R2_inst={r['R2_inst_mean_std']}  R2_past={r['R2_past_mean_std']}  "
              f"R2_lag={r['R2_lag_mean_std']}", flush=True)
        print(f"   delta_R2(memory) in-sample={r['delta_R2_memory_mean_std']}  "
              f"random-CV={r['delta_R2_memory_cv_mean_std']}", flush=True)
        print(f"   dMI/dt acf (seed {SEEDS[0]}): {r['per_seed'][0]['acf_dMIdt']}\n", flush=True)
    json.dump({"p_lags": P_LAGS, "seeds": SEEDS, "level_cols": LEVEL_COLS,
               "domains": results, "elapsed_s": round(time.time() - t0, 1)},
              open("probe1_rate_dynamics_sim_results.json", "w"), indent=2)
    print(f"[probe1_rate_dynamics_sim] wrote results ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
