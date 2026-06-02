"""
S10 item (1): is the SUBSTRATE-level signature more invariant than the
EXPRESSION-level family across the mapping-campaign knob sweep?

Session 8's mapping campaign (INFO-030) found the per-domain MI-vs-H
EXPRESSION family mutates qualitatively as one knob moves (biology
exp -> linear at high beta; chemistry linear -> ratio at low B; physics
quadratic-difference -> single-channel quadratic away from K=0.20).
That is the expression level. The frame that survived Session 9
(Base-of-Structure: a simple shared substrate, with the differences
being how each system leaves the base under perturbation) predicts the
SUBSTRATE level should be MORE invariant under the same knob.

This script tests that directly. It reuses the ORIGINAL parametrized
simulators + knob sweeps from mapping_campaign.py (verbatim, pulled from
the Session 8 branch -- not reconstructed). For each domain x knob value
x seed it builds the SAME ensemble-H operator matrix as INFO-023/030 and
extracts the SUBSTRATE signature:
  - extract_v1 null direction (6D)          -- the substrate constraint
  - KBK null-subspace rank                  -- substrate dimensionality
  - cos of the [2,3,4] projection to the (+1,+1,+2)/sqrt(6) attractor
Then per domain it measures how much the substrate DIRECTION drifts
across knob values: |cos(6D null)| between knob-value pairs (baseline
vs low, baseline vs high, low vs high), mean +/- std across seeds. High
|cos| (near 1) => substrate invariant under the knob; low |cos| =>
substrate drifts as much as the expression does.

Speaking posture (Rule C, BEFORE running): I THINK the substrate
direction may hold steadier across the knob than the expression family
did, which would support the Base-of-Structure reading. But we wait on
what the data says and where it points -- a domain whose substrate ALSO
drifts as much as its expression would be just as informative, and would
push back on "substrate = simple shared base." No verdict in advance.

Canary: 1 seed, N_ens=200, T=10, knob endpoints only.
Full:   3 seeds, N_ens=600, T=30, all 3 knob values.
"""
import sys
import json
import time
import itertools
import numpy as np

from mapping_campaign import KNOB_SWEEPS
from per_domain_kbk import build_ensemble_operator_matrix
from kbk_pipeline import extract_v1, kbk_rank_gap, project_234

CANARY = "--canary" in sys.argv
seeds = (11,) if CANARY else (11, 22, 33)
cfg = dict(N_ens=200, T=10.0, dt=0.02) if CANARY else dict(N_ens=600, T=30.0, dt=0.02)


def substrate_signature(X1, X2):
    M, _ = build_ensemble_operator_matrix(X1, X2)
    v_null, S, _, _ = extract_v1(M)
    kbk = kbk_rank_gap(S)
    _, cos_attr = project_234(v_null)
    return v_null, int(kbk["rank_null_subspace_kbk"]), float(cos_attr)


def acos_pair(u, v):
    """Sign-invariant cosine between two directions."""
    return float(abs(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-30))


def run_domain(domain, spec):
    knob = spec["knob"]
    simfn = spec["simfn"]
    baseline = spec["baseline"]
    values = spec["values"]
    if CANARY:
        values = [values[0], values[-1]]

    # per knob value: substrate direction per seed, rank, cos-to-attractor
    dirs = {v: [] for v in values}     # list of 6D null dirs across seeds
    ranks = {v: [] for v in values}
    cos_attr = {v: [] for v in values}
    for v in values:
        for s in seeds:
            X1, X2 = simfn(seed=s, **{knob: v}, **cfg)
            vec, rank, ca = substrate_signature(X1, X2)
            dirs[v].append(vec)
            ranks[v].append(rank)
            cos_attr[v].append(ca)

    # seed-stability of the substrate direction at each knob value
    seed_stab = {}
    for v in values:
        if len(seeds) > 1:
            cc = [acos_pair(a, b) for a, b in itertools.combinations(dirs[v], 2)]
            seed_stab[v] = (float(np.mean(cc)), float(np.std(cc)))
        else:
            seed_stab[v] = (None, None)

    # cross-knob drift of the substrate direction: |cos| between knob
    # pairs, matched per seed then averaged
    cross_knob = {}
    for v1, v2 in itertools.combinations(values, 2):
        cc = [acos_pair(dirs[v1][i], dirs[v2][i]) for i in range(len(seeds))]
        cross_knob[f"{v1}|{v2}"] = (float(np.mean(cc)), float(np.std(cc)))

    print(f"\n=== {domain}  (knob {knob}, baseline {baseline}) ===")
    for v in values:
        ss = seed_stab[v]
        ca = np.mean(cos_attr[v])
        rk = ranks[v]
        ss_str = f"{ss[0]:.3f}" if ss[0] is not None else "n/a"
        print(f"  {knob}={v:<6}  null_rank={rk}  cos_to_attractor={ca:+.3f}  "
              f"seed-stability|cos|={ss_str}")
    print(f"  cross-knob substrate |cos(6D null)|:")
    for k, (m, sd) in cross_knob.items():
        print(f"    {knob} {k:<14}  |cos|={m:.3f} +/- {sd:.3f}")

    return dict(domain=domain, knob=knob, baseline=baseline,
                values=[float(v) for v in values],
                ranks={str(v): ranks[v] for v in values},
                cos_to_attractor={str(v): float(np.mean(cos_attr[v])) for v in values},
                seed_stability={str(v): list(seed_stab[v]) for v in values},
                cross_knob_cos={k: list(v) for k, v in cross_knob.items()})


t0 = time.time()
print("=" * 100)
print(f"S10 SUBSTRATE-INVARIANCE (item 1) {'[CANARY]' if CANARY else '[FULL]'} "
      f"seeds={seeds} cfg={cfg}")
print("Question: does the SUBSTRATE null direction stay invariant across the same knob")
print("that mutates the EXPRESSION family (INFO-030)?  High cross-knob |cos| => invariant.")
print("=" * 100)

results = []
for domain, spec in KNOB_SWEEPS.items():
    results.append(run_domain(domain, spec))

# headline summary: min cross-knob |cos| per domain (worst-case substrate drift)
print("\n" + "=" * 100)
print("SUMMARY -- worst-case substrate direction stability across the knob sweep:")
for r in results:
    mins = min(v[0] for v in r["cross_knob_cos"].values())
    print(f"  {r['domain']:20s} min cross-knob |cos| = {mins:.3f}  "
          f"(ranks {sorted(set(sum(r['ranks'].values(), [])))})")
print("=" * 100)

out = dict(meta=dict(canary=CANARY, seeds=list(seeds), cfg=cfg,
                     elapsed_s=round(time.time() - t0, 1)),
           note="cross_knob_cos = |cos| of 6D substrate null between knob "
                "values (mean,std over seeds). Compare against INFO-030 "
                "expression-family mutation across the SAME knobs.",
           results=results)
fn = "s10_substrate_invariance_canary.json" if CANARY else "s10_substrate_invariance_results.json"
json.dump(out, open(fn, "w"), indent=2)
print(f"\nWrote {fn} ({out['meta']['elapsed_s']}s)")
