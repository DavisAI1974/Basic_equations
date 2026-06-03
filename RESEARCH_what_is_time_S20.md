# RESEARCH MAP — "What is time?" (Prong A, S20 research agents)

Two background research agents were sent out per the backlog's WHAT IS TIME / Prong A
directive (Greg, S19). Literature treated as conjecture to characterize, NOT cited as
support. What is OURS vs KNOWN is flagged throughout.

---

## PART 1 — Map of physics positions on time (+ empirical distinguishers)

| Position | Core claim | Empirical distinguisher | Data-shaped prediction? |
|---|---|---|---|
| 1. Thermodynamic arrow / entropy production | Time's direction = entropy increase; micro-arrow = asymmetry of trajectory probabilities | Forward vs time-reversed trajectory prob. ratio; KL-div >= entropy production | **YES** — fluctuation theorems <e^-S>=1, P(+S)/P(-S)=e^S on measured small-system time-series |
| 2a. Rovelli thermal/relational time | Time = modular (Tomita-Takesaki) flow of the equilibrium state; no preferred time | Physical clock singled out by the state, not chosen | MOSTLY NO — built to reproduce standard dynamics |
| 2b. Barbour timeless / shape dynamics | No fundamental time; static config space ("Platonia"); flow from records | Differs from GR only in conformal/cosmological regimes | WEAK YES — unconfirmed |
| 3. Block universe vs presentism | All times equally real vs only present real; relativity-of-simultaneity favors block | None clean — shared relativistic input | NO — metaphysical, identical predictions |
| 4. Proper time / GR | Time = invariant metric interval; clocks tick proper time | Two clocks at different potential/velocity accumulate different elapsed tau | **YES (confirmed)** — dtau/tau = dPhi/c^2 (Pound-Rebka, GPA, optical clocks, GPS) |
| 5. QM time (Page-Wootters / Pauli-POVM) | Time emergent from clock-system entanglement; no self-adjoint time op -> POVM time observables | External sees static, internal sees evolution; POVM time-of-arrival distributions | YES (proof-of-principle) — Moreva two-photon clock; arrival-time distributions |
| 6. Tolman-Ehrenfest | T*sqrt(g00)=const; equilibrium has a temperature gradient in gravity; beta ~ dtau/c | Spatial T-gradient at equilibrium in a gravitational field tracking time-dilation gradient | YES (in principle, tiny) — grad T / T = g/c^2 (~1e-16 /m on Earth, below sensitivity) |

**Strongest data-shaped positions today:** #1 (entropy production / fluctuation theorems)
and #4 (proper time). Both make quantitative, confirmed predictions on real time-series /
clock data — the cleanest "data is the output" character.

**Cross-cutting thread relevant to our operator/entropy work:** positions #1, #5, #6 all
reduce "time" to a relation between an entropy/temperature/clock degree of freedom and the
rest of a system — time as a function of marginal/relational statistics. That is the cluster
most naturally expressible in a windowed-entropy / MI operator basis.

Key refs: Lebowitz (Boltzmann entropy); Parrondo et al. arXiv:0904.1573; Seif/Hafezi/
Jarzynski Nat.Phys. 2020; Rovelli arXiv:0903.3832 (Forget Time); Connes-Rovelli thermal time;
Barbour arXiv:gr-qc/0103055; Rietdijk-Putnam (block universe); Zych/Brukner arXiv:1105.4531
(proper-time interferometric witness); Page-Wootters / Moreva et al.; Galapon arXiv:quant-ph/
9908033 (Pauli's theorem); Rovelli-Smerlak arXiv:1005.2985 (Tolman-Ehrenfest "speed of time").

---

## PART 2 — Is "time = an information / dipole expression" already published? (novelty scan)

**NOVELTY VERDICT:** The "time-from-information" *genus* is published, but our SPECIFIC
identification is NOT found explicitly.

PUBLISHED (therefore NOT ours to claim novel):
- Time emerging from a **clock<->system bipartition** — Page-Wootters (Moreva et al. PRA 89,
  052122; arXiv:1310.4691). The two-subsystem architecture is established.
- Time parametrized by a **monotonic information/entropy** quantity — entanglement entropy,
  modular flow, coarse-grained entropy ("Entropy as a Clock", IJTP 2025/2026); thermal time
  (Connes-Rovelli 1994).
- **Mutual information / transfer entropy between two subsystems** as a measure of
  correlation, synchronization, information flow — incl. conserved-MI-in-a-causal-loop
  ("Subtime", arXiv:2603.11571, 2026).

APPARENTLY OURS-TO-CLAIM-NOVEL (pending a close read of 3 near-neighbors first):
- Identifying the **rate/expression of time specifically with the two-channel mutual
  information of a dipole** — MI as *the* time observable, not as an order parameter,
  conserved quantity, or generic correlation measure.
- The **dipole / paired-channel primitive** as the carrier of the time expression (literature
  uses single-region entanglement entropy, single-algebra modular flow, or PaW pointer
  eigenvalues).

CAVEAT (read before any novelty claim): Subtime (2603.11571, 2026), Entropy as a Clock
(IJTP 2025), Entropic time endowed in quantum correlations (arXiv:1111.2855, 2011) are the
three that could collapse the gap. Skimmed via abstract/HTML only.

---

## PART 3 — Clock-rate UNIVERSALITY: the ONE law our positive findings instantiate

All three positive time findings are projections of one weak-field master relation:

    dtau/dt = 1 + Phi/c^2 - v^2/(2 c^2)      (Ashby, Living Rev. Relativity 6, 1, 2003)

| System | Observable | GR functional form | Predicted value |
|---|---|---|---|
| GPS eccentricity (RINEX) | periodic clock correction | dt_r = -(2/c^2)(r.v) = -2 sqrt(GM a) e sinE / c^2 (Ashby eqs 39-40) | coefficient -2/c^2 on (r.v) |
| Grav. redshift (Galileo eccentric-sat test) | fractional freq shift | df/f = dPhi/c^2 (Phi leg, v fixed) | coeff +1; Delva/Herrmann 2018 confirm to (0.19+/-2.48)e-5 |
| Pulsar B1913+16 Einstein delay | amplitude gamma over orbit | gamma = e (Pb/2pi)^1/3 (G/c^2)^2/3 * mc(mp+2mc)/(mp+mc)^4/3 | gamma = 0.004307(4) s (Weisberg-Huang 2016) |
| Pulsar B1913+16 orbital decay (cross-check) | Pb_dot (quadrupole, NOT clock-rate) | Pb_dot ~ -(192pi/5)(...)^5/3 f(e) | ratio meas/GR 0.9983 +/- 0.0016 |
| LIGO inspiral proper time | accumulated PN phase | integrated master relation | (chirp-mass recovery from strain) |

KNOWN: that GPS eccentricity, gravitational redshift, and pulsar Einstein delay are all
instances of dtau/dt = 1 + Phi/c^2 - v^2/(2c^2) is textbook GR.
OURS: independently *recovering* these coefficients from raw data via the OD/operator method
WITHOUT assuming the GR form, and checking the three datasets land on the single master
relation. (See probe_clock_rate_universality.py — Prong B #1 result, S20.)
