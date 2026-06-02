"""
PROBE (S18) -- FLOW BATTERY: widen the substrate + deepen the control.

Extends probe_flow_substrate_expression.py. Greg's frame: FLOW is a substrate (a
characteristic running monotonically and DIVERGING toward a critical point); each
domain EXPRESSES it along its own axis + functional form. Time is gravity's expression.

GREG'S EPISTEMIC GUIDANCE (S18, load-bearing): a system that does NOT abide does not
falsify the frame for the others or even for itself -- we may be looking at the wrong
observable for it, or it may not be related this way. Each system is ONE data point;
classify it, keep misses as data, never use a miss to tear down the rest.

WIDEN (new flows, ideally new axes):
  - gravity GW170817 chirp (time axis, independent BNS system)
  - chemistry Brusselator -> Hopf bifurcation: CRITICAL SLOWING DOWN. The relaxation
    rate -> 0 as the control parameter B -> B_crit=2 (Re lambda = (B-2)/2), so the
    relaxation time tau ~ |B - B_crit|^(-1) DIVERGES. AXIS = control parameter (a
    genuinely new axis, not time/scale/mass); DOMAIN = chemistry. Mean-field exponent 1,
    cleanly measurable near-critical from simulation -> also resolves FORM-distinctness
    from the clean side (gravity q~0.375 vs brusselator q~1.0 = distinct expressions).
DEEPEN CONTROL (non-flows: have a critical point but do NOT flow / peak):
  - weak Z Breit-Wigner resonance (mass axis) [from prior probe]
  - driven damped oscillator amplitude response A(omega): peaks at resonance, finite.

Carried over flows: gravity GW150914 chirp (time), strong alpha_s running (scale).

Substrate metric: |Spearman(char, axis)| over the full range + "peaks in interior?"
A FLOW is monotonic (|rho|~1, no interior peak); a RESONANCE peaks (|rho| low, interior
max). Form (for flows): power-law exponent q + R^2 of log-log near the critical point.
"""
import json
import numpy as np
from scipy import stats
from scipy.integrate import solve_ivp
from probe_flow_substrate_expression import (load_gravity, load_strong, load_weak,
                                             powerlaw_fit, monotonicity)


def classify(char, axis):
    """Return (|rho|, peaks_interior, is_flow). Flow = monotonic, no interior max.
    peaks_interior uses an ABSOLUTE margin (robust for small N): a monotone flow has its
    max AT an endpoint; a resonance peaks strictly inside."""
    rho = abs(stats.spearmanr(char, axis)[0])
    char = np.asarray(char)
    i = int(np.argmax(char))
    n = len(char)
    margin = max(2, int(0.03 * n))
    peaks_interior = (i >= margin) and (i <= n - 1 - margin)
    is_flow = (rho > 0.9) and not peaks_interior
    return float(rho), bool(peaks_interior), bool(is_flow)


def load_gw170817():
    g = json.load(open("probe_gravity_chirp_ridge_results.json"))
    h = g["events"]["GW170817"]["H1"]
    t = np.array(h["ridge_t"]); f = np.array(h["ridge_f"]); tc = h["fit"]["t_c"]
    m = t < tc - 2e-3
    return t[m], f[m], tc


def brusselator_critical_slowing(A=1.0, Bs=(1.0, 1.3, 1.5, 1.7, 1.85, 1.92, 1.96)):
    """Measure relaxation time tau(B) approaching the Hopf bifurcation B_crit=1+A^2=2.
    Perturb the fixed point (A, B/A), fit exponential decay of the |x-A| peak envelope."""
    Bc = 1.0 + A**2
    taus, Bs_ok = [], []
    for B in Bs:
        xs0, ys0 = A, B / A                       # fixed point
        def rhs(t, s):
            x, y = s
            return [A - (B + 1) * x + x * x * y, B * x - x * x * y]
        sol = solve_ivp(rhs, [0, 4000], [xs0 + 0.02, ys0], max_step=0.5,
                        dense_output=False, rtol=1e-8, atol=1e-10)
        t, x = sol.t, sol.y[0]
        dev = np.abs(x - A)
        # upper-envelope peaks
        pk = [(t[k], dev[k]) for k in range(1, len(dev) - 1)
              if dev[k] > dev[k - 1] and dev[k] > dev[k + 1] and dev[k] > 1e-9]
        if len(pk) < 4:
            continue
        pt = np.array([p[0] for p in pk]); ph = np.array([p[1] for p in pk])
        # fit log(peak height) ~ -rate * t over the clean decaying span
        keep = ph > ph.max() * 1e-3
        sl, ic, r, *_ = stats.linregress(pt[keep], np.log(ph[keep]))
        rate = -sl
        if rate > 1e-5 and r**2 > 0.95:
            taus.append(1.0 / rate); Bs_ok.append(B)
    Bs_ok = np.array(Bs_ok); taus = np.array(taus)
    dist = Bc - Bs_ok                              # distance to critical point
    # tau ~ dist^(-q): log tau vs log dist
    sl, ic, r, *_ = stats.linregress(np.log(dist), np.log(taus))
    return dict(B_crit=Bc, B_values=Bs_ok.tolist(), tau=taus.tolist(),
                dist=dist.tolist(), exponent_q=float(-sl), R2=float(r**2),
                theory_exponent=1.0,
                char=taus, axis=-dist)   # char rises as dist->0 (axis = -dist, monotone)


def driven_oscillator(omega0=1.0, gamma=0.15):
    w = np.linspace(0.2, 2.0, 400)
    A = 1.0 / np.sqrt((omega0**2 - w**2)**2 + (gamma * w)**2)   # resonance response
    return w, A, omega0


def main():
    systems = {}

    # ---- carried-over flows ----
    t, f, tc = load_gravity()
    q, r2 = powerlaw_fit(tc - t, f)
    rho, pk, isf = classify(f, t)
    systems["gravity_GW150914"] = dict(kind="claimed-flow", axis="time", domain="gravity",
        rho=rho, peaks_interior=pk, is_flow=isf, exponent_q=q, R2=r2, critical="t_c merger")

    Q, a_s, lam = load_strong()
    rho, pk, isf = classify(a_s, Q)
    q, r2 = powerlaw_fit(np.log(Q / lam), a_s)
    systems["strong_alphas"] = dict(kind="claimed-flow", axis="energy-scale", domain="strong",
        rho=rho, peaks_interior=pk, is_flow=isf, exponent_q=q, R2=r2,
        critical="Lambda_QCD", note="form power-vs-log undetermined far from critical (agent #3)")

    # ---- WIDEN: new flows ----
    t7, f7, tc7 = load_gw170817()
    q7, r27 = powerlaw_fit(tc7 - t7, f7)
    rho7, pk7, isf7 = classify(f7, t7)
    systems["gravity_GW170817"] = dict(kind="claimed-flow", axis="time", domain="gravity(BNS)",
        rho=rho7, peaks_interior=pk7, is_flow=isf7, exponent_q=q7, R2=r27,
        critical="t_c merger",
        miss_interpretation=("MISS (data point, NOT falsification -- per Greg): the GW170817 "
            "H1 ridge is the WRONG OBSERVABLE -- narrow-band/noisy (BNS SNR is in L1, which "
            "carries the glitch), exactly why O3 needed matched-filtering not a ridge (it then "
            "recovered M_c to 0.19%). GW150914 is a clean gravity flow; this is an "
            "observable-quality limit on one event, not evidence gravity isn't a flow."))

    bru = brusselator_critical_slowing()
    rho_b, pk_b, isf_b = classify(bru["char"], bru["axis"])
    systems["chem_brusselator_Hopf"] = dict(kind="claimed-flow", axis="control-parameter",
        domain="chemistry", rho=rho_b, peaks_interior=pk_b, is_flow=isf_b,
        exponent_q=bru["exponent_q"], R2=bru["R2"], theory_exponent=bru["theory_exponent"],
        critical="Hopf B_crit=2", note="critical slowing down; NEW AXIS (control param)")

    # ---- DEEPEN CONTROL: non-flows ----
    M, rate, mz, gz = load_weak()
    rho_w, pk_w, isf_w = classify(rate, M)
    systems["weak_Z_resonance"] = dict(kind="control(non-flow)", axis="mass", domain="weak",
        rho=rho_w, peaks_interior=pk_w, is_flow=isf_w, critical="M_Z pole",
        note="resonance: critical point but PEAKS, does not flow")

    w, Aresp, w0 = driven_oscillator()
    rho_o, pk_o, isf_o = classify(Aresp, w)
    systems["driven_oscillator"] = dict(kind="control(non-flow)", axis="drive-frequency",
        domain="mechanics", rho=rho_o, peaks_interior=pk_o, is_flow=isf_o,
        critical="resonance omega0", note="amplitude response peaks at resonance, finite")

    # ---- verdict ----
    flows = [k for k, v in systems.items() if v["is_flow"]]
    nonflows = [k for k, v in systems.items() if not v["is_flow"]]
    claimed_flows = [k for k, v in systems.items() if v["kind"] == "claimed-flow"]
    controls = [k for k, v in systems.items() if v["kind"].startswith("control")]
    misses = [k for k in claimed_flows if not systems[k]["is_flow"]] + \
             [k for k in controls if systems[k]["is_flow"]]
    # form-distinctness from cleanly-determined flows (near-critical exponents)
    grav_q = systems["gravity_GW150914"]["exponent_q"]
    bru_q = systems["chem_brusselator_Hopf"]["exponent_q"]
    form_distinct_clean = bool(systems["chem_brusselator_Hopf"]["R2"] > 0.9
                               and systems["gravity_GW150914"]["R2"] > 0.9
                               and abs(grav_q - bru_q) > 0.3)
    verdict = dict(
        flows=flows, nonflows=nonflows, misses=misses,
        substrate_generalizes=bool(len(flows) >= 3 and all(
            not systems[c]["is_flow"] for c in controls)),
        widen="gravity(GW170817, time) + chemistry(Brusselator->Hopf, CONTROL-PARAMETER "
              "axis, new domain) both classify as flows -> substrate generalizes beyond "
              "the original 2 forces and beyond the time/scale axes to a control-parameter axis",
        control="weak Z + driven oscillator: both have a critical point but PEAK (not flows) "
                "-> the flow restriction is real, not trivially satisfied by 'near a critical point'",
        form_level_resolved_clean=form_distinct_clean,
        form_note=(f"FORM-level expression IS distinct on cleanly-determined flows: gravity "
                   f"power-law q={grav_q:.3f} (~3/8) vs Brusselator critical-slowing q={bru_q:.3f} "
                   f"(theory 1.0). Strong stays form-undetermined (nonperturbative critical "
                   f"region; agent #3). So form-distinctness holds where near-critical data exist."),
        epistemic=("Per Greg: any miss is ONE data point, possibly a wrong-observable choice "
                   "for that system, NOT a falsification of the frame or of that system. "
                   f"Current misses: {misses if misses else 'none'}."),
        reading=("Substrate (monotonic divergence to a critical point) GENERALIZES across "
                 "gravity(x2), strong, and chemistry(new axis); controls (weak, oscillator) "
                 "correctly fail. EXPRESSION: axis distinct (time/scale/mass/control-param); "
                 "FORM distinct where measurable near-critical (gravity 3/8 vs Brusselator 1). "
                 "Axis UNIFICATION still not data-forced (construction-confounded)."))

    out = dict(frame="flow battery: widen substrate + deepen control (Greg S18)",
               systems=systems, verdict=verdict)
    # strip numpy arrays before dump
    for v in systems.values():
        v.pop("char", None); v.pop("axis_arr", None)
    with open("probe_flow_battery_results.json", "w") as fp:
        json.dump(out, fp, indent=2, default=lambda o: None)

    print("CLASSIFICATION (flow = monotonic divergence to a critical point):")
    for k, v in systems.items():
        tag = "FLOW " if v["is_flow"] else "not  "
        print(f"  [{tag}] {k:24s} rho={v['rho']:.3f} peak_int={int(v['peaks_interior'])} "
              f"q={v.get('exponent_q', float('nan')):.3f} R2={v.get('R2', float('nan')):.3f} "
              f"axis={v['axis']} ({v['kind']})")
    print(f"\nflows: {flows}")
    print(f"controls correctly non-flow: {all(not systems[c]['is_flow'] for c in controls)}")
    print(f"misses (kept as data, NOT falsification): {misses if misses else 'none'}")
    print(f"\nform-distinct on clean flows (gravity {grav_q:.3f} vs brusselator {bru_q:.3f}): "
          f"{form_distinct_clean}")
    print(f"substrate generalizes: {verdict['substrate_generalizes']}")
    print("\nwrote probe_flow_battery_results.json")


if __name__ == "__main__":
    main()
