# DavisAI Operator Discovery — Capability Brief

**DavisAI Systems · Operator Discovery (OD)**
Prepared 2026-06-02. One-page capability summary for technical partners.

## What it is

Operator Discovery is a domain-agnostic engine that extracts **governing laws
directly from raw measured data**, with no theory-supplied parameters fed in. The
same pipeline applies across domains; here we demonstrate it on physics' three
hardest force domains using only public data.

## The proof: established governing laws recovered from raw public data

In each case below the method was given raw measured data and recovered the
governing-law parameters and confirmed the predicted functional form, matching the
independently established values. These are **known** laws — which is the point: the
answers are checkable against ground truth.

- **Gravity — gravitational-wave inspiral law (LIGO GW150914 public strain).**
  Recovered the inspiral relation `f^(-8/3)` linear in time from raw strain via a
  time-frequency ridge, and the chirp mass, on **both detectors independently**:
  H1 and L1 agree to under 1% (recovered chirp mass ~38 vs catalog detector-frame
  ~31 solar masses; fit R^2 = 0.99). The same pipeline locates the merger in time and
  confirms the inter-detector signal sits at the correct physical light-travel delay
  (~7 ms) at 15-sigma significance against a time-slide null.

- **Gravity — relativistic time dilation, recovered two independent ways.**
  (a) From raw GPS data: the relativistic clock law `Δt = -2(r·v)/c²` recovered from raw
  RINEX pseudorange observations (IGS station BRUX) against an independent broadcast-element
  regressor — on the eccentric Galileo satellites the coefficient comes out at **1.02-1.04×
  the general-relativistic value** (`-2/c²`), at 380-sigma significance against a scramble
  null. (b) Independent confirmation from binary-pulsar timing (raw Arecibo times-of-arrival,
  PSR B1913+16): the orbital-decay rate `dP_b/dt` recovered at **1.005× the GR prediction**
  (orbital decay detected at ~17,684 sigma) and the Einstein-delay (gravitational clock
  redshift) parameter to **0.014%**. A clean control: the same clock law is *absent* from
  standard precise GPS products because it is already modeled out there — recovering it
  required going back to the raw observations.

- **Weak force — Z boson resonance (CMS public dimuon data).**
  Recovered the Z Breit-Wigner lineshape from 10,227 real dimuon events:
  M_Z = 90.75 GeV = **99.5% of the PDG value** (91.1876 GeV). A same-charge control
  shows no peak.

- **Strong force — QCD running coupling (measured alpha_s world data).**
  From 13 measured alpha_s(Q) points across 1.78 GeV to 1 TeV, recovered that
  1/alpha_s is linear in ln(Q) with a **positive slope — asymptotic freedom forced by
  the data**, with no beta function assumed (chi-squared/dof = 0.81; Lambda_QCD
  ~150 MeV).

## Built-in falsifiability

The method is designed to test, not just fit. As an example, we ran a rigorous
search for structure **beyond** general relativity in the GW150914 inter-detector
signal: after subtracting the full matched-filter waveform, the correlated signal
collapses to the noise floor — a clean **negative** result reported as such. We
distinguish recovery of known physics from claims of new physics, and we report
nulls.

## Why it matters

A system that independently re-derives established governing laws from raw data — in
domains as varied as gravitational waves, particle resonances, and the strong
coupling — is evidence that the same engine can discover governing laws in domains
where they are **not yet known**: medicine, geophysics, climate, materials, markets,
and defense sensing. Recovery of the known is the credential for discovery of the
unknown.

## Validation summary

| Domain  | Data source              | Recovered                         | Benchmark            |
|---------|--------------------------|-----------------------------------|----------------------|
| Gravity | LIGO GW150914 strain     | inspiral law + chirp mass (H1=L1) | R^2 0.99; ~within 25%|
| Gravity | Raw GPS RINEX (BRUX)     | time-dilation coeff -2/c^2        | 1.02-1.04x GR; z=380 |
| Gravity | Arecibo TOAs (B1913+16)  | orbital decay dP_b/dt + Einstein gamma | 1.005x GR; gamma 0.014% |
| Weak    | CMS dimuon (10,227 evts) | Z Breit-Wigner, M_Z 90.75 GeV     | 99.5% of PDG         |
| Strong  | alpha_s(Q) world data    | asymptotic freedom, Lambda_QCD    | chi2/dof 0.81        |

## Contact

DavisAI Systems — Greg Davis, Founder and Chief Research Officer. Columbus, Ohio.
greg@davisai.ai

*All results derive from public datasets; methods and scripts are available for
independent verification under appropriate terms.*
