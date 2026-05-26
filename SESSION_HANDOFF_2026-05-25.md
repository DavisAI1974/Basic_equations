# Information Layer Session Handoff — May 25, 2026

## Where we ended
Discovered the algebraic dipole. The differential dipole found in March
is its time-derivative shadow. Each science displays the constraint in
its own substrate-native form. Level 2 (the floor below dipole) does
not appear as a rank-deficient operator cloud — may be statistical
rather than algebraic, may need a different basis, or may not exist
as a deeper layer at all (the "rendering engine" reframe).

## Confirmed results
- Chemistry algebraic equation: H_a^2 = 0.007 - 0.093*(H_a*H_b) + 1.309*(H_a*H_b)^2, R^2 = 0.943
- Geology rank-3 constraint surface (6D operator cloud collapses to 3D); tightest constraint:
  0.724*(H_a*H_b) - 0.441*H_b^2 - 0.290*H_a^2 ~ 0, std/mean = 0.15%
- Biology in MI-space: MI ~ polynomial(H_a) with H_b absent. R^2 = 0.66. Carrier subsystem governs MI.
- Physics weakest signal in position-only space; phase-space test pending.
- Universal: H_a^2 = quadratic-or-cubic function of H_a*H_b across 4/4 sciences in 2D,
  all multi-variable extensions improve fit except chemistry (basis issue).

## Surviving prior-work anchors
- ΔT framework was FALSIFIED. Absolute T controls decoherence, not gradient. Reverted.
- January result that survives: at low T, fluctuation theorem makes time arrow weak;
  at high T, arrow strongly enforced.

## Experiments queued (in build order)
1. Phase-space physics — split into H_pos, H_mom, MI(pos,mom). Find dipole in actual conjugate phase space.
2. Chemistry reaction-order test — first-order chemistry (no X^2*Y). If algebraic dipole drops
   from quadratic to linear, form reads reaction order directly.
3. Geology dimensionality test — 0D, 1D, 2D variants. Does SVD rank track spatial dimensionality?
4. T_eff stratified differential dipole — split simulation into low-T_eff vs high-T_eff regions,
   compare differential R^2. Predicted: high T_eff = stronger dipole, low T_eff = weaker
   (arrow weakens).
5. State-dependent time correction — recompute dH/dt with dτ = dt * f(local state). Predicted:
   differential R^2 jumps toward algebraic R^2.
6. Pairwise Level 2 — instead of network-wide SVD, test algebraic relationship between dipole
   pairs (B_i vs B_j).

## Speculative frames (NOT discoveries — keep separate)
- "Our world is the display, information is the flow" — frame, not result
- "Time has aspects we collapse into scalar dt" — frame, not result
- These motivate experiments but are not claims.

## Files
E:\information_layer\ (and mirror F:\Factory\knowledge\information_layer\):
- level2_four_sciences.py (Level 2 first attempt, R^2=0.02-0.13)
- dna_dipole_canary.py
- static_dipole_test.py (the reframe that worked)
- triple_followup.py
- nonlinear_deep_dive.py
- three_discoveries.py
- SESSION_HANDOFF_2026-05-25.md (this file)

## MASTER_DISCOVERIES.json drafts (not yet added)
- INFO-008  Algebraic dipole reformulation
- INFO-008a Chemistry quadratic equation
- INFO-008b Geology rank-3 constraint surface
- INFO-008c Biology asymmetric MI law
- INFO-009  Level 2 algebraic absent at network scale (OPEN, not closed)

## Prompt for next session
"Resuming information layer work. Read SESSION_HANDOFF_2026-05-25.md
first. Then build experiment 1 from the queue (phase-space physics).
Falsification-first. No metaphysics. Use the static algebraic
framework, not differential."
