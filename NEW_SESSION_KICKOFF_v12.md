# NEW SESSION KICKOFF v12 — opens Session 12 (continues 2026-06-02 Session 11)

## Attach to the new chat
1. `CLAUDE.md` (master context, updated through Session 11)
2. `SESSION_HANDOFF_2026-06-02_v11_results.md` (read this first)
3. this file

## Branch + environment
- **Work on `claude/gravity-substrate-config-51cfp`** (the handoff-specified
  branch; `main` is fast-forwarded to it, 0/0). Confirm at start; do NOT use
  any auto-generated harness branch.
- SessionStart hook is on `main` -> PySR/Julia + numpy/scipy/sklearn
  bootstrap automatically. 29G disk free. GWOSC/PDG/PyPI reachable.
- Keep `main` synced to the work branch each session (workflow runs each
  session). No PR unless Greg asks.

## State in one breath
First real-data four-force work. LIGO method detects mergers via inter-
detector MI; the (-1,-1,+2)/sqrt6 attractor is the equal-marginal-entropy
identity, NOT coupling. PDG: three gauge forces share the running form but
do NOT meet at one point (triangle 1e13-1e17 GeV); gravity power-law,
outside the form. Per-event LIGO (3 events, no pooling) -> INFO-038: entropy
asymmetry |H_a-H_b| is an across-event noise-floor fingerprint; the merger
departure is MI-driven (orthogonal axis). Two tracks open in parallel.

## First actions (ordered)
1. **Track A — full 12-event LIGO, each separate, no pooling.**
   `python3 s11_ligo_batch.py --canary`  (2 events, N_null=8) as the canary,
   then `python3 s11_ligo_batch.py` (all 12, ~1.5GB bulk download, per-event
   off-source null -> p-value for the merger MI peak). Watch: per-event
   |H_a-H_b| spread; whether each event's ev/noise MI clears its OWN null;
   and the open INFO-038 thread (no-MI basis showed noise-OFF / event-ON
   attractor reversal on GW150914 — does it hold across the 12?).
2. **Track B — new-physics inverse problem** on the PDG couplings. Build the
   required-running-modification extractor: invert the triangle for the
   beta-coefficient shifts Delta-b_i + onset scale mu_NP that would close it
   (mechanism-agnostic, falsifiable). MAP ALTERNATIVES FIRST (Result
   Discipline): (a) nothing forces single unification [leading deflationary
   read]; (b) two-loop + thresholds may shrink the triangle (Session 10 used
   one-loop); (c) the extrapolation itself is conjecture. Then the gravity-
   into-log extraction (how violent a modification to make gravity join the
   gauge form). We can recover new physics's required FOOTPRINT, never its
   IDENTITY. Then the SM-parameter-regularity hunt (relations among coupling
   values / mass ratios / CKM-PMNS angles).
3. **Markets dipole pull** (Greg's call): mcp__claude-code-remote__list_repos
   + add_repo on DavisAI1974/Markets (and/or agent); locate the dipole result
   JSON(s) + the markets algebraic dipole eqn
   (H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2); bring where OD ingests AND write
   the "flow dipole equation" out separately as its own artifact. (No
   dipole-named result JSON exists in Basic_equations — only scripts.)

## Standing constraints
- Do NOT overwrite the Markets section in CLAUDE.md (placeholder, lines ~93).
- Medical OD stores (cardiac/cerebro) stay on their own branch, not main.
- All Operating Rules (Sessions 4-7) in force: no pre-assigned meaning,
  probe-not-falsifier, speaking posture before+after, Rule D incomplete-not-
  wrong, they-never-stacked, treat-literature-as-conjecture. Report data /
  interpretation / frame separately. No tent-widening on outliers. >=3 seeds
  + scatter for any operator-space spatial claim. Frames never grade
  themselves.
