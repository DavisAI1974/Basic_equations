# SESSION HANDOFF v11 — 2026-06-02 (Session 11): per-event LIGO + OD consolidation

Read this first, then the Session 11 note in CLAUDE.md (INFO-038). Branch:
work done on `claude/gravity-substrate-config-51cfp`; **main is now synced
to it** (fast-forward) so the SessionStart hook + all data run every
session. The harness-designated `claude/awaiting-files-TgxVq` was NOT used
(continuity precedent); Greg's promised branch-specifying file had not
arrived by session end — confirm branch at next session start.

## What happened this session

1. **Answered Greg's four-force framing questions** (no compute):
   - "Where do the grav/weak laws come from?" -> our OWN toy caricatures
     (retired Session 10), NOT literature force laws. Set 2 (PDG) USES the
     Standard Model (beta functions + gauge group are inputs), so it cannot
     derive the forces or their origins -- it re-derives known running.
     Set 1 (LIGO) is a method/measurement result, not an origins result.
   - "Will these two sets tell us what the forces are / how derived /
     origins?" -> NO. "Origins" is a why/mechanism question (outside OD
     mode), and no dataset contains the origin; the coupling values + gauge
     structure are exactly the unexplained inputs. The honest OD-shaped
     target is the INVERSE PROBLEM (below).
   - "Can we figure out the new physics if it truly calls for it?" ->
     Not its IDENTITY (3 couplings don't invert to a unique particle
     spectrum; nobody can). But YES its REQUIRED FOOTPRINT: extract the
     Delta-b_i (beta-coefficient shifts) + onset scale mu_NP that would
     close the triangle -- a mechanism-agnostic, falsifiable constraint
     surface. FIRST map alternatives (Result Discipline): (a) nothing
     forces single unification [leading deflationary read], (b) two-loop +
     thresholds may shrink the triangle (we used one-loop), (c) the
     extrapolation itself is the conjecture. Gravity is its own extraction
     (power-law -> would need to become log to join; likely confirms it is
     outside the form). THIS IS TRACK B.
   - Greg's decision: do Track A (LIGO) and Track B (new-physics inverse
     problem + SM-parameter-regularity hunt) IN PARALLEL.

2. **LIGO per-event readout, 3 events, NO pooling** (s11_ligo_perevent.py
   -> s11_ligo_perevent_results.json). Greg's directive: run each event
   separately, do not merge or average, the per-event entropy spread may
   mean something. It does:

   | event | H_a(H1) | H_b(L1) | \|H_a-H_b\| | noise cos->attr | peak-MI vs merger | ev/noise MI |
   |-------|---------|---------|-----------|-----------------|-------------------|-------------|
   | GW150914 (loud) | -1.06 | -1.68 | 0.62 | 0.980 | 16.44 vs 16.4 OK | 1.19x |
   | GW170104        | -1.32 | -2.06 | 0.74 | 0.807 | 16.56 vs 16.6 OK | 1.07x |
   | GW151226 (quiet)| -0.71 | -2.32 | 1.61 | 0.816 | 20.56 vs 16.6 NO | 1.06x |

   - Per-event entropy asymmetry |H_a-H_b| VARIES a lot (0.62 / 0.74 /
     1.61) -- a per-event detector-state signature (epoch-specific
     H1-vs-L1 noise-floor/PSD shape surviving whitening), not astrophysics.
   - It INVERSELY tracks how hard noise sits on the equal-entropy attractor
     (low asym -> high cos 0.98; high asym -> low cos 0.82). Internally
     consistent with INFO-036: the (-1,-1,+2) attractor IS the equal-marginal-
     entropy identity, so unequal per-channel entropy => further off it even
     in pure noise. The entropy reading EXPLAINS the cos reading. New ledger
     INFO-038 (isolated, 3 events, no null yet).
   - Detection (MI peak at merger) is loudness-dependent: lands at merger
     for GW150914 + GW170104, misses for quiet long-inspiral GW151226.
   - All event/noise MI ratios are modest (1.06-1.19x); MEANINGLESS without
     a null distribution (built into the batch runner below).

3. **OD consolidation onto main** (Greg: "make sure OD is updated with our
   latest jsons on dipole and 4 forces data"):
   - Brought the Session 6/7/8 four-force + per-domain result JSONs and OD
     stores (store/four_force_caricature, store/simulator_4domain) onto the
     continuity branch (purely additive) and synced main. Medical stores
     (cardiac/cerebro) deliberately LEFT on iZvY4 per Greg's call ("leave
     medical as a separate branch").
   - Four-force JSONs were verified IDENTICAL across O6ahb/iZvY4.

## DECISIONS recorded this session (act on these)

- **Dipole JSONs**: Greg -> "if the info is already in markets pull from
  there. make sure we also have the flow dipole equation listed separately
  too." ACTION (next session): use mcp__claude-code-remote__list_repos +
  add_repo to pull DavisAI1974/Markets (and/or agent), locate the dipole
  result JSON(s) + the markets algebraic dipole eqn
  (H_a^2 = a + b*(H_a*H_b) + c*(H_a*H_b)^2), bring them where OD ingests,
  AND write the "flow dipole equation" out separately as its own artifact.
  (No dipole-named result JSON exists in Basic_equations -- only scripts.)
- **Medical OD stores**: stay on their own branch, NOT main. Done.
- **Branch**: main == continuity now. Confirm the intended working branch
  at next session (Greg's branch-specifying file).

## NEXT SESSION -- ready to run

- **Track A, full 12 (built, not yet run): `s11_ligo_batch.py`.** Runs all
  12 events SEPARATELY (no pooling). Fetches one 4096s/4096Hz bulk file per
  detector per event (~130MB x2 x12 ~= 1.5GB; 29G free) and uses that SAME
  file for (a) the 32s event segment AND (b) ~100 off-source 32s NOISE
  segments -> a per-event NULL DISTRIBUTION giving the merger peak-MI a
  p-value. Canary: `python3 s11_ligo_batch.py --canary` (2 events, N_null=8)
  before the full run. Watch the per-event |H_a-H_b| spread + whether the
  ev/noise MI clears its own null. data/ligo_bulk/ is gitignored.
- **Track B (new-physics inverse problem)**: build the required-running-
  modification extractor on the PDG couplings -- invert the triangle for
  Delta-b_i + mu_NP, with alternatives (a) no-closure, (b) two-loop +
  thresholds, (c) extrapolation-is-conjecture mapped FIRST; plus the
  gravity-into-log extraction. Then the SM-parameter-regularity hunt
  (relations among coupling values / mass ratios / CKM-PMNS angles).
- **Markets dipole pull** (decision above).

All Operating Rules from Sessions 4-7 in force; no new Rule this session.
INFO-038 added (isolated). No PR.
