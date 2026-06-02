# SESSION HANDOFF — 2026-06-02 (S17) — CLAUDE.md drift fix + construction-vs-nature probe

Branch `claude/claude-md-problem-nM8Hm` (main untouched; no PR). All Operating Rules
in force. Two things this session: JOB 1 (the CLAUDE.md drift) then backlog #1 (the
construction-vs-nature probe Greg flagged "clean, decisive").

## JOB 1 — CLAUDE.md drift FIXED

The repo `CLAUDE.md` body was canonical only through S11 with a header stale at "S7";
the true master was S14 (with S15/S16 in handoffs). This branch was also cut from
`main` at the S12 kickoff, so it was missing all S13-S16 work.

- Fast-forwarded this branch to absorb the S16 work branch
  (`claude/claude-md-strategy-pwfZr`, 10 commits) — pulled in `CLAUDE_master_through_S14.md`,
  the S15 note, the S14-S15 combined handoff, the S16 handoff, the backlog, the S17
  kickoff, plus all S13-S16 scripts + result JSONs.
- Verified the S14 master is a STRICT SUPERSET of the S11 body before overwriting
  (diffed section headers + ledger entries: nothing dropped; all 44 INFO entries
  preserved). The S14 master already carried S13/S14/S15 notes.
- Rebuilt `CLAUDE.md` from the S14 master, folded in a Session 16 note + (this
  session) a Session 17 note, re-attached the START-HERE workflow block, fixed the
  header, and updated the Session Handoff Pointer.
- Kept `CLAUDE_master_through_S14.md` + the S15/S16 handoffs in-repo as historical
  artifacts. The "keep the header current" Operating Rule (added S14) is what stops
  the drift recurring.

## PROBE (backlog #1) — construction-vs-nature of the equal-entropy substrate

THE question (Greg, S16): everything we measure sits ON the equal-entropy region; is
equal-marginal-entropy a FORCED property of physical 2-channel observables (nature),
or our construction choice (bookkeeping)?

THE LEVER (`probe_construction_vs_nature.py`): scale ONE channel by a constant s,
`b -> s*b`. This is (a) MI-INVARIANT — `mi_hist_2d` uses adaptive (data-range) bin
edges, so scaling rescales the edges proportionally and the 2D histogram is
IDENTICAL, MI(a, s*b) = MI(a, b) exactly; the coupling/physics is untouched — and
(b) a pure UNITS choice — differential entropy is not scale-invariant,
`H(s*b) = H(b) + ln|s|`. So scaling injects pure construction asymmetry into the
marginal entropies while leaving the genuine coupling invariant. Reuses
`kbk_pipeline` (vasicek H + hist MI + extract_v1 + project_234) verbatim.

Design: 4 heterogeneous sim generators (OU correlated, AR(1)-Laplace, coupled
logistic chaos, sine+noise) x 5 seeds [11,22,33,44,55] x scale sweep
[0.1,0.2,0.5,1,2,5,10]; PLUS 3 real LIGO segments scaled post-whitening
(GW150914 H1/L1 aligned per s10; two independent noise-only segments V1/V2). Metric:
mean|H_a-H_b| (asymmetry, the INFO-038 knob), cos(v_null_234,(1,1,2)/sqrt6) [=
project_234, the s10/s11 substrate metric], cos to (-1,-1,2)/sqrt6 (the -(H_a-H_b)^2
identity), |MI coef in null|, mean MI.

### Results (data level)

1. **MI is EXACTLY scale-invariant.** MI_cv (coeff of var of mean MI across the
   whole scale sweep) = 1e-16 to 9e-16 (machine precision) for ALL 7 systems. The
   genuine coupling is representation-independent — a units choice changes no physics.

2. **Marginal-entropy asymmetry is a pure units knob.** asym_range (max-min of
   mean|H_a-H_b| over the scale sweep) = 2.1 to 3.8 for 6 of 7 systems, tracking
   `|H_a - H_b0 - ln s|` exactly: symmetric-construction sims (baseline H_a~=H_b)
   give a monotone `~ln s` rise; LIGO segments give the predicted V-shape with the
   minimum where `ln s` matches the baseline gap (V1 asym 0.78->0.32->1.52, min near
   s~3; GW150914 min near s~2). You can place any system on or off the
   equal-marginal-entropy point by choosing units for one channel.

3. **The substrate cos metric is NOT scale-invariant** (cos112_std over the sweep
   0.07-0.27; cosident_std 0.14-0.28) and is fragile/segment-specific. At s=1 the
   independent noise segments V1/V2 do NOT reproduce s10's GW150914-noise cos 0.984
   (they sit at 0.12 / 0.38). NOTE on the s10 "0.984": that null is DOMINATED by the
   LINEAR H_a,H_b terms (null6d [+0.569,+0.774,...]); project_234 keeps only the
   small quadratic residual, whose projection onto (1,1,2) happens to be ~0.98 for
   that particular near-constant-entropy (whitened) segment. It is an estimator-noise
   residual, not a robust geometric fact. cos generally peaks near the
   asymmetry-minimum and degrades as you rescale away.

### Reading (Result Discipline: data / interpretation / frame, alternatives mapped)

- DATA finding (INFO-051, located): MI is exactly scale-invariant; marginal-entropy
  asymmetry is set by a per-channel units choice (`H(sX)=H(X)+ln s`); the substrate
  cos metric is representation-dependent and fragile.
- INTERPRETATION (deflationary, supported): "everything falls on the equal-entropy
  substrate" is BOOKKEEPING. Systems sit there because we build/normalize channels
  to comparable scales (sims: equal noise amplitude; LIGO: whitening to unit
  variance). A physics-preserving units change (MI fixed exactly) dissolves the
  clustering. The genuinely physical, scale-invariant content in the basis is MI.
- ALTERNATIVES mapped + ruled against: (1) a hidden scale-invariant substrate
  signature — the only scale-invariant quantity in the basis is MI, which is not the
  entropy-quadratic "substrate", so no; (2) nature enforces equal scales for real
  channel pairs — for LIGO that is instrument-engineering + whitening (a construction
  choice), not a law. Both reduce to construction.
- This RESOLVES the S16 backlog-#1 open question on the construction side and is
  internally consistent with INFO-036 (attractor = equal-marginal-entropy identity)
  and INFO-038 (asymmetry moves objects off): asymmetry moves objects off BECAUSE it
  is a units knob.

### Distinct regime kept as data (don't predetermine good/bad -- Greg, S17)

`logistic_chaos` is the lone system with small asym_range (0.23). Cause: coupled
logistic maps (r=3.9, eps=0.15) ANTI-synchronize (corr -0.995) -- the two channels
become near-perfect mirror images. In that LOCKED regime the units-knob barely moves
the asymmetry (ln s shift +0.076 vs expected +2.303 at s=10), unlike every
independent system. Two things are jointly true and both are data: (1) the windowed
Vasicek estimator hits its limit on a near-deterministic signal (entropy ~ -620,
dominated by near-zero-diff windows), and (2) strongly-coupled / locked channels
genuinely respond differently to rescaling than independent ones. This is NOT garbage
to exclude -- it is a real "what happens at near-perfect coupling" data point that may
deserve its own probe (e.g. a different entropy estimator, or treating perfect
coupling as a regime of interest). Logged, not predetermined as bad.

### Caveats / open
- The quadratic-null cos metric (project_234) is estimator-noise-sensitive
  (INFO-024) and segment-specific; do not over-read individual cos values.
- Speaking posture (after): the data supports the bookkeeping reading decisively at
  the marginal-entropy level (MI exactly invariant; asymmetry = units knob). The
  fragile substrate-cos level is consistent with it but not load-bearing.

## Files this session
- `CLAUDE.md` (rebuilt canonical through S17), `SESSION_HANDOFF_2026-06-02_S17.md`,
  `BACKLOG_tests_and_probes.md` (#1 marked DONE).
- `probe_construction_vs_nature.py` + `probe_construction_vs_nature_results.json`
  (+ `_canary.json`).
- h5py confirmed in `requirements.txt` (installed at runtime; the fresh container's
  SessionStart hook had not installed it this run).

## Next (backlog, Greg picks order; clear backlog before new probes)
Remaining strong candidates: #2 construction control gravity-LIGO vs EM-HBT flow
(needs Zenodo 5113016 .rar + S13 HBT loader); #3 proper chirp recovery
(Q-transform/matched-filter); #7 SF femtoscopy-R; plus housekeeping #11 (OD /
MASTER_DISCOVERIES) and the now-mostly-done #12 (master merge).
