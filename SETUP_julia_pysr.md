# Julia + PySR per-session setup (Claude Code on the web)

How Julia + PySR come up automatically every session in this repo, how to verify
it, and how to fix the other session that's having problems.

## TL;DR
A **SessionStart hook** installs the Python stack + PySR and bootstraps Julia on
every fresh web container. The other session is almost certainly failing because
**these files are not on the branch that session starts from** (the hook only runs
if it's present on the starting branch). Copy the three files below into that repo's
branch and the problem goes away.

## The mechanism (3 files, must be on the session's STARTING branch)

1. `.claude/settings.json` -- registers the hook:
```json
{
  "hooks": {
    "SessionStart": [
      { "hooks": [ { "type": "command",
        "command": "$CLAUDE_PROJECT_DIR/.claude/hooks/session-start.sh" } ] }
    ]
  }
}
```

2. `.claude/hooks/session-start.sh` -- the bootstrap (must be `chmod +x`):
   - guards on `CLAUDE_CODE_REMOTE == "true"` (web/remote only; no-op on local).
   - `pip install --upgrade numpy scipy scikit-learn` (core, required).
   - `pip install pysr` then `python -c "import pysr; from pysr import PySRRegressor"`.
     **That import is the bootstrap** -- it triggers `juliapkg` to download Julia
     (currently 1.11.9) into `/root/.julia/environments/pyjuliapkg/` and precompile
     `SymbolicRegression.jl`. PySR is **best-effort**: a blocked Julia download warns
     but does NOT fail session start (the pure-numpy scripts still run).
   - `set -euo pipefail`, and the `CLAUDE_PROJECT_DIR`/`CLAUDE_ENV_FILE` PYTHONPATH
     append is **guarded** for unset vars (an unset var under `set -u` crashed the
     hook before S10 -- keep the guard).

3. `requirements.txt`:
```
numpy
scipy
scikit-learn
pysr
h5py
```

## Verify it worked (run these in a new session)
```bash
echo "CLAUDE_CODE_REMOTE=$CLAUDE_CODE_REMOTE"          # must be: true
ls -l .claude/hooks/session-start.sh                   # must be executable (-rwxr-x)
python3 -c "import pysr; print(pysr.__version__)"       # e.g. 1.5.10
python3 -c "from pysr import PySRRegressor; print('Julia backend reachable')"
```
The hook prints `[session-start] PySR + Julia backend ready` and
`[session-start] numeric stack ready` on success. On startup you'll also see
`[juliapkg] Installing Julia 1.11.9 ...` the first time (one-time precompile cost).

KEY GOTCHA: `julia` is NOT on the system PATH and that is NORMAL. PySR uses its own
juliapkg-managed Julia, not a system install. Do **not** `apt-get install julia` or
look for `julia` on PATH -- just check `from pysr import PySRRegressor` works.

## Manual fallback (if the hook didn't run / you're mid-session)
This is exactly what we do by hand when needed:
```bash
python3 -m pip install --upgrade numpy scipy scikit-learn pysr
python3 -c "from pysr import PySRRegressor"   # bootstraps + precompiles Julia (~1-3 min first time)
```
Then PySR is usable for the rest of the session.

## Why the OTHER session is probably failing (in likelihood order)
1. **Hook files not on the starting branch.** Most common. The SessionStart hook
   runs only if `.claude/settings.json` + the script exist on the branch the session
   is cut from. A stale/old branch (or a different repo, e.g. Markets) won't have
   them. FIX: copy the 3 files onto that branch, commit, push. (CLAUDE.md S9 already
   flagged: the Markets repos need this same hook treatment.)
2. **Script not executable.** `chmod +x .claude/hooks/session-start.sh` and commit
   the mode bit.
3. **Looking for system `julia`.** It isn't there by design (see gotcha above);
   PySR's juliapkg manages it. Test via `from pysr import PySRRegressor`, not `julia`.
4. **Julia download blocked by network policy.** The hook is best-effort and will
   warn `Julia bootstrap failed`; symbolic work is then unavailable but numpy/scipy
   scripts still run. FIX: use a session whose network policy allows the Julia S3
   download (julialang-s3.julialang.org) + the General registry; or pre-bake Julia
   into the container image.
5. **`set -u` crash on unset harness vars.** Keep the `CLAUDE_PROJECT_DIR`/
   `CLAUDE_ENV_FILE` guard in the hook (fixed S10).

## Notes
- The juliapkg-managed Julia + precompiled SymbolicRegression.jl are cached into the
  container image, so subsequent sessions on the same image are fast.
- Pydroid-3 on a phone CANNOT run PySR (no Julia on Android) -- numpy substrate
  scripts only. PySR needs the cloud container.
- Heavy extras some probes installed in-session (NOT in requirements.txt, to avoid
  breaking the hook): `pycbc` (use `--ignore-installed cryptography`, then force-
  reinstall `--no-deps scipy`), `pint-pulsar`, `hatanaka`, `uproot`/`awkward`. These
  are probe-specific, not part of the Julia/PySR path.
