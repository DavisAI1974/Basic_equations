#!/bin/bash
# SessionStart hook for the Information Layer repo.
# Installs the scientific Python stack and PySR (Julia-backed symbolic
# regression) so the substrate-extraction scripts AND the PySR
# expression-family work both run in a fresh Claude-Code-on-the-web
# container. Idempotent; safe to re-run.
set -euo pipefail

# Web/remote only -- on a local machine, skip (assume env already set up).
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# Core numeric stack (fast, reliable). Required.
python3 -m pip install --quiet --upgrade --root-user-action=ignore \
  numpy scipy scikit-learn

# PySR + its Julia backend. Best-effort: a blocked Julia download must
# NOT break session start -- the pure-numpy substrate scripts still work
# without it. The first import bootstraps Julia + precompiles
# SymbolicRegression.jl; that cost is paid once and cached into the
# container image for subsequent sessions.
if python3 -m pip install --quiet --root-user-action=ignore pysr; then
  python3 -c "import pysr; from pysr import PySRRegressor" \
    && echo "[session-start] PySR + Julia backend ready" \
    || echo "[session-start] WARN: PySR installed but Julia bootstrap failed; PySR work will be unavailable this session (numpy substrate scripts still work)"
else
  echo "[session-start] WARN: PySR pip install failed; continuing without it"
fi

# Let scripts import the repo's modules (per_domain_kbk, etc.) directly.
echo 'export PYTHONPATH="${PYTHONPATH:-}:'"$CLAUDE_PROJECT_DIR"'"' >> "$CLAUDE_ENV_FILE"

echo "[session-start] numeric stack ready"
