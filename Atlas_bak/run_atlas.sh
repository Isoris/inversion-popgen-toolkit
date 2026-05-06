#!/usr/bin/env bash
# =============================================================================
# run_atlas.sh — single-command launcher for the merged atlas server (turn 145)
# =============================================================================
# Run from the Atlas/ directory:
#
#     ./run_atlas.sh
#
# Brings up the merged server (server_turn1/popstats_server.py) with:
#   - /file + /compute subsystem      project_root = $PWD
#   - popstats subsystem               enabled iff
#                                      server_turn1/popstats_server.config.yaml exists
#
# Default bind: 127.0.0.1:8765 (localhost only). Override with --host / --port.
#
# Flags:
#   --no-popstats       skip popstats even if config file exists
#   --config <file>     use a non-default popstats config path
#   --port <n>          bind port (default 8765)
#   --host <ip>         bind host (default 127.0.0.1)
#   --check             dry run: print what would happen and exit 0
#   --help              this message
#
# See docs/SERVER_AUDIT_2026-05-05.md for the architecture.
# =============================================================================

set -euo pipefail

# Resolve Atlas/ directory (where this script lives) regardless of CWD.
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SERVER_DIR="${SCRIPT_DIR}/server_turn1"
SERVER_PY="${SERVER_DIR}/popstats_server.py"
DEFAULT_CFG="${SERVER_DIR}/popstats_server.config.yaml"

# ---- arg parsing ------------------------------------------------------------
NO_POPSTATS=0
CFG_OVERRIDE=""
HOST="127.0.0.1"
PORT="8765"
CHECK_ONLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-popstats) NO_POPSTATS=1; shift ;;
    --config)      CFG_OVERRIDE="$2"; shift 2 ;;
    --host)        HOST="$2"; shift 2 ;;
    --port)        PORT="$2"; shift 2 ;;
    --check)       CHECK_ONLY=1; shift ;;
    --help|-h)
      sed -n '2,30p' "$0" | sed 's/^# \?//'
      exit 0 ;;
    *)
      echo "unknown flag: $1" >&2
      echo "see --help" >&2
      exit 2 ;;
  esac
done

# ---- resolve config ---------------------------------------------------------
CFG_TO_USE=""
if [[ "${NO_POPSTATS}" -eq 0 ]]; then
  if [[ -n "${CFG_OVERRIDE}" ]]; then
    if [[ ! -f "${CFG_OVERRIDE}" ]]; then
      echo "error: --config file not found: ${CFG_OVERRIDE}" >&2
      exit 2
    fi
    CFG_TO_USE="${CFG_OVERRIDE}"
  elif [[ -f "${DEFAULT_CFG}" ]]; then
    CFG_TO_USE="${DEFAULT_CFG}"
  fi
fi

# ---- preflight --------------------------------------------------------------
if [[ ! -f "${SERVER_PY}" ]]; then
  echo "error: merged server not found at ${SERVER_PY}" >&2
  echo "  expected layout: <atlas_root>/server_turn1/popstats_server.py" >&2
  exit 2
fi

# Pick python3
if command -v python3 >/dev/null 2>&1; then
  PY="python3"
elif command -v python >/dev/null 2>&1; then
  PY="python"
else
  echo "error: no python3 found on PATH" >&2
  exit 2
fi

# Verify FastAPI is importable. Catches the most common new-machine error
# before the server flails.
if ! "${PY}" -c "import fastapi, uvicorn" >/dev/null 2>&1; then
  echo "error: required Python packages missing." >&2
  echo "  install with: ${PY} -m pip install fastapi uvicorn pyyaml numpy pandas pydantic" >&2
  echo "  (see ${SERVER_DIR}/requirements.txt)" >&2
  exit 3
fi

# ---- summary ----------------------------------------------------------------
echo "========================================================================="
echo " atlas server launcher"
echo "========================================================================="
echo "  Atlas root:    ${SCRIPT_DIR}"
echo "  Server:        ${SERVER_PY}"
echo "  Bind:          http://${HOST}:${PORT}"
echo "  Project root:  ${SCRIPT_DIR}        (/file + /compute subsystem)"
if [[ -n "${CFG_TO_USE}" ]]; then
  echo "  Popstats cfg:  ${CFG_TO_USE}        (popstats subsystem)"
else
  echo "  Popstats:      disabled (no config provided)"
  if [[ "${NO_POPSTATS}" -eq 0 && -z "${CFG_OVERRIDE}" ]]; then
    echo "                 to enable: cp ${SERVER_DIR}/popstats_server.config.example.yaml \\"
    echo "                              ${DEFAULT_CFG}"
    echo "                 then edit paths and re-run this launcher."
  fi
fi
echo "========================================================================="

if [[ "${CHECK_ONLY}" -eq 1 ]]; then
  echo "(--check) dry run only; not starting server"
  exit 0
fi

# ---- launch -----------------------------------------------------------------
# We exec from server_turn1/ so popstats_server's relative imports resolve.
cd "${SERVER_DIR}"

CMD=( "${PY}" popstats_server.py
      --project-root "${SCRIPT_DIR}"
      --host "${HOST}"
      --port "${PORT}" )

if [[ -n "${CFG_TO_USE}" ]]; then
  CMD+=( --config "${CFG_TO_USE}" )
fi

# exec replaces this shell with the python process so Ctrl+C goes straight
# to uvicorn — clean shutdown, no stale shell handler in the way.
exec "${CMD[@]}"
