#!/bin/bash
# =============================================================================
# scripts/registry_query.sh — quick registry interrogation from the CLI
# =============================================================================
# Usage:
#   bash scripts/registry_query.sh list_candidates
#   bash scripts/registry_query.sh list_candidates --chrom C_gar_LG28
#   bash scripts/registry_query.sh describe <cid>
#   bash scripts/registry_query.sh karyotypes <cid>
#   bash scripts/registry_query.sh fst <cid>
#   bash scripts/registry_query.sh summary
# =============================================================================
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${here}/00_config.sh"

: "${REGISTRIES:=${BASE}/inversion-popgen-toolkit/registries}"
: "${REGISTRY_API_R:=${REGISTRIES}/api/R}"

CMD="${1:-summary}"
shift || true

args=()
while [[ $# -gt 0 ]]; do
  args+=( "$1" ); shift
done

${RSCRIPT_BIN} --vanilla "${here}/R/registry_query.R" \
  --registries_root "${REGISTRIES}" \
  --registry_api_r  "${REGISTRY_API_R}" \
  --cmd "${CMD}" \
  "${args[@]}"
