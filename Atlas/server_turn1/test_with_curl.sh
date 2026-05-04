#!/bin/bash
# =============================================================================
# test_with_curl.sh
# =============================================================================
# Smoke-test driver for popstats_server. Boots the server with a synthetic
# config and exercises every endpoint with curl. No real engines required —
# the validation, cache, and routing layer can all be tested with stub
# binaries. Real engine integration tests live in the cluster smoke-test
# (turn 7 of chat A).
#
# Usage:
#   bash test_with_curl.sh                 # uses /tmp/popstats_test scratch
#   bash test_with_curl.sh --keep          # keep the scratch dir afterwards
# =============================================================================
set -euo pipefail

TEST_DIR="${TEST_DIR:-/tmp/popstats_smoketest}"
PORT="${PORT:-18765}"
KEEP=0
[[ "${1:-}" == "--keep" ]] && KEEP=1

cleanup() {
  set +e
  [[ -n "${SERVER_PID:-}" ]] && kill "${SERVER_PID}" 2>/dev/null
  wait 2>/dev/null
  if [[ "${KEEP}" -eq 0 ]]; then
    rm -rf "${TEST_DIR}"
  fi
}
trap cleanup EXIT

# ── Build synthetic test fixture ─────────────────────────────────────────────
rm -rf "${TEST_DIR}"
mkdir -p "${TEST_DIR}"/{beagle,bams,ancestry_cache,cache}

# Stub binaries that just succeed (validation tests don't need real engines)
for b in region_popstats hobs_windower angsd instant_q; do
  cat > "${TEST_DIR}/$b" <<EOF
#!/bin/bash
echo "stub $b: \$@" >&2
exit 0
EOF
  chmod +x "${TEST_DIR}/$b"
done

# Sample list with 12 individuals (3-digit padding, matches Quentin's CGA_001 convention)
seq 1 12 | awk '{printf "CGA_%03d\n", $1}' > "${TEST_DIR}/samples.ind"

# Reference .fai stub
echo -e "C_gar_LG28\t20000000\t10\t80\t81" > "${TEST_DIR}/ref.fa.fai"

cat > "${TEST_DIR}/config.yaml" <<YAML
bind: { host: "127.0.0.1", port: ${PORT} }
base: "${TEST_DIR}"
ref: "\${base}/ref.fa"
ref_fai: "\${base}/ref.fa.fai"
engines:
  region_popstats: "\${base}/region_popstats"
  hobs_windower:   "\${base}/hobs_windower"
  angsd_patched:   "\${base}/angsd"
  instant_q:       "\${base}/instant_q"
beagle_dir: "\${base}/beagle"
sample_list: "\${base}/samples.ind"
bam_dir: "\${base}/bams"
bam_suffix: "bam"
local_q_dir: "\${base}/ancestry_cache"
cache_dir: "\${base}/cache"
cache_max_bytes: 10000000
popstats_defaults:
  win_bp: 50000
  step_bp: 10000
  type: 2
  downsample: 1
  ncores: 1
hobs_defaults:
  scales:
    - "10kb:10000:2000"
angsd_defaults:
  threads: 1
  major_minor: 1
  min_maf: 0.05
  min_mapq: 20
  min_q: 20
  snp_pval: "1e-6"
  max_het_freq: 1.0
min_group_n: 3
cors_origins: ["*"]
YAML

# ── Boot server ──────────────────────────────────────────────────────────────
echo "starting server on port ${PORT}..."
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POPSTATS_CONFIG="${TEST_DIR}/config.yaml" \
  python3 -m uvicorn --app-dir "${HERE}" \
  popstats_server:app --host 127.0.0.1 --port "${PORT}" --log-level error &
SERVER_PID=$!
sleep 2

# Helper
hit() {
  local label="$1"; shift
  echo
  echo "=== ${label} ==="
  curl -sS "$@" | python3 -m json.tool || true
}

URL="http://127.0.0.1:${PORT}"

hit "GET /api/health" "${URL}/api/health"

hit "GET /api/cache/keys" "${URL}/api/cache/keys"

hit "POST /api/popstats/groupwise — empty groups → 400" \
  -X POST "${URL}/api/popstats/groupwise" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","groups":{}}'

hit "POST /api/popstats/groupwise — group too small → 400" \
  -X POST "${URL}/api/popstats/groupwise" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","groups":{"A":["CGA_001"]}}'

hit "POST /api/popstats/groupwise — unknown sample → 400" \
  -X POST "${URL}/api/popstats/groupwise" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","groups":{"A":["NOPE_001","NOPE_002","NOPE_003"]}}'

hit "POST /api/popstats/groupwise — chrom path traversal → 422" \
  -X POST "${URL}/api/popstats/groupwise" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"../../etc/passwd","groups":{"A":["CGA_001","CGA_002","CGA_003"]}}'

hit "POST /api/popstats/hobs_groupwise — bad scale label → 400" \
  -X POST "${URL}/api/popstats/hobs_groupwise" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","scales":["bogus"],"groups":{"HOM1":["CGA_001","CGA_002","CGA_003"]}}'

hit "POST /api/ancestry/groupwise_q — missing local_Q → 404" \
  -X POST "${URL}/api/ancestry/groupwise_q" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","K":8,"scale":"dense","groups":{"A":["CGA_001","CGA_002","CGA_003"]}}'

hit "POST /api/shelf_ld_test — explicit not-implemented → 501" \
  -X POST "${URL}/api/shelf_ld_test" \
  -H 'Content-Type: application/json' \
  -d '{"chrom":"C_gar_LG28","shelf":{"start_bp":15000000,"end_bp":18000000}}'

hit "GET /api/jobs/nope → 404" "${URL}/api/jobs/nope"

echo
echo "=== ALL SMOKE TESTS DONE ==="
