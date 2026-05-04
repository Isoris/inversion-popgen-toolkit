#!/usr/bin/env bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 0-01:00:00
#SBATCH -J bpv_s01
#SBATCH --account=lt200308
# =============================================================================
# STEP_A01_extract_inv_candidates.sh — v3: robust dual-caller extraction
# =============================================================================
# Extracts INV candidates from DELLY2 + Manta, applies caller-specific
# quality filters, deduplicates, and matches to Snake candidates.
#
# Self-contained: sources own config, creates dirs, validates inputs.
#
# Output columns (19):
#   inv_id, chrom, bp1_pos, bp2_pos, svlen_bp, qual, pe, sr, precise, filter,
#   caller, ct, junction_type, cipos, ciend,
#   cross_caller_match, concordance_class,
#   snake_match, snake_candidate_id
# =============================================================================
set -euo pipefail

# ── Source config (works both interactive and SLURM) ──
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
BPV_CONFIG="${SCRIPT_DIR}/00_phase3_config.sh"

if [[ ! -f "${BPV_CONFIG}" ]]; then
  echo "[FATAL] Config not found: ${BPV_CONFIG}" >&2
  exit 1
fi

# Export ALL variables so Python subprocesses see them
set -a
source "${BPV_CONFIG}"
set +a

# ── Initialize output directories ──
bpv_init_dirs

# ── Validate required inputs ──
bpv_log "================================================================"
bpv_log "  STEP 01: Extract INV candidates (v3 — robust)"
bpv_log "================================================================"
bpv_log ""
bpv_log "  Config:       ${BPV_CONFIG}"
bpv_log "  Output root:  ${BPV_ROOT}"
bpv_log "  Evidence dir: ${BPV_EVIDENCE}"
bpv_log ""

ERRORS=0

check_input() {
  local f="$1" label="$2" required="${3:-yes}"
  if [[ -f "$f" ]]; then
    local n
    if [[ "$f" == *.vcf.gz ]]; then
      n=$(bcftools view -H "$f" 2>/dev/null | wc -l)
      bpv_log "  [OK] ${label}: ${f} (${n} records)"
    else
      n=$(wc -l < "$f")
      bpv_log "  [OK] ${label}: ${f} (${n} lines)"
    fi
  elif [[ "$required" == "yes" ]]; then
    bpv_log "  [MISSING] ${label}: ${f}"
    ERRORS=$((ERRORS + 1))
  else
    bpv_log "  [SKIP] ${label}: ${f} (optional, not found)"
  fi
}

bpv_log "── Input validation ──"
check_input "${DELLY_INV_VCF}" "DELLY INV VCF" "no"
check_input "${MANTA_INV_VCF}" "Manta INV VCF" "no"
check_input "${SAMPLES_ALL}" "Sample list" "yes"

# At least one caller must exist
if [[ ! -f "${DELLY_INV_VCF}" ]] && [[ ! -f "${MANTA_INV_VCF}" ]]; then
  bpv_log "  [FATAL] Neither DELLY nor Manta INV VCF found"
  ERRORS=$((ERRORS + 1))
fi

# Check optional inputs
check_input "${SNAKE_SCORES:-missing}" "Snake scores" "no"
check_input "${DELLY_BND_VCF:-missing}" "DELLY BND VCF" "no"

if [[ ${ERRORS} -gt 0 ]]; then
  bpv_die "Input validation failed with ${ERRORS} errors"
fi

# Check BAM directory (needed by step 02, warn now)
if [[ -d "${DELLY_MARKDUP_DIR}" ]]; then
  N_BAMS=$(ls "${DELLY_MARKDUP_DIR}"/*.bam 2>/dev/null | wc -l)
  bpv_log "  [INFO] Markdup BAMs: ${N_BAMS} in ${DELLY_MARKDUP_DIR}"
else
  bpv_log "  [WARN] Markdup BAM dir not found: ${DELLY_MARKDUP_DIR} (needed for step 02)"
fi

bpv_log ""

# =============================================================================
# 1A. Parse DELLY2 INV catalog
# =============================================================================
bpv_log "── 1A: Parsing DELLY2 INV catalog ──"

DELLY_PARSED="${BPV_EVIDENCE}/delly_inv_candidates.tsv"
export DELLY_PARSED

N_DELLY=0
if [[ -f "${DELLY_INV_VCF}" ]]; then
  python3 << 'PYEOF'
import gzip, os, sys

vcf = os.environ["DELLY_INV_VCF"]
out = os.environ["DELLY_PARSED"]
min_qual = int(os.environ.get("MIN_DELLY_QUAL", "300"))
min_pe = int(os.environ.get("MIN_DELLY_PE", "3"))
max_size = int(os.environ.get("MAX_INV_SIZE_MB", "20")) * 1_000_000
min_size = int(os.environ.get("MIN_INV_SIZE_BP", "5000"))

print(f"  DELLY VCF: {vcf}")
print(f"  Filters: QUAL>={min_qual} PE>={min_pe} size=[{min_size/1e3:.0f}kb, {max_size/1e6:.0f}Mb]")

header_cols = "inv_id\tchrom\tbp1_pos\tbp2_pos\tsvlen_bp\tqual\tpe\tsr\tprecise\tfilter\tcaller\tct\tjunction_type\tcipos\tciend"

with gzip.open(vcf, 'rt') as v, open(out, 'w') as o:
    o.write(header_cols + "\n")
    n_total = 0; n_pass = 0; n_skip_qual = 0; n_skip_size = 0
    chrom_counts = {}

    for line in v:
        if line.startswith('#'): continue
        n_total += 1
        p = line.strip().split('\t')
        chrom, pos, inv_id, filt = p[0], int(p[1]), p[2], p[6]
        info = dict(kv.split('=', 1) for kv in p[7].split(';') if '=' in kv)
        flags = set(kv for kv in p[7].split(';') if '=' not in kv)

        qual = float(p[5]) if p[5] != '.' else 0
        end = int(info.get('END', pos))
        svlen = abs(end - pos)
        pe = int(info.get('PE', 0))
        sr = int(info.get('SR', 0))
        precise = 1 if 'PRECISE' in flags else 0
        ct = info.get('CT', '.')
        cipos = info.get('CIPOS', '.')
        ciend = info.get('CIEND', '.')

        # Junction type
        jtype = "."
        if ct == "3to3": jtype = "3to3_left"
        elif ct == "5to5": jtype = "5to5_right"
        elif ct == "3to5": jtype = "normal_3to5"
        elif ct == "5to3": jtype = "tandem_5to3"

        if qual < min_qual or pe < min_pe:
            n_skip_qual += 1; continue
        if svlen > max_size or svlen < min_size:
            n_skip_size += 1; continue

        n_pass += 1
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        o.write(f"{inv_id}\t{chrom}\t{pos}\t{end}\t{svlen}\t{qual:.0f}\t{pe}\t{sr}\t{precise}\t{filt}\tdelly\t{ct}\t{jtype}\t{cipos}\t{ciend}\n")

    print(f"  DELLY: {n_total} total → {n_pass} pass (skip: {n_skip_qual} qual, {n_skip_size} size)")
    print(f"  Chromosomes: {len(chrom_counts)} | Top 5: {sorted(chrom_counts.items(), key=lambda x:-x[1])[:5]}")
    n_precise = sum(1 for line in open(out) if '\t1\t' in line.split('\t')[8:9])
PYEOF

  N_DELLY=$(tail -n +2 "${DELLY_PARSED}" | wc -l)
  bpv_log "  DELLY INV candidates: ${N_DELLY}"
else
  bpv_log "  DELLY INV VCF not found — skipping"
  echo -e "inv_id\tchrom\tbp1_pos\tbp2_pos\tsvlen_bp\tqual\tpe\tsr\tprecise\tfilter\tcaller\tct\tjunction_type\tcipos\tciend" > "${DELLY_PARSED}"
fi

# =============================================================================
# 1B. Parse Manta INV catalog
# =============================================================================
bpv_log ""
bpv_log "── 1B: Parsing Manta INV catalog ──"

MANTA_PARSED="${BPV_EVIDENCE}/manta_inv_candidates.tsv"
export MANTA_PARSED

N_MANTA=0
if [[ -f "${MANTA_INV_VCF}" ]]; then
  python3 << 'PYEOF'
import gzip, os

vcf = os.environ["MANTA_INV_VCF"]
out = os.environ["MANTA_PARSED"]
min_qual = int(os.environ.get("MIN_MANTA_QUAL", "50"))
min_total_alt = int(os.environ.get("MIN_MANTA_TOTAL_ALT", "6"))
max_size = int(os.environ.get("MAX_INV_SIZE_MB", "20")) * 1_000_000
min_size = int(os.environ.get("MIN_INV_SIZE_BP", "5000"))

print(f"  Manta VCF: {vcf}")
print(f"  Filters: QUAL>={min_qual} total_alt>={min_total_alt} size=[{min_size/1e3:.0f}kb, {max_size/1e6:.0f}Mb]")

header_cols = "inv_id\tchrom\tbp1_pos\tbp2_pos\tsvlen_bp\tqual\tpe\tsr\tprecise\tfilter\tcaller\tct\tjunction_type\tcipos\tciend"

with gzip.open(vcf, 'rt') as v, open(out, 'w') as o:
    o.write(header_cols + "\n")
    sample_names = []
    n_total = 0; n_pass = 0; n_skip_qual = 0; n_skip_size = 0
    chrom_counts = {}
    n_inv3 = 0; n_inv5 = 0; n_both = 0; n_neither = 0

    for line in v:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            sample_names = line.strip().split('\t')[9:]
            continue

        n_total += 1
        p = line.strip().split('\t')
        chrom, pos, inv_id, filt = p[0], int(p[1]), p[2], p[6]
        info = dict(kv.split('=', 1) for kv in p[7].split(';') if '=' in kv)
        flags = set(kv for kv in p[7].split(';') if '=' not in kv)

        qual = float(p[5]) if p[5] != '.' else 0
        end = int(info.get('END', pos))
        svlen = abs(end - pos)
        precise = 1 if 'IMPRECISE' not in flags else 0
        cipos = info.get('CIPOS', '.')
        ciend = info.get('CIEND', '.')

        # Junction type from INV3/INV5 flags
        has_inv3 = 'INV3' in flags
        has_inv5 = 'INV5' in flags
        if has_inv3 and has_inv5:
            jtype = "both"; n_both += 1
        elif has_inv3:
            jtype = "3to3_left"; n_inv3 += 1
        elif has_inv5:
            jtype = "5to5_right"; n_inv5 += 1
        else:
            jtype = "."; n_neither += 1

        # Aggregate PR/SR across carriers
        fmt_keys = p[8].split(':') if len(p) > 8 else []
        pr_idx = fmt_keys.index('PR') if 'PR' in fmt_keys else -1
        sr_idx = fmt_keys.index('SR') if 'SR' in fmt_keys else -1
        gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else 0

        sum_pr = 0; sum_sr = 0; n_carriers = 0
        for sd in p[9:]:
            sv = sd.split(':')
            gt = sv[gt_idx] if gt_idx < len(sv) else './.'
            if gt in ('./.', '0/0', './0', '0/.', '.', '0|0'): continue
            n_carriers += 1
            if pr_idx >= 0 and pr_idx < len(sv):
                try: sum_pr += int(sv[pr_idx].split(',')[1])
                except: pass
            if sr_idx >= 0 and sr_idx < len(sv):
                try: sum_sr += int(sv[sr_idx].split(',')[1])
                except: pass

        total_alt = sum_pr + sum_sr

        if qual < min_qual or total_alt < min_total_alt:
            n_skip_qual += 1; continue
        if svlen > max_size or svlen < min_size:
            n_skip_size += 1; continue

        n_pass += 1
        chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1
        ct = "3to3" if has_inv3 else ("5to5" if has_inv5 else ".")
        o.write(f"{inv_id}\t{chrom}\t{pos}\t{end}\t{svlen}\t{qual:.0f}\t{sum_pr}\t{sum_sr}\t{precise}\t{filt}\tmanta\t{ct}\t{jtype}\t{cipos}\t{ciend}\n")

    print(f"  Manta: {n_total} total → {n_pass} pass (skip: {n_skip_qual} qual/evidence, {n_skip_size} size)")
    print(f"  Junction types: INV3={n_inv3} INV5={n_inv5} both={n_both} neither={n_neither}")
    print(f"  Chromosomes: {len(chrom_counts)} | Top 5: {sorted(chrom_counts.items(), key=lambda x:-x[1])[:5]}")
PYEOF

  N_MANTA=$(tail -n +2 "${MANTA_PARSED}" | wc -l)
  bpv_log "  Manta INV candidates: ${N_MANTA}"
else
  bpv_log "  Manta INV VCF not found — skipping"
  echo -e "inv_id\tchrom\tbp1_pos\tbp2_pos\tsvlen_bp\tqual\tpe\tsr\tprecise\tfilter\tcaller\tct\tjunction_type\tcipos\tciend" > "${MANTA_PARSED}"
fi

# Check: at least some candidates
TOTAL=$((N_DELLY + N_MANTA))
if [[ ${TOTAL} -eq 0 ]]; then
  bpv_die "No candidates from either caller — nothing to validate"
fi
bpv_log ""
bpv_log "  Total candidates before merge: ${TOTAL} (${N_DELLY} DELLY + ${N_MANTA} Manta)"

# =============================================================================
# 1C. Merge + deduplicate
# =============================================================================
bpv_log ""
bpv_log "── 1C: Merging and deduplicating ──"

MERGED_RAW="${BPV_EVIDENCE}/merged_inv_candidates_raw.tsv"
MERGED_DEDUP="${BPV_EVIDENCE}/merged_inv_candidates_dedup.tsv"
export MERGED_RAW MERGED_DEDUP

head -1 "${DELLY_PARSED}" > "${MERGED_RAW}"
tail -n +2 "${DELLY_PARSED}" >> "${MERGED_RAW}"
tail -n +2 "${MANTA_PARSED}" >> "${MERGED_RAW}"

N_RAW=$(tail -n +2 "${MERGED_RAW}" | wc -l)
bpv_log "  Raw merged: ${N_RAW}"

python3 << 'PYEOF'
import os

raw = os.environ["MERGED_RAW"]
out = os.environ["MERGED_DEDUP"]

rows = []
with open(raw) as f:
    header = f.readline().strip()
    for line in f:
        p = line.strip().split('\t')
        if len(p) < 15:
            print(f"  WARNING: short line ({len(p)} cols): {p[0] if p else '?'}")
            continue
        rows.append({
            'inv_id': p[0], 'chrom': p[1], 'bp1': int(p[2]), 'bp2': int(p[3]),
            'svlen': int(p[4]), 'qual': float(p[5]), 'pe': int(p[6]), 'sr': int(p[7]),
            'precise': p[8], 'filter': p[9], 'caller': p[10],
            'ct': p[11], 'junction_type': p[12],
            'cipos': p[13], 'ciend': p[14],
            'cross_caller_match': '.', 'concordance_class': 'single_caller',
        })

# Find overlapping DELLY-Manta pairs (≥50% reciprocal overlap)
delly_rows = [r for r in rows if r['caller'] == 'delly']
manta_rows = [r for r in rows if r['caller'] == 'manta']

matched_delly = set()
matched_manta = set()

for di, d in enumerate(delly_rows):
    for mi, m in enumerate(manta_rows):
        if d['chrom'] != m['chrom']: continue
        overlap = max(0, min(d['bp2'], m['bp2']) - max(d['bp1'], m['bp1']))
        d_len = d['bp2'] - d['bp1']
        m_len = m['bp2'] - m['bp1']
        if d_len <= 0 or m_len <= 0: continue
        if overlap / d_len >= 0.5 or overlap / m_len >= 0.5:
            d['cross_caller_match'] = m['inv_id']
            d['concordance_class'] = 'both_callers'
            m['cross_caller_match'] = d['inv_id']
            m['concordance_class'] = 'both_callers'
            matched_delly.add(di)
            matched_manta.add(mi)

n_both = len(matched_delly)
n_delly_only = len(delly_rows) - n_both
n_manta_only = len(manta_rows) - len(matched_manta)
print(f"  Cross-caller concordance:")
print(f"    Both callers:  {n_both} pairs")
print(f"    DELLY-only:    {n_delly_only}")
print(f"    Manta-only:    {n_manta_only}")

with open(out, 'w') as o:
    o.write(header + "\tcross_caller_match\tconcordance_class\n")
    for r in rows:
        vals = [r['inv_id'], r['chrom'], str(r['bp1']), str(r['bp2']),
                str(r['svlen']), f"{r['qual']:.0f}", str(r['pe']), str(r['sr']),
                r['precise'], r['filter'], r['caller'], r['ct'], r['junction_type'],
                r['cipos'], r['ciend'],
                r['cross_caller_match'], r['concordance_class']]
        o.write('\t'.join(vals) + '\n')

print(f"  Output: {out} ({len(rows)} rows, 17 columns)")
PYEOF

N_DEDUP=$(tail -n +2 "${MERGED_DEDUP}" | wc -l)
bpv_log "  After dedup tagging: ${N_DEDUP} candidates"

# =============================================================================
# 1D. Match to Snake candidates (optional — skipped if no scores file)
# =============================================================================
bpv_log ""
bpv_log "── 1D: Matching to Snake inversion candidates ──"

MATCHED="${BPV_EVIDENCE}/matched_inv_candidates.tsv"
export MATCHED

SNAKE_FILE="${SNAKE_SCORES:-}"
if [[ -n "${SNAKE_FILE}" ]] && [[ -f "${SNAKE_FILE}" ]]; then
  bpv_log "  Snake scores found: ${SNAKE_FILE}"
else
  bpv_log "  Snake scores not found (${SNAKE_FILE:-not set}) — all candidates kept without snake annotation"
  SNAKE_FILE=""
fi

python3 << 'PYEOF'
import os, gzip

dedup_file = os.environ["MERGED_DEDUP"]
snake_file = os.environ.get("SNAKE_FILE", "")
out_file = os.environ["MATCHED"]

# Load candidates
candidates = []
with open(dedup_file) as f:
    header = f.readline().strip()
    for line in f:
        p = line.strip().split('\t')
        if len(p) < 17:
            print(f"  WARNING: short line in dedup ({len(p)} cols)")
            continue
        candidates.append({
            'inv_id': p[0], 'chrom': p[1], 'bp1': int(p[2]), 'bp2': int(p[3]),
            'svlen': int(p[4]), 'qual': p[5], 'pe': p[6], 'sr': p[7],
            'precise': p[8], 'filter': p[9], 'caller': p[10],
            'ct': p[11], 'junction_type': p[12],
            'cipos': p[13], 'ciend': p[14],
            'cross_caller_match': p[15], 'concordance_class': p[16],
        })

print(f"  Candidates loaded: {len(candidates)}")

# Load Snake candidates (optional)
snake = []
if snake_file and os.path.exists(snake_file):
    opener = gzip.open if snake_file.endswith('.gz') else open
    with opener(snake_file, 'rt') as f:
        sheader = f.readline().strip().split('\t')
        for line in f:
            p = line.strip().split('\t')
            if len(p) >= len(sheader):
                snake.append(dict(zip(sheader, p)))
    print(f"  Snake candidates loaded: {len(snake)}")
else:
    print(f"  No snake scores — skipping matching")

# Match: ≥50% reciprocal overlap
n_matched = 0
out_header = (
    "inv_id\tchrom\tbp1_pos\tbp2_pos\tsvlen_bp\tqual\tpe\tsr\tprecise\tfilter\t"
    "caller\tct\tjunction_type\tcipos\tciend\t"
    "cross_caller_match\tconcordance_class\t"
    "snake_match\tsnake_candidate_id\n"
)

with open(out_file, 'w') as o:
    o.write(out_header)
    for d in candidates:
        best_match = "none"; best_cid = "."
        for s in snake:
            sc = s.get('chrom', s.get('chr', ''))
            if sc != d['chrom']: continue
            ss = int(s.get('start_bp', s.get('start', 0)))
            se = int(s.get('end_bp', s.get('end', 0)))
            overlap = max(0, min(d['bp2'], se) - max(d['bp1'], ss))
            dlen = d['bp2'] - d['bp1']; slen = se - ss
            if dlen > 0 and slen > 0:
                if overlap / dlen >= 0.5 or overlap / slen >= 0.5:
                    best_match = "overlap"
                    best_cid = s.get('candidate_id', s.get('cid', '.'))
                    n_matched += 1
                    break

        o.write(f"{d['inv_id']}\t{d['chrom']}\t{d['bp1']}\t{d['bp2']}\t{d['svlen']}\t"
                f"{d['qual']}\t{d['pe']}\t{d['sr']}\t{d['precise']}\t{d['filter']}\t"
                f"{d['caller']}\t{d['ct']}\t{d['junction_type']}\t{d['cipos']}\t{d['ciend']}\t"
                f"{d['cross_caller_match']}\t{d['concordance_class']}\t"
                f"{best_match}\t{best_cid}\n")

print(f"  Snake matches: {n_matched} / {len(candidates)}")
print(f"  Output: {out_file} ({len(candidates)} rows, 19 columns)")
PYEOF

N_MATCHED=$(grep -c 'overlap' "${MATCHED}" 2>/dev/null || echo 0)

# =============================================================================
# 1E. Summary
# =============================================================================
bpv_log ""
bpv_log "================================================================"
bpv_log "  STEP 01 COMPLETE"
bpv_log "================================================================"
bpv_log "  DELLY candidates:       ${N_DELLY}"
bpv_log "  Manta candidates:       ${N_MANTA}"
bpv_log "  Total merged:           ${N_RAW}"
bpv_log "  Cross-caller pairs:     $(grep -c 'both_callers' "${MERGED_DEDUP}" 2>/dev/null || echo 0)"
bpv_log "  Snake-matched:          ${N_MATCHED}"
bpv_log ""
bpv_log "  Output files:"
bpv_log "    DELLY parsed:     ${DELLY_PARSED}"
bpv_log "    Manta parsed:     ${MANTA_PARSED}"
bpv_log "    Merged raw:       ${MERGED_RAW}"
bpv_log "    Merged dedup:     ${MERGED_DEDUP}"
bpv_log "    Final matched:    ${MATCHED}"
bpv_log ""

# Verify outputs exist and have data
for f in "${DELLY_PARSED}" "${MANTA_PARSED}" "${MERGED_RAW}" "${MERGED_DEDUP}" "${MATCHED}"; do
  if [[ -s "$f" ]]; then
    bpv_log "    [OK] $(basename "$f"): $(tail -n +2 "$f" | wc -l) records"
  else
    bpv_log "    [WARN] $(basename "$f"): empty or missing"
  fi
done

bpv_log ""
bpv_log "  Next: bash STEP_A02_extract_breakpoint_evidence.py --candidates ${MATCHED}"
