#!/usr/bin/env bash
# =============================================================================
# annotate_population_confidence.sh — Add confidence scores to ALL PASS SVs
# =============================================================================
#
# This is NOT a filter. Every PASS SV from DELLY/Manta keeps its place.
# We ADD columns that score how confident we are about each variant:
#
#   1. pop_signal_score  — binomial -log10(P): how unlikely is this many
#      carriers showing evidence at the same spot by chance?
#
#   2. confidence_level  — LOW / MEDIUM / HIGH / VERY_HIGH based on score
#
#   3. breakpoint_precision — estimated bp-level precision from the actual
#      spread of evidence positions across carriers (narrower = more precise,
#      independent of what CIPOS says)
#
# USE CASE: When the inversion pipeline finds a boundary at some position,
# it can look up the confidence_level at that position. If there's a
# HIGH-confidence SV breakpoint within 30 bp, pick that one. If only
# LOW-confidence calls nearby, flag it as uncertain.
#
# The binomial null: at any random position, ~5% of samples show ≥1 alt
# read by chance (mapping noise, repeat artifacts). A real SV breakpoint
# will have far more carriers piling up evidence at the same narrow window.
#
# Both 226 and 81 cohorts are annotated.
# All SV types: DEL, DUP, INV, INS_small, INS_large, BND.
#
# Usage:
#   cd /scratch/.../MODULE_4H_ALL_Manta
#   bash annotate_population_confidence.sh
# =============================================================================
set -euo pipefail
source ~/.bashrc
mamba activate assembly

SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CONFIG="${SCRIPT_DIR}/00_manta_config.sh"
[[ -f "${CONFIG}" ]] || { echo "Missing config: ${CONFIG}" >&2; exit 1; }
set -a
source "${CONFIG}"
set +a

mv_log "=== Population Confidence Annotation ==="
mv_log "  Not a filter — annotates every PASS SV with confidence scores"

ANNOT_DIR="${OUTDIR}/10_population_confidence"
mkdir -p "${ANNOT_DIR}"

SV_TYPES=(DEL DUP INV INS_small INS_large BND)

# =============================================================================
# Core: annotate one VCF
# =============================================================================
annotate_vcf() {
  local vcf="$1" out_tsv="$2" svtype="$3"

  python3 - "${vcf}" "${out_tsv}" "${svtype}" << 'PYEOF'
import gzip, sys, math
from collections import Counter

vcf_path, out_path, svtype_label = sys.argv[1], sys.argv[2], sys.argv[3]

header = [
    # ── Core variant info ────────────────────────────────────────────
    "chrom", "pos", "end", "id", "svtype", "svlen_raw", "qual", "filter",
    "is_imprecise", "cipos", "ciend", "homlen", "homseq",
    "abs_svlen",

    # ── Sample counts ────────────────────────────────────────────────
    "n_samples",
    "carriers", "carrier_frac", "het_count", "hom_count", "af",

    # ── Per-sample evidence totals ───────────────────────────────────
    "sum_PR_alt", "sum_SR_alt", "sum_total_alt",
    "max_PR_alt", "max_SR_alt",
    "min_carrier_alt", "median_carrier_alt",

    # ── POPULATION CONFIDENCE (the annotation layer) ─────────────────
    "n_carriers_with_evidence",     # samples with ≥1 alt read
    "n_carriers_with_strong_ev",    # samples with ≥3 alt reads
    "carrier_evidence_frac",        # n_ev / n_samples
    "pop_signal_score",             # -log10(binomial P)
    "confidence_level",             # LOW / MEDIUM / HIGH / VERY_HIGH

    # ── Breakpoint precision estimate ────────────────────────────────
    "cipos_width",                  # from CIPOS field (caller's estimate)
    "n_split_read_carriers",        # carriers with SR_alt ≥ 1 (precise bp)
    "precision_class",              # base_pair / narrow / wide / unknown

    # ── Size classification ──────────────────────────────────────────
    "size_class",

    # ── SV-type-specific fields ──────────────────────────────────────
    "inv3", "inv5", "junction_type", "event_id",
    "left_svinsseq_len", "right_svinsseq_len",
    "mate_id", "bnd_alt", "bnd_orientation",
]

def binom_log10p(k, n, p=0.05):
    """Population signal: -log10(P) for k or more hits in n trials."""
    if n == 0 or k == 0: return 0.0
    mu = n * p
    sigma = math.sqrt(n * p * (1 - p))
    if sigma == 0: return 0.0
    z = (k - 0.5 - mu) / sigma
    if z <= 0: return 0.0
    return min(-math.log10(max(0.5 * math.erfc(z / math.sqrt(2)), 1e-300)), 300)

def confidence_from_score(score):
    """Map pop_signal_score to a human-readable confidence level."""
    if score >= 10:  return "VERY_HIGH"   # P < 1e-10
    if score >= 5:   return "HIGH"        # P < 1e-5
    if score >= 2:   return "MEDIUM"      # P < 0.01
    if score >= 1.3: return "LOW"         # P < 0.05
    return "MINIMAL"                      # P ≥ 0.05

def cipos_to_width(cipos_str):
    """Parse CIPOS=a,b and return width = b - a."""
    try:
        parts = cipos_str.split(',')
        return int(parts[1]) - int(parts[0])
    except:
        return -1

def precision_class(cipos_width, n_sr_carriers):
    """Estimate breakpoint precision from CIPOS width and split-read support."""
    if n_sr_carriers >= 3 and cipos_width >= 0 and cipos_width <= 5:
        return "base_pair"     # ≤5bp CI + multiple split-read carriers
    elif n_sr_carriers >= 1 and cipos_width >= 0 and cipos_width <= 50:
        return "narrow"        # ≤50bp CI + some split-read
    elif cipos_width >= 0 and cipos_width <= 200:
        return "moderate"      # ≤200bp CI
    elif cipos_width > 200:
        return "wide"          # >200bp CI
    else:
        return "unknown"

rows = []
with gzip.open(vcf_path, 'rt') as f:
    snames = []
    for line in f:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            snames = line.strip().split('\t')[9:]; continue
        fld = line.strip().split('\t')
        if len(fld) < 10: continue

        ch, ps, vid, alt, qs, flt = fld[0], int(fld[1]), fld[2], fld[4], fld[5], fld[6]

        # Parse INFO
        info = {}
        for it in fld[7].split(';'):
            if '=' in it:
                k, v = it.split('=', 1); info[k] = v
            else:
                info[it] = True

        svt = info.get('SVTYPE', svtype_label)
        ev = info.get('END', str(ps))
        slr = info.get('SVLEN', '.')
        imp = 'IMPRECISE' in info
        cip = info.get('CIPOS', '.')
        cie = info.get('CIEND', '.')
        hl = info.get('HOMLEN', '.')
        hs = info.get('HOMSEQ', '.')

        # Type-specific
        i3 = 'INV3' in info; i5 = 'INV5' in info
        jt = "both" if i3 and i5 else "3to3_left" if i3 else "5to5_right" if i5 else "."
        eid = info.get('EVENT', '.')
        li = info.get('LEFT_SVINSSEQ', ''); ri = info.get('RIGHT_SVINSSEQ', '')
        mid = info.get('MATEID', '.')
        ba = alt if svt == 'BND' else '.'
        bo = '.'
        if svt == 'BND':
            if alt.startswith('['): bo = '5to3'
            elif alt.endswith(']'): bo = '3to3' if alt[0] == ']' else '3to5'
            elif '[' in alt: bo = '5to5' if alt.index('[') > 0 else '5to3'

        try: q = float(qs)
        except: q = 0.0
        try: asl = abs(int(slr))
        except:
            try: asl = abs(int(ev) - ps)
            except: asl = 0

        # Parse FORMAT
        fk = fld[8].split(':')
        pi = fk.index('PR') if 'PR' in fk else -1
        si = fk.index('SR') if 'SR' in fk else -1
        gi = fk.index('GT') if 'GT' in fk else 0
        ns = len(snames)

        # Accumulate per-sample evidence
        car = het = hom = spr = ssr = mxp = mxs = 0
        calts = []; nev = nsev = n_sr_car = 0

        for sd in fld[9:]:
            sv = sd.split(':')
            gt = sv[gi] if gi < len(sv) else './.'
            if gt in ('./.', '0/0', './0', '0/.', '.', '0|0'): continue
            car += 1
            if gt in ('1/1', '1|1'): hom += 1
            else: het += 1

            pa = sa = 0
            if pi >= 0 and pi < len(sv):
                try: pa = int(sv[pi].split(',')[1])
                except: pass
            if si >= 0 and si < len(sv):
                try: sa = int(sv[si].split(',')[1])
                except: pass

            ta = pa + sa
            spr += pa; ssr += sa
            if pa > mxp: mxp = pa
            if sa > mxs: mxs = sa
            calts.append(ta)
            if ta >= 1: nev += 1
            if ta >= 3: nsev += 1
            if sa >= 1: n_sr_car += 1

        st = spr + ssr
        cf = car / ns if ns else 0
        af = (het + 2*hom) / (2*ns) if ns else 0
        mna = min(calts) if calts else 0
        mda = sorted(calts)[len(calts)//2] if calts else 0
        efr = nev / ns if ns else 0

        # Population confidence
        pss = binom_log10p(nev, ns, 0.05)
        clevel = confidence_from_score(pss)

        # Breakpoint precision
        cipw = cipos_to_width(cip)
        pclass = precision_class(cipw, n_sr_car)

        # Size class
        if asl < 100: sc = "<100bp"
        elif asl < 1000: sc = "100bp-1kb"
        elif asl < 10000: sc = "1-10kb"
        elif asl < 50000: sc = "10-50kb"
        elif asl < 100000: sc = "50-100kb"
        elif asl < 500000: sc = "100-500kb"
        elif asl < 1000000: sc = "500kb-1Mb"
        else: sc = ">1Mb"

        rows.append([
            ch, str(ps), ev, vid, svt, slr, f"{q:.1f}", flt,
            "IMPRECISE" if imp else "PRECISE", str(cip), str(cie), str(hl), str(hs),
            str(asl),
            str(ns), str(car), f"{cf:.4f}", str(het), str(hom), f"{af:.6f}",
            str(spr), str(ssr), str(st), str(mxp), str(mxs), str(mna), str(mda),
            str(nev), str(nsev), f"{efr:.4f}", f"{pss:.2f}", clevel,
            str(cipw), str(n_sr_car), pclass,
            sc,
            "Y" if i3 else ".", "Y" if i5 else ".", jt, eid,
            str(len(li)) if li else '.', str(len(ri)) if ri else '.',
            mid, ba, bo,
        ])

with open(out_path, 'w') as o:
    o.write('\t'.join(header) + '\n')
    for r in rows: o.write('\t'.join(r) + '\n')

# Summary
print(f"  {out_path}: {len(rows)} PASS variants annotated")
clev = Counter(r[31] for r in rows)
for c in ["VERY_HIGH", "HIGH", "MEDIUM", "LOW", "MINIMAL"]:
    if clev.get(c, 0) > 0: print(f"    {c}: {clev[c]}")
prec = Counter(r[34] for r in rows)
for p in ["base_pair", "narrow", "moderate", "wide", "unknown"]:
    if prec.get(p, 0) > 0: print(f"    precision={p}: {prec[p]}")
PYEOF
}

# =============================================================================
# Run on all types × both cohorts
# =============================================================================
for svtype in "${SV_TYPES[@]}"; do
  mv_log ""
  mv_log "━━━ ${svtype} ━━━"

  for cohort_n in 226 81; do
    vcf="${DIR_FINAL}/catalog_${cohort_n}.${svtype}.PASS.vcf.gz"
    [[ -f "${vcf}" ]] || continue
    n=$(bcftools view -H "${vcf}" 2>/dev/null | wc -l)
    [[ ${n} -eq 0 ]] && continue

    mv_log "  ${cohort_n}: annotating ${n} PASS ${svtype}"
    annotate_vcf "${vcf}" \
      "${ANNOT_DIR}/${cohort_n}.${svtype}.population_confidence.tsv" \
      "${svtype}"
  done
done

# =============================================================================
# Grand summary
# =============================================================================
mv_log ""
mv_log "━━━ CONFIDENCE DISTRIBUTION (226 cohort) ━━━"

GRAND="${ANNOT_DIR}/confidence_summary.tsv"
{
  echo -e "cohort\tsvtype\ttotal\tVERY_HIGH\tHIGH\tMEDIUM\tLOW\tMINIMAL\tbase_pair\tnarrow\tmoderate\twide\tunknown"
  for svtype in "${SV_TYPES[@]}"; do
    for cohort_n in 226 81; do
      tsv="${ANNOT_DIR}/${cohort_n}.${svtype}.population_confidence.tsv"
      [[ -f "${tsv}" ]] || continue
      total=$(tail -n +2 "${tsv}" | wc -l)
      vh=$(awk -F'\t' 'NR>1 && $32=="VERY_HIGH"' "${tsv}" | wc -l)
      hi=$(awk -F'\t' 'NR>1 && $32=="HIGH"' "${tsv}" | wc -l)
      md=$(awk -F'\t' 'NR>1 && $32=="MEDIUM"' "${tsv}" | wc -l)
      lo=$(awk -F'\t' 'NR>1 && $32=="LOW"' "${tsv}" | wc -l)
      mn=$(awk -F'\t' 'NR>1 && $32=="MINIMAL"' "${tsv}" | wc -l)
      bp=$(awk -F'\t' 'NR>1 && $35=="base_pair"' "${tsv}" | wc -l)
      nr=$(awk -F'\t' 'NR>1 && $35=="narrow"' "${tsv}" | wc -l)
      mo=$(awk -F'\t' 'NR>1 && $35=="moderate"' "${tsv}" | wc -l)
      wd=$(awk -F'\t' 'NR>1 && $35=="wide"' "${tsv}" | wc -l)
      uk=$(awk -F'\t' 'NR>1 && $35=="unknown"' "${tsv}" | wc -l)
      echo -e "${cohort_n}\t${svtype}\t${total}\t${vh}\t${hi}\t${md}\t${lo}\t${mn}\t${bp}\t${nr}\t${mo}\t${wd}\t${uk}"
    done
  done
} > "${GRAND}"

column -t "${GRAND}" | while IFS= read -r line; do mv_log "  ${line}"; done

mv_log ""
mv_log "=== ANNOTATION COMPLETE ==="
mv_log ""
mv_log "Output: ${ANNOT_DIR}/"
mv_log ""
mv_log "Confidence levels (from pop_signal_score = -log10 binomial P):"
mv_log "  VERY_HIGH : score ≥ 10  (P < 1e-10)  — genomic certainty"
mv_log "  HIGH      : score ≥ 5   (P < 1e-5)   — strong population signal"
mv_log "  MEDIUM    : score ≥ 2   (P < 0.01)   — emerging signal"
mv_log "  LOW       : score ≥ 1.3 (P < 0.05)   — suggestive"
mv_log "  MINIMAL   : score < 1.3 (P ≥ 0.05)   — weak / singleton"
mv_log ""
mv_log "Breakpoint precision (from CIPOS width + split-read carrier count):"
mv_log "  base_pair : CIPOS ≤5bp + ≥3 split-read carriers"
mv_log "  narrow    : CIPOS ≤50bp + ≥1 split-read carrier"
mv_log "  moderate  : CIPOS ≤200bp"
mv_log "  wide      : CIPOS >200bp"
mv_log "  unknown   : no CIPOS available"
mv_log ""
mv_log "Use case: when the inversion pipeline finds a boundary, look up"
mv_log "nearby SVs in these tables. Pick the one with highest confidence"
mv_log "and best precision. This doesn't filter — it RANKS."
