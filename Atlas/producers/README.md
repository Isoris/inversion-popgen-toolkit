# Producers — SV evidence layers for one candidate

This directory contains three independent Python producers that emit the JSON layers consumed by the atlas's SV evidence page (page `5b` in `Inversion_atlas.html`). All three write into the same per-candidate folder structure:

```
data/
└── <chrom>/
    └── candidates/
        └── <cand_id>/
            ├── sv_genotype_counts.json        (STEP_SV_GT_AGG)
            ├── sv_evidence_combinations.json  (STEP_SV_EVID_COMB)
            └── sv_support_by_sample.json      (STEP_SV_SUPPORT)
```

The atlas walks this folder via `webkitGetAsEntry()` on drag-drop and dispatches each `*.json` by its `format_version`. Each producer can ship independently — partial folders are fine; the atlas just renders fewer panels.

## STEP_SV_GT_AGG — sv_genotype_counts_v1.json

**Inputs:** merged DELLY+Manta VCF, locked karyotype labels, candidate metadata.

**Pipeline:**
1. Parse VCF → per-sample dosage 0/1/2/-1
2. **`min_carriers_per_band: 5` filter** (Quentin's hard threshold — drop SVs where no karyotype band has ≥5 carriers)
3. Per group counts: AA/AB/BB/miss
4. Fisher exact (H1/H1 vs H2/H2) + Benjamini-Hochberg FDR
5. Pattern label classification (canonical_breakpoint_marker, dominant_presence_marker, het_specific_marker, sub_haplotype_marker, internal_linked_marker, uninformative)
6. Boundary summary aggregation

**Atlas consumes:** main SV table + locus track + boundary-summary right-rail tables.

## STEP_SV_EVID_COMB — sv_evidence_combinations_v1.json

**Inputs:** BAM-evidence track from S7 / phase_8_comparative_breakpoint_fragility, DELLY+Manta GT-call status, locked labels.

**Pipeline:**
1. Per sample, build evidence vector: `{left_SA, right_SA, left_PE, right_PE, Manta_INV_GT, DELLY_INV_GT, MAPQ0_left, MAPQ0_right}` → bools
2. Group samples by identical vector
3. Sort by intersection size desc, top N (default 30)
4. Compute per-evidence-type totals (for the right-side mini-bars)

**Atlas consumes:** UpSet panel (right-rail). Click a bar → `selectedSamples` populates → locus dims non-selected SV glyphs (when sv_support layer also present).

## STEP_SV_SUPPORT — sv_support_by_sample_v1.json

**Inputs:** same VCF as STEP_SV_GT_AGG, the SV list that survived filtering, locked labels.

**Output:** per-sample × per-SV dosage matrix in compact-string row encoding (`'0'/'1'/'2'/'.'` per cell, one string per sample row). 226 × 50 SVs ≈ 11 KB.

**Atlas consumes:** dosage heatmap (right-rail compact + main-area large overlay) + the "dim non-selected SV glyphs" feature (steps 5+6 cross-talk).

## Running them together

```bash
# One candidate, three layers in one folder:
CAND=cand_LG28_15Mb
CHROM=C_gar_LG28
LEFT=15142000
RIGHT=18124000
LOCKS=locks/${CAND}.json
OUT=data/

python producers/STEP_SV_GT_AGG_aggregate_genotype_counts.py \
    --vcf merged_delly_manta.vcf.gz \
    --locks $LOCKS \
    --candidate-id $CAND --chrom $CHROM \
    --left $LEFT --right $RIGHT \
    --out-root $OUT

# Extract the surviving sv_ids from the gt_counts emitted above
jq -r '.sv_calls[].sv_id' \
    data/$CHROM/candidates/$CAND/sv_genotype_counts.json \
    > /tmp/${CAND}__sv_ids.txt

python producers/STEP_SV_SUPPORT_emit_support_by_sample.py \
    --vcf merged_delly_manta.vcf.gz \
    --sv-list /tmp/${CAND}__sv_ids.txt \
    --locks $LOCKS \
    --candidate-id $CAND --chrom $CHROM \
    --out-root $OUT

python producers/STEP_SV_EVID_COMB_emit_combinations.py \
    --bam-evidence s7_phase8/${CAND}.json \
    --vcf-status   merged_delly_manta_status.json \
    --locks $LOCKS \
    --candidate-id $CAND --chrom $CHROM \
    --left $LEFT --right $RIGHT \
    --out-root $OUT

# Drop the folder data/$CHROM/candidates/$CAND/ onto the atlas to see all three
# panels populate.
```

## Implementation status

All three are **STUBS** — the JSON-shape construction is complete and matches what the atlas validates against (locked schemas), but the **VCF parsing** and **BAM-evidence extraction** blocks are placeholders. Replace with `cyvcf2`/`pysam` (VCF) and your existing S7 BAM-evidence module's API.

Once those blocks are filled in, you can run them on LANTA against the 226-sample cohort and the per-candidate folders will populate.

## Notes on the `min_carriers_per_band: 5` filter

This is **producer-side, intentional**, and applies to STEP_SV_GT_AGG. Reasons:

- **K=3 (REF/HET/INV) on 226 samples** — bands of ~60 typically; n=5 is a robust floor for "this SV actually segregates with the band". Below that, it's noise.
- **Pre-filtering keeps the JSON small.** For ~1000 raw SV calls per candidate, the filter often drops half. The atlas's own `min_samples` UI filter then lets you dial up further.
- **The atlas defaults to `min_samples: 1`** — it doesn't second-guess the producer. The atlas trusts that anything in the JSON has been pre-vetted.

For K=6 (substructure), the same threshold works empirically — bands can be ~30-40 samples and n=5 is still a reasonable floor. If you want a stricter or looser cut for K=6, override via the `MIN_CARRIERS_PER_BAND` constant in `STEP_SV_GT_AGG_aggregate_genotype_counts.py`. Atlas-side, no change is needed.

The constant is recorded in the emitted JSON under `min_carriers_filter` so any downstream reader (or a future paper supplement) can reproduce the cut.
