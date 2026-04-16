# MODULE_1_read_prep

Read-level QC, species verification, short-read alignment, BAM merging and deduplication, population-genomics filtering, and depth QC. Produces the filtered BAMs consumed by every downstream module.

## Pipeline

```
Raw FASTQ (PE151 Illumina)
  │
  ├─ A. QC & trim ────────────────────────────────────────────────────
  │   S01  fastp trim (redo problematic pairs from EOF table)
  │   S02  scan for corrupted gzip → build redo table
  │   S03  parse fastp logs → per-sample QC table + summary
  │
  ├─ B. Species verification ─────────────────────────────────────────
  │   S04  Mash distance (k=31, s=50k) → Gar vs Mac per run-lane
  │   S05  meryl k-mer (k=21) independent verification
  │   S06  plot species assignment (scatter, hist, suspicious ranking)
  │
  ├─ C. Alignment ────────────────────────────────────────────────────
  │   S07  inventory check: fastp outputs ↔ BAM completeness
  │   S08  BAM/BAI pair consistency
  │   S09  minimap2 -ax sr → per-run sorted BAM (batched)
  │
  ├─ D. Merge & dedup ────────────────────────────────────────────────
  │   S10  merge per-run → per-sample, fixmate, markdup, clipOverlap
  │   S11  TLEN p95/p99 from merged BAM (MAPQ≥60, proper pair)
  │   S12  detailed insert size stats (raw + MAD-trimmed)
  │
  ├─ E. PopGen filter & QC ───────────────────────────────────────────
  │   S13  filter: MAPQ≥60, proper pair, same-chr, TLEN≤p99
  │   S14  mosdepth depth QC (5 region classes, per-sample, per-chr)
  │   S15  summarize BAM QC → ANGSD depth cutoff suggestions
  │   S16  provenance table: sample → raw BAM → filtered BAM
  │
  └─→ 226 filtered BAMs (~9× mean coverage)
        consumed by MODULE_2A, MODULE_3, MODULE_4A, MODULE_4B–4H
```

## Scripts

| Step | Script | Type | What it does |
|------|--------|------|-------------|
| S01 | `S01_fastp_trim_reads.slurm` | SLURM | fastp on TSV-defined read pairs. `--detect_adapter_for_pe`, `--trim_poly_g`, `-q 20`, `-f 5 -F 5 -t 5 -T 5`, `-n 0`, `-l 30`. Parallel per-sample on full node. |
| S02 | `S02_fastp_eof_redo_table.sh` | bash | Scan fastp outputs for unexpected EOF → `fastp_unexpected_eof_redo.tsv`. Searches uploadku for matching raw reads. |
| S03 | `S03_extract_fastp_qc.sh` | bash | Parse fastp logs → `fastp_qc_table.tsv` (per-sample) + `fastp_qc_summary.tsv` (min/mean/median/max). |
| S04 | `S04_species_assign_mash.slurm` | SLURM | Build Mash sketches for Gar and Mac refs (k=31, s=50000, min copies=2). Assign each run-lane by distance margin 0.002. Per-CGA majority vote. |
| S05 | `S05_species_assign_kmers.sh` | bash | Standalone Mash + meryl k-mer species assignment (k=21, log₂ ratio threshold=1.0). CLI tool, not SLURM-specific. |
| S06 | `S06_plot_species_assign.py` | python | Scatter of Mash distances, histogram of delta, suspicious-sample ranking. Outputs PNG/PDF + TSV. |
| S07 | `S07_inventory_fastp_vs_bam.sh` | bash | Cross-check: does every fastp output have a matching BAM? Reports complete/missing/inconsistent. |
| S08 | `S08_check_bam_bai_pairs.sh` | bash | Verify every `.sr.bam` has a `.sr.bam.bai` and vice versa. |
| S09 | `S09_map_minimap2.slurm` | SLURM | minimap2 `-ax sr -Y` with RG per run-lane (SM=sample, LB=sample, PU=run.lane, PL=ILLUMINA). Sort + index. |
| S10 | `S10_merge_markdup_clip.slurm` | SLURM array | Merge per-run BAMs, validate RG/SM consistency, fixmate, markdup (mark-only), BamUtil clipOverlap, flagstat + stats. |
| S11 | `S11_get_tlen_percentiles.sh` | bash | TLEN p95/p99 from properly paired, primary, non-dup reads with MAPQ≥60 and same-chr mate. |
| S12 | `S12_get_insert_size_stats.sh` | bash | Raw + MAD-trimmed (median ± 10×MAD) insert size stats. PE151 layout detection. |
| S13 | `S13_filter_bam_popgen.slurm` | SLURM array | `-f 0x2 -F 0xF0C -q 60` + same-chr mate + abs(TLEN) ≤ p99 (514 bp). Final popgen BAMs. |
| S14 | `S14_qc_depth_table.slurm` | SLURM | mosdepth across 5 region classes (whole-chr, repeat, nonrepeat, callable, noncallable). Per-sample + per-chr tables. Thresholds at 1×/5×/10×. |
| S15 | `S15_summarize_bam_qc.sh` | SLURM | Parse flagstat + stats → `all_samples.bam_qc.tsv` + `summary.txt`. Heuristic ANGSD depth cutoffs (minDepthInd, maxDepthInd). |
| S16 | `S16_make_bam_provenance.sh` | bash | Build `sample_bam_minimap2_vs_P99TLENMAPQ30.tsv` linking each CGA sample to raw + filtered BAM paths. Definitive manifest for all downstream. |

## Key outputs

| Output | Description | Used by |
|--------|-------------|---------|
| `{sample}.filtered.bam` | PopGen-filtered BAMs (226 samples) | MODULE_2A, MODULE_3, MODULE_4A, unified_ancestry |
| `{sample}.markdup.bam` | Markdup BAMs (shared, not filtered) | MODULE_4B–4H (DELLY/Manta need unfiltered markdup) |
| `sample_bam_minimap2_vs_P99TLENMAPQ30.tsv` | Provenance manifest | All modules (sample list) |
| `all_samples.coverage_context_qc.tsv` | Per-sample depth across 5 region classes | MODULE_2A (ANGSD depth filters) |
| `fastp_qc_table.tsv` | Per-sample read QC | Manuscript Table S1 |

## Parameters

| Parameter | Value | Script |
|-----------|-------|--------|
| fastp quality trim | Phred ≥ 20 | S01 |
| fastp end trim | 5 bp each end, both reads | S01 |
| fastp min length | 30 bp | S01 |
| Mash k / sketch / margin | 31 / 50000 / 0.002 | S04 |
| Aligner | minimap2 `-ax sr` | S09 |
| Dedup strategy | samtools markdup (mark-only, not removed) | S10 |
| Overlap clip | BamUtil clipOverlap `--storeOrig OC` | S10 |
| PopGen MAPQ filter | ≥ 60 | S13 |
| PopGen flag filter | `-f 0x2 -F 0xF0C` (proper pair, no unmapped/sec/dup/supp) | S13 |
| PopGen TLEN filter | ≤ sample-specific p99 (514 bp) | S13 |
| Depth QC | mosdepth, 5 region classes | S14 |
| Callable genome | 963,905,721 bp | S14–S15 |

## Dependencies

fastp (≥0.23), minimap2, samtools, BamUtil (clipOverlap), mosdepth, Mash, meryl (optional for S05), python3 + matplotlib (for S06)

## Sidecar convention

Every script writes two audit files alongside its outputs:
- **`{step}.arg`** — all parameters, tool versions, paths, exact commands (written at start)
- **`{step}.results`** — all output file paths + descriptions (written at end)

This makes every output traceable to the exact invocation that produced it.

## Methods

See [`MODULE_1_Reads_methods.md`](MODULE_1_Reads_methods.md) for manuscript-ready methods prose with all parameter values and citations.
