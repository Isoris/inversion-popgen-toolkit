# 00_markdup — Duplicate-Marked BAMs for SV Callers

## What this is

Duplicate-marked BAMs produced from the raw minimap2 alignments (MODULE_1 manifest column 2). These are the input BAMs for all structural variant callers:

- MODULE_4B (DELLY DEL)
- MODULE_4C (DELLY DUP)
- MODULE_4D (DELLY INV)
- MODULE_4E (DELLY INS)
- MODULE_4F (DELLY BND)
- MODULE_4H (Manta)

## Why separate from MODULE_1 filtered BAMs

SV callers like DELLY and Manta require **unfiltered, duplicate-marked** BAMs — not the population-genomics-filtered BAMs from MODULE_1 (which apply MAPQ≥60, proper-pair, same-chr, TLEN≤p99 filters). The aggressive filtering in MODULE_1 removes the discordant pairs and split reads that SV callers depend on for breakpoint detection.

## How they were produced

```bash
# For each sample (from MODULE_4B STEP_A01):
samtools markdup -@ 4 <raw_minimap2.bam> <sample>.markdup.bam
samtools index <sample>.markdup.bam
```

Source BAMs: MODULE_1 manifest column 2 (`sample_bam_minimap2_vs_P99TLENMAPQ30.tsv`).

Duplicate marking strategy: mark-only (not removed), matching MODULE_1's approach.

## File naming

```
00_markdup/
  CGA001.markdup.bam
  CGA001.markdup.bam.bai
  CGA002.markdup.bam
  ...
  (226 samples)
```

## HPC location

```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/delly_sv/00_markdup/
```

All MODULE_4* configs reference this path via `DIR_MARKDUP`.
