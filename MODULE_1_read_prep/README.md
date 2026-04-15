# MODULE_1_Reads

Read-level QC, species verification, alignment, BAM processing, and final filtered BAM preparation for downstream population genomics (ANGSD, PCAngsd, NGSadmix, etc.).

## Execution Order

| Step | Script | Purpose | SLURM? |
|------|--------|---------|--------|
| S01 | `S01_fastp_trim_reads.slurm` | Trim/redo fastp on problematic read pairs | Yes |
| S02 | `S02_fastp_eof_redo_table.sh` | Find corrupted fastp outputs, build redo table | No |
| S03 | `S03_extract_fastp_qc.sh` | Parse fastp logs → per-sample QC table + summary | No |
| S04 | `S04_species_assign_mash.slurm` | Mash-based species assignment (Gar vs Mac) per run-lane | Yes |
| S05 | `S05_species_assign_kmers.sh` | Combined Mash + meryl k-mer species assignment | No |
| S06 | `S06_plot_species_assign.py` | Plot species assignment results (scatter, histogram, suspicious) | No |
| S07 | `S07_inventory_fastp_vs_bam.sh` | Cross-check fastp outputs against BAM inventory | No |
| S08 | `S08_check_bam_bai_pairs.sh` | Verify BAM/BAI consistency across batches | No |
| S09 | `S09_map_minimap2.slurm` | Map fastp reads with minimap2, per-batch | Yes |
| S10 | `S10_merge_markdup_clip.slurm` | Merge per-run BAMs → per-sample, fixmate, markdup, clipOverlap, QC | Yes (array) |
| S11 | `S11_get_tlen_percentiles.sh` | Compute TLEN p95/p99 from a BAM | No |
| S12 | `S12_get_insert_size_stats.sh` | Detailed insert size statistics (raw + trimmed) | No |
| S13 | `S13_filter_bam_popgen.slurm` | Filter BAMs: MAPQ≥60, proper-pair, same-chr, TLEN≤p99 | Yes (array) |
| S14 | `S14_qc_depth_table.slurm` | Comprehensive depth/coverage QC with mosdepth (per-sample, per-chr, by region class) | Yes |
| S15 | `S15_summarize_bam_qc.sh` | Summarize BAM QC for ANGSD: depth distribution, suggested cutoffs | Yes |
| S16 | `S16_make_bam_provenance.sh` | Build provenance table linking samples to raw + filtered BAMs | No |

## Subdirectories

- `Tables/` — Reserved for output tables referenced in manuscript
- `Figures/` — Reserved for output figures referenced in manuscript
- `Manuscript/` — Methods markdown for this module

## Origin Mapping

| New script | Old script(s) | Notes |
|------------|--------------|-------|
| S01 | STEP01_00_fastp_trim_reads.slurm + .arg | Merged script+arg into single self-documenting script |
| S02 | STEP01_01_make_fastp_eof_redo_table.sh | Cleaned up |
| S03 | STEP02_00_extract_fastp_qc_tables.sh | Cleaned up |
| S04 | STEP03_00_species_assign_mash_fastp.slurm + .arg | Merged script+arg |
| S05 | STEP03_01_species_assign_kmers_mash_meryl.sh + .arg | Merged script+arg |
| S06 | STEP03_plot_mash_species_assign.py | Cleaned up |
| S07 | HELPER_STEP01_compare_fastp_and_bam_inventory.sh | Added set -euo, proper tmp cleanup |
| S08 | STEP03_00_check_sr_bam_bai_pairs.sh | Added arg/results output |
| S09 | STEP03_01_map_minimap2_batch00.slurm | Generalized batch variable |
| S10 | STEP03_02_merge_markdup_clip_qc_array.slurm | Cleaned up |
| S11 | STEP03_03_get_tlen_p95_p99.sh | Added arg/results output |
| S12 | STEP03_04_get_insert_size_stats.sh | Cleaned up |
| S13 | STEP03_04_filter_bam_popgen_pp_samechr_tlenp99.slurm | Cleaned up |
| S14 | STEP03_05_qc_depth_table_final.slurm | Cleaned up |
| S15 | STEP03_06_summarize_bam_qc_for_angsd.sh | Cleaned up |
| S16 | STEP03_07_make_bam_provenance_table.sh | Added arg/results output |

## Design Conventions

Every script that produces outputs writes two sidecar files:

1. **`<step>.arg`** — key-value TSV recording all parameters, tool versions, paths, and the exact command(s) run. Written at the start of execution.
2. **`<step>.results`** — key-value TSV listing all output file paths and a brief description. Written at the end of execution.

This makes the pipeline fully auditable: for any output file you can trace back to the exact parameters and commands that produced it.
