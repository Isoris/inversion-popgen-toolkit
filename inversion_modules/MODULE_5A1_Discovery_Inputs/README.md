# MODULE_5A1 — Discovery Inputs

## Purpose

Upstream data preparation for the inversion discovery pipeline. Builds all inputs needed before the core local-structure discovery engine can run.

## Scope

STEP01 through STEP07:

| Step | Script | What it does |
|------|--------|-------------|
| 01 | `STEP01_mask_regions_from_fasta.py` | Identify non-ACGT mask regions in reference |
| 02 | `STEP02_make_inversion_callable_angsd_mask.sh` | Build ANGSD-compatible callable mask |
| 03 | `STEP03_mask_depth_mapq_stats.py` | Depth/mapQ QC stats per sample |
| 04 | `LAUNCH_STEP04_...slurm` | ANGSD SAF estimation per chunk (SLURM array) |
| 05 | `LAUNCH_STEP05_...slurm` | Merge chunk SAFs to global |
| 06 | `LAUNCH_STEP06_...slurm` | Global folded SFS + mean pest |
| 07 | `LAUNCH_STEP07_...slurm` | ANGSD SNP calling → BEAGLE GL files |

## Main Inputs

- Reference genome (`fClaHyb_Gar_LG.fa`)
- BAM list (`bamlist_qcpass.txt`)
- Chunk RF lists for SLURM arrays

## Main Outputs

- Callable mask files (`.angsd.idx`)
- Per-chunk SAF files
- Global folded SFS
- Per-chromosome BEAGLE genotype likelihood files (`.beagle.gz`)

## Runner

```bash
bash runners/run_5A1_discovery_inputs.sh              # all steps
bash runners/run_5A1_discovery_inputs.sh --step 03    # single step
```

## Dependencies

- ANGSD, samtools in conda `assembly` environment
- `00_inversion_config.sh` sourced for all paths

## Downstream

Feeds into MODULE_5A2 (STEP08 reads the BEAGLE files produced here).
