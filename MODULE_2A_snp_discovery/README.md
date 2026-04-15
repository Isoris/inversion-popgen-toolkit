# MODULE_2B — Variants & Ancestry
## STEP_1_Variants_StructureBase

Core modern-DNA variant foundation: callable mask, ANGSD site discovery, biSNP
panels, BEAGLE generation, PCAngsd, NGSadmix, evalAdmix, ngsRelate, relatedness
pruning, and canonical best-seed-by-K selection.

STEP_1 stands alone as a complete variant + structure pipeline.
STEP_2 (Multiscale Q Extractor) is optional and consumes STEP_1 outputs.

## Public Interface (2 files)

| File | Purpose |
|------|---------|
| `run_step1.sh` | Main wrapper — the ONLY script you call directly |
| `config.sh` | All paths, ANGSD params, K range, seeds, cutoffs, palette |

## Subcommands

```bash
./run_step1.sh init               # Validate config, tools, inputs
./run_step1.sh masks_sfs          # Callable mask → RF chunks → SAF → SFS → pest
./run_step1.sh call_bisnps        # biSNP discovery per chunk using pest prior
./run_step1.sh build_panels       # Merge + thin at 6 distances + ANGSD-index
./run_step1.sh make_beagles       # BEAGLE GL per-RF + whole-genome + merge
./run_step1.sh structure_all      # PCAngsd + NGSadmix + evalAdmix + best-seed (ALL samples)
./run_step1.sh relatedness        # ngsRelate + NAToRA + pruning → pruned_samples.txt
./run_step1.sh structure_pruned   # Same structure workflow on PRUNED samples
./run_step1.sh merge_summaries    # Combined best-seed table (sample_set = all | pruned)
./run_step1.sh all                # Run everything in order
```

## Practical Run Order

1. `masks_sfs` — build masks, SAF, folded SFS, pest prior
2. `call_bisnps` — ANGSD biSNP calling per chromosome
3. `build_panels` — merge + thin (200/500/1000 bp, 5/10/25 kb)
4. `make_beagles` — BEAGLE GL generation + merging
5. `structure_all` — PCAngsd + NGSadmix + evalAdmix + best-seed on all 226 samples
6. `relatedness` — ngsRelate + NAToRA multi-cutoff + greedy pruning
7. `structure_pruned` — same as step 5 but on pruned unrelated sample set
8. `merge_summaries` — combined_best_seed_by_K.tsv with sample_set column

## Canonical STEP_1 Outputs (per sample set)

| File | Description |
|------|-------------|
| `all_seed_metrics_by_K.tsv` | Every K × seed: loglik, evalAdmix residuals, file paths |
| `best_seed_by_K.tsv` | One row per K: selected best seed + metrics |
| `best_seed_copied_files.tsv` | Provenance: which files were copied as `*_best` |
| `cluster_palette_by_K.tsv` | Stable 12-color palette per K × cluster |
| `sample_order_reference.tsv` | Canonical sample ordering |
| `sample_main_ancestry_by_K.tsv` | Per-sample: dominant cluster, full Q, color per K |
| `combined_best_seed_by_K.tsv` | Merged all + pruned with sample_set column |

## Hidden Helpers (in `helpers/`)

Orchestrators: `01_masks_sfs.sh`, `02_bisnps.sh`, `03_panels.sh`, `04_beagles.sh`,
`05_structure.sh` (the unified structure script), `06_relatedness.sh`,
`merge_beagles.sh`, `merge_structure_summaries.sh`

SLURM jobs: `slurm_saf_chunk.sh`, `slurm_merge_sfs.sh`, `slurm_bisnp_chunk.sh`,
`slurm_beagle_rf.sh`, `slurm_beagle_wg.sh`, `slurm_pcangsd.sh`,
`slurm_ngsadmix.sh`, `slurm_ngsadmix_global.sh`

Tools (preserved): `mask_regions_from_fasta.py`, `select_best_seed_by_K.R`,
`prune_first_degree_pairs.py`, `NAToRA_Public.py`, `plot_relatedness_3panel.R`

## Key Design: structure_all / relatedness / structure_pruned

`05_structure.sh` is ONE script that takes `--samples <file>`. The wrapper calls
it twice — once with all samples, once with pruned samples. Same logic, same
outputs, different sample list. The `merge_structure_summaries.sh` script then
combines both into one table with `sample_set = all | pruned`.

## ANGSD Shared Filter Parameters

| Parameter | Value |
|-----------|-------|
| -GL | 1 (SAMtools) |
| -minQ | 25 |
| -minMapQ | 25 |
| -baq | 1 |
| -C | 50 |
| -setMinDepthInd | 3 |
| -setMaxDepthInd | 57 |
| -minInd | 200 |
| -SNP_pval | 1e-6 |
| -minMaf | 0.05 |
