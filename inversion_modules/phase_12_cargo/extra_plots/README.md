# phase_12_cargo/extra_plots — supplementary visualisations

A collection of standalone R scripts producing manuscript / supplementary
figures from outputs that already exist on disk. **No new compute** — every
script is a pure read-and-plot consumer.

Conceptually adjacent to MODULE_6_Cargo but logically independent:
- The cargo module produces inversion-by-inversion biological narratives.
- These extra plots produce cohort-wide and per-candidate visualisations
  (variant landscape, sample networks, inversion sharing patterns).

## Layout

```
extra_plots/
├── README.md                                 (this file)
├── 00_extra_plots_config.sh                  (sources cargo config + extras)
├── compute/                                  (R scripts — pure read+plot)
│   ├── PLOT_01_karyotype_network.R           per-candidate sample IBS network coloured by karyotype
│   ├── PLOT_02_inversion_sharing_across_groups.R   inversion × group HOM_INV-frequency heatmap
│   ├── PLOT_03_variant_burden_dashboard.R    per-sample / per-group / per-chrom variant counts
│   ├── PLOT_04_cohort_sfs.R                  cohort-wide variant frequency spectrum by class
│   ├── PLOT_05_top_recurrent_variants.R      top-N most-shared variants (table + barchart)
│   ├── PLOT_06_private_vs_shared_by_group.R  per-group private/shared stacked bars
│   ├── PLOT_07_group_exclusivity_heatmap.R   carrier-frequency heatmap of group-specific variants
│   ├── PLOT_08_upset_variant_sharing.R       set-intersection of carrier groups
│   ├── PLOT_09_genomic_distribution.R        stacked Manhattan of private vs shared (no circlize)
│   ├── PLOT_10_context_enrichment.R          variant fold-enrichment by genomic context
│   ├── PLOT_11_sv_in_inversions.R            log2FC of SV classes inside vs outside inversions
│   └── PLOT_12_inversion_pop_frequencies.R   distribution of inversion carrier frequencies
├── tables/                                   (R scripts producing standalone TSV tables)
│   ├── TABLE_01_top_recurrent_with_group_AF.R    image 1B — adds per-group max-AF column
│   ├── TABLE_02_top_private_singletons.R         image 2C — singleton variants with sample/group
│   ├── TABLE_03_variant_categories_summary.R     image 2 footer — Private/Rare/Common/High 4-row
│   ├── TABLE_04_group_exclusivity_summary.R      image 2 footer — per-group exclusivity statistics
│   ├── TABLE_05_repeat_class_enrichment.R        image 6C — repeat class fold enrichment
│   └── TABLE_06_top_inversion_carriers.R         image 7B — clean ranked inversion table
└── launchers/
    ├── PLOT_01_karyotype_network.slurm
    ├── PLOT_02_inversion_sharing.slurm
    ├── PLOT_03_variant_burden_dashboard.slurm
    ├── run_all_extra_plots.slurm             one driver, runs all 12 plots in sequence
    └── run_all_extra_tables.slurm            one driver, runs all 6 standalone tables
```

## Two output streams

**`${EXTRAS_FIG_DIR}/`** — PDF figures (one per PLOT_* script).

**`${EXTRAS_TBL_DIR}/`** — TSV tables. Two kinds land here:
1. *Plot-companion TSVs* (named `PLOT_NN_*.tsv`) — every PLOT_NN script also emits its underlying data table, so plots are reproducible.
2. *Standalone tables* (named `TABLE_NN_*.tsv`) — manuscript-ready tables that summarize across plots or extract additional columns not visible in any single figure.

## Common conventions

- All scripts source `00_extra_plots_config.sh` (which itself sources
  `00_cargo_config.sh` → `00_inversion_config.sh`).
- Output PDFs land under `${EXTRAS_FIG_DIR}` (default `${INVDIR}/30_cargo/extras/`).
- Per-candidate scripts take `[candidate_id|all]` as a CLI argument.
- Each script is fail-safe: if its required upstream input is missing, it
  emits a `[skip]` message and exits 0 — so `run_all_extra_plots.slurm`
  doesn't break when only some inputs exist.

## What's deliberately NOT here

- **MODULE_6_Cargo content** (per-arrangement burden, GO/family enrichment,
  configuration spectrum) lives in the cargo module itself; this folder
  doesn't duplicate it.
- **Regional burden Manhattan / new windowed scans** — that's new compute,
  not a plot. If you want it, it's a separate module.
- **`circlize` circos plots** — replaced with stacked Manhattan
  (`PLOT_09_genomic_distribution.R`). Same information, no extra package
  dependency, scales better at 226 samples × 28 chromosomes.

## Required inputs

Per script, the upstream files consumed:

| Plot | Reads |
|---|---|
| 01 | `${DOSAGE_DIR}/<chr>.dosage.tsv.gz`, `${INVDIR}/<cand_dir>/data/sample_karyotypes.tsv` |
| 02 | `${SAMPLE_REGISTRY}/groups/inv_*_HOM_INV.txt`, ancestry Q at canonical K |
| 03 | `${VARIANT_MASTER}` — sample/chrom carrier columns |
| 04 | `${VARIANT_MASTER}` — carrier_count column |
| 05 | `${VARIANT_MASTER}` — sort by n_carriers |
| 06 | `${VARIANT_MASTER}` + ancestry group memberships |
| 07 | Same as 06 |
| 08 | Same as 06 |
| 09 | `${VARIANT_MASTER}` + `${REF_FAI}` |
| 10 | `${VARIANT_MASTER}` + `${GENE_BED}` + `${REPEAT_BED}` (optional) + ROH BED (optional) |
| 11 | DELLY / Manta SV outputs at `${BASE}/sv_calling/...` + `${SNAKE_CAND_FILE}` |
| 12 | `${SNAKE_CAND_FILE}` + per-candidate karyotype counts |
