# MODULE_6_Cargo (phase_7_cargo) — integration notes

## Where to drop this in the repo

```
inversion_modules/
├── 00_inversion_config.sh
├── phase_2_discovery/
├── phase_3_refine/
├── phase_4_postprocessing/
├── phase_5_followup/
├── phase_6_secondary/
└── phase_7_cargo/   ← THIS MODULE
    ├── README.md
    ├── 00_cargo_config.sh
    ├── STEP_60_build_gene_intervals.sh
    ├── STEP_61_eggnog_annotate.sh
    ├── compute/
    ├── analysis/
    ├── plot/
    └── launchers/
```

## Required patch before first run

Apply `STAT_ENUM_cargo_v1.patch` (shipped separately under `registries/patches/`)
to extend the `put_interval_summary` stat enum:

```bash
cd ${BASE}
patch -p1 < registries/patches/STAT_ENUM_cargo_v1.patch
```

Without the patch, MODULE_6 will still run — it falls back to writing files
directly under `${CARGO_DIR}/per_arrangement_burden/<cid>/` and emits a WARN
for each rejected `put_interval_summary` call. Apply the patch to get proper
registry rows.

## Required upstream state

Before STEP_C61 runs, these must already exist:

1. **Sample groups for arrangements** registered by `STEP_C01i_d_seal.R`:
   `inv_<cid>_HOM_REF`, `inv_<cid>_HET`, `inv_<cid>_HOM_INV` for each
   candidate.

2. **Per-chromosome cohort VCFs** at `${MODCONS}/03_variants/normalized/${CHR}.clair3.norm.vcf.gz`
   (produced by MODULE_CONSERVATION step 3).

3. **Variant master table** at `${MODCONS}/16_merged_variant_tables/variant_master_scored.tsv`
   with columns `var_key, chr, pos, ref, alt, gene_id, snpeff_annotation,
   snpeff_impact, bcsq_consequence, bcsq_gene, ...`.

4. **VESM scores** at `${MODCONS}/05B_vesm/vesm_variant_scores.tsv`
   with columns `var_key, vesm_llr, ...`.

5. **Candidate intervals** at `${SNAKE_CAND_FILE}` (already wired by
   `00_inversion_config.sh`).

## Optional inputs

- `${REPEAT_BED}` — RepeatMasker / EarlGrey BED of repeats. If unset,
  `repeat_overlap_frac` column in `genes.tsv` is empty.
- eggNOG-mapper output (built by STEP_61) — required for GO/KEGG enrichment.
  If not present, family-only enrichment runs.

## Known column-name dependencies

The cargo scripts assume `variant_master_scored.tsv` has these columns:

| Column | Used by | Notes |
|---|---|---|
| `var_key` | C61, C62 | "chr:pos:ref:alt" key, joined to VCF |
| `chr`, `pos` | C60, C61, C62 | filtering by interval |
| `gene_id` | C60, C61, C62 | falls back to `bcsq_gene` if missing |
| `snpeff_annotation` | C60, C61, C62 | classified into missense / synonymous / LoF |
| `bcsq_consequence` | C60, C61, C62 | secondary annotation source |

If MODULE_CONSERVATION's column names differ, the simplest fix is to add an
alias map at the top of each script (search `v.get('chr', '')` and friends).

## Output layout under `${INVDIR}/30_cargo/`

```
30_cargo/
├── gene_intervals.bed.gz, gene_intervals.bed.gz.tbi
├── gene_lengths.tsv
├── transcript_gene_map.tsv
├── eggnog_annotations.tsv      (optional)
├── gene_function.tsv           (per-gene resolved annotation)
├── inventory/
│   ├── diagnostic_table.tsv    ← THE GATING TABLE
│   └── <cid>/
│       ├── genes.tsv
│       └── paralog_summary.tsv
├── per_arrangement_burden/<cid>/   (only if registry stat-enum patch not applied)
│   ├── per_gene_summary__inv_<cid>_HOM_REF.tsv
│   ├── per_gene_summary__inv_<cid>_HOM_INV.tsv
│   ├── sfs__inv_<cid>_HOM_REF.tsv
│   ├── sfs__inv_<cid>_HOM_INV.tsv
│   ├── fst_per_site.tsv
│   └── burden_diff.tsv
├── config_spectrum/<cid>/      (only if registry stat-enum patch not applied)
│   ├── config_mi__inv_<cid>_HOM_REF.tsv
│   └── config_mi__inv_<cid>_HOM_INV.tsv
├── functional_enrichment/<cid>/
│   ├── go_enrichment.tsv
│   ├── kegg_enrichment.tsv
│   └── family_enrichment.tsv
├── cross_inversion/
│   ├── family_convergence.tsv
│   ├── go_convergence.tsv
│   ├── top_candidate_genes.tsv
│   └── cross_inversion_summary.tsv
└── figures/
    └── <cid>__profile.pdf
```

When the registry patch is applied, the per_arrangement_burden / config_spectrum
files instead live under `${RESULTS_REGISTRY_DIR}/interval/<chrom>/...` and
`${RESULTS_REGISTRY_DIR}/pairwise/<stat>/<chrom>/...`, with manifest rows
linking back to candidate_id and group_id.

## Testing recommendation

For initial validation, pick a single high-confidence candidate (e.g. one
from your SAFR set with clean three-band PCA), run:

```bash
sbatch STEP_60_build_gene_intervals.sh
# wait for completion
sbatch compute/STEP_C60_cargo_inventory.slurm --candidate-id <CID>
# wait, inspect ${CARGO_INVENTORY_DIR}/diagnostic_table.tsv
sbatch --array=1-1 launchers/run_cargo_all.slurm   # just one task
```

This exercises every code path on one candidate before unleashing the array
across all 28 chromosomes.

## Things to think about for the manuscript

The cargo module produces per-inversion narratives at three depths:

1. **Inventory + diagnostics** — a one-sentence-per-inversion description:
   "candidate LG12_17 spans 1.8 Mb on LG12, contains 47 protein-coding genes
   at density 26/Mb, 23 missense and 3 LoF segregating sites; 22 HOM_REF and
   18 HOM_INV homozygotes available for arrangement comparison."

2. **Per-arrangement signatures** — the beetle-paper analog upgraded:
   "Of the 23 segregating missense sites, 14 segregate at higher frequency on
   HOM_INV (mean Δfreq = +0.21), with 4 sites at FST > 0.5. Per-gene VESM
   burden is asymmetric for 6 genes (perm p < 0.05); the top-ranked is
   GENE_X with Δburden = +12.4 LLR units carried preferentially on HOM_INV."

3. **Configuration spectrum** — the beetle paper does NOT do this:
   "Within the HOM_INV arrangement, 3 genes show pairwise missense MI Z >
   3 indicating coupled (haplotype-resolved) variants; on HOM_REF, only 1
   gene reaches that threshold, consistent with the inversion locking
   together evolved haplotype configurations."

4. **Convergence** — synthesis across inversions:
   "Family X is captured by 6 of the 14 analyzable inversions (null 95th
   percentile = 2), suggesting recurrent biological selection rather than
   chance gene-density coincidence."

The synthesis tables (cross_inversion/) feed directly into a single
multi-panel manuscript figure showing convergence + per-inversion top genes.
