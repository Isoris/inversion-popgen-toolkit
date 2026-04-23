# MODULE_6_Cargo — Inversion gene content & per-arrangement evolution

Final analytical module for the inversion manuscript. Consumes results from
all earlier phases (karyotypes from phase 4, candidate intervals from phase 2,
deleterious annotations from MODULE_CONSERVATION) and resolves them into a
per-inversion biological narrative: **what genes does each inversion contain,
and how do those genes evolve differently between arrangements?**

This module is a *capstone integration* layer. It does **no new variant calling
and no new annotation**. Every score reuses an annotation already produced by
MODULE_CONSERVATION (SnpEff impact, csq haplotype-aware effect, VESM LLR, SIFT)
and re-stratifies it by karyotype.

## What the module produces

For every inversion candidate that passes the diagnostic filter:

**Inventory (6A)** — gene list with density, repeat content, breakpoint
distance, gene-family / paralog assignment per gene.

**Functional enrichment (6B)** — GO / KEGG / family enrichment of inversion
gene set vs three matched backgrounds (genome-wide, matched collinear regions,
within-locus across arrangements).

**Per-arrangement evolutionary signatures (6C)** — for each gene inside each
candidate, separately on `inv_<cid>_HOM_REF` and `inv_<cid>_HOM_INV`
chromosomes:

  - `cargo_burden` — sum and mean VESM LLR per gene per arrangement
  - `cargo_lof_count` — count of HIGH-impact (stop_gained, frameshift,
    splice_acceptor/donor, start_lost) variants per gene per arrangement
  - `cargo_nonsyn_synon` — segregating NS and SS counts per arrangement
    (the beetle paper's Table S4 split by arrangement, the upgrade)
  - `cargo_sfs` — joint missense SFS as a 2D histogram across arrangements
  - `cargo_fst_per_site` — Hudson FST per missense site between arrangements
    (the beetle paper's Table S5 as a continuous distribution, not bins)
  - `cargo_config_mi` — pairwise mutual information between missense sites
    inside each arrangement, with Z-score against a frequency-preserving
    permutation null (Level-2 configuration spectrum)

**Cross-inversion synthesis (6E)** — convergence analysis: do multiple
independent inversions repeatedly capture the same gene families? Per-inversion
two-page profile figures.

## What this module deliberately does NOT do

- **No CNV** — DELLY's small-DUP floor (~300 bp practical) misses most
  single-exon duplications. Either skip CNV (current default) or repurpose
  the PAV mosdepth machinery as a future addition.
- **No MK / dN/dS** — requires outgroup codon alignment via the Fish10K
  Cactus output. Deferred to the *C. macrocephalus* manuscript where the
  outgroup design is natural.
- **No re-annotation** — VESM, SIFT, SnpEff, csq are run once per variant
  in MODULE_CONSERVATION. This module only re-aggregates by karyotype subset.

## Directory layout

```
phase_7_cargo/
├── README.md                              # this file
├── 00_cargo_config.sh                     # paths + cargo-specific settings
├── STEP_60_build_gene_intervals.sh        # one-time: GFF3 → gene_intervals.bed.gz
├── STEP_61_eggnog_annotate.sh             # one-time: REF_PEP → eggnog_annotations.tsv
├── compute/
│   ├── STEP_C60_cargo_inventory.py        # 6A: inventory + diagnostics
│   ├── STEP_C61_per_arrangement_burden.py # 6C Level 1: load, NS/SS, LoF, FST, SFS
│   └── STEP_C62_config_spectrum.py        # 6C Level 2: pairwise MI + null
├── analysis/
│   ├── STEP_C63_functional_enrichment.R   # 6B: enrichment vs 3 backgrounds
│   └── STEP_C64_cross_inversion_synthesis.R # 6E: convergence
├── plot/
│   └── STEP_C65_per_inversion_profile.R   # 6E: per-inversion profile figure
└── launchers/
    └── run_cargo_all.slurm                # SLURM array driver, one task per candidate
```

## Inputs (all already exist after upstream pipelines run)

| Source | Path | Used by |
|---|---|---|
| Candidate intervals | `${SNAKE_CAND_FILE}` | 6A, 6C, 6E |
| Sample karyotypes per candidate | `${INVDIR}/<cand_dir>/data/sample_karyotypes.tsv` | 6A diagnostics |
| Karyotype sample groups | `sample_registry`: `inv_<cid>_HOM_REF/HET/HOM_INV` | 6C, 6E |
| Variant master (annotated) | `${MODCONS}/16_merged_variant_tables/variant_master_scored.tsv` | 6C |
| VESM LLR per variant | `${MODCONS}/05B_vesm/vesm_variant_scores.tsv` | 6C burden |
| Per-chrom cohort VCFs | `${MODCONS}/03_variants/normalized/${CHR}.clair3.norm.vcf.gz` | 6C genotypes |
| Reference GFF3 | `${REF_GFF3}` | STEP_60 |
| Reference FASTA index | `${REF_FAI}` | STEP_60 |
| Reference protein FASTA | `${REF_PEP}` | STEP_61 |
| RepeatMasker BED (if available) | `${REPEAT_BED}` | 6A (optional) |

## Outputs (all wired through `results_registry`)

Per-arrangement per-gene results land as `interval_summary` rows under the
`inv_<cid>_HOM_REF` and `inv_<cid>_HOM_INV` group_ids. A-vs-B comparisons
land as `pairwise` rows. New stat values introduced by this module:

| stat | kind | who | Description |
|---|---|---|---|
| `cargo_burden` | interval_summary | one arrangement group | per-gene VESM sum/mean |
| `cargo_lof_count` | interval_summary | one arrangement group | per-gene HIGH-impact count |
| `cargo_nonsyn_synon` | interval_summary | one arrangement group | per-gene NS/SS counts |
| `cargo_sfs` | interval_summary | one arrangement group | per-gene missense SFS |
| `cargo_config_mi` | interval_summary | one arrangement group | per-gene MI Z-score |
| `cargo_fst_per_site` | pairwise | HOM_REF vs HOM_INV | per-site FST table |
| `cargo_burden_diff` | pairwise | HOM_REF vs HOM_INV | per-gene Δburden + permutation p |

These names are added to the stat enum in `registry_loader.py` (see patch
`registries/patches/STAT_ENUM_cargo_v1.patch` shipped with this module).

## Diagnostic gate

A candidate is analyzable at each level only if the homozygote counts and
missense density are sufficient. STEP_C60 produces a diagnostic table:

| Column | Meaning |
|---|---|
| `n_HOM_REF`, `n_HOM_INV`, `n_HET` | sample counts per karyotype |
| `n_genes_inside` | genes overlapping the candidate interval |
| `median_NS_per_gene`, `max_NS_per_gene` | missense density |
| `level_1_ok` | `n_HOM_REF >= 5 AND n_HOM_INV >= 5` (descriptive per-arrangement) |
| `level_2_ok` | `n_HOM_REF >= 10 AND n_HOM_INV >= 10 AND max_NS_per_gene >= 5` (configuration spectrum) |

Candidates failing `level_1_ok` get inventory only. Candidates passing
`level_1_ok` but failing `level_2_ok` get burden + NS/SS + FST + SFS.
Candidates passing `level_2_ok` additionally get the configuration spectrum.

## Run order

```bash
# One-time setup (idempotent — skips if outputs exist)
sbatch STEP_60_build_gene_intervals.sh
sbatch STEP_61_eggnog_annotate.sh   # optional, enables 6B family enrichment

# Per-candidate compute (SLURM array)
sbatch launchers/run_cargo_all.slurm

# Cross-candidate analysis (single node)
sbatch analysis/STEP_C63_functional_enrichment.slurm
sbatch analysis/STEP_C64_cross_inversion_synthesis.slurm

# Per-candidate figures
sbatch plot/STEP_C65_per_inversion_profile.slurm
```
