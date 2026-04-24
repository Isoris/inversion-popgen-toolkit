# inversion-popgen-toolkit

Analytical pipelines for low-coverage whole-genome population genomics in
aquaculture species, developed for haplotype-resolved North African
catfish (*Clarias gariepinus*) but designed to be broadly applicable.

## What this repo does

End-to-end detection and characterization of chromosomal inversions in
a 226-sample pure *C. gariepinus* broodstock cohort sampled from a
commercial Thai hatchery, at ~9× Illumina short-read WGS coverage,
mapped against the Gar subgenome reference (`fClaHyb_Gar_LG.fa`, 28
pseudochromosomes, ~964 Mb) extracted from the haplotype-resolved F₁
hybrid assembly (§1 of the companion manuscript). Pipeline outputs
feed the `MS_Inversions_North_african_catfish` manuscript.

## Repo layout — one-paragraph version

Two parallel module trees. `Modules/` (MODULE_1..4) handles
preprocessing + SV-caller outputs: read QC, SNP discovery, population
structure, heterozygosity, Clair3 SNP/INDEL, DELLY, Manta, conservation
annotation. `inversion_modules/` (phase_1..7) is the inversion-discovery
pipeline itself: input prep → phase-2 discovery (local PCA, MDS, GHSL,
staircase detector, landscape blocks) → phase-3 bp-resolution
refinement → phase-4 classification (5 axes) → phase-5 per-candidate
deep analysis → phase-6 secondary confirmation (LD, Fst, Hobs) →
phase-7 gene content. For a full map with manuscript-section
mapping, see **[`docs/MODULE_MAP.md`](docs/MODULE_MAP.md)**.

## Key docs

| File | What it is |
|---|---|
| [`docs/MODULE_MAP.md`](docs/MODULE_MAP.md) | Map of both module trees → manuscript sections |
| [`docs/OPEN_TODOS.md`](docs/OPEN_TODOS.md) | Live work queue (laptop-doable / HPC-blocked / backlog) |
| [`docs/ENGINES.md`](docs/ENGINES.md) | Engine inventory (instant_q, GHSL, popstats, …) |
| [`docs/toolkit_audit.md`](docs/toolkit_audit.md) | Structural reorg plan + known gaps |
| [`docs/kde_mode_estimation_review.md`](docs/kde_mode_estimation_review.md) | Methodological note on KDE mode estimation |
| `inversion_modules/README.md` | Inversion-discovery pipeline root |
| `inversion_modules/phase_5_qc_triage/README.md` | QC diagnosis substrate |

Every module / phase has its own README. Sub-blocks under
`phase_2_discovery/` (2a–2e) each have READMEs, as does each of the
twelve top-level phases.

## Dependencies

- Bash 4+, SLURM on LANTA
- Python 3.10+ (pandas, pysam, numpy, scipy)
- R 4.2+ (data.table, RSpectra, jsonlite, ggplot2, ComplexHeatmap, modeest)
- ANGSD (fork: `github.com/Isoris/angsd_fixed_HWE`), BEAGLE, PCAngsd,
  NGSadmix, evalAdmix, ngsRelate, ngsF-HMM
- DELLY2, Manta, Clair3, WhatsHap
- bcftools, samtools, mosdepth, SIFT4G, VESM, bcftools csq
- MASH (for species-identity verification vs both parental subgenomes)

## Running on LANTA

Expected project root:
`/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`
with this repo cloned inside. SLURM account `lt200308`, conda env
`assembly`. Install/update helper:
`tools/install_scripts/install_update_2026-04-20.sh`.

## Dev utilities

Under `tools/`:
- `_schema_check.py` — validate the three BK JSON schemas
- `_code_field_check.py` — check each `keys_extracted` path maps to a source field
- `_bk_rename.py` — apply the BK key rename pass
- `_rcheck.py` — lightweight R syntax balance checker
- `build_inventory.py` — repo file inventory
- `standardize_batch.py` — batch header/echo standardization helper

Run all from the repo root: `python3 tools/_schema_check.py`, etc.

## R packages (vendored)

Under `R_packages/`:
- `figrid/` — panel-aware figure composition
- `exactplot/` — plotting helpers

## License

MIT — see `LICENSE`.

## Manuscript

This repo accompanies the manuscript
`MS_Inversions_North_african_catfish` (in preparation, target:
Nature Communications / Genome Research tier). A separate companion
manuscript describes the haplotype-resolved F₁ hybrid assembly
(*C. gariepinus* ♀ × *C. macrocephalus* ♂) used as the reference
here.
