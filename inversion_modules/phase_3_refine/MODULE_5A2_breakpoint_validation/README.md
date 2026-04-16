# MODULE_5A2: Breakpoint Validation & Seeded Deconvolution

## What this module does

Takes DELLY2 inversion calls and validates them against per-sample
read-level breakpoint evidence (discordant pairs, split reads, soft clips),
cross-references with Snake 2 band assignments (REF/HET/INV), runs
Fisher's exact + Cochran-Armitage trend tests, and produces
**breakpoint-validated seed sets** for C01i seeded deconvolution.

## Pipeline

```
STEP 01 — Extract DELLY2 INV calls, match to Snake candidates
STEP 02 — Per-sample BAM evidence extraction (pysam)
           + group assignment from Snake2 or local PCA
STEP 03 — Fisher/CA/chi-square tests + concordance + seed generation
STEP 04 — Publication-quality validation plots
```

## Key outputs

| File | Description |
|------|-------------|
| `01_per_sample_evidence/{inv_id}_evidence.tsv` | Per-sample support/no-support + evidence breakdown |
| `03_statistical_tests/all_candidates_tests.tsv` | All test results across candidates |
| `03_statistical_tests/manuscript_sentences.txt` | Copy-paste text for Methods/Results |
| `04_deconvolution_seeds/{inv_id}_seeds.tsv` | Confirmed REF/HET/INV sample seeds for C01i |
| `04_deconvolution_seeds/seed_summary.tsv` | Which candidates qualify for seeding |
| `05_plots/{inv_id}_validation.pdf` | Per-candidate figure |
| `05_plots/genome_wide_validation_summary.pdf` | Forest plot + trend across all candidates |

## Dependencies

- Python 3 with pysam, numpy, matplotlib
- DELLY2 INV catalog (from delly_sv_INV/)
- Markdup BAMs (from delly_sv/00_markdup/)
- Snake 2 band assignments (optional, from 06_mds_candidates/snake2_community/)
- Dosage files (optional, for local PCA fallback)

## How seeds feed into C01i

The seed files contain confirmed sample-to-karyotype assignments:

```
sample    seed_class       evidence_score
CGA045    INV_confirmed    12
CGA078    HET_confirmed    6
CGA102    REF_confirmed    0
```

C01i can use these as initialization for marker co-segregation clustering,
replacing the unsupervised discovery with prior-informed decomposition.

## Usage

```bash
cd /scratch/.../inversion_codebase_v8.5/MODULE_5A2_breakpoint_validation
bash run_breakpoint_validation.sh
```
