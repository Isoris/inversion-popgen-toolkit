# MODULE_5A3 — Discovery Postprocessing

## Purpose

Downstream annotation and summary after the core discovery engine. Overlaps candidate regions with population-level statistics, runs regional PCA with adaptive-k grouping, and produces combined diagnostic panels.

## Scope

STEP11 through STEP13 + supplementary helpers:

| Step | Script | What it does |
|------|--------|-------------|
| 11 | `STEP11_overlap_candidate_regions_with_theta_and_het.R` | Overlap candidates with population theta (tP) |
| 12 | `STEP12_candidate_region_pca_groups_and_plotC.R` | Regional PCA + BIC adaptive-k grouping (k=2..5) |
| 13 | `STEP13_make_combined_panels_ABC.R` | Combined diagnostic panels A+B+C |

### Supplementary helpers

| Script | What it does |
|--------|-------------|
| `STEP_SUPP_inversion_candidate_counts.sh` | Count candidates per chromosome |
| `STEP_SUPP_inversion_candidate_master_table.py` | Build master candidate summary table |
| `STEP_SUPP_inversion_notes.sh` | Diagnostic notes |
| `STEP_SUPP_inversion_step_manifest.sh` | Step file manifest |
| `STEP_SUPP_theta_overlap_manifest.sh` | Theta overlap manifest |

## Main Inputs

- Candidate regions table from 5A2 (`candidate_regions.tsv.gz`)
- Per-chr dosage from 5A2
- Population theta windows from STEP10b
- Sample metadata

## Main Outputs

- `07_het_overlap/` — theta-overlapped candidate annotations
- `08_regional_pca/` — per-candidate regional PCA groups + plot C data
- `09_combined_plots/` — combined diagnostic panel figures

## Runner

```bash
bash runners/run_5A3_discovery_postprocessing.sh           # all STEP11–13
bash runners/run_5A3_discovery_postprocessing.sh --step 12 # STEP12 only
```

## Dependencies

- Requires 5A2 outputs (candidate table, dosage, theta bridge)
- R: data.table, ggplot2, mclust
- `00_inversion_config.sh`

## Downstream

Feeds into MODULE_5B (followup), MODULE_5C (LD), MODULE_5D (FST), MODULE_5E (HOBS).
