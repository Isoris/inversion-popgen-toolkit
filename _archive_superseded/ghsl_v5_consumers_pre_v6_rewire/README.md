# ghsl_v5_consumers_pre_v6_rewire

Pre-patch copies of files that were rewired from GHSL v5 to GHSL v6
in chat 14 (2026-04-18).

## Why this archive exists

Chat 13 swapped `phase_2_discovery/2e_ghsl/` from v5 to v6 but did
not update the downstream consumers. Those still read `<chr>.ghsl_v5.annot.rds`
and `<chr>.ghsl_v5.karyotypes.rds` and expected v5 column names
(`ghsl_v5_score`, `ghsl_rank_stability`, `ghsl_div_contrast_z`,
`ghsl_div_bimodal`). v6 did not emit those paths, and its column
names differ. Net effect: phase-4b Tier-3 karyotype confirmation was
silently returning UNINFORMATIVE for every sample. Flagged in the
chat-13 2e_ghsl/README.md as the wiring gap for chat 14/15.

## Files in this directory

| File | Superseded reason |
|------|-------------------|
| `lib_ghsl_confirmation.R` | v5 path/column reads. Chat 14 rewires to v6 and adds panel-based APIs. |
| `run_all.R` | v5 path/column reads at lines 453–525. Chat 14 rewires to v6. |
| `STEP_C01d_candidate_scoring_wired_25_v934_registry.R` | Reads `iv$ghsl_v5_score_max`. Chat 14 renames to v6. |
| `STEP_D05_ghsl_stability.R` | Header doc referenced `STEP_C04_snake3_ghsl_v5.R`. Chat 14 updates doc (code is a stub, not used). |
| `STEP_C04b_snake3_ghsl_classify.R` | v6 classifier pre-patch. Chat 14 adds per-chrom annot / karyotype / per-sample-panel RDS emit. |
| `STEP_C04_snake3_ghsl_v6.R` | v6 heavy engine pre-patch. Chat 14 adds chunking fix (single-chrom runs load only the active chrom) and finer default scales. |
| `registry_loader.R` | Chat 14 extends `load_compute_api()` with `ghsl_*` methods so registries can query GHSL panel data on demand (same pattern as `pairwise_stat`). |

## Rollback

If chat 14's changes need to be reverted, these files are the pre-change
state. None should be restored to the live tree as-is unless you also
revert the paired changes in the classifier (the new per-chrom RDS emit)
and the registry_loader compute API extension.

Archived: 2026-04-18 (chat 14).
