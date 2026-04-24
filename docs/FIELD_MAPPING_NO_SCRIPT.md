# Field-mapping analysis — NO_SCRIPT_DECLARED schemas

For the 6 schemas where `source_script` points at missing files, or the declared writer doesn't exist. I did a broader filesystem search to find any candidate writer that could plausibly produce the schema's keys.

Outcomes fall into 3 categories:

- **Drift**: A live script exists but its name/path differs from `source_script` (e.g., `C01n` → `STEP_D09n`). Fix is to update the schema's `source_script` and/or the script's output field names.
- **Candidate exists but in unexpected phase**: e.g., burden writers live in phase_12_cargo instead of phase_7 validation. Possibly a deliberate architecture choice — confirm with Quentin.
- **Schema may be vestigial**: e.g., `inv_detect_blocks`, `synteny_dollo` — no candidate writer anywhere, and another schema (e.g., `synteny_v6`) may supersede it.

## `burden`

- **Candidate writer(s) found by broader search:**
  - `Modules/MODULE_CONSERVATION_CORE/STEP_15_burden_tables.sh`
  - `inversion_modules/phase_12_cargo/compute/STEP_C61_per_arrangement_burden.py`
  - `inversion_modules/phase_9_classification/figures/fig4e_cumulative_burden_per_group.R`
- **Note:** Schema source_script mentions `STEP_C01f_c_burden_regression.R` which doesn't exist. The closest live writers compute burden differently: STEP_15 is a bash orchestrator, STEP_C61 is Python phase_12 burden, fig4e is a figure script. Need to decide which one OR build the declared C01f_c writer. Also phase_12_cargo territory — check with Quentin before wiring.

### Candidate writer: `Modules/MODULE_CONSERVATION_CORE/STEP_15_burden_tables.sh`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_burden_beta` | `beta_hom_inv_vs_hom_ref` | **MISSING** | — |
| `q6_burden_p` | `p_value` | **MISSING** | — |
| `q6_burden_ci_lower` | `ci_lower` | **MISSING** | — |
| `q6_burden_ci_upper` | `ci_upper` | **MISSING** | — |
| `q6_n_del_sites` | `n_del_sites_tested` | **MISSING** | — |
| `q6_selection_signal` | `selection_signal` | **MISSING** | — |

### Candidate writer: `inversion_modules/phase_12_cargo/compute/STEP_C61_per_arrangement_burden.py`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_burden_beta` | `beta_hom_inv_vs_hom_ref` | **MISSING** | — |
| `q6_burden_p` | `p_value` | **MISSING** | — |
| `q6_burden_ci_lower` | `ci_lower` | **MISSING** | — |
| `q6_burden_ci_upper` | `ci_upper` | **MISSING** | — |
| `q6_n_del_sites` | `n_del_sites_tested` | **MISSING** | — |
| `q6_selection_signal` | `selection_signal` | **MISSING** | — |

### Candidate writer: `inversion_modules/phase_9_classification/figures/fig4e_cumulative_burden_per_group.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q6_burden_beta` | `beta_hom_inv_vs_hom_ref` | **MISSING** | — |
| `q6_burden_p` | `p_value` | **MISSING** | — |
| `q6_burden_ci_lower` | `ci_lower` | **MISSING** | — |
| `q6_burden_ci_upper` | `ci_upper` | **MISSING** | — |
| `q6_n_del_sites` | `n_del_sites_tested` | **MISSING** | — |
| `q6_selection_signal` | `selection_signal` | **MISSING** | — |

## `carrier_reconciliation`

- **Candidate writer(s) found by broader search:**
  - `inversion_modules/phase_12_cargo/extra_plots/tables/TABLE_06_top_inversion_carriers.R`
  - `inversion_modules/phase_12_cargo/extra_plots/compute/_lib_group_carriership.R`
- **Note:** Schema says `audit_cores_vs_sv (new)` + `carrier_audit.R` — neither exists. phase_12 has carrier-related scripts but they're downstream (visualization / breeding implications). The schema covers Q3+Q6+Q7 reconciliation of carrier lists — probably needs a fresh writer in phase_7 or phase_9, not phase_12.

### Candidate writer: `inversion_modules/phase_12_cargo/extra_plots/tables/TABLE_06_top_inversion_carriers.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_n_reconciliation_sources` | `sources_used` | **MISSING** | — |
| `q7_reconciliation_agreement` | `agreement_fraction` | **MISSING** | — |
| `q7_n_disagreeing_samples` | `n_disagreeing` | **MISSING** | — |

### Candidate writer: `inversion_modules/phase_12_cargo/extra_plots/compute/_lib_group_carriership.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_n_reconciliation_sources` | `sources_used` | **MISSING** | — |
| `q7_reconciliation_agreement` | `agreement_fraction` | **MISSING** | — |
| `q7_n_disagreeing_samples` | `n_disagreeing` | **MISSING** | — |

## `gene_conversion_tracts`

- **Candidate writer(s) found by broader search:**
  - `inversion_modules/phase_7_karyotype_groups/proposal/gene_conversion_detector.R`
- **Note:** Schema source_script is literally blank (None). Writer exists — `gene_conversion_detector.R` in phase_7 proposal. 487 lines, no registry calls, no file outputs, no CLI entry. Probably a function library waiting on a runner. Status mirrors the Q4 mechanism cheats (SCRIPT_IS_LIB pattern) — maybe Quentin is managing this one too?

### Candidate writer: `inversion_modules/phase_7_karyotype_groups/proposal/gene_conversion_detector.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q2_gc_total_tracts` | `total_tracts` | **MATCH** | L440 |
| `q2_gc_total_samples_with_gc` | `total_samples_with_gc` | **MATCH** | L441 |
| `q2_gc_n_snps_pass_qc_diagnostic` | `snp_qc.n_snps_pass_qc_diagnostic` | **DRIFT** | `max_flagged_snps`, `n_dropped_non_diagnostic`, `n_snps_pass_qc`, `n_snps_pass_qc_diagnostic`, `n_snps_total` |

## `inv_detect_blocks`

- **Candidate writers:** *none found in live tree*
- **Note:** Schema source_script is `inv_detect_v9.3 + 12_bridge_to_codebase.R`. Neither file exists in live tree. Searched entire tree including _archive — no match. This is likely a vestigial schema from an older architecture where `inv_detect_v9.3` was a standalone pipeline. Check if it still has a role post-pass-15 rename.

Declared keys (no candidate writer to diff against):

- `q1_staircase_score` (from `staircase_score`)
- `q1_staircase_rank` (from `staircase_rank`)
- `q1_n_scales_detected` (from `n_scales_detected`)
- `q1_consensus_strength` (from `consensus_strength`)
- `q3_regime_change_left` (from `regime_change_left_bp`)
- `q3_regime_change_right` (from `regime_change_right_bp`)

## `peel_diagnostic`

- **Candidate writer(s) found by broader search:**
  - `inversion_modules/phase_2_discovery/2d_candidate_detection/STEP_D09n_peeling_diagnostic.R`
- **Note:** Schema source_script says `C01n peel / inv_detect Phase 9` — this is a renamed-script case. Live script is `STEP_D09n_peeling_diagnostic.R` (naming drift: C01n → D09n). Script computes `l1b_peel` results but not the exact fields the schema expects (`l1b_score`, `l2_score`, `peel_combined_score`, `peel_survives`). Writer exists but output structure drifted.

### Candidate writer: `inversion_modules/phase_2_discovery/2d_candidate_detection/STEP_D09n_peeling_diagnostic.R`

| Key (schema) | `from` field | Status | Similar names in writer |
|---|---|---|---|
| `q7_peel_l1b_score` | `l1b_score` | **MISSING** | — |
| `q7_peel_l2_score` | `l2_score` | **MISSING** | — |
| `q7_peel_combined` | `peel_combined_score` | **DRIFT** | `max_peel_frac`, `peel`, `peel_mode`, `peel_samples`, `v_i_peeled` |
| `q7_peel_survives` | `peel_survives` | **DRIFT** | `max_peel_frac`, `peel`, `peel_mode`, `peel_samples`, `v_i_peeled` |

## `synteny_dollo`

- **Candidate writers:** *none found in live tree*
- **Note:** Schema source_script is `synteny_inversion_detection.sh + dollo_parsimony_inversions.py + mashmap`. None of those files exist in the live tree. `synteny_v6.schema.json` (OK_DIRECT) likely supersedes this. `synteny_dollo` may be archiveable.

Declared keys (no candidate writer to diff against):

- `q5_conservation_class` (from `conservation_class`)
- `q5_present_c_gariepinus` (from `present_in_c_gariepinus`)
- `q5_present_c_macrocephalus` (from `present_in_c_macrocephalus`)
- `q5_outgroup_state` (from `outgroup_state`)
- `q5_dollo_state` (from `dollo_inferred_state`)
