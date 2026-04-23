# =============================================================================
# STRUCTURED RESULT BLOCK SCHEMAS (Tier 2 Storage)
# =============================================================================
#
# 18 block types, one per major analysis.
# Each block is a JSON file in per_candidate/<cid>/structured/
#
# Every block has a COMMON HEADER:
#   {
#     "block_type": "fisher_test",
#     "candidate_id": "LG12_inv17",
#     "chrom": "C_gar_LG12",
#     "region": [8204500, 12510000],
#     "source_script": "STEP03_statistical_tests_and_seeds.py",
#     "source_version": "v2",
#     "computed_at": "2026-04-16T14:32:00",
#     "groups_used": {
#       "source": "C01i_v2",
#       "version": "2026-04-15",
#       "validation_level": "VALIDATED",
#       "counts": {"HOM_STD": 89, "HET": 98, "HOM_INV": 32, "REC": 7}
#     },
#     "keys_produced": ["q7_layer_d_fisher_or", "q7_layer_d_fisher_p", ...],
#     "raw_files": ["STEP03_contingency_table.tsv"],
#     ... block-specific fields ...
#   }
#
# In R, write with: jsonlite::write_json(block, path, auto_unbox=TRUE, pretty=TRUE)
# In R, read with:  jsonlite::read_json(path)
# In Python:        json.dump(block, f, indent=2) / json.load(f)
# On terminal:      cat structured/fisher_test.json  (human-readable)
#
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 1: existence_layer_a.json (Local PCA evidence)
# Source: C01a precomp + C01b cores + C01d scoring
# Question: Q1, Q7
# Group requirement: NONE
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "existence_layer_a",
#   "detected": true,
#   "inv_likeness_mean": 0.7435,
#   "inv_likeness_max": 0.9213,
#   "inv_likeness_sd": 0.0821,
#   "family_likeness_mean": 0.1203,
#   "beta_pval": 1.2e-09,
#   "beta_qval": 8.9e-08,
#   "core_family": "1L",
#   "pa_pattern": "SML",
#   "n_windows": 86,
#   "shape_class": "strong_square",
#   "squareness": 0.8912,
#   "nn_birth": 12,
#   "nn_survives_40": true,
#   "nn_survives_80": true,
#   "dimensions": {
#     "D1": 0.842, "D2": 0.780, "D3": 0.710, "D4": 0.623,
#     "D5": 0.801, "D6": 0.750, "D7": 0.890, "D8": 0.920,
#     "D9": 0.683, "D10": 0.712, "D11": 0.850, "D12": 0.780
#   },
#   "composite_score": 0.798,
#   "dim_positive": 11,
#   "keys_produced": ["q7_layer_a_detected", "q7_layer_a_inv_likeness",
#                      "q7_layer_a_beta_pval", "q7_layer_a_beta_qval",
#                      "q7_layer_a_core_family", "q7_layer_a_pa_pattern",
#                      "q1_d01..q1_d12", "q1_composite_score", "q1_dim_positive"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 2: existence_layer_b.json (SV caller evidence)
# Source: C00 flashlight + MODULE_4D/4E/4G
# Question: Q7
# Group requirement: NONE
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "existence_layer_b",
#   "detected": true,
#   "delly_inv": true,
#   "manta_inv": true,
#   "bnd_triangulated": false,
#   "delly_id": "DEL00012345",
#   "manta_id": "MantaINV:0:1:2:0",
#   "bp1_pos": 8204832,
#   "bp2_pos": 12509700,
#   "bp1_cipos": [-300, 300],
#   "bp2_ciend": [-350, 350],
#   "pe_support": 23,
#   "sr_support": 12,
#   "n_carriers_sv": 103,
#   "sv_audit": {
#     "delly_raw_carriers": 137,
#     "delly_strict_carriers": 103,
#     "delly_carrier_loss": 34,
#     "delly_site_qual": 1250,
#     "delly_site_pe": 23,
#     "delly_passes_strict": true,
#     "manta_raw_carriers": 128,
#     "manta_pass_carriers": 118,
#     "manta_carrier_loss": 10,
#     "expected_dropout_pct": 34.2,
#     "observed_dropout_pct": 20.8,
#     "pop_prior_applied": false,
#     "pop_prior_n_rescued": 0,
#     "bnd_3to3_present": true,
#     "bnd_5to5_present": true,
#     "bnd_rescued": false
#   },
#   "keys_produced": ["q7_layer_b_detected", "q7_layer_b_delly", "q7_layer_b_manta",
#                      "q7_layer_b_n_carriers", "q7_layer_b_pe_support",
#                      "q7_layer_b_sr_support", "q7_layer_b_cipos", "q7_layer_b_ciend",
#                      "q7b_delly_raw_carriers", "q7b_delly_strict_carriers", ...]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 3: existence_layer_c.json (GHSL haplotype contrast)
# Source: C04 GHSL (Layer C — haplotype contrast)
# Question: Q7
# Group requirement: NONE (GHSL produces its own partition)
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "existence_layer_c",
#   "detected": true,
#   "ghsl_version": "v6",
#   "ghsl_contrast_mean": 0.142,
#   "ghsl_n_pass": 72,
#   "ghsl_n_total_windows": 86,
#   "ghsl_pct_pass": 0.837,
#   "ghsl_quality": "HIGH",
#   "partition_stable": true,
#   "partition_entropy": 0.34,
#   "partition_groups": 3,
#   "pca_ghsl_concordance": 0.91,
#   "keys_produced": ["q7_layer_c_ghsl_detected", "q7_layer_c_ghsl_contrast",
#                      "q7_layer_c_ghsl_n_pass", "q7_layer_c_ghsl_pct_pass",
#                      "q7_layer_c_ghsl_quality", "q7_layer_c_partition_stable"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 4: existence_layer_d.json (Genotype-breakpoint association)
# Source: STEP03 statistical tests
# Question: Q7
# Group requirement: uses groups but IS the group validation
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "existence_layer_d",
#   "tested": true,
#   "contingency_table": [[45, 8], [12, 161]],
#   "table_labels": {"rows": ["carrier", "non_carrier"],
#                     "cols": ["support_yes", "support_no"]},
#   "fisher_or": 23.7,
#   "fisher_ci": [8.9, 62.1],
#   "fisher_p": 4.98e-06,
#   "armitage_z": 5.12,
#   "armitage_p": 3.01e-07,
#   "concordance": 0.89,
#   "concordance_description": "89% of PCA-INV have breakpoint support, 7% of PCA-REF do",
#   "inv_support_frac": 0.849,
#   "het_support_frac": 0.612,
#   "ref_support_frac": 0.069,
#   "samples_inv_with_support": ["CGA009", "CGA045", "..."],
#   "samples_inv_no_support": ["CGA078", "CGA201"],
#   "seed_qualified": true,
#   "seed_reason": "Fisher p < 0.05, inv_support 85%, concordance 89%, n_seeds 45",
#   "keys_produced": ["q7_layer_d_tested", "q7_layer_d_fisher_or",
#                      "q7_layer_d_fisher_p", "q7_layer_d_armitage_z",
#                      "q7_layer_d_armitage_p", "q7_layer_d_concordance"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 5: boundary_left.json (Left breakpoint characterization)
# Source: C01g boundary catalog
# Question: Q3
# Group requirement: NONE
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "boundary_left",
#   "position_bp": 8204832,
#   "position_sources": {
#     "pca_boundary": 8200000,
#     "sv_breakpoint": 8204832,
#     "staircase_vote": 8200000,
#     "delta12_changepoint": 8210000,
#     "fst_step": 8205000
#   },
#   "concordance_kb": 10.0,
#   "hardness": 0.85,
#   "sharpness": 0.12,
#   "type": "STRUCTURAL_SHARP",
#   "verdict": "confirmed_structural",
#   "cheats_supporting": {
#     "cheat4_fst_step": {"value": 0.12, "significant": true},
#     "cheat10_depth_ratio": {"value": 1.12, "anomalous": false},
#     "cheat11_clip_count": {"value": 23, "significant": true},
#     "cheat17_fossil": false,
#     "cheat21_te_density": {"value": 0.34, "enrichment": 1.8, "status": "ENRICHED"}
#   },
#   "n_cheats_supporting": 5,
#   "fst_profile": {
#     "distances_kb": [-200, -100, -50, 0, 50, 100, 200],
#     "values": [0.02, 0.03, 0.05, 0.12, 0.14, 0.13, 0.12],
#     "step_magnitude": 0.09
#   },
#   "hobs_profile": {
#     "distances_kb": [-100, 0, 100],
#     "values": [0.31, 0.38, 0.35],
#     "step_magnitude": 0.04
#   },
#   "repeat_context": {
#     "density_10kb": 0.34,
#     "dominant_class": "LINE",
#     "gc_content": 0.41,
#     "enrichment_vs_genome": 1.8
#   },
#   "mask_offset": {
#     "pca_to_sv_distance_bp": 4832,
#     "bp_in_repeat": true,
#     "repeat_class_at_bp": "LINE/L2"
#   },
#   "keys_produced": ["q3_left_bp", "q3_left_sv_bp", "q3_left_hardness",
#                      "q3_left_sharpness", "q3_left_type", "q3_left_verdict",
#                      "q3_left_n_cheats", "q3_left_clip_count",
#                      "q3_left_depth_anomaly", "q3_left_is_fossil",
#                      "q3_left_fst_*", "q3_left_hobs_*", "q3_left_repeat_*"]
# }
# NOTE: boundary_right.json has identical schema with "right" prefix


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 6: carrier_reconciliation.json
# Source: audit_cores_vs_sv.R + population_regenotype.py
# Question: Q3 (carrier keys), Q6, Q7
# Group requirement: N/A (this IS about the groups)
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "carrier_reconciliation",
#   "pca_carriers": {
#     "n": 130,
#     "source": "C01i_decomposition",
#     "includes_het": true,
#     "includes_hom_inv": true,
#     "sample_ids": ["CGA009", "CGA010", "..."]
#   },
#   "sv_carriers": {
#     "n": 103,
#     "source": "DELLY2_strict",
#     "genotypes": {"0/1": 78, "1/1": 25},
#     "sample_ids": ["CGA009", "CGA033", "..."]
#   },
#   "ghsl_carriers": {
#     "n": 118,
#     "source": "GHSL_v6_partition",
#     "sample_ids": ["CGA009", "CGA010", "..."]
#   },
#   "reconciliation": {
#     "all_three_agree": 95,
#     "pca_sv_agree": 103,
#     "pca_ghsl_agree": 115,
#     "pca_only": 15,
#     "pca_carrier_sv_dropout": 27,
#     "sv_carrier_pca_miss": 0,
#     "concordance_pca_sv": 0.79,
#     "concordance_pca_ghsl": 0.88,
#     "dropout_rate": 0.21
#   },
#   "group_validation": {
#     "level": "VALIDATED",
#     "reason": "OR_p=4.98e-06 + GHSL_concordant + jackknife_robust"
#   },
#   "per_sample_detail_file": "carrier_reconciliation_per_sample.tsv.gz",
#   "keys_produced": ["q3_n_carriers_pca", "q3_n_carriers_sv",
#                      "q3_carrier_concordance", "q3_n_pca_carrier_sv_ref",
#                      "q3_dropout_rate", "q3_population_prior_applied",
#                      "q3_n_rescued_by_prior"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 7: internal_dynamics.json
# Source: C01i decomposition + decomp_stats.R + cheat24
# Question: Q2
# Group requirement: UNCERTAIN (self-validating via switching pattern)
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "internal_dynamics",
#   "n_recombinant": 7,
#   "n_gene_conversion": 5,
#   "n_double_crossover": 2,
#   "recombinant_detail": [
#     {"sample": "CGA009", "crossover_bp": 9450000, "left_class": "HET",
#      "right_class": "HOM_STD", "type": "gene_conversion", "near_boundary": "left"},
#     {"sample": "CGA045", "crossover_bp": 11200000, "left_class": "HOM_INV",
#      "right_class": "HET", "type": "gene_conversion", "near_boundary": "right"}
#   ],
#   "position_distribution": {"left_20pct": 3, "center_60pct": 2, "right_20pct": 2},
#   "switching_rates": {
#     "HOM_STD": 0.0012, "HET": 0.0089, "HOM_INV": 0.0015, "REC": 0.1823,
#     "kruskal_wallis_p": 3.2e-06
#   },
#   "class_stability_pct": 96.9,
#   "class_entropy_mean": 0.312,
#   "diversity_gradient": {"r2": 0.412, "shape": "shallow_U"},
#   "phase_blocks": {
#     "n50_inside": 312000, "n50_flank": 89000, "ratio": 3.51,
#     "switch_rate_inside": 0.0012, "switch_rate_flank": 0.0089
#   },
#   "indel_concordance": {"n_classes": 3, "concordance": 0.823, "ashman_D": 2.14},
#   "founder_haplotypes": {
#     "chao1_std": 8.2, "chao1_inv": 5.4,
#     "saturated_std": true, "saturated_inv": false,
#     "n_founders_std": 7, "n_founders_inv": 4
#   },
#   "keys_produced": ["q2_n_recombinant", "q2_n_gene_conversion", ..., "q2_chao1_*"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 8: mechanism.json
# Source: cheat14 (SDs) + cheat27 (BISER2) + cheat29 (junction) + cheat21 (TE)
# Question: Q4
# Group requirement: NONE
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "mechanism",
#   "sd_evidence": {
#     "has_inverted_sd": true,
#     "left": {"start": 8190000, "end": 8205000, "length": 15000},
#     "right": {"start": 12505000, "end": 12520000, "length": 15000},
#     "identity_pct": 97.3,
#     "orientation": "inverted",
#     "alignment_length": 14200,
#     "n_mismatches": 383,
#     "divergence_pct": 2.7,
#     "biser2_concordance": "agree_nahr"
#   },
#   "junction_evidence": {
#     "left": {
#       "type": "extended_mh",
#       "mh_length": 23,
#       "mh_sequence": "ATCGATCGATCGATCGATCGATC",
#       "flanking_50bp": "ACGT...50bp...ACGT",
#       "te_name": null,
#       "te_family": null
#     },
#     "right": {
#       "type": "microhomology",
#       "mh_length": 8,
#       "mh_sequence": "GCTAAGCT",
#       "flanking_50bp": "ACGT...50bp...ACGT",
#       "te_name": null,
#       "te_family": null
#     }
#   },
#   "te_evidence": {
#     "enrichment": "NEUTRAL",
#     "enrichment_fold": 1.1,
#     "same_family": false,
#     "same_orientation": "NA",
#     "tr_enrichment": "NEUTRAL"
#   },
#   "gene_content": {
#     "n_genes_inside": 42,
#     "n_genes_spanning_bp": 1,
#     "gene_names_at_bp": ["slc6a15"],
#     "gene_density_inside": 9.8,
#     "gene_density_genome": 8.2,
#     "go_enrichment": null
#   },
#   "classification": {
#     "mechanism": "NAHR",
#     "support": "sd+junction+biser2",
#     "confidence": "3/4 = likely",
#     "decision_tree_path": "has_SD → identity>90% → inverted → BISER2_agree → NAHR"
#   },
#   "keys_produced": ["q4_has_inverted_sd", "q4_sd_identity_pct", ...,
#                      "q4_mechanism_confidence", "q4_decision_tree_path"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 9: age_evidence.json
# Source: cheat30 (GDS) + popstats (Fst/dXY/theta) + cheat19 (diversity)
# Question: Q5
# Group requirement: SUPPORTED (for Fst/dXY); NONE (for GDS)
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "age_evidence",
#   "gds_based": {
#     "gds_gap": 0.092,
#     "gds_gap_percentile": 88,
#     "gds_within_std": 0.034,
#     "gds_within_inv": 0.041,
#     "gds_between": 0.126,
#     "gds_het_pattern": "intermediate",
#     "origin_dip_p": 0.23,
#     "origin_class": "recurrent",
#     "requires_groups": false
#   },
#   "fst_based": {
#     "fst_b1b3": 0.142,
#     "requires_groups": true,
#     "groups_validation": "VALIDATED"
#   },
#   "divergence_based": {
#     "dxy_inside": 0.0089,
#     "dxy_flanking": 0.0034,
#     "dxy_ratio": 2.62,
#     "da_net": 0.0051,
#     "theta_pi_std": 0.0028,
#     "theta_pi_inv": 0.0031,
#     "theta_pi_het": 0.0045,
#     "tajima_D_std": 0.23,
#     "tajima_D_inv": 0.45,
#     "tajima_D_pooled": 2.14,
#     "segregating_sites_std": 423,
#     "segregating_sites_inv": 389,
#     "fixed_differences": 12,
#     "shared_polymorphisms": 847,
#     "requires_groups": true,
#     "groups_validation": "VALIDATED"
#   },
#   "diversity_profile": {
#     "shape": "shallow_U",
#     "r2": 0.412,
#     "requires_groups": false
#   },
#   "cross_species": {
#     "dollo_node": "not_tested",
#     "dollo_mya": null,
#     "n_species_sharing": null
#   },
#   "cross_validation": {
#     "gds_fst_spearman_rho": 0.73,
#     "gds_fst_spearman_p": 0.002,
#     "age_proxies_agree": true
#   },
#   "conclusion": {
#     "age_class": "old",
#     "evidence_for": ["GDS_88pctl", "Fst_top_third", "U_shape", "dxy_ratio>2", "12_fixed_diffs"],
#     "evidence_against": [],
#     "confidence": "ANSWERED"
#   },
#   "keys_produced": ["q5_gds_gap", "q5_gds_gap_percentile", "q5_fst_b1b3", ...,
#                      "q5_fixed_differences", "q5_shared_polymorphisms"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 10: frequency.json
# Source: C01i decomposition + C01f jackknife
# Question: Q6
# Group requirement: SUPPORTED
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "frequency",
#   "freq_inv": 0.381,
#   "genotype_counts": {
#     "HOM_STD": 89, "HET": 98, "HOM_INV": 32, "REC": 7, "UNCLASS": 0, "total": 226
#   },
#   "hwe_assessment": {
#     "expected_het": 0.4718,
#     "observed_expected_ratio": 0.919,
#     "deviation": "hwe_consistent",
#     "chisq_p": 0.312
#   },
#   "per_qgroup_freq": {
#     "Q1": 0.42, "Q2": 0.35, "Q3": 0.51, "Q4": 0.38,
#     "Q5": 0.31, "Q6": 0.40, "Q7": 0.34, "Q8": 0.29
#   },
#   "freq_cv_across_qgroups": 0.178,
#   "family_linkage": "multi_family",
#   "jackknife": {
#     "max_delta": 0.0089,
#     "sensitive_qgroup": "Q3",
#     "status": "robust_multi_family"
#   },
#   "selection_signal": {
#     "tajima_D_inside": 2.14,
#     "tajima_D_flanking": 0.23,
#     "pattern": "three_layer_consistent",
#     "caveat": "hatchery Ne~20, cannot distinguish selection from drift"
#   },
#   "keys_produced": ["q6_freq_inv", "q6_n_HOM_STD", ..., "q6_selection_pattern"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 11: hypothesis_verdict.json
# Source: C01f hypothesis tests
# Question: Q7
# Group requirement: varies by test
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "hypothesis_verdict",
#   "verdict": "confirmed_real_inversion",
#   "verdict_confidence": "high",
#   "fdr_qval": 0.0012,
#   "layers_tested": 4,
#   "layers_passed": 4,
#   "independence_class": "confirmed",
#   "supporting_tests": {
#     "t8_clair3_concordance": 0.823,
#     "t9_jackknife_status": "robust_multi_family",
#     "t9_max_delta": 0.0089,
#     "t10_theta_concordance": 0.781
#   },
#   "hypothesis_history": [
#     {"version": 1, "date": "2026-04-10", "verdict": "likely_real_inversion"},
#     {"version": 2, "date": "2026-04-12", "verdict": "confirmed_real_inversion",
#      "change_reason": "Clair3 concordance added (T8)"}
#   ],
#   "keys_produced": ["q7_verdict", "q7_verdict_confidence", "q7_fdr_qval",
#                      "q7_n_layers_passed", "q7_independence_class", "q7_t8_*", ...]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 12: morphology.json
# Source: C01a precomp inv_likeness components
# Question: Q1
# Group requirement: NONE
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "morphology",
#   "flat_inv_score": 0.8234,
#   "spiky_inv_score": 0.1102,
#   "fragmentation_score": 0.0891,
#   "s_het_mean": 0.7412,
#   "s_pve_mean": 0.6023,
#   "s_dip_mean": 0.8891,
#   "inv_likeness_mean": 0.7435,
#   "inv_likeness_max": 0.9213,
#   "inv_likeness_sd": 0.0821,
#   "family_likeness_mean": 0.1203,
#   "dosage_het_cv_mean": 0.4512,
#   "pve1_excess_mean": 0.3891,
#   "dip_pval_median": 0.0012,
#   "keys_produced": ["q1_flat_inv_score", "q1_spiky_inv_score", ..., "q1_dip_pval_median"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 13: band_composition.json
# Source: C01i decomposition + instant_q
# Question: Q1
# Group requirement: UNCERTAIN (bands come from decomposition)
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "band_composition",
#   "band_qgroup_fracs": {
#     "HOM_STD": {"Q1":0.18,"Q2":0.15,"Q3":0.08,"Q4":0.12,"Q5":0.14,"Q6":0.11,"Q7":0.13,"Q8":0.09},
#     "HET":     {"Q1":0.14,"Q2":0.12,"Q3":0.10,"Q4":0.13,"Q5":0.15,"Q6":0.12,"Q7":0.11,"Q8":0.13},
#     "HOM_INV": {"Q1":0.08,"Q2":0.10,"Q3":0.22,"Q4":0.15,"Q5":0.09,"Q6":0.14,"Q7":0.12,"Q8":0.10}
#   },
#   "dominant_qgroup_std": "Q1",
#   "dominant_qgroup_inv": "Q3",
#   "q_group_overlap": 0.72,
#   "keys_produced": ["q1_band_std_qgroup_fracs", "q1_band_het_qgroup_fracs",
#                      "q1_band_inv_qgroup_fracs", "q1_dominant_qgroup_std",
#                      "q1_dominant_qgroup_inv", "q1_q_group_overlap"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCK 14: burden.json
# Source: STEP_C01f_c_burden_regression.R
# Question: Q6 (selection_pattern)
# Group requirement: VALIDATED
# ─────────────────────────────────────────────────────────────────────────────
#
# {
#   "block_type": "burden",
#   "burden_mean_ref": 12.3,
#   "burden_mean_het": 14.1,
#   "burden_mean_inv": 18.7,
#   "burden_fold_inv_ref": 1.52,
#   "kw_p": 0.003,
#   "kw_qval": 0.012,
#   "wt_inv_vs_ref_p": 0.001,
#   "inside_outside_fold": 2.3,
#   "verdict": "INV_ELEVATED",
#   "keys_produced": ["burden_mean_ref", "burden_mean_inv", "burden_fold_inv_ref",
#                      "burden_kw_qval", "burden_verdict"]
# }


# ─────────────────────────────────────────────────────────────────────────────
# BLOCKS 15-18: Additional schemas (lighter)
# ─────────────────────────────────────────────────────────────────────────────
#
# BLOCK 15: block_detect.json
#   Source: PHASE_01C_block_detect.R
#   Contains: boundary classifications from ±80 sim_mat analysis
#   Fields: boundaries[], blue_cross_verdicts[], block_concordance_matrix
#
# BLOCK 16: triangle_insulation.json
#   Source: C01c
#   Contains: interval classifications, squareness, insulation depths,
#             sub-regime assignments, off-diagonal bridges
#
# BLOCK 17: peel_diagnostic.json
#   Source: C01f or C01n
#   Contains: L1b_effect, L2_effect, peel_delta, samples_removed, samples_gained
#
# BLOCK 18: synteny_dollo.json
#   Source: synteny_inversion_detection.sh + dollo_parsimony_inversions.py
#   Contains: mashmap breaks, dollo node, species sharing, min_age_mya


# =============================================================================
# SUMMARY: 18 BLOCK TYPES
# =============================================================================
#
#  Block                    Source          Question  Groups needed
#  ─────────────────────────────────────────────────────────────────
#  1  existence_layer_a     C01a/b/d        Q1, Q7    NONE
#  2  existence_layer_b     C00/4D/4E/4G    Q7        NONE
#  3  existence_layer_c     C04 GHSL        Q7        NONE
#  4  existence_layer_d     STEP03          Q7        IS validation
#  5  boundary_left         C01g            Q3        NONE
#  6  carrier_reconciliation audit+regeno   Q3,Q6,Q7  N/A
#  7  internal_dynamics     C01i+stats      Q2        UNCERTAIN
#  8  mechanism             cheat14/27/29   Q4        NONE
#  9  age_evidence          cheat30+popstat Q5        SUPPORTED*
#  10 frequency             C01i+C01f       Q6        SUPPORTED
#  11 hypothesis_verdict    C01f            Q7        varies
#  12 morphology            C01a            Q1        NONE
#  13 band_composition      C01i+Q          Q1        UNCERTAIN
#  14 burden                burden_regr     Q6        VALIDATED
#  15 block_detect          01C             Q1,Q3     NONE
#  16 triangle_insulation   C01c            Q1,Q3     NONE
#  17 peel_diagnostic       C01f/C01n       Q7        NONE
#  18 synteny_dollo         mashmap+dollo   Q5        NONE
#
#  * age_evidence: GDS sub-block needs NONE; Fst/dXY sub-block needs SUPPORTED
#
# =============================================================================
