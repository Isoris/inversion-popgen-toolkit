# MODULE_MAP.md — how the two module trees fit together

*Last updated: 2026-04-24*

The repo has **two** parallel module trees with overlapping purposes.
This document explains which is which and maps both to the manuscript
(`MS_Inversions_North_african_catfish_v19.docx`).

Renumbering is deliberately **not** proposed — too disruptive. The goal
here is to let a new reader (or future-Quentin) navigate without
tripping on the duplication.

---

## Two trees, two purposes

### `Modules/` — preprocessing + caller outputs (MODULE_1..4)

Wet-lab → BAM → SNP panel → per-sample SV caller VCFs. Every module is
self-contained, SLURM-driven, and produces a well-defined artefact
consumed by `inversion_modules/` or by the manuscript directly.

### `inversion_modules/` — the inversion-discovery pipeline (phase_1..7)

Takes the outputs of `Modules/` and runs the actual inversion analysis:
discovery → breakpoint refinement → postprocessing → gene-content
analysis. Each `phase_N_*/` directory is **one logical stage of the
paper's §3 results**.

### Why not one tree?

Historical: `Modules/` was built first as general-purpose preprocessing
for the whole lab pipeline (SNP discovery, structure, SV calling).
`inversion_modules/` was spun up later when the inversion work grew
beyond one module and needed its own multi-phase pipeline. The
old-tree / new-tree split is documented here rather than merged because
(a) `Modules/` is reused by other planned papers (pure *C. macrocephalus*
cohort, hybrid assembly) and (b) merging would shuffle ~100 file paths.

---

## `Modules/` — map

| Dir | Role | Manuscript section |
|---|---|---|
| `MODULE_1_read_prep` | Read QC, species verification (Mash vs both parental subgenomes), alignment, dedup, depth QC | §2.1 |
| `MODULE_2A_snp_discovery` | ANGSD biSNP panel, BEAGLE GLs, PCAngsd, NGSadmix, evalAdmix, ngsRelate pruning | §2.2 |
| `MODULE_2B_structure` | Extended K=2–20 ancestry, per-chromosome structure, kinship figures, sample registry | §2.3 |
| `MODULE_3_heterozygosity_roh` | Per-sample H, multiscale theta, ngsF-HMM ROH, FROH | §2.4 |
| `MODULE_4A_clair3_snp_indel` | Clair3 hard-genotype SNP/INDEL catalog (separate from 2A GL-based panel) | §2.5 |
| `MODULE_4B_delly_del` | DELLY2 deletions | §3.1 (SV catalogs) |
| `MODULE_4C_delly_dup` | DELLY2 duplications | §3.1 |
| `MODULE_4D_delly_inv` | DELLY2 inversions (primary SV-caller evidence for inversions) | §3.1 |
| `MODULE_4E_delly_bnd` | DELLY2 breakends | §3.1 |
| `MODULE_4F_delly_tra` | DELLY2 translocations | §3.1 |
| `MODULE_4G_manta_sv` | Manta (assembly-based) SV caller — independent second caller | §3.1 |
| `MODULE_CONSERVATION_CORE` | Conservation/deleterious variant annotation: bcftools csq + SIFT4G + VESM | §4 (burden analyses) |
| `MODULE_CON_conservation_deleterious` | Three-track conservation pipeline: trackA/B/C | §4 |
| `Others/MODULE_6_founder_packs` | **Placeholder** (.gitkeep only) — downstream founder-diversity analysis | future |
| `Others/MODULE_PAV_presence_absence` | **Placeholder** (.gitkeep only) — PA-Roary-like presence/absence tool | future |

Note on numbering holes: there is no `MODULE_5*` in `Modules/`.
MODULE_5A1/5A2/5A3/5B/5C/5D/5E all migrated into `inversion_modules/`
under phase-named directories (see next section). The `5*` numbers
survive in filenames and comments as provenance pointers.

---

## `inversion_modules/` — map

Each phase is one stage in the inversion pipeline. `_archive/` under
this root holds legacy v8.5 HPC code kept for diff reference.
`_archive_superseded/` holds dead-end directions dropped from the
pipeline.

| Dir | Legacy name | Role | Manuscript section |
|---|---|---|---|
| `phase_1_inputs/` | MODULE_5A1 | Callable masks, ANGSD SAF/SFS, BEAGLE inputs | §3.2 |
| `phase_2_discovery/` | MODULE_5A (discovery portions) | Genome scan for inversion-like regions — **5 sub-blocks (2a–2e)** | §3.3 |
| `phase_3_refine/` | MODULE_5A2_breakpoint_validation | Bring SV-caller evidence into the registry for phase_7/validation's VALIDATED gate | §3.4 |
| `phase_4_catalog/` | (old 4a_existence_layers) | Catalog birth: C01d scoring + C01e blocks + C01g boundary-unification | §3.5 |
| `phase_5_qc_triage/` | (old 4b_qc_triage, formerly phase_qc_shelf) | Data-quality QC per candidate; emits soft `q_qc_shelf_flag` | §3.5 |
| `phase_6_breakpoint_refinement/` | (old 4c_breakpoint_refinement, formerly breakpoint_pipeline) | bp-resolution refinement → `candidate_breakpoints_consensus.tsv` | §3.5 |
| `phase_7_karyotype_groups/` | (old 4d_group_proposal + 4e_group_validation) | `proposal/` (C01i) + `validation/` (C01f); emits SUPPORTED/VALIDATED/SUSPECT gate | §3.5 |
| `phase_8_evidence_biology/` | (old 4f_group_dependent) | Mechanism (Q4) + age/origin (Q5) + existence audit (Q7) + bp_bridge | §3.5 |
| `phase_9_classification/` | (old 4g_final_classification) | characterize_candidate + compute_candidate_status + axis 5 + 7-Q tier | §3.5 |
| `phase_10_followup/` | MODULE_5A_Discovery_Visualisations | Per-candidate deep analysis (dosage rasters, region-grow plots, diagnostic figures) | §3.6 |
| `phase_11_secondary/` | MODULE_5C/5D | LD / Fst secondary analyses (5E archived 2026-04-24; still under legacy names inside) | §3.7 |
| `phase_12_cargo/` | MODULE_6_Cargo | Gene content + per-arrangement evolution (breeding implications) | §4 |
| `utils/` | — | Shared helpers (ancestry_bridge, sample_map, theme plate) | — |

> **Layout note 2026-04-24 (pass 15).** `phase_4_postprocessing/` was
> flattened into 6 dedicated top-level phases (phase_4_catalog through
> phase_9_classification) to surface the catalog → QC → breakpoint →
> group_proposal → group_validation → evidence_biology → classification
> chain. Downstream phases (followup/secondary/cargo) shifted to 10/11/12.
> See `TWELVE_PHASE_RENAME_PROPOSAL.md` for the rationale and
> `SESSION_AUDIT_2026-04-24_chat-3-end.md` for prior-state history.

### `phase_2_discovery/` sub-blocks

The five 2*x* blocks form a DAG:

| Block | Role | README |
|---|---|---|
| `2a_local_pca/` | BEAGLE → dosage → per-chr local PCA windows | ✅ (added 2026-04-24) |
| `2b_mds/` | lostruct distances + per-chr MDS → candidate region seeds | ✅ (added 2026-04-24) |
| `2c_precomp/` | SV-evidence prior + window precompute + seeded regions + landscape detector | ✅ (existing, 411 lines) |
| `2d_candidate_detection/` | Staircase detector — primary boundary track | ✅ (existing, 326 lines) |
| `2e_ghsl/` | GHSL haplotype-contrast (Layer C in the 4-layer framework) | ✅ (existing, 175 lines) |

Data flow:
`2a` → `2b` → `2c` → `2d` + `2e` → `phase_3_refine/` → `phase_4_catalog/` → … → `phase_9_classification/`

### `phase_3_refine/` scripts

A 6-step chain (+ config + orchestrator). Each step is tagged with the
evidence Layer it **feeds** (Layer A/B/C/D in the 4-layer framework):

| Step | Role | Layer |
|---|---|---|
| `00_phase3_config.sh` | Module config (sources master `00_inversion_config.sh`) | — |
| `STEP_A01_extract_inv_candidates.sh` | DELLY+Manta INV extraction, dedup, match to phase-2d candidates | A upstream |
| `STEP_A02_extract_breakpoint_evidence.py` | Per-sample pysam BAM evidence + REF/HET/INV group assignment | A upstream (→ feeds D) |
| `STEP_D03_statistical_tests_and_seeds.py` | Fisher / χ² / Armitage tests + Layer D registry writes + C01i seeds | **D** |
| `STEP_D04_validation_plots.py` | Publication figures (OR forest, contingency, evidence heatmap) | D diagnostic |
| `STEP_B05_delly_manta_concordance.py` | DELLY × Manta concordance report + manuscript sentences | **B** |
| `STEP_B06_bnd_rescue.py` | Orphan BND pair rescue + `existence_layer_b_bnd_rescue` writes | **B** (supplementary) |
| `run_breakpoint_validation.sh` | SLURM orchestrator chaining 01→06 | — |
| `annotate_population_confidence.sh` | Side-tool: add confidence scores to ALL PASS SVs (not a filter) | — |
| `breakpoint_validator_standalone.py` | Single-candidate standalone validator | — |

No `STEP_C*` script exists in phase 3 because Layer C (GHSL) is
produced in `phase_2_discovery/2e_ghsl/`, not here.

### `phase_4_catalog/` → `phase_9_classification/` (post pass 15)

After pass 15, the old `phase_4_postprocessing/` umbrella no longer
exists; the seven-block DAG is now six top-level phases (4→9) with
each phase's README giving its own role and producer/consumer contract.
The strict sequential DAG is preserved — each stage still consumes the
output of the previous and writes into the shared registries:

| Phase | Role |
|---|---|
| `phase_4_catalog/` | Catalog birth: C01d scoring + C01e blocks + C01g boundary-unification |
| `phase_5_qc_triage/` | Data-quality QC per candidate: SNP density, BEAGLE uncertainty, coverage CV, per-group Hobs; emits `q_qc_shelf_flag` (clean / low_snp / high_uncertain / coverage_artifact / messy). Soft gate — messy candidates still pass through. |
| `phase_6_breakpoint_refinement/` | bp-resolution breakpoint refinement via dosage signal + per-carrier ancestral fragment distribution → `candidate_breakpoints_consensus.tsv` with `final_left_bp`, `final_right_bp`, CI bounds. |
| `phase_7_karyotype_groups/proposal/` | C01i decompose / multi_recomb / nested_comp / seal (writes UNCERTAIN) |
| `phase_7_karyotype_groups/validation/` | C01f hypothesis tests + gate; reads Layer D from phase_3 for VALIDATED promotion |
| `phase_8_evidence_biology/` | Q4 mechanism + Q5 age/origin + Q7 existence audit + bp_bridge (scripts grouped by question in `q4_mechanism/`, `q5_age_and_origin/`, `q7_existence_audit/`, `bp_bridge/`) |
| `phase_9_classification/` | characterize_candidate + compute_candidate_status + axis 5 (wired via `V7_FINAL_DIR` env) + q_qc_shelf_* reader; final tier + 7-question characterization |

Cross-phase material lives in `registries/schemas/` (structured-block
schemas + specs) and at the repo-level `docs/` (architecture notes,
design rationale).

---

## Quick lookup — "which module does X?"

| Question | Where |
|---|---|
| "Where are my filtered BAMs?" | `Modules/MODULE_1_read_prep/` |
| "Where is the SNP panel?" | `Modules/MODULE_2A_snp_discovery/` |
| "Where is global Q (K=8)?" | `Modules/MODULE_2B_structure/` |
| "Where are DELLY inversion calls?" | `Modules/MODULE_4D_delly_inv/` |
| "Where is the inversion candidate catalog?" | `inversion_modules/phase_4_catalog/` |
| "Where does bp-resolution refinement happen?" | `inversion_modules/phase_3_refine/` + `inversion_modules/phase_6_breakpoint_refinement/` |
| "Where are per-candidate figures?" | `inversion_modules/phase_10_followup/` + `phase_12_cargo/plot/` |
| "Where is the Hobs secondary confirmation?" | `inversion_modules/phase_5_qc_triage/STEP_Q07b_hobs_per_group.sh` + `STEP_Q07c_hobs_windower.sh` (per-group Hobs; MODULE_5E archived 2026-04-24) |
| "Where is the burden analysis?" | `Modules/MODULE_CONSERVATION_CORE/` + `phase_8_evidence_biology/` |

---

## Known numbering oddities (not fixed — documented)

1. **`Modules/MODULE_5*` missing.** Migrated to `inversion_modules/`.
   Filenames like `STEP_5A2_*.R` still carry the old number for
   provenance.
2. **`inversion_modules/phase_11_secondary/MODULE_5C` and `MODULE_5D`.**
   Still carry `MODULE_5*` folder names internally. Consistent with the
   legacy naming; renaming would break several self-relative paths in
   their SLURM launchers. Tracked for the HPC-coordinated Cat 3
   RENAMING pass (see `phase_2_discovery/2c_precomp/RENAMING.md`).
   `MODULE_5E_Inversion_HOBS` was archived on 2026-04-24 to
   `_archive_superseded/MODULE_5E_Inversion_HOBS_superseded_by_Q07b/`
   — Hobs confirmation is now done per-karyotype-group by
   `phase_5_qc_triage/STEP_Q07b + STEP_Q07c`.
3. **`Modules/Others/` is placeholders only.** `MODULE_6_founder_packs`
   and `MODULE_PAV_presence_absence` contain only `.gitkeep`. The real
   work for these (if/when done) will likely land in
   `inversion_modules/phase_12_cargo/` (founder packs) or a new
   standalone module (PAV).
4. **Phase restructuring history.** (a) Pass 12 (2026-04-24) folded the standalone `breakpoint_pipeline/` and `phase_qc_shelf/` modules into `phase_4_postprocessing/` as sub-blocks 4c and 4b respectively. (b) Pass 15 (2026-04-24) then flattened the entire `phase_4_postprocessing/` umbrella into 6 dedicated top-level phases (phase_4_catalog through phase_9_classification) and shifted downstream phases (followup/secondary/cargo) to 10/11/12. The data-dependency chain is now enforced by folder numbering.
5. **Double-tree for phase 2 discovery**: see
   `inversion_modules/phase_2_discovery/` vs
   `inversion_modules/_archive/MODULE_5A2_Discovery_Core/`. The archive
   is the frozen v7-era pre-refactor copy. Live work is in `phase_*`.

---

## Missing READMEs

None as of 2026-04-24. Every live directory in the two module trees
has a README. When adding new subfolders, add a short README at the
same time.
