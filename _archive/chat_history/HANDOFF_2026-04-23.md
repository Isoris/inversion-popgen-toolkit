# HANDOFF — Next chat (continuing from 2026-04-23 packaging session)

**User:** Quentin Andres (Kasetsart Univ Bangkok, PhD computational genomics)
**Previous session date:** 2026-04-23, evening Bangkok time
**HPC status:** LANTA offline until ~2026-04-27 Monday
**Working channel:** laptop + GitHub Desktop (no git CLI on HPC, no git on laptop that the user uses — everything via GitHub Desktop GUI)
**User's mood at session end:** energized enough to continue wiring, but recognized this chat was getting long and staleness issues were creeping in

---

## 🚨 CRITICAL — FILE STALENESS PROBLEM

**The end-of-session discovery**: user uploaded `compute_candidate_status.R`
at the START of session (705 lines). At the END of session, user re-uploaded
the SAME file — but it was now **977 lines** (272 lines bigger, different md5).

This means the user's repo has been evolving during the session — so
ANY R or Python file I pulled from tarballs earlier may be stale against
the current repo state.

**The axis5 tarball I built just before this handoff was based on the 705-line
version. I deleted it because dropping it would have destroyed 272 lines of
newer code.**

**ALWAYS RE-ASK THE USER TO UPLOAD THE CURRENT VERSION before patching any
R or Python file.** Do not assume what was uploaded earlier in a previous
chat is still current.

Files to treat as potentially stale (re-ask for current versions):
- `compute_candidate_status.R` — known stale (977 vs 705)
- `characterize_candidate.R` — unknown, probably stale
- `BREEDING_A_broodstock_compatibility.R`, `BREEDING_C_*`, `BREEDING_D_*` — unknown
- `cheat27_sd_nahr_substrate.R`, `cheat28_tandem_repeat_context.R`,
  `cheat29_junction_forensics.R`, `cheat30_gds_by_genotype.R`,
  `cheat30_ibs_by_genotype.R` — unknown
- `02_ancestral_fragments.R` — was updated mid-session from v1.0 → v1.1,
  probably current now
- All scripts in phase_7_cargo/ and catfish-synteny-toolkit/ — unknown

Files that are DEFINITELY current (Claude shipped them this session):
- `assign_structural_class_v7.py`, `bnd_sided_support.py`,
  `bp_pipeline_bridge.py`, `cheat29b_assembled_junction.py`,
  `cross_species_bridge_v6.py` — Claude shipped these in the phase4 bundle
- `_axis5_final_label.R`, `_lib_final_labels.R` — Claude wrote these this session
- All 5 new JSON schemas in `registries/schemas/structured_block_schemas/`
- `run_phase4_v5_new_evidence.sh`, `run_phase4_v7_wiring.sh` — Claude
  patched/wrote these this session

---

## WHAT WAS PUSHED TO GITHUB THIS SESSION (all confirmed by user)

### Repo 1: `Isoris/inversion-popgen-toolkit` (existing)

**Commit `ce8949c` — Updated_phase_qc_shelf_v3.1** (41 files changed, +6167 -307)
- phase_qc_shelf v2.1 + v3.0 + v3.1 layers
- Engine H (Hobs per-group) layer — shipped in same tarball but untested on LANTA
- Bug 3 (QC_OUT default path) fixed in `00_config.sh`
- Bug 1 (install_update LIVE_MODULE path) fixed

**Commit (next, after v3.1) — phase4 v7 wiring** (68 new/modified files)
- `breakpoint_pipeline/02_ancestral_fragments.R` — v1.1 registry-wired + HSM
  patched (582 lines; HSM from Bickel & Frühwirth 2006, needs `modeest` CRAN pkg)
- `inversion_modules/phase_4_postprocessing/4d_group_dependent/`:
  - cheat29b_assembled_junction.py
  - bnd_sided_support.py
  - cross_species_bridge_v6.py
  - bp_pipeline_bridge.py
  - run_phase4_v5_new_evidence.sh (4 bugs fixed vs raw bundle: v6/v7 script
    names; dropped --query_species/--sister_species; added STEP 3b
    bp_pipeline_bridge; added DOLLO_TSV env var)
- `inversion_modules/phase_4_postprocessing/4e_final_classification/`:
  - assign_structural_class_v7.py
  - _axis5_final_label.R (standalone helper)
  - run_phase4_v7_wiring.sh (master orchestrator)
- `phase_7_cargo/` (48 files) — FULL folder integrated as subfolder of the
  main repo. Includes breeding/BREEDING_A/C/D, all extra_plots, tables,
  launchers, compute, _lib_final_labels.R helper.
- `registries/schemas/structured_block_schemas/`: 5 new schemas
- `docs/`: 5 new phase4 docs (HANDOFF.md, SESSION_AUDIT_append.md,
  kde_mode_estimation_review.md, phase4_framework_v7_correction.md,
  toolkit_audit.md)
- `LOCAL_INSTALL_README.md` at repo root

### Repo 2: `Isoris/catfish-synteny-toolkit` (NEW)

Created this session. 28 files. Includes:
- README.md with tool table (wfmash, miniprot, TRF, IRF, minimap2, etc.)
- HANDOFF.md, DEPENDENCIES.md, MANUSCRIPT_SKELETON.md
- LICENSE (MIT), .gitignore
- run_pipeline.sh
- config/ (3 files: species_manifest, ncbi_accessions, species_tree.nwk)
- scripts/ (18 STEP_ scripts)

Description field on GitHub:
> Graph-based structural rearrangement detection and validation across catfish genomes

Loose coupling with inversion-popgen-toolkit via 3 output files:
- `output/step_02/events.bed`
- `output/step_09c/flank_coherence.tsv`
- `output/step_11/polarized.tsv`

These are consumed by `cross_species_bridge_v6.py` in the other repo.

---

## WHAT'S STILL PENDING (session scope, not yet done)

### Step 6a — Axis 5 in compute_candidate_status.R 🚨 REDO NEEDED
Claude built an axis5_activation_local.tar.gz based on the OLD 705-line
version. DELETED because user uploaded 977-line current version too late.
**Next chat: re-ask for current compute_candidate_status.R, re-apply the
Axis 5 patch to THAT version.**

Design of the Axis 5 patch:
- Insert between `status_dt <- rbindlist(...)` and `fwrite(status_dt, ...)`
- Read `V7_FINAL_DIR` env var
- If set + dir exists: source `_axis5_final_label.R` helper (already in
  repo at `inversion_modules/phase_4_postprocessing/4e_final_classification/`)
- Call `append_axis5_to_rows(status_dt, v7_final_dir)` — adds 4 columns:
  q_overall_structural_class, axis5_weakest_component, axis5_justification,
  axis5_source
- If V7_FINAL_DIR unset → no-op, original behavior preserved
- Helper resolution: V7_AXIS5_LIB env var → 4 relative-path fallbacks

### Step 6b — Gap 2 of HANDOFF step 8
Add GDS consistency check in `characterize_q5()`:
```r
gds_ok <- keys$q5_hom_inv_gds_clean == "clean"
if (isTRUE(gds_ok) && !is.na(age_point_estimate)) {
  age_class <- "answered"
} else {
  age_class <- "measured"
}
```
Requires current `characterize_candidate.R` from user.

### Step 7 — BREEDING_A/C/D filter by v7 label
Add 8-line snippet near top of each BREEDING script (after `cands <- fread(...)`):
```r
V7_FINAL_DIR <- Sys.getenv("V7_FINAL_DIR", "")
if (nzchar(V7_FINAL_DIR) && dir.exists(V7_FINAL_DIR)) {
  source("phase_7_cargo/extra_plots/compute/_lib_final_labels.R")
  labels <- load_final_labels(V7_FINAL_DIR)
  usable_cids <- filter_usable_candidates(labels)
  before <- nrow(cands)
  cands <- cands[candidate_id %in% usable_cids]
  cat(sprintf("[v7 filter] %d → %d candidates\n", before, nrow(cands)))
}
```
The helper `_lib_final_labels.R` is already shipped. Just need to edit
3 BREEDING files. **Re-ask for current versions.**

### Step 8 — Close 4 wiring gaps from toolkit_audit.md
- Gap 1: Q7B keys from breakpoint_evidence_audit.py — extend q7 key list
  in `build_key_spec()` of compute_candidate_status.R (~10 lines)
- Gap 2: cheat30 GDS consistency → covered in Step 6b above
- Gap 3: Axis 5 → covered in Step 6a above
- Gap 4: gene_conversion_detector into characterize_q2 — similar key list
  extension

All of these need the CURRENT version of the file being edited.

### Step 2 + Step 9 — HPC-ONLY (wait until 2026-04-27)
- Step 2: Run `breakpoint_pipeline` on LG28 to generate 3 TSV outputs
  (consensus, fragments, per_method) that bp_pipeline_bridge.py consumes
- Step 9: Full 226-cohort v7 end-to-end run via `sbatch run_phase4_v7_wiring.sh`

### Step 10 — Tag v7.0 after HPC validates

---

## MONDAY HPC CHECKLIST (for user when LANTA returns)

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit
git pull                                     # pull v3.1 + phase4

# Install modeest for HSM
Rscript -e 'if (!require(modeest)) install.packages("modeest", repos="https://cloud.r-project.org")'

# Also clone the new catfish-synteny repo
cd ..
git clone https://github.com/Isoris/catfish-synteny-toolkit.git

# HANDOFF step 2: run breakpoint_pipeline on LG28 to make TSVs
cd inversion-popgen-toolkit/breakpoint_pipeline
bash run_pipeline.sh --chrom C_gar_LG28                      # adjust per actual CLI
# verify HSM columns appear
head -1 results/candidate_breakpoints_consensus.tsv | tr '\t' '\n' | \
  grep -E "mode_(hsm|kde|agreement|method)"

# HANDOFF step 9: full phase4 v7 run
cd ..
BASE=/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04 \
BP_OUT=${BASE}/breakpoint_pipeline_output \
  sbatch inversion_modules/phase_4_postprocessing/4e_final_classification/run_phase4_v7_wiring.sh

# After completion
awk -F'\t' 'NR>1 {print $2}' ${BASE}/phase4_v7_blocks/final/final_catalog.tsv | sort | uniq -c | sort -rn
```

Expected distribution: most candidates in `supported_balanced_inversion_*`,
small fraction in `complex_rearrangement_out_of_scope`, small fraction in
`family_confounded_locus`, some in `supported_shared_between_species_inversion`.

If `complex_rearrangement_out_of_scope` > 50% → threshold too aggressive,
drop rule 1 threshold in assign_structural_class_v7.py from 0.70 → 0.80.

---

## KEY FILE PATHS

On HPC:
- Project root: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`
- Main repo: `${BASE}/inversion-popgen-toolkit/`
- Synteny repo (to clone Monday): `${BASE}/catfish-synteny-toolkit/`
- conda env: `assembly`
- SLURM account: `lt200308`

On GitHub:
- https://github.com/Isoris/inversion-popgen-toolkit (v3.1 + phase4 pushed)
- https://github.com/Isoris/catfish-synteny-toolkit (NEW, pushed)

User's workflow:
- Edits on laptop via text editor / GitHub Desktop drag-drop
- Pushes via GitHub Desktop GUI
- Pulls on HPC via git CLI (user doesn't use git CLI locally)
- **No git access on HPC except pull/fetch; no push from HPC**

---

## SCIENTIFIC CONTEXT (for new Claude to understand)

**Project**: `MS_Inversions_North_african_catfish` — population genomics
study detecting chromosomal inversions in a **226-sample pure Clarias
gariepinus hatchery broodstock cohort**, sequenced at ~9× Illumina coverage
across 28 linkage groups.

**K clusters** in this cohort = hatchery broodline structure, **NOT hybrid
or population structure**. Repeat: this is a pure-species cohort.

**Three cohorts that MUST NEVER be conflated:**
1. F1 hybrid (C. gariepinus × C. macrocephalus) — ONLY for genome assembly paper
2. 226-sample pure C. gariepinus hatchery — current inversion paper ← THIS ONE
3. Pure C. macrocephalus wild cohort — future paper

**Target venue**: Nature Communications / Genome Research tier
**Flagship finding**: LG28 inversion at C_gar_LG28:15,115,000–18,005,000,
~2.89 Mb sub-telomeric, karyotype 60/106/60, shelf Fst 0.308 vs flanking
~0.032–0.055. But **LG28 is a test case, not the flagship**. Full-cohort
v7 output will identify cleaner candidates for the manuscript centerpiece.

**Key architectural principle**: "Wiring, not new code." v6 tried to add
a new interior-structure diagnostic; v7 correctly restructured to wire
the user's existing breakpoint_pipeline outputs instead of re-deriving
them. **DO NOT write new analysis modules; only connect existing ones.**

---

## USER PREFERENCES & COMMUNICATION STYLE

- **Pushes back productively** when Claude is wrong — saves repeated bugs.
  Examples this session:
  - "its already v1.1 but not kde fixed" → caught stale 02_ancestral_fragments.R
  - "where can I get the phase 7 cargo?" → caught incomplete phase_7_cargo bundle
  - "I FEEL THAT YOU HAVE THE WRONG VERSIONS OF FILES" → caught 705 vs 977 line staleness
- Uses GitHub Desktop GUI, **not git CLI on laptop**
- Prefers concise, direct answers; doesn't want lots of hand-holding
- Is juggling 4 projects at once ("teacher asks 1 student to do like 4 projects at once")
- Is fresh/energetic for this sync work, treats it as housekeeping before harder science
- Will say "idk" or "lol" when unsure — friendly register, keep tone matched
- Sometimes mistypes ("HSA" meant "yes"); don't over-interpret

---

## THREE-SNAKE CONTINUITY FRAMEWORK (context for analysis)

The inversion-popgen-toolkit uses a three-snake framework for candidate
validation. Brief primer for new Claude:

- **Snake 1**: Local PCA/MDS continuity across candidates
- **Snake 2**: Sample-grouping continuity (karyotype-consistent grouping)
- **Snake 3**: GHSL haplotype-contrast within-sample divergence (v11+)

Engines:
- **Engine A**: NGSadmix + PCAngsd
- **Engine B**: C++ fixed-F EM (instant_q)
- **Engine H**: Hobs per-group (new in v3.1)

The phase_qc_shelf is the QC substrate that feeds these engines.

---

## WHAT TO DO IN THE NEXT CHAT (PRIORITIZED)

### Priority 1 — Re-do Axis 5 patch against CURRENT file (977 lines)
Ask user to upload current compute_candidate_status.R. Apply the 45-line
Axis 5 block between `compl_dt <- rbindlist(...)` and
`fwrite(status_dt, ...)`. Deliver as drag-drop tarball.

### Priority 2 — Step 7 (BREEDING filter)
Ask user to upload current BREEDING_A/C/D. Add 8-line snippet to each.
One tarball for all three.

### Priority 3 — Step 8 Gap 1 (Q7B audit keys) + Gap 4 (gene conversion)
Ask user to upload current compute_candidate_status.R (same one from
Priority 1). Extend `build_key_spec()` q7 and q2 lists. Consolidate with
Priority 1 into ONE commit if possible.

### Priority 4 (only if user wants) — Step 6b (Gap 2, GDS consistency)
Ask user to upload current characterize_candidate.R. Add ~6 lines in
characterize_q5().

### Wait for HPC Monday
Steps 2 and 9 can't be done on laptop.

---

## PRODUCT / TOOLING NOTES FOR NEXT CLAUDE

- User uses GitHub Desktop. Drag-drop bundles ("staging" folders packaged
  as .tar.gz) are the preferred delivery format.
- Always PUT THE EDITED FILE IN THE EXACT REPO PATH in the staging folder
  so drag-drop auto-merges correctly.
- Always generate a short DROP_README.md or LOCAL_INSTALL_README.md with
  commit message pre-written.
- Always include a bash -n / python ast.parse / backtick-aware R balance
  check before shipping.
- User has R installed locally. Code must be valid R.
- The user knows his stack very well. When he pushes back, investigate
  immediately — don't explain why I might be right.

---

## DELIVERABLES STATE AT CHAT END

Files in `/mnt/user-data/outputs/` at end of session (may expire):
- `phase_qc_shelf_v3.1_local.tar.gz` (122 KB) — DONE, pushed
- `phase4_v7_wiring_local.tar.gz` (134 KB) — DONE, pushed
- `catfish_synteny_toolkit_local.tar.gz` (65 KB) — DONE, pushed
- ~~`axis5_activation_local.tar.gz`~~ — DELETED, was based on stale input

User's current state:
- ✅ All three repos synced with latest session work
- ⏳ Axis 5 wiring deferred to next chat (needs current 977-line file)
- ⏳ Step 7 BREEDING filter deferred to next chat
- ⏳ Step 8 gaps deferred to next chat
- ⏳ HPC validation deferred to Monday 2026-04-27
