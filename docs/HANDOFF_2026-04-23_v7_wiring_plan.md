# Handoff — installing and wiring `phase4_chat_final/`

**Audience:** next session (user working on LANTA, HPC back Monday 2026-04-27).
**Goal:** take the 6 scripts + 5 schemas + 1 launcher from this bundle
and wire them into `inversion-popgen-toolkit/` so everything is
registry-backed, reproducible, and callable from SLURM.

**Non-goal:** do not write new analysis code. v6 already almost broke
this rule (interior-structure diagnostic was a new module when user
asked for wiring — v7 corrected it). Stay disciplined.

---

## Prerequisites checklist

Before starting, confirm on LANTA:

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/

# 1. Toolkit repo is cloned and current
cd inversion-popgen-toolkit/
git status                                   # should be clean or expected WIP
git log --oneline -5                         # confirm current HEAD

# 2. breakpoint_pipeline is at v1.0 or later with scripts 01-07 present
ls breakpoint_pipeline/01_dosage_signal.R \
   breakpoint_pipeline/02_ancestral_fragments.R \
   breakpoint_pipeline/03_consensus_merge.R \
   breakpoint_pipeline/04_diagnostic_figure.R \
   breakpoint_pipeline/05_four_encoding_diagnostic.R \
   breakpoint_pipeline/06_regime_stream_graph.R \
   breakpoint_pipeline/07_breakpoint_evidence_pileup.py
# All 7 should exist

# 3. Phase 4 structure is intact
ls inversion_modules/phase_4_postprocessing/{4a_existence_layers,4b_group_proposal,4c_group_validation,4d_group_dependent,4e_final_classification}

# 4. Registry API is importable
python3 -c "
import sys
sys.path.insert(0, 'registries/api/python')
from registry_loader import load_registry
reg = load_registry()
print('registry OK:', type(reg).__name__)
"

# 5. conda env 'assembly' is active and has modeest
conda activate assembly
Rscript -e 'if(!require(modeest)) install.packages("modeest", repos="https://cloud.r-project.org"); library(modeest); cat("modeest OK\n")'
```

If any of these fail, stop and fix them before proceeding. The wiring
assumes all of these.

---

## Step 1 — Install the HSM patch into `02_ancestral_fragments.R`

**Why first:** the `bp_pipeline_bridge.py` (step 5) reads the TSVs
produced by this script. Installing the patch first means all
downstream consumers see the improved mode estimates.

### 1a. Back up the current script

```bash
cd breakpoint_pipeline/
cp 02_ancestral_fragments.R 02_ancestral_fragments.R.bak_pre_hsm
```

### 1b. Add the `library(modeest)` line

Open `02_ancestral_fragments.R`. Find the `suppressPackageStartupMessages`
block (near the top, around line 15–20 in the current version):

```r
suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  ...
})
```

Add `library(modeest)` as a new line inside the block. Order doesn't
matter.

### 1c. Replace `summarize_fragment_side()`

The current function spans lines 132–206 (start: `summarize_fragment_side <- function(`;
end: the closing `}` of the outer function before the next `# =====`
divider).

Delete those lines. Paste in the contents of
`scripts/02_ancestral_fragments_hsm_patch.R` **starting from the
`summarize_fragment_side <- function(` line (line 25 of the patch,
after the header comment).**

Do NOT paste the patch's header comment block — just the function
body. The existing file already has enough comments; the header of
the patch file is for standalone readability only.

### 1d. Sanity-check R syntax

```bash
Rscript -e 'source("02_ancestral_fragments.R")' 2>&1 | head -5
# Should either print nothing or an error pointing to a specific line.
# If no error, the parse is clean.
```

If there's a syntax error, the most likely cause is a mismatched
brace when pasting. Diff against the backup:

```bash
diff 02_ancestral_fragments.R.bak_pre_hsm 02_ancestral_fragments.R | head -30
```

### 1e. Smoke-test on a small synthetic run

The breakpoint_pipeline repo has `run_pipeline.sh` — run it on ONE
candidate (not all of LG28) to verify the patch doesn't crash:

```bash
# Use whatever tiny test case the pipeline has; if none, skip this
# and go straight to full LG28 in step 2 — the patch is a drop-in
# replacement and the R syntax check above is usually enough.
```

### 1f. Commit the patch

```bash
cd ..
git add breakpoint_pipeline/02_ancestral_fragments.R
git commit -m "Replace Silverman-KDE mode with HSM (Bickel & Frühwirth 2006)

Silverman's rule oversmooths sharp-peaked right-skewed distributions,
which is the shape of per-carrier fragment boundary distributions
(sharp peak at true breakpoint + recombinant tail).

Replace summarize_fragment_side() to use modeest::hsm() as primary
mode estimator. Keep KDE as secondary for sensitivity check. Add
three diagnostic fields: mode_hsm_bp, mode_kde_bp, mode_agreement_kb.

Existing mode_bp, ci_low, ci_high semantics preserved; 03_consensus_merge
requires no changes.

Adds CRAN dependency: modeest."
```

---

## Step 2 — Run the full breakpoint_pipeline on LG28

This produces the three TSVs that the bridge script consumes.

```bash
cd breakpoint_pipeline/
bash run_pipeline.sh --chrom C_gar_LG28 --candidates CANDIDATE_TABLE.tsv
# or whatever the actual invocation is per its README
```

Expected output files (in whatever outdir run_pipeline.sh uses):
- `candidate_breakpoints_consensus.tsv`
- `candidate_ancestral_fragments.tsv`
- `candidate_breakpoints_per_method.tsv`

**Sanity check the new HSM output:**

```bash
# New columns should be present
head -1 candidate_breakpoints_consensus.tsv | tr '\t' '\n' | grep -E "mode_(hsm|kde)|mode_agreement|mode_method"
# Should print: mode_hsm_bp, mode_kde_bp, mode_agreement_kb, mode_method
```

```bash
# For LG28, predicted outcome:
# - mode_hsm_bp within 5-15 kb of 15115000 (left) and 18005000 (right)
# - mode_agreement_kb < 20 kb typically
# - mode_method = "hsm"
awk -F'\t' '$1 ~ /LG28/ {print $0}' candidate_breakpoints_consensus.tsv
```

If mode_agreement_kb > 20 kb, that's a signal the KDE was substantially
oversmoothed on this distribution — keep the HSM value and note the
disagreement in the diagnostic PDF.

---

## Step 3 — Install the 5 schemas

All schemas go into `registries/schemas/structured_block_schemas/`:

```bash
cd inversion-popgen-toolkit/
REG_SCHEMAS=registries/schemas/structured_block_schemas/

cp phase4_chat_final/schemas/mechanism_assembled.schema.json   $REG_SCHEMAS
cp phase4_chat_final/schemas/bnd_sided_support.schema.json     $REG_SCHEMAS
cp phase4_chat_final/schemas/flank_coherence.schema.json       $REG_SCHEMAS
cp phase4_chat_final/schemas/synteny_v6.schema.json            $REG_SCHEMAS
cp phase4_chat_final/schemas/fragment_distribution.schema.json $REG_SCHEMAS
```

Verify all 5 are registered by the registry:

```bash
python3 -c "
import sys, json
sys.path.insert(0, 'registries/api/python')
from registry_loader import load_registry
reg = load_registry()
want = ['mechanism_assembled', 'bnd_sided_support', 'flank_coherence',
        'synteny_v6', 'fragment_distribution']
for bt in want:
    try:
        schema = reg.schemas.get(bt)
        print(f'OK   {bt}: {len(schema.get(\"properties\", {}))} properties')
    except Exception as e:
        print(f'FAIL {bt}: {e}')
"
```

All five should print OK. If any fails, check that the
`block_type` field inside the JSON matches the filename.

### 3a. Deprecate superseded schemas (optional)

The following schemas are superseded but kept for backwards compat:

- `synteny_dollo.schema.json` → replaced by `synteny_v6.schema.json`.
  Edit the JSON to add a `"deprecated": true` key and a `"replaced_by":
  "synteny_v6"` key. Don't delete it — existing blocks written under
  this schema remain valid.

- `interior_structure.schema.json` (if it was ever deployed from v6)
  → replaced by `fragment_distribution.schema.json`. Same treatment.

---

## Step 4 — Install the 6 scripts

Placement follows the role-based directory layout (see
`docs/toolkit_audit.md` — "don't reorganise the 4a/4b/4c/4d/4e
structure").

```bash
P4=inversion_modules/phase_4_postprocessing

# 4d_group_dependent/ — scripts that consume group assignments
cp phase4_chat_final/scripts/cheat29b_assembled_junction.py    $P4/4d_group_dependent/
cp phase4_chat_final/scripts/bnd_sided_support.py              $P4/4d_group_dependent/
cp phase4_chat_final/scripts/cross_species_bridge_v6.py        $P4/4d_group_dependent/
cp phase4_chat_final/scripts/bp_pipeline_bridge.py             $P4/4d_group_dependent/

# 4e_final_classification/ — final label synthesis
cp phase4_chat_final/scripts/assign_structural_class_v7.py     $P4/4e_final_classification/
```

Python syntax check:

```bash
for f in $P4/4d_group_dependent/cheat29b_assembled_junction.py \
         $P4/4d_group_dependent/bnd_sided_support.py \
         $P4/4d_group_dependent/cross_species_bridge_v6.py \
         $P4/4d_group_dependent/bp_pipeline_bridge.py \
         $P4/4e_final_classification/assign_structural_class_v7.py; do
  python3 -c "import ast; ast.parse(open('$f').read())" && echo "OK $f"
done
```

All five should print OK.

### 4a. Remove superseded scripts

If any of these are already in the repo from prior sessions, delete them:

```bash
rm -f $P4/4d_group_dependent/cross_species_bridge.py    # v5, replaced by v6
rm -f $P4/4e_final_classification/assign_structural_class_v5.py  # v5
rm -f $P4/4e_final_classification/assign_structural_class_v6.py  # v6
rm -f $P4/4b_group_proposal/interior_structure_diagnostic.py    # v6 mistake

# Also the superseded schemas if they were deployed:
# (confirm first with git log before deleting)
```

---

## Step 5 — Install the launcher

The v5 launcher orchestrates the new-evidence layer. It lives in
`launchers/` at the toolkit root (or in the 4d directory — follow
whatever convention the repo uses for launchers; check a few existing
ones).

```bash
cp phase4_chat_final/launchers/run_phase4_v5_new_evidence.sh \
   inversion_modules/phase_4_postprocessing/4d_group_dependent/

chmod +x inversion_modules/phase_4_postprocessing/4d_group_dependent/run_phase4_v5_new_evidence.sh
```

Inspect it and update any hardcoded paths to match the toolkit's
environment-variable conventions (BASE, PROJECT_DIR, etc.). The
current launcher uses positional args; if your repo pattern is
`source 00_module_config.sh` first, wrap accordingly.

### 5a. Create `run_phase4_v7_wiring.sh` — the master orchestrator

Create a new master launcher that calls everything in order. Suggested
content:

```bash
#!/bin/bash
#SBATCH --job-name=phase4_v7_wiring
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/phase4_v7_%j.out

set -euo pipefail

# Source project config
source /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit/00_inversion_config.sh

# Paths
P4=${BASE}/inversion-popgen-toolkit/inversion_modules/phase_4_postprocessing
BP_OUT=${BASE}/breakpoint_pipeline_output
V7_OUT=${BASE}/phase4_v7_blocks
REG=${BASE}/inversion-popgen-toolkit/registries

mkdir -p ${V7_OUT} logs

# 1. Run v5 new-evidence layer (cheat29b + bnd_sided_support +
#    original cross_species_bridge has been superseded, skip)
bash ${P4}/4d_group_dependent/run_phase4_v5_new_evidence.sh

# 2. Run v6 cross-species bridge (Kuang scheme)
python3 ${P4}/4d_group_dependent/cross_species_bridge_v6.py \
  --candidates ${BASE}/CANDIDATE_TABLE.tsv \
  --between_species_bed ${BASE}/catfish-synteny-toolkit/output/step_02/events.bed \
  --flank_coherence_tsv ${BASE}/catfish-synteny-toolkit/output/step_09c/flank_coherence.tsv \
  --polarized_tsv ${BASE}/catfish-synteny-toolkit/output/step_11/polarized.tsv \
  --dollo_tsv ${BASE}/dollo_output.tsv \
  --outdir ${V7_OUT} \
  --registries_root ${REG}

# 3. Run bp_pipeline bridge (reads breakpoint_pipeline TSVs)
python3 ${P4}/4d_group_dependent/bp_pipeline_bridge.py \
  --consensus_tsv ${BP_OUT}/candidate_breakpoints_consensus.tsv \
  --fragments_tsv ${BP_OUT}/candidate_ancestral_fragments.tsv \
  --per_method_tsv ${BP_OUT}/candidate_breakpoints_per_method.tsv \
  --outdir ${V7_OUT} \
  --registries_root ${REG}

# 4. Run final synthesis
python3 ${P4}/4e_final_classification/assign_structural_class_v7.py \
  --candidates ${BASE}/CANDIDATE_TABLE.tsv \
  --keys_dir ${V7_OUT} \
  --v5_blocks_dir ${V7_OUT} \
  --outdir ${V7_OUT}/final

echo "phase4 v7 complete; final catalog at ${V7_OUT}/final/final_catalog.tsv"
```

Save as `${P4}/4e_final_classification/run_phase4_v7_wiring.sh` and
`chmod +x`.

Adjust paths to your actual project layout — the templates above
assume the user's documented paths from SESSION_AUDIT.

---

## Step 6 — Wire into `compute_candidate_status.R`

This is where v7's final label becomes part of the Axis system.

### 6a. Read `final_label.json` as Axis 5

Edit `${P4}/4e_final_classification/compute_candidate_status.R`:

```r
# Add a helper function
read_final_label <- function(cid, v7_blocks_dir) {
  p <- file.path(v7_blocks_dir, "final", cid, "final_label.json")
  if (!file.exists(p)) {
    return(list(
      q_overall_structural_class = NA_character_,
      weakest_component = NA_character_,
      justification = NA_character_
    ))
  }
  j <- jsonlite::fromJSON(p, simplifyVector = TRUE)
  list(
    q_overall_structural_class = j$q_overall_structural_class,
    weakest_component = j$weakest_component,
    justification = paste(j$justification, collapse = "; ")
  )
}

# In the axis assembly section, add:
axis_5 <- read_final_label(cid, V7_BLOCKS_DIR)
```

Then write axis_5 into the per-candidate output alongside axes 1-4.

### 6b. Add new wired keys to `build_key_spec()`

The Spec vs Reality audit (`SPEC_VS_REALITY.md`) tracks which keys
are wired vs aspirational. The new blocks from this bundle contribute:

| Block | New wired keys | Count |
|---|---|---|
| `fragment_distribution` | q2_interior_class, q2_fragment_interior_recomb_fraction, q2_n_fragment_carriers, q3_final_left_bp, q3_final_right_bp, q3_left_ci_width_kb, q3_right_ci_width_kb, q3_breakpoint_precision_class, q3_n_methods_agreeing_left, q3_n_methods_agreeing_right | 10 |
| `synteny_v6` | q5_bs_event_overlap, q5_bs_event_type, q5_tree_polarization_direction, q5_tree_polarization_confidence, q5_dollo_vs_tree_concordance, q5_conservation_class | 6 |
| `mechanism_assembled` | q4b_asm_precise_record_available, q4b_asm_junction_class, q4b_asm_homlen, q4b_asm_source | 4 |
| `bnd_sided_support` | q7b_bnd_single_sided_left, q7b_bnd_single_sided_right, q7b_bnd_support_score | 3 |

In `build_key_spec()`, move these keys from `*_aspir` to `*_wired`.

Expected spec update:
- Q1 49/49 (unchanged)
- Q2 55+3/55 (fragment keys are Q2)
- Q3 65+5/74 (consensus keys are Q3) → 70/74
- Q4 28+4/47 → 32/47
- Q5 28+6/39 → 34/39
- Q6 12/28 (unchanged)
- Q7 34+3/75 → 37/75

Total: 271+21 = **292/367 wired** (up from 271/367 = 73.8% → 79.6%).

Update `SPEC_VS_REALITY.md` summary table accordingly.

---

## Step 7 — Wire `phase_7_cargo` consumers

The BREEDING_A/C/D scripts need `final_label.json` to work.

### 7a. Add final_label reader to `phase_7_cargo/compute/_lib_group_carriership.R`

Add a helper (or in whatever lib file is the shared helpers for these
scripts):

```r
load_final_labels <- function(v7_final_dir) {
  label_files <- list.files(v7_final_dir,
                             pattern = "final_label.json$",
                             recursive = TRUE, full.names = TRUE)
  rbindlist(lapply(label_files, function(p) {
    j <- jsonlite::fromJSON(p, simplifyVector = TRUE)
    data.table(
      candidate_id = j$candidate_id,
      q_overall_structural_class = j$q_overall_structural_class,
      weakest_component = j$weakest_component
    )
  }))
}
```

### 7b. Filter breeding inputs by label

In `BREEDING_A_broodstock_compatibility.R`, at the start:

```r
labels <- load_final_labels(V7_FINAL_DIR)

# Exclude out-of-scope and family-confounded candidates from breeding logic
valid_labels <- c(
  "supported_balanced_inversion_simple",
  "supported_balanced_inversion_NAHR_like_supported_by_assembly",
  "supported_balanced_inversion_NAHR_like_hypothesis",
  "supported_balanced_inversion_NHEJ_like_supported_by_assembly",
  "supported_balanced_inversion_NHEJ_like_hypothesis",
  "supported_balanced_inversion_with_substrate_mechanism_unresolved",
  # Allow suffixed variants
  "supported_balanced_inversion_with_edge_recombinants",
  "supported_shared_between_species_inversion",
  "supported_nested_inversion"
)
# Pattern match on "supported_balanced_inversion_" prefix also OK
usable_candidates <- labels[
  grepl("^supported_(balanced_inversion|shared|nested)", q_overall_structural_class)
]$candidate_id
```

Same pattern in BREEDING_C and BREEDING_D.

### 7c. Document in phase_7_cargo/README.md

Add a section "Inputs from Phase 4":

> BREEDING_A/C/D require the Phase 4 final labels. Ensure
> `run_phase4_v7_wiring.sh` has run and produced
> `phase4_v7_blocks/final/final_catalog.tsv` before running these.
> Labels are loaded from `V7_FINAL_DIR` (set via env var
> `V7_FINAL_DIR=${BASE}/phase4_v7_blocks/final` before sourcing
> the breeding config).

---

## Step 8 — Close the 4 wiring gaps from `toolkit_audit.md`

From §D.2 of `docs/toolkit_audit.md`, four small wiring gaps exist
beyond this bundle. Close them in the same wiring PR:

### Gap 1 — Q7B keys from `breakpoint_evidence_audit.py`

In `build_key_spec()`:

```r
q7_wired <- c(q7_wired,
  "q7b_pca_carriers_strong_sv",
  "q7b_pca_carriers_weak_sv",
  "q7b_pca_carriers_no_sv",
  "q7b_sv_only_carriers",
  "q7b_audit_status",
  "q7b_audit_concordance",
  "q7b_audit_score",
  "q7b_audit_coverage_fraction"
)
```

Remove corresponding entries from `q7_aspir`. ~10 lines.

### Gap 2 — cheat30 GDS keys into characterize_q5

In `characterize_q5()`:

```r
# Before computing age, check GDS consistency
gds_ok <- keys$q5_hom_inv_gds_clean == "clean"
if (isTRUE(gds_ok) && !is.na(age_point_estimate)) {
  age_class <- "answered"
} else {
  age_class <- "measured"
}
```

### Gap 3 — Axis 5 (done in Step 6a above)

### Gap 4 — gene_conversion_detector into characterize_q2

In `characterize_q2()`:

```r
gc <- load_block(cand_dir, "gene_conversion_tracts")
if (!is.null(gc)) {
  keys$q2_gc_tract_count <- gc$n_tracts
  keys$q2_gc_tract_total_bp <- gc$total_bp
}
```

Add to `q2_wired`:
```r
q2_wired <- c(q2_wired, "q2_gc_tract_count", "q2_gc_tract_total_bp")
```

---

## Step 9 — End-to-end validation run

Run on the full 226-sample cohort and verify:

```bash
sbatch ${P4}/4e_final_classification/run_phase4_v7_wiring.sh
# Wait for completion (~1-4 hours depending on candidate count)

# Check the final catalog
cat ${V7_OUT}/final/final_catalog.tsv | head

# Distribution of labels
awk -F'\t' 'NR>1 {print $2}' ${V7_OUT}/final/final_catalog.tsv | sort | uniq -c | sort -rn
```

Expected qualitative result:
- Most candidates: `supported_balanced_inversion_*` variants
- Small fraction: `complex_rearrangement_out_of_scope` (the gate working)
- Small fraction: `family_confounded_locus`
- Some: `supported_shared_between_species_inversion` (Kuang scheme working)

If `complex_rearrangement_out_of_scope` is >50% of candidates, the
recombinant-fraction threshold is too aggressive — drop from 0.70 to
0.80 in `assign_structural_class_v7.py::interior_class` rule 1.

---

## Step 10 — Commit and tag

```bash
cd inversion-popgen-toolkit/
git add -A
git status  # review before committing

git commit -m "Phase 4 v7 wiring: bridge to breakpoint_pipeline + HSM mode estimator

- Install HSM patch for fragment-boundary mode estimation
  (Bickel & Frühwirth 2006; replaces Silverman-bandwidth KDE)
- Add bp_pipeline_bridge.py — derives Phase 4 keys from
  breakpoint_pipeline outputs
- Add cross_species_bridge_v6.py — Kuang 2026 style annotation
  with Dollo as cross-check only
- Add v7 assigner with complex_rearrangement_out_of_scope gate
- Add 5 new schemas, 21 new wired keys, 1 new axis
- Close 4 wiring gaps from toolkit_audit.md

Spec coverage: 271/367 → 292/367 wired (73.8% → 79.6%)
"

git tag -a v7.0 -m "Phase 4 v7 wiring complete"
```

---

## What this session must NOT do

Stay out of these traps (from session audit §D19 + §D22 + §D25):

1. **Don't write new analysis modules.** The bundle is a wiring
   package. Everything novel has been done in prior sessions. Only
   mechanical connection work.

2. **Don't merge/refactor yet.** The toolkit_audit.md lists 4 merge
   candidates and 3 simplifications. They are not urgent. Ship wiring
   first; do the clean-up PR in a separate session.

3. **Don't improve DUP-TRP / INV-DUP gating.** The gate threshold
   of 0.70 is a first pass. If it misbehaves on real data, ADJUST
   the threshold, don't replace the approach.

4. **Don't replace HSM with something else mid-session.** If HSM
   disagrees with KDE by a lot on real LG28, that's the DIAGNOSTIC
   working (mode_agreement_kb is supposed to flag that). HSM is
   the new primary; KDE is the secondary. Don't switch back.

5. **Don't make LG28 the flagship.** User was explicit that LG28 is
   a test case. Use the full-cohort v7 output to identify a cleaner
   candidate (`q2_interior_class == "clean_simple"` AND
   `q5_conservation_class != "unclassified_complex"`) for the
   manuscript centrepiece.

---

## Quick-reference glossary of where everything lives after wiring

| File | Path | Consumed by |
|---|---|---|
| `02_ancestral_fragments.R` (patched) | `breakpoint_pipeline/` | `03_consensus_merge.R`; produces 3 TSVs |
| Three breakpoint TSVs | wherever run_pipeline.sh writes | `bp_pipeline_bridge.py` |
| `cheat29b_assembled_junction.py` | `phase_4/4d_group_dependent/` | `assign_structural_class_v7.py` via `mechanism_assembled.json` |
| `bnd_sided_support.py` | `phase_4/4d_group_dependent/` | `assign_structural_class_v7.py` via `bnd_sided_support.json` |
| `cross_species_bridge_v6.py` | `phase_4/4d_group_dependent/` | `assign_structural_class_v7.py` via `synteny_v6.json` |
| `bp_pipeline_bridge.py` | `phase_4/4d_group_dependent/` | `assign_structural_class_v7.py` via `fragment_distribution.json` |
| `assign_structural_class_v7.py` | `phase_4/4e_final_classification/` | `compute_candidate_status.R` via `final_label.json` |
| `final_label.json` (per candidate) | `phase4_v7_blocks/final/<cid>/` | `compute_candidate_status.R` Axis 5; phase_7_cargo BREEDING_A/C/D |

---

## Timeline estimate

- Step 1 (HSM patch): 20 min
- Step 2 (run breakpoint_pipeline on LG28): 30-60 min compute
- Step 3 (schemas): 10 min
- Step 4 (scripts): 15 min
- Step 5 (launchers): 30 min
- Step 6 (compute_candidate_status wiring): 30 min
- Step 7 (phase_7_cargo wiring): 45 min
- Step 8 (4 wiring gaps): 60 min
- Step 9 (end-to-end run): 2-4 hours compute + 30 min review
- Step 10 (commit + tag): 10 min

**Total active work: ~4 hours + 3-5 hours compute.** One focused day.

Do not deviate from this plan. If anything in the pipeline surprises
you, document it in a new SESSION_AUDIT entry and continue — do not
rewrite mid-session.
