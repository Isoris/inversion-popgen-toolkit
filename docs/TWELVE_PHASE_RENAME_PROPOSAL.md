# Design proposal: 12-phase rename (Option 1 from chat-3-end discussion)

**Status:** design doc only — no moves in this pass. Execute in a fresh chat.
**Originated:** 2026-04-24, end of chat 3 (passes 9, 10, 11-partial, 12, 13 shipped)
**Supersedes:** `docs/PHASE4_RENUMBER_PROPOSAL.md` (pass 12's internal phase_4 renumber — that proposal only restructured *inside* phase_4; this one breaks phase_4 open entirely)

---

## Problem statement

`phase_4_postprocessing/` currently holds 7 sub-blocks (4a–4g) doing at least three orthogonal things:

| Current sub-block | What it produces | Is it "postprocessing"? |
|---|---|---|
| `4a_existence_layers/` | The candidate catalog (C01d scoring, C01g boundary unification) | Catalog construction |
| `4b_qc_triage/` | Data-quality flag per candidate (clean/messy) | Triage, not postprocessing |
| `4c_breakpoint_refinement/` | bp-resolution breakpoints + CI via ancestral fragments | **Real analysis** — manuscript figures |
| `4d_group_proposal/` | Karyotype-group assignments (Hom1/Het/Hom2) | Analysis |
| `4e_group_validation/` | Hypothesis tests + VALIDATED gate | Analysis |
| `4f_group_dependent/` | Mechanism (NAHR/NHEJ/MMBIR), age, cross-species conservation | **Real evolutionary biology** |
| `4g_final_classification/` | Tier assignment + 7-question characterization | Final classification |

Quentin's concern, verbatim:

> "the biology is not only post process, its also real analysis. so idk i know its strange but find a way bc now its nesting too much"

And:

> "phase 4 post processing is a strange name its too general. maybe it must be flattened and have many more phases?"

The folder name `phase_4_postprocessing` misrepresents what lives there. A reader of the manuscript's Methods section who greps the repo will find the mechanism / age / cross-species work under a folder labeled `postprocessing`, which undersells what it is.

---

## Proposed layout

Flatten `phase_4` into 6 dedicated phases, shift downstream phases accordingly:

```
inversion_modules/
├── phase_1_inputs/                     (unchanged)
├── phase_2_discovery/                  (unchanged)
├── phase_3_refine/                     (unchanged — SV caller evidence → registry)
│
├── phase_4_catalog/                    ← was phase_4/4a_existence_layers/
│   PRODUCES: the candidate catalog
│   CONTAINS: C01d scoring, C01e blocks, C01g boundary unification
│
├── phase_5_qc_triage/                  ← was phase_4/4b_qc_triage/
│   PRODUCES: data-quality flag per candidate
│   CONTAINS: Q01–Q10 shelf QC, hobs per-group, bridge_from_phase4
│
├── phase_6_breakpoint_refinement/      ← was phase_4/4c_breakpoint_refinement/
│   PRODUCES: bp-resolution breakpoints + CI
│   CONTAINS: dosage signal, ancestral fragments, consensus merge
│
├── phase_7_karyotype_groups/           ← was phase_4/4d + 4e merged
│   PRODUCES: per-candidate Hom1/Het/Hom2 assignments with validation level
│   CONTAINS: 4d_group_proposal/ + 4e_group_validation/ as sub-blocks
│   (OR: keep them as separate phase_7 + phase_8 if merge feels wrong)
│
├── phase_8_evidence_biology/           ← was phase_4/4f_group_dependent/
│   PRODUCES: mechanism + age + cross-species + SV audit keys per candidate
│   CONTAINS: cheats, bridges, breakpoint_evidence_audit, etc.
│   (4f sub-organization — see section below)
│
├── phase_9_classification/             ← was phase_4/4g_final_classification/
│   PRODUCES: final tier + 7-question characterization
│   CONTAINS: compute_candidate_status.R, characterize_candidate.R,
│             assign_structural_class_v7.py, _axis5_final_label.R,
│             _qc_shelf_reader.R (from pass 13)
│
├── phase_10_followup/                  ← was phase_5_followup
│   (per-candidate deep analysis / dosage rasters)
│
├── phase_11_secondary/                 ← was phase_6_secondary
│   (LD + Fst secondary confirmations)
│
├── phase_12_cargo/                     ← was phase_7_cargo
│   (gene content + breeding implications)
│
└── utils/                              (unchanged)
```

**Net effect:** 7 phases → 12 phases. Each phase has one clear purpose in its name. The word "postprocessing" disappears.

### Open question — phase_7_karyotype_groups merge vs split

4d_group_proposal and 4e_group_validation are currently separate because
the group-validation gate (SUPPORTED / VALIDATED / UNCERTAIN / SUSPECT)
is what guards downstream analyses. Merging them into one `phase_7_karyotype_groups/`
folder with `proposal/` and `validation/` subfolders keeps the gate
visible while reducing the top-level phase count.

Alternative: keep them as `phase_7_group_proposal/` + `phase_8_group_validation/`,
then shift everything else by one (evidence_biology → 9, classification → 10,
followup → 11, secondary → 12, cargo → 13). **13 phases feels like a lot.** Recommend the merge.

### phase_8_evidence_biology internal organization

The current flat 4f_group_dependent/ (14 files) should be regrouped by question:

```
phase_8_evidence_biology/
├── README.md
├── _archive/
├── q4_mechanism/
│   ├── cheat27_sd_nahr_substrate.R
│   ├── cheat28_tandem_repeat_context.R
│   ├── cheat29_junction_forensics.R
│   └── cheat29b_assembled_junction.py
├── q5_age_and_origin/
│   ├── cheat30_gds_by_genotype.R
│   └── cross_species_bridge_v6.py
├── q7_existence_audit/
│   ├── cheat6_ancestry_jackknife_v934.R
│   ├── bnd_sided_support.py
│   └── breakpoint_evidence_audit.py
├── bp_bridge/
│   └── bp_pipeline_bridge.py       (Q2+Q3 — reads from phase_6)
└── run_phase4_v5_new_evidence.sh   (orchestrator, rename to run_evidence_biology.sh)
```

**Also pull out of phase_8_evidence_biology entirely:**
- `population_regenotype.py` → `Modules/MODULE_4D_delly_inv/utils/`
- `SLURM_A03b_population_regenotype.sh` → `Modules/MODULE_4D_delly_inv/slurm/`

(reason: the author's declared pipeline slot is "STEP_A03b after bcftools merge in MODULE_4D/4E", not phase 4; it's SV-caller infrastructure misplaced in 4f. The SLURM_A03b launcher already has a fallback lookup at `../utils/population_regenotype.py` showing the original intent.)

---

## Blast radius

Enumerated against the current repo state on 2026-04-24:

| Category | File count | Notes |
|---|---|---|
| Shell scripts (.sh, .slurm) | 77 | Orchestrators, LAUNCHers, SLURM arrays |
| R scripts (.R) | 22 | Source() refs, helper paths |
| Python (.py) | 7 | Import paths, argparse defaults |
| Markdown (.md) | 50 | READMEs, docs, session audits |
| JSON (.json) | 3 | Registry schemas with `source_script` fields |
| **TOTAL unique files** | **159** | |

By top-level directory:
- 52 files in `inversion_modules/phase_4_postprocessing/` (internal)
- 59 files in `Modules/MODULE_4*` (external — DELLY/Manta launchers that source phase_4 launchers)
- 7 files in `inversion_modules/phase_7_cargo/`
- 11 files in `docs/`
- Remaining: scattered in tools/, utils/, registries/, other phases

**Comparison to previous passes:**
- Pass 9 (phase_3 Layer rename): 18 files
- Pass 12 (phase_4 sub-block renumber + fold in qc_shelf + breakpoint_pipeline): 66 files
- **Pass 15 (this proposal): ~159 files** — about 2.4× pass 12

---

## Execution plan

### Stage 1: checkpoint tarball

Capture the current state. If anything goes wrong, we restore from this.

### Stage 2: directory moves (in correct order to avoid collisions)

**Inside phase_4_postprocessing** (6 moves from sub-blocks to siblings):

```bash
cd inversion_modules/phase_4_postprocessing

# Move sub-blocks OUT to become top-level phases. Do this BEFORE shifting
# phase_5/6/7 to avoid name collisions at every step.
mv 4a_existence_layers        ../phase_4_catalog
mv 4b_qc_triage               ../phase_5_qc_triage
mv 4c_breakpoint_refinement   ../phase_6_breakpoint_refinement
# 4d + 4e merge path: create phase_7 folder, move both in
cd ..
mkdir -p phase_7_karyotype_groups
mv phase_4_postprocessing/4d_group_proposal    phase_7_karyotype_groups/proposal
mv phase_4_postprocessing/4e_group_validation  phase_7_karyotype_groups/validation
mv phase_4_postprocessing/4f_group_dependent   phase_8_evidence_biology
mv phase_4_postprocessing/4g_final_classification  phase_9_classification

# At this point phase_4_postprocessing/ holds only README.md, docs/,
# orchestrator/, patches/, schemas/, specs/, tests/.
# These need distribution:
#   orchestrator/ → keep at phase_9_classification/orchestrator/ (it orchestrates the chain)
#   schemas/      → move to registries/schemas/ (they're cross-phase)
#   specs/        → move to registries/schemas/ (same — they describe registry blocks)
#   tests/        → distribute per-phase OR keep as phase_9_classification/tests/
#   patches/      → keep at phase_9_classification/patches/
#   docs/         → distribute per-phase based on content
#   README.md     → delete (replaced by per-phase READMEs)
```

**Then shift downstream phases** (work from highest-numbered to lowest to avoid collisions):

```bash
# 7_cargo -> 12_cargo
mv phase_7_cargo      phase_12_cargo
# 6_secondary -> 11_secondary
mv phase_6_secondary  phase_11_secondary
# 5_followup -> 10_followup
mv phase_5_followup   phase_10_followup

# Finally delete the now-empty phase_4_postprocessing shell
rmdir phase_4_postprocessing   # after verifying nothing remains
```

### Stage 3: evidence_biology internal regroup + regenotype move

```bash
cd inversion_modules/phase_8_evidence_biology

mkdir q4_mechanism q5_age_and_origin q7_existence_audit bp_bridge
mv cheat27_* cheat28_* cheat29_* cheat29b_* q4_mechanism/
mv cheat30_* cross_species_bridge_v6.py     q5_age_and_origin/
mv cheat6_*  bnd_sided_support.py breakpoint_evidence_audit.py  q7_existence_audit/
mv bp_pipeline_bridge.py                     bp_bridge/
mv run_phase4_v5_new_evidence.sh             run_evidence_biology.sh

# Move SV-caller regenotyping out of phase 4 entirely
mv population_regenotype.py                  ../../Modules/MODULE_4D_delly_inv/utils/
mv SLURM_A03b_population_regenotype.sh       ../../Modules/MODULE_4D_delly_inv/slurm/
```

### Stage 4: sed pass — 5 rename tiers in strict order

Each tier must complete entirely before the next one starts, otherwise
source and destination names collide.

```bash
# Build FILE LIST from scratch BEFORE each tier, because earlier tiers
# modify which files contain which old strings.

# Tier 1: downstream phases FIRST (high-to-low, no collisions)
#   phase_7_cargo   -> phase_12_cargo
#   phase_6_secondary -> phase_11_secondary
#   phase_5_followup  -> phase_10_followup
FILES=$(grep -rlE 'phase_(5_followup|6_secondary|7_cargo)' --include='*.sh' --include='*.R' --include='*.py' --include='*.md' --include='*.json' --include='*.slurm' . | grep -v _archive | grep -v .git/)
for f in $FILES; do
  sed -i 's|phase_7_cargo|phase_12_cargo|g;
          s|phase_6_secondary|phase_11_secondary|g;
          s|phase_5_followup|phase_10_followup|g' "$f"
done

# Tier 2: phase_4 sub-blocks -> new top-level phases
FILES=$(grep -rlE 'phase_4_postprocessing/4[a-g]_' --include='*.sh' --include='*.R' --include='*.py' --include='*.md' --include='*.json' --include='*.slurm' . | grep -v _archive | grep -v .git/)
for f in $FILES; do
  sed -i 's|phase_4_postprocessing/4g_final_classification|phase_9_classification|g;
          s|phase_4_postprocessing/4f_group_dependent|phase_8_evidence_biology|g;
          s|phase_4_postprocessing/4e_group_validation|phase_7_karyotype_groups/validation|g;
          s|phase_4_postprocessing/4d_group_proposal|phase_7_karyotype_groups/proposal|g;
          s|phase_4_postprocessing/4c_breakpoint_refinement|phase_6_breakpoint_refinement|g;
          s|phase_4_postprocessing/4b_qc_triage|phase_5_qc_triage|g;
          s|phase_4_postprocessing/4a_existence_layers|phase_4_catalog|g' "$f"
done

# Tier 3: bare phase_4_postprocessing refs that don't have a sub-block
# (these now point nowhere — most are in docs; a few may be schemas/
# or specs/ moves)
FILES=$(grep -rlE 'phase_4_postprocessing' --include='*.sh' --include='*.R' --include='*.py' --include='*.md' --include='*.json' --include='*.slurm' . | grep -v _archive | grep -v .git/)
# Manual inspection required for each — no safe mass sed here.

# Tier 4: bare 4[a-g]_ sub-block refs without parent path
# (e.g., "../../4d_group_proposal/lib_decompose_helpers.R" from inside
# phase_7_cargo). These now need new relative paths.
# Manual inspection required.

# Tier 5: the regenotype move
FILES=$(grep -rlE 'population_regenotype|SLURM_A03b_population_regenotype' --include='*.sh' --include='*.R' --include='*.py' --include='*.md' . | grep -v _archive | grep -v .git/)
for f in $FILES; do
  sed -i 's|phase_4_postprocessing/4f_group_dependent/population_regenotype|Modules/MODULE_4D_delly_inv/utils/population_regenotype|g;
          s|phase_4_postprocessing/4f_group_dependent/SLURM_A03b_population_regenotype|Modules/MODULE_4D_delly_inv/slurm/SLURM_A03b_population_regenotype|g' "$f"
done

# Tier 6: the evidence_biology internal subfolder refs
# Orchestrator LAUNCH_group_cheats.sh has ${PHASE4D_DIR}/cheat27_*, etc.
# Need updating to ${PHASE8_DIR}/q4_mechanism/cheat27_*, etc.
# Manual inspection required — not a simple sed.
```

### Stage 5: doc + README rewrites

- **New README for each of the 12 phases** — each phase gets a dedicated README explaining what it produces and what it consumes.
- **`inversion_modules/README.md`** — tree diagram updated to show all 12 phases
- **`docs/MODULE_MAP.md`** — phase table rewritten
- **`phase_4_postprocessing/README.md`** — deleted (replaced by phase-specific READMEs)

### Stage 6: verification pass

For every one of the 159 files:
1. `bash -n` on shell scripts
2. `python3 -m ast.parse` on Python
3. R parse check on R scripts
4. JSON parse on schemas
5. Markdown fence balance
6. `grep -rn` for any remaining old strings

Plus functional tests of the pass-8 axis-5 wiring, pass-13 q_qc_shelf
reader, and pass-11 bridge script — all three use env vars
(`V7_FINAL_DIR`, `QC_SHELF_EVIDENCE_DIR`) or `$here`-relative paths, so
should survive the folder moves — but verify.

---

## Risks

1. **HPC tree divergence during transition.** The laptop rename happens in
   one commit; HPC needs to pull (or rsync) to match. Jobs running on HPC
   with old paths in sbatch scripts will fail silently. Coordinate the push.

2. **Output directories under `FOLLOWUP_DIR`, `CATALOG_DIR`, etc.** These are
   not under `MODULE_DIR`, so they're safe. But if any launcher uses
   `${MODULE_DIR}/phase_4_postprocessing/<something>/output/`, that's broken.
   Verify with `grep MODULE_DIR` on the 77 shell scripts.

3. **Registry `source_script` fields.** 3 JSON schemas have
   `source_script: "phase_4_postprocessing/..."`. These are metadata strings,
   not paths — updating them is for human-readability, not functionality.

4. **Registry `method=` tags.** `method="phase_qc_shelf"` in existing
   interval_registry rows stays valid — it's a string tag, not a path.
   Do NOT touch it. (Same principle as pass 12's Q10 METHOD_TAG preservation.)

5. **Pass 8's axis 5, pass 13's q_qc_shelf reader, pass 11's bridge.** All
   three pieces are env-gated (`V7_FINAL_DIR`, `QC_SHELF_EVIDENCE_DIR`) or
   use `$here`-relative resolution. Folder moves shouldn't break them —
   verify in stage 6.

6. **The pass 12 DROP_README still talks about "4b_qc_triage sits inside
   phase_4_postprocessing."** Needs updating or archiving as historical.

7. **Claude context limit.** Pass 12 was 66 files and took most of a chat.
   This pass is 159 files. Plan for 2 chats, ideally:
   - Chat A: stages 1–4 (moves + sed)
   - Chat B: stages 5–6 (doc rewrites + verification)

---

## Don't-do list

Things NOT to do in this pass:

- Don't attempt the 4f biology subfolder regroup AND the phase renumber in
  the same sed pass — do the phase moves first, then regroup.
- Don't re-tag registries (`METHOD_TAG` stays `phase_qc_shelf` forever —
  changing it orphans existing rows).
- Don't rename the registry block_types (`existence_layer_d`, `age_evidence`,
  etc.) — these are stable public API.
- Don't rename pass 13's `q_qc_shelf_*` keys (same reason).
- Don't renumber phases 1/2/3 — they're stable and haven't changed.

---

## Suggested opening prompt for fresh chat

```
I'm continuing from chat-3-end (passes 9, 10, 11-partial, 12, 13
shipped). Option 1 (12-phase rename) was chosen. Read
docs/TWELVE_PHASE_RENAME_PROPOSAL.md for the full plan. Start with
stage 1 (checkpoint tarball), then stages 2–3 (directory moves).
Verify each stage before proceeding.

Cohort identity reminder: 226 pure C. gariepinus hatchery broodstock.
Quentin Andres is the user. Registry method tags stay as
"phase_qc_shelf" — string tags, not paths.
```

---

## Alternative: don't do this

Option 4 from the same discussion was "keep folders, fix the README to
explain what lives where." Zero code changes, only a clarifying table
at the top of `phase_4_postprocessing/README.md`. That option is still
available and may genuinely be right if:

- Manuscript submission deadline is near (this pass is 2 chats of time
  that don't produce new biology)
- HPC coordination is expensive right now
- You'd rather spend the time on OPEN_TODOS items (T2/T3/T7/T12) or
  on threshold tuning for the q_qc_shelf reader

The proposal above is the **structural right answer**. Whether it's the
**right answer for right now** is a scheduling question.

---

## Checkpoint: passes already shipped (do not redo)

- ✅ Pass 9: phase_3_refine Layer-tagged rename (STEP_{A,B,D}NN_)
- ✅ Pass 10: phase_6 MODULE_5E archive
- 🟡 Pass 11 (partial): phase_qc_shelf → 4b_qc_triage bridge script
- ✅ Pass 12: phase_4 restructured from 5 sub-blocks (4a–4e) to 7 (4a–4g);
   phase_qc_shelf folded in as 4b, breakpoint_pipeline folded in as 4c
- ✅ Pass 13: q_qc_shelf_* reader wired into 4g_final_classification
   (the file now lives at 4g_final_classification/_qc_shelf_reader.R,
   will move to phase_9_classification/_qc_shelf_reader.R after this pass)

---

## Pending tasks (not in this proposal)

From OPEN_TODOS + previous handoffs:

- Pass 11 task 2: rewrite `4b_qc_triage/README.md` with 2-mode structure
- Pass 11 task 3: polish `docs/MODULE_MAP.md`
- OPEN_TODOS T2/T3: BREEDING + Q7B wiring (needs prior-chat spec)
- OPEN_TODOS T7: RENAMING Cat 2 sed pass
- OPEN_TODOS T12: manuscript MAPQ inconsistency
- HPC-blocked: T5 full v7 run, T6 Engine H validation

These can interleave with this pass, but not inside the same chat.
