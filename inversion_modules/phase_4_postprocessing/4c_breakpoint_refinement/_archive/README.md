# _archive — other files from this session

These are files produced during this session that are NOT part of the core breakpoint pipeline. They're kept here so nothing is lost, but they don't fit into the `01 → 02 → 03 → 04` flow.

## What's here

### `heatmap_ordering_upgrade/` — the heatmap row-ordering fix

From the first half of the session. The problem: STEP28/STEP29/STEP35/q08 were sorting heatmap rows using `setorder(ri, coarse_group_refined, u)` with no within-karyotype tiebreaker, which collapsed sub-structure for composite candidates like LG28.

**The fix** was a shared helper `within_group_ordering.R` that does u-primary sorting with a configurable secondary axis (hclust, v, residual_pc1, error_profile, or none). STEP36 was upgraded to render 3-version (V1/V2/V3) comparison heatmaps so you can diagnose which ordering is right per candidate. STEP41 and q08b were upgraded to consume the helper.

**Status**: functional, brace-balance-checked, never runtime-tested. Ready to drop into your codebase.

**Naming collision warning**: the scripts in here are named `STEP36_…`, `STEP40_…`, `STEP41_…` — but your actual codebase has existing scripts with those numbers doing different things (Clair3 local signatures, internal coherence, membership trajectories). **These need to be renumbered before use** — likely `STEP42`, `STEP43`, `STEP44`. This is flagged in `docs/SESSION_AUDIT.md` section "Things the next session should double-check."

### `architecture_classifier/` — the classifier I was asked to stop building

From the first third of the session. An architecture-class classifier with 10 features and 7 classes (SIMPLE_STRONG, SIMPLE_WEAK, COMPOSITE_INTERNAL, COMPOSITE_OVERLAP, DEMOGRAPHIC_BLOCK, TECHNICAL, UNRESOLVED) and 4 tiers. It produces per-sample subgroup labels via a "separate first, correct second" partition.

**Why it's archived, not in the pipeline**:
1. You explicitly told me to stop building more classifiers partway through the session.
2. It has a critical-review note from me: the classifier is probably noob-looking-smart. Simple vs composite is maybe one real distinction; 7 classes and 4 tiers is granularity for granularity's sake. A two-class (simple / composite) heuristic probably does the same job.
3. It also has the name-collision problem (your STEP40 exists).

**Status**: brace-balance-checked. Includes a Python synthetic-data self-test (`test_classifier_selftest.py`) with 5 scenarios — 3 passed, 2 failed due to test-harness RNG issues, not classifier bugs.

**Recommendation**: leave it here. If you later decide you actually need an architecture classifier, start from scratch with a simpler two-class version. This one is a reference for what NOT to do.

### `original_workflow_diagrams/` — superseded diagrams

`WORKFLOW_DIAGRAM.md` was the breakpoint-only tallbar from mid-session. `WORKFLOW_DIAGRAM_BP05.md` was a separate tallbar for the 4-encoding module.

**Why archived**: both are replaced by `../PIPELINE_DIAGRAM.md` in the main folder, which is one unified view covering everything (01 through 06) in one place.

**Keep because**: if you ever want to cite just the breakpoint sub-workflow in the manuscript (without the diagnostics), the standalone `WORKFLOW_DIAGRAM.md` is a good source figure.

### `STEP28_patch.txt`, `STEP35_patch.txt` (inside `architecture_classifier/`)

Notes-only patches describing changes to make to your existing STEP28 and STEP35 scripts to use the within_group_ordering helper. These are documentation, not executable code — read them if/when you apply the heatmap ordering upgrade.

---

## If you want to use anything here

1. **Heatmap ordering upgrade**: most likely to be useful. Rename the three colliding scripts (STEP36, STEP40, STEP41 → STEP42/43/44 or whatever convention you prefer), drop them into your codebase alongside `within_group_ordering.R`, and you have the V1/V2/V3 comparison output that shows whether your heatmap ordering is exposing real substructure.

2. **Architecture classifier**: probably don't. See critique in `docs/SESSION_AUDIT.md`.

3. **Original workflow diagrams**: already superseded. Keep as reference.
