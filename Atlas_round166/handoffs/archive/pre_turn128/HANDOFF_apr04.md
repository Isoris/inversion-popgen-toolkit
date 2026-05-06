# HANDOFF — start of next session

**Last session ended:** April 29, 2026
**State:** documentation complete, cluster pipeline never run, atlas never tested with real data.

---

## In one sentence

Run LG28 on the cluster, look at the JSON, look at the atlas. Everything else waits.

---

## First action — LG28 cluster dry-run (~1 hour)

Open `Atlas/_scripts/RUNBOOK_produce_phase2_jsons.md`, follow Section 1 step-by-step. Produces:

- `Atlas/json/LG28_ghsl.json` (~5–30 MB, 8 layers)
- `Atlas/json/LG28_theta.json` (~30 MB, 5 layers)

If something fails, the runbook's troubleshooting section names the 5 most likely failure modes (field-name mismatches in upstream RDSes, D17 script paths, JSON shape questions). Most failures are 1-line fixes once seen.

## Second action — atlas browser test (~5 min)

Open `Atlas/Inversion_atlas.html`. Load existing dosage precomp first. Drag-drop `LG28_ghsl.json`. Click page 15 tab.

**Expected:** sim_mat heatmap renders. Four other panels show TODO placeholders.

**If sim_mat looks wrong** (uniform, all NaN, blocked-out): field-name mismatch in `export_ghsl_to_json_v3.R`. Look at console (F12), report back.

## Third action — bring back evidence

Whichever of these is easiest to send when picking the next session up:

- A screenshot of page 3 with real GHSL data loaded
- A structural dump (`length(layers)`, `dim(sim_mat)`, sample sizes) from one JSON
- Console error log if anything throws

That evidence unblocks the next iteration.

---

## Active todo, ranked

| # | Task | Status | Blocking on |
|---|---|---|---|
| 1 | LG28 dry-run on cluster | not started | you (cluster access) |
| 2 | Browser smoke test with real data | blocked | task 1 |
| 3 | Replace 4 stub renderers (drawGhZ, drawGhLines, drawGhPCA, drawGhMDS) | blocked | task 2 (need real shapes) |
| 4 | STEP_TR_B v5 retrofit (sim_mat + MDS for θπ) | not started | sim_mat sizing decision (ADR-7) |
| 5 | Page 12 panel renderers | blocked | task 4 |
| 6 | Page-3-bis K-stripe heatmap view | not started | not blocking — can do anytime |
| 7 | Synthetic test JSON | not started | optional pre-task-1 sanity check |
| 8 | Manuscript methods section | premature | tasks 1–5 (need real results) |

**Critical path:** 1 → 2 → 3, in order. Everything else is dependent or parallel.

---

## Don't redo this work

These are decisions already made. Don't re-derive them next session unless real data forces a reversal:

- **Detector hierarchy is asymmetric** (ADR-2). GHSL's PASS-runs primary; θπ's |Z|-threshold + D17 both primary. Don't unify.
- **D17 wrappers, not modifications** (ADR-3). Don't touch D17 source.
- **D17 is O(N²), not O(N³)** (ADR-4). Don't re-investigate.
- **K=3 GHSL palette mirrors dosage K=3** (ADR-5). Don't pick different colors.
- **Sim_mat as upper-triangle Float32 flat array** (ADR-6). Don't switch to dense or banded for GHSL.
- **Atlas patches are one-shot, not source of truth** (ADR-8). Don't re-run patchers unless rebuilding atlas.html from scratch.
- **Default to ship cheap optional outputs** (ADR-12). Don't remove fields "for safety."

Full reasoning in `Atlas/_scripts/ARCHITECTURE_DECISIONS.md`.

---

## Don't trust without verification

Five claims in last session's deliverables that haven't been validated against real data and could be wrong:

1. **STEP_C04 v6 RDS field names.** Exporter probes `annot$dt`, `annot$annot`, bare data.frame for the annot RDS, and `$div_mat` + `$n_phased_het_mat` for matrices. If actual names differ → exporter fails immediately, easy fix.
2. **`ghsl_v6_status` column name.** Exporter tries `ghsl_v6_status`, `ghsl_status`, `status` — first match wins. Wrong match → empty `ghsl_envelopes` layer.
3. **STEP_TR_B `pc1_loadings` shape** — flat or nested. STEP_TR_C handles both via fallback parsing, but might not.
4. **D17 wall times** were extrapolated from numpy benchmark, not measured in R. Could be off by 2–3×.
5. **STEP_TR_B writes to `JSON_OUT_DIR/$CHROM/...`** per `00_theta_config.sh`. If the actual output path is different, STEP_TR_C/D won't find the JSON.

If something in the runbook fails, suspect one of these first.

---

## Open questions to resolve eventually

- **θπ sim_mat sizing.** Three options on the table (ADR-7): coarse grid / banded / int8 quantized. Decision deferred until you see real fine-grid θπ sim_mats. STEP_TR_B v5 retrofit can't ship without picking one.
- **Page-3 arrow-key navigation step size.** Dosage units (cross-page muscle memory) or GHSL units (matches what's on screen)? Default in current code is dosage; Alt+→ for native GHSL. Confirm at first browser test.
- **Tab double-click vs explicit toggle for page-3-bis.** Lean toward explicit "view: default | bis" toggle for discoverability.

---

## Critical reminders

- **226 pure C. gariepinus hatchery cohort ONLY.** Not F1 hybrid. Not C. macrocephalus.
- **9× coverage, hatchery F1 violates HWE** — folded SFS, doSaf 5, doMajorMinor 1.
- **K clusters = hatchery broodline structure**, not geographic populations.
- **GHSL grid is 5kb fixed.** θπ at win10000.step2000 (~16,500 windows on LG28).
- **NAToRA pruning yields 81 unrelated** for downstream analyses needing independence.
- **BAM manifest:** `pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv` (col3 = filtered BAM).
- **Phase blocks from WhatsHap on short reads:** 44–819 bp.
- **Rscript binary:** `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`.

---

## tl;dr

1. Run LG28 from the runbook
2. Check the atlas in browser
3. Send back what you saw

Everything else waits on those three things.

---

*Use this file as the entry point. SESSION_MANIFEST.md has the comprehensive view; this is the action-oriented short version.*
