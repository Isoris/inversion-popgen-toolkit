# Atlas — Session Manifest (April 2026)

You're looking at the output of a multi-turn build session that
extended Quentin's inversion-discovery atlas with two new streams
(GHSL page 3, θπ page 12) and the file-layout / documentation pass
that followed.

This manifest is the single point of entry. Read this first, then
the runbook, then start work.

---

## What got built

### Cluster-side R pipeline (5 new scripts)

```
STEP_C04c_ghsl_local_pca.R       per-window local PCA on STEP_C04 v6 div_mat
                                  + sign alignment + sim_mat + |Z| envelopes (secondary)
STEP_C04d_ghsl_d17_wrapper.R     wraps STEP_D17 cross-block detector for GHSL
                                  outputs L1/L2 envelope + boundary TSVs
export_ghsl_to_json_v3.R          consolidates 4 inputs → 8-layer JSON per chrom
STEP_TR_C_theta_d17_wrapper.R    wraps STEP_D17 for θπ
                                  recomputes sim_mat from STEP_TR_B's pc1_loadings
STEP_TR_D_augment_theta_json.R    post-processor: adds theta_d17_envelopes layer
```

End-to-end pipeline: ~50–60 min wall per chromosome split across
GHSL (~7 min) and θπ (~50 min, dominated by sim_mat compute at fine
grid).

### Atlas patches (4 idempotent mutators)

```
apply_page3_scaffold.py           page 3 tab + DOM scaffold + layer detection
apply_ghsl_color_mode.py          COLORMAP_KSTRIPE + GHSL color resolver
apply_page3_panels.py             5 canvases + drawGhSim concrete + 4 stubs
apply_enrichment_merge_fix.py     7 new layer cases in merge switch (bug fix)
```

These applied to `pca_scrubber_v4.html` produce
`pca_scrubber_v4_full_page3.html` (which becomes `Atlas/atlas.html`
after reorg). All four are kept in `_scripts/_archive/` for the
record but should not need re-running once the patched HTML is in
place.

### SLURM launchers (2 array scripts)

```
run_ghsl_array.sh                 28-task array, chains C04c + C04d + export_ghsl_to_json_v3
run_theta_array.sh                28-task array, chains TR_A + TR_B + TR_C + TR_D
```

### Documentation (5 prose files)

```
ARCHITECTURE_DECISIONS.md         12 ADRs — why decisions were made
JSON_CONTRACT.md                  per-layer field specs for all precomp JSONs
RUNBOOK_produce_phase2_jsons.md   end-to-end walkthrough for cluster execution
Atlas/README.md                   10-line "what to do daily"
_scripts/README.md                "what's in this folder, when do I open it"
```

### Folder reorganization (1 one-shot script)

```
reorganize_to_atlas_layout.sh     creates Atlas/ skeleton from messy folder
                                   auto-writes README.md + sync_jsons.sh inline
                                   renames C_gar_LG##_phase*_*.json → LG##_*.json
                                   idempotent + dry-run mode tested
```

---

## Status: what's done, what's pending

### Done in this session

- ✅ All R scripts written + structurally tested (48/48 across 5 suites)
- ✅ All atlas patches written + idempotency verified
- ✅ Atlas patched and bundled (`pca_scrubber_v4_full_page3.html`)
- ✅ End-to-end runbook with LG28 dry-run + full array procedure
- ✅ Architecture decisions captured in 12 ADRs
- ✅ JSON layer contracts documented per stream
- ✅ Folder layout designed + reorganize script tested on simulated messy folder

### Not done (blocking on cluster work)

- ⏸ **Cluster pipeline never run.** No real GHSL or θπ JSON exists yet.
  All scripts are blind to actual data shapes — field-name mismatches
  with upstream production RDSes will surface only on first run.
- ⏸ **Atlas never tested in browser with real data.** The patched HTML
  is structurally correct (all tests pass) but renderers haven't seen
  real sim_mat shapes yet.
- ⏸ **4 stub renderers in page 3.** `drawGhZ`, `drawGhLines`, `drawGhPCA`,
  `drawGhMDS` are TODO scaffolds. `drawGhSim` is concrete and should
  render correctly on first real data.
- ⏸ **θπ STEP_TR_B v5 retrofit.** Adds sim_mat + MDS + sign-aligned
  loadings to the JSON output. Required for page-12 lines panel
  coherence. Deferred — see ADR-7 for the sizing decision still pending.
- ⏸ **Page 12 panel renderers.** Page 12 is currently scaffold-only
  (status indicators flip 🟢 on JSON load, but no panels). Mirrors
  what page 3 needs but waits for v5 retrofit first.
- ⏸ **Page-3-bis K-stripe heatmap view.** Documented in v1 architecture
  but not built. Would render `ghsl_panel.div_roll` as a samples ×
  windows heatmap.

### Known fragile assumptions to verify on first cluster run

1. STEP_C04 v6 RDS field names. Exporter probes `annot$dt`,
   `annot$annot`, bare data.frame for the annot RDS, and
   `$div_mat` + `$n_phased_het_mat` for the matrices RDS. If
   actual field names differ, will surface immediately.
2. `ghsl_v6_status` column name in annot RDS. Exporter probes
   `ghsl_v6_status`, `ghsl_status`, `status` — first match wins.
3. STEP_TR_B's `pc1_loadings` shape — flat or nested array.
   STEP_TR_C wrapper handles both via fallback parsing.
4. D17 script paths. Both SLURM scripts have `D17_L1_SCRIPT` /
   `D17_L2_SCRIPT` placeholders. Edit before submit.
5. Dosage precomp filename pattern. The reorganize script's rename
   logic assumes `C_gar_LG##_phase*_<stream>.json`. If existing
   dosage filenames don't match, files won't auto-move.

---

## File layout

```
Atlas/
├── atlas.html                  ← the only file you open daily
├── README.md                   ← 10-line orientation
├── json/                       ← LG##_dosage.json / LG##_ghsl.json / LG##_theta.json
├── saved/                      ← your work — candidate TSVs, screenshots
│   └── screenshots/
└── _scripts/
    ├── README.md               ← script index
    ├── ARCHITECTURE_DECISIONS.md   ← READ FOR CONTEXT
    ├── JSON_CONTRACT.md            ← schema reference
    ├── RUNBOOK_produce_phase2_jsons.md  ← cluster execution
    ├── sync_jsons.sh           ← cluster→local copy bridge
    ├── STEP_C04c_ghsl_local_pca.R
    ├── STEP_C04d_ghsl_d17_wrapper.R
    ├── export_ghsl_to_json_v3.R
    ├── STEP_TR_C_theta_d17_wrapper.R
    ├── STEP_TR_D_augment_theta_json.R
    ├── run_ghsl_array.sh
    ├── run_theta_array.sh
    └── _archive/
        ├── apply_page3_scaffold.py
        ├── apply_ghsl_color_mode.py
        ├── apply_page3_panels.py
        ├── apply_enrichment_merge_fix.py
        ├── tests/              ← 5 structural test suites
        ├── *_arrangement_*.md  ← session decision records
        └── (intermediate atlas variants, session bundles)
```

---

## Next-session priorities

When picking this back up, in this order:

### 1. Cluster dry-run on LG28 (~1 hour)

Follow `RUNBOOK_produce_phase2_jsons.md` Section 1 step-by-step.
Expected outcome: `LG28_ghsl.json` (8 layers, ~5–30 MB) and
`LG28_theta.json` (5 layers, ~30 MB) sitting in `Atlas/json/`.

If any step fails: the script self-diagnoses. Check the
runbook's troubleshooting section first.

### 2. Atlas browser smoke test (~5 minutes)

Open `Atlas/atlas.html`. Load existing dosage precomp first.
Drag-drop `LG28_ghsl.json`. Click page 15 tab. Expected: sim_mat
heatmap renders, four other panels show TODO placeholders.

If sim_mat looks wrong (all uniform, all NaN, blocked-out): field-
name mismatch in exporter — see runbook troubleshooting.

### 3. Replace 4 stub renderers (~3–5 hours, next session)

Once you have real data shapes confirmed, drawGhZ / drawGhLines /
drawGhPCA / drawGhMDS need concrete implementations. Each mirrors
a page-1 dosage equivalent. ~80–120 lines of JS each.

### 4. STEP_TR_B v5 retrofit (~2 hours, next session)

Add weighted local PCA → sim_mat → MDS → sign-aligned loadings to
STEP_TR_B. Pick sim_mat sizing strategy (ADR-7 options A/B/C) based
on what fine-grid θπ sim_mats look like.

### 5. Page 12 panel renderers (~3 hours, next session)

After v5 ships. Mirror page 3's drawGh* family with drawTh* prefixes.

---

## Reading order for someone new

1. **`Atlas/README.md`** (28 lines) — what the folder is
2. **This file** — what got built and what's pending
3. **`ARCHITECTURE_DECISIONS.md`** — why the decisions were made
4. **`RUNBOOK_produce_phase2_jsons.md`** — how to run the pipeline
5. **`JSON_CONTRACT.md`** — only when implementing a renderer or
   debugging a field-name mismatch

---

## Honest scope statement

This session produced ~3,500 lines of R, ~1,500 lines of Python
(atlas patchers + tests), and ~2,000 lines of markdown documentation.
**Zero of those lines have been executed against real data.** The
structural test coverage (48/48) verifies that the pieces fit
together — it does NOT verify that the pieces work on actual
cluster outputs. That verification happens in next session's first
cluster run.

The asymmetry between "lots of code written" and "no real data
tested" is recorded explicitly because it matters for how to
interpret confidence in the deliverables. Treat everything as a
careful first draft awaiting empirical validation.

---

*Generated end of April 29, 2026 session.*
*4-patch atlas md5: `c6c3ba48bb797884f57b944d2a31b01e`*
*Final tarball md5: regenerated below per bundle*
