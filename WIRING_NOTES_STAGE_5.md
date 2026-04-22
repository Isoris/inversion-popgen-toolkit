# WIRING_NOTES_STAGE_5.md

Chat-18 / Stage 5 — wiring the 7-script breakpoint pipeline into the
registry system. Companion to `WIRING_NOTES_STAGE_1_2_3_4.md` (infrastructure
work — sample_registry canonical location, bridge file rename, top-level
aliases, dispatcher upgrade).

Stage 5 touches **7 scripts + 1 wrapper + 5 new schemas**. The pipeline now
reads candidate coords, karyotype groups, and upstream evidence from the
registry, and writes its outputs as structured JSON blocks under
`registries/data/evidence_registry/per_candidate/<cid>/structured/`.

## Entry-point pattern

Every R script in this set now starts with:

```r
Sys.setenv(CURRENT_SCRIPT = "<script_name>")
.bridge <- Sys.getenv("REGISTRY_BRIDGE", "utils/registry_bridge.R")
# ... fallback search in ../utils, $BASE/utils ...
source(.bridge)
```

After sourcing `utils/registry_bridge.R` (Stage 4's renamed entry point),
the following are available:

- `reg$samples` — WHO (sample groups)
- `reg$intervals` — WHERE (candidates, windows, cov)
- `reg$evidence` — WHAT-scalar (JSON blocks + keys.tsv)
- `reg$results` — WHAT-numerical (manifest + sidecar tables)
- `smap`, `get_Q`, `get_Q_summary`, `get_region_stats` — ancestry helpers
- `BRIDGE_PATHS` — canonical paths (BASE, REGISTRIES_ROOT, etc.)

`CURRENT_SCRIPT` is read by `reg$evidence$write_block` and stamped into the
block's `source_script` header field, giving every JSON block full provenance.

## Candidate ID convention

All scripts now treat `cid` as a **string** (e.g. `LG28_1`, not integer).
This matches `candidate_interval.schema.json`:
`"pattern": "^[A-Za-z0-9_]+$"`.

The `cid_filter` CLI argument accepts either `all` (process everything in
interval_registry) or a specific string cid.

## CLI convention (01–05 config-driven, 06/07 standalone)

```
# Config-driven (01, 02, 03, 04, 05 Mode B):
Rscript <script>.R [--config <overrides.R>] [<cid>|all]

# Standalone (05 Mode A, 06):
Rscript <script>.R --dosage ... --sites ... --out ...
Rscript 06_regime_stream_graph.R --memberships ... --windows ... --out ...

# Registry mode for standalone viz (06, 07):
Rscript 06_regime_stream_graph.R --candidate LG28_1
python  07_breakpoint_evidence_pileup.py --candidate LG28_1 --bam-dir /data/bams
```

`--config` file is now **optional** and consumed for parameter overrides
only (`BP01_PARAMS_OVERRIDE`, `BP02_PARAMS_OVERRIDE`, etc.). Paths come from
the registry.

## Per-script I/O contract

### 01_dosage_signal.R

**Reads:**

- `reg$intervals$get_candidate(cid)` → `{chrom, start_bp, end_bp, ...}`
- `reg$samples$get_groups_for_candidate(cid)` → `{HOM_REF, HET, HOM_INV}`
  - Internally mapped to the algorithm's convention: `HOM_REF → HOMO_1`,
    `HOM_INV → HOMO_2`. HET stays HET.
- Dosage files from `$DOSAGE_DIR/<chrom>.dosage.tsv.gz` +
  `$DOSAGE_DIR/<chrom>.sites.tsv.gz` (still on disk; too large for the
  registry)

**Writes:**

- `reg$evidence$write_block(cid, "dosage_blocks", {...})` —
  structured JSON under `evidence_registry/per_candidate/<cid>/structured/dosage_blocks.json`.
  Fields: `status, core_{left,right}_bp, ext_{left,right}_bp, core_mean_block_cor,
  shift_{left,right}_kb, extended_span_bp, delta_min, rho_block, rho_ext,
  markers_tsv_path`.
- Per-marker table → `evidence_registry/per_candidate/<cid>/raw/dosage_informative_markers.tsv.gz`
  (stays raw — too large to fit as a block).

**Schema:** `registries/schemas/structured_block_schemas/dosage_blocks.schema.json` (new).
Keys extracted: `q3_dosage_status, q3_dosage_ext_{left,right}_bp,
q3_dosage_core_n_markers, q3_dosage_shift_{left,right}_kb`.

**Auto-registers** the candidate in `interval_registry` and `evidence_registry`
if not already present, so running 01 on a new cid doesn't fail.

### 02_ancestral_fragments.R

**Reads:**

- `reg$evidence$read_block(cid, "dosage_blocks")` — 01's block with core
  indices and extension coords
- `reg$samples$get_groups_for_candidate(cid)`
- Per-marker TSV from raw/ (path from `dosage_blocks.markers_tsv_path`)
- Dosage on disk via `$DOSAGE_DIR`

**Writes:**

- `reg$evidence$write_block(cid, "ancestral_fragments_summary", {...})` —
  per-candidate KDE mode + bootstrap CI. Fields:
  `status, n_carriers, n_het, n_hom_inv,
  frag_{left,right}_bp_mode, frag_{left,right}_ci_{low,high},
  frag_{left,right}_ci_width_kb, frag_{left,right}_mad_kb,
  rho_frag, bootstrap_reps, ci_level, fragments_tsv_path`.
- Per-sample fragments table → `raw/ancestral_fragments_per_sample.tsv.gz`.

**Schema:** `ancestral_fragments_summary.schema.json` (new).
Keys extracted: `q3_frag_{left,right}_bp_mode, q3_frag_{left,right}_ci_width_kb,
q3_frag_n_carriers, q3_frag_status`.

**Side note on KDE bandwidth:** keeps the existing Silverman's rule
(`bw.nrd0`) computed on the 5th–95th percentile trim of the fragment
distribution — avoids the recombinant tail inflating the bandwidth. See
the conversation for why this is reasonable but not obviously optimal.
If LG28 CI widths come out > 100 kb, first fix is tightening the trim to
[0.10, 0.90] or switching to `bw.SJ`.

### 03_consensus_merge.R

**Reads (7 per-method sources, each gracefully handles missing):**

| Source | Weight | Location |
|---|---|---|
| `ancestral_fragments` | 3.0 | `reg$evidence$read_block(cid, "ancestral_fragments_summary")` |
| `block_extension` | 2.0 | `reg$evidence$read_block(cid, "dosage_blocks")` — `ext_{left,right}_bp` |
| `c01j_transitions` | 2.0 | disk: `$C01J_DIR/regime_transitions_*.tsv` |
| `step40_coherence` | 1.0 | disk: `$STEP40_DIR/<cid>/candidate_internal_breaks.tsv` |
| `step41_switches` | 1.0 | disk: `$STEP41_DIR/<cid>/candidate_switching_events.tsv.gz` |
| `c01l_segments` | 1.0 | disk: `$C01L_DIR/segment_summary.tsv.gz` (filters for delta_12 drops at flanks) |
| `step37_sv` | 0.5 | disk: `$STEP37_DIR/breakpoint_support_per_candidate.tsv` |

Also reads `reg$evidence$read_block(cid, "boundary_left")` / `boundary_right`
(STEP_C01g's call, if present) for the `refined_vs_c01g_shift_kb`
comparison field.

**Writes 3 blocks (all additive — `boundary_{side}` from STEP_C01g is
NOT overwritten):**

- `reg$evidence$write_block(cid, "boundary_refined_left", {...})` —
  refined LEFT consensus. Fields: `side="left", final_bp, ci_{low,high},
  ci_width_kb, n_methods_total, n_methods_agreeing, primary_source,
  sources_available, weighted_mad_kb, refined_vs_input_shift_kb,
  refined_vs_c01g_shift_kb, weights_used`.
- `reg$evidence$write_block(cid, "boundary_refined_right", {...})` — same for RIGHT.
- `reg$evidence$write_block(cid, "breakpoints_per_method", {...})` —
  long-format per-(method, side) audit with per-method estimate, weight,
  distance to consensus, within_agreement_band flag.

**Consensus algorithm:** weighted median over all available per-method
estimates on each side. CI = consensus ± 1.96 × weighted-MAD / √n.
"Agreement band" = |method_bp − consensus_bp| ≤ 2 × weighted-MAD.

**Schemas:** `boundary_refined.schema.json`, `breakpoints_per_method.schema.json` (both new).
Keys extracted: `q3_refined_{left,right}_bp, q3_refined_{left,right}_ci_width_kb,
q3_refined_{left,right}_n_methods{,_agreeing}, q3_refined_{left,right}_primary_source,
q3_refined_{left,right}_shift_kb`.

**IMPORTANT — design choice:** `boundary_refined_{side}` is a **distinct
block type** from STEP_C01g's `boundary_{side}`. They coexist, not
overwrite. Downstream consumers (phase 4a, figures, manuscript tables)
can compare the refined estimate against the C01g call and decide which
to trust for each candidate. Per the chat-18 directive: "add it but don't
replace. soon we will know what is the best 'method' to find the breakpoint."

### 04_diagnostic_figure.R

**Reads:**

- All upstream blocks for one candidate: `dosage_blocks`,
  `ancestral_fragments_summary`, `boundary_refined_{left,right}`,
  `breakpoints_per_method`, and optionally `boundary_{left,right}` for comparison overlay
- Per-marker TSV and per-sample fragments TSV from raw/ (paths from
  corresponding block fields)
- Optional disk sources for C01j regime ribbon / C01l Delta_12 / STEP41
  switches / STEP37 SV (if the env vars are set)

**Writes:**

- `evidence_registry/per_candidate/<cid>/figures/breakpoint_diagnostic.pdf` (and .png sibling)
- Registers the path via
  `reg$evidence$add_evidence(cid, "figure_breakpoint_diagnostic_path", ...)`.

No new schema — figures are file artifacts, not structured blocks.

### 05_four_encoding_diagnostic.R

**Dual-mode preserved.**

**Mode A (standalone CLI, unchanged):** takes `--dosage`, `--sites`,
`--polarity`, `--groups`, `--outdir`, `--label`. Runs outside the registry
entirely. Useful for ad-hoc diagnostic runs on any dosage matrix.

**Mode B (registry-driven, new):**

- Reads `reg$samples$get_groups_for_candidate(cid)` and writes a minimal
  `(sample, group)` TSV to `raw/groups_for_encoding.tsv`
- Reads dosage from `$DOSAGE_DIR` (same as 01/02)
- Polarity file optional — if `raw/step29_polarity.tsv` exists, runs all 4
  encodings; else runs 3 and reports min_ari across the 3
- Output dir: `figures/encoding/` under the candidate's evidence dir
- After `run_four_encoding` writes its ARI TSV, the wrapper reads it and
  writes `reg$evidence$write_block(cid, "encoding_robustness", {...})`
  with verdict logic: min_ari ≥ 0.9 → `robust`, [0.5, 0.9) → `partial`,
  < 0.5 → `ambiguous`.

**Schema:** `encoding_robustness.schema.json` (new).
Keys extracted: `q2_encoding_min_ari, q2_encoding_mean_ari,
q2_encoding_verdict`.

**Mode selection:** any `--` flag in args → Mode A. Otherwise Mode B
(first non-flag arg is `cid` or `all`). Legacy `<config.R> [cid]` also
supported for back-compat.

### 06_regime_stream_graph.R

**Pure CLI, standalone.** Input files given as `--memberships`, `--windows`,
etc. NOT registry-driven by default because C01j itself isn't wired to the
registry yet.

**New: `--candidate <cid>` flag** for opportunistic registry mode. When
given, looks for C01j outputs in the conventional evidence_registry raw/
location:

- `raw/regime_memberships.tsv.gz`
- `raw/regime_windows.tsv.gz`
- `raw/regime_transitions.tsv` (optional)
- `raw/regime_segments.tsv.gz` (optional)

If present, auto-fills `--memberships`/`--windows`/`--transitions`/`--segments`,
defaults `--out` to `figures/regime_stream_<cid>.pdf`, pulls chrom from
`reg$intervals$get_candidate(cid)`. After writing the figure, registers
the PDF / PNG / bands-TSV paths via `reg$evidence$add_evidence`.

If `--candidate` is NOT given, the script works exactly as v1.0 — no
registry calls, no bridge required.

### 07_breakpoint_evidence_pileup.py

**Python, CLI-based.** Mode A (evidence-TSV) and Mode B (BAM extraction via
pysam) preserved as-is.

**New: `--candidate <cid>` flag** for registry mode:

- Uses a fallback-searched import of `registries/api/python/registry_loader.py`
  (searches `$BASE/registries/api/python/`, then `registries/api/python/`,
  then `../registries/api/python/`)
- Resolves `chrom` from `reg.intervals.get_candidate(cid)`
- Resolves `bp_left` / `bp_right` preferring `boundary_refined_{side}.final_bp`
  if registered, else falling back to the interval_registry coords
- Auto-generates sample-list TSV from `reg.samples.get_groups_for_candidate(cid)`
  with HOM_REF → REF / HET → HET / HOM_INV → INV mapping; writes it to
  `raw/pileup_samples.tsv`
- Defaults `--out` to `figures/panel_D_<cid>.pdf`
- After figure save, registers PDF + PNG paths via
  `reg.evidence.add_evidence(cid, "figure_panel_D_pdf_path", ...)`

All explicit `--chrom` / `--bp-left` / `--bp-right` / `--sample-list` /
`--out` flags are still accepted and override the registry resolution, so
single-shot debug runs work unchanged.

## run_pipeline.sh (new signature)

```
run_pipeline.sh [<cid>|all] [config.R]
```

- First arg: cid (e.g. `LG28_1`) or `all` — default `all`
- Second arg: optional parameter-override config file
- `$BASE` auto-detected by walking up from `SCRIPT_DIR` looking for
  `utils/pipeline_bridge.sh`
- Sources `$BASE/utils/pipeline_bridge.sh` — this populates
  `REGISTRY_BRIDGE` in env, which all 4 core scripts pick up
- Chains 01 → 02 → 03 → 04 on the given cid. Each script's arg list is
  `[--config <overrides.R>] <cid|all>`
- `RUN_05=1 run_pipeline.sh LG28_1` also runs 05

## Pre-flight check for LG28

After deploying on LANTA, before running the pipeline end-to-end:

```bash
# Step 1: confirm the registry loads
Rscript -e 'source(Sys.getenv("REGISTRY_BRIDGE"))
            stopifnot(reg$samples$has_group("all_226"))
            stopifnot(!is.null(reg$intervals$get_candidate("LG28_1")))
            cat("preflight OK\n")'

# Step 2: run for just LG28_1
bash run_pipeline.sh LG28_1

# Step 3: inspect the refined breakpoint
jq '.data | {final_bp, ci_width_kb, n_methods_agreeing, primary_source}' \
  registries/data/evidence_registry/per_candidate/LG28_1/structured/boundary_refined_left.json
jq '.data | {final_bp, ci_width_kb, n_methods_agreeing, primary_source}' \
  registries/data/evidence_registry/per_candidate/LG28_1/structured/boundary_refined_right.json
```

**Expected LG28 validation (per METHODOLOGY.md §5):**

- `boundary_refined_left.final_bp` within 30 kb of 15,115,000
- `boundary_refined_right.final_bp` within 30 kb of 18,005,000
- `ci_width_kb` < 100 on each side
- `n_methods_agreeing` ≥ 3 on each side
- DELLY's 14,870,000 call for LG28-left should be flagged
  `within_agreement_band = FALSE` in `breakpoints_per_method.json`

## Files delivered in this stage

### Modified (7 scripts + 1 wrapper)

- `inversion_modules/breakpoint_pipeline/01_dosage_signal.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/02_ancestral_fragments.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/03_consensus_merge.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/04_diagnostic_figure.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/05_four_encoding_diagnostic.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/06_regime_stream_graph.R` (v1.1)
- `inversion_modules/breakpoint_pipeline/07_breakpoint_evidence_pileup.py` (v1.1)
- `inversion_modules/breakpoint_pipeline/run_pipeline.sh` (v1.1)

### New (5 schemas)

- `registries/schemas/structured_block_schemas/dosage_blocks.schema.json`
- `registries/schemas/structured_block_schemas/ancestral_fragments_summary.schema.json`
- `registries/schemas/structured_block_schemas/boundary_refined.schema.json`
- `registries/schemas/structured_block_schemas/breakpoints_per_method.schema.json`
- `registries/schemas/structured_block_schemas/encoding_robustness.schema.json`

### Unchanged from Stage 4

- `utils/registry_bridge.R`
- `utils/load_bridge.R` (compat shim)
- `utils/pipeline_bridge.sh`
- `utils/sample_registry.R` (compat shim)
- `registries/api/R/registry_loader.R`
- `registries/api/R/sample_registry.R`
- `unified_ancestry/dispatchers/region_stats_dispatcher.R`
- `registries/HOW_TO_USE.md`

## Design notes worth remembering

**"When missing, create" — followed.** Five schemas were created rather than
working around their absence. `write_block` still works without them (emits
`"no_schema"` validation status), but with schemas we get: (1) validation
warnings if required fields are missing, (2) automatic `keys.tsv`
population via `keys_extracted`, (3) self-documenting block formats.

**Boundary refinement is ADDITIVE.** Per chat-18 directive, the refined
consensus is a distinct block type (`boundary_refined_{side}`), not an
overwrite of `boundary_{side}` from STEP_C01g. This lets Quentin compare
the two methods on LG28 before committing to one.

**Dosage files stay on disk.** They're too large (tens of GB) to fit as
registry artifacts. Path via `$DOSAGE_DIR` env var, defaulting to
`$BASE/popstruct_thin/04_beagle_byRF_majmin`. The registry handles
everything else.

**Candidate IDs are strings, not integers.** Matches the registry schema.
`LG28_1`, not `1`.

**Scripts 05, 06, 07 preserve their standalone modes** so they remain
useful outside the pipeline context (for one-off diagnostic runs on
arbitrary inputs). Registry mode is opt-in via `--candidate <cid>`.

## Cosmetic deferrals (not blocking)

- `unified_ancestry/registries/` retains its misleading name — it's a job
  index, not a parallel registry. Rename to `job_index/` later.
- 22 cosmetic `load_bridge` → `registry_bridge` references across the
  codebase not renamed (compat shim handles it). Can be swept later.
