# S6 — Dosage heatmap streaming viewer (bridge + sampling modes)

**Status:** spec only. No code this turn.
**Scope:** consolidates four user prompts (atlas bridge, viewer architecture, sampling modes, minimal local server) into one design.
**Supersedes / extends:** patch `P4_1_dosage_chunk_endpoint.md` — that patch is the *minimum-viable* slice of this spec (Mode B with `mode=raw` only). Implement the bridge first against the existing atlas, then grow into the standalone viewer.
**Used by:** `candidateDosageHeatmapHtml(c)` and the per-candidate panel in `Inversion_atlas.html`; future standalone `dosage_viewer/` app.

---

## TL;DR

The atlas already has a working dosage heatmap renderer (`renderDosageHeatmap` + `selectTopMarkers` + `_wireCandidateDosageHeatmap`). What's missing is the **data bridge** that fills `state.data.dosage_chunks.chunks` from `04_dosage_by_chr/<chrom>.{sites,dosage}.tsv.gz`. We do not redesign the renderer.

This spec covers two deliverables on the same backend:

1. **Bridge mode (P4.1, must-ship):** server endpoint that hands the existing renderer a chunk in the shape it already expects. This unblocks the candidate panel.
2. **Standalone viewer (future):** a separate `dosage_viewer/` app for free-form region/tile exploration — region mode, candidate mode, breakpoint mode, tile mode, and six sampling modes (`raw / even / random / variance / hybrid / aggregate`).

The bridge and the viewer share the same indexed backend store, the same sampling logic, and the same metadata envelope. The atlas just consumes one shape of response; the viewer consumes a richer one.

---

## Why both deliverables on one spec

Splitting them invites schema drift. If the bridge endpoint and the viewer endpoint each read the gzipped TSVs independently, we end up with two parsers, two cap policies, two missing-value conventions, and divergent metadata. One backend, two consumers.

The bridge is a **subset** of the viewer endpoint, not a sibling. Implement the shared core first.

---

## Reverse-engineered renderer contract (authoritative)

Read directly from the live atlas, not invented:

- **Index** (what `state.data.dosage_chunks` holds): `{ chrom, chunks: [{ start_bp, end_bp, url, ... }, ...] }`. Detected by `Array.isArray(data.dosage_chunks.chunks)` (see `detectSchemaAndLayers`).
- **Chunk** (what each chunk URL must return as JSON):
  ```jsonc
  {
    "chrom": "C_gar_LG28",
    "start_bp": 15000000,
    "end_bp":   18300000,
    "samples": ["CGA_001", "CGA_002", ...],            // length = n_samples; STRING IDs
    "markers": [
      { "pos_bp": 15010234, "missingness": 0.012, "diagnostic_score": null },
      ...
    ],                                                 // length = n_sites
    "dosage": [                                        // n_sites rows × n_samples cols
      [0, 1, 2, -1, ...],                              // missing encoded as -1 (NOT null)
      ...
    ]
  }
  ```

  Confirmed reads in the live file:
  - `chunk.markers[mi].pos_bp` / `.missingness` / `.diagnostic_score` (`selectTopMarkers`, line ~14566).
  - `chunk.dosage[mi]` is an array of per-sample dosage values, with `-1` (or `< 0`) treated as NA (line ~14602: `if (v == null || !Number.isFinite(v) || v < 0) continue;   // NA encoded as -1`).
  - `chunk.samples[i]` is a string sample ID, resolved against `state.data.samples[].id|cga|ind|sample` via `_buildSampleLookups` (line ~15128).
  - `_findCoveringChunk` matches by `ch.start_bp <= region.start_bp && ch.end_bp >= region.end_bp`, falling back to any overlap (line ~15084).

- **Renderer expectations** (do not change):
  - Orientation: `dosage` is **sites × samples**, not samples × sites. The renderer transposes for drawing.
  - Missing encoding: **`-1`** (the `null` in the original P4.1 doc was wrong — `selectTopMarkers` filters on `v < 0`, so a `null` becomes `NaN` via `Number.isFinite` and is also dropped, but `-1` is the documented convention).
  - Chunk index keyed by URL: same URL must return the same chunk (the cache uses URL/region as key in `_findCoveringChunk` and `_fetchAndCacheChunk`).

These four facts are the contract. Anything in the original P4.1 patch that contradicts them is wrong; this spec is the source of truth.

## Anchors in the live atlas

| Function | Approx. line | Role |
|---|---|---|
| `selectTopMarkers(chunk, region, capN, opts)` | 14553 | top-N selection inside a chunk; uses `diagnostic_score` if available on every marker, else variance fallback |
| `_findCoveringChunk(region)` | 15075 | resolves region → chunk via `state.data.dosage_chunks.chunks` |
| `_fetchAndCacheChunk(chunkRef)` | 15097 | fetches `chunkRef.url` and caches |
| `_buildSampleLookups(chunk, cand)` | 15128 | maps `chunk.samples[i]` (string) to cohort index via `state.data.samples` |
| `renderDosageHeatmap(canvas, region, cand, capN)` | 15306 | top-level renderer, calls `_findCoveringChunk` then `selectTopMarkers` |
| `candidateDosageHeatmapHtml(c)` | 15600 | static panel HTML; reads `state.data.dosage_chunks.chunks.length` for the loaded badge |
| `_wireCandidateDosageHeatmap(c)` | 15632 | wires the panel post-render; calls `renderDosageHeatmap` then re-selects via `selectTopMarkers` |
| Schema detector for `dosage_chunks` | 41196 | `ok = !!data.dosage_chunks && Array.isArray(data.dosage_chunks.chunks)` |

These are the strings to grep for. If line numbers shift, the function names won't.

---

## Source data layout

```
${base}/04_dosage_by_chr/
   C_gar_LG01.sites.tsv.gz
   C_gar_LG01.dosage.tsv.gz
   C_gar_LG02.sites.tsv.gz
   C_gar_LG02.dosage.tsv.gz
   ...
   C_gar_LG28.sites.tsv.gz
   C_gar_LG28.dosage.tsv.gz
```

(On the user's WSL path: `/mnt/e/01-catfish_assembly_manuscript_CGA/04_dosage_by_chr/`. On LANTA: under `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`.)

**`<chrom>.sites.tsv.gz`** — one row per site, sorted by `pos`:

| col | type | content |
|---|---|---|
| 1 | int | `pos` (1-based bp) |
| 2 | str | `ref` |
| 3 | str | `alt` |
| 4 | float | `missingness` (optional; `NA` allowed) |
| 5 | float | `diagnostic_score` (optional; `NA` allowed) |
| 6 | str  | `site_id` (optional; falls back to `chrom:pos:ref>alt`) |

**`<chrom>.dosage.tsv.gz`** — same row count as sites, one row per site:

- Each row is tab-separated, `n_samples` cells.
- Cells are `0`, `1`, `2`, or one of `NA / . / -1 / ''` for missing → normalized to **`-1`** on the wire.

**Sample IDs**: come from a per-base `samples.tsv` (one ID per line) or from the dosage header row if present. Detection rule: if line 1 of `<chrom>.dosage.tsv.gz` parses as integers, it's data; otherwise it's a header with sample IDs. Document this in `01_prepare_dosage_store.py` (see below).

---

## Architecture (shared backend, two consumers)

```
                    +--------------------------+
                    |  04_dosage_by_chr/       |   raw TSV.GZ (canonical)
                    |    *.sites.tsv.gz        |
                    |    *.dosage.tsv.gz       |
                    +-----------+--------------+
                                |
                                | 01_prepare_dosage_store.py
                                v
                    +--------------------------+
                    |  dosage_store/           |   indexed backend
                    |    <chrom>/              |   (Parquet partitioned by chrom)
                    |      sites.parquet       |
                    |      dosage.parquet      |
                    |    samples.tsv           |
                    |    manifest.json         |
                    +-----------+--------------+
                                |
                                | 02_run_server.py (FastAPI)
                                v
            +-------------------+-------------------+
            |                                       |
            v                                       v
   +------------------+                   +-------------------+
   |  Atlas bridge    |                   | Standalone viewer |
   |  /api/dosage/    |                   | dosage_viewer/    |
   |    chunk         |                   |   static/index    |
   |  (P4.1 shape)    |                   |   region / tile   |
   |                  |                   |   /candidate /    |
   |                  |                   |   breakpoint      |
   +------------------+                   +-------------------+
            |                                       |
            v                                       v
   Inversion_atlas.html                  standalone HTML page
   renderDosageHeatmap                   Canvas/WebGL renderer
```

One preprocessing step. One on-disk store. One server, two endpoint families.

---

## Storage format choice: Parquet

Pick **Parquet partitioned by chromosome**. Reasons:

- `pyarrow` is already in the popstats server's stack.
- Row-group filter pushdown gives O(log n) random access by `pos` if the file is sorted.
- Columnar — reading a region pulls only the row-group(s) that overlap.
- Portable; no HDF5 driver headaches on LANTA, no Zarr metadata sprawl.

Reject HDF5 (driver fragility), Zarr (good for cloud, overkill for one workstation), Arrow IPC (no built-in row-group filter), per-chrom NumPy `.npy` (no metadata, no filters).

If a future chromosome has >5 M sites and Parquet row-group seek becomes the bottleneck, revisit and add a positional `.tbi`-style index or move to Zarr. Not now.

### Layout on disk

```
dosage_store/
  manifest.json                   { schema_version, n_samples, chroms: [...], created_utc }
  samples.tsv                     # one sample ID per line, in column order
  C_gar_LG01/
    sites.parquet                 # cols: pos, ref, alt, missingness, diagnostic_score, site_id
                                  # row group size: 50_000 sites
                                  # sorted by pos; min/max stats enable predicate pushdown
    dosage.parquet                # cols: pos, s0, s1, ..., s225 (int8: 0/1/2/-1)
                                  # same row order as sites.parquet
                                  # row group size: 50_000 sites
  C_gar_LG02/
    ...
```

Why two parquets per chromosome instead of one wide file: site metadata is hot (every region query reads it); the dosage matrix is cold (only fetched after sites have been picked). Splitting them halves the bytes read for sample queries and for the `aggregate` mode.

`int8` for dosage saves 4× over int32. `-1` fits.

---

## Sampling modes (six)

Implemented once in the server. Both the bridge and the viewer endpoints accept `mode=`. The bridge defaults to `mode=raw` for backward compatibility with the existing renderer (which already does its own top-N selection via `selectTopMarkers`). The viewer defaults to `mode=hybrid`.

| mode | what it does | when to use | bias warning |
|---|---|---|---|
| `raw` | return all sites in region; if `n_sites > max_sites`, return error with a suggested mode | small regions where every site matters | none |
| `even` | pick `max_sites` sites evenly across `[start, end]` by genomic position | gene-conversion-like tracts, breakpoint transitions, local ancestry | none — positionally unbiased |
| `random` | reproducible random subset using `seed` (default 1); always re-sort by `pos` before returning | unbiased visual sample | none |
| `variance` | compute per-site dosage variance across samples; return top `max_sites`; re-sort by `pos` | inversion-band visualization, structural-haplotype separation | enriches family-structure markers, repeats, paralogy — not a neutral view |
| `hybrid` | 70 % `even` + 30 % `variance`, dedup, sort by `pos` | **default** for exploration | mild — same as `variance` for the 30 % share |
| `aggregate` | bin region into `n_bins` equal-width bp bins; per bin per sample, return mean dosage (NA-aware) | chromosome-scale or whole-genome overview | smooths real bands; not for tier-3 evidence |

### Critical contract (applies to all modes)

> **After selection, always re-sort returned sites by genomic position.** The heatmap x-axis is positional. `variance` and `hybrid` selections produce out-of-order indices; sort before sending.

This rule is duplicated in the renderer (`selectTopMarkers` step 5 already does it) and must also hold server-side because the standalone viewer doesn't go through `selectTopMarkers`.

### `raw` cap policy

Don't auto-switch modes silently. If `mode=raw` and `n_sites > max_sites`, return:

```json
{
  "error": "raw_cap_exceeded",
  "n_sites_total": 12453,
  "max_sites": 1000,
  "suggestion": "Use mode=even, mode=hybrid, or mode=aggregate, or raise max_sites (HARD_CAP=20000)."
}
```

with HTTP 400. Forcing the user to choose keeps the discipline rule "raw means all" honest.

### Hard caps

```
DEFAULT_MAX_SITES   = 1000
DEFAULT_HARD_CAP    = 20000        # absolute ceiling regardless of mode
DEFAULT_MAX_REGION_BP = 50_000_000 # > 50 Mb requires mode=aggregate
```

Document in `manifest.json` so frontends can render the right warning.

---

## Endpoints

All endpoints under one FastAPI app, served on `127.0.0.1:8767` (different port from popstats `8766` to keep deployments independent — the user already runs the popstats server).

### `GET /api/manifest`

Returns the contents of `manifest.json` plus runtime info:

```json
{
  "schema_version": "1.0",
  "n_samples": 226,
  "chroms": [
    { "name": "C_gar_LG01", "n_sites": 412034, "length_bp": 41203400 },
    ...
  ],
  "samples": ["CGA_001", ...],
  "candidate_regions_loaded": true,
  "sample_order_loaded": false,
  "limits": { "default_max_sites": 1000, "hard_cap": 20000, "max_region_bp": 50000000 }
}
```

### `GET /api/region`

The general-purpose endpoint. Used by the standalone viewer.

```
GET /api/region?
    chrom=C_gar_LG28
   &start=15000000
   &end=18300000
   &max_sites=1000
   &mode=hybrid
   &seed=1
```

Response (rich envelope, used by viewer):

```json
{
  "chrom": "C_gar_LG28",
  "start_bp": 15000000,
  "end_bp": 18300000,
  "sampling_mode": "hybrid",
  "n_sites_total": 12453,
  "n_sites_returned": 1000,
  "downsampled": true,
  "warning": "hybrid mode includes 30% top-variance markers — biased toward informative sites",
  "samples": ["CGA_001", ...],
  "site_ids": ["LG28:15010234:A>G", ...],
  "positions": [15010234, 15013011, ...],
  "missingness": [0.012, 0.030, ...],
  "diagnostic_score": [null, null, ...],
  "matrix": [[0,1,2,-1,...], ...],          // n_sites × n_samples, int8 with -1 = NA
  "_meta": { "ms_io": 312, "ms_select": 14, "row_groups_read": 4 }
}
```

### `GET /api/dosage/chunk` — atlas bridge endpoint

Translates the rich envelope into the **renderer-compatible chunk shape** documented above. This is the only endpoint the atlas needs.

```
GET /api/dosage/chunk?
    chrom=C_gar_LG28
   &start=15000000
   &end=18300000
   &cap=500          # mapped to max_sites
   &mode=raw         # default for the bridge; renderer does its own selection via selectTopMarkers
```

Response: the chunk shape (`chrom`, `start_bp`, `end_bp`, `samples` (string IDs), `markers[]` with `pos_bp` / `missingness` / `diagnostic_score`, `dosage[][]` with `-1` for NA).

The bridge does **not** invent its own top-N logic — at `mode=raw` it returns up to `cap` sites in positional order and lets `selectTopMarkers` do the rest. This preserves the existing diagnostic-score / variance fallback chain.

If `cap` is too small to fit all sites, the bridge falls back to `mode=even` (because the renderer assumes positional coverage in the chunk's `[start_bp, end_bp]` window) and adds `_meta.fallback_mode: "even"` so the atlas can flag it.

### `GET /api/candidate/{candidate_id}`

Wraps `/api/region` using `candidate_regions.tsv`. Same envelope.

```
GET /api/candidate/INV_LG28_001?max_sites=1500&mode=variance
```

### `GET /api/breakpoint/{candidate_id}/{side}`

`side ∈ {left, right}`. Returns ±`window` bp around the breakpoint.

```
GET /api/breakpoint/INV_LG28_001/left?window=500000&max_sites=1000&mode=even
```

### `GET /api/tile`

Whole-chromosome overview. Builds a `width × height` grid by `aggregate` mode (forced) — no other mode makes sense at this resolution.

```
GET /api/tile?chrom=C_gar_LG28&width=500&height=226
```

Returns the same envelope as `/api/region` with `sampling_mode: "aggregate"` and `n_sites_returned == width`.

`zoom` and `tile` parameters from the original prompt are reserved for a future tile-server (Mapbox-style) — out of scope for v1.

### Atlas chunk-index injection

The bridge alone isn't enough — the atlas reads `state.data.dosage_chunks.chunks` first. Two options (already analyzed in P4.1, but stated here for completeness):

- **Option A — synthesized index** (preferred): when the atlas detects the dosage server is reachable, inject one synthetic chunk per chromosome with a templated URL containing `__START__` / `__END__` / `__CAP__`. The chunk fetcher substitutes these before `fetch()`. The renderer is none the wiser.
- **Option B — direct fetch from renderer**: bypass the index entirely. Worse because it forks the codepath (JSON-loaded chunks vs server chunks behave differently).

Pick A. Implement in a new function `_atlasInstallServerDosageBridge()` called after the server status badge flips to `available`.

---

## Sample order

- **Default**: order from `dosage_store/samples.tsv` (= column order in the parquet).
- **Override**: optional `config/sample_order.tsv` (one ID per line). If present, the server reorders the returned matrix columns to match.
- **Per-request override**: `?sample_order=group_by_pc1` etc. — out of scope for v1; only built-in orderings later.

The **atlas bridge endpoint always returns samples in the canonical store order** so cohort-index lookups in `_buildSampleLookups` stay deterministic. Custom orderings are viewer-only.

---

## File layout

### Atlas bridge (P4.1 — minimum slice)

Lives inside the popstats server (no new server process):

```
Atlas/server_turn1/
  popstats_server.py          # add /api/dosage/chunk handler
  test_dosage_endpoint.py     # new test file
```

Reads `04_dosage_by_chr/` directly with the streaming gzip parser. **No preprocessing required for the bridge.** The full-fat preprocessing is for the standalone viewer.

This is what unblocks the atlas panel today. Everything below (`dosage_viewer/`) is the future work.

### Standalone viewer (future)

```
dosage_viewer/
  README.md
  01_prepare_dosage_store.py        # TSV.GZ → Parquet store
  02_run_server.py                  # FastAPI server, all endpoints
  config/
    sample_order.tsv                # optional
    candidate_regions.tsv           # optional, schema below
  static/
    index.html                      # minimal viewer UI
    app.js                          # canvas/WebGL renderer + fetch
    style.css
  output/
    dosage_store/                   # Parquet store (built by step 01)
    logs/
  tests/
    test_prepare_dosage_store.py
    test_server_endpoints.py
    test_sampling_modes.py
    SIM_test_data/                  # tiny simulated dataset, prefix SIM_
```

**`candidate_regions.tsv`** schema:

| col | type | content |
|---|---|---|
| `candidate_id` | str | e.g. `INV_LG28_001` |
| `chrom` | str | matches store |
| `start` | int | 1-based bp |
| `end` | int | 1-based bp |
| `left_breakpoint` | int | optional, bp |
| `right_breakpoint` | int | optional, bp |
| `notes` | str | free text |

---

## Preprocessing — `01_prepare_dosage_store.py`

Required behavior:

1. Read `<chrom>.sites.tsv.gz` and `<chrom>.dosage.tsv.gz` for each chromosome listed (or auto-discover from a base path).
2. **Detect dosage orientation**: if line 1 has the same count of fields as `samples.tsv` and parses to integers/NA, treat as headerless sites×samples (the canonical case). If line 1 has string values that match a known sample-ID pattern, it's a header — capture sample IDs and skip line 1.
3. **Validate**:
   - `n_rows(dosage) == n_rows(sites)` for each chrom.
   - All cell values ∈ `{0, 1, 2, NA}` (configurable).
   - `pos` strictly increasing in `sites.tsv.gz` (sort if not).
4. Convert dosage cells to `int8` with `-1` = NA.
5. Write `dosage_store/<chrom>/sites.parquet` with row group size 50,000.
6. Write `dosage_store/<chrom>/dosage.parquet` with row group size 50,000, columns `pos, s0, s1, ..., s{n-1}`.
7. Write `dosage_store/samples.tsv` (one ID per line, canonical order).
8. Write `dosage_store/manifest.json`.

CLI:

```
python 01_prepare_dosage_store.py \
    --base /mnt/e/01-catfish_assembly_manuscript_CGA/04_dosage_by_chr \
    --out  dosage_viewer/output/dosage_store \
    --samples /mnt/e/01-catfish_assembly_manuscript_CGA/samples.tsv \
    --chroms C_gar_LG01,C_gar_LG02,...,C_gar_LG28
```

Idempotent: skip a chrom if `<chrom>/manifest.json` exists and source mtime is older.

---

## Frontend — minimal UI (standalone viewer only; not the atlas)

Single-page `static/index.html`:

- top bar: chrom selector (from `/api/manifest`), start, end, max_sites, mode dropdown, seed, **Render** button
- candidate selector (populated from `/api/manifest.candidate_regions_loaded` if present)
- below: a 1024 × N canvas (N = `n_samples`); rows = samples, cols = sites
- metadata strip: `n_sites_total`, `n_sites_returned`, `downsampled`, `mode`, `warning`
- legend strip: dosage 0 / 1 / 2 / NA color squares
- no SVG, no `<table>`, no React for v1 — vanilla JS + Canvas2D

Color scheme: same as atlas (FIG_C08-style). Document the exact hex values in `static/app.js` so they don't drift from the atlas.

WebGL is out of scope for v1. Canvas2D handles 1000 × 226 = 226,000 cells without breaking a sweat.

### Performance rules (frontend)

- one canvas, redraw on Render only
- never request more than `manifest.limits.hard_cap` sites
- recent-region cache via `Map<query_key, response>` with LRU size 16
- no auto-rerender on slider drag — Render-button-only until/unless we add debouncing

---

## Validation against contract (renderer-side)

After implementing the bridge:

1. `state.data.dosage_chunks` is `{ chrom: 'C_gar_LG28', chunks: [{ start_bp: 0, end_bp: 41203400, url: '...?__START__=...' }], source: 'atlas_server' }`.
2. `_findCoveringChunk({ start_bp: 17_000_000, end_bp: 17_300_000 })` resolves to that synthetic chunk; URL placeholders are substituted in `_resolveChunkUrl` before `fetch()`.
3. Fetched chunk has `markers.length === dosage.length`.
4. `chunk.dosage[i].length === chunk.samples.length` for every `i`.
5. Missing values are `-1`, not `null` (so `selectTopMarkers`' `v < 0` filter catches them).
6. Cohort-side lookup: every `chunk.samples[i]` matches one `state.data.samples[j].id` (or `.cga` / `.ind` / `.sample`). Spot-check 10 samples manually first time.
7. Heatmap renders at `~17 Mb` LG28 candidate: 226 rows × ≤500 cols, three bands visible (REF / HET / INV).

If any of those fail, the contract is wrong, not the renderer. Fix the bridge.

---

## Test plan

### Bridge (P4.1)

- `test_dosage_endpoint.py::test_dosage_lg28_window` — 226 × ≤500, cell types = `{0,1,2,-1}`.
- `test_dosage_endpoint.py::test_dosage_empty_region` — `start=1, end=2` returns `n_sites: 0`.
- `test_dosage_endpoint.py::test_dosage_chunk_shape_roundtrip` — fetch chunk, run mock `selectTopMarkers` against it, assert `selected_indices.length > 0` for a known-active region.
- `test_dosage_endpoint.py::test_missing_encoding_negative_one` — assert no `null` cells in matrix.

### Standalone viewer (future)

- `test_prepare_dosage_store.py` — synthetic 100-site × 10-sample TSV.GZ → Parquet round-trip; values preserved; NAs → `-1`.
- `test_sampling_modes.py::test_raw_cap_exceeded_returns_400` — 12k sites, `max_sites=1000`, `mode=raw` → 400 with `error: raw_cap_exceeded`.
- `test_sampling_modes.py::test_even_positional_continuity` — even mode picks roughly equally-spaced positions (KS test against uniform).
- `test_sampling_modes.py::test_random_seed_reproducibility` — same `seed` → same selection on two calls.
- `test_sampling_modes.py::test_variance_top_n` — variance mode returns top-N by per-site variance, then sorted by `pos`.
- `test_sampling_modes.py::test_hybrid_70_30_split` — 700 even + 300 variance after dedup (or close, given dedup eats some).
- `test_sampling_modes.py::test_aggregate_bin_count` — `n_sites_returned == width` regardless of `n_sites_total`.
- `test_sampling_modes.py::test_post_selection_sort_invariant` — for every mode, `positions` array is monotonically increasing.
- `test_server_endpoints.py::test_candidate_endpoint_uses_regions_tsv` — `INV_LG28_001` → coords from candidate_regions.tsv.
- `test_server_endpoints.py::test_breakpoint_left_window` — `±500kb` around `left_breakpoint`.

All tests use `SIM_*` fixtures — the simulated dataset is committed.

---

## Decision points (need user confirmation)

1. **Port for the dosage server.** `8767` (proposed) so popstats `8766` is undisturbed. Or fold into popstats server as a second router (less isolated, fewer processes)? I lean **separate process** for the standalone viewer; the bridge endpoint lives on popstats `8766` either way.

2. **`mode=raw` cap policy** — error vs auto-switch to `even`? I propose **error with suggestion**, because silent auto-switching breaks the discipline ("raw means all") and the rich envelope doesn't have a great way to flag the substitution to a tired user at 2 am.

3. **Storage format final call** — Parquet (proposed) vs HDF5 vs Zarr. Open to override, but I'd want a concrete reason given the pyarrow-already-in-stack argument.

4. **Sample-order config file location** — `dosage_viewer/config/sample_order.tsv` (proposed) vs the existing `Atlas/registries/` location? The viewer is in its own folder, but if the same sample order should drive popstats too, registries makes more sense.

5. **Defaults** — `DEFAULT_MAX_SITES = 1000`, `HARD_CAP = 20000`, `MAX_REGION_BP = 50_000_000`. Comfortable, or too tight / too loose?

---

## What this spec deliberately does NOT include

- WebGL renderer. Canvas2D is fine for 226k cells; revisit only if perf bites.
- Auth. The server is `127.0.0.1`-only.
- Tile pyramid (Mapbox-style multi-zoom). The `tile` endpoint is a single-level overview.
- Sample clustering / PCA-based reordering. Future work.
- Streaming responses (`Transfer-Encoding: chunked`). Responses are <2 MB even at hard cap; no need.
- Compression negotiation. Let nginx / Caddy handle it if anyone ever puts this behind a reverse proxy.
- An `/api/sites` discovery endpoint that returns positions only. The atlas doesn't need it; the viewer can re-query `/api/region` with `mode=aggregate, max_sites=2000` if it ever needs a "where are sites" overview.

---

## Implementation order

1. **P4.1 bridge first** — `popstats_server.py::dosage_chunk()` + atlas-side `_atlasInstallServerDosageBridge()`. Unblocks the candidate panel. 1 session.
2. **Validate** the bridge against three real candidates on LG28 (the ~17 Mb shelf, plus two more). Confirm contract.
3. **Preprocessing script** `01_prepare_dosage_store.py`. 1 session including tests against `SIM_*` fixtures.
4. **Standalone server** `02_run_server.py` with `/api/region`, `/api/manifest`, `/api/candidate`, `/api/breakpoint`. Sampling modes layered in via one shared `select_sites(positions, variances, mode, max_sites, seed)` function. 2 sessions.
5. **`/api/tile`** + frontend `static/`. 1 session.
6. **Hook the bridge to the new store** instead of raw TSV.GZ — once the store exists, the bridge can read parquet too, faster than gzipped TSV. Optional optimization.

The bridge is independent. Step 1 ships value alone if steps 3+ are deferred.
