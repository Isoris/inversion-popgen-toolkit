# Dosage heatmap bridge — operator README

**Status:** P4.1 slice shipped (atlas bridge endpoint + installer). The
full standalone viewer (Parquet store, sampling-modes endpoint, frontend
canvas) is the next slice.

**Spec:** `specs_todo/from_turn129/S6_dosage_heatmap_streaming_viewer.md`.

## What this is

The atlas's chromosome-wide dosage heatmap (`renderDosageHeatmap` in
`Inversion_atlas.html`) reads from `state.data.dosage_chunks.chunks[i].url`.
Without a server, the only way to populate that index is to bake static
JSON chunks per chromosome ahead of time. **The bridge replaces those
static chunks with a live server endpoint** — the user runs the popstats
server, the atlas detects it, and the heatmap reads on-demand.

## Architecture

```
  +------------------------+
  |  renderDosageHeatmap   |   atlas (Inversion_atlas.html, line ~15306)
  |  selectTopMarkers      |
  |  _findCoveringChunk    |   reads state.data.dosage_chunks
  |  _fetchAndCacheChunk   |   substitutes __PLACEHOLDERS__ in templated URLs
  +-----------+------------+
              |
              | fetch(http://localhost:8765/api/dosage/chunk?...)
              v
  +------------------------+
  |  popstats_server.py    |   server_turn1/popstats_server.py
  |  GET /api/dosage/chunk |
  |  GET /api/dosage/manifest
  +-----------+------------+
              |
              | dosage_bridge.py  (handle_dosage_chunk)
              v
  +------------------------+
  |   <base>/04_dosage_by_chr/   |   gzipped TSV, one chrom per pair
  |     C_gar_LG01.sites.tsv.gz  |   pos, ref, alt, missingness, diag, site_id
  |     C_gar_LG01.dosage.tsv.gz |   n_sites x n_samples cells; -1=NA
  |     ...                      |
  +------------------------+
```

Three pieces:

1. **Server side** (`server_turn1/dosage_bridge.py` + popstats wiring)
   - `GET /api/dosage/manifest` — lists available chromosomes + length_bp
   - `GET /api/dosage/chunk?chrom=...&start=...&end=...&cap=...&mode=...`
   - Reads gzipped TSV directly; no preprocessing required for P4.1.
   - Sampling modes: `raw` (default, errors on cap exceeded with structured
     hint) and `even` (positionally-uniform site selection).

2. **Client side** (`js/atlas_dosage_bridge.js`)
   - `AtlasDosageBridge.install({ baseUrl, activeChrom, cap })` probes the
     manifest and writes a synthetic chunk index into `state.data.dosage_chunks`.
   - Each chunk has full-chromosome bounds and a templated URL (one chunk
     per chrom, one URL pattern that gets parameterised at fetch time).
   - `setActiveChrom(chrom)` switches the chrom badge without re-probing.
   - `uninstall()` restores the previous index (static or null).

3. **Two surgical edits to `Inversion_atlas.html`** (counted as one change):
   - `_resolveChunkUrl(chunkRef, region)` substitutes `__CHROM__` /
     `__START__` / `__END__` / `__CAP__` placeholders in templated URLs
     at fetch time. Static URLs pass through unchanged.
   - `_fetchAndCacheChunk(chunkRef, region)` now takes a region and uses
     it for cache keying when the URL is templated (so different regions
     don't collide on the same cached payload). Static URLs keep their
     previous region-agnostic cache key — no behavior change for the
     pre-existing pre-baked-JSON workflow.

## How to set it up (operator)

### 1. On LANTA — confirm the dosage TSVs exist

```
ls $BASE/04_dosage_by_chr/
# Expected: C_gar_LG01.sites.tsv.gz  C_gar_LG01.dosage.tsv.gz  ...
```

If absent, build them with the standard pipeline. The schema is:

- `*.sites.tsv.gz`: `pos<TAB>ref<TAB>alt<TAB>missingness<TAB>diagnostic_score<TAB>site_id`
  (last 3 columns optional; "NA" allowed for missing/diag)
- `*.dosage.tsv.gz`: one row per site, n_samples cells per row, values in
  `{0, 1, 2, NA, ., -1}`. Missing → -1 on the wire.
- An optional first row of sample IDs is auto-detected and skipped.

### 2. Update the popstats config

Add to `popstats_server.config.yaml`:

```yaml
dosage_dir: "${base}/04_dosage_by_chr"
# Optional override; defaults to sample_list:
# dosage_samples: "${base}/samples.tsv"
dosage_max_region_bp: 50000000
```

Restart the popstats server.

### 3. On the atlas side — install the bridge

The bridge installer can be wired into the atlas's existing
`atlasServer` health-check loop. Minimal usage from the dev console:

```js
// Probe + install (uses window.atlasServer.url and the manifest):
AtlasDosageBridge.install().then(r => console.log('bridge:', r));
// -> { installed: true, chroms: 28, baseUrl: 'http://localhost:8765',
//      activeChrom: 'C_gar_LG01' }
```

To wire it into the atlas startup, add a script tag near the
`atlas_sv_evidence.js` include:

```html
<script src="js/atlas_dosage_bridge.js"></script>
<script>
  if (typeof atlasServer !== 'undefined') {
    atlasServer.isAvailable().then(ok => {
      if (ok) AtlasDosageBridge.install({ activeChrom: state.data && state.data.chrom });
    });
  }
</script>
```

### 4. Verify in the browser

Open the dev console and run:

```js
AtlasDosageBridge.getStatus()
// { installed: true, chroms: 28, baseUrl: '...', activeChrom: '...', ts: ... }

state.data.dosage_chunks
// { source: 'atlas_server', chunks: [...], chrom: 'C_gar_LG28', ... }
```

Then trigger the candidate-panel dosage heatmap. The renderer's
`_findCoveringChunk` resolves the active candidate region to the synthetic
chrom-wide chunk, `_fetchAndCacheChunk` substitutes placeholders to the
real region bounds, the server returns a region-scoped chunk, the renderer
draws.

## Testing

```
# Server side (pure unit + integration via synthesized fixtures)
python3 tests/dosage_bridge/test_dosage_bridge.py             # 79 tests

# Atlas-side installer (Node, no DOM required)
node tests/dosage_bridge/test_atlas_dosage_bridge.js          # 62 tests

# Surgical-edit verification (extracts edited functions from the HTML
# and exercises them in a sandbox; catches regressions in the URL resolver
# and the cache-key change)
node tests/dosage_bridge/test_html_surgical_edits.js          # 29 tests
```

## What's deferred (next slice)

The full viewer spec lays out 4 more slices:

1. **`01_prepare_dosage_store.py`** — TSV.GZ → Parquet store (faster random
   access; the bridge currently reads gzipped TSV linearly).
2. **`02_run_server.py`** — standalone FastAPI app on port 8767 with the
   richer `/api/region` endpoint (6 sampling modes: raw, even, random,
   variance, hybrid, aggregate) plus `/api/candidate`, `/api/breakpoint`,
   `/api/tile`.
3. **`static/`** — the standalone HTML viewer (canvas2D, chrom selector,
   region inputs, sampling-mode dropdown).
4. **Optional bridge-store hookup** — point P4.1's bridge at the parquet
   store instead of raw TSV.GZ for speed.

These can ship independently. The P4.1 slice already unblocks the candidate
panel; the rest is performance + standalone-viewer work.

## Constants you may want to tweak

| Constant | Default | Where | Notes |
|----------|---------|-------|-------|
| `DEFAULT_MAX_SITES` | 1000 | `dosage_bridge.py` | per-request site cap |
| `DEFAULT_HARD_CAP` | 20000 | `dosage_bridge.py` | absolute ceiling on `cap=` |
| `DEFAULT_MAX_REGION_BP` | 50_000_000 | `dosage_bridge.py` | rejects chrom-scale requests |
| `DEFAULT_BASE_URL` | `http://localhost:8765` | `atlas_dosage_bridge.js` | popstats server URL |
| `DEFAULT_CAP_DEFAULT` | 1000 | `atlas_dosage_bridge.js` | falls back into __CAP__ when not overridden |
| `DEFAULT_TIMEOUT_MS` | 1500 | `atlas_dosage_bridge.js` | manifest probe timeout |

## File map

```
server_turn1/
  dosage_bridge.py                              ← handler module (no fastapi-app of its own)
  popstats_server.py                            ← imports + wires GET /api/dosage/{chunk,manifest}
  popstats_server.config.example.yaml           ← documents the dosage_dir / dosage_samples keys

js/
  atlas_dosage_bridge.js                        ← installer + URL builder + setActiveChrom

Inversion_atlas.html                            ← +27 lines: _resolveChunkUrl + region in cache key

tests/dosage_bridge/
  test_dosage_bridge.py                         ← 79 tests; gz fixtures + handler round-trip
  test_atlas_dosage_bridge.js                   ← 62 tests; install + setActiveChrom + uninstall
  test_html_surgical_edits.js                   ← 29 tests; verifies the HTML edit didn't break static

specs_todo/from_turn129/
  S6_dosage_heatmap_streaming_viewer.md         ← full viewer spec (this slice = §P4.1 only)
```
