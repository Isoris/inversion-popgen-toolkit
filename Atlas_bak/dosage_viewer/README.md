# Dosage viewer — standalone canvas viewer

**Spec:** `specs_todo/from_turn129/S6_dosage_heatmap_streaming_viewer.md`

This is the **standalone viewer** — separate from the atlas-side P4.1 bridge.
The bridge (`server_turn1/dosage_bridge.py` + `js/atlas_dosage_bridge.js`)
hooks into the atlas's existing dosage heatmap renderer; the viewer is a
self-contained webapp for region exploration with all six sampling modes.

## Pieces

```
dosage_viewer/
  01_prepare_dosage_store.py    # TSV.GZ → Parquet store
  02_run_server.py              # FastAPI server (port 8767)
  region_handler.py             # endpoint orchestrators
  sampling.py                   # 6 sampling modes (raw/even/random/variance/hybrid/aggregate)
  store_reader.py               # parquet read helpers
  static/
    index.html                  # top bar, canvas, legend, recent-history
    style.css                   # FIG_C08-aware light + dark theme
    app.js                      # bootstrap + render flow + LRU + hover

tests/dosage_viewer/
  test_prepare_dosage_store.py  # 40 tests
  test_sampling_modes.py        # 68 tests
  test_store_reader.py          # 32 tests
  test_region_handler.py        # 70 tests
  test_run_server.py            # 61 tests (FastAPI TestClient)
  test_frontend_pure.js         # 21 tests (colour ramp, LRU)
                                # = 292 viewer tests total
```

## End-to-end setup

### 1. Build the parquet store

```bash
python3 dosage_viewer/01_prepare_dosage_store.py \
    --base    /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/04_dosage_by_chr \
    --out     /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/dosage_store \
    --samples /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/het_roh/01_inputs_check/samples.ind
```

Idempotent — re-run safely; only chroms whose source mtime changed since
the last build are rebuilt. Use `--force` to rebuild everything.

Outputs:

```
dosage_store/
  manifest.json                 # store-level: schema_version, n_samples, chroms[]
  samples.tsv                   # canonical sample order
  C_gar_LG01/
    sites.parquet               # pos, ref, alt, missingness, diagnostic_score, site_id
    dosage.parquet              # pos + 226 int8 sample columns; -1 = NA
    _chrom_manifest.json        # per-chrom mtime stamps for skip-if-done
  C_gar_LG02/
    ...
```

Storage cost: at 226 samples × ~5 M sites × int8 = **~1.1 GB per chrom**
uncompressed (~250 MB compressed with zstd row groups). Sites parquet is
~6 cols × 10 bytes/site = ~50 MB/chrom uncompressed (~10 MB compressed).

### 2. Optional — write a candidate config

```
config/candidate_regions.tsv:
candidate_id    chrom         start     end       left_breakpoint    right_breakpoint    notes
INV_LG28_001    C_gar_LG28    14000000  18000000  15750000           17890000            shelf inversion
INV_LG12_001    C_gar_LG12    11200000  12100000                                         small candidate
```

### 3. Start the server

```bash
python3 dosage_viewer/02_run_server.py \
    --store      /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/dosage_store \
    --candidates /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/config/candidate_regions.tsv \
    --static     dosage_viewer/static \
    --port       8767
```

Open `http://localhost:8767/` in a browser.

## Sampling modes

| mode | what | when to use |
|------|------|-------------|
| `raw` | every site in window | default for small windows; errors with `raw_cap_exceeded` if n > max_sites |
| `even` | bp-uniform: bin into max_sites bins, pick closest site to each midpoint | the neutral large-window view |
| `random` | reproducible random subset (seed) | sanity-check that a pattern isn't an artefact of a particular site selection |
| `variance` | top max_sites by per-site dosage variance | "where's the signal" — biased toward informative sites |
| `hybrid` | 70% even + 30% variance, deduped | DEFAULT for arbitrary windows; balances coverage with informativeness |
| `aggregate` | bin into n_bins, mean dosage per bin per sample | chrom-scale tile view; the only mode that makes sense at >50 Mb |

All non-`raw` modes return positions sorted ascending. Aggregate matrices
have `null` cells for empty bins (vs `-1` for per-site NA cells in other
modes).

## API quick reference

| endpoint | params | returns |
|----------|--------|---------|
| `GET /api/manifest` | — | store metadata + chrom list |
| `GET /api/region` | chrom, start, end, max_sites=1000, mode=hybrid, seed=1 | rich envelope with `matrix` |
| `GET /api/dosage/chunk` | chrom, start, end, cap=1000, mode=raw | renderer-compat shape (markers + dosage) |
| `GET /api/candidate/{id}` | max_sites, mode, seed | wraps /api/region using candidate_regions.tsv |
| `GET /api/breakpoint/{id}/{side}` | window_bp=500000, max_sites, mode, seed | ± window around L/R breakpoint |
| `GET /api/tile` | chrom, width=500 | whole-chrom aggregate overview |
| `GET /api/health` | — | service status |

Error envelope (status 400/404):

```json
{
  "detail": {
    "error": "raw_cap_exceeded",
    "n_sites_total": 53129,
    "max_sites": 1000,
    "suggestion": "Use mode=even/random/hybrid/variance/aggregate, or raise max_sites."
  }
}
```

## Frontend behaviour notes

- **No auto-rerender on slider drag.** The Render button is the only
  thing that triggers a fetch. Spec §UX explicitly forbids debounced-
  drag fetching because chrom-scale aggregate calls are too expensive.
- **LRU recent-region cap = 16.** Click any entry to refill the inputs
  and re-render. Cache key includes the mode (so `raw` and `even` of the
  same window are separate entries).
- **Theme toggle** swaps light/dark via `[data-theme]` on `<body>`.
  Canvas is redrawn (NA pattern colours change with theme).
- **Hover** on a tile shows `chrom:pos · sample · dosage = [0|1|2|NA]`.
  Aggregate cells show `dosage = 1.34` (2 decimal places).
- **Escape** aborts an in-flight fetch.

## Test commands

```bash
# Pure unit tests
python3 tests/dosage_viewer/test_sampling_modes.py        # 68
python3 tests/dosage_viewer/test_prepare_dosage_store.py  # 40
python3 tests/dosage_viewer/test_store_reader.py          # 32
python3 tests/dosage_viewer/test_region_handler.py        # 70

# FastAPI integration (uses TestClient — no real socket)
python3 tests/dosage_viewer/test_run_server.py            # 61

# Frontend pure logic (extracts functions from app.js, runs in vm sandbox)
node tests/dosage_viewer/test_frontend_pure.js            # 21
```

## What's NOT in this slice

- **No server-side caching.** Each request re-reads the parquet. The atlas-
  side bridge installs a per-region LRU on the client; the standalone
  viewer's frontend has the same. If server-side caching is needed, mirror
  popstats_server's `DiskCache` class — keyed by `(chrom, start, end, max_sites, mode, seed)`.
- **No predicate-pushdown on parquet reads.** `read_sites_in_window` reads
  the whole sites column then slices. Move to `pq.read_table(..., filters=...)`
  if region scans on chrom-scale parquets get slow (they shouldn't for
  typical 1–5 Mb windows).
- **No /api/candidates list endpoint.** Frontend can't populate the
  candidate dropdown automatically yet — spec defers that to the next
  slice. For now, the operator types `INV_LG28_001` etc. into the
  candidate select, OR uses the URL like `/api/candidate/INV_LG28_001`
  directly. (Easy add: emit `candidate_ids` array in /api/manifest when
  candidates are loaded.)

## File map

```
server_turn1/dosage_bridge.py                      # P4.1 — atlas-side bridge handler
server_turn1/popstats_server.py                    # +/api/dosage/chunk + /api/dosage/manifest
js/atlas_dosage_bridge.js                          # P4.1 — atlas installer
Inversion_atlas.html                               # P4.1 — surgical _resolveChunkUrl edit

dosage_viewer/01_prepare_dosage_store.py           # P4 — Parquet builder
dosage_viewer/02_run_server.py                     # P4 — standalone FastAPI on port 8767
dosage_viewer/sampling.py                          # P4 — 6 sampling modes
dosage_viewer/store_reader.py                      # P4 — Parquet read helpers
dosage_viewer/region_handler.py                    # P4 — endpoint orchestrators
dosage_viewer/static/index.html                    # P4 — frontend
dosage_viewer/static/style.css                     # P4 — CSS-variable theme
dosage_viewer/static/app.js                        # P4 — viewer logic

tests/dosage_bridge/                               # P4.1 tests (170 assertions)
tests/dosage_viewer/                               # P4 tests (292 assertions)
```
