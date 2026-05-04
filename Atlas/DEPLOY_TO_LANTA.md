# DEPLOY_TO_LANTA.md — server-side deployment for the Inversion Atlas

End-to-end deployment guide for the live-popstats server (turn 1) +
fast_ld LD endpoint (turn 11c). Run this on LANTA.

## What's in this tarball

```
server_bundle/
├── server_turn1/                  Live-popstats FastAPI server
│   ├── popstats_server.py         The server (~1670 lines)
│   ├── popstats_server.config.example.yaml
│   ├── requirements.txt           fastapi, pydantic, pyyaml, numpy, pandas
│   ├── test_units.py              14 unit tests
│   ├── test_with_curl.sh          10 curl integration tests
│   └── README.md
│
├── server_turn11c_ld_fast/        fast_ld LD endpoint
│   ├── fast_ld_endpoint.py        FastAPI handler
│   ├── lazy_windows_json.py       Per-chrom windows.json cache
│   ├── test_fast_ld_endpoint.py   11 endpoint tests
│   └── SPLICE_POINTS.md           Source of truth for the patch
│
└── engine_fast_ld/                Standalone fast LD C engine
    ├── fast_ld.c                  ~1100 LOC C source
    ├── Makefile                   make / make test / make static
    ├── fast_ld_wrapper.py         Python adapter (compute_ld(req))
    ├── build_windows_json.py      Preprocessor: atlas JSON → windows JSON
    ├── test_fast_ld.py            8 engine tests
    ├── test_wrapper.py            2 wrapper tests
    ├── bench.py                   Realistic-scale benchmark
    └── README.md
```

The atlas-side bundle (HTML + JS modules) ships separately as
`inversion_atlas_bundle.tar.gz`.

## Step 1 — Build the fast_ld engine

```bash
cd engine_fast_ld
make            # builds ./fast_ld
make test       # 10 tests should pass (8 engine + 2 wrapper)
```

Toolchain: `gcc` (or `clang`), `zlib-dev`, OpenMP (`libgomp` on
Ubuntu, `libomp` on macOS). On LANTA these are in the default module
environment.

For a portable, statically-linked binary suitable for distribution:

```bash
make static
```

Copy the binary to your LANTA tools directory. Recommended location:

```
/scratch/lt200308-agbsci/.../tools/fast_ld
```

Save the absolute path — you need it in step 4.

## Step 2 — Stage the Python files

Create `popstats_server/` deployment dir and populate:

```bash
mkdir -p popstats_server
cp server_turn1/popstats_server.py                  popstats_server/
cp server_turn1/popstats_server.config.example.yaml popstats_server/
cp server_turn1/requirements.txt                    popstats_server/
cp server_turn1/test_units.py                       popstats_server/

cp server_turn11c_ld_fast/fast_ld_endpoint.py       popstats_server/
cp server_turn11c_ld_fast/lazy_windows_json.py      popstats_server/
cp server_turn11c_ld_fast/test_fast_ld_endpoint.py  popstats_server/

cp engine_fast_ld/fast_ld_wrapper.py                popstats_server/
cp engine_fast_ld/build_windows_json.py             popstats_server/
```

Install Python deps (project uses pip; all deps pip-installable):

```bash
cd popstats_server
pip install --user -r requirements.txt
```

## Step 3 — Patch popstats_server.py

Two edits. Both are surgical — see
`server_turn11c_ld_fast/SPLICE_POINTS.md` for the canonical version.

**Edit A — imports near the top of `popstats_server.py`:**

```python
from pathlib import Path
from fast_ld_endpoint import FastLDReq, handle_split_heatmap, fast_ld_engine_hash
from lazy_windows_json import WindowsJsonCache
```

**Edit B — server bootstrap + route handler, near the other endpoint
mounts:**

```python
# Bootstrap: next to where engines are loaded
WINDOWS_CACHE = WindowsJsonCache(
    cache_dir=Path(CFG["windows_cache_dir"]),
    atlas_json_dir=Path(CFG["atlas_json_dir"]),
    sites_dir=Path(CFG["sites_dir"]),
    builder_script=Path(CFG["build_windows_json_path"]),
)
FAST_LD_BIN  = Path(ENGINES.path("fast_ld"))
FAST_LD_HASH = fast_ld_engine_hash(FAST_LD_BIN)

@app.post("/api/ld/split_heatmap")
async def ld_split_heatmap(req: FastLDReq) -> JSONResponse:
    _ensure_ready()
    def cache_get(k):    return CACHE.get_json(k, "ld")
    def cache_put(k, v): CACHE.put_json(k, "ld", v)
    payload = await handle_split_heatmap(
        req,
        fast_ld_bin=FAST_LD_BIN,
        fast_ld_engine_hash=FAST_LD_HASH,
        dosage_dir=Path(CFG["dosage_dir"]),
        windows_cache=WINDOWS_CACHE,
        sample_filter=SAMPLES.filter_known,
        cache_get=cache_get, cache_put=cache_put,
    )
    return JSONResponse(payload)
```

This is the canonical `/api/ld/split_heatmap` route. If you also have
the older turn-11a ngsLD-based route mounted and you want both for
side-by-side validation, rename one — e.g.
`/api/ld/split_heatmap_ngsld` for the old one, leave
`/api/ld/split_heatmap` for fast_ld.

## Step 4 — Configure popstats_server_config.yaml

Copy the example and edit paths:

```bash
cd popstats_server
cp popstats_server.config.example.yaml popstats_server_config.yaml
```

Add the following **new** keys (existing turn-1 keys remain
untouched):

```yaml
# === Turn 11c (fast_ld) — new keys ===
dosage_dir:              /scratch/lt200308-agbsci/.../dosage_per_chrom
atlas_json_dir:          /scratch/lt200308-agbsci/.../atlas_data
sites_dir:               /scratch/lt200308-agbsci/.../sites_per_chrom
windows_cache_dir:       /scratch/lt200308-agbsci/.../cache/windows_json
build_windows_json_path: /scratch/lt200308-agbsci/.../tools/build_windows_json.py

# === Engines block — add fast_ld; keep existing entries ===
engines:
  region_popstats: <existing>
  hobs_windower:   <existing>
  angsd_patched:   <existing>
  instant_q:       <existing>
  fast_ld:         /scratch/lt200308-agbsci/.../tools/fast_ld
```

Notes:
- `dosage_dir` should contain `<chrom>.dosage.tsv.gz` per chromosome
  (BEAGLE-imputed; same files the atlas reads).
- `atlas_json_dir` is the directory holding the atlas chromosome
  JSON (e.g. `LG28_v8.json`).
- `sites_dir` should contain `<chrom>.sites.tsv` (or `.tsv.gz`) —
  per-SNP `chrom\tpos\tmaj\tmin` for the same SNPs as the dosage
  matrix.
- `windows_cache_dir` is writable scratch — created if missing.
- `build_windows_json_path` is the absolute path to the
  `build_windows_json.py` you copied into `popstats_server/` (or
  wherever you stage it).

## Step 5 — Run unit tests

```bash
cd popstats_server
python3 test_units.py                  # 14 tests (turn 1)
python3 test_fast_ld_endpoint.py       # 11 tests (turn 11c, real binary E2E)
```

Both must be green before serving.

## Step 6 — Start the server

```bash
cd popstats_server
python3 -m popstats_server
# binds 127.0.0.1:8765 by default
```

(Run under `tmux` / `screen` / `systemd-run --user` as appropriate
for your LANTA login session.)

## Step 7 — Verify

Health check — confirms fast_ld is registered with a content hash:

```bash
curl -s http://127.0.0.1:8765/api/health | jq '.engines.fast_ld'
# → "abcd1234ef..."   (~16-char content hash)
```

Run a real LD call (replace the chrom + sample lists with valid ones
for your cohort):

```bash
curl -X POST http://127.0.0.1:8765/api/ld/split_heatmap \
  -H 'Content-Type: application/json' \
  -d '{
    "chrom": "C_gar_LG28",
    "window_range": [2150, 2200],
    "groups": {
      "ALL":     ["S001", "S002", "..."],
      "HOM_INV": ["S017", "S041", "..."]
    },
    "shelf_bp": [15200000, 18100000]
  }' | jq '.timing, .cache_state, .matrices.HOM_INV.shelf_ratio'
```

First call: builds `LG28.windows.json` (~1 s for LG28-sized chroms),
runs fast_ld for both groups, returns ~`0.4 s` total. `cache_state`
will be `"miss"`.

Second identical call: `cache_state` will be `"hit"`, total wallclock
should be sub-millisecond.

## Step 8 — Point the atlas at the server

The atlas-side `atlas_ld.js` defaults to `POST /api/ld/split_heatmap`
on the same origin as the page. If you serve the atlas HTML from a
different host (e.g. open it locally and the server runs on
`http://lanta.example:8765`), set up CORS or use a reverse proxy.

Same applies to all other live-popstats endpoints
(`/api/popstats/groupwise`, `/api/popstats/hobs_groupwise`,
`/api/ancestry/groupwise_q`, `/api/popstats/hwe`).

## Smoke-test checklist after deployment

```
[ ] make + make test in engine_fast_ld/  → 10 tests green
[ ] python3 test_units.py                → 14 tests green
[ ] python3 test_fast_ld_endpoint.py     → 11 tests green
[ ] /api/health returns engines.fast_ld content hash
[ ] First /api/ld/split_heatmap call: cache_state="miss", < 1 s
[ ] Second identical call: cache_state="hit", < 50 ms
[ ] Atlas page 3 LD panel renders heatmap on ▶ run with focal candidate
```

## Troubleshooting

**"engines.fast_ld is null"** — the `fast_ld` path in the config
points to a missing or non-executable file. `chmod +x`, verify path.

**"window_range out of bounds"** — the request's `[start_w, end_w]`
exceeds the chromosome's window count in the atlas JSON. Check that
`atlas_json_dir/<chrom>.json` has at least `end_w + 1` windows.

**Slow first call (>5 s)** — the lazy windows-JSON build is reading
a very large sites file. Consider pre-building (run
`build_windows_json.py` once per chrom and place output in
`windows_cache_dir`).

**Cache hits don't appear** — the cache key includes the engine
content hash, so any rebuild of `fast_ld` invalidates all keys. This
is by design.

**Atlas LD panel shows "error: fetch ..."** — server unreachable
from the browser. Test with the curl verification command first; if
that works, it's a CORS / origin issue between the HTML and the
server.

## What this deployment does NOT do

- Validate fast_ld against ngsLD on real LANTA data — that's a
  separate cross-check task. Math has been validated against
  `numpy.corrcoef²` on synthetic data (within 1/255 quantization
  across 1.4M pairs). For real-data validation, mount the older
  turn-11a `ld_endpoint.py` at `/api/ld/split_heatmap_ngsld` in
  parallel and compare shelf-ratio values for the same candidate.
- Provide whole-chromosome live LD — endpoint hard-caps
  `window_range` to keep latency reasonable. Whole-chromosome views
  use the pre-baked `STEP_LD_00b` SLURM cache.
- GHSL-stratified live mode (parallel server endpoint over phased
  haplotypes). Tracked as future work.
