# popstats_server

Live-engine HTTP wrapper around the inversion popstats engines
(`region_popstats`, `hobs_windower`, `angsd_fixed_HWE`, `instant_q`).
Sits on the cluster (LANTA) where the BEAGLE / BAM / Engine-B caches live;
the atlas (running in the user's browser) reaches it over an SSH tunnel.
Single-file Python service, ~1.5 k lines, 5 dependencies.

## What it is for

The atlas does most popgen analysis live in the browser using precomputed
data layers (dosage chunks, ancestry Q-matrix, per-sample θπ). But three
classes of computation can't run client-side:

- **Hudson Fst / dXY** — quadratic over BEAGLE GLs the browser doesn't have
- **HoverE per group** — needs Mérot's per-group ANGSD HWE re-run
- **θπ per group with proper variance** — same BEAGLE-GL dependency

This server wraps the cluster-side C engines that *do* compute these,
exposes them as HTTP endpoints, and caches results so a re-request of the
same group composition returns in <50 ms instead of re-running the engine.

The atlas decides what a group **is** (by lassoing samples, applying
K-means, intersecting families with regimes, etc.). The server only knows
how to compute statistics given a list of sample IDs per group. This
separation is deliberate: the atlas can ship rich group-composition UI
without requiring server changes, and the server stays a stateless
engine-runner.

## Endpoints

| Method | Path                              | Purpose |
|--------|-----------------------------------|---------|
| GET    | `/api/health`                     | engine versions + cache stats |
| POST   | `/api/popstats/groupwise`         | wraps `region_popstats` (Engine F) |
| POST   | `/api/popstats/hobs_groupwise`    | wraps `angsd_fixed_HWE` + `hobs_windower` |
| POST   | `/api/ancestry/groupwise_q`       | reads `instant_q` cache, aggregates per group |
| POST   | `/api/shelf_ld_test`              | 501 — atlas does this in JS |
| GET    | `/api/cache/keys`                 | debug: list cache hashes |
| DELETE | `/api/cache/keys/<hash>`          | debug: drop one entry |
| GET    | `/api/jobs/<id>`                  | (reserved for future async runs) |

All `POST` bodies are JSON. Group definitions are explicit member lists:

```json
{
  "chrom": "C_gar_LG28",
  "region": { "start_bp": 15000000, "end_bp": 18000000 },
  "groups": {
    "HOM1": ["CGA_001", "CGA_007", ...],
    "HET":  ["CGA_002", "CGA_011", ...],
    "HOM2": ["CGA_003", "CGA_017", ...]
  },
  "metrics": ["fst", "dxy", "theta_pi"],
  "win_bp": 50000, "step_bp": 10000
}
```

Group names match the regex `[A-Za-z0-9_]+`. Member IDs are validated
against the canonical sample list (`SAMPLE_LIST_POPSTATS`); unknowns are
rejected with a clear 400. Up to 10 groups per request (Engine F's
`MAX_GROUPS`).

## Install (cluster side)

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
git clone <atlas-repo>     # or copy the popstats_server/ dir into place
cd popstats_server

# Conda / mamba env (assumes 'assembly' env from 00_ancestry_config.sh)
conda activate assembly
pip install -r requirements.txt

# Compile engines if not already done. Each is independent.
make -C ../inversion-popgen-toolkit/unified_ancestry/engines/fst_dxy
make -C ../inversion-popgen-toolkit/unified_ancestry/engines/hobs_hwe/scripts hobs_windower
# angsd_fixed_HWE is the patched ANGSD binary — see "About the ANGSD patch" below.

# Copy and edit the config
cp popstats_server.config.example.yaml popstats_server.config.yaml
${EDITOR:-nano} popstats_server.config.yaml
```

## Run

```bash
python popstats_server.py --config popstats_server.config.yaml
```

The server logs a startup banner showing the resolved paths and the
engine binary hashes:

```
2026-05-02 12:00:48 [INFO] popstats: config loaded from ...
2026-05-02 12:00:48 [INFO] popstats: cache: loaded 0 entries (0.0 MB) from /scratch/.../popstats_server_cache
2026-05-02 12:00:48 [INFO] popstats: sample_list: 226 ids from /scratch/.../samples.ind
2026-05-02 12:00:48 [INFO] popstats: startup complete: 226 samples, 4 engines, 0 cached entries
2026-05-02 12:00:48 [INFO] popstats: listening on http://127.0.0.1:8765
```

Verify with curl from the same node:

```bash
curl -s http://127.0.0.1:8765/api/health | python3 -m json.tool
```

## Connecting the atlas (SSH tunnel)

The server binds to `127.0.0.1` only. The atlas runs in your local
browser. Bridge them with an SSH local-port-forward:

```bash
# On your laptop
ssh -L 8765:127.0.0.1:8765 lanta
# (leave open in a tmux pane; the tunnel lives as long as ssh does)
```

Then in the browser, the atlas's `atlasServer.url` defaults to
`http://localhost:8765`, which now reaches the cluster. The atlas's
status-badge in the toolbar should flip to green within 30 seconds.

If you're running both an atlas launcher (`atlas.sh`) AND the popstats
server on the same machine, give them different ports. The launcher
defaults to 8765; if you keep that, run popstats_server on **8766** and
update `atlas.localStorage['inversion_atlas.serverUrl']` to
`http://localhost:8766` (or click the badge → "set URL").

## About the ANGSD patch

`angsd_fixed_HWE` is **stock ANGSD with a one-line bugfix from Claire
Mérot's group** that corrects the per-site `F` estimator under the model
where samples within a group are not all in HWE (the case for grouped
inversion karyotypes). The CLI is unchanged from upstream ANGSD —
this server invokes it as a normal subprocess with the standard
`-bam -ref -out -GL 1 -doMajorMinor -doMaf 1 -SNP_pval 1e-6 -doHWE 1
-maxHetFreq 1.0 -minMaf -minMapQ -minQ -r <chr>: -nThreads N` argv.

Compile from source on the cluster:

```bash
cd <patch-source-dir>
make clean && make -j8
cp angsd /scratch/.../angsd_patched/angsd
```

Point the config's `engines.angsd_patched` at the resulting binary.
The server hashes the binary's bytes and includes that hash in every
HWE cache key, so a recompile auto-invalidates downstream cached results.

## Cache layout and invalidation

```
${cache_dir}/
  popstats/<hash>.json    region_popstats results (Engine F)
  hobs/<hash>.json        Q07b+Q07c chained results
  ancestry/<hash>.json    instant_q-cache aggregations
  hwe/<hwe_hash>.hwe.gz   raw ANGSD HWE per-group output
  shelf_ld/               (reserved)
  index.jsonl             append-only log; the filesystem is authoritative
```

Cache keys are content-addressable:

```
sha256(chrom + region + sorted(group_name → sorted(members))
       + metric_set + win + step + binary_hash)[:32]
```

Sorting at every level guarantees that the atlas can send groups in any
UI order with identical hits. Recompiling any engine binary changes its
content hash → all dependent cache entries become unreachable on lookup
(eventually evicted by LRU at the configurable byte cap).

The cache directory is **safe to delete** at any time. The server will
rebuild the index from disk on the next start and silently re-run engines
on cache miss.

## Testing

```bash
# Unit tests (no engines needed, no server boot)
python test_units.py

# Smoke tests (boots a server with stub binaries + synthetic data)
bash test_with_curl.sh
```

Both are dependency-free against the cluster. Real-engine integration
testing happens against the live LANTA installation — see
`docs/lg28_demo.md` (turn 7 of chat A) for the end-to-end walkthrough.

## Operational notes

- **Single-user assumption**: the job manager and ANGSD scheduling
  assume one user at a time. Multi-user scaling (queue with backpressure,
  per-user cache namespaces) is out of scope for chat A.
- **HWE cost**: first request for a novel group composition triggers an
  ANGSD run (~5–15 s per group). Three groups in parallel ≈ same wall
  time. Subsequent requests with the same composition return in <100 ms.
- **`region_popstats` cost**: ~50–500 ms per chromosome depending on
  region restriction. Use the `region` field in requests to keep
  per-request latency interactive.
- **`local_Q_samples.tsv.gz` cost**: pure file I/O + pandas group-by;
  ~50 ms uncached, <10 ms cached.

## License

Same as the inversion-popgen-toolkit umbrella project.
