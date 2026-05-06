# P4.1 — Dosage chunk bridge: server endpoint + atlas wiring

## Risk: medium
## Lines changed: ~120 server + ~40 atlas
## Depends on: P1.5 stable server startup
## Verification: candidate dosage heatmap renders for an LG28 candidate

---

## What you said

> Dosage heatmap status:
> The atlas already has dosage heatmap UI/renderer code, but the
> data bridge is missing/incomplete.
>
> The panel checks `state.data.dosage_chunks`. If missing, it shows:
> "no dosage_chunks index — load enrichment JSON to enable"
>
> Preferred bridge: Add local server endpoint because popstats
> server already runs.
>
> Suggested endpoint:
>   GET /api/dosage/chunk?chrom=C_gar_LG28&start=15000000&end=18300000&cap=500
>
> It should read:
>   ${base}/04_dosage_by_chr/C_gar_LG28.sites.tsv.gz
>   ${base}/04_dosage_by_chr/C_gar_LG28.dosage.tsv.gz
>
> Filter rows by start/end. Rank/thin markers if needed. Return JSON
> chunk for the existing renderer.
>
> Main dosage heatmap meaning:
>   rows = samples/fish
>   columns = SNPs/markers across candidate region
>   color = dosage 0/1/2/missing
>
> Do not collapse this into one band for the main heatmap.

## Reverse-engineer the existing renderer

Anchor: `function renderDosageHeatmap` and `function _wireCandidateDosageHeatmap`
in your atlas. Read what they expect from `state.data.dosage_chunks`.

In my snapshot (turn 128), the renderer reads from
`state.data.dosage_chunks.chunks` — an array of `{chrom, start_bp,
end_bp, url}`. The candidate code calls something like
`fetch(chunk.url)` and expects JSON with:

```jsonc
// the chunk JSON, not the index
{
  "chrom": "C_gar_LG28",
  "start_bp": 15000000,
  "end_bp": 18300000,
  "n_samples": 226,
  "samples": [{"cga": "...", "ind": "..."}, ...],   // sample order
  "n_sites": 487,
  "sites": [
    {"pos": 15010234, "ref": "A", "alt": "G"},
    {"pos": 15013011, "ref": "C", "alt": "T"},
    // ...
  ],
  "dosage_matrix": [
    // n_samples rows × n_sites cols, each cell is 0/1/2 or null/-1 for missing
    [0, 1, 2, 0, ...],
    [1, 1, 2, 1, ...],
    // ...
  ]
}
```

**Verify this shape against your actual renderer before implementing
the endpoint.** Search the renderer for accesses like
`chunk.dosage_matrix` / `chunk.sites[i].pos` / `chunk.samples[si]`
to confirm.

## Server endpoint

In `popstats_server.py`, add:

```python
import gzip
import csv
from typing import Optional

@app.get("/api/dosage/chunk")
async def dosage_chunk(
    chrom: str,
    start: int,
    end: int,
    cap: int = 500,
    sort_by_pos: bool = True,
):
    """Read a dosage matrix slice from disk and return as JSON for
    the atlas dosage heatmap renderer.

    Inputs:
      ${base}/04_dosage_by_chr/<chrom>.sites.tsv.gz
        column 1: pos (bp, 1-based)
        column 2: ref
        column 3: alt
      ${base}/04_dosage_by_chr/<chrom>.dosage.tsv.gz
        same line count as sites; each line is a tab-separated
        vector of n_samples dosage values (0/1/2 or 'NA' / '-1')

    Output JSON shape — see P4.1 patch doc for the schema the
    renderer expects.

    Cap behaviour: if more than `cap` sites fall in [start, end],
    we evenly subsample (every k-th row). Future: use markrank
    column if present.

    Path is ALWAYS resolved relative to ctx.base — no user-controlled
    paths.
    """
    if _ctx is None:
        return _server_not_initialized()
    base = _ctx.base
    sites_path = os.path.join(base, '04_dosage_by_chr',
                              f'{chrom}.sites.tsv.gz')
    dosage_path = os.path.join(base, '04_dosage_by_chr',
                               f'{chrom}.dosage.tsv.gz')
    if not os.path.isfile(sites_path):
        return {"error": f"sites file not found: {sites_path}"}
    if not os.path.isfile(dosage_path):
        return {"error": f"dosage file not found: {dosage_path}"}

    # First pass: collect site indices in [start, end]
    keep_indices = []
    site_records = []
    with gzip.open(sites_path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if not row or row[0].startswith('#'): continue
            try:
                pos = int(row[0])
            except (ValueError, IndexError):
                continue
            if pos < start: continue
            if pos > end:   break  # files are pos-sorted
            ref = row[1] if len(row) > 1 else ''
            alt = row[2] if len(row) > 2 else ''
            keep_indices.append(i)
            site_records.append({"pos": pos, "ref": ref, "alt": alt})

    n_in_region = len(keep_indices)
    if n_in_region == 0:
        return {
            "chrom": chrom, "start_bp": start, "end_bp": end,
            "n_samples": _ctx.n_samples,
            "samples": _ctx.samples,
            "n_sites": 0, "sites": [],
            "dosage_matrix": [],
            "_note": "no sites in region",
        }

    # Cap by even subsampling
    if n_in_region > cap:
        step = n_in_region / cap
        wanted = set(int(i * step) for i in range(cap))
        keep_indices = [keep_indices[i] for i in sorted(wanted)]
        site_records = [site_records[i] for i in sorted(wanted)]

    keep_set = set(keep_indices)

    # Second pass: read dosage rows for those sites
    # Dosage file is row-per-site (n_sites rows × n_samples cols).
    dosage_rows = []
    with gzip.open(dosage_path, 'rt') as f:
        for i, line in enumerate(f):
            if i not in keep_set: continue
            cells = line.rstrip('\n').split('\t')
            row = []
            for c in cells:
                c = c.strip()
                if c in ('NA', '-1', '', '.'):
                    row.append(None)
                else:
                    try: row.append(int(c))
                    except ValueError: row.append(None)
                    if isinstance(row[-1], int) and row[-1] not in (0, 1, 2):
                        row[-1] = None
            dosage_rows.append(row)

    # Transpose: renderer expects n_samples × n_sites
    n_samples = _ctx.n_samples
    if not dosage_rows:
        matrix = []
    else:
        n_sites = len(dosage_rows)
        matrix = [[dosage_rows[si][s] if s < len(dosage_rows[si]) else None
                   for si in range(len(dosage_rows))]
                  for s in range(n_samples)]
        # The above transposes correctly only if every row has
        # n_samples columns. Defensive check:
        if any(len(r) != n_sites for r in matrix):
            matrix = [[None]*n_sites for _ in range(n_samples)]

    return {
        "chrom": chrom,
        "start_bp": start,
        "end_bp": end,
        "n_samples": n_samples,
        "samples": _ctx.samples,
        "n_sites": len(site_records),
        "sites": site_records,
        "dosage_matrix": matrix,
        "_meta": {
            "n_in_region": n_in_region,
            "cap": cap,
            "subsampled": n_in_region > cap,
        },
    }
```

Notes:
- The endpoint never writes anywhere.
- File paths are always under `ctx.base` — no path injection.
- Subsampling is even (every k-th); future improvement: use a
  marker-rank column if your sites file has one.
- `None` for missing dosage in the JSON, which the renderer should
  paint as the "missing" colour.

## Atlas wiring

The atlas dosage heatmap currently checks
`state.data.dosage_chunks`. We have two options:

### Option A — virtual dosage_chunks index pointing at the server

When the atlas detects a server is available, it injects a
`dosage_chunks` index that points at server URLs:

```js
// turn 129 P4.1: when the server is reachable, synthesize a
// dosage_chunks index that delegates fetches to the server
// endpoint. The renderer doesn't need to know the source —
// it just sees an index of chunks with URLs.
function _atlasInstallServerDosageBridge() {
  if (state.data.dosage_chunks && state.data.dosage_chunks.chunks) {
    return; // already have an index from JSON enrichment
  }
  if (!atlasServer || !atlasServer.isAvailable) return;
  if (atlasServer.status !== 'available') return;

  const chrom = state.data.chrom;
  const chrEnd = state.data.chrom_length_bp || 1e9;
  // One synthetic chunk per chromosome — the endpoint handles slicing.
  state.data.dosage_chunks = {
    chunks: [{
      chrom: chrom,
      start_bp: 0,
      end_bp: chrEnd,
      url: atlasServer.url + '/api/dosage/chunk' +
           '?chrom=' + encodeURIComponent(chrom) +
           '&start=__START__&end=__END__&cap=__CAP__',
      _server_synthesized: true,
    }],
    source: 'atlas_server',
  };
}
```

Then in the existing dosage chunk fetcher, replace
`__START__`, `__END__`, `__CAP__` with the actual values:

```js
function _resolveChunkUrl(chunk, start, end, cap) {
  if (!chunk._server_synthesized) return chunk.url;
  return chunk.url
    .replace('__START__', start)
    .replace('__END__', end)
    .replace('__CAP__', cap || 500);
}
```

Hook `_atlasInstallServerDosageBridge` into the existing
"server became available" callback or call it from
`_atlasServerInitBadge` after status flips to `available`.

### Option B — direct fetch from renderer

Less elegant but simpler: when `state.data.dosage_chunks` is null
AND the server is available, the renderer fetches the endpoint
directly:

```js
async function _fetchDosageFromServer(cand) {
  if (!atlasServer || atlasServer.status !== 'available') return null;
  const chrom = state.data.chrom;
  const start = cand.start_bp;
  const end   = cand.end_bp;
  const url = atlasServer.url + '/api/dosage/chunk' +
              '?chrom=' + encodeURIComponent(chrom) +
              '&start=' + start +
              '&end=' + end +
              '&cap=500';
  try {
    const r = await fetch(url);
    if (!r.ok) return null;
    return await r.json();
  } catch (e) {
    console.warn('[dosage bridge]', e);
    return null;
  }
}
```

Pick **Option A** — it preserves the renderer's existing code path
and means the same renderer works for JSON-loaded chunks (R-side
emit) and server-loaded chunks (live).

## Verification

1. Apply patch to server. Restart server.
2. Apply atlas patch.
3. Reload atlas. Server status badge goes green.
4. Activate a candidate on LG28 (e.g. the ~17 Mb shelf).
5. Open the dosage heatmap panel.
6. Network tab shows
   `GET /api/dosage/chunk?chrom=C_gar_LG28&start=...&end=...&cap=500`.
7. Heatmap renders: 226 rows × <500 columns, cells coloured by
   dosage 0/1/2.

## Performance check

```
time curl -s "http://127.0.0.1:8766/api/dosage/chunk?chrom=C_gar_LG28&start=15000000&end=18300000&cap=500" | wc -c
```

Should complete in <2s for typical candidate widths. If it's slow,
the bottleneck is the second pass through the gzipped dosage file —
all rows have to be read and only kept ones decoded. For a 28 Mb
chromosome with ~500K sites this is ~1-3s. If unacceptable, add a
positional index file (next time).

## Test (server-side)

```python
# server_turn1/test_dosage_endpoint.py
from fastapi.testclient import TestClient
import os, sys
os.environ.setdefault('POPSTATS_CONFIG',
    os.path.join(os.path.dirname(__file__), 'popstats_server.local.yaml'))
sys.path.insert(0, os.path.dirname(__file__))
from popstats_server import app

client = TestClient(app)

def test_dosage_lg28_window():
    r = client.get('/api/dosage/chunk',
                   params={'chrom':'C_gar_LG28','start':15000000,'end':18300000,'cap':500})
    assert r.status_code == 200
    body = r.json()
    assert body['chrom'] == 'C_gar_LG28'
    assert body['n_samples'] == 226
    assert body['n_sites'] <= 500
    assert len(body['dosage_matrix']) == 226   # rows = samples
    if body['n_sites'] > 0:
        assert len(body['dosage_matrix'][0]) == body['n_sites']
        # Sanity: cells are 0/1/2/None
        for row in body['dosage_matrix'][:3]:
            for v in row[:10]:
                assert v in (0, 1, 2, None)

def test_dosage_empty_region():
    r = client.get('/api/dosage/chunk',
                   params={'chrom':'C_gar_LG28','start':1,'end':2,'cap':500})
    assert r.status_code == 200
    assert r.json()['n_sites'] == 0
```
