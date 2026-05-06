# Server architecture audit — 2026-05-05 (turn 145)

**One server. One command. All local.**

This doc replaces the ambiguity around "atlas_server vs popstats_server vs
ld_endpoint vs dosage_bridge vs dosage_viewer". Read it once, never wonder again.

---

## TL;DR

You run **one** Python process from the Atlas/ directory:

```bash
./run_atlas.sh
```

That's it. It serves everything the atlas frontend needs:

- `/health` — status probe
- `/file/<path>`, `/compute/<name>` — local file IO + named computes (was `atlas_server.py`)
- `/api/popstats/groupwise`, `/api/popstats/hobs_groupwise` — Engine F + ANGSD wrappers
- `/api/ancestry/groupwise_q` — local Q reads
- `/api/shelf_ld_test`, `/api/ld/split_heatmap` — LD endpoints
- `/api/dosage/chunk`, `/api/dosage/manifest` — dosage chunks for the heatmap
- `/api/cache/keys`, `/api/cache/keys/{key}`, `/api/jobs/{job_id}` — debug/cache

Default bind: `127.0.0.1:8765`. Localhost only. No auth, no HTTPS, no network exposure.

---

## What used to be confusing

Pre-turn-145 there were five Python files that looked like servers, plus a
duplicate folder. After audit:

| File / directory | What it is | Status |
|---|---|---|
| `atlas_server.py` | stdlib `http.server` shim. `/health`, `/file/*`, `/compute/*`. ~260 LOC. | **Folded into popstats_server.py** (turn 145). Kept as a thin compatibility shim that re-exports the same surface for back-compat. |
| `server_turn1/popstats_server.py` | FastAPI app. Real backend. ~1810 LOC pre-turn-145, ~2050 LOC post. | **Active.** This is the unified server. |
| `server_turn1/dosage_bridge.py` | Handler module. Imported by popstats_server. NOT a separate server. | Module of popstats. |
| `server_turn1/ld_endpoint.py` | Handler module. Imported by popstats_server. NOT a separate server. | Module of popstats. |
| `server_turn11c_ld_fast/` | Drop-in spike folder; copies were merged into `server_turn1/` long ago. | **Deleted (turn 145)**. Test relocated to `server_turn1/test_ld_endpoint.py`. |
| `dosage_viewer/02_run_server.py` | Standalone webapp. Has its own static frontend (`dosage_viewer/static/`), parquet store, six sampling modes. The atlas HTML never calls it. | **Standalone product.** Out of scope for the atlas's "one server" model. Run separately if you want to use the standalone viewer. |

The "LANTA cluster" framing in the original `popstats_server.py` top-of-file
comment was an error from earlier handoffs. Everything is local. Drop your
processed data (dosage TSV.GZ, SVs, FASTA, RDS precomps, popstats binaries)
into the project folder, point the server at it via the config file, run
the launcher.

---

## Local layout

```
Atlas/
├── run_atlas.sh                          single launcher (POSIX shells)
├── run_atlas.py                          single launcher (Windows / no bash)
├── Inversion_atlas.html                  the atlas frontend
├── atlas_server.py                       compat shim → calls into server_turn1
├── server_turn1/
│   ├── popstats_server.py                the unified FastAPI server
│   ├── popstats_server.config.yaml       site-local config (paths, bind, engines)
│   ├── ld_endpoint.py                    LD route handler module
│   ├── dosage_bridge.py                  dosage route handler module
│   ├── lazy_windows_json.py              per-chrom window cache
│   ├── requirements.txt                  python deps
│   ├── test_units.py                     popstats unit tests
│   ├── test_ld_endpoint.py               LD endpoint tests (relocated turn 145)
│   └── test_with_curl.sh                 curl smoke tests
├── engine_fast_ld/                       C engine for LD (separate from server)
└── …everything else
```

---

## Bootstrap model

The merged server has **two independent subsystems**, each usable on its own:

### Subsystem 1: `/file` + `/compute` (was atlas_server.py)

- **What it needs:** `--project-root <dir>` (default: current working directory)
- **What you can do without it:** call `/file/*` to read/write any file under
  the project root, call `/compute/<name>` to run a registered Python compute
  function. No external binaries needed.
- **When the subsystem is unavailable:** never. If you launched the server,
  this subsystem is live.

### Subsystem 2: popstats / LD / dosage / ancestry (was popstats_server.py)

- **What it needs:** `--config <yaml>` pointing at `popstats_server.config.yaml`
  with valid paths to `engines:`, `sample_list:`, `beagle_dir:`, etc.
- **What you can do without it:** nothing — these endpoints return HTTP 503 with
  `{"detail": "popstats subsystem not configured. Pass --config to enable."}`
- **When the subsystem is unavailable:** when no `--config` was passed OR the
  config references files that don't exist. The `/health` endpoint shows
  `subsystems.popstats.ready: false` with a `reason` field explaining why.

The launcher script tries to start with both subsystems active:

1. If `server_turn1/popstats_server.config.yaml` exists, popstats subsystem
   is enabled (engines, sample list, etc. loaded at startup)
2. The `/file` + `/compute` subsystem is always enabled with `project-root`
   defaulting to the Atlas/ directory

If you don't have the popstats binaries / BEAGLE / sample list yet, the
server still starts — popstats endpoints just return 503. The atlas
frontend handles this gracefully (the server-status button shows "running,
popstats subsystem unavailable").

---

## Decision matrix

| What I want to do | Run this | Need config? |
|---|---|---|
| Just open the atlas, load JSON files locally, no popstats | `./run_atlas.sh` | No |
| Run popstats live on local data | `./run_atlas.sh` | Yes — copy and edit `server_turn1/popstats_server.config.example.yaml` |
| Use the standalone dosage viewer (separate product) | `python3 dosage_viewer/02_run_server.py --store … --port 8767` | Per its own README |
| Verify what's reachable | Click the **● server** button in the atlas top bar | — |

---

## Status button (top bar, turn 145)

The `● server` chip used to live inside the `data ▾` dropdown. Turn 145
lifts it out to its own top-bar button. States:

- **● server** with green dot — server reachable, all probed endpoints OK
- **● server** with amber dot — reachable but one or more subsystems unavailable
- **● server** with grey dot — not reachable

Click → small popup showing:

1. **Start command** — paste in a terminal at Atlas/ root: `./run_atlas.sh`
2. **Config** — where to copy/edit the example config to enable popstats
3. **Recheck** button — re-probes status without page reload

Nothing more. If you need verbose diagnostics, the popstats `/health`
endpoint returns the full subsystem matrix as JSON.

---

## What's NOT in this audit's scope

- **Folding the standalone dosage_viewer** (`dosage_viewer/02_run_server.py`)
  into the atlas server. It's a different product with its own frontend
  and a parquet store that the atlas doesn't use. Folding it in adds
  ~1500 LOC and one heavy dependency (pyarrow) for a workflow the atlas
  itself doesn't need. Run it separately if you want to use the
  standalone region-explorer viewer.

- **Cluster deployment** (`DEPLOY_TO_LANTA.md`). That document is a
  historical artifact from when there was a LANTA-cluster path. The
  atlas is now local-only. If you do want to run the server on a remote
  machine and tunnel back to localhost for the atlas, the same launcher
  + an SSH tunnel works (`ssh -L 8765:127.0.0.1:8765 lanta`), but
  there's no atlas-side config for remote URLs — just keep the atlas
  pointing at `localhost:8765`.

---

## Migration notes

If you had any local scripts that hit `atlas_server.py` directly via
`python3 atlas_server.py …`, they keep working — `atlas_server.py` is
now a thin shim that delegates to the merged server. The same CLI
arguments (`<project_root>`, `--port`, `--host`) are accepted.

If you were running both `atlas_server.py` AND `popstats_server.py` on
port 8765, that conflict no longer exists — the merged server is one
process and there's only one launcher.
