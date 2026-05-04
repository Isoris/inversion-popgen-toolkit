# Chat A — live popstats system · Bundle index + LG28 demo

This is the close-out bundle for chat A. The system lets you compose
arbitrary sample groupings in the atlas, push them to a thin server
wrapping `region_popstats` / `hobs_windower` / `angsd_fixed_HWE`, and
get window-aligned per-group population-genetic statistics back into
the popstats stack — instantly rerendering, no cluster reruns.

## What's in this bundle

| file | what |
|---|---|
| `LG28_DEMO.md` | The click sequence with manuscript context |
| `lg28_demo_dashboard.html` | Standalone interactive checklist (open in any browser) |
| `README.md` | This file |

## The full chat-A delivery

Eight modules across two halves:

### Server side (1 module, deploys to LANTA)

| bundle | purpose |
|---|---|
| `popstats_server_turn1.tar.gz` | Python FastAPI server. 6 endpoints, DiskCache LRU, engine-binary SHA-256 versioning, parallel ANGSD. 24 tests (14 unit + 10 curl). Runs on the LANTA login node. Bind 127.0.0.1:8765. |

### Atlas side (7 modules, drop into `Inversion_atlas.html`)

| bundle | purpose | tests |
|---|---|---|
| `atlas_groups_turn2_v2.tar.gz` | Group composition engine — dimensions, expressions, slots, sets, lasso history, snapshots. UMD `window.popgen`. | 70 |
| `atlas_request_layer_turn3.tar.gz` | Server adapter — IndexedDB cache, debounce, dedup, abort, content-addressed key matching server byte-for-byte. UMD `window.popgenLive`. | 31 |
| `atlas_renderers_turn4.tar.gz` | Canvas renderers + adapters — multi-line, bars, heatmap, categorical strip, mode chips. UMD `window.popgenRenderers`. | 63 |
| `atlas_floating_dock_turn5a.tar.gz` | Floating dock shell — drag, resize, launcher, three-tab bar, keybinds. UMD `window.popgenDock`. | 50 |
| `atlas_dock_tabs_turn5b.tar.gz` | Three tab UIs — Groups (dimensions/draft/slots/Compute), Lasso, Snapshots. UMD `window.popgenDockTabs`. | 41 |
| `atlas_page6_wiring_turn6.tar.gz` | Page 6 wiring — getData() resolvers (static/live/split), tooltip, slot-combine UI. UMD `window.popgenPage6`. | 40 |
| `atlas_track_gallery_turn6_5.tar.gz` | Dynamic track discovery + view presets. UMD `window.popgenGallery`. | 78 |
| `atlas_turn7.tar.gz` | Q09b shelf-LD test in-browser + page-7 stratified Δ12. UMD `window.popgenTurn7`. | 43 |

**Atlas-side total: 416 tests. Server-side total: 24 tests. Grand total: 440 tests, all green.**

## Load order

In `Inversion_atlas.html`, place these `<script>` tags after the existing
atlas script. Order matters — turn N depends on turns 1..N-1:

```html
<script src="atlas_group_engine.js"></script>      <!-- turn 2 -->
<script src="atlas_request_layer.js"></script>     <!-- turn 3 -->
<script src="atlas_renderers_turn4.js"></script>   <!-- turn 4 -->
<script src="atlas_floating_dock.js"></script>     <!-- turn 5a -->
<script src="atlas_dock_tabs.js"></script>         <!-- turn 5b -->
<script src="atlas_page6_wiring.js"></script>      <!-- turn 6 -->
<script src="atlas_track_gallery.js"></script>     <!-- turn 6.5 -->
<script src="atlas_turn7.js"></script>             <!-- turn 7 -->
```

Each bundle includes its own `SPLICE_POINTS.md` with the specific
patches into the existing atlas code. The total splice across all
seven turns is roughly 60–80 lines of edits in `Inversion_atlas.html`
plus the 8 `<script>` tags.

## Quick start

1. **Server side (LANTA login node)**:
   ```bash
   cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/popstats_server/
   tar xzf popstats_server_turn1.tar.gz
   cd popstats_server/
   pip install -r requirements.txt
   python -m popstats_server  # binds 127.0.0.1:8765
   ```

2. **Tunnel from your laptop**:
   ```bash
   ssh -L 8765:127.0.0.1:8765 lanta
   ```
   Verify: `curl -s http://127.0.0.1:8765/api/health | jq .ok` → `true`.

3. **Atlas side**: extract all 7 atlas bundles into the same directory
   as `Inversion_atlas.html`. Apply the splices from each bundle's
   `SPLICE_POINTS.md`. Open the atlas in a browser.

4. **Run the LG28 demo**: open `lg28_demo_dashboard.html` in a second
   tab and follow the click sequence. Every step is reproducible and
   the dashboard tracks your progress.

## Configuration

The server reads its paths from `popstats_server_config.yaml`. Adjust
to your LANTA layout:

```yaml
base_dir:           /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
beagle_dir:         {base_dir}/<existing BEAGLE path>
local_q_dir:        {base_dir}/<existing per-K Q matrices>
region_popstats_bin: {base_dir}/inversion-popgen-toolkit/.../region_popstats
ref_fai:            {base_dir}/<reference index>.fa.fai
sample_master_list: {base_dir}/<sample master>.txt
cache_dir:          {base_dir}/popstats_cache/
log_dir:            {base_dir}/popstats_logs/
```

These mirror the values from `inversion-popgen-toolkit/unified_ancestry/00_ancestry_config.sh`.

## Verification

After splicing, the LG28 demo dashboard's pre-requisite checklist
walks you through verifying every integration point. If all six
prereqs check, the click sequence will work.

If anything fails:
1. Open browser dev tools → console. The atlas modules log a startup
   line each. Missing module = missing `<script>` tag or wrong load
   order.
2. The dock keybind G is suppressed when an input is focused — make
   sure you're not in a text field.
3. The Compute button returns `err: groups` if no slots are populated.
   Add at least one slot first.
4. If the popstats stack stays empty after Compute, check `state.popstatsLive.lastResponse`
   in the console — should be a real object with `windows` array.

## Manuscript anchor

This system underwrites Result 3 / Figure 3 of MS_Inversions. The
Methods text drop-ins are in both `LG28_DEMO.md` and the demo
dashboard. The reproducibility footprint (440 tests across 8 modules)
is suitable for inclusion in the manuscript supplement.

## What's next (not chat A)

Out of scope but tracked for future work:
- GHSL stratified live-mode (parallel server endpoint over phased
  haplotype data)
- L2 boundary scan integration (turn 9 idea — spreadsheet-style
  candidate overview with per-row Compute on click)
- `compareGroups(A, B)` helper exposing |A∩B|, jaccard, overlap
  (turn 10 idea, used by step 10 of the LG28 demo)
- Auto-snapshot on every Compute (recovery from accidental dock
  state loss)

These are nice-to-haves. The system as shipped produces every panel
of Figure 3 with reviewer-grade reproducibility.
