# SPEC — Multi-chromosome JSON load orchestrator

**Status**: drafted turn 130 final session. Not yet implemented.
Estimated ~0.5–1 turn.

**Trigger** (Quentin):
> *"When we load all JSONs at once."*
>
> *"Soon its all ready then upload all chromosomes and it will be all
> automatic."*

The atlas today loads one chromosome JSON at a time. To support the
full-automatic vision (all 28 chroms → all candidates auto-promoted →
genome-wide ideogram), we need an orchestrator that handles loading
many JSONs without OOM-ing the browser.

---

## 1. The constraints

- 28 chromosome JSONs, ~50–300 MB each (varies with what layers are
  bundled). Worst case: 28 × 300 MB = 8.4 GB. Browser tabs typically
  cap around 4 GB heap on desktop, less on mobile.
- IndexedDB cache exists today (turn 95-ish from `07759823` chat) for
  per-chrom persistence. Per-origin quota varies (some browsers cap at
  200 MB, others at 60% of disk).
- The expensive computes (lineage, sweep, sliding-window) are
  per-chromosome and don't need cross-chrom data simultaneously.

## 2. Architecture

### 2.1 Lazy load with active-chrom focus

Default mode: **only the active chromosome's full JSON is in memory**.
Other chroms are reduced to a lightweight summary blob:

```
chromSummary[chrom] = {
  n_l2_envelopes,
  n_candidates,
  layersPresent,
  hasLineageCache,
  hasInheritanceCache,
  // small enough to keep in memory for all 28 chroms
}
```

When the user switches active chrom, the previous chrom's full data
goes back to IndexedDB; the new chrom's full data is loaded from
IndexedDB into memory.

### 2.2 Background pre-compute scheduling

When the user loads the multi-chrom bundle, the orchestrator schedules
the per-chrom expensive computes (lineage, sweep) via
`requestIdleCallback` chains:

```
queue = chrom_list.copy()
function nextChromCompute() {
  if (queue.empty) return;
  const chrom = queue.shift();
  loadChromIntoMemory(chrom);
  runLineageCompute({ chrom });
  runL2SweepInheritance({ chrom });
  cacheToIDB(chrom);
  unloadChromFromMemory(chrom);
  requestIdleCallback(nextChromCompute);
}
```

So the user can browse one chrom while the others are computing in the
background. Progress strip in the sidebar: "12 / 28 chromosomes
processed."

### 2.3 Bulk upload UI

New page or sidebar widget: "Bulk load 28 JSONs." Drag-drop a folder
or zip file. Atlas detects which files are chromosome JSONs (by
filename pattern + JSON `chrom` field), persists each to IDB, builds
chromSummary, schedules background compute.

## 3. State additions

```
state.chromSummary: { [chrom]: ChromSummary }
state.activeChrom: string
state.bulkLoadProgress: { total, processed, current }
state.computeQueue: chrom[]
```

## 4. Implementation slices

### Slice 1 — chromSummary infrastructure (~0.5 turn)
- Build summary blob on chrom load
- Persist to IDB alongside the full JSON
- Lazy-load full JSON only when chrom becomes active

### Slice 2 — bulk-load UI (~0.5 turn)
- Drag-drop zone for folder/zip
- Filename pattern detection
- Progress feedback during ingestion

### Slice 3 — background compute scheduler (~0.5 turn)
- requestIdleCallback queue
- Per-chrom progress
- Skip chroms already computed (cache hit)

## 5. Open questions

1. **What goes in chromSummary?** Minimum to support genome-wide
   ideogram + cross-chrom lineage matching. Probably:
   `n_candidates`, `n_auto_promoted`, `n_lineages`, `chromosome_size_bp`,
   per-candidate `(start_bp, end_bp, tier, source)`.
2. **What if the user wants two chroms in memory simultaneously?**
   E.g., comparing LG12 and LG28 side-by-side. Edge case; defer.
3. **What about Save-Session?** The session export today is per-chrom
   data + state. Multi-chrom session adds chromSummary + IDB pointers.
   May need a session format bump.

## 6. Tests

- Bulk-load 3 synthetic JSONs → all 3 land in IDB, chromSummary
  populated, only the active one is in `state.data`.
- Switch active chrom → previous chrom's data unloaded, new chrom
  loaded.
- Background compute: schedule 3 chroms, verify each gets lineage +
  sweep results in IDB.
- Memory budget: synthetic 28 × 200MB JSONs → browser doesn't OOM
  (only one in memory at a time).

## 7. Dependencies

- Existing IDB cache (already shipped).
- Lineage compute (shipped turn 130).
- L2-sweep auto-promote (`SPEC_l2_sweep_inheritance.md`, queued).
- Cross-chrom lineage aggregator (`SPEC_cross_chromosome_lineages.md`,
  queued).
