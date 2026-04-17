# SPEC — deferred scan and classification features (post-chat-16)

This doc specifies features that are **designed but not implemented**,
for chat 17+ once LANTA data is in hand. Each has a clear input/output
contract, a name for where it should live, and acceptance criteria.

---

## 1. `reg$compute$scan_with_ancestry_resolution(chrom, start, end, ...)`

**Status:** spec only — requires ancestry cache (`LOCAL_Q_DIR/K08/`) to
be populated. Blocked until chat-17 Step 2 (canonical-K precompute)
succeeds on LANTA.

**Purpose.** Like `scan_pairwise`, but at each window the "group" is
dynamically defined by which ancestry cluster each sample belongs to
at that window (not at cohort scale). Reveals locus-specific ancestry
switches that happen inside an inversion — particularly useful when an
inversion region shows introgressed ancestry only in one karyotype.

**Input contract.**
```r
scan_with_ancestry_resolution(
  chrom,
  start_bp, end_bp,
  window_size = 10000, step = NULL,
  K = NULL,                    # ancestry K level (default canonical)
  karyotype_filter = "HOM_INV",# restrict to this karyotype
  min_group_n = 5,
  stat = "fst",                # fst between ancestry clusters
  persist = TRUE
)
```

**Behavior per window.**
1. Call `reg$compute$ancestry_q_vector(chrom, pos - ws/2, pos + ws/2, K)`
   to get per-sample ancestry at this window.
2. Intersect with the karyotype group (e.g. `inv_<cid>_HOM_INV` resolved
   via `resolve_smallest_at`).
3. Bucket samples by their assigned_pop at this window → dynamic groups.
4. Compute pairwise FST between the top-two ancestry clusters within
   the karyotype group.
5. Persist as `interval_summary` with `stat="ancestry_q_mean"` carrying
   the per-window ancestry composition.

**Acceptance criteria.**
- Handles windows where all samples belong to one ancestry cluster
  (output NA, not error)
- Handles windows where karyotype × ancestry intersection has <
  min_group_n samples (falls back to parent karyotype group and logs)
- Output shape: data.table with `window_mid`, `candidate_id`,
  `cluster_1`, `cluster_2`, `n_1`, `n_2`, `fst`, `used_fallback`

---

## 2. `scan_all_candidates_driver.R` — the full automatic sweeper

**Status:** spec only — ~200 lines of R, depends on segment registration
being routine (chat-17 validation of `STEP_C01i_b_multi_recomb.R`
segment-writer).

**Purpose.** One script that iterates every registered candidate,
calls `scan_candidate_full()` for each, writes the full population
genetics panel (FST between all karyotype pairs + Δ12/entropy/ENA
within each karyotype) to the results_registry. This is the "click
and leave it running for an hour, come back to a fully populated
database" step.

**Input contract.**
```bash
Rscript scan_all_candidates_driver.R \
  --tier_max 2              `# only scan tier 1 + 2 candidates` \
  --window_size 10000 \
  --flank_kb 100 \
  --min_group_n 3 \
  --karyotypes HOM_REF,HET,HOM_INV \
  --stats fst,dxy \
  --region_stats delta12,entropy,ena \
  --resume                  `# skip cids already cached` \
  --parallel 4              `# mcmapply workers`
```

**Behavior.**
1. `cids <- reg$evidence$list_candidates()` (or tier-filtered via
   `registry_list_candidates_by_tier`)
2. For each cid:
   - Check `reg$results$ask_what_for_candidate(cid)` — if already
     covered for the requested stats × karyotype pairs, skip (`--resume`)
   - Call `reg$compute$scan_candidate_full(cid, ...)`
   - Log `reg$results$session_summary()` per candidate
3. After loop: `reg$results$integrity_check()` — fail exit if any
   FK violations or stale group_versions

**Acceptance criteria.**
- Fully automatic — no per-candidate configuration
- Idempotent — second run with `--resume` writes nothing
- Logs per-candidate wall time + n_rows written
- Exits 0 only if integrity_check passes

**Deferred decisions** (should be made after LG12 is scanned once):
- Default window size (10 kb? 20 kb? 50 kb?)
- Default flank_kb (depends on local recombination rate at 9× coverage)
- min_group_n threshold for scan fallback
- Whether to scan `RECOMBINANT*` carriers separately or just fold into
  parent HOM_INV

---

## 3. `stale_segment_gc.R` — recombinant-map re-run cleanup

**Status:** spec only — needed once segmentation is routine. Not a
chat-17 blocker (a stale segment doesn't break scans, just causes
incorrect group membership for affected samples).

**Purpose.** When `STEP_C01i_b_multi_recomb.R` is re-run on a candidate
(e.g. after a threshold change or a sample-level QC fix), the previous
segment candidates and their groups become stale. This script deletes
invalidated segments and their downstream results_registry rows.

**Input contract.**
```bash
Rscript stale_segment_gc.R --cid LG12_17 [--dry-run]
```

**Behavior.**
1. Load the current recombinant_map block for cid
2. Enumerate existing segments `reg$intervals$get_children(cid)` with
   scale ∈ {`seg`, `seg_dco`}
3. Compare against the current block's breakpoint set
4. For any segment no longer in the breakpoint set:
   - Delete its candidate from interval_registry
   - Delete all `inv_<seg_cid>_*` groups from sample_registry
   - Delete all manifest rows mentioning the seg_cid
   - Move any files to `results_registry/_archive/<timestamp>/`
5. Re-register new segments via the existing multi_recomb writer

**Acceptance criteria.**
- `--dry-run` produces a "would delete N candidates, M groups, P manifest rows" report
- Archive files are preserved so historical runs can still be audited
- Post-run `reg$results$integrity_check()` passes

---

## 4. `STEP_C01i_e_subcluster.R` — recursive PCA sub-cluster detector

**Status:** spec only — feature for the 10% biology case (like your
cand46 figure). Needed once HPC run shows candidates with visible
sub-structure inside karyotypes.

**Purpose.** For each candidate where PCA inside a karyotype shows
dbscan-separable sub-clusters (as in FIG_C14), register those
sub-clusters as groups with the naming convention
`inv_<cid>_<KARYO>__sub<N>` / `__noise`.

**Input contract.**
```bash
Rscript STEP_C01i_e_subcluster.R \
  --candidates <c01d>/candidates.tsv \
  --decomp_dir <decompose_out>/ \
  --karyotypes HOM_REF,HET,HOM_INV \
  --dbscan_eps 5.0          `# in PC units` \
  --dbscan_minpts 5 \
  --min_cluster_size 5 \
  --outdir <out>
```

**Behavior per (cid, karyotype).**
1. Load per_window_class.rds to identify samples in `inv_<cid>_<karyo>`
2. Run dbscan on the candidate-wide PC1/PC2 of those samples
3. For each cluster with ≥ min_cluster_size members:
   - `reg$samples$add_group(sprintf("inv_%s_%s__sub%d", cid, karyo, N), ...)`
     with `dimension="karyotype_subcluster"`
4. For noise points (cluster label -1):
   - `reg$samples$add_group(sprintf("inv_%s_%s__noise", cid, karyo), ...)`
5. Write a `subcluster_report` Tier-2 block documenting the n_clusters
   found, silhouette scores, and sample-to-cluster assignments

**Acceptance criteria.**
- Noise-only karyotypes (no dbscan clusters) register no __sub groups
  (just __noise if nonempty)
- Re-running is idempotent: existing groups are updated not duplicated
- The `subcluster_report` block's `keys_extracted` adds
  `q2_n_subclusters_<KARYO>` to the BK manifest

---

## 5. `classification_loop.R` — little-by-little candidate review

**Status:** workflow spec — not code; describes your "stop on each
inversion, understand it, upgrade if needed" iteration loop.

**Per-candidate review checklist.**
1. `have <- reg$results$ask_what_for_candidate(cid)` — should have:
   - 3 pairwise FST rows (HOM_REF vs HET, HOM_REF vs HOM_INV, HET vs HOM_INV)
   - 3 or 9 interval_summary rows (Δ12/entropy/ENA × n karyotypes × n segments)
   - 1 candidate_q + 1 candidate_f
2. If all present → run `reg$results$integrity_check()` → pass → move on
3. If something missing or flagged:
   - Open the relevant scan tsv.gz and inspect windows
   - Check `q6_group_validation` — if UNCERTAIN, the sample grouping
     is what needs work, not the scan
   - Check `n_1 / n_2` columns — if too many `used_fallback=TRUE`, the
     segmentation is degenerate; consider relaxing min_group_n or
     widening flank
   - If scan shows zero signal where expected, the candidate is likely
     spurious — lower its tier or demote to SUSPECT
   - If scan shows signal only in ONE segment, the rest of the candidate
     may be a false extension — consider splitting into two candidates

**When a specific inversion can't be explained by the current registry
semantics**, update DATABASE_DESIGN.md with a new case study in §
"Known complex cases" and extend the schema or conventions
incrementally. Never special-case in code without documenting.

**Classification output format.** One TSV per candidate in
`classification_results/<cid>/status.tsv` with:
- `cid`, `chrom`, `start_bp`, `end_bp`, `size_kb`
- `tier` (final), `validation_level` (final)
- `n_hom_ref`, `n_het`, `n_hom_inv`, `n_recombinant*`
- `mean_fst_HOM_REF_vs_HOM_INV` (from scan, within parent interval)
- `peak_delta12_HOM_INV` (from region scan)
- `is_introgressed` (bool, from ancestry scan)
- `is_under_selection` (bool, from Tajima's D or similar — future)
- `n_subclusters_<KARYO>` (if sub-cluster detector has run)
- `n_segments` (from recombinant_map)
- `notes` — free text for manual review

---

## 6. `sweep_figures_gallery.R` + `build_inversion_atlas.R` — reviewable outputs

**Status:** spec only — chat 18 or later, once the first per-candidate
figures exist. Not a code blocker for HPC; purely a workflow
convenience.

**Purpose.** Solve the "200 folders on Windows" browsing problem.
Canonical per-candidate storage is great for organization but terrible
for review sessions. Generate two derived outputs from the canonical
storage:

1. **Flat gallery** — `figures_gallery/<cid>__summary.png` per
   candidate. Browseable with left/right arrows in Windows Photos or
   any file manager.
2. **PDF atlas** — `inversion_atlas_<date>.pdf`, one page per candidate.
   Searchable, shareable, bookmarkable.

### `sweep_figures_gallery.R` contract

```bash
Rscript sweep_figures_gallery.R \
  --root        <BASE>/registries/data/evidence_registry/per_candidate \
  --gallery     <BASE>/figures_gallery \
  --kinds       summary,dbscan,heatmap,scan_fst \
  --sort_by     tier,chrom,start_bp \
  --symlink     FALSE       # TRUE on Linux-only; FALSE copies (Windows-safe)
```

**Behavior.**
1. Walk `per_candidate/<cid>/figures/<kind>.png` for each registered
   candidate and each requested kind
2. Copy (or symlink) to `gallery/<kind>/<NNN>_<cid>.png` where NNN is
   the sort index (so file manager sort-by-name = sort by whatever
   criterion was chosen)
3. Write `gallery/INDEX.tsv` — the sort key mapping, so a reviewer can
   see "which candidate is page 47?"
4. Idempotent: skip if source is older than target

**Acceptance.**
- Works over SMB/SSHFS to Windows (copy mode, not symlink)
- Sort order survives alphabetic file-manager sort (zero-padded indices)
- Missing figures don't crash the sweep (log warning, continue)

### `build_inversion_atlas.R` contract

```bash
Rscript build_inversion_atlas.R \
  --root        <BASE>/registries/data/evidence_registry/per_candidate \
  --output      <BASE>/figures_atlas/inversion_atlas_2026-MM-DD.pdf \
  --layout      composite   # composite|single_panel
  --sort_by     tier,chrom,start_bp \
  --filter      "tier<=2"   # optional R expression filter
```

**Behavior per page.**
1. Header: cid, chrom, start-end bp, span, tier, validation level
2. Multi-panel composite (via figrid) — PCA, heterozygosity, DBSCAN,
   scan panels, marker heatmap
3. Caption with key BK numbers (n_HOM_REF, n_HOM_INV, n_recombinants,
   mean FST, etc.)
4. Bookmarks by chromosome for jump navigation

**Sort options.**
- `tier,chrom,start_bp` — default, good for review
- `tier,confidence_desc` — most-likely-real first, for triage
- `complexity_desc` — candidates with sub-clusters or segments first,
  for "interesting cases" review

**Deferred decisions.**
- Exact composite layout (which panels on which grid positions)
- Whether to include scan tracks inline or as appendix pages
- Page size (A4 letter vs US letter; probably A4 given Thai affiliation)
- Whether to include a summary table as page 1

### `reg$query$candidate_figure_paths(cid, kinds)`

Helper needed by both scripts above. Returns a named list of canonical
paths for a candidate's figures, e.g.:

```r
reg$query$candidate_figure_paths("LG12_17",
  kinds = c("summary", "dbscan", "heatmap"))
#> $summary  "evidence_registry/per_candidate/LG12_17/figures/summary.png"
#> $dbscan   "evidence_registry/per_candidate/LG12_17/figures/dbscan.png"
#> $heatmap  "evidence_registry/per_candidate/LG12_17/figures/heatmap.png"
```

Returns NA for kinds not found. Makes gallery / atlas scripts portable
across directory structure changes.

---

## Order of work chat 17 onwards

1. **Chat 17** — LANTA canonical-K precompute + BK extraction validation.
   No new code; just run.
2. **Chat 17 (end)** — pick one candidate from the sweep, run
   `scan_candidate_full` on it manually, verify the manifest shape and
   segment resolution. If segments don't exist because multi_recomb
   hasn't emitted any recombinants for the test candidate, confirm
   parent-scale scan works.
3. **Chat 18** — iterate the `classification_loop` on LG12 candidates.
   One at a time. Note which ones work cleanly (update TRACK_GOOD.md)
   and which don't (update TRACK_ODD.md with the specific failure).
4. **Chat 19** — build `scan_all_candidates_driver.R` with thresholds
   informed by chat-18 iteration. This is the "leave it running, come
   back to a full database" script.
5. **Chat 20+** — `STEP_C01i_e_subcluster.R` if chat-18 surfaced
   sub-structure that needs it. `scan_with_ancestry_resolution` if
   ancestry-switching inversions are seen. `stale_segment_gc.R` if any
   multi_recomb re-runs happen.

---

## What NOT to build pre-emptively

- A generic pan-stat scanner that handles FST + dxy + Tajima's D + piN/piS
  + XP-EHH in one framework. `scan_pairwise` supports adding new `stat`
  values one at a time; bulk-add creates lock-in on an API that may not
  match the stat library used (angsd, moments, dadi, etc.)
- A GUI for candidate review. Classification is more productive as a
  checklist in a markdown file than as clicky software.
- Pre-emptive sub-cluster schemas for every possible biology case.
  Extend when a specific candidate demands it.
