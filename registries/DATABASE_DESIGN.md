# Registry database — design reference

One coherent database. Four tables. Explicit foreign keys. One query plane.

This document is the single source of truth for the registry system. If the
code disagrees with this document, the code is wrong.

**Status:** chat 16 rewrite of the chat-15 `stats_cache` sketch. Previous
state: four sub-APIs (`reg$samples` / `reg$intervals` / `reg$evidence` /
`reg$stats`) with no explicit join relationships and an opaque 6-char sha1
tag used to identify sample sets in filenames. Current state: four tables
with FK discipline, a `results_registry` replacing `stats_cache`, named
sample groups everywhere, and one `reg$ask()` query plane.

---

## The four tables

| Table | Physical location | Primary key | What it stores |
|---|---|---|---|
| **sample_registry**   | `data/sample_registry/`   | `group_id` | WHO: sample master, groups |
| **interval_registry** | `data/interval_registry/` | `candidate_id`, `(chrom, start_bp, end_bp)` | WHERE: windows, candidate intervals, cov paths |
| **evidence_registry** | `data/evidence_registry/` | `(candidate_id, key)` | WHAT-SCALAR: per-candidate Tier-2 JSON blocks + flat keys.tsv |
| **results_registry**  | `data/results_registry/manifest.tsv` | `row_id` (uuid) | WHAT-NUMERICAL: persisted FST tracks, Q/F matrices, ancestry summaries |

No table stores data that logically belongs in another table. No row is
written to two tables. Results go in `results_registry`; per-candidate
scalar facts go in `evidence_registry`; group membership goes in
`sample_registry`; coordinates go in `interval_registry`.

---

## Foreign keys — how the tables join

All FKs reference by string identifier, never by hash or content fingerprint.
The referenced identifier must already exist in the target table at write
time; `put_*` methods reject writes with unresolved FKs.

```
                        results_registry.manifest.tsv
                        ─────────────────────────────
                            row_id (PK)
       ┌──── where.candidate_id ──► interval_registry.candidate_intervals.candidate_id
       │                                        │
       │                                        └──► (chrom, start_bp, end_bp)
       │
       ├──── where.chrom (no FK; just a string)
       │
       ├──── who_1.group_id  ──────► sample_registry.sample_groups.group_id
       ├──── who_2.group_id  ──────► sample_registry.sample_groups.group_id
       │
       └──── what.stat       (enum; no FK)
             what.K          (nullable integer; no FK)

                        evidence_registry.per_candidate/<cid>/
                        ──────────────────────────────────────
                            candidate_id ──► interval_registry.candidate_id
                            structured/<block>.json
                              ├─ validated against
                              │  schemas/structured_block_schemas/<block>.schema.json
                              └─ keys.tsv  (one row per flat key)

                        sample_registry.sample_groups
                        ─────────────────────────────
                            group_id (PK)
                            members_file ──► groups/<group_id>.txt
                            created (timestamp ≈ version)
```

---

## Group versioning

Every `sample_groups.tsv` row has a `created` timestamp. We treat this as the
group's version. When a `results_registry` row is written, the current
`created` value of the referenced group(s) is copied into the manifest as
`who_1.group_version` / `who_2.group_version`.

On read, the results_registry compares the stored `group_version` to the
current `created` value. A mismatch means the group's membership has
changed since the result was computed — the result is flagged stale. This
replaces chat-15's scheme of hashing sample lists to detect collisions: a
registered group's version is already a single source of truth, so no
second source (hash) is needed.

**What group versioning catches.** Someone adds a sample to `all_226` (it
becomes `all_227`, or `all_226` gets overwritten). Any result computed
against the old definition is automatically flagged stale by
`integrity_check()`. Without this, the old results would silently be mixed
with new computations and reviewers would have no way to tell.

**What it does not catch.** Sample master updates (e.g. a sample's metadata
row changes but it stays in the same group). That's OK — metadata changes
don't invalidate numerical results; group membership changes do.

---

## The manifest row — results_registry's schema in plain English

Every numerical artifact in the pipeline is a single row in
`data/results_registry/manifest.tsv`. The row declares four things:

1. **Where** — which region of the genome (chrom + coords or candidate_id)
2. **Who** — which sample group(s) this computation was over
3. **What** — which statistic / metric
4. **How we got here** — the provenance block (script, engine, config hash)

plus the physical location of the data file on disk and its shape.

The full schema is in `schemas/registry_schemas/result_row.schema.json`.

### Four `kind` values

| kind | Populated fields | Example |
|---|---|---|
| `pairwise` | chrom, who_1, who_2, stat | FST between `ancestry_K8_Q3` and `ancestry_K8_Q5` on LG12 |
| `candidate_q` | candidate_id, who_1, K | Per-sample Q matrix for candidate LG12_17 at K=8 |
| `candidate_f` | candidate_id, K | Global F matrix at K=8 (no sample-group dimension) |
| `interval_summary` | chrom, start_bp, end_bp, who_1, stat, K | Ancestry window summary for LG12 1Mb–1.1Mb on `all_226` |

---

## File naming convention

No hashes in filenames. All identifiers are either registered group names or
candidate IDs.

```
results_registry/
├── manifest.tsv
├── pairwise/
│   └── <stat>/                     e.g. fst, dxy, theta_pi
│       └── <chrom>/
│           └── <group1>__vs__<group2>.tsv.gz
├── candidate/
│   └── <candidate_id>/
│       ├── Q_K<NN>.<group>.tsv.gz
│       ├── F_K<NN>.tsv.gz
│       └── meta.tsv
└── interval/
    └── <chrom>/
        └── <start_bp>_<end_bp>.<group>.<stat>.K<NN>.tsv.gz
```

`<group1>__vs__<group2>` is canonically sorted so `A vs B` and `B vs A`
land in the same file. No duplicate storage possible.

---

## The query plane — `reg$ask()`

One method. Four filter dimensions. Any subset can be NULL.

```r
reg$ask(
  where   = list(chrom = NULL, start_bp = NULL, end_bp = NULL, candidate_id = NULL),
  who     = NULL,                  # a group_id, or c(g1, g2) for pairwise pairs
  what    = NULL,                  # stat enum: "fst" | "dxy" | "Q" | "F" | ...
  kind    = NULL,                  # "pairwise" | "candidate_q" | "candidate_f" | "interval_summary"
  K       = NULL,                  # filter by K level
  overlap = "any"                  # "any" | "contains" | "contained_by"
)
# returns data.table with all manifest rows matching the filters
```

### Convenience shortcuts

| Shortcut | Calls `ask()` with |
|---|---|
| `reg$ask_what_at(chrom, start, end)` | `where=list(chrom, start_bp, end_bp)` |
| `reg$ask_what_for_candidate(cid)` | `where=list(candidate_id=cid)` |
| `reg$ask_what_for_group(gid)` | `who=gid` |
| `reg$ask_provenance(row_id)` | returns the full provenance block for one row |

This is the "database feel" — same shape as a SQL SELECT with WHERE clauses.

---

## Integrity check — the one command that proves the database is consistent

```r
reg$integrity_check()
```

Checks, in order:

1. Every `who_1.group_id` / `who_2.group_id` in the manifest exists in sample_registry
2. Every `where.candidate_id` in the manifest exists in interval_registry
3. Every `group_version` in the manifest matches the current `created` value for that group (flags stale, doesn't fail)
4. Every `file` in the manifest exists on disk
5. Every `sha256` in the manifest matches the file's actual sha256 (optional; off by default because it's expensive)
6. Every file under `results_registry/pairwise|candidate|interval/` is listed in the manifest (no orphans)
7. Every `candidate_id` mentioned in evidence_registry exists in interval_registry

Returns a data.table with one row per check, columns `check / pass / n_violations / details`. A clean run prints `all 7 checks passed`. A dirty run lists the specific rows / files at fault.

This is the command Quentin runs before submitting the tarball to a journal. And it's the command a reviewer runs to verify reproducibility claims.

---

## Three-language bindings

| Layer | R (`api/R/registry_loader.R`) | Python (`api/python/registry_loader.py`) | Bash (`api/bash/registry_loader.sh`) |
|---|---|---|---|
| Read manifest | full | full | path only (`registry_results_path`) |
| Query (`ask` + shortcuts) | full | read-only (`ask()` only) | — |
| Write (atomic `put_*`) | full | full (for scripts that write results in Python, e.g. `STEP_C01i_c_nested_composition.py`) | — |
| Schema validation | full | full | — |
| Integrity check | full | — | — |

Bash stays atomic — same philosophy as the existing three registries. R is the only full-API binding.

---

## Migration notes (chat-15 → chat-16)

- `reg$stats` → `reg$results` (R). `reg$stats` kept as a deprecated alias for one chat cycle; emits `.Deprecated` on call.
- `data/stats_cache/` → `data/results_registry/`. Loader auto-moves on init if the old dir exists and the new dir does not.
- Filename hashes (`N226_fe01a5` etc.) → group names. Back-compat reader tolerates the old pattern and logs a deprecation message.
- `get_sample_set_tag()` → `get_active_sample_group()`. The active group is read from `$SAMPLE_GROUP` env (default `all_226`) and must be a registered group or a clear error is thrown.
- Launcher `SAMPLE_SET_TAG` → `SAMPLE_GROUP`. No hashing in bash.

All changes are covered by the test suite in `registries/tests/test_results_registry.R`.

---

## Sample group naming convention (including sub-clusters)

The database is agnostic about what a sample group *means* — it only knows
that `group_id` is a registered string in `sample_registry`. The meaning
lives in the naming convention. For the 10% of candidates that have real
sub-structure inside karyotypes (multiple sub-clusters visible in the
candidate-PCA, driven by ancestry / nested haplotypes / GHSL sub-bands),
this convention lets both pipeline writers and classification readers
agree on group naming without extra schema fields.

### Top-level karyotype groups (always present)

```
inv_<cid>_<KARYOTYPE>       — KARYOTYPE in {HOM_REF, HET, HOM_INV,
                                 RECOMBINANT, RECOMBINANT_GC, RECOMBINANT_DCO,
                                 HOM_STD (v9 alias for HOM_REF)}
```

Examples: `inv_LG12_17_HOM_REF`, `inv_LG12_17_HET`, `inv_LG12_17_HOM_INV`.

These are the 90% case. Three (sometimes six with recombinants) groups per
candidate. Written by `STEP_C01i_d_seal.R` at group-validation seal time.

### Sub-cluster groups (when sub-structure is detected inside a karyotype)

When the PCA inside a karyotype shows multiple distinct clusters (see e.g.
Figure [cand46] — `HOM_INV` splits into 3–4 vertical sub-groups), each
sub-cluster gets its own registered group. Naming follows a two-segment
suffix: `__<dimension><N>` or `__<dimension>_<label>`.

The sub-cluster integer `<N>` is the dbscan cluster id (or equivalent
clustering algorithm output). It is stable within a run but carries no
geometric meaning (it is not "top" or "bottom"); it is just a group-by
bucket. Samples that don't cluster cleanly are registered in the
`__noise` bucket for that karyotype.

```
inv_<cid>_<KARYOTYPE>__sub<N>             — dbscan cluster member (N = dbscan id)
inv_<cid>_<KARYOTYPE>__noise              — dbscan "noise" (didn't cluster)
inv_<cid>_<KARYOTYPE>__ancestry<K>_<k>    — split by ancestry at K clusters
inv_<cid>_<KARYOTYPE>__ghsl_<band>        — split by GHSL haplotype band
inv_<cid>_<KARYOTYPE>__family_<fid>       — split by family_id
```

Examples:
- `inv_LG12_17_HOM_INV__sub1`, `inv_LG12_17_HOM_INV__sub2` — two dbscan sub-clusters found inside HOM_INV
- `inv_LG12_17_HOM_INV__noise` — samples inside HOM_INV that didn't assign to any sub-cluster
- `inv_LG12_17_HOM_INV__ancestry8_3` — samples that are in HOM_INV *and* ancestry_K8_Q3
- `inv_LG12_17_HET__ghsl_B2` — samples that are HET *and* in GHSL band B2

The key property: **every sub-cluster group is a subset of its parent
karyotype group**. An integrity check can verify this:
`get_group("inv_<cid>_<K>__sub<N>") ⊆ get_group("inv_<cid>_<K>")`.

### Dimension convention

Every sample-group write can optionally set the `dimension` field in
`sample_groups.tsv`. Current values:

| dimension | meaning |
|---|---|
| `.` (default) | top-level (cohort, pruned, karyotype) |
| `karyotype` | top-level karyotype group for a candidate |
| `karyotype_subcluster` | sub-cluster inside a karyotype (from recursive PCA) |
| `ancestry` | ancestry cluster (from NGSadmix / local Q) |
| `ghsl` | GHSL haplotype band |
| `family` | pedigree family |
| `intersect` | intersection of two dimensions (e.g. karyotype × ancestry) |

The dimension is metadata only — it doesn't change how `put_pairwise` or
`ask()` behave. It's there so that classification scripts can filter with
e.g. `reg$samples$list_groups()[dimension == "karyotype_subcluster"]`.

### Querying for sub-clusters

A classification script iterating over a candidate asks the database
"what sub-structure has been detected?" with one call. The natural shape
is a new helper (will be added in a follow-up): `reg$samples$get_subgroups(cid)`
returning a data.table of `(group_id, karyotype, dimension, parent_group, n)`
so the script can see all sub-clusters at once.

Until that helper exists, the same question answers in a few lines:

```r
all_groups <- reg$samples$list_groups()
cand_prefix <- sprintf("inv_%s_", cid)
cand_groups <- all_groups[grepl(paste0("^", cand_prefix), group_id)]
# Sub-clusters have a __ separator after the karyotype
cand_groups$is_subcluster <- grepl("__", cand_groups$group_id)
```

### What this means for the 10% case in your image

For candidate 46 with 3-sub-cluster HOM_INV + 3-sub-cluster HET:

1. `STEP_C01i_d_seal.R` writes the top-level karyotype groups as usual
2. A follow-up step (`STEP_C01i_e_subcluster.R`, not yet implemented)
   runs recursive PCA inside each karyotype and calls
   `reg$samples$add_group("inv_cand46_HOM_INV__sub1", ...)` for each
   detected sub-cluster
3. Classification / population-genetics scripts query the sub-clusters,
   compute pairwise FST between them (sub1 vs sub2 inside HOM_INV), and
   put the results into the manifest
4. The manifest row for each FST carries both sub-cluster names as
   `group_1` / `group_2`, so the provenance is explicit without the
   database needing a "sub-cluster" concept baked in

**No schema changes, no special-casing in put_*, no new manifest fields.**
The sub-cluster story is purely a naming convention + an optional metadata
column. The detection logic (when, how, how many sub-clusters) is a
pipeline concern, and can evolve independently of the database.

### The remaining 1% — nested inversions, chromosome-scale complications

Nested inversions, translocations, or candidates where even sub-clusters
have further structure are handled by the same pattern: deeper suffixes
(`__sub1_a`, `__sub1_b`) or a different dimension tag. The database
doesn't limit the depth. If a future candidate genuinely breaks the
convention, the manifest still records exactly which group-pair was
compared, so reproducibility is preserved even when the biology is
messy.

### Multi-scale karyotype — same sample, different label at different intervals

Rare but real: a sample is `HOM_INV` for a parent candidate (say a 7.5 Mb
inversion `LG12_17`) but is effectively `HOM_REF` inside a nested
double-crossover interval within that parent (say a 50 kb DCO tract
`LG12_17_DCO_3`). Both labels are true simultaneously, at different
coordinate resolutions. The sample's "karyotype" depends on which
interval you ask about.

**The database already supports this natively via nested candidates.**
`interval_registry.candidate_intervals` has a `parent_id` column (chat
11 addition). A DCO sub-interval is registered as its own candidate
with `parent_id = LG12_17`:

```r
# Parent candidate
reg$intervals$add_candidate(
  candidate_id = "LG12_17",
  chrom = "C_gar_LG12", start_bp = 12510000, end_bp = 20030000,
  scale = "100")

# Nested DCO sub-interval
reg$intervals$add_candidate(
  candidate_id = "LG12_17_DCO_3",
  chrom = "C_gar_LG12", start_bp = 14480000, end_bp = 14530000,
  scale = "50kb",
  parent_id = "LG12_17")

# The same sample belongs to BOTH karyotype groups — one per candidate,
# no ambiguity because each group_id is namespaced by its candidate_id:
reg$samples$get_group("inv_LG12_17_HOM_INV")         # includes CGA042
reg$samples$get_group("inv_LG12_17_DCO_3_HOM_REF")   # also includes CGA042
```

No schema changes. No new fields. No alternate-label column. Each
(sample, candidate) pair has its own karyotype label independently;
the coordinate hierarchy in interval_registry (via `parent_id`) tells
you how the candidates relate.

### Segmented karyotypes for recombinant candidates

The DCO example above has ONE nested sub-interval, but the correct
full picture for a recombinant carrier is a **complete segmentation**
of the parent interval: every genomic position in the parent maps to
exactly one segment, and each segment has its own karyotype group per
sample. This is the standard "piecewise genotype" representation used
in linkage and ancestry-painting workflows.

Worked example — sample CGA042 is HOM_INV across `LG12_17` (12.5–20.0 Mb)
except for a 50 kb double-crossover window (14.48–14.53 Mb) where the
inversion is reverted. The correct registration:

```r
# The parent candidate — ALWAYS registered (cohort-level truth)
reg$intervals$add_candidate("LG12_17",
  chrom = "C_gar_LG12", start_bp = 12500000, end_bp = 20030000,
  scale = "100", parent_id = ".")

# Piecewise segmentation — three child intervals covering the parent.
# Scale tag "seg" identifies these as recombinant-analysis segments so
# consumers can distinguish them from other nested inversions.
reg$intervals$add_candidate("LG12_17_seg_L",
  chrom = "C_gar_LG12", start_bp = 12500000, end_bp = 14480000,
  scale = "seg", parent_id = "LG12_17")
reg$intervals$add_candidate("LG12_17_seg_DCO1",
  chrom = "C_gar_LG12", start_bp = 14480000, end_bp = 14530000,
  scale = "seg_dco", parent_id = "LG12_17")
reg$intervals$add_candidate("LG12_17_seg_R",
  chrom = "C_gar_LG12", start_bp = 14530000, end_bp = 20030000,
  scale = "seg", parent_id = "LG12_17")

# Karyotype memberships — ALL of these simultaneously for CGA042:
# 1. Cohort level (parent truth — sample carries the inversion overall)
reg$samples$add_group("inv_LG12_17_HOM_INV",              c(..., "CGA042"))
# 2. Also flagged as a recombinant at the parent scale (chat-12 R∧G gate)
reg$samples$add_group("inv_LG12_17_RECOMBINANT_DCO",       c(..., "CGA042"))
# 3. Segment-level truth (local arrangement in each piece)
reg$samples$add_group("inv_LG12_17_seg_L_HOM_INV",         c(..., "CGA042"))
reg$samples$add_group("inv_LG12_17_seg_DCO1_HOM_REF",      c(..., "CGA042"))
reg$samples$add_group("inv_LG12_17_seg_R_HOM_INV",         c(..., "CGA042"))
```

CGA042 is in **all five** groups. No contradiction: membership in the
parent HOM_INV group means "carries the inversion across the cohort-level
candidate"; membership in a segment HOM_REF group means "the local
arrangement in this segment is reference." Both are true.

### Group count is bounded — this is not an explosion

For a candidate with `n_seg` segments (number of segments = number of
crossover breakpoints + 1), the maximum group count per sample is:

| recombination pattern | segments | max carrier groups |
|---|---|---|
| no recombinants (90% of candidates) | 1 (= parent) | 3 (HOM_REF / HET / HOM_INV) |
| single crossover | 2 | 3 parent + 2×3 = 9 |
| double crossover (your example) | 3 | 3 parent + 3×3 = 12 |
| triple crossover (rare) | 4 | 3 parent + 4×3 = 15 |

Plus the `RECOMBINANT` / `RECOMBINANT_GC` / `RECOMBINANT_DCO` parent-level
markers (3 more). In the worst case a candidate has ~20 groups. That's
not an explosion — that's well within the scan budget of `list_groups()`.

### Which segment scheme to use — do we always segment, or only when needed?

The pipeline decides. Two reasonable conventions:

- **Always register the parent groups.** Every candidate gets
  `inv_<cid>_<KARYO>` groups at parent scale (chat-12 seal output).
  This is the cohort-level truth and is the primary input to
  classification.
- **Register segment groups only for samples with detected recombinant
  tracts.** If no recombinant breakpoints were detected for candidate
  LG12_17 (90% case), skip the segmentation — it would add nothing.
  Only the ~10% of candidates with recombinant carriers get
  `_seg_L` / `_seg_DCO<N>` / `_seg_R` sub-intervals, and only for those
  specific samples.

This second rule keeps `list_groups()` small for the common case while
preserving full resolution for the cases that need it. The recombinant
detector (`lib_recomb_combination.R::derive_R_from_regime`) already
identifies which samples have crossover tracts; the segmentation step
calls `add_candidate` + `add_group` only for those.

**For classification: "what is this sample at this coordinate?"**
Walk up the candidate hierarchy from the smallest containing interval.
The helper `reg$query$effective_karyotype(sample_id, chrom, pos)`
returns the most specific label plus the parent chain, e.g.
`list(candidate_id="LG12_17_seg_DCO1", karyotype="HOM_REF",
      parent_candidate="LG12_17", parent_karyotype="HOM_INV")`.
The caller decides which scale they need. If no segments are registered
(90% case), the helper just returns the parent-level karyotype.

### Multi-window scans (FST / Δ12 / ENA / expH across an inversion)

For a window-by-window scan across an inversion with recombinant
segments, the correct group at each window is the karyotype group
defined at the SMALLEST candidate containing that window. For windows
outside any recombinant segment, that's the parent candidate. For
windows inside a DCO tract, that's the segment candidate — which
automatically moves CGA042 from "HOM_INV side" to "HOM_REF side" for
just those windows.

**Primitive for per-window candidate resolution:**
```r
for (pos in window_centers) {
  cid <- reg$intervals$resolve_smallest_at("C_gar_LG12", pos)
  g1  <- reg$samples$get_group(sprintf("inv_%s_HOM_REF", cid))
  g2  <- reg$samples$get_group(sprintf("inv_%s_HOM_INV", cid))
  r   <- reg$compute$pairwise_stat(pos = pos, flank_size = 5000,
                                    group1 = g1, group2 = g2,
                                    stat = "fst", chrom = "C_gar_LG12",
                                    cache = TRUE)
  # manifest row lands with candidate_id=cid, so downstream queries
  # know exactly which scale the FST was computed at
}
```

`resolve_smallest_at(chrom, pos)` returns the candidate_id of the
smallest registered interval containing pos (e.g. `LG12_17_seg_DCO1`
for windows inside the DCO, `LG12_17_seg_L` for windows in the left
segment, or `LG12_17` itself for windows in non-segmented zones). The
scan driver doesn't need to know about segment boundaries — the geometry
is handled by interval_registry.

A higher-level convenience `reg$compute$scan_pairwise(chrom, s, e,
window_size, step, stat, karyo_pair)` that wraps this loop (with
empty-group fallback, output shaped as a plotting-ready data.table,
manifest rows stamped per window) is not yet built — see chat-17+
backlog. Until then, the 5-line loop above is the standard pattern.

**Why this is better than an "alternate label" column.** An `alt_karyotype`
field on sample_groups would make ambiguity a property of the group
row, conflating the thing being labeled (sample × interval) with the
label itself. Keeping multi-scale labels as independent group
memberships per candidate preserves the invariant that sample_groups is
a clean `(group_id, sample_id)` relation. It also scales naturally:
a triply-nested event (parent → child → grandchild) is just three
candidates with three chains of karyotype groups; no special casing.

**Uncertainty is separate.** Whether we're *confident* about a sample's
karyotype is a different concern from what the karyotype is. Confidence
lives in `evidence_registry` keys (`q6_group_validation` ∈ {NONE,
UNCERTAIN, SUPPORTED, VALIDATED}; `q2_pca_cluster_separation_flag` ∈
{clean, noisy}). A sample labeled `HOM_INV` with
`q6_group_validation = UNCERTAIN` is registered normally in
`inv_LG12_17_HOM_INV`, but downstream scripts filtering by validation
level can exclude it. Two separate axes: the label and the confidence
in the label.

---

## Per-candidate folder layout — flat filesystem, tree on demand

The filesystem layout for per-candidate artefacts (Tier-2 evidence
blocks, figures, per-candidate numerical caches) is **flat**:

```
evidence_registry/per_candidate/
├── LG12_17/                      ← parent inversion
│   ├── internal_dynamics.json
│   ├── recombinant_map.json
│   └── figures/
├── LG12_17_seg_L/                ← segment, sibling in filesystem
│   └── ...
├── LG12_17_seg_DCO1/
│   └── ...
└── LG12_17_seg_R/
    └── ...

results_registry/candidate/
├── LG12_17/
│   ├── Q_K08.all_226.tsv.gz
│   └── F_K08.tsv.gz
├── LG12_17_seg_L/
│   └── ...
└── ...
```

Parent candidates, nested children, and segment children are all
siblings at the top level of `per_candidate/` and `candidate/`. The
hierarchy is **defined by `parent_id` in interval_registry**, not by
filesystem nesting.

### Why flat filesystem + tree-on-demand

An earlier design considered nesting children inside their parent's
directory (`LG12_17/children/LG12_17_seg_DCO1/`), which matches the
visual hierarchy but creates two real problems:

1. **Re-classification requires moving files.** If later analysis
   decides "LG12_17 is actually three independent systems, each with
   four inversions inside", the 1→3→12 restructuring would require
   physically moving 12 directories. Every hard-coded path in
   downstream scripts would break.
2. **Parent lookup per access.** Every `per_candidate/<cid>/` read
   would need a `parent_id` walk first to construct the path. One
   extra registry call per file access.

Under the flat convention, re-classification is a single edit: change
the 12 rows' `parent_id` values in `candidate_intervals.tsv`. Files
don't move. Plot assembly code doesn't change. Manifest rows stay
valid.

### Querying the tree

Tree views are built on demand by `reg$intervals$get_tree(root_cid)`:

```r
tree <- reg$intervals$get_tree("LG12_17")
tree$root                   # "LG12_17"
tree$all_cids               # ["LG12_17", "LG12_17_seg_L", "LG12_17_seg_DCO1", "LG12_17_seg_R"]
tree$by_depth[[1]]          # ["LG12_17"]
tree$by_depth[[2]]          # the three segments
tree$edges                  # parent / child columns for drawing the tree
```

For the "1 → 3 systems × 4 inversions each" case, after re-classification
`get_tree("LG12_17")` returns a 3-level tree with 16 candidates total.
No file moves — just the tree semantics change.

### Figure assembly under this convention

Plot scripts that need a composite figure for a candidate system query
the tree, iterate `all_cids`, and pull per-candidate panels from the
flat directories:

```r
tree <- reg$intervals$get_tree("LG12_17")
panels <- lapply(tree$all_cids, function(cid) {
  readRDS(file.path(evidence_dir, "per_candidate", cid, "panel.rds"))
})
# Assemble with figrid using tree$by_depth to decide rows/columns
```

If the biology view changes (a candidate gets promoted to a
three-system parent, or demoted to a leaf), only the `parent_id` in
interval_registry changes; the panel RDS files, the results_registry
rows, and the plot assembly code all stay valid.

### What if a candidate disappears entirely?

If a candidate is ruled out (e.g. demoted to spurious after more
analysis), the recommended action is **not** to delete its directory.
Instead:

1. Set the candidate's tier to the spurious band in its evidence block
2. Let `reg$results$integrity_check()` continue to pass (the files are
   still valid caches of prior computation)
3. Downstream classification scripts filter by tier before rendering

Physical deletion should wait until a periodic archive sweep, not
happen mid-analysis. This matches the posture on segment groups —
register defensively, filter at read time.

### Reviewing 200 inversions without clicking into 200 folders

Practical workflow problem: canonical storage is
`per_candidate/<cid>/figures/` (one folder per candidate), but for a
review pass you want to flip through all 200 candidates with left/right
arrow keys in Windows Photos or a PDF viewer, without opening each
folder.

**Solution: store organized, browse consolidated.** Keep the per-candidate
folders as canonical storage; generate browsing-friendly consolidated
outputs alongside. Two complementary outputs:

1. **Flat figure gallery** (`figures_gallery/`) — one PNG per candidate,
   copied (not symlinked — works over SMB/SSHFS to Windows), named
   `<cid>__summary.png`. Browseable with left/right arrows in Windows
   Photos or any file manager. Regenerated by a sweep script when new
   figures are produced.

2. **Consolidated PDF atlas** (`inversion_atlas_<date>.pdf`) — one page
   per candidate, multi-panel layout per page (PCA + heterozygosity +
   DBSCAN + marker heatmap + FST scan). Searchable by cid, bookmarkable,
   shareable as one attachment. Typical use: review session with
   advisor, supplementary material for manuscript.

Both are **derived outputs**, not canonical. The canonical per-candidate
folder stays the single source of truth; gallery and atlas are
regenerated from it.

```
evidence_registry/per_candidate/                ← canonical, per-cid
├── LG12_17/figures/summary.png
├── LG12_17_seg_DCO1/figures/summary.png
└── ... (200 candidate folders)

figures_gallery/                                ← derived, flat
├── LG12_17__summary.png                        ← copy of per_candidate/LG12_17/figures/summary.png
├── LG12_17_seg_DCO1__summary.png
└── ... (200 PNGs in one folder, arrow-browseable)

figures_atlas/
└── inversion_atlas_2026-04-20.pdf             ← derived, 200 pages
```

Implementation is one sweep script (~30 lines + one ggplot pipeline):

```r
# sweep_figures_gallery.R
tree <- reg$evidence$list_candidates()
for (cid in tree) {
  src <- file.path(per_candidate_dir, cid, "figures", "summary.png")
  if (file.exists(src)) {
    file.copy(src, file.path(gallery_dir, sprintf("%s__summary.png", cid)),
              overwrite = TRUE)
  }
}
```

The PDF atlas uses `figrid` (already in the toolkit) to compose each
page from per-candidate panels, then `pdf()` / `ggsave(device="pdf",
onefile=TRUE)` to concatenate pages. Sort order is configurable —
typically by chrom then start_bp, but a tier-sorted or "needs review
first" atlas is a one-line sort change.

**Helper for gallery scripts:** `reg$query$candidate_figure_paths(cid,
kind = "summary")` — returns canonical paths for standard figure kinds
(summary, dbscan, heatmap, scan_fst, scan_delta12). Makes the sweep
script portable across changes in directory structure.

This is not yet implemented — see `SPEC_DEFERRED.md` § "Figure gallery
+ atlas" for the full contract. It's a chat-18 feature (needs real
figures in hand to decide the page layout).

---

## Where TSVs live — three categories, three homes

Not all tsv output lives in the registry. There are three distinct
categories, each with its own root directory:

### Category 1: precomp artefacts (per-chromosome, per-method)

Written by precomp scripts into their own output directories. Not
registered — they are intermediate lookup tables that downstream
scripts read directly.

```
<BASE>/inversion_localpca_v7/06_mds_candidates/snake_regions_multiscale/precomp/
├── inv_likeness.tsv.gz                    ← STEP_C01a
├── windows_master.tsv.gz                  ← STEP_C01a (window registry)
├── seeded_regions_windows_<chr>.tsv.gz    ← STEP_C01b
├── landscape/
│   ├── block_registry_<chr>.tsv.gz        ← PHASE_01C_block_detect
│   ├── boundary_catalog_<chr>.tsv.gz
│   ├── blue_cross_verdicts_<chr>.tsv.gz
│   ├── block_concordance_<chr>.tsv.gz
│   └── 01C_window_pa.tsv.gz               ← per-window block membership
```

**None of these are in the registry.** They feed into C01d, which is
the candidate-birth step that registers the surviving candidates.

### Category 2: ancestry cache (per-chromosome, per-K)

```
<BASE>/ancestry_cache/
├── K02/
│   ├── C_gar_LG01.all_226.local_Q_summary.tsv.gz
│   ├── C_gar_LG01.all_226.local_Q_samples.tsv.gz
│   └── ... × 28 chromosomes
├── K03/
├── ... (K02 through K20)
```

Written by `LAUNCH_instant_q_precompute.slurm`. The K=8 canonical set
is read into precomp RDS by `STEP_C01a_precompute.R`. Intermediate —
not in the registry directly, but accessed by `reg$compute$ancestry_*`
methods on demand.

### Category 3: registry data

```
<BASE>/registries/data/
├── sample_registry/    ← WHO (226 samples + groups)
├── interval_registry/  ← WHERE (candidates + segments, with parent_id)
├── evidence_registry/  ← WHAT-SCALAR (Tier-2 blocks + keys.tsv)
└── results_registry/   ← WHAT-NUMERICAL (FST tracks, Q/F matrices)
```

Everything queryable through `reg$*` lives here. Intermediate
artefacts stay in categories 1 and 2 and are read by registry compute
methods when needed (paths stored in evidence blocks).

**Which scripts write to interval_registry?** Only two (audited chat 16):
- `STEP_C01d_catalog_birth.R` — the main candidate set
- `STEP_C01i_b_multi_recomb.R` — segments for recombinant candidates
  (chat 16 addition, with sanity guard at MAX_N_SEGMENTS_HARD=10)

No other script registers candidates. This keeps interval_registry
bounded to `n_main_candidates + Σ n_segments_per_recombinant` — a few
hundred rows total for a 28-chrom cohort.

---

## Why this is different from chat-15's `stats_cache`

Chat 15 treated pairwise FSTs and candidate Q/F matrices as "cached
numerical results that used to evaporate." That framing got us a flat
folder with opaque filenames and ad-hoc sample-set identifiers. The design
was correct at the file-content level (Q/F tables have the right shape)
but incorrect at the database level (no schema for the manifest, no FKs,
no query plane, no integrity check).

This rewrite promotes the cache to a **first-class registry**:

- It has a schema (`result_row.schema.json`) that any language can validate against.
- Its identifiers are FKs into the existing three registries, not ad-hoc hashes.
- It has a query plane (`reg$ask()`) and an integrity check (`reg$integrity_check()`).
- Its files are named using human-readable, registered identifiers — so reviewers can audit by reading the file tree, not by decoding hashes.

The result: when the pipeline writes an FST track, that FST track is
visible in one SELECT-style query; its provenance is one row; its
integrity is verifiable in one command. That's what "coherent database"
means for this paper.
