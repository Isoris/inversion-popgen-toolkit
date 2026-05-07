# Atlas — merged tree, 2026-05-06 evening

This is the canonical merged Atlas tree as of 2026-05-06 evening. It
combines four separate tarballs from the morning's parallel work
streams plus the evening's band-tracking validation and arrangement
spec.

**Read first**: this README, then `HANDOFF_2026-05-06_evening.md`,
then `HANDOFF_MERGE.md` (the sub-atlas merge work, still pending).

## Test status

```
FOUNDATION:    455 / 455 ✓   (shared/ + modular smoke + legacy parity)
BATCHES:       194 / 194 ✓   (per-page modules from the 5-batch split)
BAND_TRACKING:  88 /  88 ✓   (this environment's subset; the full count
                              is higher when upstream tarballs are layered)
─────────────────────────────────
TOTAL HERE:    737 / 737 ✓
```

## What this tree contains (and what it doesn't)

### What's present

```
.
├── HANDOFF_2026-05-06_evening.md        ← latest (this chat) — the 60-sec orientation
├── HANDOFF_MERGE.md                     ← the page-split merge work plan (327 lines)
├── HANDOFF_SUPPLEMENT.md                ← corrections to HANDOFF_MERGE (page5 isn't multi-species, etc.)
├── README.md                            ← this file
│
├── Atlas/
│   ├── inversion_review.html            ← top-level review entry
│   ├── inversion_review/                ← review-stage pages (4, 6, 7, 11, sv_evidence)
│   ├── inversion_discovery/             ← discovery-stage pages (1, 2, 8, 12, 15, 19)
│   ├── inversion_catalogue/             ← catalogue-stage pages (3, 9, 10, 17, 18, 21, overview)
│   ├── inversion_comparative/           ← comparative-stage pages (5*, 16, 16b)  *=help, see SUPPLEMENT
│   │
│   ├── shared/                          ← MODULAR CORE (locked / read-only per HANDOFF_MERGE)
│   │   ├── contingency.js, hungarian.js, het_rate.js, kmeans.js
│   │   ├── per_l2_cluster.js, state.js, state_io.js
│   │   └── band_tracking/               ← BAND-TRACKING LAYER (incomplete — see below)
│   │       ├── vote_evidence.js         ← *PATCHED 2026-05-06: opts.band_groups
│   │       ├── partition_consensus.js   ← *PATCHED 2026-05-06: best-anchored multi-layer
│   │       ├── band_voters.js           ← Layer A (per-band view)
│   │       ├── partition_enumerate.js   ← bruteforce partition enumeration + scoring
│   │       ├── index.js                 ← production re-export surface
│   │       │
│   │       ├── projection.js            ← STUB (PATTERN_CLASS enum only) — gets overwritten when
│   │       │                              you layer Atlas_band_layers_2026-05-05.tar.gz on top
│   │       ├── projection_STUB.js       ← duplicate of projection.js, kept as a label-of-warning
│   │       ├── index_min.js             ← driver-only re-exports (skips the upstream-tarball modules)
│   │       └── _kmeans_imported.js      ← copy of ../kmeans.js for standalone driver use
│   │
│   ├── tests/                           ← per-page + foundation + band_tracking tests
│   │   ├── test_shared_*.js (7)         ← READ-ONLY foundation
│   │   ├── test_modular_smoke.js, test_legacy_parity.js  ← READ-ONLY
│   │   ├── test_<phase>_page<N>.js (21) ← per-page from the 5 batches
│   │   ├── test_band_consensus.js       ← *UPDATED 2026-05-06: +F11, +F12
│   │   └── fixtures/                    ← NEW 2026-05-06
│   │       ├── synth_LG28_kmeans.json   ← K=3 LG28-shape synthetic (60/106/60)
│   │       ├── synth_LG28_projections.json   ← K=2 reference shape
│   │       └── build_synth_kmeans.mjs
│   │
│   ├── drivers/                         ← NEW 2026-05-06 — standalone runners
│   │   ├── run_consensus.mjs            ← generic 3-mode driver (voterecords/projections/kmeans)
│   │   └── run_consensus_precomp.mjs    ← reads precomp JSON directly
│   │
│   ├── data/
│   │   ├── precomp/                     ← real LANTA precomp JSONs
│   │   │   ├── LG28.json                ← 4302 windows × 226 samples × 48 L2 envelopes
│   │   │   ├── catfish_226_relatedness.json
│   │   │   ├── cohort_diversity_v1.json
│   │   │   ├── cs_breakpoints_v1.json
│   │   │   ├── LG28.repeat_density.scrubber_windows.json
│   │   │   └── sv_genotype_counts/
│   │   │
│   │   └── arrangement_calls/           ← NEW 2026-05-06 — real-LG28 consensus outputs
│   │       └── lg28_2026-05-06_run/
│   │           ├── CLEAN_envs_02-03-04_consensus.json                ← K=3 recovered, score 0.867
│   │           ├── CLEAN_envs_03-04-05_acrossanomaly_consensus.json  ← spans the 91/67/68 envelope
│   │           ├── CLEAN_envs_04-05-06_consensus.json
│   │           └── LABEL_AMBIGUOUS_envs_01-02_consensus.json         ← F12 case
│   │
│   ├── specs_todo/
│   │   ├── SPEC_arrangement_color_mode_and_arrangement_calls_v1.md   ← NEW 2026-05-06 (main)
│   │   ├── SPEC_busco_4d_age_brackets.md                              ← from morning batch
│   │   └── SPEC_inversion_age_atlas_surface_AMENDMENT.md              ← from morning batch
│   │
│   ├── findings/                        ← NEW 2026-05-06
│   │   └── FINDING_2026-05-06_both_excluded_bug.md
│   │
│   ├── docs/
│   │   ├── MIGRATION_INVENTORY.md
│   │   ├── MIGRATION_LOG.md
│   │   └── merge_inputs/                ← TODO inventories for the page-split merge work
│   │       ├── TODO_MISSING_inventory.txt
│   │       ├── TODO_PROMOTE_inventory.txt
│   │       ├── TODO_SLOTS_inventory.txt
│   │       └── all_todo_locations.txt
│   │
│   ├── dist/
│   │   ├── inversion_review_flat.html   ← pre-built artifact
│   │   └── Atlas-inversion-2026-05.tar.gz
│   │
│   ├── build/
│   │   ├── flatten.py
│   │   └── package_workflow.py
│   │
│   ├── README.md                        ← original Atlas/ README
│   └── (empty: catalogue/, comparative/, review/, scrubber/)
│
├── handoff_docs/                        ← prior handoffs from morning chats
│   ├── HANDOFF_2026-05-05_turn164_lasso_linkage.md
│   ├── HANDOFF_2026-05-05_turn165_g_panel_auto_tab.md
│   ├── HANDOFF_2026-05-05_turn166_spec_extension_round2.md
│   ├── HANDOFF_2026-05-06_morning_age_specs.md      ← from the spec-batch chat
│   └── SPEC_band_track_extraction_and_l3_single_band_rows.md
│
└── legacy/
    └── Inversion_atlas.html             ← turn-166 monolith, 75617 lines, 3.5 MB
                                          kept for parity testing only
```

### What's MISSING (Quentin has these on his end)

`Atlas/shared/band_tracking/` is **incomplete**. The five files
present here come from `Atlas_band_consensus_2026-05-05.tar.gz`
(the 5th of 5 additive tarballs in the morning handoff's
five-tarball progression). The four earlier additive tarballs ship
the upstream modules that `index.js` re-exports from but that aren't
in this tree:

| Tarball | Files it ships |
|---|---|
| `Atlas_FINAL_2026-05-05.tar.gz` (#1) | `single_band.js`, `het.js`, `hom.js`, `iv.js` |
| `Atlas_band_tracking_2026-05-05.tar.gz` (#2) | (additions on top of #1: single-band tracking + het skeleton + iv calls) |
| `Atlas_band_layers_2026-05-05.tar.gz` (#3) | `trajectory.js`, **`projection.js`** (real, not stub), `karyotype_model.js` |
| `Atlas_band_consensus_2026-05-05.tar.gz` (#4) | `vote_evidence.js`, `band_voters.js`, `partition_enumerate.js`, `partition_consensus.js`, `index.js` ← present in this tree |

To get a complete `Atlas/shared/band_tracking/`, layer the four
prior tarballs in order, **then** unpack this merged tree on top.
Since this tree's evening-patched `vote_evidence.js` and
`partition_consensus.js` need to win, the order matters: prior
tarballs first, this tree last.

### Two cleanups after layering

1. The driver-only files exist because this tree's `index.js`
   re-exports from upstream modules that aren't here. Once layered:
   - **Delete** `projection_STUB.js`, `index_min.js`,
     `_kmeans_imported.js` (the stub `projection.js` is overwritten
     by the real one from tarball #3 — that's correct, leave the
     replaced file alone).
   - **Update** `Atlas/drivers/run_consensus*.mjs` import lines from
     `./Atlas/shared/band_tracking/index_min.js` to
     `./Atlas/shared/band_tracking/index.js`.

2. The page-split merge work (per `HANDOFF_MERGE.md`) hasn't started
   yet. There are ~197 TODO markers across page modules referencing
   functions still in the legacy monolith, plus 31 missing state
   slots, plus 2 promote-to-shared candidates. The merge protocol in
   HANDOFF_MERGE.md walks through resolving them. This work is
   **independent** of the band-tracking work and can proceed in
   parallel.

## Setup checklist for a working LANTA / desktop tree

1. Extract this tree to your work directory.
2. Layer the four prior `Atlas_*_2026-05-05.tar.gz` band_tracking
   tarballs (in numeric order, oldest first), each over the previous.
3. Unpack THIS tree on top — the `vote_evidence.js`,
   `partition_consensus.js`, `test_band_consensus.js` overwrites
   land here.
4. Delete the three driver-only helper files (see above).
5. Update the driver imports (see above).
6. Run all three test groups:
   ```bash
   cd Atlas
   # foundation
   for t in tests/test_shared_*.js tests/test_modular_smoke.js tests/test_legacy_parity.js; do
     LEGACY_ATLAS=$PWD/../legacy/Inversion_atlas.html node "$t"
   done
   # batches
   for t in tests/test_discovery_*.js tests/test_review_*.js tests/test_catalogue_*.js tests/test_comparative_*.js; do
     node "$t"
   done
   # band-tracking (full — requires upstream modules to be present)
   node tests/test_band_consensus.js
   ```
   Expected: 455 + 194 + (≥88, exact number depends on how many
   fixture-asserts the full upstream layer adds — morning handoff
   cited 991 for the cumulative band_tracking suite).

## The 2026-05-06 evening patches (what changed)

Two patches landed in the band_tracking layer this afternoon:

### `vote_evidence.js`
- New optional `opts.band_groups` parameter on `extract_votes` —
  integer group id per band.
- `build_coassociation_matrix` reads it: when two bands belong to
  the same group, the "both excluded → together" rule is suppressed
  (counted as `apart`). Fixes the K_arrangements > 2 over-counting
  bug.
- Backwards compatible: omit `band_groups` → original behavior.
- See `Atlas/findings/FINDING_2026-05-06_both_excluded_bug.md`.

### `partition_consensus.js`
- Multi-layer detection now compares each picked partition against
  the **best** partition (testing whether refinement holds), not
  arbitrary pairs of siblings. Pairs of mutually-incompatible
  siblings that are both refinements of the best are no longer
  multi-layer evidence — they're noise neighbors of one truth.
- Adds a reasoning trace line documenting whether `band_groups` was
  applied.

### `test_band_consensus.js`
- **F11** — K=3 arrangements with `band_groups` correction
  (regression test for the bug fix). Without `band_groups` →
  score ≤ 0.78, MULTI_LAYER_STRUCTURE false positive. With
  `band_groups` → CLEAN_PARTITION at score ≥ 0.86, recovers the
  planted arrangement structure.
- **F12** — N=2 windows × K=3 arrangements label-ambiguity case.
  Pinned as a known limitation. Score=1.0, residual=0, but
  classifies MULTI_LAYER_STRUCTURE because the data is genuinely
  insufficient to disambiguate without a third anchor probe.

## Real LG28 results (the empirical confirmation)

LG28 prototype, K=3 K-means per envelope, default thresholds, with
patches in place:

| Probes | n_envs | Result |
|---|---|---|
| 2 | 2 | LABEL_AMBIGUOUS (F12 case): score 1.0, residual 0, but 6 valid label permutations tied |
| **3** | **3** | **CLEAN_PARTITION at 0.8667**, top = `[[0,3,6],[1,4,7],[2,5,8]]` ✓ |
| 4–6 | 4–6 | NO_CLEAN_CONSENSUS (K_locus exceeds bruteforce hard_fail; needs greedy fallback) |

Per-envelope cluster sizes in the prototype interval (registered
prototype is 60/106/60):

| Envelope | bp range | n_per_cluster | note |
|---|---|---|---|
| _0010_01 | 15.087–15.208 Mb | 62/102/62 | ✓ |
| _0010_02 | 15.195–16.222 Mb | 61/103/62 | ✓ |
| _0010_03 | 16.209–16.580 Mb | 62/104/60 | ✓ |
| _0010_04 | 16.569–17.473 Mb | 71/95/60  | ≈ |
| _0010_05 | 17.460–17.937 Mb | **91/67/68** | ✗ — possible breakpoint regime change |
| _0010_06 | 17.926–18.330 Mb | 76/96/54  | ≈ |

The K=3 prototype karyotype is recovered correctly across **every**
3-envelope grouping that includes envelopes 02–04 (the cleanly-shaped
core of the prototype) AND across groupings that span the
`_0010_05` 91/67/68 anomaly. The size mismatch in that envelope
reflects PC1 cluster boundary placement, not sample-identity
scrambling.

## Three open work streams for the next chat

| # | Stream | Status | Next action |
|---|---|---|---|
| 1 | **Arrangement-mode spec implementation** | spec drafted (`SPEC_arrangement_color_mode_and_arrangement_calls_v1.md`), naming decisions still open (`[NAME-1/2/3]`) | pick names → implement decoder + runner + color-mode |
| 2 | **Page-split merge work** | post-batches, pre-merge: 197 TODOs + 31 slots + 2 promotes | follow `HANDOFF_MERGE.md` step-by-step protocol |
| 3 | **Age-layer atlas patch** | spec drafted (Method 1+2+3 + Row D), `STEP_C01f_e/f` LANTA scripts not yet written | independent; can proceed in parallel |

Stream 1 (arrangement viz) is the spec deliverable from this
afternoon. Stream 2 (merge work) is queued from this morning's
parallel batches. Stream 3 (age layer) is queued from this
morning's spec batch. None blocks the others.

## clair3 status

clair3 is running on Quentin's LANTA jobs (per his message during
this chat). Do not block any of these work streams on its outputs.

End of README.
