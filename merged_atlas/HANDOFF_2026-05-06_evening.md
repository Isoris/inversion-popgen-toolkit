# HANDOFF — 2026-05-06 — band_tracking real-data validation + arrangement-mode spec

**Date**: 2026-05-06 (afternoon, post LG28 real-data run)
**Continues from**: HANDOFF_2026-05-06_age_specs_and_band_tracking.md
(this morning) — the band_tracking spec batch.
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

---

## TL;DR (60-second orientation)

This chat picked Option B from the morning handoff: run band_tracking
on real LG28. Three things happened:

1. **A bug was found and fixed.** The K=3 case (HOM_REF/HET/HOM_INV)
   classified MULTI_LAYER_STRUCTURE on real LG28 data despite being
   a clean inversion. Root cause: `build_coassociation_matrix` rule
   "both excluded → together" overcounts when K_arrangements > 2.
   Two patches landed in `vote_evidence.js` and
   `partition_consensus.js`:
   - **Patch 1** (`vote_evidence.js`): optional `band_groups` array
     — when two bands belong to the same group (e.g. two K-means
     clusters of one window/L2 envelope), "both excluded" is
     correctly counted as `apart`, not `together`. Backwards
     compatible (omit `band_groups` → original behavior).
   - **Patch 2** (`partition_consensus.js`): multi-layer detection
     compares each picked partition against the **best** partition,
     not arbitrary pairs of siblings. Pairs of mutually-incompatible
     siblings that are both refinements of the best are no longer
     flagged as multi-layer.
   Test suite went from 10 fixtures (88 checks) to 12 fixtures with
   2 new ones pinning the fix (F11) and the known label-ambiguity
   limitation (F12). 88/88 pass on this environment's subset.

2. **Real LG28 results.** With both patches, 3-envelope candidate
   regions on LG28 prototype (15.115–18.005 Mb) classify
   CLEAN_PARTITION at score 0.8667 with the correct K=3 partition
   `[[0,3,6],[1,4,7],[2,5,8]]` — recovering the registered 60/106/60
   karyotype across multiple envelope groupings. The 91/67/68
   anomaly at envelope `_0010_05` (17.46–17.94 Mb) does not break
   the call. See §3 for the table.

3. **A spec was drafted, not built.** Per Quentin's request: don't
   build a viz prototype yet, write the spec. Result:
   `SPEC_arrangement_color_mode_and_arrangement_calls_v1.md` — adds
   "arrangement" as a 9th color mode in the existing per-sample-lines
   panel infrastructure (the `_resolveSampleColorByMode` dispatch
   that's already wired into `Inversion_atlas.html`), defines the
   `arrangement_calls_v1.json` wire format, and explicitly handles
   Quentin's flag that the analytic unit isn't always L2-locked
   (per-window n=5/10/etc, custom bp ranges, half-L2 splits all
   supported via the `probe` abstraction).

**Estimated cost of implementing the spec**: ~510 LOC atlas-side +
~200 LOC tests. About 60% the size of the parallel age-layer patch.

---

## What's in this tarball

```
Atlas/
├── shared/band_tracking/
│   ├── vote_evidence.js              ← PATCHED (band_groups)
│   ├── partition_consensus.js        ← PATCHED (multi-layer detection vs best)
│   ├── band_voters.js                ← unchanged (from Atlas_band_consensus_2026-05-05)
│   ├── partition_enumerate.js        ← unchanged
│   ├── index.js                      ← unchanged (production index)
│   ├── index_min.js                  ← driver-only (skips not-shipped-here modules)
│   ├── projection_STUB.js            ← rename from projection.js — local stub only,
│   │                                    DO NOT ship over your production projection.js
│   │                                    (which is in Atlas_band_layers_2026-05-05.tar.gz)
│   └── _kmeans_imported.js           ← copy of your Atlas/shared/kmeans.js for driver use
├── tests/
│   └── test_band_consensus.js        ← +F11 (band_groups regression) +F12 (label-ambiguity pin)
└── specs_todo/
    └── SPEC_arrangement_color_mode_and_arrangement_calls_v1.md   ← MAIN DELIVERABLE

drivers/
├── run_consensus.mjs                 ← generic 3-mode driver (voterecords/projections/kmeans)
└── run_consensus_precomp.mjs         ← 4th mode: reads precomp JSON directly
                                        (chrom + windows[].pc1[] + l2_envelopes[])

findings/
└── FINDING_2026-05-06_both_excluded_bug.md   ← the bug writeup, with
                                                 synthetic + real-data confirmation

synthetic_fixtures/
├── synth_LG28_kmeans.json            ← 226 samples, 60/106/60 K=3 karyotype, 3 windows
├── synth_LG28_projections.json       ← K=2 reference (the existing test-suite shape)
└── build_synth_kmeans.mjs            ← script that built the kmeans fixture

lg28_real_data_runs/
├── CLEAN_envs_02-03-04_consensus.json       ← THE result: K=3 recovered, score 0.867
├── CLEAN_envs_03-04-05_acrossanomaly_consensus.json  ← K=3 recovered across the 91/67/68
├── CLEAN_envs_04-05-06_consensus.json       ← K=3 recovered at the breakpoint side
└── LABEL_AMBIGUOUS_envs_01-02_consensus.json ← F12 case: score=1.0 but label-ambiguous

HANDOFF_2026-05-06.md   ← this file
```

NOT in this tarball (you have these on your end, from prior tarballs):

- `Atlas/shared/band_tracking/single_band.js`, `het.js`, `hom.js`,
  `iv.js`, `trajectory.js`, `karyotype_model.js`, **`projection.js`**.
  These ship in `Atlas_FINAL_2026-05-05.tar.gz`,
  `Atlas_band_tracking_2026-05-05.tar.gz`, and
  `Atlas_band_layers_2026-05-05.tar.gz`. The `index.js` re-exports
  from them. The driver (`run_consensus*.mjs`) uses `index_min.js`
  to avoid pulling them in.
- `Atlas/shared/per_l2_cluster.js`, `kmeans.js`, `state.js`,
  `state_io.js`, `contingency.js`, etc. — your full
  `Atlas/shared/` tree from `Atlas_handoff_2026-05-05_tar.gz`.
- `Inversion_atlas.html` — the full atlas (you have it as
  `legacy/Inversion_atlas.html` and as a turn-166 version inside
  `original_uploads/Atlas_turn166_round2_2026-05-05_tar-1.gz`).

---

## Naming decisions still open (block the spec rename, but not impl)

The spec body uses provisional names tagged `[NAME-N]`. Three pending:

| Tag | Provisional | Means |
|---|---|---|
| `[NAME-1]` | **arr-cluster** | one K-means cluster within one probe |
| `[NAME-2]` | **candidate-region** | the genomic span being analyzed |
| `[NAME-3]` | **probe** | the unit of clustering (1-window / 5-window / L2 / custom range) |

`[NAME-3]` matters most: it's the new abstraction Quentin flagged
("not only per L2"). The spec uses "probe" everywhere because it
captures the role ("we're probing the sample partition at this
position") without colliding with L1/L2 envelope, locus, candidate,
band, cluster, or arrangement.

If "probe" doesn't fit your codebase's vocabulary — alternatives:
**clustering unit**, **slice**, **window-group**, **probe-region**.
One global find/replace.

---

## Real LG28 results table (key numbers for the manuscript)

All runs against `Atlas/json/LG28.json` from
`Atlas_turn166_round2_2026-05-05_tar-1.gz` (4302 windows, 226 samples,
48 L2 envelopes, schema_version 3). K=3 K-means per probe, default
thresholds.

| Probes | n_envs | result |
|---|---|---|
| `_0010_01,_0010_02` | 2 | LABEL_AMBIGUOUS (F12 case): score 1.0, residual 0, but 6 valid label permutations tied. |
| `_0010_01,_0010_02,_0010_03` | 3 | **CLEAN_PARTITION** at score 0.8667, top = `[[0,3,6],[1,4,7],[2,5,8]]` ✓ |
| `_0010_02,_0010_03,_0010_04` | 3 | **CLEAN_PARTITION** at score 0.8667, top = `[[0,3,6],[1,4,7],[2,5,8]]` ✓ |
| `_0010_03,_0010_04,_0010_05` | 3 | **CLEAN_PARTITION** at score 0.8667 (across the 91/67/68 anomaly) ✓ |
| `_0010_04,_0010_05,_0010_06` | 3 | **CLEAN_PARTITION** at score 0.8667 ✓ |
| `_0010_01..._0010_04` | 4 | NO_CLEAN_CONSENSUS (K_locus=12 exceeds bruteforce hard_fail=13) |
| `_0010_01..._0010_06` (full prototype) | 6 | NO_CLEAN_CONSENSUS (K_locus=18 > hard_fail). Per-band view still runs and shows correct cross-envelope structure. |

**Per-envelope K-means cluster sizes** (registered prototype is 60/106/60):

| Envelope | bp range | n_per_cluster | note |
|---|---|---|---|
| _0010_01 | 15.087–15.208 Mb | 62/102/62 | ✓ |
| _0010_02 | 15.195–16.222 Mb | 61/103/62 | ✓ |
| _0010_03 | 16.209–16.580 Mb | 62/104/60 | ✓ |
| _0010_04 | 16.569–17.473 Mb | 71/95/60  | ≈ |
| _0010_05 | 17.460–17.937 Mb | **91/67/68** | ✗ — possible breakpoint regime change. Sample identity preserved in band-tracking. |
| _0010_06 | 17.926–18.330 Mb | 76/96/54  | ≈ |

---

## What the next chat should do

### Priority 1 — Implement the spec

The spec is the main deliverable. Implementation order (from §10 of
the spec):

1. **Naming decisions** (`[NAME-1]`, `[NAME-2]`, `[NAME-3]`) → global
   rename across spec, FINDING, drivers, tests. ~50 LOC of churn.
2. `arrangement_decode.js` — decoder from top-partition →
   per-sample arrangement IDs. ~80 LOC.
3. `arrangement_call_runner.js` — wraps `clusterL2` per probe +
   voteRecords assembly + `consensus_partition` + decode. ~150 LOC.
4. `_LINES_COLOR_MODES` registration + `_arrangementColor` +
   `_arrangementScopeColor` + state slots. ~80 LOC into atlas.
5. Probe-strategy expander UI (page 3). ~100 LOC.
6. QC badges under candidate panel header. ~60 LOC.
7. Label-ambiguity UX (F12 case). ~40 LOC.
8. Tests. ~200 LOC.
9. Optional: LANTA-side batch step `STEP_C01f_g`. ~150 Python.

Critical path: 1→2→3→4 yields a working "arrangement" mode. 5–9 are
follow-ups.

### Priority 2 — Wire the greedy fallback

`kt_infer_macro_band_groups` is in the band_layers tarball #3 (you
have it; this chat doesn't). Threading it into `consensus_partition`
when K_locus > 10 enables the full 6-envelope LG28 prototype to
classify, not just 3-envelope subsets. Independent of the spec.

### Priority 3 — Then Option A (the age-layer patch)

The morning handoff queued the ~880-LOC age-layer patch. It's
independent of band_tracking and the spec; can ship in parallel.
Confirm whether `Atlas_turn164_2026-05-05_tar.gz` is still your
production baseline before patching.

### Priority 4 — Investigate `_0010_05` 91/67/68 anomaly

Manuscript-level question: why does K-means find 91/67/68 instead of
~60/106/60 in the 17.46–17.94 Mb sub-region? Possibilities: real
breakpoint signal (registered prototype ends at 18.005 Mb, this is
~80 kb upstream), PC1 shape change near the breakpoint, or genuine
sub-arrangement structure. Band_tracking confirms sample identity is
preserved (the 91 samples at this envelope's cluster 0 are the same
samples that are at cluster 0 in adjacent envelopes), so the
mismatch is in cluster boundaries, not assignments. Worth a focused
look during manuscript figure prep.

---

## Open decisions, ranked by urgency

### Now (block the spec implementation)

1. **Naming**: `[NAME-1]`, `[NAME-2]`, `[NAME-3]`. See top of spec.
2. **Probe-strategy default for single-L2 candidates** (spec §11.2):
   auto-subdivide or require user override?
3. **Color palette match with K-means** (spec §11.4): match first
   3 colors, or deliberately distinct?

### Soon (block atlas integration polish)

4. **LABEL_AMBIGUOUS as 7th classifier class?** (spec §11.1) — UX
   handles it without algorithm change; algorithm-level decision can
   defer.
5. **Probe-strategy persistence** (spec §11.3) — recommendation in
   spec is to record in `arrangement_calls_v1.json` itself; confirm.
6. **Behavior when consensus_class != CLEAN** (spec §11.5) — emit
   the call anyway with QC warning, or grey out? Recommendation:
   emit with QC warning.

### Later (block manuscript)

7. **Where does Method 3 BUSCO 4D pipeline run on LANTA?** (from
   morning handoff, parallel age-spec work, not affected by this
   chat).
8. **`_0010_05` anomaly investigation** (Priority 4 above).

---

## Constants worth carrying forward

(Same as morning handoff. Repeated here so the next chat doesn't have
to bounce between docs.)

### Three catfish cohorts — never conflate

1. F₁ hybrid (*C. gariepinus* × *C. macrocephalus*) — genome assembly
   paper only.
2. **226-sample pure *C. gariepinus* hatchery cohort on LANTA — current
   inversion work.** K clusters reflect hatchery broodline structure,
   NOT species admixture.
3. Pure *C. macrocephalus* wild cohort — future paper.

### KaryoConstants

- NAToRA-pruned 81 unrelated samples
- NGSadmix K=8 on pruned set
- `mergeThr 0.85`, `minNGroup 5`, `minNWin 5`, `alpha 0.001`
- `state.k 3`, `purity 0.80`, `ambiguous 0.5`
- KING/Manichaikul: 1st≥0.177, 2nd≥0.0884, 3rd≥0.0442
- Reference: `fClaHyb_Gar_LG.fa` (28 pseudochromosomes, ~964 Mb)

### LG28 prototype

- 15.115 – 18.005 Mb, 28-band assembly
- karyotype 60/106/60 (HOM_REF / HET / HOM_INV)
- shelf Fst ~0.308 vs flanks ~0.04
- spans L2 envelopes `C_gar_LG28_d17L2_0010_01` through `_0010_06`

### Mutation rate (from morning handoff, locked)

- μ_year = 3×10⁻⁹/site/year (Liu 2023 Siluriformes)
- Method 3 brackets with μ_low=1×10⁻⁹, μ_high=9×10⁻⁹

### Band-consensus algorithm thresholds (defaults)

- trajectory τ = 0.7
- projection split / merge: 0.4 / 0.85
- primary partner threshold = 0.7
- secondary partner = [0.4, 0.7]
- ambiguous-band threshold = 0.5
- δ for top-N adaptive = 0.10
- K_cap for bruteforce = 10
- hard-fail K = 13

**Real-LG28 evidence (added this chat)**: at K=3, default thresholds
work as designed once the band_groups correction is in. No
recalibration needed for the hatchery cohort.

---

## Logistical

- **clair3 is still running** on Quentin's LANTA jobs (from his
  message earlier). Do not block on its outputs for this work.
- **Inversion_atlas.html turn baseline**: Quentin had two copies in
  the prior tarball — `legacy/Inversion_atlas.html` (75617 lines)
  and a turn-166 version inside the truncated inner tarball
  `Atlas_turn166_round2_2026-05-05_tar-1.gz`. The next chat should
  confirm which is the patching target. The morning handoff said
  turn 164; this chat saw turn 166 in the tarball name.
- **Test count reconciliation**: morning handoff said "991/991"
  for the synthetic suite. This chat's local subset (using
  `index_min.js`) ran 88 checks across 12 fixtures and they all
  pass. The 991 number includes assertions across the full
  `band_tracking` + upstream modules; need to re-run on a tree with
  all tarballs unpacked to confirm 991+.

End of handoff.
