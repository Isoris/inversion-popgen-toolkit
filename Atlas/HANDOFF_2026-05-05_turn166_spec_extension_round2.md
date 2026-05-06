# HANDOFF — Spec extension: het-anchored band-track skeleton

**Date**: 2026-05-05
**Scope**: spec-only turn — extending
`specs_todo/SPEC_band_track_extraction_and_l3_single_band_rows.md`
**Atlas main file**: not modified this turn (still 76,407 LOC from
turn 165 close).
**Project**: `MS_Inversions_North_african_catfish`

This was a spec-design conversation, not a code-shipping turn. The
deliverable is a fully-fleshed-out spec that subsequent turns can
implement against.

---

## What this turn produced

The original turn-prior spec covered:
- Per-band rows under K×K contingency tables in the L3 panel
- Band-track extraction (cohort-tracking through window-to-window
  K-means transitions)
- Multi-segment / A-B-A track detection
- Four supplementary TSV outputs

Across three back-and-forth messages this turn, Quentin redirected
the spec in three substantial ways. Each redirect was additive — the
prior content stayed, new sections layered on top:

### Message 1 — *"its like the per sample lines but with lines colored by their 'intervals'"*

→ added §1.5, §7.4, §7.5

The lines panel becomes the band-track view: each fish's line
colored by which track it belongs to at each window. K-means raw
labels rotate arbitrarily, so coloring lines by raw label is
useless — coloring by track membership (which is identity-stable
across windows because tracks follow cohorts) makes the lines
panel readable.

Implements as a new value of `state.linesColorMode`:
`'band_track_all'` / `'band_track_het'` /
`'band_track_principal_only'` / `'band_track_with_het_overlay'`.
Hooks into the existing v3.99 turn 14e+ stub
`_resolveSampleColorByMode` (which has been waiting for resolvers
on its non-`family` branches).

### Message 2 — *"track the het band from dosage bc its the skeleton for our candidate inversion interval"*

→ added §4A (entire het-anchored mode), §6.4 het skeleton TSV,
§6.5 het carrier edges TSV

The het band is not just one of the tracks — it's the **principal
track**. The inversion's biological backbone IS the chain of
windows where the same heterozygous carriers persist. Symmetric
tracks (every band gets equal treatment) is now a fallback mode;
het-anchored is the default when `dosage_chunks` is loaded.

§4A.1 identifies the het band(s) per pane via a hierarchy:
`_computeHetRateForL2` (turn 129 primitive) → mean dosage →
PC1-middle heuristic. Allows zero, one, or multiple het bands per
pane (avoids forcing exactly one when complex).

§4A.3 introduces dropout-skip lookahead: forward-walk can skip up
to `MAX_DROPOUT_GAP` panes when downstream carriers re-establish
high overlap with the current band. This is the
"dropout vs boundary" distinction made operational.

§4A.4 distinguishes het-skeleton outputs from symmetric-track
outputs in `state.bandTracks[cand.id]`: a `het_skeleton` field
sits alongside the symmetric `tracks` array, with skeleton-specific
fields (consensus_carriers, carrier_stability_score,
dropout_panes, dropout_diagnosis).

§4A.5 records why this lives in the same spec rather than a
separate one — every primitive in §3-§5 generalises with one
extra metric (overlap coefficient) and one extra step
(dropout-skip).

### Message 3 — *"add this to our single band ... track the het band from dosage bc its the skeleton ... be careful with double crossovers"*

→ added §3.2, §4A.6, §4A.7, expanded §2 vocabulary

§3.2 maps the HPC-conceptual input columns (chr/start/end/sample_id/
band_label/PC1/PC2/dosage/het count/missingness/local ancestry Q)
to atlas primitives. Most are already shipped; only local-ancestry
Q is gated on a future JSON layer. The mapping table records the
equivalence so a future R-module mirror has a stable input
contract.

§4A.6 records the seven possible causes of het-signal dropout
(recombination, gene conversion, marker dropout, local low
diversity, shared ancestral haplotype, ancestry mosaicism, noisy
K-means, lack of informative SNPs). Three biological + four
technical-or-statistical. The algorithm cannot distinguish them.
The spec's operational rule: "het signal stops" is a signal
property; "inversion stops" is a biological claim.

§4A.7 records the canonical manuscript slogans:
**"dropout does not necessarily mean boundary"** (for
skeleton_with_dropout cases) and **"different downstream
carriers = different system"** (for genuine boundaries). Use in
tooltips, methods sections, figure captions.

§2 vocabulary expanded with:
- het signal stops vs inversion stops (operational vs biological)
- The two slogans above
- carrier stability score (now formally defined)

---

## Spec final shape

`specs_todo/SPEC_band_track_extraction_and_l3_single_band_rows.md`
1,769 lines, 12 top-level sections:

```
1. The problem this solves (with §1.4-1.5 het+lines integration)
2. Vocabulary discipline (with operational slogans)
3. The data flow
   3.1 Band cohorts
   3.2 Input columns (HPC ↔ atlas mapping)
4. The connection rule (per Quentin's design)
   4.1 Why not Hungarian alone
4A. Het-anchored mode (the principal-track variant)
   4A.1 Identifying the het band(s) per pane
   4A.2 Carrier-continuity edges
   4A.3 Skeleton extraction with internal-dropout tolerance
   4A.4 Distinguished outputs for het skeletons
   4A.5 Why this lives inside the symmetric-mode spec
   4A.6 Dropout causes (why "het signal stops" ≠ "inversion stops")
   4A.7 The dropout-vs-boundary slogan
5. Building band tracks from edges
   5.1 Path construction
   5.2 Continuous track summary
   5.3 Split-track / A-B-A detection
   5.4 System classification
6. Output objects (in-session + optional JSON/TSV)
   6.1 In-session state
   6.2 Output files (six TSVs A-F)
   6.3 TSV writer naming + location
7. UI surfaces in the L3 panel
   7.1 Single-band rows under each K×K table (slice 1)
   7.2 Track strip below the L3 strip (slice 2)
   7.3 System classification chip in the L3 toolbar (slice 3)
   7.4 Lines-panel coloring by track membership (slice 2)
   7.5 Het skeleton overlay on candidate strip + lines panel (slice 2)
   7.6 What this is NOT
8. Implementation slices
   Slice 1 — single-band rows + het-band identifier (~1 turn)
   Slice 2 — band-track extraction + skeleton + lines coloring + track strip (~1.5 turns)
   Slice 3 — TSV exports + classification chip + mode chip (~0.5 turn)
9. Open questions
10. Tests (slice 1, 2, 3 test plans)
11. Cross-references
12. Manuscript figure / supplementary integration
```

Total implementation estimate: 2.5–3.5 turns split across three
slices.

---

## Why this is one spec, not two

Quentin asked at one point whether to split into "atlas-side
spec" vs "HPC R-module spec." The decision recorded in §3.2 +
§4A.5: **one spec, two homes**. The atlas implementation is the
canonical algorithm in JS; an R-module mirror can consume the
same column contract (§3.2) and produce comparable outputs.
Differences only arise from threshold choices documented in TSV
headers. A separate R-module spec, if needed, would mirror this
spec rather than redefine the algorithm.

---

## What's NEXT

Implementation in atlas-turn order:

**Turn 166 (slice 1, ~1 turn):**
- `_buildBandEdgesForPanes(panes, K)`
- `_classifyHetBandsForPane(l2idx, K, opts)` — uses turn 129's
  `_computeHetRateForL2` primitive
- `_ssBandRowHtml(M, rowIdx, KA, KB, edgeStats, hetClass)`
- `_ssBandRowsForPaneHtml(ct, edgesForPair, leftK, hetClasses)`
- `state.l3SingleBandRows` toggle + localStorage
- L3 toolbar checkbox
- `renderL3Panel()` integration

**Turn 167 + 168 (slice 2, ~1.5 turns):**
- `_buildBandTracksForCandidate(cand, opts)`
- `_extractHetSkeleton(panes, K, hetClassesByPane)`
- `_detectSplitTracks(tracks, opts)`
- `_classifyBandSystem(tracks, hetSkeleton, panes)`
- `_resolveSampleColorByMode('band_track', ...)` — fills the
  v3.99 turn 14e+ stub
- Lines-panel color-mode picker gains four `band_track_*` modes
- `_drawBandTrackStrip(canvas, opts)`
- Het skeleton overlays
- `_invalidateBandTracksCache` hooks

**Turn 169 (slice 3, ~0.5 turn):**
- Six serializers (band_tracks, band_track_edges, system_band_summary,
  sample_band_paths, het_skeletons, het_carrier_edges)
- Classification chip in the L3 toolbar
- Mode chip (symmetric / het-anchored)
- TSV-export button (`📊 tracks`)

Real-LG28 review by Quentin between turns to calibrate thresholds:
- `JACCARD_MIN`, `RETENTION_MIN`, `HET_JACCARD_MIN`,
  `HET_OVERLAP_MIN`, `HET_RETENTION_MIN`
- `MAX_DROPOUT_GAP`, `DROPOUT_OVERLAP_MIN`
- `MAX_GAP_BP`, `GAP_OVERLAP_MIN`, `MIN_CORE_OVERLAP_N`

All defaults documented in the spec are starting points; Quentin
should review them on real LG28 data before considering them
locked for the manuscript.

---

## Interaction with already-shipped work

This spec extends and consumes:

- **Turn 129 (`_computeHetRateForL2`)** — the per-sample het-rate
  primitive that powers §4A.1's het-band identifier
- **Turn 130 (`_isAutoCandidate`, lineage compute)** — review-surface
  framing reused; lineage compute is the dual representation
- **Turn 141 (`_paintCandidateBands`)** — lines panel candidate-bands
  Slice 2 alpha intervals naturally compose with §7.5 het skeleton
  overlay
- **Turn 152 (G-panel inheritance tab)** — cross-candidate matrix
  consumes `_gatherActiveCandidatesForInheritance`'s filter; band
  tracks are intra-candidate, complementary
- **Turn 156 (V-shape modal)** — the modal pattern §7's
  classification chip popover should mirror
- **Turns 160-164 (band-trace pipeline)** — the §7.4 lines coloring
  reuses `setBandTraceFishSet` as the click target for track rows;
  band-trace is the single-cohort interactive, this spec is the
  multi-cohort enumerative
- **Turn 165 (G-panel auto tab)** — review surfaces extend with
  band-track classification once both ship

---

## Files in this handoff

- `specs_todo/SPEC_band_track_extraction_and_l3_single_band_rows.md` —
  the fully-extended spec (1,769 lines)
- This handoff note

No atlas changes this turn. No new tests. The spec is the
deliverable.

End of handoff.
