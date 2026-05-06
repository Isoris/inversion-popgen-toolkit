# HANDOFF — turn 166 — band-track + het-skeleton spec extension

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (76,407 LOC; **unchanged** this turn)
**Spec target**: `specs_todo/SPEC_band_track_extraction_and_l3_single_band_rows.md` (1612 lines, +640 LOC from the turn 165 close baseline)
**Working dir**: `/home/claude/work/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**What this turn shipped**: spec writing only — no atlas code changes,
no test runs needed. The spec from the prior session was extended with
het-anchored mode (the het band as the principal track + lines-panel
coloring + skeleton overlay) so the next implementer has a single
unified spec covering everything from per-band rows in the L3 K×K
view to the full het-skeleton-driven candidate skeleton.

**No tarball regression risk.** Atlas binary is byte-identical to the
turn 165 close. The 3652/0 sweep result still holds.

---

## 0. Why this is a spec turn, not a code turn

Two earlier interactions in this session established the design:

1. The first message was the original "band track extraction" spec
   request (§3 onward of the existing spec — Jaccard + retention,
   path construction, A-B-A double-crossover detection, four TSV
   outputs).
2. A clarifying message: *"its like the per sample lines but with
   lines colored by their 'intervals'"* — meaning band tracks
   should color the per-sample-lines panel, making each fish's
   line a stable color tied to its track membership.
3. A second clarifying message: *"add this to our single band. Its
   basically we track the het band from dosage bc its the skeleton
   for our candidate inversion interval"* + the long het-carrier
   continuity proposal. This made the het band the **principal
   track** of the spec, not just one band among K.

Quentin's final instruction: *"Its the spec to extend our contingency
and so on we have mode het and so on. Idkkk but try your best. Then
bundle and handoff for next chat."*

So this turn integrated the het-anchored mode into the existing
spec rather than writing a separate doc. Single spec, one place to
read, one set of cross-references.

---

## 1. What the spec now covers

**Symmetric mode** (the original spec, unchanged):
- Per-band rows under each K×K contingency table in the L3 panel
- Band-track extraction by greedy path construction over connection edges (§4 connection rule: Jaccard ≥ 0.65 OR retention pair ≥ 0.75)
- A-B-A / split-track detection across gaps (§5.3) with the cross-segment Jaccard guard against accidental merges
- System-level classification: simple_single_inversion / nested / split / fragmented / complex
- Four TSVs (band_tracks, band_track_edges, system_band_summary, sample_band_paths)
- Three UI surfaces: per-band rows in K×K (slice 1), track strip + classification chip (slice 2/3)

**Het-anchored mode** (new this turn — the extension):
- §1.4 frames the het band as the **principal track** — biologically privileged, the inversion's skeleton/backbone
- §1.4 establishes the dropout-vs-boundary distinction with explicit wording discipline ("het signal stops here," not "inversion ends here")
- §1.5 frames the per-sample-lines panel as the band-track display via the existing `state.linesColorMode` slot
- §4A.1 het-band identifier per pane: het-enrichment score with three signal-source fallbacks (mean het rate via `_computeHetRateForL2` from turn 129 → mean dosage → middle-PC1 heuristic)
- §4A.1 supports zero / one / multiple het bands per pane (no forced uniqueness)
- §4A.2 carrier-continuity edges: extends the §4 connection predicate with overlap coefficient (`O = |A∩B|/min(|A|,|B|)`) so subset/dropout cases pass when Jaccard alone would fail
- §4A.3 skeleton extraction with **internal-dropout tolerance**: forward walk can skip up to MAX_DROPOUT_GAP panes when downstream het bands re-establish high carrier overlap (Quentin's exact A-B-C example: don't split when carriers vanish briefly and resume the same)
- §4A.4 het-skeleton output shape: principal flag, dropout panes list, dropout diagnosis, consensus carriers, carrier stability score, classification (`clean_skeleton` / `skeleton_with_dropout` / `reappearing_skeleton` / `mosaic` / `fragmented`)
- §4A.5 explains why this lives inside the symmetric-mode spec rather than a separate doc (it's a specialisation, reusing every primitive from §3-§5 with one extra metric and one extra walk-step)
- §6 adds TSVs E (`het_skeletons`) and F (`het_carrier_edges`) to the existing four
- §7.4 adds four `band_track_*` sub-modes for `state.linesColorMode`: all / het / principal_only / with_het_overlay; the resolver fills in the v3.99 turn 14e+ stub
- §7.5 adds three skeleton overlays (candidate strip warm-tint, lines panel warm-tint band, L3-strip pane-header dots)
- §8 slices updated: slice 1 includes the het-band identifier; slice 2 grows to ~1.5 turns covering symmetric path construction + het skeleton + lines coloring + skeleton overlay; slice 3 grows from four to six TSVs
- §10 tests extended with het-anchored coverage (≈12 new test items spanning identifier, skeleton extraction, dropout-skip, lines-panel resolver)

**Vocabulary discipline (§2)** — the code and TSV outputs follow:

- "het signal stops here" / "carrier-set continuity breaks here" — operational, used in code/output
- "the inversion ends here" — biological, NOT used by the algorithm; reserved for Quentin's manuscript text
- "this is a double crossover" — NOT used; outputs use `class: split_track_reappearing` + `possible_interpretation: double_crossover_like` as separate fields
- "track A is the same as track B" — NOT used; algorithm reports a candidate, user decides

---

## 2. Files touched this turn

```
specs_todo/SPEC_band_track_extraction_and_l3_single_band_rows.md   +640 LOC
  - Status header updated (estimate 2.5–3.5 turns; "Extended this session" note)
  - Trigger expanded with the two clarifying messages
  - §1.4 het band as principal track (NEW, 56 lines)
  - §1.5 lines panel as band-track view (NEW, 44 lines)
  - §2 vocabulary extended (+8 new terms: het band, het skeleton, het carrier, carrier set, dropout, "detectable divergent haplotype signal", "the inversion ends here" forbidden)
  - §4A het-anchored mode (NEW, 5 subsections, 196 lines)
    - §4A.1 het-band identifier with 3 signal-source fallbacks
    - §4A.2 carrier-continuity edges with overlap coefficient
    - §4A.3 skeleton extraction with internal-dropout tolerance
    - §4A.4 distinguished outputs (HetSkeleton shape, classifications)
    - §4A.5 why this lives in the symmetric-mode spec
  - §6.1 in-session state shape extended with het_skeleton fields
  - §6.2 TSVs E and F added (het_skeletons, het_carrier_edges)
  - §6.3 filename list expanded to six files
  - §7.3 mode chip note added
  - §7.4 lines-panel `band_track_*` color sub-modes (NEW)
  - §7.5 het skeleton overlay (NEW)
  - §7.6 (former §7.4 "What this is NOT") renumbered
  - §8 slices: slice 1 +het-band identifier item; slice 2 +skeleton extraction +lines resolver +overlay; slice 3 +het TSVs +mode chip
  - §10 tests: ≈12 new test items in slice 1 + slice 2 + slice 3 sections
  - §11 cross-references: added _computeHetRateForL2, _HET_RAMP, state.linesColorMode

specs_todo/README.md                                              (unchanged this turn — already added in prior turn)

Inversion_atlas.html                                              (unchanged this turn)
tests/                                                             (unchanged this turn)
```

No code changes. No test runs needed. The atlas binary at the end of
this turn is byte-identical to the close of turn 165.

---

## 3. Bundle contents

```
Atlas_turn166_2026-05-05_tar.gz
  Atlas/
    Inversion_atlas.html                           76,407 LOC (= turn 165)
    specs_todo/
      SPEC_band_track_extraction_and_l3_single_band_rows.md   1612 LOC (NEW for this turn — extended)
      README.md                                    (already mentions this spec from prior turn)
      [all other specs unchanged]
    tests/                                         (unchanged from turn 165)
    handoffs/
      current/HANDOFF_2026-05-05_turn165_g_panel_auto_tab.md   (turn 165 work)
      current/HANDOFF_2026-05-05_turn166_band_track_spec.md    (THIS handoff)
    HANDOFF_2026-05-05_turn166_band_track_spec.md  (also at root)
```

---

## 4. What the next chat should do

**Read order for whoever picks this up**:

1. This handoff (you're reading it)
2. `SPEC_band_track_extraction_and_l3_single_band_rows.md` end-to-end
   — it's long (~1600 lines) but covers everything; §1.4-1.5 +
   §4A are the new substance
3. `SPEC_distant_band_concordance_fish_trajectory.md` (turn 130, partly shipped)
   — the dual representation, useful background
4. `_computeHetRateForL2` definition in `Inversion_atlas.html` line 16735
   — this is the primitive §4A.1 builds on
5. `renderL3Panel()` in `Inversion_atlas.html` line 48685 — this is
   where slice 1's single-band rows render

**Suggested next turn (turn 167)**: implement Slice 1.

The deliverables (per spec §8 slice 1):
- `_buildBandEdgesForPanes(panes, K)` — pure compute over the L3
  strip's panes; returns the edge bundle. Use `getL2Cluster(l2idx).fixedKLabels`
  per pane.
- `_classifyHetBandsForPane(l2idx, K, opts)` — per-pane het-band
  identifier with the three signal-source fallback chain.
- `_ssBandRowHtml`, `_ssBandRowsForPaneHtml` — renderers in the
  style of `_ssContingencyTableHtml`.
- `state.l3SingleBandRows` toggle + checkbox in the L3 head.
- `renderL3Panel()` call site update.
- Cache invalidation on candidate / K / mode change.
- Test: `test_turn167_l3_single_band_rows_and_het_identifier.js`
  with the sections enumerated in §10 slice 1 tests (~80 assertions
  expected based on prior-turn density).

**After Slice 1, suggested order**:

- Slice 2 (turn 168, probably 1.5 turns split into 168A + 168B):
  - 168A: symmetric path construction + A-B-A detection + system
    classification + track-strip canvas
  - 168B: het-anchored mode (skeleton extraction with dropout-skip)
    + lines-panel `band_track_*` resolver + skeleton overlays
- Slice 3 (turn 169): TSV exports (six serializers) + classification
  chip + mode chip in the L3 toolbar

**Why slice 1 first**: the per-band rows under K×K tables are the
visible diagnostic that explains every downstream decision. Once
the user can see "g0 (n=56) → c0 with 56 fish, persists; g2 (n=30)
→ {c2: 12, c3: 18}, splits" with het-class tints, the rest of the
work is predictable. Slice 1 is also self-contained (no skeleton
extraction, no lines-panel changes) so it can ship cleanly.

---

## 5. Open design questions still in the spec

These are flagged in §9 and should be revisited when slice 2
starts:

1. **K=3 vs K=6 default for compute.** Recommendation in spec:
   compute at `state.k` (whatever the user has set); document in
   the chip that K=6 produces denser tracks.
2. **Pane source.** L3 panes = `state.candidate.l2_indices` for
   now; free-form pane lists deferred.
3. **Cross-system band tracks.** Defer to fish-trajectory spec.
4. **Auto-label proposals.** Out of scope; future spec hook.
5. **Persistence in session save/load.** Recommendation in spec:
   no — derived state, recompute on load.
6. **Render cost.** Spec estimates negligible (<300 rects on
   track strip, <200 row fragments per L3 strip). Confirm in
   slice 2 implementation if Quentin reports lag on real LG28.

New questions raised by the het-anchored extension:

7. **HET_JACCARD_MIN default = 0.5 vs §4's JACCARD_MIN = 0.65.**
   The het mode is intentionally looser because dropout/marker
   coverage causes legitimate carrier-set shrinkage between
   panes. Confirm with Quentin on real LG28; consider exposing
   both as user-tunable in the L3 toolbar if the defaults need
   per-region tuning.
8. **MAX_DROPOUT_GAP = 3 panes.** Empirical default; should be
   tuned on real LG28 data. The skeleton class output makes this
   visible (a user-observed false-bridge would show up as an
   incorrectly-merged skeleton with low cross-gap Jaccard,
   easy to flag).
9. **Lines-panel sub-mode default.** When `state.linesColorMode
   = 'band_track'` is selected for the first time, which sub-mode
   should be active? Recommendation: `band_track_with_het_overlay`
   — every track gets a unique color but the het skeleton pops
   visually. Confirm in slice 2 review.

---

## 6. Honest risks for the next session

**Reading-comprehension risk for the next implementer.** The spec
is now 1612 lines. §1.4-1.5 + §4A are dense. Whoever picks this
up should budget ≈30 minutes for a thorough read before
implementing — especially the §4A.3 skeleton-extraction algorithm
with dropout-skip, which has subtle ordering (forward walk +
lookahead + dropout-classification post-step).

**Implementation effort estimate could be off.** Slice 2 is
estimated at ~1.5 turns. If the het-anchored mode's skeleton
extraction interacts badly with the symmetric path construction
(e.g. a synthetic A-B-A fixture gets bridged by the dropout-skip
when it shouldn't, requiring the §4A.3 classification step to
reach into the §5.3 cross-segment-Jaccard guard), slice 2 might
need to split as 168A / 168B as suggested above. Plan for that.

**Het signal availability.** The het-band identifier requires
`dosage_chunks` to be loaded for the primary `het_rate` source.
On chromosomes where this layer is missing (or hasn't been
emitted by the HPC pipeline), the algorithm falls back to
PC1-middle heuristic and emits a warning. The slice 2 implementer
should make sure the fallback path doesn't silently produce
poor-quality skeletons — the warning must surface visibly in the
mode chip's tooltip and the het-skeleton TSV's `het_signal_source`
column.

**`state.linesColorMode` resolver scaffolding.** The spec assumes
the v3.99 turn 14e+ stub still exists with the right signature.
Slice 2 implementer should verify before coding by grepping
`_resolveSampleColorByMode` / `_resolveSampleScopeColor` in the
atlas; the spec cited those names but I didn't re-check live this
turn (no atlas changes shipped). If the function names have
drifted, slice 2 needs a small adjustment to the resolver
integration plan.

**No code regressions because no code changed.** This turn was
spec-only, so the 3652/0 turn-test sweep from turn 165 close still
holds. The atlas tarball is structurally identical to turn 165
except for the spec file diff.

---

## 7. What's NEXT after the band-track stack lands

(Reproduced from turn 165 §6 with this spec inserted in priority order):

1. **Run the auto-promote tab on real LG28** (turn 165 deliverable
   under Quentin's review)
2. **Implement Slice 1 of the band-track spec** (turn 167) —
   per-band rows under K×K + het-band identifier
3. **Implement Slice 2 of the band-track spec** (turn 168A/B) —
   symmetric path construction + het skeleton + lines coloring
4. **Implement Slice 3 of the band-track spec** (turn 169) —
   six TSVs + classification chip + mode chip
5. **Run band-track on real LG28** — calibrate thresholds
6. **`_inhTooltip` DPR-scaling fix** — latent bug from turn 162
7. **Slice 3 of `SPEC_review_surfaces_auto_and_lineages.md`**
   (lineages tab) — gated on Quentin running lineage compute on
   real data
8. **Run brackets / combinatorial enumeration** (band-trace turn
   161 §6) — may become redundant once the band-track spec ships
   since band tracks generalise the band-trace concept

---

End of handoff.
