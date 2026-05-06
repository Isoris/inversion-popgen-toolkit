# SPEC — Fully automatic genome-wide inversion identification (umbrella vision)

**Status**: drafted turn 130 final session. **Umbrella spec** that ties
together every piece below into the "click a button, get all inversions
on all chromosomes" pipeline Quentin has been asking for since the
beginning.

**Trigger** (Quentin, multiple chats):
> *"Soon its all ready then upload all chromosomes and it will be all
> automatic."*
>
> *"You think that when we will activate the select candidate button in
> detailed per window mode, it will use all 3 pages CUSUMS + contingency
> tables + dosage. And identify every inversion for every chromosome all
> automatic when we load all JSONs at once. Bc I'm lazy person."*

This is the master vision. Every other spec in this folder serves it.

---

## 1. The end-state user experience

**Day 1.** Quentin uploads all 28 chromosome JSONs (or starts a chrom-by-chrom
upload). Each JSON contains: precomp, L1/L2 envelopes, K-means labels,
GHSL panel, dosage chunks, θπ tracks.

**Day 1, 5 minutes later.** The atlas has automatically computed:

- L2-sweep inheritance per chromosome (auto-promoted candidates)
- Fish-trajectory lineages per chromosome
- Cross-candidate inheritance Jaccard (where ≥2 candidates exist)
- Sliding-window inheritance at 5w/10w resolution
- Per-candidate karyotype assignments (REF/HET/INV)
- Per-candidate boundary refinement (CUSUM consensus)
- Per-candidate validation gates (silhouette, agreement, sample sizes)

**Day 1, 10 minutes later.** Quentin reviews:

- A single "auto-promoted" tab in the G-panel showing every chromosome's
  candidates, ranked by confidence
- Dashed-outline rows in the candidate strip + page-2 list (already
  shipped Slice 0)
- A genome-wide ideogram showing all auto-promoted candidates colored
  by tier
- Per-candidate breeding-readiness cards (Atlas 5)

**Day 1, 30 minutes later.** Quentin clicks "export manuscript bundle."
Atlas emits a Markdown file with: per-candidate tables, boilerplate
text, citation keys, figure references, supplementary atlas links.
Quentin polishes the prose by hand.

## 2. What's already built

| Component | Status | Spec / location |
|---|---|---|
| Local PCA precomp + L1/L2 staircase | shipped (LANTA-side) | `STEP_D17_multipass_*` |
| GHSL panel + K-means | shipped (LANTA-side) | `STEP_C04_*` |
| θπ tracks | shipped (LANTA-side) | `STEP_TR_*` |
| Atlas viewer (per chrom) | shipped | `Inversion_atlas.html` |
| Candidate registry + bridge | shipped turn 129 | atlas-side |
| Inheritance Jaccard (across confirmed candidates) | shipped turn 117 | atlas-side |
| Fish-trajectory lineage compute | shipped turn 130 Slice 1 | `SPEC_distant_band_concordance_fish_trajectory.md` |
| Lineage UI (color mode + strip) | shipped turn 130 Slice 2 | same spec |
| Auto-review dashed-row infrastructure | shipped turn 130 follow-up | `SPEC_review_surfaces_auto_and_lineages.md` Slice 0 |

## 3. What's NOT yet built — the gap to "fully automatic"

All ordered by dependency. Each is a spec in this folder.

### Tier 1 — producers (the critical path)

1. **L2-sweep auto-promote** — `SPEC_l2_sweep_inheritance.md` Slice 1.
   Runs the same Jaccard pipeline on every L2 envelope (not just
   user-promoted), auto-promotes candidates when 5 gates pass. **The
   single biggest unlock**: removes the "user must manually promote 5+
   candidates" cold-start gate.

2. **Sliding-window inheritance** — `SPEC_sliding_window_inheritance.md`.
   Same Jaccard, but on fixed-width tile windows (5w/10w/Nw). Catches
   inheritance signal between L2 boundaries.

3. **Cross-chromosome lineage aggregator** —
   `SPEC_cross_chromosome_lineages.md` (queued in this batch).
   Match lineages across chromosomes via shared fish memberships.
   Tells the user "fish 12 and 47 are in lineage A on LG12 and lineage
   A on LG28" — they share inheritance backgrounds.

4. **Multi-chrom JSON load orchestrator** —
   `SPEC_multichrom_load_orchestrator.md` (queued).
   Loading 28 JSONs in one go. Memory budget, lazy load, per-chrom
   compute scheduling.

### Tier 2 — review surfaces

5. **G-panel scaffold** — `SPEC_g_panel_unified_groups.md` Slice 1.
   Popup with karyotype / inheritance / manual / auto / lineages tabs.
   Single review surface for everything algorithm-proposed.

6. **G-panel `auto` tab** — `SPEC_review_surfaces_auto_and_lineages.md`
   Slice 2. Bulk confirm / dismiss / inspect for auto-promoted
   candidates.

7. **G-panel `lineages` tab** — same spec, Slice 3. Track / pin / dismiss
   for fish-trajectory lineages.

8. **Genome-wide ideogram** —
   `SPEC_genome_wide_ideogram.md` (queued). All 28 chromosomes,
   confirmed + auto candidates colored by tier. Click → jump to
   candidate.

### Tier 3 — manuscript export

9. **Manuscript bundle export** —
   `SPEC_manuscript_bundle_export.md` (queued). Markdown + TSV + JSON
   from the candidate catalogue. Boilerplate per-candidate prose +
   tables. Already partially supported by existing TSV/MD export on
   page 3.

10. **Per-candidate breeding-readiness cards** — Atlas 5, already
    spec'd in `c03fc41e` plan-arrangement chat. Drives marker-design +
    pairing recommendations.

### Tier 4 — diagnostic / decision-support

11. **Inheritance diagnostic protocol** —
    `SPEC_inheritance_diagnostic_protocol.md` (already drafted).
    What to look at on real LG28 data to calibrate thresholds.

12. **Inheritance unification decision** —
    `SPEC_inheritance_unification.md` (already drafted, decision-gated).
    Whether to collapse all four inheritance modes (candidate-only /
    L2-sweep / sliding / fish-trajectory) into one dispatcher.

## 4. Estimated total work

Tier 1 (producers): ~3-4 turns total.
Tier 2 (review surfaces): ~2-3 turns.
Tier 3 (manuscript export): ~1-2 turns.
Tier 4 (diagnostic): no code, runs against real data.

**Total: ~6-9 turns to ship the fully-automatic vision.**

## 5. Honest acknowledgments

The phrase "click a button → manuscript text" is misleading. What the
pipeline can realistically produce is **boilerplate** per candidate:

> *"On chromosome LG28, candidate at 15.115–17.985 Mb shows 60 REF, 106
> HET, 60 INV samples (Cochran-Armitage trend p < 1e-8). Inheritance
> Cramér's V with adjacent LG28 candidate at 22.4 Mb is 0.71, suggesting
> shared chromosome ancestry. Fish-trajectory lineage analysis recovers
> 4 lineages (n=47, 38, 76, 65)."*

The interpretive paragraphs — what these inversions mean for the
hatchery, what selection might be doing, what recurrent breakpoints
suggest — those are **always Quentin's job**. No pipeline writes those.

The role of the atlas is: produce numbers Quentin trusts, present them
with enough evidence to defend at peer review, and emit them in a form
that drops cleanly into the manuscript skeleton.

## 6. What this spec is NOT

- **Not a replacement for human judgment** on candidate confirmation.
  The auto-promote produces proposals; the user confirms.
- **Not a publication automaton.** It writes per-candidate facts,
  not the discussion section.
- **Not a black box.** Every auto-promote, every lineage call, every
  boundary refinement is inspectable in the atlas via existing surfaces.

## 7. Cross-references

- `SPEC_distant_band_concordance_fish_trajectory.md` — Slice 1 + 2
  shipped turn 130.
- `SPEC_l2_sweep_inheritance.md` — auto-promote producer.
- `SPEC_sliding_window_inheritance.md` — fine-resolution complement.
- `SPEC_review_surfaces_auto_and_lineages.md` — G-panel + dashed UI.
- `SPEC_g_panel_unified_groups.md` — popup scaffold.
- `SPEC_inheritance_diagnostic_protocol.md` — calibration protocol.
- `SPEC_inheritance_unification.md` — decision-gated refactor.
- New specs queued in this batch:
  - `SPEC_cross_chromosome_lineages.md`
  - `SPEC_multichrom_load_orchestrator.md`
  - `SPEC_genome_wide_ideogram.md`
  - `SPEC_manuscript_bundle_export.md`
  - `SPEC_recombinant_dosage_changepoint_detector.md`
  - `SPEC_per_candidate_breeding_readiness_card.md`
  - `SPEC_marker_panel_design_atlas.md`
  - `SPEC_boundary_consensus_aggregator.md`
  - `SPEC_metric_overlay_priors.md`
  - `SPEC_karyotype_per_interval_intersection.md`
