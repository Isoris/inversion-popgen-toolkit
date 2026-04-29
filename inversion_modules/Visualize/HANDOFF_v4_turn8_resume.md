# HANDOFF — resume in next chat

This is the master handoff for the v4 scrubber + 226-sample C. gariepinus
inversion-detection work as of end of this chat session.

The current chat made significant progress on UX polish, then hit a long
biological-design discussion that uncovered three big asks the next chat
needs to implement.

---

## Latest staged file

`/mnt/user-data/outputs/pca_scrubber_v3.html` — **1,517,307 bytes**

Last shipped at v4 turn 7 (the prior turn before this one). All test suites
green at staging:

```
test_v399.js                              754/754
test_v4_diversity_atlas_runtime.js         51/51
test_population_atlas_v1_runtime.js        75/75
test_v41_metrics_runtime.js                68/68
test_v41_recluster_rebuild_runtime.js      38/38
test_v41_uv_modes_runtime.js               38/38
test_v399_theta_pi_scrubber_runtime.js     32/32
test_v4_pca_lasso_runtime.js               20/20
test_v399_winsum_color_runtime.js          33/33
test_v399_blocks_runtime.js                62/62
test_v399_candidate_bar_runtime.js         25/25
test_v399_lines_color_mode_runtime.js      12/12
test_v399_chips_inspector_runtime.js       34/34
test_v399_karyotype_tier_runtime.js       140/140
test_v399_catalogue_export_runtime.js      32/32
test_v399_registry_interop_runtime.js      31/31
                                       --------
                                       1,445 / 1,445 ✅
```

(The v4 turn 7 handoff and earlier copies of this handoff said "1,485"
— that was an arithmetic error I propagated forward. Re-summed this turn
when bundling: actual total is 1,445. No tests were lost; only the
displayed total was wrong.)

`/tmp/w2/pca_scrubber_v3.html` matches the staged file (verified, same byte count).

---

## What shipped this session (v4 turns 1–7)

A long UX-polish arc. Brief recap:

**Turns 1–3 (UX polish round 1):**
- Tracked-samples aside redesigned compactly with K-cycle button + 6 K-colored
  band buttons + lasso checkbox + dark Clear button
- PC1/PC2 axis labels show λ_k variance %
- L3 toolbar hover styling on inline-styled buttons
- L3 column overlap fix (`.ct-content { flex: 0 0 auto }`)
- theta pi → θπ rename
- Stats-to-right-of-table layout in L3 contingency
- L3 toolbar wrap fix
- |Z| label below L2 in windows-page strip
- L1/L2 fragment dark-red separators
- PCA lasso wired (Shift-drag = manual group; checkbox = tracked-samples)

**Turn 6:**
- Cross-pane sample spotlight: click any L3 mini-PCA dot to label that sample
  in all panes + light up its contingency cell. Show tracked button +
  Clear spotlight button complete the loop.

**Turn 7:**
- 1w/5w/10w/Nw button styling matched to L2 styling (consistent min-width)
- Promote-to-candidate green flash on both promote buttons
- Reconciled with parallel work from a separate chat session: kept the
  Save/Load JSON sessions and ★ promote button additions, layered turn-7
  changes on top

Three handoffs already in `/mnt/user-data/outputs/`:
- `HANDOFF_v4_population_atlas_v1.md` (sibling HTML scaffold work)
- `HANDOFF_v4_turn4.md`
- `HANDOFF_v4_turn7.md`

Read those for the full UX-polish arc context.

---

## What's queued from prior turns (NOT YET SHIPPED)

Items that were proposed or partially designed in this chat but did NOT
make it into a staged build:

### 1. Ancestry-confounding test (designed last turn, no code shipped)

**What it is:** A new panel on the Candidate Focus page (page 3) that
computes Cramér's V + max |Pearson correlation| between K-band
assignments and genome-wide Q ancestry. Three-tier verdict:
`BAND_TRACKS_ANCESTRY` (suspect candidate, V≥0.6 or max|corr|≥0.7),
`BAND_INDEPENDENT_OF_ANCESTRY` (real local signal, V<0.3 and max|corr|<0.4),
`INCONCLUSIVE` (otherwise). Empty state when q_proportions layer absent.

**Status:** complete code design written in chat, NOT applied to file
because tools were intermittently down at the time. The full patch is in
the prior turn's response in this chat — three blocks of code (helpers,
HTML scaffold, renderer + wiring). Should drop in cleanly.

**Next step:** apply patch + verify q_proportions layer detection finds
the right field name in the precomp JSON. May require schema confirmation
from MODULE_2B's ngsadmix output spec.

**Estimated effort to finish:** ~30 min (it's already written, just needs
applying + adjusting the layer detector cascade + writing the test).

### 2. Move chrom axis ticks below per-sample lines panel

Originally Quentin's request: "The axis text and ticks for the genomic
coordinates we add them below the Per sample lines instead of below the
robust Z so that the per sample lines and the tracked samples PCA are at
the same level in pixels."

Invasive layout change. Not started. ~30–45 min.

### 3. Drag handle JS for `#pcaAsideResize`

HTML+CSS already in place from a much earlier turn; pointer-handler JS
never wired. ~10 min.

### 4. Free-mode polish

Three sub-items: tab raising, |Z| drag handle, merge button placement.
Each ~20–30 min. Low priority.

---

## What's NEW this turn (the BIG asks)

This chat ended on a long biological-design discussion that produced three
substantial requirements. Schema entries §24/§25/§26 are in
`/mnt/user-data/outputs/SCHEMA_v2_25_DRAFT.md`. The implementation spec is in
`/mnt/user-data/outputs/SPEC_v4_window_row_overlapping_regimes.md`.

**Summary:**

### A. Per-window row in candidate strip (smallest, highest priority)

L2 envelopes are too coarse for some inversions. User needs to be able to
draft a candidate at window resolution. Adds a third zone-bar row to the
candidate strip, visible only in candidate mode. Lets user click + drag
at window precision. Existing candidate `start_w`/`end_w` already supports
window ranges — this is rendering + click-handling, not data model.

State change: `state.l3Draft` extended with `resolution: 'L2' | 'W'`,
`start_w`, `end_w`, `l3_cuts`. All optional; old shape continues to work.

**~250 lines, single turn.** Highest-value, smallest-scope deliverable.

### B. Overlapping candidates (medium)

Currently the candidate-bar renders all candidates in one flat row. With
overlap, they z-stack poorly. Refactor to lane-assign overlapping candidates
into vertical lanes inside the bar. Click-handling becomes 2D. Catalogue
gets a "regimes" column (empty until C ships).

No data-model changes — just rendering. **~400 lines, single turn.**

### C. Regime registry (biggest)

Many-to-many candidate ↔ regime mapping. New `state.regimeList[]` per chrom,
new persistence key, full CRUD UI on candidate page. Each regime has an
`axis_topology` enum (one_axis_3band, two_axes_independent, etc. — see
SCHEMA §26 for the full vocabulary).

**Critical design decision:** `axis_topology` is user-chosen, NOT auto-
computed by the scrubber. The architecture-call thresholds aren't validated
yet, and committing to defaults would be misleading. The scrubber records
the user's interpretation. The R-side `STEP_R45_arrangement_tree.R` (planned)
will eventually populate this field automatically.

**~700 lines, likely 2 turns.**

### D. Scale-stability test (toggleable L3 mode) — NEW THIS TURN

A toggle on the L3 toolbar that swaps the existing L3 contingency-table
panel for a scale-stability panel. Same focal region, K-means recomputed
at three zoom levels (single window / slab / full candidate), rectangular
contingency tables show how bands at scale-A map to bands at scale-B.

Verdict cascade: `STABLE_3BAND` (classic single inversion) /
`STABLE_6BAND` (stable multiband, needs band-tree to decide single vs
compound) / `NESTED_3IN6` (3-band inversion with substructure at fine
scale) / `OVERLAP_BREAKS_3` (possible two overlapping inversions in the
middle of the candidate) / `UNSTABLE` (artifact-suspect).

Uses existing `getSlabClusterAt(s, e, K)` — no new clustering primitive
needed. Strided sampling for large candidates keeps compute under 200 ms.

**Critical caveat:** verdict thresholds (ARI ≥ 0.85, mapping rate ≥ 0.80)
are heuristic defaults, NOT validated against ground truth. Documented as
adjustable, NOT publication-ready calls.

**~600 lines, likely 2 turns.** Independent of B/C — could ship right after A.

### Order

**A → D → B → C**, with real-data testing checkpoints between each.
Why insert D before B/C: D is the test that distinguishes "real 6-band"
from "3-band-with-substructure," which directly informs whether the
many-to-many regime mapping (C) is even needed. If most candidates come
back NESTED_3IN6 / STABLE_3BAND, C can be scoped down. If STABLE_6BAND
or OVERLAP_BREAKS_3 dominate, C is justified.

If only one deliverable can ship before the next data freeze: A.
If two: A + D.

---

## Critical context the next chat needs

### Cohort isolation (NEVER conflate)

- 226-sample pure C. gariepinus hatchery cohort on LANTA → "MS_Inversions_
  North_african_catfish" → the inversion scrubber + this work
- F1 hybrid (C. gariepinus × C. macrocephalus) → genome assembly paper, separate
- Pure C. macrocephalus wild cohort → future paper

### Working style

- Compact non-fragmented pipelines
- Original tested code used VERBATIM, never reimplemented
- Single central config per module (`00_inversion_config.sh` style)
- Targeted fixes over rewrites
- Deliverables as files in `/mnt/user-data/outputs/`
- Multi-session parallel chats — ALWAYS reconcile script versions before
  assuming state. This session reconciled twice with parallel work
  (v4 turn 5 added Save/Load JSON + ★ promote button, kept those changes).
- Concise/abbreviated English from Quentin
- Manuscript-deliverable angle: minimize headache, dump results in 3 atlases,
  send to teacher, finish.

### Scrubber architecture invariants

- 3-atlas plan: pca_scrubber_v3.html (Inversions) + population_atlas_v1.html
  (Population, sibling) + variants_atlas_v1.html (Variants, planned)
- Schema currently at v2.24, going to v2.25 with this turn's additions
- `state.candidate` = active focal candidate; `state.candidateList[]` =
  saved registry per chrom; `state.regimeList[]` = NEW from §25, default []
- L1/L2 are algorithmic envelope levels (from similarity-matrix clustering).
  L3 is NOT a third hierarchy — see §24. L3 cuts are per-candidate UI
  subdivisions only.

### Don'ts (the ones that matter)

- Do NOT mass-rename step files
- Do NOT hardcode thresholds; configurable in central config
- Do NOT auto-assign `axis_topology` in the scrubber — keep user-editable
- Do NOT introduce L3 as a hierarchical analysis level; L3 cuts are
  per-candidate UI only
- Do NOT block A on B/C — they're independent features
- Do NOT migrate old candidates eagerly — read-side migration only

---

## Recommended next-chat plan

**Turn 1:** Apply the queued ancestry-confound deliverable (~30 min). Run
on LG28's existing K=6 candidates. The verdict distribution is informative
on its own — if most candidates come back BAND_TRACKS_ANCESTRY, the K=6
calls are likely artifact. If most are INDEPENDENT, they're likely real.

**Turn 2:** Implement Deliverable A (per-window row). Real-data test on
3–5 candidates where the L2 is much larger than the inversion peak.

**Turn 3–4:** Implement Deliverable D (scale-stability test). Run on
LG28's full candidate catalogue. Verdict distribution decides whether
C is justified at full scope or scoped down.

**Turn 5:** Implement Deliverable B (overlapping candidates, lane stacking).

**Turn 6–7:** Implement Deliverable C (regime registry).

**Turn 8+:** Queue items 2–4 (chrom axis ticks, drag handle JS, free-mode
polish) when there's appetite for return to UX polish.

---

## Files in /mnt/user-data/outputs/

Already present from prior session:
- `pca_scrubber_v3.html` (1,517,307 bytes — v4 turn 7 staged)
- `population_atlas_v1.html` (54,532 bytes — sibling scaffold)
- All test suite files (`test_*.js`)
- `HANDOFF_v4_population_atlas_v1.md`, `HANDOFF_v4_turn4.md`, `HANDOFF_v4_turn7.md`
- `SCHEMA_V2.md` (354 KB, schema through v2.24)
- Plus all earlier v3.95–v3.99 handoffs

NEW from this turn (no scrubber file change):
- `SCHEMA_v2_25_DRAFT.md` — schema additions §24 (window-row + L3 cuts),
  §25 (overlapping candidates + regime mapping), §26 (architecture-call
  vocabulary), §27 (scale-stability verdict vocabulary)
- `SPEC_v4_window_row_overlapping_regimes.md` — implementation spec for the
  four new deliverables A/B/C/D (D = scale-stability test, added this turn)
- `HANDOFF_next_chat_resume.md` — this document

Apply when next chat starts:
1. Read this handoff
2. Reconcile `/tmp/w2/pca_scrubber_v3.html` against the staged file (should
   match at 1,517,307 bytes; if larger, there's parallel-chat work to
   integrate)
3. Run all test suites to confirm 1,445/1,445 baseline
4. Pick the next deliverable and proceed
