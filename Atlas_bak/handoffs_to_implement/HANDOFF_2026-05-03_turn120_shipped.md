# HANDOFF — 2026-05-03 atlas session, turns 117 → 120

**Author**: prior Claude (this session).
**Recipient**: future Claude (next session) and Quentin.
**Atlas line count**: **59,033** (started session at 57,556 partial → +1,477 over 5 turns).
**Tests**: **218/218** green across 5 suites. No regression at any point in the session.

---

## TL;DR — what was shipped in this session

| Turn | Description | Atlas Δ | Tests | Status |
|---|---|---|---|---|
| 117 | Focal-vs-bg widget — full integration on pages 11 + 16 | +223 | 15 | ✅ |
| 117b | Dot plot pure-hover model (drop click-to-pin) | 0 atlas / -30 module | 14 | ✅ |
| 118 | Manuscript bundle Tier-1 enrichment (6 new blocks + 23 TSV cols) | +552 | 82 | ✅ |
| 119 | Age & origin atlas panel (cheat30 GDS surface) | +512 | 74 | ✅ |
| 120 | PCA scree-plot inset on tracked-sample plane | +190 | 33 | ✅ |
| **Total** | | **+1,477** | **218** | |

Three audit/spec deliverables also shipped (no atlas changes):

- `AUDIT_confirmed_candidate_export_gaps.md` (528 lines) — the audit that drove turn 118
- `SPEC_inversion_age_atlas_surface.md` (559 lines) — Task A (between-inversion ranking) + Task B (within-group ranking) using electric-catfish μ_year anchor
- `SPEC_age_origin_panel.md` (430 lines) — turn 119 spec
- `SPEC_pca_scree_inset.md` (217 lines) — turn 120 spec
- `PATCH_export_precomp_lam_top_k.md` (145 lines) — R-side patch for turn 120

---

## What the next session should pick up

Listed in priority order based on what's spec'd, what's blocking, and
what's most manuscript-relevant.

### Priority 1 — Tier-2 manuscript bundle enrichment (~150 atlas lines, 1 turn)

The audit (`AUDIT_confirmed_candidate_export_gaps.md`) identified seven
Tier-2 gaps after the Tier-1 ship in turn 118. Architecture is cleanly
extensible — each new enrichment helper plugs into
`_bundleEnrichmentForCandidate` dispatcher with no other changes:

1. **Repeat density at boundaries** — read `state.repeatDensity[chrom]`,
   filter windows around `c.boundary_left.zone_*` and
   `c.boundary_right.zone_*`. Add to bundle as
   `_bundleRepeatDensityBlock(c)`.
2. **ncRNA density at boundaries** — same shape as repeat density,
   reading `state.ncRNADensity[chrom]`.
3. **Stats profile rows** (page 14) — distance to telomere/centromere,
   gene density, GO/KEGG cargo, deleterious burden Δ. Reads from the
   `stats_profile.json` overlay.
4. **Diamond detector status** (page 1) — presence/absence + splitting
   band info per candidate.
5. **Regime breadth + transition boundaries** — from
   `state.candidateList[i]` directly (already there per turn 48 work).
6. **Two-track regime IDs and per-track active bands** — already
   partially in bundle (turn 50); double-check completeness in
   `_bundleEnrichmentForCandidate` output.
7. **Inversion age** (Task A + Task B) — blocked on the
   `inversion_age_v1.json` HPC emitter shipping. See
   `SPEC_inversion_age_atlas_surface.md`.

Effort: ~20 atlas lines per helper, ~150 total. Tests follow the
turn-118 pattern (sandbox + behavioural + scalar round-trip).

### Priority 2 — LD decay overlay (sibling to focal-vs-bg)

`SPEC_ld_decay_overlay.md` is drafted. Lives below the focal-vs-bg
widget on pages 11 and 16. Reuses the existing `fast_ld` stack. The
natural next visualisation after focal-vs-bg.

### Priority 3 — Implement `SPEC_inversion_age_atlas_surface.md`

Two display surfaces: (a) page 5 catalogue gains an `age_my_year`
sortable column (Task A, between-inversion ranking via dXY); (b) page
3 candidate focus gains a per-group within-class age bar chart (Task
B, within-inversion ranking via π).

Blocked on:
1. HPC emitter `STEP_C01f_d_emit_age_json.py` (~80 lines, spec'd in §X
   of the SPEC doc) producing `inversion_age_v1.json`.
2. Modify the existing `STEP_C01f_c_burden_regression.R` v9.7 §4b to
   emit `age_my_year` alongside the existing `age_mya` (one-line
   change: `age_my_year = dXY / (2 * mu_year)` where `mu_year =
   3e-9`).

Once those two HPC pieces are done, atlas-side work is ~400 lines.

### Priority 4 — Cheat30 (turn 119) HPC integration

Atlas-side is done. Quentin needs to:
1. Run `cheat30_gds_by_genotype.R` on per-chrom BEAGLE dosage matrices.
2. Wrap the per-candidate output as `cheat30_gds_results_<chrom>.json`
   per the wire format in `SPEC_age_origin_panel.md` §2.
3. Pre-compute the KDE densities for the ridgeline plot HPC-side (the
   `pair_density` field in the JSON) — atlas does NOT redo KDE in JS.

Once that JSON is loaded, the page-3 panel populates automatically.
Manuscript bundle's `_bundleAgeOriginBlock` also activates.

### Priority 5 — R-side patch for turn 120 (PCA scree)

`PATCH_export_precomp_lam_top_k.md` documents the change. Quentin
applies to `export_precomp_to_json_v3.R` on LANTA, regenerates one
chrom (LG28 is the canonical proving ground), drops the new JSON in
the atlas, sees 5 bars instead of 2 + hint. Then regenerates the rest
at his convenience.

No urgency — atlas-side fallback (2 bars + ratio) is already
informative. The 5-bar version adds elbow shape.

### Priority 6 — Q-completion chips (cross-cutting feature)

Quentin sketched this in turn 119: small `Q5 ⓘ` chip in panel headers,
hover for the 39 Q5 keys with resolved/missing/aspirational status.

Needs the registry-key TSV loader from
`compute_candidate_status.R` (Quentin pasted the source in turn 119)
wired into the atlas first. The 352-key registry has hard-coded
"aspirational" lists per question that distinguish "ought to have"
from "ought to have eventually" for completion accounting.

Once the loader exists, the chip is a small wrapper that reads the
registry keys for the candidate-under-cursor and renders a tooltip.

This is its own spec doc, not yet drafted.

### Priority 7 — Strategic / multi-species direction

Five spec docs from earlier sessions in `specs_todo/` describe a
larger architectural direction:

- `SPEC_OVERVIEW_multispecies_architecture.md` — the umbrella
- `SPEC_hypothesis_registry_and_multispecies.md` — Q-question registry
  extended to multi-species claims
- `SPEC_comparative_te_breakpoint_fragility.md` — TE-mediated fragility
  as a manuscript hook
- `SPEC_busco_anchors.md` — BUSCO anchor markers
- `SPEC_phylogenetic_tree_integration.md` — adding a tree view

Strategic shift, not tactical. Worth reading before the next major
direction change.

### Priority 8 — Card flip + bottom-right inset (UX polish)

Two parked items from the session:

- Card flip on candidate page (turn 119): Quentin sketched, said "this
  is a detail". Wait until age-origin panel is in production and we
  see if the page is actually crowded.
- Bottom-right scree-inset slot (turn 120): originally meant for DA
  eigenvalues from the screenshot, but DAPC isn't applicable to local
  PCA. Slot stays empty. Could host a per-band silhouette diagnostic
  ("how good is the K-means clustering at this window?") if Quentin
  asks for it.

---

## Atlas state slots added across this session (for context)

These are the new top-level `state.*` fields that didn't exist before:

| Field | Turn | Purpose |
|---|---|---|
| `state._focalVsBg` | 117 | `{ csRadiusBp: 2_000_000 }` — focal-vs-bg widget radius for cs_breakpoints |
| `state.cheat30Results` | 119 | Per-chrom dict of cheat30 GDS results for the age & origin panel |
| `state.screePlotEnabled` | 120 | Boolean for the PCA scree inset toggle |

These are in addition to fields already documented in the prior
session's reference (turn 117 inherits from turns 1–116, see
`ATLAS_REFERENCE_for_phase4b_synthesis.md`).

---

## Test inventory

Five test suites, all in `/home/claude/work/build/`. Each is a
standalone Node script; run with `node test_<name>.js`.

| Suite | Assertions | Covers |
|---|---|---|
| `test_turn117_integration.js` | 15 | focal-vs-bg widget integration |
| `test_dotplot_hover_only.js` | 14 | dot plot pure-hover model |
| `test_turn118_enriched_bundle.js` | 82 | manuscript bundle Tier-1 enrichment |
| `test_turn119_age_origin.js` | 74 | age & origin panel + cheat30 enrichment |
| `test_turn120_scree_inset.js` | 33 | PCA scree-plot inset |
| **Total** | **218** | |

All green. No regression at any point in the session.

The vm-sandbox helper-extraction pattern (used in turns 118/119/120)
works well for behavioural tests: pull the helper source out of the
atlas with a brace-balanced regex, exec in a Node `vm` context with a
fake state, assert against the output. Future tests should follow this
pattern — it doesn't need a real DOM, runs fast, and validates the
actual atlas source rather than a re-implementation.

---

## Files in the bundle

**Core deliverables** (the actual work):
- `Inversion_atlas.html` — the canonical atlas, **59,033 lines**
- `atlas_dotplot.js` — 814 lines, post-117b cleanup (no pin state)
- `atlas_focal_vs_bg.js` — 976 lines, unchanged from prior session

**Specs in `specs_todo/`** (work plan for future turns):
- `README.md`
- `SPEC_age_origin_panel.md` — turn 119
- `SPEC_busco_anchors.md`
- `SPEC_comparative_te_breakpoint_fragility.md`
- `SPEC_cross_species_dotplot.md` — already implemented, kept for reference
- `SPEC_focal_vs_background_widget.md` — already implemented (turn 117), kept for reference
- `SPEC_hypothesis_registry_and_multispecies.md`
- `SPEC_inversion_age_atlas_surface.md`
- `SPEC_ld_decay_overlay.md`
- `SPEC_ncrna_density_layer.md` — already implemented, kept for reference
- `SPEC_OVERVIEW_multispecies_architecture.md`
- `SPEC_pca_scree_inset.md` — turn 120
- `SPEC_phylogenetic_tree_integration.md`
- `later/README.md`
- `later/SPEC_xpehh_track.md`

**Audits / patches**:
- `AUDIT_confirmed_candidate_export_gaps.md`
- `PATCH_export_precomp_lam_top_k.md`

**Changelogs** (one per turn):
- `CHANGELOG_turn115.md` (carried forward from prior session)
- `CHANGELOG_turn116.md` (carried forward)
- `CHANGELOG_turn117.md`
- `CHANGELOG_turn117b_dotplot_hover.md`
- `CHANGELOG_turn118_enriched_bundle.md`
- `CHANGELOG_turn119_age_origin.md`
- `CHANGELOG_turn120_scree_inset.md`

**Tests** (all in repo root):
- `test_turn117_integration.js`
- `test_dotplot_hover_only.js`
- `test_turn118_enriched_bundle.js`
- `test_turn119_age_origin.js`
- `test_turn120_scree_inset.js`

**Reference docs** (carried forward):
- `ATLAS_REFERENCE_for_phase4b_synthesis.md`
- `aggregate_ncrna_density.R`
- `sample_ncrna_density_LG28.json`
- `HANDOFF_2026-05-03_session_summary.md` (prior session, kept for context)

---

## Operating notes for the next session

A few patterns established this session that future Claude should
follow:

### 1. "Read first, propose second"

This was a Quentin preference from earlier sessions and it worked well
across this session. Before making any change to the atlas, always
read the relevant existing code (`grep` for the function, `view` the
file at the right line range, understand the pattern) BEFORE writing
the patch. Several times in this session I avoided silly bugs by
checking first what state slot or DOM container already existed.

### 2. Push back when the user's intuition is right but the framing is off

Two examples in this session:

- Turn 120: Quentin asked for *"the same plots as on the DAPC reference
  image"*. The PCA-eigenvalues inset was a great idea; the
  DA-eigenvalues inset was inappropriate (the atlas doesn't run DAPC).
  Worth explaining the difference and dropping the second inset.
- Turn 119: Quentin floated *"maybe a card flip"* but also *"this is
  a detail"*. Right call to acknowledge it as a future-state idea, ship
  the panel without the flip, and revisit if crowding is actually a
  problem.

### 3. Graceful degradation on missing data

Every new enrichment helper in turns 117/118/119/120 returns
`null`/empty-string when the relevant state isn't loaded. The bundle
gracefully degrades. The panel shows an empty-state hint with the
right next-step instruction (load this JSON, run that pipeline).

This pattern means atlas changes ship NOW even when HPC pipeline work
is months away. The bundle test in turn 118 (Test 1: empty state → 0
blocks emitted, full back-compat) is the canonical assertion of this
contract.

### 4. localStorage persistence pattern

For any new toggle (turn 117 widget radii, turn 119 panel state, turn
120 scree toggle): use `inversion_atlas.<key>` localStorage namespace,
restore on init with try/catch wrapped in a fail-soft default,
persist on change. Pattern is consistent across all turns.

### 5. Keep tests close to source

The vm-sandbox extraction pattern (turns 118/119/120) lets us test
atlas-source helpers without a real DOM. Brace-balanced function
extraction handles even the trickier multi-line literals. Each new
helper should get behavioural tests in this style. Don't shortcut to
"it compiles, ship it" — the regression suite is what gives Quentin
confidence that follow-up turns won't break what we shipped today.

### 6. Spec before code on anything > ~100 lines

Turn 119 (cheat30 panel) and turn 120 (scree inset) both spec'd before
implementing. The specs caught design decisions (per-window vs per-L2
scope, fallback behavior, what NOT to ship) before any code was
written. Saves rework. Future-Claude should keep this pattern for any
major surface.

### 7. Honest uncertainty over false confidence

Several places in this session I told Quentin "I haven't checked X yet"
or "this depends on Y which I don't have visibility into". That's
correct behavior — better to flag the uncertainty than guess and
cement the wrong answer. Quentin appreciated the candor.

---

End of handoff. Atlas is at a clean state: 59,033 lines, 218/218 tests
green, no regressions, all changes documented. Next session can pick up
any of the priority items above.
