# HANDOFF — turn 159 — propose-all karyotype groups + auto-jump to inheritance

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (74,104 lines, +331 LOC from 73,773)
**Working dir**: `/home/claude/work/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

**Closes**: Quentin's request to make the G-panel a one-shot grouping
flow — "push G when on my candidate then it auto propose the groups for
the karyotype, push enter, it auto fills it and merge tabs and look at
inheritance and so on."

**Picked up from**: post-turn-158 reorganized working tree
(`HANDOFF_2026-05-05_reorganization.md`).

---

## 0. What this turn ships

A unified UX layer on top of the existing turn-136 karyotype tab. Before
this turn, the user had to click "+ manual" once per band to promote
each band to the manual group system, then manually navigate to the
inheritance tab to see the cross-candidate I·g view. Now:

1. **G** opens the panel on the focal candidate (existing behavior, turn 135).
2. **🎯 propose all groups** button → generates per-band proposals and
   shows them in a strip with classification chips + counts.
3. **Enter** accepts → addToManualGroup ×N → auto-jumps to inheritance tab.
4. **Esc** cancels the proposal strip.

Three core changes to `Inversion_atlas.html`:

### `_gpKaryoProposeAll(candidate)` — pure helper

Returns:
```js
{
  candidateId, K,
  n_proposals, n_empty_bands, n_ungrouped,
  classifier_ran,                 // true iff _classifyHLabelBands ran
  proposals: [{
    name,                         // <candId>_<bandLabel>[_<classification>]
    bandIdx, bandLabel,
    classification,               // 'HET' | 'HOM' | 'AMBIGUOUS' | 'NO_DOSAGE' | null
    color, sampleIdxList, n,
  }, ...],
}
```

Naming convention: `<candTag>_<bandLabel>` + optional classification
suffix. The suffix is appended ONLY when the dosage classifier
disagrees with the band's nominal label (e.g. `I3_HOM_REF_HET_dosage`
when a band labeled HOM_REF is classified HET by dosage; this surfaces
the conflict in the group name rather than burying it). When dosage
agrees with the label, the suffix is omitted to keep names short.

`AMBIGUOUS` and `NO_DOSAGE` always get a suffix because they're
provisional regardless of label.

### `_gpKaryoAcceptProposals(bundle)` — applies proposals + tab jump

Calls `addToManualGroup` per non-empty proposal. Returns:
```js
{
  accepted, skipped,
  statuses: [{ name, ok, group_id?, n_members?, reason?, error? }, ...]
}
```

Then, if `state.gPanelOpen`, switches `state.gPanelTab = 'inheritance'`
and re-renders the modal. Inheritance compute itself is unaffected
(manual groups are not inputs to inheritance — inheritance reads
candidate `locked_labels`); the jump is pure UX, putting the user in
front of the I·g cross-candidate view immediately after promoting.

### UI: button + proposal strip + Enter/Esc handler

- **Button** in karyo tab footer, next to existing 🔒/⬇/🔍/📊 actions.
  Bordered green to make the primary action obvious.
- **Proposal strip** rendered in karyo tab body when
  `state._gpKaryoProposals` is set and matches the focal candidate ID.
  Shows per-band rows with color swatch, classification chip, proposed
  name, and `n=N` count, plus accept/cancel buttons.
- **Enter/Esc IIFE** `_wireGPanelKaryoProposalKeys` — keydown handler
  that fires only when (a) panel open, (b) on karyotype tab,
  (c) `_gpKaryoProposals` set, (d) target is not an input field.
  Independent of the existing `_wireGPanelHotkey` IIFE so binding
  scope stays narrow.

---

## 1. What this does NOT do (deliberate scope discipline)

- Does NOT modify `addToManualGroup` itself. Existing exclusivity rules
  (a CGA in any new group is removed from prior groups) apply
  unchanged. If you propose all bands of candidate I3 *and* I3 has
  CGAs that were already in a manual group from candidate I2, those
  CGAs migrate. This matches the existing per-band "+ manual" semantics.
- Does NOT change inheritance compute. `runInheritanceCompute` reads
  `c.locked_labels`, not `state.manualGroups`. The inheritance tab
  jumps to show what's already there; it does not re-trigger compute.
- Does NOT touch K-means clusters (concept 1 in the
  `SPEC_g_panel_unified_groups.md` taxonomy) — those remain a per-L2
  geometric primitive.
- Does NOT add a "collect ALL candidates' groups in one shot" pass.
  This turn solves the per-candidate one-shot flow. Cross-candidate
  bulk collection would be a separate turn (matches the gap I flagged
  earlier — `_collectAllGroups()` accessor).

---

## 2. Files touched

```
Inversion_atlas.html                                +331 LOC
  - _gpKaryoProposeAll                              new function
  - _gpKaryoAcceptProposals                         new function
  - _gPanelRenderTabKaryotype                       +button +proposal strip
  - karyo tab event-listener block                  +3 button handlers
  - _wireGPanelKaryoProposalKeys                    new IIFE (Enter/Esc)
  - window.* exports                                +2 exports

tests/test_turn159_propose_all_karyotype_groups.js  new (45 assertions)
```

No existing tests modified. No tests left in a broken state.

---

## 3. Test results

**Single test**: 45 / 0 across 8 sections:
1. Source-pattern checks (12 assertions)
2. `_gpKaryoProposeAll` behaviour on synthetic 6-sample K=3 candidate (16)
3. Empty-bands skipped (2)
4. `_gpKaryoAcceptProposals` calls addToManualGroup ×N (6)
5. Tab auto-jump to inheritance only when panel open (2)
6. Empty/null inputs handled gracefully (3)
7. Naming sanity (2)
8. Re-render hook fires after accept (1)

**Full sweep**: **3036 / 0** (was 2991 / 0 at bundle close).
Zero regressions. JS syntax clean. HTML parser 0 errors.

The 10 environment-broken pre-turn-128 tests (missing fixture files at
`/home/claude/work/build/Inversion_atlas.html` etc) remain unchanged
— same as bundle baseline.

---

## 4. What Quentin should exercise

1. **The flow**: load a chrom, focus a candidate (page 1 → promote &
   lock), press **g** → karyo tab → click **🎯 propose all groups**
   → strip appears with one row per non-empty band, each with a
   classification chip → press **Enter** → groups are added, panel
   jumps to inheritance tab.
2. **Naming check**: if dosage classifier ran (you ran "compute stripe
   quality" on the dosage heatmap), check the names. A clean K=3
   with no surprises should give `<candId>_HOM_REF`, `<candId>_HET`,
   `<candId>_HOM_INV`. A band where dosage disagrees with the karyo
   label should pick up `_HET_dosage` or `_HOM_dosage` suffix.
   `AMBIGUOUS` / `NO_DOSAGE` bands always get a suffix.
3. **Cancel**: hit Esc instead of Enter while the strip is showing.
   Strip clears, no groups added.
4. **No-classifier case**: do the flow on a candidate before loading
   dosage data. `bundle.classifier_ran === false`, names get only
   the bandLabel (no classification suffix). Names still work — the
   acceptance still works.
5. **Empty-bands edge case**: propose-all on a degenerate candidate
   where K=3 but only one band is populated. Strip shows 1 row.
   `n_empty_bands = 2` reported in header.

---

## 5. What's NEXT (relative to the broader queue)

After this turn, the per-candidate one-shot flow works. The remaining
grouping work in priority order:

1. **`_collectAllGroups()` accessor** — programmatic snapshot for ALL
   candidates × all fish × all group systems (karyotype + inheritance
   + manual). Unblocks bulk-popstats, Save Session round-trip of
   manual + inheritance state, manuscript Supplementary Data exports.
   ~120 LOC + ~80 LOC tests. This is the "collect all inversions in
   one time" half of Quentin's intent.
2. **Status update for `SPEC_g_panel_unified_groups.md`** — that spec
   is in `specs_done/` but the turn 159 propose-all flow is a strict
   superset of what shipped in 135–136. Either add a "Slice 4
   (turn 159)" section to that spec or write a small standalone spec
   capturing the propose-all UX design as a turn-159 reference doc.
3. **Save Session payload extension** — `_buildSessionPayload` doesn't
   currently include `state.manualGroups`. After 159 the user is much
   more likely to make many manual groups in a sitting; losing them
   on F5 is bad UX. ~30 LOC.
4. **Apply propose-all to inheritance tab** — the inheritance tab has
   per-group "+ make manual" but no "promote all groups at once"
   button. Same pattern as 159, ~80 LOC.

---

## 6. Honest framing

**What's solid:**
- 3036 / 0 sweep. Clean JS + HTML.
- Pure helpers separable from UI; sandbox-tested with stubs.
- The accept path goes through the existing `addToManualGroup` so
  exclusivity, persistence, and re-render contracts are all unchanged.
- Tab jump only fires when panel is open — no side effect when called
  programmatically from elsewhere (e.g. a future bulk path).

**What's risky:**
- Naming: when the user has multiple candidates with overlapping
  karyo labels (e.g. two K=3 candidates both have HOM_REF/HET/HOM_INV),
  proposed names will be `I1_HOM_REF`, `I2_HOM_REF`, etc — distinct
  by candidate prefix. But if the user re-runs propose-all on the
  same candidate after editing labels, the same names will be
  reused, which means addToManualGroup will EXTEND the existing
  groups (its idempotent contract). This is probably what Quentin
  wants but worth flagging — there's no "are you sure you want to
  re-propose" guard.
- Tab jump assumes inheritance compute already ran or will run on
  next render. If the user is on a chromosome with only the focal
  candidate (no other locked_labels), the inheritance tab will show
  the empty-state "Need ≥2 candidates" message after the jump.
  That's correct UX — they just promoted the first batch — but it
  looks empty until they promote a second candidate.
- The Enter handler is a global keydown listener. The guard chain is
  `gPanelOpen && gPanelTab==='karyotype' && _gpKaryoProposals && !inputField`
  but if any future code sets `_gpKaryoProposals` outside the karyo
  tab flow, Enter would fire. Safe today, narrow blast radius.

End of handoff.
