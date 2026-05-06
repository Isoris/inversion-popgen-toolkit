# SPEC — Unified G-panel: karyotype groups + inheritance groups + manual groups

**Status**: SHIPPED in turns 135 (slice 1) + 136 (karyotype slice 2). See
`tests/test_turn135_g_panel_slice1.js` and `tests/test_turn136_g_panel_karyotype_slice2.js`.
Drafted turn 128d. Replaces and
absorbs the implicit "G panel" Quentin referred to as the manual-groups
sidebar section.

**Trigger**: Quentin's request (turn 128d):
> "Our G panel of groups didn't load groups of inversion karyotypes and I
> think both specs are somehow related. Because we could have 1 tab in the
> G popup group panel for karyotype of candidate and one for inheritance
> groups too. ?"

The user wants karyotype-derived groups (HOM_REF / HET / HOM_INV per
candidate) to appear in the same G-panel that today only holds manual
groups, and asks whether inheritance groups should also live there as a
sibling tab.

---

## 1. The five "group" concepts in the atlas — distinguish before unifying

Five distinct concepts share the word "group". A single popup that
unifies them naively will confuse the user. The spec must keep them
crisp:

| # | Concept | Source | Scope | Persistence |
|---|---|---|---|---|
| 1 | **K-means cluster** | `getL2Cluster(l2idx)` | per-L2 | precomp + cache |
| 2 | **Karyotype regime** | derived from K=3 + PC1 ordering | per-candidate | live UI |
| 3 | **Inheritance group** | `inheritanceGroupClustering()` | cross-candidate | live UI |
| 4 | **Family group** | ngsRelate JSON `family_id` | cohort-level | data layer |
| 5 | **Manual group** | user-curated in sidebar | cohort-level | localStorage |

### What goes in the G-panel

- ✅ **Karyotype regimes** (concept 2): per-candidate, per-band group of
  samples with their HOM_REF / HET / HOM_INV assignment.
- ✅ **Inheritance groups** (concept 3): cross-candidate clusters of
  shared-membership patterns.
- ✅ **Manual groups** (concept 5): the existing in-sidebar set.

### What does NOT go in the G-panel

- ❌ **K-means clusters** (concept 1) are a geometric primitive, not a
  user-facing group. They surface elsewhere (color-mode bars, K-cycle
  buttons, contingency tables).
- ❌ **Family groups** (concept 4) are kinship, not phenotype. They have
  a separate ontology (related vs unrelated, hub vs singleton) and a
  separate color-mode (`'family'`). Mixing them with karyotype groups
  in one tab strip would conflate "this fish has the same kinship as X"
  with "this fish has the same karyotype as Y" — they're orthogonal
  questions.

---

## 2. Why Quentin's intuition is correct (mostly)

The karyotype + inheritance + manual concepts are all "**phenotype
groupings of the cohort that the user can inspect, recolor by, or
export from**." They share:

- A list of named groups
- Each group has a sample-membership set
- The user wants to colorize / filter / export by group
- Group definitions are stable across L2 navigation (unlike K-means,
  which redefines per-L2)

So a unified popup with three tabs is right for these three. Family
groups operate similarly but on kinship — separate enough that a
fourth tab is at most an inert link, not a peer surface.

---

## 3. Where the panel lives today + where it should live

### Today

The sidebar has a `manualGroupsList` div (line ~5258) and the compact
panel mirrors it (`manualGroupsListCompact`, line ~6196). There's no
karyotype-group surface in the sidebar — karyotype assignments are
visible only on page 4 (karyotype/tier).

### Proposed

A single floating popup (`#groupsDialogPanel`, mirroring `regimeDialogPanel`)
opened from a sidebar button or hotkey. Three tabs:

```
┌────────────────────────────────────────────────────┐
│  Groups [karyotype] [inheritance] [manual]    [✕]  │
├────────────────────────────────────────────────────┤
│  Tab content                                       │
│  …                                                  │
└────────────────────────────────────────────────────┘
```

Why a popup (vs sidebar): three tabs each with a sample list of up to
226 entries don't fit in the sidebar's 280px width without wrapping or
scrolling that fights the rest of the sidebar's vertical real estate.
A popup centred over the page is the natural form factor.

The existing sidebar `#manualGroupsList` stays — it's a quick-glance
preview. The popup becomes the editing surface and the cross-tab
viewer.

---

## 4. Tab 1 — Karyotype groups

### Trigger and scope

When the user has a focused candidate (`state.candidate` non-null),
the karyotype tab shows that candidate's per-band sample memberships
as named groups:

- `g0 / HOM_1` (PC1-low band)
- `g1 / HET` (PC1-mid band, K=3)
- `g2 / HOM_2` (PC1-high band)

For K=6 candidates, six bands; the regime mapping is fuzzier and the
column shows the K=6 band with a tooltip "no biological regime
assigned at K=6".

### What the tab shows

Per group:

- The group name (editable — defaults to `g0/HOM_1`, etc.)
- A swatch in the K-cluster color
- Sample count badge (`N=84`)
- A collapsed sample list — click to expand to see the CGAs.
- Action buttons: **Make manual group from this** (copies the band to
  Tab 3 as a manual group, snapshot at the moment of click), **Color
  by this** (sets `state.colorMode='cluster'` and locks via
  `state.lockedRefL2 = candidate.ref_l2`), **Export TSV** (one CGA
  per line).

### When no candidate is focused

Empty-state message: "Select a candidate in the catalogue or pin one
on page 1 to see its karyotype groups here."

### When candidate is multi-track (rare)

Karyotype tab shows the focal track by default; a small `track 1 / 2`
toggle picks the other.

---

## 5. Tab 2 — Inheritance groups

### Trigger and scope

When `state.confirmed.length >= 2` AND the user has run inheritance
clustering (today: implicit, runs on the I·g labels strip in page 1).
Surfaces the cross-candidate Jaccard cluster cuts as named groups.

Each group represents "samples that share band-membership patterns
across multiple confirmed candidates" — the manuscript-relevant
inheritance signal.

### What the tab shows

Per group:

- Group ID (`I1`, `I2`, …)
- Group name (editable)
- A swatch in the inheritance palette (Tol qualitative)
- Per-candidate band breakdown: `cand_42 → g1 (n=80)` rows, one per
  candidate the group touches.
- Sample count badge (total across candidates).
- Actions: same as karyotype (make manual / color by / export).

### Settings

- **Threshold slider**: cosine-distance cutoff (`_IGC_DEFAULT_COSINE_DIST_THRESHOLD`).
  Re-runs `inheritanceGroupClustering()` and refreshes the tab. This
  is the same parameter as on the I·g strip; the popup is the editing
  surface.
- **Compute button**: explicit `[ Compute ]` action so the user
  controls when the (potentially expensive) clustering runs.

### When fewer than 2 confirmed candidates

Empty-state: "Confirm at least 2 candidates on page 2 to compute
inheritance groups."

---

## 6. Tab 3 — Manual groups

The existing `manualGroupsList` content but in a wider format:

- Sample list per group with CGA + family info
- Drag-drop sample reassignment between manual groups (deferred)
- Add / rename / delete (existing actions)
- Import/export TSV
- "From tracked" creation (existing)

This is mostly a re-host of the sidebar UI in a roomier popup, with
better sample-list rendering.

---

## 7. State design

```js
// Existing
state.manualGroups   = [...];        // unchanged
state.candidate      = {...};        // unchanged
state.confirmed      = [...];        // unchanged
state.familyPalette  = {};           // unchanged

// New
state.gPanelOpen     = false;        // popup visibility
state.gPanelTab      = 'karyotype';  // 'karyotype' | 'inheritance' | 'manual'
state.gPanelGroups   = {             // ephemeral cache of derived groups
  karyotype:    null,                // populated when popup opens with a candidate
  inheritance:  null,                // populated when popup opens with confirmed >= 2
};
```

The karyotype + inheritance groups are **derived**, not stored. The
popup re-derives them when it opens. The user's "Make manual group
from this" action snapshots a derived group into `state.manualGroups`
where it persists via the existing localStorage path.

---

## 8. UI affordances

### Opening the popup

- A small `[ G ]` button in the sidebar header (above the Atlas tools
  folder), or in the per-sample-lines panel header where I·g pills
  live now.
- Hotkey `g` (lowercase) opens the popup. Capital `G` is unused.

### Closing

- ✕ button, click outside the panel, Esc key. Mirror the `regimeDialogPanel`
  pattern.

### Tab switching

- Click a tab header. Persist `state.gPanelTab` to localStorage so
  reopen restores last view.

### Color-by-this consistency

When the user clicks "Color by this" inside the popup, the popup stays
open — but a small toast/banner in the popup head reads "PCA recolored
by [group name]. Click again to revert."

---

## 9. Why karyotype groups didn't load before

The current manual-groups system has no mechanism to ingest derived
groups. Each manual group is a static `{ name, sample_ids, color }`
created either by user action ("from tracked") or import. Karyotype
assignments are computed live from `getL2Cluster(candidate.ref_l2)` and
were never plumbed into the groups system.

The fix is **not** to auto-populate manualGroups with karyotype groups
(that would mix snapshots with live derivations). Instead the karyotype
tab reads `state.candidate.ref_l2 → getL2Cluster() → groupLabels[]` on
every render and shows the live groups. The user explicitly snapshots
into manualGroups via "Make manual group from this" when they want
persistence.

---

## 10. Implementation slices

### Slice 1 (smallest viable, ~1 turn)
- Add `state.gPanelOpen / state.gPanelTab` slots + localStorage restore.
- Build `#groupsDialogPanel` overlay following the `regimeDialogPanel`
  template.
- Tab 3 (manual): host the existing manualGroups UI in the new popup.
  The sidebar manualGroupsList stays as a preview.
- Sidebar `[ G ]` button opens the popup; hotkey `g` too.
- Tabs 1 + 2 are stub "Coming next turn" empty-states.

### Slice 2 (~1 turn)
- Tab 1 (karyotype): live render from `state.candidate`. "Make manual
  group from this" action.
- Empty state when no candidate focused.

### Slice 3 (~1 turn)
- Tab 2 (inheritance): live render from `inheritanceGroupClustering()`.
  Threshold slider + Compute button.
- Empty state when < 2 confirmed candidates.

### Slice 4 (~0.5 turn) — polish
- Per-tab Export TSV.
- "Color by this" toast/banner.
- Badge count on tab headers (e.g. `karyotype (3)` if focal candidate
  has 3 bands).

---

## 11. Open design questions (for Quentin to resolve before Slice 1)

1. **Sidebar button location**: above Atlas tools, in lines header, or
   both?
2. **Hotkey**: `g` or `G`? Currently neither bound.
3. **K=6 karyotype labels**: today the per-K6-band → biological-regime
   mapping is unsettled. Should the karyotype tab show 6 unnamed bands
   at K=6, or refuse and message "K=3 only"?
4. **Family groups**: include a 4th tab as a read-only family directory?
   Or keep family entirely separate (current default; the family color
   mode + sidebar dropdown handle it).

---

## 12. Tests

- Source: panel element + button + hotkey wiring.
- Tab switching persists `state.gPanelTab` to localStorage.
- Karyotype tab: with a synthetic candidate + getL2Cluster stub,
  derives 3 groups with correct sample counts.
- Inheritance tab: threshold slider triggers `inheritanceGroupClustering()`
  re-run; group count updates.
- "Make manual group from this": adds a new entry to
  `state.manualGroups` with the snapshot membership.
- Manual tab parity with the sidebar's existing manual-groups list.

---

## 13. Why this is medium-sized, not small

- Three derivation pipelines (karyotype lookup, inheritance compute,
  manual snapshot) need their own tab renderers.
- The popup component itself is moderately complex (overlay + tabs +
  per-tab body + actions).
- Inheritance compute is potentially expensive — needs a Compute
  button rather than reactive recompute on every threshold drag.
- Cross-tab consistency (e.g. snapshot from karyotype creates a
  manual group that should appear immediately in tab 3 without a
  reopen) needs deliberate state plumbing.

Estimated total: **~3 turns** across 4 slices. Slice 1 alone is a
~1-turn win — it adds the popup scaffold and re-hosts the existing
manual groups, which is already a UX improvement (more space for the
list). Slices 2 + 3 deliver the karyotype + inheritance content.

---

## 14. Summary for Quentin

> Yes, your intuition is right: karyotype groups + inheritance groups +
> manual groups all belong in one G popup with three tabs. Family groups
> stay separate (different ontology — kinship not phenotype).
>
> The reason karyotype groups didn't load into the existing manual-groups
> system is that manual groups are static snapshots while karyotype
> groups are live derivations from the focal candidate. The fix is a new
> popup with a live karyotype tab, not auto-populating the manual list.
>
> Slice 1 (popup scaffold + manual tab re-host) is ~1 turn. Karyotype
> and inheritance tabs come in subsequent slices.

---

## 15. Late-turn diagnosis: why the inheritance/grouping panel is currently broken

Quentin flagged after this spec was drafted:
> "The grouping panel is bugged — it doesn't auto select or compute when
> we click it for the candidate. It doesn't find the 'groups' also
> shouldn't it be related or linked to the groups in candidate page as
> well? Or its a different grouping?"

**Short answer**: it IS the same grouping as the candidate page (both
derive from `candidate.locked_labels`). The bug is a registry mismatch
in the data pipeline, not a conceptual difference.

**Long answer**:

The inheritance compute pipeline reads from a different registry than
the one the user actually populates:

| Registry | Type | Populated by | Read by |
|---|---|---|---|
| `state.candidateList` | Array | user actions: promote, lock-and-promote, draft confirm | catalogue, page 2, candidate strip, exports |
| `state.candidates` | Dict (id→cand) | Save-Session import ONLY | `_gatherActiveCandidatesForInheritance` (line 34362) → `runInheritanceCompute` → I·g pills |
| `state.candidates_detailed` | Dict | mode-switch detailed pipeline only | same gatherer, detailed branch |

So when the user saves 5 candidates interactively, `state.candidateList`
has 5 entries but `state.candidates` is still empty. The I·g compute
silently finds < 2 items and returns null. **The grouping panel is
reading from the wrong slot.**

The candidate-page band-composition cards work fine because they read
`state.candidate.locked_labels` directly from a single candidate; the
inheritance pipeline tries to traverse `state.candidates` and finds it
empty. Same biological signal, divergent reader paths.

### Fix options (this spec absorbs option c)

a. **Bridge approach** (cheap): make `_gatherActiveCandidatesForInheritance`
   read from `state.candidateList`. One function change. Risk: anywhere
   downstream that mutates `state.candidates` (Save-Session import) loses
   its effect — would need to merge into candidateList on import.

b. **Single source of truth refactor** (clean): converge on
   `state.candidateList` everywhere; deprecate `state.candidates` /
   `state.candidates_detailed` dicts; rebuild `state.candidate /
   .candidate_detailed / .activeMode` dispatch on top of array filters.
   2–3 turns of careful work; touches many sites.

c. **Auto-sync** (pragmatic middle ground, RECOMMENDED): on every
   mutation of `state.candidateList` (add / remove / restore from
   localStorage), rebuild `state.candidates` as a derived index. Slim
   helper `_rebuildCandidateRegistries()` called from the existing
   4–5 mutation entry points. ~1 turn.

### Pre-Slice for the G-panel

Adopting option (c) means the G-panel Slice 1 should **NOT** ship the
popup until the registry bridge is in place; otherwise the karyotype
+ inheritance tabs will be empty-state on first load even with saved
candidates.

**Updated ordering for the next chat**:

| Order | Step | Effort |
|---|---|---|
| Pre-Slice | `_rebuildCandidateRegistries()` bridge, hook into candidateList mutations, tests | ~0.3 turn |
| Slice 1 | G-panel popup scaffold + manual tab re-host | ~1 turn |
| Slice 2 | Karyotype tab (works because bridge exists) | ~1 turn |
| Slice 3 | Inheritance tab (computes correctly because bridge exists) | ~1 turn |
| Slice 4 | Polish — Export TSV, Color-by-this toast, badge counts | ~0.5 turn |

### Tests for the Pre-Slice

- **Empty list**: `state.candidateList = []` → `runInheritanceCompute()`
  returns null.
- **One candidate**: add one to candidateList → `state.candidates` now
  has a matching entry → `runInheritanceCompute()` returns null
  (needs ≥ 2).
- **Two candidates**: add second → `state.candidates` has two entries
  → `runInheritanceCompute()` returns a result with > 0 groups.
- **Remove a candidate**: removing from candidateList also removes from
  `state.candidates`.
- **Save-Session round-trip**: import candidates via Save-Session (the
  legacy path that populates `state.candidates` directly) — bridge does
  the inverse, copying back into candidateList for consistency.
