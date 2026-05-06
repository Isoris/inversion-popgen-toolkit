# SPEC — Review surfaces for auto-promoted candidates + lineage blocs

**Status**: Slice 0 (infrastructure) SHIPPED in turn 130 follow-up
(see `tests/test_turn130_auto_review_infra.js`,
`tests/test_turn130_lineage_compute.js`, `tests/test_turn130_lineage_ui.js`).
The L2-sweep auto-promote producer that this spec was waiting on has since
landed (turns 133–134); the remaining UI slices (lineage-bloc grouping tab in
the G-panel, review queue) are now unblocked. Drafted turn 130.

**Trigger** (Quentin, turn 130 follow-up):
> *"They could (autopromoted inversion candidate blocs) eventually be
> added as a grouping tab in the grouping G menu. So that we can rapidly
> review it and it would let us confirm them to add them to catalogue
> because otherwise it will spam the catalogue or add them to the
> catalogue but last rows on the list and outlines dashed."*
>
> *"Same for inheritance blocs / regimes."*

---

## 1. The problem this solves

L2-sweep + lineage compute will produce two parallel streams of
**proposed structure** that the user hasn't confirmed:

1. **Auto-promoted candidates** from `runL2SweepInheritance()` —
   inversion-candidate blocs the algorithm thinks are real, gated by
   silhouette + ≥2-group + dedupe rules. Output of
   `SPEC_l2_sweep_inheritance.md` §3.
2. **Lineage blocs / regimes** from `runLineageCompute()` — fish-trajectory
   lineages. Each lineage is a coherent group of fish that share a
   chromosome-wide K-means band-trajectory. Output of `SPEC_distant_band_concordance_fish_trajectory.md`.

Without a review surface, both will either:

- **Spam the existing catalogue** if dumped in raw → user can't tell
  algorithm-proposed from user-confirmed at a glance, the catalogue
  becomes noisy.
- **Stay invisible** if hidden → user can't act on them, defeats the
  purpose of computing them.

The fix: **two complementary review surfaces** that work together.

## 2. Surface A — dashed-outline rows in existing UI

### 2.1 What it is

`cand.source === 'auto_l2_sweep'` (or `'auto_lineage'`) candidates render
in the existing candidate strip + page-2 candidate list with a
**distinct visual treatment**:

- **Dashed border** (1px dashed instead of solid, slight opacity reduction)
- **Robot icon** `🤖` prefix on the candidate name
- **Sorted to the bottom** of any sort that includes `confirmed` —
  unconfirmed auto-promotes never crowd out user-confirmed candidates
- **No I·g pill participation** until confirmed (don't pollute the
  inheritance compute with unreviewed proposals)

When the user clicks "confirm" on an auto-promoted candidate, the
`source` flips to `'auto_l2_sweep_confirmed'` (or just removes the auto
prefix) and the candidate joins the standard catalogue.

### 2.2 Where it shows up

| Surface | Treatment for auto candidates |
|---|---|
| Candidate strip (page 1 top) | Dashed border, 🤖 prefix, sorted to right end |
| Page-2 candidate list | Dashed row, 🤖 prefix, sorted below confirmed |
| Catalogue (page 3) | NOT shown unless catalogue mode = `auto` (new view) |
| Cross-candidate matrix | NOT included until confirmed (data integrity) |
| I·g pills | NOT painted (don't claim inheritance for unreviewed) |
| Lines panel candidate-bands (§ SPEC_lines_panel_candidate_bands) | NOT painted (would overlay wrong color on regions the user hasn't accepted) |

### 2.3 State / persistence

- Existing `cand.source` field carries the discriminator
  (`'manual'` | `'auto_l2_sweep'` | `'auto_lineage'` |
  `'auto_l2_sweep_confirmed'` etc.)
- Existing `cand.confirmed` boolean flips to `true` on user confirm
- New per-chrom `state.dismissedAutoIds: Set` in localStorage at
  `pca_scrubber_v3.dismissedAutoIds.<chrom>` so dismissed auto-promotes
  don't reappear on the next sweep. Already in
  `SPEC_l2_sweep_inheritance.md` §3.

## 3. Surface B — G-panel review tabs

### 3.1 What it is

The G-panel popup (per `SPEC_g_panel_unified_groups.md`) gains two new
tabs at the END of its tab strip:

- **`auto`** — review queue for auto-promoted candidates from L2-sweep.
- **`lineages`** — review queue for fish-trajectory lineages.

These are not "groupings" in the same sense as the existing karyotype /
inheritance / manual tabs; they're **triage queues**. Distinct visual
treatment (e.g., dashed tab border, "review" subtitle) makes the
distinction clear.

### 3.2 The `auto` tab

For each `cand` in `state.candidateList` where
`cand.source && cand.source.startsWith('auto_')` and not yet confirmed:

```
[🤖] L2:13  · 4.2–4.7 Mb  · K=3, silh=0.42, n_groups=2
     Bands: g0 (n=98)  g1 (n=84)  g2 (n=44)
     [Confirm]  [Dismiss]  [Inspect]
```

- **Confirm** → flips `cand.confirmed = true`, drops `auto_` prefix from
  source, removes from this tab, appears in the regular candidate strip
  in solid (not dashed) form.
- **Dismiss** → adds `cand.id` to `state.dismissedAutoIds`, removes from
  list. Next L2-sweep won't re-suggest the same L2.
- **Inspect** → navigates to that L2 on page 1 and opens the L3 mini-PCA
  panel — same flow as clicking a candidate in the strip.

Bulk actions header: **Confirm all** / **Dismiss all** with confirmation.

### 3.3 The `lineages` tab

For each lineage in `state.lineageResult.lineage_id_per_sample`:

```
[Lineage 0]  n_fish = 47   chromosome-wide
     ▮▮▮▮▮▮▮ (hue swatch)
     [Track]  [Pin as group]  [Dismiss]
```

- **Track** → set `state.tracked = [those fish]` (replaces selection).
  User sees the lineage in the per-sample-lines panel.
- **Pin as group** → create a manual group with these fish. Adds to the
  manual-groups list (which the existing G-panel manual tab manages).
- **Dismiss** → flag the lineage as not interesting (no practical effect
  on the compute; just hides it from this list).

The lineages tab has no "confirm" action — lineages aren't candidates.
They feed into the manual-groups system if the user wants to keep them.

### 3.4 Tab visibility

- `auto` tab visible ONLY when `state.candidateList` contains at least
  one unconfirmed auto-promoted candidate. Otherwise hidden.
- `lineages` tab visible ONLY when `state.lineageResult` exists and has
  ≥2 lineages.

Empty G-panel never has these tabs cluttering the strip.

## 4. Why both surfaces, not just one

The G-panel is **modal** (popup). It blocks the user from looking at the
chromosome while reviewing. Good for batch triage; bad for "I'm working
on candidate I3 and I want to know what auto-promotes are nearby."

The dashed-outline strip rows are **inline**. They sit alongside the
user's confirmed candidates as they work. Good for ambient awareness; bad
for batch triage (can't bulk-confirm 12 auto-promotes from the strip).

Both serve the same data, different workflows. Cheap to maintain because
the underlying state (`cand.source`, `cand.confirmed`,
`state.dismissedAutoIds`) is the single source of truth.

## 5. Implementation slices

### Slice 0 — infrastructure (~0.5 turn)  ← SHIPPED turn 130 follow-up
- [x] `_isAutoCandidate(cand)` central predicate (returns true when
      `cand.source.startsWith('auto_')` AND `!cand.confirmed`)
- [x] `_gatherActiveCandidatesForInheritance` skips auto candidates
      (so I·g pills + cross-candidate matrix exclude unreviewed proposals)
- [x] Defensive fallback when `_isAutoCandidate` not in scope (sandbox
      tests still work)
- [x] CSS rule `.cand-list-item.is-auto` (dashed border, opacity 0.85)
- [x] CSS rule `.cli-auto-prefix` for the 🤖 prefix
- [x] `refreshCandidateListUI` sort modifier — auto candidates last
- [x] `refreshCandidateListUI` emits `is-auto` class + 🤖 prefix
- [x] Tests for predicate edge cases, gather-skip behavior, sort order

### Slice 1 — full dashed-outline visual treatment (~0.5 turn, after L2-sweep ships)
- [ ] Filter modifier in cross-candidate matrix to skip auto candidates
- [ ] Lines panel candidate-bands renderer skips auto (when that ships)
- [ ] Confirm action flips source + drops dashed treatment

### Slice 2 — `auto` G-panel tab (~0.5 turn)
- [ ] Tab markup + visibility gate
- [ ] Per-row card layout with Confirm / Dismiss / Inspect
- [ ] Bulk Confirm-all / Dismiss-all buttons with confirmation prompt
- [ ] Tests: tab visibility, action handlers.

### Slice 3 — `lineages` G-panel tab (~0.5 turn)
- [ ] Tab markup + visibility gate
- [ ] Per-row card with hue swatch + Track / Pin as group / Dismiss
- [ ] Tests: visibility on lineageResult presence, Track→state.tracked,
      Pin→manualGroups.

## 6. Open design questions

1. **Default new-state on auto-promotes**: `cand.confirmed = false`
   (as the L2-sweep spec already says) means the auto-promotes don't
   appear in the confirmed-only walkthrough on page 8. That's correct
   default behavior. Confirm this is what Quentin wants.

2. **Dismiss vs delete**: dismissing an auto-promote should NOT trigger
   `removeCandidateFully()` because that purges it from `state.candidateList`
   entirely. Instead, dismiss should: drop from candidateList AND add
   to `dismissedAutoIds` so it doesn't come back. Same end state but
   the dismissed-set is what blocks re-promotion.

3. **Lineage "Pin as group"**: the existing manual-groups system stores
   group definitions. Pinning a lineage means freezing the current
   `state.lineageResult.lineage_id_per_sample` for that lineage_id into
   a static fish-id list. After the freeze, recomputes won't update
   the pinned group — it's a snapshot.

4. **Should the catalogue (page 3) get a parallel auto view?** The
   current spec says no — catalogue is an L2-row registry, not a
   candidate registry, so it doesn't naturally have "user-confirmed
   vs auto" semantics. But if Quentin wants auto-promotes catalogued,
   add a `mode = 'auto'` to the existing L2 / L1 / fav / L3 selector.

5. **Lineage blocs as inversion candidates**: a lineage could in
   principle become a candidate (a region where these fish co-segregate
   is one). The lineages tab does NOT auto-promote candidates from
   lineages — that's a separate spec if Quentin wants it. The Pin-as-
   group action is the bridge.

## 7. What this is NOT

- **Not a replacement for the existing candidate strip / page-2 list.**
  Auto candidates show up there too, just dashed.
- **Not a refactor of the catalogue.** The catalogue is L2-row-based
  and stays that way. Auto candidates aren't catalogue rows; they're
  candidateList entries with a special source tag.
- **Not new compute.** Both review tabs read existing
  `state.candidateList` + `state.lineageResult`. No new algorithm.
- **Not the G-panel base implementation.** That's
  `SPEC_g_panel_unified_groups.md` Slice 1. The auto + lineages tabs
  in this spec are slotted as Slices 4 + 5 of the G-panel work.

## 8. Tests

### Slice 1 (dashed outline)
- Auto candidates render with dashed border CSS class.
- Auto candidates sort below confirmed in candidate strip + page-2 list.
- I·g pills compute excludes unconfirmed auto candidates.
- Lines panel candidate-bands renderer excludes unconfirmed auto.
- Confirm action flips source + drops dashed treatment.

### Slice 2 (auto tab)
- Tab visible only when ≥1 unconfirmed auto candidate exists.
- Confirm-row → `cand.confirmed = true`, source loses `auto_` prefix.
- Dismiss-row → removed from candidateList, added to dismissedAutoIds,
  persisted to localStorage.
- Bulk Confirm-all / Dismiss-all run on the visible list with
  confirmation prompt.

### Slice 3 (lineages tab)
- Tab visible only when `state.lineageResult.n_lineages >= 2`.
- Track → state.tracked = lineage members, lines panel re-renders.
- Pin → manualGroup created with the lineage's fish, persists.
- Dismiss → lineage hidden from this list (no compute side-effects).

## 9. Dependencies

This spec depends on:

- `SPEC_l2_sweep_inheritance.md` Slice 1 (auto-promote mechanic) — must
  ship first.
- `SPEC_g_panel_unified_groups.md` Slice 1 (G-panel scaffold) — must
  ship before Slices 2 and 3 here.
- `SPEC_distant_band_concordance_fish_trajectory.md` Slice 1 (compute)
  — already shipped turn 130.

The dashed-outline Slice 1 can ship independently after the L2-sweep
auto-promote lands; the G-panel tabs need both prerequisites.
