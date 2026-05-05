# SPEC — Genome-wide ideogram (28-chrom inversion atlas page)

**Status**: drafted turn 130 final session. Not yet implemented.
Estimated ~1 turn.

**Trigger** (Quentin, multiple chats):
> *"Genome-wide ideogram with 28 chr, colored by tier."*
>
> *"Figure 1 of the paper: chromosomes as horizontal bars, candidates
> as colored segments, sized proportional to span."*

A page in the atlas that shows all 28 chromosomes as horizontal bars
with confirmed + auto-promoted candidates overlaid. Click a candidate
→ jump to that candidate's focal page.

---

## 1. Visual contract

```
[chrom name]  ▬▬▬▮▮▬▬▬▬▬▬▬▬▮▮▮▮▬▬▬▬▬▬▬▮▬▬   length bar (Mb scale)
     LG01     ───────[==]──────────[====]─────────[]──────────  220 Mb
     LG02     ───[==]────────────[==]─────────────                180 Mb
     ...

Legend:
  [==] confirmed Tier 1 (red)
  [==] confirmed Tier 2 (orange)
  [==] confirmed Tier 3 (yellow)
  [==] auto-promoted (dashed border, grey-on-color)
```

Sort: by chromosome name (LG01 → LG28) by default; by candidate count or
total candidate span as alternative orderings.

## 2. Data sources

Read from `state.chromSummary` (per `SPEC_multichrom_load_orchestrator.md`):

```
chromSummary[chrom] = {
  size_bp,
  candidates: [
    { id, start_bp, end_bp, tier, source, confirmed }
  ]
}
```

If chromSummary doesn't exist (single-chrom mode), only the active
chromosome shows up — graceful degradation.

## 3. Interactions

- **Hover** on a candidate segment → tooltip with id, span, tier, source,
  K, n_groups.
- **Click** → load that chromosome (if not active) + jump to candidate
  focal page.
- **Filter strip** at top: tier 1 / 2 / 3, source (manual / auto), confirmed
  (yes / no).
- **Search box** → highlight candidates matching id pattern.

## 4. Implementation slices

### Slice 1 — basic ideogram (~0.5 turn)
- Canvas with one row per chromosome
- Read from chromSummary
- Render confirmed candidates as solid colored boxes; auto as dashed
- Hover tooltip
- Click handler that loads the chrom + navigates

### Slice 2 — filter / search (~0.3 turn)
- Filter strip
- Search highlight
- Live re-render

### Slice 3 — orderings + zoom (~0.3 turn)
- Sort dropdown
- Zoom by candidate count / span
- Per-chrom mini-strip toggle

## 5. Manuscript figure export

The ideogram should export as SVG with:
- Manuscript typography (Arial, 7pt body)
- Configurable color scheme (matches Atlas 1 vs Atlas 5 figure context)
- Caption template

## 6. Open questions

1. **Where in the atlas does this page live?** New top-level page, or
   sidebar widget? Probably a new page given how visually-heavy it is.
2. **What about candidate-system clustering across chroms?** If two
   candidates on different chroms share inheritance (Cramér's V high),
   should they be visually linked? Probably yes, as a secondary
   overlay (curved lines connecting linked candidates) — but defer to
   Slice 4.
3. **Atlas vs manuscript figure**: should the same code emit both?
   Yes — a "PDF export" button on the page.

## 7. Tests

- Synthetic chromSummary with 3 chroms × 5 candidates each → all 15
  render at correct positions.
- Click a candidate → state.activeChrom switches, focal page opens.
- Filter "tier 1 only" → only tier-1 candidates remain visible.
- Search "LG12_inv5" → that one candidate highlights.

## 8. Dependencies

- `SPEC_multichrom_load_orchestrator.md` (chromSummary infrastructure).
- `SPEC_l2_sweep_inheritance.md` (auto-promoted candidates).
- Existing candidate page navigation.
