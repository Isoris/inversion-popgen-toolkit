# Handoff — Population Atlas v1.1 + Diversity Atlas v2.2

## What shipped

### Population_atlas.html (v1.1)

Pages 4 / 5 / 6 (Heterozygosity / Diversity / Inbreeding) converted from
"scaffold + planned panels" to **summary + cross-link landing cards**.

- New tab labels with `↗` indicator and updated tooltips
- Green callout button on each page with prominent "Open Diversity Atlas → ..." link
  - Page 4 → `Diversity_atlas.html#page1` (Samples, where H lives)
  - Page 5 → `Diversity_atlas.html#page2` (Chromosomes, where θπ lives)
  - Page 6 → `Diversity_atlas.html#page5` (ROH, where F_ROH composition lives)
- 4-cell summary grid on each page with cohort headline numbers
  (mean H, median H, ρ(H, F_ROH), KW H × K=8 — analogous on the other pages)
- "What's in the Diversity Atlas for this topic" depth-pointer card listing the
  exact tabs and tables a reader can drill into
- Caveats / bin-schema preserved (folded SFS, ROH dual-bin design)
- About page (page 7) updated to four-atlas family per ADR-14, with a new
  card explaining the scope decision

Title and header meta updated to reflect the new role:
*Population Atlas — Cohort Summary & Cross-Atlas Hub*

### Diversity_atlas.html (v2.2)

Three small but important updates:

1. **Hash routing** added: `Diversity_atlas.html#page2` now opens directly
   on the Chromosomes tab (was: opened on last-visited tab from
   localStorage). Required for the Population Atlas's green callouts to
   work as expected. Also responds to `hashchange` for in-page navigation.
2. **Reverse cross-link buttons** on the About page (tab 7): yellow ←
   Population Atlas button + blue → Inversion Atlas button, mirroring the
   green callouts on the other side. Discoverable bidirectional graph.
3. **Roadmap card** on tab 8 updated to mark items 3 / 4 / 5 as resolved
   (with strikethrough), preserving the history. Two new items added
   for next session: kinship heatmap surfacing and per-sample θπ ribbons.

## Validated end-to-end

JSDOM smoke test on both files: 0 errors, 0 warnings. Cross-link hrefs
verified:
- Population → Diversity_atlas.html (×3 hash anchors), Inversion, Genome
- Diversity → Population, Inversion, Genome (Genome via dropdown only)

Hash routing tested across all 8 page anchors plus a bogus hash
(falls back to localStorage / page1). hashchange events also fire
correctly.

## Cross-atlas relationship after this session

```
              ┌─────────────────────────┐
              │   Population Atlas      │  ← cohort QC, families, breeding-program lens
              │   (v1.1, scaffold +     │     • pages 1-3: Samples / QC / Families
              │    summary cards)       │     • pages 4-6: summary + green callout to Diversity
              └────────┬────────────────┘     • page 7: about + 4-atlas map
                       │
                       │ green callouts (3 buttons, hash-anchored)
                       ▼
              ┌─────────────────────────┐
              │   Diversity Atlas       │  ← deep popgen-reviewer view
              │   (v2.2, shipped)       │     • 8 tabs, 26 tables, 22 plots
              │                         │     • 226 strip plots overlaid on K=8 boxes
              │                         │     • all-226 / top-100 stack composition toggle
              │                         │     • hash routing for inbound deep-links
              └────────┬────────────────┘
                       │ ← back-link button on About tab
                       ▼
              ┌─────────────────────────┐
              │   Population + Inversion (peer atlases via shared dropdown)
              └─────────────────────────┘
```

Two atlases, two depths, one source of truth — the MODULE_3
supplementary tables. The Population Atlas's pages 4/5/6 are now the
"breeding-program collaborator's quick lookup"; the Diversity Atlas
is the "population-genetics reviewer's full deliverable".

## Files in this turn

- `Diversity_atlas.html` (1.4 MB) — adds hash routing + reverse cross-links + updated roadmap
- `Population_atlas.html` (63 KB) — pages 4/5/6 rewritten as summary cards; About updated to 4-atlas family
- `manuscript_paragraphs_module3.md` (already shipped earlier this session, unchanged) — Methods + Results paragraphs

## What's still open

- **Inversion Atlas not touched this session** — it still has the green/yellow/orange/blue dropdown widget but doesn't yet have a callout from any of its tabs back to the Diversity Atlas's hotspot tab. Easy 5-minute addition next session: on Inversion Atlas's per-candidate page, add a "θπ at this locus →" button to `Diversity_atlas.html#page3` filtered to the relevant chromosome.
- **Kinship heatmap** — ngsRelate output (J1 estimator) is available but not yet JSON-ified. Adding it as a Diversity Atlas Tab 4 panel would let the kinship structure (the actual *cause* of the all226 → pruned81 KW collapse) be visible in the same atlas as the F_ROH × K=8 boxes. Probably the highest-value next addition.
- **Per-sample θπ ribbons** — need cluster-side `thetaStat do_stat` batch run. Each sample's per-window θπ → 226 × 1,895-window matrix → could be surfaced as a click-to-detail expansion on Tab 1 sample rows.
- **Genome Atlas** — still un-built. The fourth atlas in ADR-14, scheduled for chromosome-scale assembly stats, synteny, gene tracks, and TE landscape. No urgency unless the assembly paper resurfaces.
