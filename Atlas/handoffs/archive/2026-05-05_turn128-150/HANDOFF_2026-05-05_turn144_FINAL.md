# HANDOFF — turn 144 — Breeding-readiness card, Turn C (render)

**Date**: 2026-05-05
**Atlas main file**: `Inversion_atlas.html` (70,687 lines)
**Working dir**: `/home/claude/Atlas/Atlas/`
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort on LANTA HPC.
**Supersedes**: turn 143 Turn B handoff.

This is **Turn C** of the 4-turn build of
`SPEC_per_candidate_breeding_readiness_card.md`. Turn A (turn 142)
shipped the cohort_diversity loader. Turn B (turn 143) shipped the
pure data builder + Wilcoxon + advice rule engine. **Turn C ships the
render layer + a help-page Tutorials placeholder** Quentin asked for
this turn. Turn D (next) wires bulk catalogue export.

---

## 0. Cohort discipline (NEVER conflate)

1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — assembly paper only.
2. **226-sample pure C. gariepinus hatchery** — current inversion work.
   K=8 clusters = hatchery broodline structure, NOT species admixture.
3. **C. macrocephalus wild** — future paper.

User: **Quentin Andres** (Kasetsart University Bangkok).

---

## 1. What turn C ships

### 1.1 Six render helpers + DOM dispatcher (~830 LOC, ~line 21326–~22155)

All inserted right after Turn B's `_buildBreedingCard` window export.
All but two are PURE (no DOM, no state mutation), all exposed on
`window.*`.

| Helper | Pure? | Purpose |
|---|---|---|
| `_brHtmlEsc(s)` | ✓ | HTML-safe escape for the renderer's interpolated strings. Mirrors `_svgEsc` (line ~67943) for HTML — separate fn because it's defined earlier in the file. |
| `_renderBreedingCardFROHBarsSVG(burden, opts)` | ✓ | Inline SVG bar chart of per-arrangement F_ROH means with min–max whiskers, median tick, n-labels, Wilcoxon p annotation. Defaults to 520×170; opts override `width` / `height`. Empty-state placeholder when burden unavailable. |
| `_renderBreedingCardHTML(card, opts)` | ✓ | The main body builder. Returns one self-contained `<div class="brc-card">…</div>` with the SPEC §1 sections in order: header strip · karyotype chips · F_ROH (SVG + summary table + Wilcoxon line) · K8 × karyotype contingency · MAF · pairing-advice cards · footer. Defensive against missing fields — null cards / no-h_class / no-cohort all render gracefully. |
| `_breedingCardPrintCSS()` | ✓ | The print stylesheet as a string. `@page A4` with 14mm margins, 11pt body, `page-break-inside: avoid` on `.brc-section` and `.brc-advice-item`, `.brc-no-print` hidden under `@media print`. |
| `_breedingCardPrintHTML(card, opts)` | ✓ | Wraps the body in a complete `<!doctype>` document with the print stylesheet, an optional auto-print script, and an optional dismissable hint. Opts: `auto_print` (default true), `show_hint` (default true). |
| `openBreedingCardPrintWindow(cand, opts)` | ✗ DOM | Top-level dispatcher: builds card → wraps as print HTML → `window.open` popup → writes HTML → auto-print fires. Mirrors the proven `_galleryExportPDF` pattern (line ~67925). Returns the opened window or `false` on failure (popup blocked, null cand, builder error). |
| `_wireCandKaryoPrintCardBtn(cand)` | ✗ DOM | Idempotent injector: adds the "🧬 print breeding card" button to the candidate-focus karyotype-tab toolbar (`#candKaryoPane .ck-toolbar`), inserted before `#ckInfo`. Guards against double-insert. |

Constants:

```js
_BR_KARYO_COLORS       // REF blue / HET teal / INV orange / MID violet / HOM grey
_BR_SEVERITY_COLORS    // strong / warn / info palette for advice items
```

### 1.2 Render contract (SPEC §1 visual hierarchy)

The card body emits these CSS classes / structural anchors so Turn D
can reuse the renderer for bulk catalogue export without changes:

```
.brc-card
  .brc-header             — title strip + meta (chrom · coords · span · K · tier · confidence)
  .brc-section.brc-karyo
    .brc-chip-row         — color-coded count chips (REF / HET / INV / MID / HOM_AMBIG)
    .brc-extra            — classification_summary, regime, multi-haplotype tag, ambig tag
  .brc-section.brc-froh
    .brc-svg-wrap         — _renderBreedingCardFROHBarsSVG output
    table.brc-table       — per-arrangement n / mean / median / sd / IQR / min–max
    .brc-wilcoxon         — REF vs INV (or REF vs all_carriers fallback)
  .brc-section.brc-ancestry
    table.brc-table       — K8 × karyotype contingency, with frac% in parens, totals row
  .brc-section.brc-maf
    .brc-maf-row          — p(REF), p(INV), MAF, minor_allele, n_alleles
  .brc-section.brc-advice
    .brc-advice-item × N  — one per advice item, severity-colored left border + tag pill
    .brc-disclaimer       — SPEC §2 "auto-generated; final breeding decisions require human review"
  .brc-footer             — schema · generated_at · cohort size · atlas reviewer URL
```

### 1.3 The "Print breeding card" button

Wired into `renderCandidateKaryotype` right after `expBtn` (export-TSV)
gets its handler. Lives in the candidate-focus karyotype-tab toolbar.
Click → `openBreedingCardPrintWindow(currentCandidate)` → popup with
the rendered card → browser print dialog auto-opens → user picks
"Save as PDF" as destination.

### 1.4 Help-page Tutorials section (placeholder)

Quentin asked: *"You must add the command lines tutorials and some
tutorials in the help page, we can click and it will open the
tutorials same as how we did for the 30 seconds to get figure 3."*

Status: **placeholder section shipped this turn** — eight cards, all
marked `data-status="pending"`, displaying as `coming soon`. The
authoring framework is in place; tutorial *content* is not. This is
deliberate — Quentin said "I don't know" when I asked which CLI
tutorials to author and in what order, so I shipped the catalogue
stub and queued the authoring as a future turn.

The eight cards seeded:

**Atlas walkthroughs:**
- `discover_first_inversion` — load chrom JSON → scrub → lock K=3 → promote
- `boundary_refinement` — refine [start_bp, end_bp] on page 4
- `breeding_card` — the new Turn C deliverable: print breeding-readiness card
- `manuscript_bundle` — page 5 catalogue → 📝 manuscript bundle

**Command-line / pipeline walkthroughs:**
- `cli_inversion_toolkit` — end-to-end inversion-popgen-toolkit run on LANTA
- `cli_te_density` — `RUN_RD_TE_FULL_BOTH_SPECIES.sh`
- `cli_step_d17_boundary` — STEP_D17 L1/L2 boundary scan
- `cli_module_2b_popstruct` — NGSadmix → PCAngsd → evalAdmix → ngsRelate

To author one, swap the card's `data-status` to `available` and wire
its `data-tutorial-id` to the tutorial launcher (which itself is
queued — needs to mirror the diversity atlas's "30 seconds to figure 3"
overlay pattern).

---

## 2. What was NOT done — deferred to follow-up turns

Quentin uploaded four atlas screenshots this turn flagging additional
issues. After triage I scoped them and asked which to do; he said
*"finish the atlas first but at least make a small part of help page
that is a placeholder for tutorials."* So I shipped Turn C + the
tutorials placeholder and parked the rest. Each of the four is a
small focused turn:

### 2.1 L3 contingency-table toolbar consolidation (page 1)

The toolbar currently sprawls across **three rows** before the slab
panels appear (visible in screenshot 1):
- Row 1: layout buttons + color modes + K=3/K=6/K=3+6 + het/L2-sweep + recluster dropdown
- Row 2: L2/1w/5w/10w/Nw + b.continuity + off/scale + draft + promote/detailed/etc.
- Row 3: live dosage + edit (F/B/C)

That's ~120px of vertical real estate on top of the actual content.
Fix: pure CSS/HTML restructure to one or two rows. ~30 LOC, low risk.

### 2.2 "Windows (N)" button shows no engagement state

The arrow-step-size cycler button at top-right (`📊 Windows (N)`) —
when N ≠ 1, the button is engaged but visually identical to default.
Quentin's expectation: when the user has clicked it and step-mode is
no longer the default 1, the button should reflect that visually.

Constraint to preserve: the existing comment at line ~5418 is
explicit that the button must NOT look like a zoom-mode toggle (so
genome / L1 zoom / L2 zoom can't be confused with it). Solution:
introduce a third state — a subtle highlight different from the zoom
buttons' active style (e.g., a colored border-bottom or accent dot)
that indicates "engaged but not a viewmode". ~15 LOC.

### 2.3 1w / 5w / 10w / Nw mode L3 panel uses old slab styling

This is the biggest of the four. Comparing screenshots 3 (L2 mode,
rich) vs 4 (1w mode, downgraded):

- L2 mode renders rich pane heads with per-pane recluster dropdown +
  layout buttons + H-system legend with `(n=N)` counts
- 1w/Nw mode has bare slab panels, no per-pane controls, no legend

Root cause: `renderL3Panel` (line 42706) is the rich L2-mode
renderer; `renderL3PanelSlab` (line 43230) is the slab-mode renderer.
They diverged. To bring slab mode to L2-mode parity is a focused
rewrite of `renderL3PanelSlab` (and its `slabFocalContentHtml` /
`compareSlabPair` helpers) to emit the same pane structure. Estimated
~150–250 LOC + ~50 tests. Should be its own turn.

### 2.4 G-key panel ↔ popgen page merge

Currently:
- The `popgen` page tab (visible top-right in screenshot 2) is one page
- The G hotkey opens a modal-overlay popup `_gPanelToggle()` with three
  tabs: karyotype / inheritance / manual

Quentin's ask: merge the modal into the popgen page as a second row
of tabs (different background shade), so the popgen page has a single
tab strip per row — popgen tabs in row 1, unified-groups tabs in row
2. The G hotkey then jumps to the popgen page and focuses the
unified-groups row.

Architecturally invasive: removes/hides the overlay, adds a tab strip
to the popgen page, reroutes `_gPanelRenderTab*` to render into
popgen-page DOM. ~400 LOC + ~50 tests. Its own turn.

---

## 3. Tests

`tests/test_turn144_breeding_card_render.js` — **162 / 0** across
13 sections:

1. Source-level definitions (10) — all functions + constants + Turn C banner present
2. Window exports (6) — every renderer hooked to `window.*`
3. `_brHtmlEsc` behavior via output observation (2) — `<script>` is escaped
4. `_renderBreedingCardHTML` structural shape (33) — every SPEC §1 section + sub-element
5. Fallback states (10) — null card, no h_class, no cohort
6. `_renderBreedingCardFROHBarsSVG` shape (15) — viewBox, ticks, colors, n-labels, Wilcoxon, K=4 MID bar adapt, custom dims, null-burden placeholder
7. `_breedingCardPrintCSS` contents (10) — `@page A4`, 14mm margin, every key class, font in pt, `page-break-inside: avoid`, `@media print` hides `.brc-no-print`
8. `_breedingCardPrintHTML` wrapping (14) — !doctype, charset, title, embedded `<style>`, auto_print on/off, show_hint on/off, null-card path
9. `openBreedingCardPrintWindow` mocked DOM (12) — success path: window.open / document.open/write/close fire in order with !doctype HTML; null cand → alert + false; popup blocked → alert + false
10. `_wireCandKaryoPrintCardBtn` DOM idempotence (7) — button created with right id/class/title/text, inserted before `#ckInfo`, second call doesn't duplicate, no-toolbar doesn't throw
11. Real-fixture render (7) — 226-sample synthetic 60/106/60 split renders; embeds 60/106/60 counts; full doc has @page A4
12. Severity styling on advice items (3) — strong/info severities both rendered, color border-left applied
13. Atlas-side wiring + tutorials placeholder (16) — `renderCandidateKaryotype` wires the print button after expBtn; tutorials section heading; grid container; all 8 cards present; all 8 default to `data-status="pending"`; "coming soon" pills; authoring note

Adjacent suites unchanged:
- turn 143 breeding-card compute: **196 / 0**
- turn 142 cohort_diversity loader: **133 / 0**
- turn 141 candidate-bands: **62 / 0**
- turn 140 H-label chip: **45 / 0**
- turn 139 H-label classifier: **143 / 0**

Full turn-numbered suite: **1955 / 0** (up from 1793 / 0 baseline at
turn 143 start; +162 from the new test file).

---

## 4. Atlas state

| | LOC | Tests | Files |
|---|---|---|---|
| Pre-turn (turn 143 baseline) | 69,665 | 1793 | 42 |
| Post-turn (this) | 70,687 | 1955 | 43 |
| Δ | +1,022 | +162 | +1 |

Bundle includes one new test file (`test_turn144_breeding_card_render.js`)
and a smoke-render output (`tests/breeding_card_smoke.html`, 13.4 KB)
that renders the synthetic 226-sample LG28 candidate as a real
print-ready breeding card. The smoke file is also copied to
`/mnt/user-data/outputs/` for direct inspection.

---

## 5. Smoke-render verification

A real render against the cohort_diversity_v1.json fixture with a
synthetic 60/106/60 split (matching the actual reported LG28 karyotype
shape) produced:

```
candidate.id:                LG28_smoke
karyotype counts:            60 / 106 / 60   (REF / HET / INV)
burden.n_resolved:           226             (100% F_ROH lookup hit)
Wilcoxon REF vs INV p:       0.5619          (not significant, expected for synthetic)
delta_mean_inv_minus_ref:    +0.0039         (tiny, expected — F_ROH values uniform)
ancestry rows (K8 clusters): 8               (all K8 clusters represented)
maf:                         0.500           (perfectly balanced — also expected)
advice items:                3               (default + 2 data_pending)
document HTML length:        13,415 bytes
```

The Wilcoxon p of 0.56 here is the synthetic null. On real LG28 data
with the actual planted F_ROH asymmetry, this number will be the
manuscript-quotable p-value.

---

## 6. Backups

```
Inversion_atlas.html.bak_pre_breeding_readiness_turnC   (this turn — pre-Turn-C)
Inversion_atlas.html.bak_pre_breeding_readiness_turnB   (turn 143 baseline)
Inversion_atlas.html.bak_pre_breeding_readiness_turnA   (turn 142 baseline)
Inversion_atlas.html.bak_pre_lines_cand_bands           (turn 141 baseline)
```

`.bak_*` files NOT in bundle. Re-derivable.

---

## 7. What's NOT done — Turn D scope (unchanged from turn 143 handoff)

After this Turn C, Turn D adds:

- "Generate breeding cards" bulk button on page 5 catalogue
- Per-tier filter (default Tier 1 + Tier 2)
- Per-card download as standalone HTML (using `_breedingCardPrintHTML`)
- Bulk export (zip of all cards, or one combined HTML)

The Turn C renderer is fully reusable for Turn D — its output is a
self-contained HTML doc string with no external dependencies (CSS is
inlined, no external image refs, no JS dependencies). Turn D can call
`_breedingCardPrintHTML(card, { auto_print: false })` per candidate
and bundle the resulting strings without any changes to Turn C's code.

Estimated Turn D scope: ~150 LOC + ~30 tests.

---

## 8. Where to start the next chat

### Option 8a — Turn D (RECOMMENDED, finishes the breeding-card story)

Build the bulk export per §7. The catalogue page already has the
"📝 manuscript bundle" button; "Generate breeding cards" can sit
alongside it. Estimated ~150 LOC + ~30 tests. After Turn D, Atlas 5
Part B is fully shipped — chat-`c03fc41e` framing ("converts the paper
from a population genomics study into a hatchery management resource")
has its complete deliverable.

### Option 8b — Tackle Quentin's queued atlas fixes (§2)

Four small focused turns, each independent:

1. **L3 toolbar consolidation** (~30 LOC, low risk, big visible win)
2. **Windows-(N) button engagement state** (~15 LOC)
3. **1w/Nw L3 panel parity with L2 mode** (~150–250 LOC, biggest)
4. **G-popup → popgen page merge** (~400 LOC, most architectural)

Recommended order if going this route: 1, 2, 3, 4 — increasing scope,
each delivers visible improvement on its own.

### Option 8c — Author the first tutorial

Pick one of the eight `data-status="pending"` cards in the help-page
Tutorials section, write its step-by-step content, build the tutorial
launcher overlay (mirroring the diversity-atlas "30 seconds to figure
3" pattern), wire that one card to `data-status="available"`. The
breeding-card tutorial is the obvious first candidate since the
underlying functionality just shipped this turn.

### My recommendation

**8a.** Turn D is the closest thing to "done" — one more 150-LOC turn
finishes the four-turn build, hatchery managers have something
concrete to scan, and the manuscript Atlas-5-Part-B section becomes
fully demonstrable. The four atlas fixes (8b) are real but they're
quality-of-life polish, not critical-path; they slot in after Turn D
without needing any of its work as a prerequisite.

---

## 9. Honest framing

**What turn 144 actually delivered:**

- Six render helpers + a DOM dispatcher + a toolbar wirer, all
  individually unit-tested
- A complete print-ready HTML/CSS render of the SPEC §1 visual contract
  — every section in the mockup is implemented
- A reusable inline SVG bar chart that adapts column count to the
  candidate's haplotype regime (3 bars for K=3 clean, 4 bars when
  HOM_MID present)
- An `@page A4` print stylesheet with proper page-break controls so
  long advice / ancestry tables don't split mid-section
- A DOM dispatcher that mirrors the proven `_galleryExportPDF` popup
  pattern — same keyboard shortcut behavior the user is already
  trained on for SVG/PNG/PDF gallery export
- An idempotent toolbar button injection that survives multiple
  `renderCandidateKaryotype` calls without duplicating
- A help-page Tutorials placeholder section with eight named cards
  covering atlas walkthroughs and CLI/pipeline walkthroughs — the
  authoring framework is ready for content

**What it deliberately didn't deliver:**

- Tutorial *content*. Quentin said "I don't know" when I asked which
  ones to author and in what order. Cards are placeholders, all
  marked `coming soon` so users see what's coming without being
  promised content that doesn't exist.
- The tutorial launcher overlay. Cards are `cursor: not-allowed` and
  `pointer-events`-passive — clicking does nothing. Authoring the
  launcher mirrors the diversity-atlas pattern; it's a focused next-
  turn job once we know which tutorial gets authored first.
- The four atlas fixes Quentin flagged from the screenshots (toolbar
  sprawl, Windows-(N) styling, slab-mode parity, G-popup merge).
  Each is its own focused turn — see §2.

**Manuscript impact:**

Page 7 (karyotype tab, candidate-focus) now has a "🧬 print breeding
card" button that produces a manuscript-quality, reviewer-scannable
print PDF for any promoted candidate with a locked karyotype + a
loaded cohort_diversity_v1.json. The Wilcoxon p-values surfaced are
manuscript-quotable. Hatchery managers / aquaculture-genomics
reviewers who don't read the full paper can grab one PDF per
inversion and decide whether to use it for breeding decisions. That
was the original chat-`c03fc41e` framing for Atlas 5 Part B.

Walk the map carefully, respect cohort discipline, don't break the
test suite. Turn C's render layer is structurally sound, fully
covered, and produces real output verifiable by visual inspection of
`tests/breeding_card_smoke.html`.
