# Batch 3 — extraction notes

**Pages shipped**: page3, page9, page17, page18, page21
**Sub-atlas**: `inversion_catalogue/`
**Foundation regression**: 455/455 tests still pass, `shared/` byte-identical to handoff.
**New tests**: 5 files, 65 assertions, all passing.

---

## Discrepancies vs. HANDOFF_BATCH_3.md

The handoff doc's **"What it does"** column was inaccurate for two pages.
Line ranges (which the doc also gives) are authoritative; descriptions are not.
The actual page contents per the line ranges:

| Page  | Lines       | Doc said                  | Actually is                                         |
|-------|-------------|---------------------------|-----------------------------------------------------|
| page3 |  7261–7368  | "5 catalogue" — table     | ✓ correct (catalogue toolbar + table shell)         |
| page9 |  7782–7812  | "11 confirmed" carousel   | ✓ correct (carousel shell only)                     |
| page17|  8110–8127  | Marker readiness panel    | ✗ **Statistical profile** (`spBody`)                |
| page18|  8134–8157  | Genome-wide linkage table | ✗ **Marker readiness panel** (`mpBody`)             |
| page21|  7822–7859  | Manual karyotype groups   | ✗ **Annotation cockpit** (canvas-driven)            |

The two "marker / linkage" descriptions in the doc were swapped, and page21's
description of "manual karyotype groups list" does not match the legacy HTML
at that line range — that region is the annotation cockpit. Extraction
followed the line ranges.

---

## What's in each module

Each page's `pageN.js` follows the same shape:
- A doc header block (legacy line refs, public entry, dependency list)
- A `state` shim (reads `window.state` if present, else `{}`)
- The verbatim legacy code body wrapped between BEGIN/END markers
- Public ES-module exports at the bottom

### page17.js (stats profile)
- Verbatim legacy lines 28445–29306 (~860 lines).
- All 24 `_sp*` helpers + `SP_DEFAULT_ROWS` + the `_renderStatsProfilePage`
  entry are in one file.
- Public entry: `renderStatsProfilePage()`.
- Internal helpers re-exported for tests + merge chat.

### page18.js (marker readiness panel)
- Verbatim legacy lines 29307–30160 (~850 lines).
- All `_mp*` helpers + Tier-1/2/3/4 logic + AF scoring (`_mpScoreVariantAf`).
- Public entry: `renderMarkerPanelPage()`.
- Note: `_mpDeriveAutoPanel` is read by sibling page17 — exported.

### page21.js (annotation cockpit)
- Verbatim legacy lines 46970–47616 (~640 lines).
- Canvas drawing (`_annoCockpitDraw`), keyboard nav (`_annoCockpitOnKey`),
  click-to-jump (`_annoCockpitOnClick`), readouts, linkage panel, hap panel.
- Public entry: `refreshAnnotationCockpit()`.
- The legacy block ended with `window.foo = foo` re-exports; those are
  replaced with ES-module `export { ... }` at the bottom. Otherwise the
  body is unchanged.

### page9.js (confirmed carousel)
- **No verbatim legacy code exists.** The HTML shell (lines 7782–7812)
  has #confirmedNavBar / #confirmedNavPrev / #confirmedNavNext etc., but
  `grep -nE 'confirmedNav' legacy/Inversion_atlas.html` returns only the
  HTML lines — no JS handlers, no renderer. The page is a stub even in
  the legacy build.
- Shipped as a thin shell that shows the empty-state message and reports
  the count of confirmed candidates.

### page3.js (catalogue)
- **Same situation.** `renderCatalogue` is referenced 9 times in legacy
  (`if (typeof renderCatalogue === 'function') renderCatalogue();`) but
  never defined. Also confirmed:
  ```
  $ grep -nE 'function renderCatalogue|renderCatalogue\s*=|window\.renderCatalogue' \
        legacy/Inversion_atlas.html
  (no matches)
  ```
- The toolbar HTML (~110 lines) is fully wired with IDs, but no
  JavaScript populates `#catBody` / `#catHead` or wires the filter,
  sort, export, regime, or view-mode buttons.
- The ONE piece of catalogue toolbar JS that DOES exist in legacy is
  `_wireCatalogueBreedingExportBtns` at lines 23668–23715, plus its
  ~480-line dependency block (Turn 146 breeding-card export, lines
  ~23236–23715). That block is too large to extract into page3 without
  overstepping batch scope and is left as `TODO_MISSING(_wireCatalogue...)`
  for the merge chat.
- Shipped as a thin shell that shows the legacy "Load a JSON to populate
  the catalogue." empty-state message.

---

## TODO_MISSING markers — frequency table

44 total markers across the 5 modules. Grouped here by likely owner so the
merge chat can route them.

### Likely shared/ promotion candidates (used by 2+ pages)
| Function                                    | Pages    | Note |
|---------------------------------------------|----------|------|
| `_esc`                                      | 17, 18   | HTML-escape helper. High frequency across all batches probably. Strong promotion candidate. |

### Likely owned by other batches
| Function                                    | Owner    | Note |
|---------------------------------------------|----------|------|
| `_csGetSyntenyBlocks`                       | Batch 5  | Cross-species synteny lookup |
| `_csPermutationTest`                        | Batch 5  | Cross-species permutation test |
| `renderCandidateFocus`                      | Batch 1  | page2 candidate-focus renderer (page9 carousel reuses it) |
| `_gatherActiveCandidatesForInheritance`     | shared?  | Legacy line 41196; gathers active candidates for current chrom |
| `computeTrackedLinkageProjection`           | shared?  | Legacy line ~46912; tracked-fish band purity |
| `_wireCandidateHaplotypeAnnotations`        | Batch 1  | Band-annotation UI (probably page2) |
| `candidateHaplotypeAnnotationsHtml`         | Batch 1  | HTML builder for haplotype-annotation panel |

### Internal to this batch (cross-module)
| Function                                    | Owner            | Note |
|---------------------------------------------|------------------|------|
| `_mpDeriveAutoPanel`                        | page18 (this)    | Already exported from page18.js — merge chat just wires the import in page17.js |

### Page3 catalogue rendering — entire renderer is missing in legacy
Listed for completeness; none of these exist anywhere in `legacy/Inversion_atlas.html`.
Implementing them is its own future task, not a merge-chat fix:
- `_buildCatalogueRows`
- `_filterCatalogueRows`
- `_sortCatalogueRows`
- `_paintCatalogueRow`
- `_wireCatalogueToolbar`
- `_exportCatalogueTSV` / `_exportCatalogueMarkdown` / `_exportCatalogueJSON`
- `_exportCatalogueGallerySVG` (and PNG/PDF siblings)
- `_openRegimeRegistryDialog` / `_assignSelectedToRegime`
- `_promoteSelectedAsCandidate`
- `makeShelfLDPanel` (legacy comment says: "atlas_turn7.js's makeShelfLDPanel()")
- `makeLDSplitPanel` (legacy comment says: "atlas_ld.js's makeLDSplitPanel()")

### Page3 breeding-card export — exists in legacy, just not extracted in this batch
Lives at legacy lines ~23236–23715. The merge chat or a follow-up batch should
extract it as a sibling helper (probably `inversion_catalogue/breeding_export.js`):
- `_wireCatalogueBreedingExportBtns` (legacy 23668)
- `_exportBreedingCardsHTML` / `_exportBreedingCardsJSON`
- `_filterCandsForBreedingExport`
- `_buildBreedingCardsCombinedHTML` / `_buildBreedingCardsJSONBundle`
- `_BREEDING_EXPORT_TIER_MODES`

### Page9 carousel rendering — entire renderer is missing in legacy
- `_renderConfirmedCarousel`
- `_wireConfirmedCarouselNav`

---

## Decisions

1. **Followed line ranges, not the doc's "What it does" column**, when the two
   disagreed. The doc's labels for page17/18/21 were wrong; line ranges were
   right.

2. **Extracted page17, page18, page21 verbatim** with a thin ES-module wrapper
   (header + `state` shim + verbatim body + bottom exports). No edits to
   legacy logic. The merge chat decides what becomes shared and what stays
   per-module.

3. **Shipped page3 and page9 as shells** matching legacy reality. Both were
   already shells in legacy (the renderers don't exist there either),
   so faithful extraction meant shipping shells. Documented this loudly in
   the module headers and in this notes file so the merge chat doesn't
   spend time looking for renderers that don't exist.

4. **`window.foo = foo` re-export blocks at the bottom of the legacy
   annotation-cockpit section** were replaced with ES-module
   `export { ... }`. Behavior preserved (modules can still see these
   helpers; the API surface is just `import` instead of `window`).

5. **Tests exercise pure helpers where possible.** page17 / page18 validators
   and page21's `_ackBandColor` / `_annoCockpitCandidateAtCursor` are
   pure functions and have direct assertions. The DOM-dependent renderers
   are only checked for "doesn't throw in Node-no-DOM environment".

6. **One bug found in test assumptions, fixed before shipping.** Initial
   tests asserted `_spIsValidProfileJson(null) === false` but the legacy
   validator returns the `null` operand from short-circuit evaluation
   (`return parsed && ...`), not strict `false`. Tests rewritten to use
   truthy/falsy checks, which is what the legacy code's call sites do too.

---

## Files shipped

```
Atlas/inversion_catalogue/
├── BATCH_3_NOTES.md       (this file)
├── page3.html             (108 lines, sed -n '7261,7368p')
├── page3.js               (stub; renderCataloguePage + initCataloguePage)
├── page9.html             ( 31 lines, sed -n '7782,7812p')
├── page9.js               (stub; refreshConfirmedCarousel + initConfirmedCarousel)
├── page17.html            ( 18 lines, sed -n '8110,8127p')
├── page17.js              (verbatim from legacy 28445–29306)
├── page18.html            ( 24 lines, sed -n '8134,8157p')
├── page18.js              (verbatim from legacy 29307–30160)
├── page21.html            ( 38 lines, sed -n '7822,7859p')
└── page21.js              (verbatim from legacy 46970–47616)

Atlas/tests/
├── test_catalogue_page3.js   ( 5 assertions)
├── test_catalogue_page9.js   ( 6 assertions)
├── test_catalogue_page17.js  (17 assertions)
├── test_catalogue_page18.js  (19 assertions)
└── test_catalogue_page21.js  (18 assertions)
                              ─────
                              65 assertions, all passing
```

Foundation: 455/455 still pass. `shared/` and the 8 foundation tests
(`test_shared_*.js`, `test_modular_smoke.js`, `test_legacy_parity.js`)
are byte-identical to the input handoff bundle.
