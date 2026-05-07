# HANDOFF SUPPLEMENT — corrections + lessons from the 5 batch chats

This document captures things that came up during the parallel extraction
that the original `HANDOFF_MERGE.md` would not have known. Read both.

---

## Correction 1: page5 is the HELP page, not multi-species

**The original `HANDOFF_BATCH_5.md` was wrong** about what page5 is.

### What I told Batch 5

> page5 (lines 8159–9316): "Multi-species comparison page. **Big** (1158 lines of HTML). Pairwise dotplots, synteny graph, breakpoint annotations."

### What's actually at lines 8159–9316

A static **help/quick-reference page** (`<button data-page="page5" data-stage="help">` at legacy line 5142). Pure declarative HTML — tables of tabs, hotkeys, schema definitions, pipeline diagrams. No JS render function in the legacy file.

Batch 5 caught the discrepancy and handled it correctly:
- Extracted lines 8159–9316 as `page5.html` (1158 lines, the help content)
- Made `page5.js` a stub with explanatory comments
- Tagged `PAGE5_META.stage = 'help'`, `static: true`
- Documented in `BATCH_5_NOTES.md`

### Where the multi-species page actually lives

**`page16b`** (legacy lines 8033–8100) is the multi-species classification cockpit.
The renderer is `_renderMultiSpeciesPage` at legacy line 27497. **This is correctly extracted** in `inversion_comparative/page16b.js` (2,403 LOC).

### What this means for you (merge chat)

When building `inversion_comparative.html`:
- `page5` is the help page → it goes in `inversion_catalogue.html` actually, NOT comparative. It's stage `help`. Move it.
- Update the comparative HTML's tab bar to NOT include page5.
- Update the catalogue HTML's tab bar to INCLUDE page5 with stage `help`.

Or alternatively: since page5 is just static help content, you could put it in any sub-atlas (or all of them). Quentin's call. Default: catalogue.

---

## Correction 2: stage filtering is more nuanced than my original mapping

The `data-stage` attribute on legacy tab buttons doesn't map cleanly 1-1 to sub-atlases. Here's the actual stage distribution I see in the legacy:

```
discovery       → page1, page12, page15, page2, page19, page3
refinement      → page11, page4, page8
classification  → page21, page6, page7
synthesis       → page9, page17, page18
compare         → page16, page16b
help            → page5, page10, page_overview, page_sv_evidence (mixed?)
```

vs. the original sub-atlas mapping I wrote:

```
inversion_discovery   → page1, page12, page15, page2, page8, page19
inversion_review      → page11, page_sv_evidence, page4, page7, page6
inversion_catalogue   → page3, page9, page21, page17, page18, page10, page_overview
inversion_comparative → page16, page16b, page5
```

**The merge chat decision**: don't try to perfectly align stages with sub-atlases. The sub-atlas split was driven by *workflow*, not stage. Each sub-atlas's HTML should include only the buttons whose pages it owns:

```
inversion_discovery.html  buttons: page1, page12, page15, page2, page8, page19
inversion_review.html     buttons: page11, page_sv_evidence, page4, page7, page6
inversion_catalogue.html  buttons: page3, page9, page21, page17, page18, page10, page_overview, page5  (+ help)
inversion_comparative.html buttons: page16, page16b
```

Just literally include the per-page tab `<button>` HTML for the pages that sub-atlas owns. Don't filter by `data-stage`.

---

## Correction 3: the "9 confirmed" page is page9 in catalogue, but it consumes confirmed candidates produced by the discovery scrubber

Page9 (catalogue) is downstream of page1 (discovery) — it shows confirmed candidates. They're connected via `state.candidateList` (a `cross_atlas` slot in `SLOT_REGISTRY`).

This means: **the cross-atlas state contract is doing real work** — Quentin promotes a candidate on page2 (discovery), and page9 (catalogue) reads that from review/inversion/sessions.

The merge chat doesn't need to do anything special here — the foundation already supports it. Just be aware that the user workflow crosses sub-atlases.

---

## Lessons learned

### What worked

1. **The "mark TODOs and ship" rule** kept extraction time short. The 5 chats did real work (~14,800 LOC of JS, ~3,700 LOC of HTML extracted) but didn't get bogged down chasing dependencies.

2. **The non-overlapping filename convention** (`inversion_<phase>/page<N>.{js,html}`) eliminated all merge conflicts. Five tarballs untar into the same directory cleanly.

3. **The "shared/ is read-only" rule** held perfectly. None of the 5 chats touched the foundation. All 455 foundation tests pass on the merged tree.

4. **Per-page minimal tests** (just "module imports without throwing") gave fast feedback that extraction was syntactically valid. 194 pass on the batches.

### What didn't work as well

1. **The handoff docs had errors** that the chats had to detect and document (the page5/multi-species mixup). Future workflow: have batch chats explicitly verify line-range descriptions BEFORE extracting.

2. **Some pages were over-described.** `BATCH_3` had 7 pages assigned but the descriptions in the handoff were thin — chat had to do a lot of legacy spelunking to figure out what each page actually does.

3. **The TODO_MISSING tail is long** (180 single-reference markers). Most of these are tiny utility functions — `xToPx`, `yToPx`, `withAlpha`. Resolving them is mechanical but voluminous. Could have been faster if the handoff had said "extract these helper categories upfront as `shared/dom.js`, `shared/geom.js`, etc."

### What to do better next time

If you ever do this kind of refactor again:

- **Have one prep chat extract ALL utility helpers first** (`_esc`, `xToPx`, `yToPx`, `withAlpha`, `formatTrackVal`, etc.) into `shared/dom.js`, `shared/geom.js`, `shared/format.js` BEFORE splitting into batches. Cuts the tail of single-reference TODOs by ~50%.
- **Auto-validate handoff line ranges** by asking each batch chat to grep for `data-page="<assigned_page>"` and confirm the surrounding HTML matches what the handoff describes. Batch 5 caught this manually but wasted a little time.
- **Have batch chats expose a uniform `renderPage<N>(state)` entry point** (already de-facto convention). Document in the handoff explicitly so the merge chat's dispatcher is trivial.

---

## Final notes for the merge chat

You have everything you need:
- Foundation: clean, tested, locked
- Batch outputs: clean, tested, all imports OK
- TODO inventories: pre-computed in `Atlas/docs/merge_inputs/`
- Original handoffs: in `handoff_docs/` for reference

The work is bounded. Most of it is mechanical. Keep your scope tight (don't try to fix everything; defer when in doubt). Foundation tests stay green throughout — that's the contract.

Ship `Atlas_FINAL_2026-05-05.tar.gz` when done.
