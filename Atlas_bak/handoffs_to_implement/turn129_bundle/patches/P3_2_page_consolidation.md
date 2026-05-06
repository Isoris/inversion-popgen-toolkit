# P3.2 — Page consolidation

## Risk: medium
## Lines changed: ~50 (rename + redirect; new content is in P5.1 + P6.1)
## Depends on: P3.1 nav must be in place
## Verification: clicking the consolidated tabs lands on the merged page

---

## What's wrong

Two duplications make the nav heavier than it should be:

1. **SV evidence as a separate page** — SV records are evidence
   for boundary location. They belong on the Boundaries page, not
   on their own.
2. **Catalogue separate from karyotype/tier** — the catalogue is
   the final candidate table: candidate ID + karyotype +
   tier + confidence + boundary status + notes. Karyotype/tier is
   already showing those things per-candidate. They want to
   become one page.

## The two consolidations

### Consolidation A: SV evidence into Boundaries

Atlas currently has `page11` (boundaries). Rename to
"Boundaries + SV evidence" semantically; the page-id stays
`page11` so existing handlers don't break.

Sub-button label in `SUB_BUTTONS_BY_MAIN.refine`:

```js
{ page: 'page11', label: 'Boundaries + SV' }
```

Inside `<div class="page" id="page11">`, the page becomes a
multi-section layout. P5.1 specifies the section list and skeleton.
This patch just renames the tab.

If there's an existing standalone SV page (e.g.
`<div class="page" id="pageSV">` or similar), absorb its content
under a new section heading inside page11, then delete the
standalone page from the DOM. Anchor: search for `'pageSV'`,
`'sv_evidence'`, `'sv_panel'`. If hits exist, list them and we'll
plan the merge in a follow-up.

### Consolidation B: Catalogue into Karyotype / Tier

Atlas currently has `page4` (karyotype/tier) and `page3`
(catalogue). Merge.

The merged page should have **3 sub-views inside it** (not nav
sub-buttons — internal toggle):

- **Karyotype** — per-promoted-candidate sample-level regime
  breakdown (HOMO_1 / HET / HOMO_2 from K=3).
- **Tier** — 14-axis classification view (existence layers,
  boundary quality, etc.).
- **Catalogue** — the sortable / filterable candidate table.

The catalogue table itself probably already lives in page3's
DOM. Move it inside page4 as a new section, keep its existing
event handlers. Then delete page3.

Sub-button label in `SUB_BUTTONS_BY_MAIN.refine`:

```js
{ page: 'page4', label: 'Karyotype / Tier / Catalogue' }
```

P6.1 specifies the merged page's internal toggle and layout.
This patch just renames the tab and moves DOM.

## DOM moves

### Page 11 absorbs SV

Anchor: search for `<div class="page" id="page11">` and
`<div class="page" id="pageSV">` (if exists).

Move pageSV's children inside page11 wrapped in
`<section class="boundary-section sv-section">`. Delete the
empty pageSV div.

If page11 doesn't have a section structure yet, add one:

```html
<div class="page" id="page11">
  <header class="page-header">
    <h2>Boundaries + SV evidence</h2>
  </header>
  <section class="boundary-section" id="boundary-coords">
    <!-- existing boundary coordinate track -->
  </section>
  <section class="boundary-section" id="boundary-zones">
    <!-- left/right statistical boundary zones -->
  </section>
  <!-- ... 8 more sections per P5.1 ... -->
</div>
```

P5.1 gives the full skeleton.

### Page 4 absorbs Catalogue

Anchor: search for `<div class="page" id="page4">` and
`<div class="page" id="page3">`.

Move page3's catalogue table into page4 inside a new section.
Delete the empty page3 div.

```html
<div class="page" id="page4">
  <header class="page-header">
    <h2>Karyotype / Tier / Catalogue</h2>
    <nav class="internal-toggle">
      <button data-view="karyotype" class="active">Karyotype</button>
      <button data-view="tier">Tier</button>
      <button data-view="catalogue">Catalogue</button>
    </nav>
  </header>
  <section data-view-content="karyotype" class="active">
    <!-- existing karyotype content -->
  </section>
  <section data-view-content="tier">
    <!-- existing tier content -->
  </section>
  <section data-view-content="catalogue">
    <!-- moved from page3 -->
  </section>
</div>
```

P6.1 specifies the internal-toggle JS.

## Reference cleanup

After the consolidation, anywhere code says
`document.getElementById('page3')` or
`document.getElementById('pageSV')` will return null. Find all
such references and update.

Common spots:
- The old tab-click handler may have `if (target === 'page3')`
  branches — remove or redirect to `page4`.
- Any URL/hash routing that mentions these page ids.

```bash
# Audit
grep -rn "'page3'\|\"page3\"\|page === 'page3'" Atlas/Inversion_atlas.html
grep -rn "'pageSV'\|\"pageSV\"" Atlas/Inversion_atlas.html
```

Replace `page3` → `page4` (with `?view=catalogue` query if you want
deep-linking). Replace `pageSV` → `page11`.

## Verification

1. Apply P3.1 first (nav must be there).
2. Apply P3.2.
3. Reload.
4. Refine row shows: "Boundaries + SV", "Karyotype / Tier / Catalogue",
   "Windows". (3 subs, not 5.)
5. Click "Karyotype / Tier / Catalogue" → page4 active. Internal
   toggle shows Karyotype | Tier | Catalogue. Click Catalogue → the
   table that used to live on page3 renders here.
6. Click "Boundaries + SV" → page11 active. SV evidence appears
   inline below the boundary panels.

## Test

```js
// tests/test_p3_2_consolidation.js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

// Page 3 (catalogue) and pageSV are gone as standalone pages
ok('page3 removed as standalone',
   !/<div class="page" id="page3">/.test(html));
ok('pageSV removed as standalone',
   !/<div class="page" id="pageSV">/.test(html));

// Page 4 has the three internal views
ok('page4 has karyotype view',
   /data-view-content="karyotype"/.test(html));
ok('page4 has tier view',
   /data-view-content="tier"/.test(html));
ok('page4 has catalogue view',
   /data-view-content="catalogue"/.test(html));

// Sub-button labels
ok('Refine sub: Boundaries + SV',
   /'page11'[^}]*'Boundaries \+ SV'/.test(html));
ok('Refine sub: Karyotype / Tier / Catalogue',
   /'page4'[^}]*'Karyotype \/ Tier \/ Catalogue'/.test(html));

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## Risk note

If `page3`'s catalogue or `pageSV`'s SV viewer have their own JS
modules that fire on `DOMContentLoaded` and look up element IDs, they
may need updating to look inside the merged page. Audit `js/` for:

```
grep -rn "getElementById('page3')\|getElementById('pageSV')" Atlas/js/
```

Each hit needs review. The fix may be as simple as updating the
selector, or as involved as scoping the module to the new section.
