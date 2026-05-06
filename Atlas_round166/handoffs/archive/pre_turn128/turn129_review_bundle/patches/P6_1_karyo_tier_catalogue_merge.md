# P6.1 — Karyotype / Tier / Catalogue merged page

## Risk: medium
## Lines changed: ~150
## Depends on: P3.2 (page3 must be marked for absorption)
## Verification: page4 has 3-way internal toggle; catalogue table shows up

---

## What you said

> Catalogue should be combined with Karyotype / Tier page.
> Catalogue is basically the final candidate table:
>   candidate ID, karyotype structure, tier, confidence, boundary
>   status, notes.

## Page structure

`page4` (currently "Karyotype / Tier") absorbs `page3`'s catalogue.
Inside page4, three internal views — they're an in-page toggle, NOT
sub-buttons in row 2 (which would defeat the purpose of consolidating).

```
┌─────────────────────────────────────────────────────────────────┐
│  Karyotype / Tier / Catalogue                                   │
│  ┌─────────────┬───────────┬──────────────┐                     │
│  │ Karyotype   │ Tier      │ Catalogue    │  ← internal toggle  │
│  └─────────────┴───────────┴──────────────┘                     │
│                                                                 │
│  <active view's content>                                        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

The Catalogue view IS the master candidate table:

| candidate_id | chrom | start_bp | end_bp | karyotype | tier | confidence | boundary_status | notes |
|---|---|---|---|---|---|---|---|---|
| LG28-cand-001 | C_gar_LG28 | 14.8M | 17.7M | H1/H1: 61, H1/H2: 103, H2/H2: 62 | A1 | high | resolved | shelf |

## DOM

```html
<!-- turn 129 P6.1: page4 absorbs page3 (catalogue) and exposes a
     3-way internal view toggle. The toggle is in-page, not in row 2,
     because adding 3 more sub-buttons defeats the consolidation. -->
<div class="page" id="page4">

  <header class="page-header">
    <h2>Karyotype / Tier / Catalogue</h2>
    <nav class="internal-toggle" id="page4-view-toggle">
      <button data-view="karyotype" class="active">Karyotype</button>
      <button data-view="tier">Tier</button>
      <button data-view="catalogue">Catalogue</button>
    </nav>
  </header>

  <section data-view-content="karyotype" class="active" id="page4-karyotype">
    <!-- existing karyotype content from old page4 lives here -->
  </section>

  <section data-view-content="tier" id="page4-tier">
    <!-- existing tier (14-axis classification) content from old page4 -->
  </section>

  <section data-view-content="catalogue" id="page4-catalogue">
    <!-- moved from old page3 -->
    <div class="catalogue-toolbar">
      <input type="search" id="cat-filter" placeholder="filter candidates...">
      <button id="cat-export-tsv">⬇ TSV</button>
      <button id="cat-export-md">⬇ Markdown</button>
    </div>
    <div class="catalogue-table-wrap">
      <table id="cat-table">
        <thead>
          <tr>
            <th data-sort="id">candidate_id</th>
            <th data-sort="chrom">chrom</th>
            <th data-sort="start">start_bp</th>
            <th data-sort="end">end_bp</th>
            <th data-sort="width">width</th>
            <th data-sort="karyotype">karyotype</th>
            <th data-sort="tier">tier</th>
            <th data-sort="confidence">confidence</th>
            <th data-sort="boundary_status">boundary_status</th>
            <th>notes</th>
          </tr>
        </thead>
        <tbody id="cat-tbody">
          <!-- rows rendered by _renderCatalogueRows() -->
        </tbody>
      </table>
    </div>
    <div id="cat-empty" class="empty-hint">
      No candidates yet. Promote candidates on page 1 (local PCA |z|)
      or page 3 (candidate focus).
    </div>
  </section>

</div>
```

## CSS

```css
.internal-toggle {
  display: flex; gap: 0;
  border-bottom: 1px solid var(--rule);
  margin-top: 6px;
}
.internal-toggle button {
  background: transparent; color: var(--ink-dim);
  border: 0; padding: 6px 18px;
  border-bottom: 2px solid transparent;
  font-family: var(--mono); font-size: 11.5px;
  cursor: pointer;
  transition: color 0.1s, border-color 0.1s;
}
.internal-toggle button:hover { color: var(--ink); }
.internal-toggle button.active {
  color: var(--ink);
  border-bottom-color: var(--accent);
}

[data-view-content] { display: none; }
[data-view-content].active { display: block; }

#page4 {
  display: flex; flex-direction: column;
}

.catalogue-toolbar {
  padding: 8px 22px; display: flex; gap: 8px; align-items: center;
  border-bottom: 1px solid var(--rule);
}
.catalogue-toolbar input[type="search"] {
  flex: 1; max-width: 320px;
  background: var(--panel-2); border: 1px solid var(--rule);
  color: var(--ink); padding: 4px 8px;
  font-family: var(--mono); font-size: 11px;
  border-radius: 3px;
}
.catalogue-toolbar button {
  background: var(--panel-2); color: var(--ink);
  border: 1px solid var(--rule); border-radius: 3px;
  padding: 4px 10px; font-family: var(--mono); font-size: 10.5px;
  cursor: pointer;
}
.catalogue-toolbar button:hover {
  background: var(--panel-3); border-color: var(--accent);
}

.catalogue-table-wrap {
  overflow: auto; padding: 0 22px 16px;
}
#cat-table {
  width: 100%; border-collapse: collapse;
  font-family: var(--mono); font-size: 11px;
}
#cat-table thead th {
  background: var(--panel-2);
  border-bottom: 1px solid var(--rule);
  padding: 6px 10px; text-align: left;
  font-weight: 500; cursor: pointer; user-select: none;
}
#cat-table thead th:hover { background: var(--panel-3); }
#cat-table tbody td {
  padding: 4px 10px;
  border-bottom: 1px solid rgba(255,255,255,0.04);
}
#cat-table tbody tr:hover { background: var(--panel-2); }
```

## JS

```js
// turn 129 P6.1: internal-toggle controller for page4
function _atlasInitPage4Toggle() {
  const toggle = document.getElementById('page4-view-toggle');
  if (!toggle) return;
  toggle.addEventListener('click', e => {
    const btn = e.target.closest('button[data-view]');
    if (!btn) return;
    const view = btn.dataset.view;
    // Activate the clicked button
    toggle.querySelectorAll('button').forEach(b => {
      b.classList.toggle('active', b.dataset.view === view);
    });
    // Activate the matching content section
    document.querySelectorAll('#page4 [data-view-content]').forEach(sec => {
      sec.classList.toggle('active', sec.dataset.viewContent === view);
    });
    // Persist
    try { localStorage.setItem('inversion_atlas.page4_view', view); } catch (_) {}
  });
  // Restore last view
  try {
    const last = localStorage.getItem('inversion_atlas.page4_view');
    if (last) {
      const btn = toggle.querySelector(`button[data-view="${last}"]`);
      if (btn) btn.click();
    }
  } catch (_) {}
}
_atlasInitPage4Toggle();

// Render the catalogue table from state.candidateList
function _renderCatalogueRows() {
  const tbody = document.getElementById('cat-tbody');
  const empty = document.getElementById('cat-empty');
  if (!tbody) return;
  const list = state.candidateList || [];
  if (list.length === 0) {
    tbody.innerHTML = '';
    if (empty) empty.style.display = '';
    return;
  }
  if (empty) empty.style.display = 'none';

  // Filter
  const filter = (document.getElementById('cat-filter') &&
                  document.getElementById('cat-filter').value || '').toLowerCase();

  const rows = list
    .filter(c => !filter || JSON.stringify(c).toLowerCase().includes(filter))
    .map(c => {
      const w = (c.end_bp - c.start_bp);
      const wMb = (w / 1e6).toFixed(2) + ' Mb';
      const ll = c.locked_labels || [];
      const counts = [0, 0, 0, 0]; // h1h1, h1h2, h2h2, ambig
      for (const r of ll) {
        if (r === 0) counts[0]++;
        else if (r === 1) counts[1]++;
        else if (r === 2) counts[2]++;
        else counts[3]++;
      }
      const karyo = `H1/H1: ${counts[0]}, H1/H2: ${counts[1]}, H2/H2: ${counts[2]}` +
                    (counts[3] ? `, ambig: ${counts[3]}` : '');
      const tier = c.tier || '—';
      const conf = c.confidence != null ? c.confidence : '—';
      const bs = c.boundary_status || '—';
      const notes = c.notes || '';
      return `<tr data-cand-id="${c.id}">
        <td>${_esc(c.id || '—')}</td>
        <td>${_esc(c.chrom || state.data.chrom || '—')}</td>
        <td>${c.start_bp != null ? c.start_bp.toLocaleString() : '—'}</td>
        <td>${c.end_bp != null ? c.end_bp.toLocaleString() : '—'}</td>
        <td>${wMb}</td>
        <td>${_esc(karyo)}</td>
        <td>${_esc(tier)}</td>
        <td>${_esc(String(conf))}</td>
        <td>${_esc(bs)}</td>
        <td>${_esc(notes)}</td>
      </tr>`;
    })
    .join('');
  tbody.innerHTML = rows;

  // Click row to activate that candidate
  tbody.querySelectorAll('tr').forEach(tr => {
    tr.addEventListener('click', () => {
      const cid = tr.dataset.candId;
      const c = state.candidateList.find(x => x.id === cid);
      if (c && typeof activateCandidate === 'function') activateCandidate(c);
    });
  });
}

// Filter live
const _catFilterInput = document.getElementById('cat-filter');
if (_catFilterInput) {
  _catFilterInput.addEventListener('input',
    () => _renderCatalogueRows());
}

// Export buttons
const _catExportTsv = document.getElementById('cat-export-tsv');
if (_catExportTsv) {
  _catExportTsv.addEventListener('click', () => _exportCatalogueTSV());
}
const _catExportMd = document.getElementById('cat-export-md');
if (_catExportMd) {
  _catExportMd.addEventListener('click', () => _exportCatalogueMarkdown());
}

// Hook: re-render the catalogue whenever the candidate list changes.
// Anchor: find existing function that fires on candidate list mutations
// (e.g. _emitCandidateListChanged or similar) and call _renderCatalogueRows
// from there. If no such hook exists, call after every state mutation
// you know about.
```

## Migration

The old `page3` JS (catalogue rendering, sort, filter, export) gets
**moved** here, not rewritten. Anchor: search the atlas for
`'cat-tbody'` / `'catalogue-table'` / `'#page3'` references and
re-target them at the new IDs.

If the existing renderer uses different IDs, choose one approach:
- (a) Rename the old IDs in JS to match the new DOM (recommended).
- (b) Rename the new DOM IDs to match the old JS (if existing JS is
  large and you don't want to risk it).

Pick (a). The IDs above are reasonable and explicit.

## Verification

1. P3.2 must be applied (page3 absorbed → page4).
2. Apply P6.1.
3. Reload. Click Refine → Karyotype / Tier / Catalogue.
4. Internal toggle shows: Karyotype | Tier | Catalogue.
5. Karyotype is active by default.
6. Click Catalogue. Table renders with all candidates.
7. Type in filter box → table filters live.
8. Click "⬇ TSV" → download a .tsv with all rows.
9. Reload. The toggle remembers the last view (Catalogue, in this case).

## Test

```js
// tests/test_p6_1_karyo_tier_catalogue.js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

ok('page4 has internal toggle nav',
   /<nav class="internal-toggle" id="page4-view-toggle">/.test(html));
ok('toggle has karyotype button',
   /data-view="karyotype"[^>]*class="active"/.test(html));
ok('toggle has tier button',
   /data-view="tier"/.test(html));
ok('toggle has catalogue button',
   /data-view="catalogue"/.test(html));
ok('section content for karyotype',
   /data-view-content="karyotype"[^>]*class="active"/.test(html));
ok('section content for tier',
   /data-view-content="tier"/.test(html));
ok('section content for catalogue',
   /data-view-content="catalogue"/.test(html));
ok('catalogue table present',
   /<table id="cat-table">/.test(html));
ok('catalogue toolbar with filter',
   /<input[^>]*id="cat-filter"/.test(html));
ok('toggle controller defined',
   /function _atlasInitPage4Toggle\(\)/.test(html));
ok('catalogue render fn defined',
   /function _renderCatalogueRows\(\)/.test(html));
ok('view persistence to localStorage',
   /'inversion_atlas\.page4_view'/.test(html));

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```
