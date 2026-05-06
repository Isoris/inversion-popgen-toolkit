# P3.1 — Two-row nav: main row + sub row

## Risk: medium
## Lines changed: ~150
## Depends on: nothing (replaces turn 128 single-row pills)
## Verification: clicking a main button updates Row 2, doesn't append next to clicked button

---

## What's wrong

Turn 128 shipped pills with in-row fold/unfold. You rejected this:

> When clicking a top button like "Refine", its child buttons appeared
> to the right in the same row. This is bad.

You want **two true rows**:

```
Row 1 (#atlas_main_nav):  Overview  Discovery  Refine  Evidence  Compare  Output  Help
Row 2 (#atlas_sub_nav):   <subs of activeMain>
```

Click Refine → Row 1 unchanged, Row 2 becomes:
`Boundaries + SV evidence  |  Karyotype / Tier / Catalogue  |  Windows`

Click Discovery → Row 1 unchanged, Row 2 becomes:
`local PCA |z|  |  local PCA θπ  |  local PCA GHSL  |  candidate focus  |  negative regions`

## Architecture

### State

```js
let activeMain = 'discovery';   // current main section
let activePage = 'page1';       // current actual page (within active main)
```

`activeMain` decides which sub-buttons render.
`activePage` decides which content area is visible.

When user clicks a main button: `activeMain` updates, sub-row
re-renders, but `activePage` does NOT change (the content area stays
on its current page until the user picks a sub-button).

When user clicks a sub-button: `activePage` updates, content area
switches.

When user navigates programmatically (e.g. cs-bp jumps to page16),
both `activeMain` and `activePage` should update so the user lands
in the right main section automatically.

### Mapping

```js
const SUB_BUTTONS_BY_MAIN = {
  overview: [],   // landing — no subs, just one content area
  discovery: [
    { page: 'page1',  label: 'local PCA |z|' },
    { page: 'page12', label: 'local PCA θπ' },
    { page: 'page15', label: 'local PCA GHSL' },
    { page: 'page2',  label: 'candidate focus' },
    { page: 'page19', label: 'negative regions' },
  ],
  refine: [
    { page: 'page11', label: 'Boundaries + SV' },
    { page: 'page4',  label: 'Karyotype / Tier / Catalogue' },
    { page: 'page8',  label: 'Windows' },
  ],
  evidence: [
    { page: 'page6',  label: 'Popstats' },
    { page: 'page7',  label: 'Ancestry' },
  ],
  compare: [
    { page: 'page16',  label: 'Cross-species' },
    { page: 'page16b', label: 'Multi-species' },
  ],
  output: [
    { page: 'page9',          label: 'Confirmed' },
    { page: 'page21',         label: 'Annotation' },
    { page: 'page10',         label: 'Markers' },
    { page: 'page17',         label: 'Stats profile' },
    { page: 'page18',         label: 'Marker panel' },
    { page: 'page_overview',  label: 'Overview spreadsheet' },
  ],
  help: [],
};

// Reverse lookup: given a page, which main owns it?
const MAIN_FOR_PAGE = {};
for (const [main, subs] of Object.entries(SUB_BUTTONS_BY_MAIN)) {
  for (const s of subs) MAIN_FOR_PAGE[s.page] = main;
}
MAIN_FOR_PAGE.page5 = 'help';   // help itself
```

## Markup

Replace the entire turn 128 tab bar (with pills + tab buttons in
one row) with two separate `<div>` containers:

```html
<!-- turn 129 P3.1: two-row nav. Row 1 = main sections only. Row 2 =
     subpages of the active main section. Replaces the single-row pill
     fold/unfold from turn 128. -->
<div id="atlas_main_nav">
  <button class="main-btn" data-main="overview"   data-color="grey">Overview</button>
  <button class="main-btn" data-main="discovery"  data-color="blue"   data-active="1">Discovery</button>
  <button class="main-btn" data-main="refine"     data-color="cyan">Refine</button>
  <button class="main-btn" data-main="evidence"   data-color="green">Evidence</button>
  <button class="main-btn" data-main="compare"    data-color="orange">Compare</button>
  <button class="main-btn" data-main="output"     data-color="yellow">Output</button>
  <button class="main-btn" data-main="help"       data-color="grey">Help</button>
</div>

<div id="atlas_sub_nav">
  <!-- populated by _renderSubNav() based on activeMain -->
</div>
```

## CSS

```css
#atlas_main_nav {
  display: flex; gap: 6px; padding: 6px 22px;
  background: var(--panel); border-bottom: 1px solid var(--rule);
  align-items: center; flex-wrap: nowrap;
}
#atlas_sub_nav {
  display: flex; gap: 4px; padding: 4px 22px;
  background: var(--panel-2); border-bottom: 1px solid var(--rule);
  align-items: center; min-height: 30px;
  flex-wrap: nowrap; overflow-x: auto;
}

.main-btn {
  background: var(--panel-2); color: var(--ink-dim);
  border: 1px solid var(--rule); border-radius: 3px;
  padding: 5px 14px; font-family: var(--mono); font-size: 11.5px;
  letter-spacing: 0.04em; cursor: pointer;
  white-space: nowrap; flex-shrink: 0;
  transition: background 0.1s, color 0.1s, border-color 0.1s;
}
.main-btn:hover {
  background: var(--panel-3); color: var(--ink);
}

/* Color accents per stage. Match turn 128 palette so users don't
   relearn. */
.main-btn[data-color="blue"]   { border-left: 3px solid rgba(79, 163, 255, 0.55); }
.main-btn[data-color="cyan"]   { border-left: 3px solid rgba(122, 211, 219, 0.55); }
.main-btn[data-color="green"]  { border-left: 3px solid rgba(60, 192, 138, 0.55); }
.main-btn[data-color="orange"] { border-left: 3px solid rgba(228, 138, 60, 0.55); }
.main-btn[data-color="yellow"] { border-left: 3px solid rgba(232, 196, 76, 0.55); }
.main-btn[data-color="grey"]   { border-left: 3px solid rgba(138, 148, 163, 0.50); }

/* Active state: the section whose row 2 is showing */
.main-btn[data-active="1"]                    { color: var(--ink); }
.main-btn[data-active="1"][data-color="blue"]   { background: rgba(79, 163, 255, 0.20); border-color: rgba(79, 163, 255, 0.65); }
.main-btn[data-active="1"][data-color="cyan"]   { background: rgba(122, 211, 219, 0.20); border-color: rgba(122, 211, 219, 0.65); }
.main-btn[data-active="1"][data-color="green"]  { background: rgba(60, 192, 138, 0.20); border-color: rgba(60, 192, 138, 0.65); }
.main-btn[data-active="1"][data-color="orange"] { background: rgba(228, 138, 60, 0.22); border-color: rgba(228, 138, 60, 0.70); }
.main-btn[data-active="1"][data-color="yellow"] { background: rgba(232, 196, 76, 0.22); border-color: rgba(232, 196, 76, 0.70); }
.main-btn[data-active="1"][data-color="grey"]   { background: rgba(138, 148, 163, 0.18); border-color: rgba(138, 148, 163, 0.55); }

.sub-btn {
  background: transparent; color: var(--ink-dim);
  border: 0; border-bottom: 2px solid transparent;
  padding: 5px 14px; font-family: var(--mono); font-size: 11px;
  letter-spacing: 0.04em; cursor: pointer;
  white-space: nowrap; flex-shrink: 0;
  transition: color 0.1s, border-color 0.1s;
}
.sub-btn:hover { color: var(--ink); }
.sub-btn[data-active="1"] {
  color: var(--ink);
  border-bottom-color: var(--accent);
}

/* When activeMain has no subs, hide the sub row entirely (overview, help) */
#atlas_sub_nav:empty { display: none; }
```

## JS

```js
// turn 129 P3.1: two-row nav controller. Replaces the turn 128
// in-row pill fold/unfold. Row 1 = main section buttons; Row 2 =
// subpages of the active main section.
let _atlasActiveMain = 'discovery';
let _atlasActivePage = 'page1';

function _renderMainNav() {
  document.querySelectorAll('#atlas_main_nav .main-btn').forEach(btn => {
    btn.dataset.active = (btn.dataset.main === _atlasActiveMain) ? '1' : '0';
  });
}

function _renderSubNav() {
  const subRow = document.getElementById('atlas_sub_nav');
  if (!subRow) return;
  const subs = SUB_BUTTONS_BY_MAIN[_atlasActiveMain] || [];
  if (subs.length === 0) {
    subRow.innerHTML = '';   // CSS :empty hides the row
    return;
  }
  subRow.innerHTML = subs.map(s =>
    `<button class="sub-btn" data-page="${s.page}" data-active="${s.page === _atlasActivePage ? 1 : 0}">` +
    s.label +
    `</button>`
  ).join('');
  // Wire clicks
  subRow.querySelectorAll('.sub-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      _atlasSwitchToPage(btn.dataset.page);
    });
  });
}

function _atlasSwitchToMain(main) {
  if (!SUB_BUTTONS_BY_MAIN[main]) return;
  _atlasActiveMain = main;
  _renderMainNav();
  _renderSubNav();
  // Persist
  try { localStorage.setItem('inversion_atlas.active_main', main); } catch (_) {}
  // Note: do NOT auto-switch the page. User clicks a main → sees row 2,
  // then picks a sub. The content area stays on the current page until
  // they pick a sub.
}

function _atlasSwitchToPage(page) {
  _atlasActivePage = page;
  // Auto-update activeMain so the row reflects where we are
  const owner = MAIN_FOR_PAGE[page];
  if (owner && owner !== _atlasActiveMain) {
    _atlasActiveMain = owner;
  }
  _renderMainNav();
  _renderSubNav();
  // Switch content
  document.querySelectorAll('.page').forEach(p => p.classList.remove('active'));
  const targetEl = document.getElementById(page);
  if (targetEl) targetEl.classList.add('active');
  // Persist
  try {
    localStorage.setItem('inversion_atlas.active_main', _atlasActiveMain);
    localStorage.setItem('inversion_atlas.active_page', page);
  } catch (_) {}
}

function _atlasInitNav() {
  // Wire main-nav click handlers
  document.querySelectorAll('#atlas_main_nav .main-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      _atlasSwitchToMain(btn.dataset.main);
    });
  });
  // Restore from localStorage if available
  try {
    const m = localStorage.getItem('inversion_atlas.active_main');
    const p = localStorage.getItem('inversion_atlas.active_page');
    if (m && SUB_BUTTONS_BY_MAIN[m]) _atlasActiveMain = m;
    if (p) _atlasActivePage = p;
  } catch (_) {}
  _renderMainNav();
  _renderSubNav();
  // Initial page activation (mirrors what the previous tab handler did
  // on first load).
  if (_atlasActivePage) {
    document.querySelectorAll('.page').forEach(el => el.classList.remove('active'));
    const tgt = document.getElementById(_atlasActivePage);
    if (tgt) tgt.classList.add('active');
  }
}

// Expose for programmatic page switches (cs-bp jumps, etc.)
if (typeof window !== 'undefined') {
  window._atlasSwitchToMain = _atlasSwitchToMain;
  window._atlasSwitchToPage = _atlasSwitchToPage;
}

// Run after DOM ready (script lives at end of HTML)
_atlasInitNav();
```

## Migration from turn 128

1. **Remove turn 128 markup.** Delete the row of `tab-stage-pill`
   buttons and the `<button data-page="..." data-stage="...">` flat
   list inside `<nav id="tabBar">`. Keep the `⚙` gear button — move
   it into Row 1 as the leftmost element.

2. **Remove turn 128 CSS.** Delete:
   - `#tabBar .tab-stage-pill ...` rules
   - `#tabBar[data-active-stage="..."]` folding selectors
   - The data-stage band rules can stay; they don't hurt anything.

3. **Remove turn 128 JS.** Delete `_setActiveStage`,
   `_initTabGroupPills`. Replace with the new `_atlasInitNav` and
   helpers above.

4. **Update existing tab.click() callers.** Anywhere code calls
   `btn.click()` on a `data-page="..."` button, replace with
   `_atlasSwitchToPage('pageX')`.

   Anchors:
   - `document.querySelector('#tabBar button[data-page="page2"]').click()`
   - `document.querySelector('#tabBar button[data-page="page16"]').click()`
   - etc.

   Search for `'#tabBar button[data-page'` — there are about a dozen
   such call sites. Each becomes `_atlasSwitchToPage('pageX')`.

5. **Update turn 128 test files.** They assert on `.tab-stage-pill`
   markup. Replace with assertions on `.main-btn` / `.sub-btn`. The
   turn 128 test file should be deleted or rewritten.

## Verification

1. Reload atlas. Page 1 (local PCA |z|) is active by default.
2. Row 1 highlights "Discovery" with blue tint.
3. Row 2 shows: |z|, θπ, GHSL, candidate focus, negative regions.
4. The |z| sub-button is highlighted (active page).
5. Click "Refine". Row 1 highlights Refine (cyan). Row 2 changes to:
   Boundaries + SV, Karyotype/Tier/Catalogue, Windows.
   Content area still shows Page 1 (|z|) — user hasn't picked a sub yet.
6. Click "Boundaries + SV". Content area switches to page11.
7. Reload. Row 1 still on Refine. Row 2 still shows the three subs.
   Boundaries + SV is still highlighted.

## Test

```js
// tests/test_p3_1_two_row_nav.js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

ok('row 1 container exists',  /<div id="atlas_main_nav">/.test(html));
ok('row 2 container exists',  /<div id="atlas_sub_nav">/.test(html));
ok('main button for discovery', /class="main-btn" data-main="discovery"/.test(html));
ok('main button for refine',    /class="main-btn" data-main="refine"/.test(html));
ok('main button for evidence',  /class="main-btn" data-main="evidence"/.test(html));
ok('main button for compare',   /class="main-btn" data-main="compare"/.test(html));
ok('main button for output',    /class="main-btn" data-main="output"/.test(html));
ok('SUB_BUTTONS_BY_MAIN defined',
   /SUB_BUTTONS_BY_MAIN\s*=\s*\{/.test(html));
ok('refine has page11 (boundaries) sub',
   /refine:\s*\[[\s\S]*?page:\s*'page11'/.test(html));
ok('refine has page4 (karyotype/tier/catalogue) sub',
   /refine:\s*\[[\s\S]*?page:\s*'page4'/.test(html));
ok('compare has page16b sub',
   /compare:\s*\[[\s\S]*?page:\s*'page16b'/.test(html));
ok('_atlasSwitchToMain defined',
   /function _atlasSwitchToMain\(main\)/.test(html));
ok('_atlasSwitchToPage defined',
   /function _atlasSwitchToPage\(page\)/.test(html));
ok('turn 128 pills removed',
   !/class="tab-stage-pill"/.test(html));

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## Open question

**Where does the gear button (`⚙ globalSettingsBtn`) go?**

Currently it's inside `<nav id="tabBar">` as the leftmost child. With
the two-row nav it could:
- (a) sit at the leftmost slot of Row 1 (looks like a main button but
  isn't — opens sidebar)
- (b) sit in the header bar (turn 127 header) — fits with other
  global controls

Recommend (a) for now — minimum disruption.
