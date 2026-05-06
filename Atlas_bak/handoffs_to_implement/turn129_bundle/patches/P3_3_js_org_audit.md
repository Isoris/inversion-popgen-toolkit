# P3.3 — JS file organization audit

## Risk: low
## Lines changed: <10 in HTML script tags
## Depends on: nothing
## Verification: console shows no module-not-found errors; tests aren't loaded into production

---

## What you said

> JS file organization:
> We had confusion between production JS modules and test JS files.
>
> Decision:
> - production atlas modules belong under: Atlas/js/
> - tests belong under: Atlas/tests/
>
> If Inversion_atlas.html imports external modules, it should import
> them from js/, for example:
>
> <script src="js/atlas_group_engine.js"></script>
> <script src="js/atlas_request_layer.js"></script>
> <script src="js/atlas_renderers.js"></script>
>
> Do not import files from tests/ into production HTML.
>
> There was also a checklist saying:
> "All eight atlas-side <script> tags loaded in correct order."

## What to audit

### 1. Inventory script tags

```bash
cd Atlas
grep -nE '<script[^>]*src=' Inversion_atlas.html
```

For each `<script src="...">` tag, verify:

- [ ] Path starts with `js/` (not `./` or `../tests/`)
- [ ] Target file exists in `Atlas/js/`
- [ ] Target file is a production module, not a test
- [ ] Target file does NOT live in `tests/`
- [ ] No script reference to `server_turn1/` or
  `server_turn11c_ld_fast/` (those are server-side, can't load in
  browser anyway)
- [ ] No `<script src="">` with empty src (silent 404)

### 2. List `Atlas/js/` contents

```bash
ls -la Atlas/js/
```

For each `.js` file in `js/`, decide:
- Is it imported by `Inversion_atlas.html`?
- Or by `Diversity_atlas.html` / `Genome_atlas.html` /
  `Population_atlas.html`?
- Or is it dead code from an earlier turn?

If dead, move to `Atlas/_audits/dead_js/`.

### 3. List `Atlas/tests/` contents

```bash
ls -la Atlas/tests/
```

Tests must NOT be imported by any HTML file. Confirm:

```bash
grep -rn 'src="tests/\|src="../tests/\|src="\./tests/' Atlas/*.html
# Should be empty.
```

### 4. Dependency order

If module B references symbols from module A, A must load before B.
Common ordering for atlas modules (rough guess based on typical naming):

```html
<!-- Order: foundations → request → renderers → engines → page wires -->
<script src="js/atlas_constants.js"></script>      <!-- if exists -->
<script src="js/atlas_request_layer.js"></script>  <!-- HTTP client -->
<script src="js/atlas_renderers.js"></script>       <!-- DOM renderers -->
<script src="js/atlas_group_engine.js"></script>    <!-- group dock -->
<script src="js/atlas_dotplot.js"></script>         <!-- visualizations -->
<script src="js/atlas_focal_vs_bg.js"></script>
<!-- ... -->
<!-- Last: any page-specific wiring that depends on all of the above -->
<script src="js/atlas_init.js"></script>            <!-- if exists -->
```

The user's quoted line "All eight atlas-side `<script>` tags loaded
in correct order" implies there are exactly 8. Verify.

### 5. Check the console

After fixes:

```
F12 → Console
Reload
```

Look for:
- `Failed to load resource: ...js/<name>.js` → broken path
- `Uncaught ReferenceError: <symbol> is not defined` → wrong order

## Document the inventory

Once audited, write a short README inside `Atlas/js/`:

```md
# Atlas JS modules

Production modules loaded by Inversion_atlas.html in this order:

1. atlas_request_layer.js — HTTP client, talks to popstats_server
2. atlas_renderers.js     — generic DOM render helpers
3. atlas_group_engine.js  — group dock state + dim builder
4. atlas_dotplot.js       — multi-species dotplot widget
5. atlas_focal_vs_bg.js   — focal-vs-background panel
6. ...
7. ...
8. atlas_init.js          — wires page handlers, runs last

Tests live in ../tests/. Do not import test files into HTML.
```

This README is the single source of truth for "what's a module".

## Sample fixed `<script>` block

```html
<!-- turn 129 P3.3: production modules in js/, tests are NEVER loaded
     here. Order matters: low-level dependencies first, page wiring
     last. Audit done in patches/P3_3. -->
<script src="js/atlas_request_layer.js"></script>
<script src="js/atlas_renderers.js"></script>
<script src="js/atlas_group_engine.js"></script>
<script src="js/atlas_dotplot.js"></script>
<script src="js/atlas_focal_vs_bg.js"></script>
<script src="js/atlas_popstats_panel.js"></script>
<script src="js/atlas_dosage_heatmap.js"></script>
<script src="js/atlas_init.js"></script>
```

(Names are illustrative — match what's actually in your `js/`.)

## Test

```js
// tests/test_p3_3_script_paths.js
const fs = require('fs');
const path = require('path');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

const scripts = [...html.matchAll(/<script[^>]*src="([^"]+)"/g)].map(m => m[1]);

ok('all scripts under js/',
   scripts.every(s => s.startsWith('js/') || s.startsWith('https://')));
ok('no scripts from tests/',
   !scripts.some(s => /tests\//.test(s)));
ok('no scripts from server_turn',
   !scripts.some(s => /server_turn/.test(s)));

// Each src exists on disk
const root = path.dirname(__filename) + '/..';
for (const s of scripts) {
  if (s.startsWith('https://')) continue;
  const exists = fs.existsSync(path.join(root, s));
  ok(`script exists: ${s}`, exists);
}

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## What NOT to do

- Don't blindly add `<script>` tags for every file in `js/`. Some may
  be deliberately not loaded by `Inversion_atlas.html` (e.g.
  Diversity-only modules). Check each one before adding.
- Don't move files between `js/` and `tests/` without grepping for
  imports first.
- Don't introduce ES module syntax (`import` / `export`) without
  switching all scripts to `type="module"` — partial conversion
  breaks loading.
