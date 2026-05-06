# P1.1 — `const out` reassignment fix in `collectPopstatsTracks`

## Risk: trivial
## Lines changed: 1
## Depends on: nothing
## Verification: page loads without console error

---

## What's wrong

Browser console:
```
Uncaught TypeError: invalid assignment to const 'out'
   at collectPopstatsTracks
   Inversion_atlas.html ~line 57495
```

Pattern: code declares `const out = ...` then later reassigns
`out = ...`. JS won't let you reassign a `const`.

## Anchor strings

To find the right place in your `Inversion_atlas.html`:

```
function collectPopstatsTracks
```

If that exact name doesn't exist, try:
```
collectPopstatsTracks
```

The function is somewhere around line 57495 in the version where
this bug manifests. Look for the pattern:

```js
function collectPopstatsTracks(...) {
  ...
  const out = ...;
  ...
  out = ...;     // <-- this line is the bug
  ...
  return out;
}
```

## The fix

Pick one of two equally-correct options.

### Option A — change to `let` (smallest patch, recommended)

Change:
```js
  const out = ...;
```
to:
```js
  let out = ...;
```

That's it. One word.

### Option B — restructure to never reassign

If the reassignment is doing something like `out = out.filter(...)`
or `out = [...out, ...more]`, replace with mutating ops:

```js
  // before
  const out = [];
  for (...) { ... out = out.concat(other); }

  // after
  const out = [];
  for (...) { ... out.push(...other); }
```

Use Option A unless there's a clear reason to refactor.

## Why this matters beyond the immediate error

The thrown TypeError aborts the rest of `collectPopstatsTracks`
mid-execution. Anything that depended on its return value (popgen
panel render, popstats track list, popstatsLive wiring) silently
falls through to its empty-state branch. So the user-visible
symptom is "popstats panel won't populate" even after the data
arrives — same as missing data, but it's actually a runtime crash.

## Verification

1. Apply the fix.
2. Reload the atlas page.
3. Open browser console. No `invalid assignment to const 'out'`.
4. Open the popstats panel. If P2.1 is also applied, the panel
   populates with charts after Compute.

## Test (optional, source-level)

Add to `tests/test_const_out_fix.js` (a new tiny file):

```js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

// Find the function body
const m = html.match(/function collectPopstatsTracks[\s\S]*?\n\}/);
if (!m) {
  console.log('FAIL: collectPopstatsTracks not found');
  process.exit(1);
}
const body = m[0];

// `const out = ...` followed by `out = ...` (without `let|var|const`)
// is the bug pattern.
const decl = /\bconst\s+out\s*=/.test(body);
const reassign = /(?<!const\s|let\s|var\s)\bout\s*=(?!=)/.test(body);

if (decl && reassign) {
  console.log('FAIL: const out is still reassigned in body');
  process.exit(1);
}
console.log('PASS: const out reassignment fix in place');
```

Run with `node tests/test_const_out_fix.js`. Should print PASS.
