# P1.2 — Server URL: canonical localStorage key with fallback

## Risk: low
## Lines changed: ~30
## Depends on: nothing
## Verification: localStorage probe shows one canonical key, atlas talks to server on 8766

---

## What's wrong

You currently set three different localStorage keys to point to the
local server:

```js
localStorage.setItem("inversion_atlas.serverUrl",  "http://127.0.0.1:8766");
localStorage.setItem("inversion_atlas.server_url", "http://127.0.0.1:8766");
localStorage.setItem("atlasServer.url",            "http://127.0.0.1:8766");
```

Three keys means three places in code reading from different keys —
or worse, racing to read whichever was set last. Symptom: setting
one key doesn't always switch the URL because something reads the
other key.

## The fix

Pick **one** canonical key and centralize all reads through one
helper. Other keys become deprecated fallbacks (read but never
written).

### Canonical

```
inversion_atlas.server_url
```

(Snake_case to match other `inversion_atlas.*` keys in your
codebase. Easy to grep for.)

### Helper function

Anchor: search for the existing `atlasServer` object (turn 127
shipped this around line 48050 in my snapshot). Add this method
or rewrite the existing URL getter:

```js
const ATLAS_SERVER_URL_KEYS = [
  'inversion_atlas.server_url',   // canonical
  'inversion_atlas.serverUrl',    // legacy camelCase, read-only
  'atlasServer.url',              // legacy alt-namespace, read-only
];
const ATLAS_SERVER_URL_DEFAULT = 'http://127.0.0.1:8766';

function _atlasGetServerUrl() {
  for (const k of ATLAS_SERVER_URL_KEYS) {
    try {
      const v = localStorage.getItem(k);
      if (v && typeof v === 'string' && v.trim()) {
        return v.trim();
      }
    } catch (_) { /* localStorage may be disabled */ }
  }
  return ATLAS_SERVER_URL_DEFAULT;
}

function _atlasSetServerUrl(url) {
  try {
    localStorage.setItem(ATLAS_SERVER_URL_KEYS[0], url);  // canonical only
    // optional: clean up legacy keys
    for (let i = 1; i < ATLAS_SERVER_URL_KEYS.length; i++) {
      try { localStorage.removeItem(ATLAS_SERVER_URL_KEYS[i]); } catch (_) {}
    }
  } catch (_) {}
}
```

### Migrate every read

Anchor: grep for `localStorage.getItem` calls related to server URL.
Replace each with `_atlasGetServerUrl()`.

Common patterns to find and replace:

```js
// before
const url = localStorage.getItem('atlasServer.url') || 'http://...';
// after
const url = _atlasGetServerUrl();
```

```js
// before
const url = localStorage.getItem('inversion_atlas.serverUrl');
// after
const url = _atlasGetServerUrl();
```

### Migrate every write

Anchor: grep for `localStorage.setItem(...server`.
Replace each with `_atlasSetServerUrl(url)`.

There may be a settings UI that lets the user paste a URL — that
input handler should also call `_atlasSetServerUrl()`.

### Default port 8766

The existing `atlasServer` object likely has a hardcoded default
URL — change it to `http://127.0.0.1:8766` to match
`popstats_server.local.yaml`.

Anchor:
```
http://127.0.0.1:8765
```
or
```
const _ATLAS_SERVER_URL_DEFAULT
```

Replace 8765 → 8766 in any default. (Confirm by grep; there's
likely just one or two hits.)

## Verification

```js
// Run in console after applying the patch
console.log('canonical:', localStorage.getItem('inversion_atlas.server_url'));
console.log('legacy 1: ', localStorage.getItem('inversion_atlas.serverUrl'));  // null after _atlasSetServerUrl
console.log('legacy 2: ', localStorage.getItem('atlasServer.url'));            // null after _atlasSetServerUrl
console.log('resolved: ', _atlasGetServerUrl());                               // 'http://127.0.0.1:8766'
```

Then click the server-status badge in the yellow header folder
(turn 127). It should ping `http://127.0.0.1:8766/health` and flip
to green `● server`.

## Test

```js
// tests/test_server_url_canonical.js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

ok('canonical key constant defined',
   /ATLAS_SERVER_URL_KEYS\s*=\s*\[[^\]]*'inversion_atlas\.server_url'/.test(html));
ok('default port 8766 (not 8765)',
   /127\.0\.0\.1:8766/.test(html) && !/127\.0\.0\.1:8765/.test(html));
ok('helper _atlasGetServerUrl defined',
   /function _atlasGetServerUrl\(\)/.test(html));
ok('helper _atlasSetServerUrl defined',
   /function _atlasSetServerUrl\(/.test(html));
ok('legacy keys still readable',
   /'inversion_atlas\.serverUrl'/.test(html) &&
   /'atlasServer\.url'/.test(html));

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## Risk note

The fallback chain reads ALL legacy keys before the canonical. If a
user has stale data in `atlasServer.url` from before, AND the
canonical key is empty, the legacy will win. To avoid this, the
order in `ATLAS_SERVER_URL_KEYS` must be **canonical first**. The
helper above does this correctly.
