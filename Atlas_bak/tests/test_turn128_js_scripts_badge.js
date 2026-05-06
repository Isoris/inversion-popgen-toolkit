// =============================================================================
// turn 128 integration test — JS scripts loaded badge (sibling to schemaBadge)
//
// Verifies, in order:
//   1. Source-level: helpers, registry, CSS classes, click handler wired,
//      badge + overlay div present in DOM markup.
//   2. Registry shape: every entry has the required fields, every domain has
//      ≥1 script, no duplicate file paths.
//   3. _detectScriptTagsInDOM: correctly extracts file paths from a fake
//      document.scripts collection (handles plain src + ?query + ./prefix).
//   4. _computeScriptLoadState: returns 'global' / 'tag_only' / 'missing'
//      under the right conditions.
//   5. _summarizeJsRegistry: aggregates correctly across the registry.
//   6. _renderJsScriptsBadge: writes correct text and class to the badge
//      element under 'all green' / 'partial' scenarios.
//   7. CSS contract: jsScriptsBadge style rules exist alongside schemaBadge.
//   8. The literal "</script>" footgun is avoided in the registry source
//      (a regression-guard for the bug discovered + fixed during dev).
// =============================================================================

const fs = require('fs');
const path = require('path');
const vm = require('vm');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// =============================================================================
// Source-level checks
// =============================================================================
console.log('\n=== Source-level checks ===');

ok('_JS_REGISTRY constant defined',
   /const _JS_REGISTRY = \[/.test(html));
ok('_detectScriptTagsInDOM defined',
   /function _detectScriptTagsInDOM\(/.test(html));
ok('_computeScriptLoadState defined',
   /function _computeScriptLoadState\(/.test(html));
ok('_summarizeJsRegistry defined',
   /function _summarizeJsRegistry\(/.test(html));
ok('_renderJsScriptsBadge defined',
   /function _renderJsScriptsBadge\(/.test(html));
ok('_renderJsScriptsRegistryModal defined',
   /function _renderJsScriptsRegistryModal\(/.test(html));
ok('_wireJsScriptsBadge IIFE present',
   /\(function _wireJsScriptsBadge\(\) \{/.test(html));
ok('jsScriptsBadge click handler wired',
   /badge\.addEventListener\('click', _renderJsScriptsRegistryModal\)/.test(html));

// =============================================================================
// HTML markup
// =============================================================================
console.log('\n=== HTML markup ===');

ok('jsScriptsBadge span in header',
   /<span id="jsScriptsBadge"/.test(html));
ok('jsScriptsRegistryOverlay div in header',
   /<div id="jsScriptsRegistryOverlay"/.test(html));
ok('jsScriptsBadge starts hidden (display:none)',
   /<span id="jsScriptsBadge"[^>]*style="display:none/.test(html));
ok('jsScriptsBadge has cursor:pointer',
   /<span id="jsScriptsBadge"[^>]*cursor: pointer/.test(html));
ok('jsScriptsBadge sits AFTER schemaBadge in header markup',
   html.indexOf('<span id="schemaBadge"') < html.indexOf('<span id="jsScriptsBadge"'));

// =============================================================================
// CSS contract
// =============================================================================
console.log('\n=== CSS contract ===');

ok('CSS rule for #jsScriptsBadge exists',
   /header #jsScriptsBadge\b/.test(html));
ok('CSS rule for #jsScriptsBadge.v2 exists',
   /header #jsScriptsBadge\.v2/.test(html));
ok('CSS rule for #jsScriptsBadge.partial exists',
   /header #jsScriptsBadge\.partial/.test(html));
ok('jsScriptsBadge shares base style block with schemaBadge',
   /header #schemaBadge,\s*\n\s*header #jsScriptsBadge \{/.test(html));

// =============================================================================
// </script> footgun guard
// =============================================================================
console.log('\n=== </script> footgun regression guard ===');

// The registry source must NOT contain literal "</script>" inside a JS string.
// If it did, the browser would prematurely terminate the inline <script> block
// at parse time, breaking everything below. We use the canonical
// "<' + '/script>" idiom to embed a script-end tag in a JS string.
//
// Detection: walk line-by-line; flag a line if it contains "</script>"
// AND that line is clearly inside a JS string context. We treat a line as
// "inside JS context" if it does NOT start with HTML markup (i.e. is not
// the opening or closing line of a real <script src=...> tag pair, and is
// not the standalone </script> end-tag of an inline block).
const dangerLines = [];
const allLines = html.split('\n');
for (let i = 0; i < allLines.length; i++) {
  const line = allLines[i];
  if (line.indexOf('</script>') < 0) continue;
  // Real <script src=...></script> single-line: starts with <script src
  if (/^\s*<script\b[^>]*\bsrc=/.test(line)) continue;
  // Real standalone </script> end-tag for an inline block: starts with </script>
  if (/^\s*<\/script>\s*$/.test(line)) continue;
  dangerLines.push((i + 1) + ': ' + line.trim());
}
ok('no literal </script> inside JS string literals',
   dangerLines.length === 0,
   dangerLines.length ? 'lines: ' + dangerLines.slice(0, 3).join(' | ') : '');
ok('canonical "<\' + \'/script>" idiom present',
   html.indexOf("<' + '/script>") >= 0);

// =============================================================================
// Registry shape
// =============================================================================
console.log('\n=== Registry shape ===');

// Pull _JS_REGISTRY into a sandbox to inspect it as data.
function pullVarBlock(src, varName) {
  const re = new RegExp('const\\s+' + varName + '\\s*=\\s*\\[', 'g');
  const m = re.exec(src);
  if (!m) return null;
  let i = m.index + m[0].length - 1;   // position of opening '['
  let depth = 0;
  for (let j = i; j < src.length; j++) {
    const ch = src[j];
    if (ch === '[') depth++;
    else if (ch === ']') {
      depth--;
      if (depth === 0) {
        // Continue past the closing bracket and trailing semicolon
        let end = j + 1;
        while (end < src.length && /\s/.test(src[end])) end++;
        if (src[end] === ';') end++;
        return src.substring(m.index, end);
      }
    } else if (ch === '"' || ch === "'" || ch === '`') {
      const q = ch;
      j++;
      while (j < src.length) {
        if (src[j] === '\\') { j += 2; continue; }
        if (src[j] === q) break;
        j++;
      }
    } else if (ch === '/' && src[j+1] === '/') {
      while (j < src.length && src[j] !== '\n') j++;
    } else if (ch === '/' && src[j+1] === '*') {
      j += 2;
      while (j < src.length - 1 && !(src[j] === '*' && src[j+1] === '/')) j++;
      j++;
    }
  }
  return null;
}

const regBlock = pullVarBlock(html, '_JS_REGISTRY');
ok('_JS_REGISTRY block extractable', regBlock !== null);

let registry = null;
if (regBlock) {
  const sandbox = { module: { exports: {} } };
  vm.createContext(sandbox);
  // Wrap "const _JS_REGISTRY = [...];" then expose
  vm.runInContext(regBlock + '\nmodule.exports = _JS_REGISTRY;', sandbox);
  registry = sandbox.module.exports;
}

ok('_JS_REGISTRY parses as array', Array.isArray(registry));
ok('_JS_REGISTRY non-empty', registry && registry.length > 0);

if (registry) {
  let allDomainsValid = true;
  let allEntriesValid = true;
  const seen = new Set();
  for (const dom of registry) {
    if (!dom || typeof dom.domain !== 'string' || !Array.isArray(dom.scripts)) {
      allDomainsValid = false;
      continue;
    }
    if (dom.scripts.length < 1) allDomainsValid = false;
    for (const e of dom.scripts) {
      if (!e || typeof e.file !== 'string' || typeof e.global !== 'string' ||
          typeof e.desc !== 'string' || typeof e.turn !== 'string') {
        allEntriesValid = false;
        continue;
      }
      if (seen.has(e.file)) {
        allEntriesValid = false;   // duplicate file path
        console.log('    duplicate file: ' + e.file);
      }
      seen.add(e.file);
    }
  }
  ok('every domain has {domain, scripts[]}', allDomainsValid);
  ok('every script has {file, global, desc, turn} and unique file', allEntriesValid);

  // Specific known scripts must be in the registry
  const allFiles = registry.flatMap(d => d.scripts.map(e => e.file));
  ok('registry includes js/atlas_group_engine.js',
     allFiles.indexOf('js/atlas_group_engine.js') >= 0);
  ok('registry includes js/atlas_sv_evidence.js',
     allFiles.indexOf('js/atlas_sv_evidence.js') >= 0);
  ok('registry includes js/atlas_dosage_bridge.js',
     allFiles.indexOf('js/atlas_dosage_bridge.js') >= 0);
}

// =============================================================================
// Behavioural tests — pure helpers extracted into a sandbox
// =============================================================================
console.log('\n=== Behavioural tests (sandboxed) ===');

function pullFunction(src, fnName) {
  const startRegex = new RegExp(
    '^function\\s+' + fnName.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\$&') + '\\s*\\(', 'm');
  const m = src.match(startRegex);
  if (!m) return null;
  const start = m.index;
  const open = src.indexOf('{', start);
  if (open < 0) return null;
  let depth = 1, i = open + 1;
  while (i < src.length && depth > 0) {
    const ch = src[i];
    if (ch === '{') depth++;
    else if (ch === '}') depth--;
    else if (ch === '"' || ch === "'" || ch === '`') {
      const quote = ch;
      i++;
      while (i < src.length) {
        if (src[i] === '\\') { i += 2; continue; }
        if (src[i] === quote) break;
        i++;
      }
    } else if (ch === '/' && src[i+1] === '/') {
      while (i < src.length && src[i] !== '\n') i++;
    } else if (ch === '/' && src[i+1] === '*') {
      i += 2;
      while (i < src.length - 1 && !(src[i] === '*' && src[i+1] === '/')) i++;
      i++;
    }
    i++;
  }
  return src.substring(start, i);
}

const fnDetect    = pullFunction(html, '_detectScriptTagsInDOM');
const fnLoadState = pullFunction(html, '_computeScriptLoadState');
const fnSummarize = pullFunction(html, '_summarizeJsRegistry');
const fnRenderBadge = pullFunction(html, '_renderJsScriptsBadge');

ok('_detectScriptTagsInDOM extractable', !!fnDetect);
ok('_computeScriptLoadState extractable', !!fnLoadState);
ok('_summarizeJsRegistry extractable', !!fnSummarize);
ok('_renderJsScriptsBadge extractable', !!fnRenderBadge);

// Build a minimal sandbox that gives the helpers what they need.
function makeSandbox(stateOverlay) {
  const opts = stateOverlay || {};

  // Fake DOM that exposes document.scripts + getElementById('jsScriptsBadge').
  const fakeBadge = {
    textContent: '',
    className: '',
    style: { display: 'none' },
    title: '',
  };
  const fakeDoc = {
    scripts: opts.scripts || [],
    getElementById: function (id) {
      if (id === 'jsScriptsBadge') return fakeBadge;
      return null;
    },
  };
  const fakeWindow = opts.windowGlobals || {};

  const sandbox = {
    document: fakeDoc,
    window:   fakeWindow,
    console:  console,
    setTimeout, clearTimeout,
    Number, Array, Math, JSON, Object, String,
    Set, Map,
    isNaN, isFinite,
    _badge: fakeBadge,
  };
  vm.createContext(sandbox);

  // Run the registry first so the helpers can reference it.
  vm.runInContext(regBlock, sandbox);
  vm.runInContext(fnDetect,    sandbox);
  vm.runInContext(fnLoadState, sandbox);
  vm.runInContext(fnSummarize, sandbox);
  vm.runInContext(fnRenderBadge, sandbox);
  return sandbox;
}

// --- Test: _detectScriptTagsInDOM
console.log('\nTest: _detectScriptTagsInDOM');
{
  // Build a mock document.scripts with various src forms
  const mockScripts = [
    { getAttribute: (n) => n === 'src' ? 'js/atlas_group_engine.js' : null },
    { getAttribute: (n) => n === 'src' ? './js/atlas_request_layer.js' : null },
    { getAttribute: (n) => n === 'src' ? 'js/atlas_renderers_turn4.js?v=2' : null },
    { getAttribute: (n) => n === 'src' ? 'js/atlas_floating_dock.js#x' : null },
    { getAttribute: (n) => null },   // inline <script> (no src)
  ];
  const sb = makeSandbox({ scripts: mockScripts });
  vm.runInContext('var __r = _detectScriptTagsInDOM();', sb);
  const r = sb.__r;
  ok('extracts plain src',           r.has('js/atlas_group_engine.js'));
  ok('strips leading "./"',          r.has('js/atlas_request_layer.js'));
  ok('strips ?query',                r.has('js/atlas_renderers_turn4.js'));
  ok('strips #hash',                 r.has('js/atlas_floating_dock.js'));
  ok('also adds basename for tolerance', r.has('atlas_group_engine.js'));
  ok('inline <script> without src skipped',
     !Array.from(r).some(v => v === '' || v == null));
}

// --- Test: _computeScriptLoadState
console.log('\nTest: _computeScriptLoadState');
{
  const sb = makeSandbox({});
  // 'global' → tag present AND global registered
  vm.runInContext(
    'var __r = _computeScriptLoadState({file:"js/foo.js", global:"FooApi"}, ' +
    'new Set(["js/foo.js"]), function(name){ return name === "FooApi"; });', sb);
  ok('global path: tag + global → "global"', sb.__r === 'global');

  vm.runInContext(
    'var __r = _computeScriptLoadState({file:"js/foo.js", global:"FooApi"}, ' +
    'new Set(["js/foo.js"]), function(name){ return false; });', sb);
  ok('tag_only path: tag without global → "tag_only"', sb.__r === 'tag_only');

  vm.runInContext(
    'var __r = _computeScriptLoadState({file:"js/foo.js", global:"FooApi"}, ' +
    'new Set(), function(name){ return false; });', sb);
  ok('missing path: no tag → "missing"', sb.__r === 'missing');

  // Falls back to basename match if path-form is missing
  vm.runInContext(
    'var __r = _computeScriptLoadState({file:"js/foo.js", global:"FooApi"}, ' +
    'new Set(["foo.js"]), function(name){ return false; });', sb);
  ok('basename-only tag still detected as "tag_only"', sb.__r === 'tag_only');

  // Global present but no tag is still "global" (the global is what matters)
  vm.runInContext(
    'var __r = _computeScriptLoadState({file:"js/foo.js", global:"FooApi"}, ' +
    'new Set(), function(name){ return name === "FooApi"; });', sb);
  ok('global wins over missing tag → "global"', sb.__r === 'global');
}

// --- Test: _summarizeJsRegistry
console.log('\nTest: _summarizeJsRegistry');
{
  // All globals registered
  const sb = makeSandbox({});
  vm.runInContext(
    'var __r = _summarizeJsRegistry(new Set(), function(name){ return true; });', sb);
  ok('all-global: nGlobal === nTotal',
     sb.__r.nGlobal === sb.__r.nTotal && sb.__r.nGlobal > 0,
     'got ' + JSON.stringify(sb.__r));
  ok('all-global: nMissing === 0', sb.__r.nMissing === 0);
  ok('all-global: nTagOnly === 0', sb.__r.nTagOnly === 0);

  // No globals, no tags
  vm.runInContext(
    'var __r = _summarizeJsRegistry(new Set(), function(name){ return false; });', sb);
  ok('all-missing: nGlobal === 0',
     sb.__r.nGlobal === 0 && sb.__r.nMissing === sb.__r.nTotal);

  // Tag present for one specific entry, no globals
  vm.runInContext(
    'var __r = _summarizeJsRegistry(new Set(["js/atlas_group_engine.js"]), ' +
    'function(name){ return false; });', sb);
  ok('one-tag-only: nTagOnly === 1', sb.__r.nTagOnly === 1);
  ok('one-tag-only: nMissing === nTotal - 1',
     sb.__r.nMissing === sb.__r.nTotal - 1);
}

// --- Test: _renderJsScriptsBadge
console.log('\nTest: _renderJsScriptsBadge — happy path');
{
  // Pretend every expected global is registered.
  const wins = {};
  // Pull the registry to know what globals to register
  for (const dom of registry) {
    for (const e of dom.scripts) wins[e.global] = {};
  }
  const mockScripts = registry.flatMap(d => d.scripts).map(e => ({
    getAttribute: (n) => n === 'src' ? e.file : null,
  }));
  const sb = makeSandbox({ scripts: mockScripts, windowGlobals: wins });
  vm.runInContext('_renderJsScriptsBadge();', sb);
  ok('badge text matches "JS · N scripts"',
     /^JS · \d+ scripts?$/.test(sb._badge.textContent),
     'got "' + sb._badge.textContent + '"');
  ok('badge className === "v2" when all green',
     sb._badge.className === 'v2',
     'got "' + sb._badge.className + '"');
  ok('badge made visible',
     sb._badge.style.display === '');
  ok('badge title includes "expected globals registered"',
     sb._badge.title.indexOf('expected globals registered') >= 0);
}

console.log('\nTest: _renderJsScriptsBadge — partial');
{
  // Only register one global; that should land us in "partial".
  const mockScripts = registry.flatMap(d => d.scripts).map(e => ({
    getAttribute: (n) => n === 'src' ? e.file : null,
  }));
  const wins = { popgen: {} };   // just one
  const sb = makeSandbox({ scripts: mockScripts, windowGlobals: wins });
  vm.runInContext('_renderJsScriptsBadge();', sb);
  ok('badge className === "partial" when not all green',
     sb._badge.className === 'partial',
     'got "' + sb._badge.className + '"');
  ok('badge text reports just the registered count (1)',
     /^JS · 1 script$/.test(sb._badge.textContent),
     'got "' + sb._badge.textContent + '"');
}

console.log('\nTest: _renderJsScriptsBadge — singular vs plural');
{
  // 1 global → "JS · 1 script"; 0 → "JS · 0 scripts"
  const sb = makeSandbox({ scripts: [], windowGlobals: {} });
  vm.runInContext('_renderJsScriptsBadge();', sb);
  ok('zero: plural form',
     sb._badge.textContent === 'JS · 0 scripts',
     'got "' + sb._badge.textContent + '"');
}

// =============================================================================
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail > 0 ? 1 : 0);
