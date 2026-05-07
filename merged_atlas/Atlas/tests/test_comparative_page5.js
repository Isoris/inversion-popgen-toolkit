// tests/test_comparative_page5.js
// Smoke test for inversion_comparative/page5.js (static help / quick-reference).
//
// Page5 is purely declarative content (the help tab). The JS module is a
// stub that exists so the page-loader can dispatch to it the same way it
// dispatches to the other pages (uniform routing). This test verifies:
//   1. The module parses + imports as an ES module.
//   2. Both expected exports (renderPage5, PAGE5_META) are present.
//   3. PAGE5_META declares the page is static.
//   4. The HTML fragment is non-trivial (~1158 lines of help content).
//
// Run: node tests/test_comparative_page5.js

import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '..');
const PAGE_JS  = resolve(ROOT, 'inversion_comparative/page5.js');
const PAGE_HTM = resolve(ROOT, 'inversion_comparative/page5.html');

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log('--- page5.js: parse + import ---');

const chk = spawnSync('node', ['--check', PAGE_JS], { encoding: 'utf8' });
check('node --check passes', chk.status === 0, chk.stderr.trim() || '');

let mod = null;
try { mod = await import(PAGE_JS); }
catch (e) { console.log('    import threw:', e.message.split('\n')[0]); }
check('module imports without throwing', mod !== null);

if (mod) {
  check('exports renderPage5',          typeof mod.renderPage5 === 'function');
  check('exports PAGE5_META',           typeof mod.PAGE5_META === 'object');
  if (mod.PAGE5_META) {
    check('PAGE5_META.id === page5',    mod.PAGE5_META.id === 'page5');
    check('PAGE5_META.stage === help',  mod.PAGE5_META.stage === 'help');
    check('PAGE5_META.static === true', mod.PAGE5_META.static === true);
  }
  check('renderPage5() does not throw', (() => {
    try { mod.renderPage5({}); return true; } catch { return false; }
  })());
}

// HTML fragment sanity
const html = readFileSync(PAGE_HTM, 'utf8');
check('html is >800 lines',                  html.split('\n').length > 800, `${html.split('\n').length} lines`);
check('html contains <div id="page5"',       html.includes('<div id="page5"'));
check('html contains "Quick reference"',     html.includes('Quick reference'));

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail ? 1 : 0);
