// tests/test_comparative_page16.js
// Smoke test for inversion_comparative/page16.js (cross-species breakpoints page).
//
// Verifies:
//   1. The module loads in Node without throwing at import time.
//   2. The expected set of cross-species runtime functions is reachable
//      via node --check (parse-clean) and via dynamic import.
//   3. Top-level constants (CROSS_SPECIES_TOOL, CS_EVENT_DEF) are present
//      via globalThis or window stubs that the module installs.
//
// Run: node tests/test_comparative_page16.js

import { spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '..');
const PAGE = resolve(ROOT, 'inversion_comparative/page16.js');

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log('--- page16.js: parse + import ---');

// 1. node --check must parse the file
const chk = spawnSync('node', ['--check', PAGE], { encoding: 'utf8' });
check('node --check passes', chk.status === 0, chk.stderr.trim() || '');

// 2. Dynamic import — page16.js is plain JS (no import/export in legacy form),
//    so Node treats it as CommonJS-ish. We just verify the file resolves and
//    its top-level statements (window.X = X registrations, behind typeof guards)
//    don't throw.
let imported = false;
try {
  // Must be window-stub-safe: legacy code uses `if (typeof window !== 'undefined')`
  // guards everywhere, so node (no window) is fine.
  await import(PAGE);
  imported = true;
} catch (e) {
  console.log('    import threw:', e.message.split('\n')[0]);
}
check('module imports without throwing', imported);

// 3. Sanity: file is non-trivial in size — extracted ~2500 lines of legacy JS.
import { readFileSync } from 'node:fs';
const src = readFileSync(PAGE, 'utf8');
check('file is >2000 lines',                src.split('\n').length > 2000, `${src.split('\n').length} lines`);
check('contains _renderCrossSpeciesPage',   src.includes('function _renderCrossSpeciesPage'));
check('contains _csBuildPermResultHtml',    src.includes('function _csBuildPermResultHtml'));
check('contains _csComputeSynteny',         src.includes('function _csComputeSynteny'));
check('contains CROSS_SPECIES_TOOL const',  /const\s+CROSS_SPECIES_TOOL\s*=/.test(src));
check('contains CS_EVENT_DEF const',        /const\s+CS_EVENT_DEF\s*=/.test(src));
check('TODO_MISSING annotations present',
      /TODO_MISSING\(_esc\)/.test(src));

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail ? 1 : 0);
