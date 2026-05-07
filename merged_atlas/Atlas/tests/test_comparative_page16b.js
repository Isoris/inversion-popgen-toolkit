// tests/test_comparative_page16b.js
// Smoke test for inversion_comparative/page16b.js (multi-species cockpit).
//
// Verifies parse-clean + dynamic import + presence of the major rendering
// entry points and IO loaders for the four JSON layers this page consumes:
//   dotplot_mashmap_v1, synteny_multispecies_v1, phylo_tree_v1,
//   dxy_per_inversion_v1, comparative_te_breakpoint_fragility_v1,
//   karyotype_lineage_v1
//
// Run: node tests/test_comparative_page16b.js

import { spawnSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, resolve } from 'node:path';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '..');
const PAGE = resolve(ROOT, 'inversion_comparative/page16b.js');

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log('--- page16b.js: parse + import ---');

const chk = spawnSync('node', ['--check', PAGE], { encoding: 'utf8' });
check('node --check passes', chk.status === 0, chk.stderr.trim() || '');

let imported = false;
try { await import(PAGE); imported = true; }
catch (e) { console.log('    import threw:', e.message.split('\n')[0]); }
check('module imports without throwing', imported);

const src = readFileSync(PAGE, 'utf8');
check('file is >2000 lines', src.split('\n').length > 2000, `${src.split('\n').length} lines`);

// Core renderers
check('contains _renderMultiSpeciesPage',       src.includes('function _renderMultiSpeciesPage'));
check('contains _msRenderTreeSvg',              src.includes('function _msRenderTreeSvg'));
check('contains _msRenderActiveHeader',         src.includes('function _msRenderActiveHeader'));
check('contains _msRenderCenterColumn',         src.includes('function _msRenderCenterColumn'));
check('contains _msRenderDetailColumn',         src.includes('function _msRenderDetailColumn'));

// Auto-suggest
check('contains _msAutoSuggestArchitecture',    src.includes('function _msAutoSuggestArchitecture'));
check('contains _msAutoSuggestAgeModel',        src.includes('function _msAutoSuggestAgeModel'));

// Six layer-loaders
const loaders = [
  '_isDotplotMashmapJSON',           '_storeDotplotMashmap',
  '_isSyntenyMultispeciesJSON',      '_storeSyntenyMultispecies',
  '_isPhyloTreeJSON',                '_storePhyloTree',
  '_isDxyPerInversionJSON',          '_storeDxyPerInversion',
  '_isCompTEFragilityJSON',          '_storeCompTEFragility',
  '_isKaryotypeLineageJSON',         '_storeKaryotypeLineage',
];
for (const fn of loaders) {
  check(`contains ${fn}`, src.includes(`function ${fn}`));
}

// Constants (re-anchored prelude)
check('contains DOTPLOT_MASHMAP_TOOL const',
      /const\s+DOTPLOT_MASHMAP_TOOL\s*=/.test(src));
check('contains DOTPLOT_MASHMAP_LS_KEY const',
      /const\s+DOTPLOT_MASHMAP_LS_KEY\s*=/.test(src));

// Manifest is present
check('TODO_MISSING annotations present', /TODO_MISSING\(_esc\)/.test(src));

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail ? 1 : 0);
