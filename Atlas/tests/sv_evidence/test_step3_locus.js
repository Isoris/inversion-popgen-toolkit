// =============================================================================
// tests/sv_evidence/test_step3_locus.js
// =============================================================================
// Step 3 regression test for atlas_sv_evidence.js. Verifies:
//   - LOCUS_TRACKS constant exposed and ordered correctly.
//   - _drawTrackZones produces 5 colour-coded segments matching the layer.
//   - _drawTrackAxis emits one tick per 500 kb across the visible window.
//   - _drawTrackSVCalls emits one glyph per SV call, in the right lane,
//     with FDR-coloured halo and click-target data attributes.
//   - Filtered-out SVs render at low opacity (0.18); kept SVs at 1.0.
//   - Highlighted SV gets a stroke ring.
//   - Placeholder tracks for missing state-bound layers (dosage / pca /
//     genes / repeats) render with appropriate hint text.
//
// Run:
//   cd Atlas && node tests/sv_evidence/test_step3_locus.js
// =============================================================================

'use strict';

const path = require('path');
const fs   = require('fs');

const MODULE_PATH  = path.resolve(__dirname, '../../js/atlas_sv_evidence.js');
const FIXTURE_PATH = path.resolve(__dirname, 'fixture_sv_genotype_counts_v1.json');

let pass = 0, fail = 0;
function assert(cond, msg) {
  if (cond) { pass++; console.log('  ✓ ' + msg); }
  else      { fail++; console.error('  ✗ ' + msg); }
}
function section(name) { console.log('\n— ' + name + ' —'); }

const mod = require(MODULE_PATH);
const C   = mod._const;
const h   = mod._internals;
const layer = JSON.parse(fs.readFileSync(FIXTURE_PATH, 'utf-8'));

// Build the bp → x mapper the renderer would use
const z = layer.zone_definitions_bp;
const lo = z.left_flank[0];
const hi = z.right_flank[1];
const drawW = 800;
const x = bp => ((bp - lo) / (hi - lo)) * drawW;

// ---------------------------------------------------------------------------
section('LOCUS_TRACKS constant');
assert(Array.isArray(C.LOCUS_TRACKS) && C.LOCUS_TRACKS.length === 7,
  'LOCUS_TRACKS has 7 entries (per spec §3.4 7-track stack)');

const trackIds = C.LOCUS_TRACKS.map(t => t.id);
const expectedOrder = ['zones','axis','dosage','pca','sv_calls','genes','repeats'];
expectedOrder.forEach((id, i) => {
  assert(trackIds[i] === id, 'track #' + i + ' is ' + id);
});

// Self-contained vs state-bound classification
const selfTracks  = C.LOCUS_TRACKS.filter(t => t.kind === 'self').map(t => t.id);
const stateTracks = C.LOCUS_TRACKS.filter(t => t.kind === 'state').map(t => t.id);
assert(selfTracks.indexOf('zones')    !== -1, 'zones is self-contained');
assert(selfTracks.indexOf('axis')     !== -1, 'axis is self-contained');
assert(selfTracks.indexOf('sv_calls') !== -1, 'sv_calls is self-contained');
assert(stateTracks.indexOf('dosage') !== -1, 'dosage requires atlas state');
assert(stateTracks.indexOf('pca')    !== -1, 'pca requires atlas state');
assert(stateTracks.indexOf('genes')  !== -1, 'genes requires atlas state');
assert(stateTracks.indexOf('repeats')!== -1, 'repeats requires atlas state');

// ---------------------------------------------------------------------------
section('_drawTrackZones — 5 colour-coded segments');
const zonesSvg = h._drawTrackZones(layer, lo, hi, drawW, 22, x);
assert(typeof zonesSvg === 'string' && zonesSvg.startsWith('<svg'),
  'returns a <svg> string');
const rectMatches = zonesSvg.match(/<rect /g) || [];
assert(rectMatches.length === 5, '5 <rect> elements (one per zone)');

// Spot-check each zone colour appears in the SVG
['#5a6472','#4fa3ff','#f5a524','#bf616a'].forEach(col => {
  assert(zonesSvg.indexOf(col) !== -1, 'zone colour present: ' + col);
});

// Boundary labels should appear above the boundary segments
['Left boundary','Right boundary','Inversion body'].forEach(lbl => {
  assert(zonesSvg.indexOf(lbl) !== -1, 'zone label rendered: "' + lbl + '"');
});

// Label colour should be the zone fill (text uses fill="..."  )
assert(/<text[^>]*fill="#4fa3ff"[^>]*>Left boundary<\/text>/.test(zonesSvg),
  'Left boundary label inherits its zone fill colour');

// ---------------------------------------------------------------------------
section('_drawTrackAxis — 500 kb ticks');
const axisSvg = h._drawTrackAxis(layer, lo, hi, drawW, 22, x);
assert(axisSvg.startsWith('<svg'), 'axis returns SVG');

// Window is 14.14 → 19.12 Mb → ~10 tick marks at 500 kb spacing.
const tickLines  = (axisSvg.match(/<line [^>]*y2="\d+(\.\d+)?"/g) || []).length;
const tickLabels = (axisSvg.match(/<text [^>]*>\d+\.\d+ Mb<\/text>/g) || []).length;
assert(tickLines >= 10,
  'axis emits ≥ 10 tick line elements (got ' + tickLines + ')');
assert(tickLabels >= 9,
  'axis emits ≥ 9 tick label texts (got ' + tickLabels + ')');

// Spot-check a couple of expected labels
assert(axisSvg.indexOf('15.00 Mb') !== -1, 'axis labels include 15.00 Mb');
assert(axisSvg.indexOf('18.00 Mb') !== -1, 'axis labels include 18.00 Mb');

// ---------------------------------------------------------------------------
section('_drawTrackSVCalls — glyph per call');
// Stub the global window/state for the renderer. _drawTrackSVCalls reads
// _state.filters and _state.highlightedSvId from module scope; modify
// _state directly via the public _state ref.
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
mod._state.highlightedSvId = null;

const callsSvg = h._drawTrackSVCalls(layer, lo, hi, drawW, 60, x);
assert(callsSvg.startsWith('<svg'), 'calls track returns SVG');

const glyphMatches = callsSvg.match(/data-sv-glyph="(SV\d+)"/g) || [];
assert(glyphMatches.length === layer.sv_calls.length,
  'one data-sv-glyph element per SV call (got ' + glyphMatches.length +
  ', expected ' + layer.sv_calls.length + ')');

// Each call must have a polygon (the triangle glyph)
const polygonMatches = callsSvg.match(/<polygon /g) || [];
assert(polygonMatches.length === layer.sv_calls.length,
  'one <polygon> per SV call');

// Each call must have an FDR halo <circle>
const circleMatches = callsSvg.match(/<circle /g) || [];
assert(circleMatches.length === layer.sv_calls.length,
  'one FDR-halo <circle> per SV call');

// Lane assignment: BND glyphs should appear before INV before DEL etc.
// We can't easily test y coordinates from the regex; instead sample SV001
// (BND) which should be in lane 0 (laneH = 12, yMid ≈ 6).
const sv001Match = callsSvg.match(/data-sv-glyph="SV001"[\s\S]*?<polygon points="([^"]+)"/);
assert(sv001Match, 'SV001 glyph is renderable');
if (sv001Match) {
  const yCoords = sv001Match[1].split(' ').map(p => parseFloat(p.split(',')[1]));
  const yMid = (Math.min(...yCoords) + Math.max(...yCoords)) / 2;
  assert(yMid < 60 / 5, 'SV001 (BND) sits in lane 0 (top, y < ' + (60/5) + ')');
}
const sv006Match = callsSvg.match(/data-sv-glyph="SV006"[\s\S]*?<polygon points="([^"]+)"/);
if (sv006Match) {
  // SV006 is DUP → lane 3 (out of 5)
  const yCoords = sv006Match[1].split(' ').map(p => parseFloat(p.split(',')[1]));
  const yMid = (Math.min(...yCoords) + Math.max(...yCoords)) / 2;
  assert(yMid > 60 * 2/5 && yMid < 60 * 4/5,
    'SV006 (DUP) sits in lane 3 (yMid in [24, 48], got ' + yMid.toFixed(1) + ')');
}

// FDR-halo colour: SV001 has fdr=2.1e-9 → red #bf616a
assert(/data-sv-glyph="SV001"[\s\S]*?<circle [^>]*fill="#bf616a"/.test(callsSvg),
  'SV001 (FDR < 1e-6) has red halo');
// SV006 has fdr=0.18 → grey #888888
assert(/data-sv-glyph="SV006"[\s\S]*?<circle [^>]*fill="#888888"/.test(callsSvg),
  'SV006 (FDR > 0.05) has grey halo');

// Filtered-out SVs (with show_only_associated=true) render at low opacity
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: true,
                       fdr_threshold: 0.05 };
const callsSvg2 = h._drawTrackSVCalls(layer, lo, hi, drawW, 60, x);
// SV007 has fdr=0.49 → filtered out at threshold 0.05
const sv007 = callsSvg2.match(/<g class="sv-locus-glyph" data-sv-glyph="SV007"[^>]*opacity="([^"]+)"/);
assert(sv007 && parseFloat(sv007[1]) < 0.5,
  'SV007 (FDR > 0.05) is dimmed when show_only_associated is on (opacity ' +
  (sv007 ? sv007[1] : '?') + ')');
// SV001 has fdr=2.1e-9 → kept
const sv001 = callsSvg2.match(/<g class="sv-locus-glyph" data-sv-glyph="SV001"[^>]*opacity="([^"]+)"/);
assert(sv001 && parseFloat(sv001[1]) > 0.9,
  'SV001 (FDR ≪ 0.05) is full opacity when show_only_associated is on');

// Highlight ring: when _state.highlightedSvId is set, that glyph gets a stroke
mod._state.highlightedSvId = 'SV004';
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
const callsSvg3 = h._drawTrackSVCalls(layer, lo, hi, drawW, 60, x);
assert(/data-sv-glyph="SV004"[\s\S]*?<polygon [^>]*stroke="var\(--accent-2/.test(callsSvg3),
  'SV004 polygon gets accent-2 stroke when highlighted');

// Cleanup module state
mod._state.filters = { ...C.DEFAULT_FILTERS };
mod._state.highlightedSvId = null;

// ---------------------------------------------------------------------------
section('_drawTrackPlaceholder — labelled empty state');
['dosage','pca','genes','repeats'].forEach(id => {
  const t = C.LOCUS_TRACKS.find(x => x.id === id);
  const svg = h._drawTrackPlaceholder(t, drawW, t.h);
  assert(svg.startsWith('<svg'), id + ' placeholder returns SVG');
  assert(svg.indexOf('font-style="italic"') !== -1,
    id + ' placeholder text is italicised');
  assert(svg.indexOf('Drop') !== -1 || svg.indexOf('Load') !== -1,
    id + ' placeholder includes a load hint');
});

// ---------------------------------------------------------------------------
section('edge cases');
// Empty sv_calls: track renders a hint text instead of throwing
const emptyLayer = { ...layer, sv_calls: [] };
const emptySvg = h._drawTrackSVCalls(emptyLayer, lo, hi, drawW, 60, x);
assert(emptySvg.indexOf('no SV calls in layer') !== -1,
  'empty sv_calls renders "no SV calls in layer" hint');

// Unknown SV type falls into the Other lane
const weirdLayer = { ...layer,
  sv_calls: [{ ...layer.sv_calls[0], sv_id: 'SVXX', sv_type: 'WEIRD',
               position_bp: 16000000,
               fisher: { fdr_bh: 0.5, p_value: 0.5, odds_ratio: 1, comparison: '' } }] };
mod._state.filters = { ...C.DEFAULT_FILTERS, zone: 'all', show_only_associated: false };
const weirdSvg = h._drawTrackSVCalls(weirdLayer, lo, hi, drawW, 60, x);
assert(weirdSvg.indexOf('data-sv-glyph="SVXX"') !== -1,
  'unknown SV type still renders (falls into Other lane)');

// Off-domain calls are culled
const farLayer = { ...layer,
  sv_calls: [{ ...layer.sv_calls[0], sv_id: 'SVFAR', position_bp: 99999999,
               fisher: { fdr_bh: 0.5, p_value: 0.5, odds_ratio: 1, comparison: '' } }] };
const farSvg = h._drawTrackSVCalls(farLayer, lo, hi, drawW, 60, x);
assert(farSvg.indexOf('SVFAR') === -1,
  'off-domain SV is culled from rendering');

// ---------------------------------------------------------------------------
console.log('\n=========================================');
console.log('  PASS: ' + pass + '   FAIL: ' + fail);
console.log('=========================================');
process.exit(fail === 0 ? 0 : 1);
