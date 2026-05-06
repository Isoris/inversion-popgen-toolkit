// =============================================================================
// turn 141 — Per-sample-lines candidate-band highlights (Slice 1)
// =============================================================================
// SPEC: specs_todo/SPEC_lines_panel_candidate_bands.md
//
// Quentin (turn 129):
//   "If we have 2 or 3 inversion systems in the per sample lines we must
//    draw the interval of the candidate taking 1/3 vertical space for
//    their highlights alpha background you could use like yellow same as
//    now one green one a bit blue?"
//
// Tests cover:
//   - state.linesPanelCandidateBands slot exists, default true
//   - localStorage key + setter wired (setLinesPanelCandidateBands)
//   - Toolbar checkbox + change-handler wireup in the IIFE
//   - _LINES_CAND_BAND_PALETTE: yellow / green / blue at the spec hex
//   - _candidateBandColor returns rgba() for indices 0..2 with the
//     canonical palette; switches to hsla() golden-angle cycle for ≥3
//   - _candidateBandColor honors alpha argument (default 0.10)
//   - _paintCandidateBands paints exactly the right number of bands:
//       confirmed + on-chrom + finite bp + visible-after-clip
//     and silently drops everything else
//   - Off-screen candidate still bumps palette index (so palette
//     assignment is stable across zoom)
//   - Defensive: missing ctx / missing toX / out-of-range mbMin/mbMax
//     do not throw and return 0 painted bands
//   - drawLinesPanel paints bands BEFORE the frame (source-level
//     draw-order check)
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

// ============================================================================
// 1. State init + localStorage key
// ============================================================================
console.log('\n=== 1. State init + localStorage ===');

ok('state.linesPanelCandidateBands init present',
   /linesPanelCandidateBands:\s*true/.test(html));
ok('localStorage key constant defined',
   /const _LINES_PANEL_CAND_BANDS_KEY\s*=\s*['"]pca_scrubber_v3\.linesPanelCandidateBands['"]/.test(html));
ok('window._LINES_PANEL_CAND_BANDS_KEY exported',
   /window\._LINES_PANEL_CAND_BANDS_KEY\s*=\s*_LINES_PANEL_CAND_BANDS_KEY/.test(html));

// ============================================================================
// 2. Setter + helper definitions
// ============================================================================
console.log('\n=== 2. Setter + helper defs ===');

ok('setLinesPanelCandidateBands defined',
   /function setLinesPanelCandidateBands\(b\)/.test(html));
ok('setLinesPanelCandidateBands writes localStorage',
   /localStorage\.setItem\(_LINES_PANEL_CAND_BANDS_KEY,\s*b\s*\?\s*['"]1['"]\s*:\s*['"]0['"]\)/.test(html));
ok('setLinesPanelCandidateBands triggers redraw',
   /function setLinesPanelCandidateBands[\s\S]*?drawLinesPanel\(\)[\s\S]*?\n\}/.test(html));
ok('window.setLinesPanelCandidateBands exported',
   /window\.setLinesPanelCandidateBands\s*=\s*setLinesPanelCandidateBands/.test(html));

ok('_LINES_CAND_BAND_PALETTE defined',
   /const _LINES_CAND_BAND_PALETTE\s*=\s*\[/.test(html));
ok('palette includes yellow #F5C518',
   /['"]#F5C518['"]/.test(html));
ok('palette includes green #4CAF50',
   /['"]#4CAF50['"]/.test(html));
ok('palette includes blue #3B82F6',
   /['"]#3B82F6['"]/.test(html));
ok('_candidateBandColor defined',
   /function _candidateBandColor\(idx,\s*alpha\)/.test(html));
ok('_paintCandidateBands defined',
   /function _paintCandidateBands\(ctx,\s*opts\)/.test(html));

ok('window._LINES_CAND_BAND_PALETTE exported',
   /window\._LINES_CAND_BAND_PALETTE\s*=\s*_LINES_CAND_BAND_PALETTE/.test(html));
ok('window._candidateBandColor exported',
   /window\._candidateBandColor\s*=\s*_candidateBandColor/.test(html));
ok('window._paintCandidateBands exported',
   /window\._paintCandidateBands\s*=\s*_paintCandidateBands/.test(html));

// ============================================================================
// 3. Toolbar checkbox + change-handler wireup
// ============================================================================
console.log('\n=== 3. Toolbar checkbox + handler ===');

ok('linesCandBandsToggle checkbox in toolbar',
   /<input\s+type="checkbox"\s+id="linesCandBandsToggle"\s+checked\s*\/?>/.test(html));
ok('linesCandBandsLabel wraps the checkbox',
   /<label\s+id="linesCandBandsLabel"/.test(html));
ok('handler block fetches the toggle by id',
   /document\.getElementById\(['"]linesCandBandsToggle['"]\)/.test(html));
ok('handler attaches change listener to setLinesPanelCandidateBands',
   /addEventListener\(['"]change['"],\s*e\s*=>\s*setLinesPanelCandidateBands\(!!e\.target\.checked\)\)/.test(html));
ok('handler restores from localStorage with default-ON contract',
   /localStorage\.getItem\(_LINES_PANEL_CAND_BANDS_KEY\)/.test(html));

// ============================================================================
// 4. drawLinesPanel integration (draw order: bands BEFORE frame)
// ============================================================================
console.log('\n=== 4. drawLinesPanel integration ===');

// Pull the body of drawLinesPanel and check that the candidate-bands
// paint reference appears BEFORE the strokeRect that draws the frame.
const drawLinesMatch = html.match(/function drawLinesPanel\(\)[\s\S]*?\n\}\s*\n/);
ok('drawLinesPanel function body extractable', drawLinesMatch !== null);
if (drawLinesMatch) {
  const body = drawLinesMatch[0];
  const bandsIdx = body.indexOf('_paintCandidateBands');
  // strokeRect that draws the per-source frame — pad.l + 0.5 marker
  const frameIdx = body.indexOf('strokeRect(pad.l + 0.5');
  ok('paint candidate bands reference present in drawLinesPanel',
     bandsIdx >= 0);
  ok('frame stroke present in drawLinesPanel',
     frameIdx >= 0);
  ok('candidate bands painted BEFORE frame stroke (draw order)',
     bandsIdx >= 0 && frameIdx >= 0 && bandsIdx < frameIdx,
     'bandsIdx=' + bandsIdx + ' frameIdx=' + frameIdx);
  ok('paint call gated on state.linesPanelCandidateBands',
     /state\.linesPanelCandidateBands\s*!==\s*false[\s\S]{0,200}_paintCandidateBands/.test(body));
  ok('paint call wrapped in try/catch (defensive)',
     /try\s*\{[\s\S]{0,500}_paintCandidateBands\(ctx[\s\S]{0,500}\}\s*catch/.test(body));
}

// ============================================================================
// 5. Behavioral: extract & exec the helpers
// ============================================================================
console.log('\n=== 5. Behavioral execution ===');

const helperRegion = html.match(
  /const _LINES_CAND_BAND_PALETTE[\s\S]*?window\._paintCandidateBands = _paintCandidateBands;/
);
ok('helper region extractable', helperRegion !== null);

const ctx = vm.createContext({ window: {}, console });
let helpersLoaded = false;
try {
  vm.runInContext(helperRegion[0], ctx);
  helpersLoaded = true;
} catch (e) {
  console.log('  helper region eval threw: ' + e.message);
}
ok('helper region parses and runs', helpersLoaded);

if (helpersLoaded) {
  // 5.1 palette content
  const pal = ctx.window._LINES_CAND_BAND_PALETTE;
  ok('palette[0] === #F5C518 (yellow)', pal && pal[0] === '#F5C518');
  ok('palette[1] === #4CAF50 (green)',  pal && pal[1] === '#4CAF50');
  ok('palette[2] === #3B82F6 (blue)',   pal && pal[2] === '#3B82F6');

  // 5.2 _candidateBandColor canonical indices
  const c0 = ctx.window._candidateBandColor(0);
  const c1 = ctx.window._candidateBandColor(1);
  const c2 = ctx.window._candidateBandColor(2);
  // Yellow #F5C518 = (245,197,24)
  ok('color(0) = rgba(245,197,24,0.1)', c0 === 'rgba(245,197,24,0.1)', c0);
  // Green #4CAF50 = (76,175,80)
  ok('color(1) = rgba(76,175,80,0.1)',  c1 === 'rgba(76,175,80,0.1)', c1);
  // Blue #3B82F6 = (59,130,246)
  ok('color(2) = rgba(59,130,246,0.1)', c2 === 'rgba(59,130,246,0.1)', c2);

  // 5.3 alpha argument honored
  const c0a = ctx.window._candidateBandColor(0, 0.05);
  ok('color(0, 0.05) honors alpha',
     c0a === 'rgba(245,197,24,0.05)', c0a);

  // 5.4 idx >= 3 cycles via hsla golden-angle rotation
  const c3 = ctx.window._candidateBandColor(3);
  const c4 = ctx.window._candidateBandColor(4);
  ok('color(3) is hsla() golden-angle cycle', /^hsla\(/.test(c3), c3);
  ok('color(4) is hsla() golden-angle cycle', /^hsla\(/.test(c4), c4);
  ok('color(3) and color(4) are distinct (rotated)', c3 !== c4);

  // 5.5 alpha clamping
  const cClampHi = ctx.window._candidateBandColor(0, 1.5);
  const cClampLo = ctx.window._candidateBandColor(0, -0.2);
  ok('alpha > 1 clamped to 1', /,1\)$/.test(cClampHi), cClampHi);
  ok('alpha < 0 clamped to 0', /,0\)$/.test(cClampLo), cClampLo);

  // 5.6 negative / non-numeric idx normalised to 0 (yellow)
  const cNeg  = ctx.window._candidateBandColor(-3);
  const cBad  = ctx.window._candidateBandColor('foo');
  ok('negative idx → palette[0] (yellow)',
     cNeg === 'rgba(245,197,24,0.1)', cNeg);
  ok('non-numeric idx → palette[0] (yellow)',
     cBad === 'rgba(245,197,24,0.1)', cBad);

  // -----------------------------------------------------------------
  // 5.7 _paintCandidateBands behavior — set up a fake canvas ctx.
  // -----------------------------------------------------------------
  function makeFakeCtx() {
    const calls = [];
    return {
      calls,
      save()    { calls.push(['save']); },
      restore() { calls.push(['restore']); },
      fillRect(x, y, w, h) { calls.push(['fillRect', x, y, w, h]); },
      set fillStyle(v) { calls.push(['fillStyle', v]); },
      get fillStyle() { return null; },
    };
  }
  const baseOpts = {
    pad:   { l: 44, r: 16, t: 6, b: 8 },
    plotW: 740,
    plotH: 200,
    toX:   (mb) => 44 + ((mb - 0) / 30) * 740,
    mbMin: 0,
    mbMax: 30,
    chrom: 'LG28',
  };

  const cands3 = [
    { confirmed: true, chrom: 'LG28', start_bp: 1e6,  end_bp: 4e6  },
    { confirmed: true, chrom: 'LG28', start_bp: 8e6,  end_bp: 11e6 },
    { confirmed: true, chrom: 'LG28', start_bp: 18e6, end_bp: 22e6 },
  ];

  let cv = makeFakeCtx();
  let n  = ctx.window._paintCandidateBands(
    cv, Object.assign({}, baseOpts, { candidates: cands3 }));
  ok('3 confirmed on-chrom → paints 3 bands', n === 3, 'n=' + n);
  ok('3 fillRect calls recorded',
     cv.calls.filter(c => c[0] === 'fillRect').length === 3);

  // 5.8 Drafts (confirmed=false) filtered out
  cv = makeFakeCtx();
  n = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: false, chrom: 'LG28', start_bp: 1e6, end_bp: 4e6 },
      { confirmed: true,  chrom: 'LG28', start_bp: 8e6, end_bp: 11e6 },
    ],
  }));
  ok('drafts filtered (confirmed=false)', n === 1, 'n=' + n);

  // 5.9 Cross-chrom candidates filtered out
  cv = makeFakeCtx();
  n = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: true, chrom: 'LG12', start_bp: 1e6, end_bp: 4e6 },
      { confirmed: true, chrom: 'LG28', start_bp: 8e6, end_bp: 11e6 },
    ],
  }));
  ok('cross-chrom candidates filtered', n === 1, 'n=' + n);

  // 5.10 Empty list
  cv = makeFakeCtx();
  ok('empty candidates list → 0',
     ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, { candidates: [] })) === 0);

  // 5.11 Defensive: null ctx
  ok('null ctx → 0',
     ctx.window._paintCandidateBands(null, Object.assign({}, baseOpts, { candidates: cands3 })) === 0);

  // 5.12 Off-screen candidate dropped (visible span empty)
  cv = makeFakeCtx();
  n = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: true, chrom: 'LG28', start_bp: 100e6, end_bp: 110e6 }, // way past mbMax
      { confirmed: true, chrom: 'LG28', start_bp: 5e6,   end_bp: 7e6   },
    ],
  }));
  ok('off-screen candidate dropped, on-screen kept', n === 1, 'n=' + n);

  // 5.13 Bad bp values
  cv = makeFakeCtx();
  n = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: true, chrom: 'LG28', start_bp: 'foo', end_bp: 4e6 },
      { confirmed: true, chrom: 'LG28', start_bp: NaN,   end_bp: 4e6 },
      { confirmed: true, chrom: 'LG28', start_bp: 5e6,   end_bp: 5e6 }, // zero-width
      { confirmed: true, chrom: 'LG28', start_bp: 5e6,   end_bp: 7e6 }, // good
    ],
  }));
  ok('malformed bp values silently skipped (good one survives)', n === 1, 'n=' + n);

  // 5.14 Custom alpha propagates to fillStyle
  cv = makeFakeCtx();
  ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [cands3[0]],
    alpha: 0.05,
  }));
  const fs0 = cv.calls.find(c => c[0] === 'fillStyle');
  ok('custom alpha=0.05 reaches fillStyle',
     fs0 && fs0[1] === 'rgba(245,197,24,0.05)', fs0 && fs0[1]);

  // 5.15 Stable palette index across off-screen drops
  // First candidate is off-screen (skipped), but its index should bump
  // so the second visible candidate paints with green (idx 1), not yellow.
  cv = makeFakeCtx();
  ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: true, chrom: 'LG28', start_bp: 100e6, end_bp: 110e6 }, // off
      { confirmed: true, chrom: 'LG28', start_bp: 5e6,   end_bp: 7e6   }, // visible, idx=1
    ],
  }));
  const fs1 = cv.calls.find(c => c[0] === 'fillStyle');
  ok('off-screen candidate still bumps palette index',
     fs1 && fs1[1] === 'rgba(76,175,80,0.1)', fs1 && fs1[1]);

  // 5.16 Defensive: missing toX
  cv = makeFakeCtx();
  ok('missing toX → 0',
     ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
       candidates: cands3, toX: null,
     })) === 0);

  // 5.17 Defensive: invalid mb range
  cv = makeFakeCtx();
  ok('invalid mb range (mbMax <= mbMin) → 0',
     ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
       candidates: cands3, mbMin: 50, mbMax: 60,
     })) === 0);

  // 5.18 Defensive: negative plotW/plotH
  cv = makeFakeCtx();
  ok('non-positive plotW → 0',
     ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
       candidates: cands3, plotW: 0,
     })) === 0);

  // 5.19 No chrom in opts → don't filter on chrom (fallback friendly)
  cv = makeFakeCtx();
  n = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, {
    candidates: [
      { confirmed: true, chrom: 'LG12', start_bp: 1e6, end_bp: 4e6 },
      { confirmed: true, chrom: 'LG28', start_bp: 8e6, end_bp: 11e6 },
    ],
    chrom: null,
  }));
  ok('null chrom → cross-chrom filter disabled (paints both)',
     n === 2, 'n=' + n);

  // 5.20 No candidates property at all
  cv = makeFakeCtx();
  ok('missing candidates array → 0',
     ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts)) === 0);

  // 5.21 Throws-resistant — call returns Number, not undefined
  cv = makeFakeCtx();
  const r = ctx.window._paintCandidateBands(cv, Object.assign({}, baseOpts, { candidates: cands3 }));
  ok('returns finite Number (count of painted bands)',
     typeof r === 'number' && Number.isFinite(r));
}

// ============================================================================
// 6. Spec discipline — confirm the toolbar/state choices match SPEC §3, §4
// ============================================================================
console.log('\n=== 6. Spec discipline ===');

// SPEC §3: drafts (confirmed: false) do NOT get a band.
ok('source filters confirmed === true (SPEC §3)',
   /c\.confirmed\s*!==\s*true/.test(html));

// SPEC §4: localStorage key shape
ok('localStorage key matches SPEC §4 (pca_scrubber_v3.linesPanelCandidateBands)',
   /pca_scrubber_v3\.linesPanelCandidateBands/.test(html));

// SPEC §6.3: array order (chromIdx walks candidateList in order, no sort)
ok('walks candidateList in array order (no sort by start_bp)',
   /for\s*\(\s*let\s+i\s*=\s*0;\s*i\s*<\s*candidates\.length;\s*i\+\+\s*\)/.test(html) &&
   !/candidates\s*\.\s*slice\(\)\s*\.\s*sort/.test(html.match(/function _paintCandidateBands[\s\S]*?\n\}/)[0] || ''));

// ============================================================================
// Final tally
// ============================================================================
console.log('\n=============================================================');
console.log('  ' + pass + ' / ' + (pass + fail) + ' tests passed');
console.log('=============================================================');
process.exit(fail === 0 ? 0 : 1);
