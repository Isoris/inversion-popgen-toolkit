// =============================================================================
// turn 140 — H-label classification chip on PCA per-band legend
// =============================================================================
// Quentin (turn 138):
//   "a mini round chip icon to classify the per sample lines [...] background
//    color tells you the system, age slot reserved for future age data."
//
// Tests cover:
//   - _hlabelChipColor(classification) → correct hex for HOM/HET/AMBIGUOUS/
//     NO_DOSAGE; falls back to grey for unknown values
//   - _HLABEL_CHIP_COLOR exported
//   - _hlabelDrawChip writes to the canvas context (fill + stroke + arc calls)
//   - drawPCA legend block: chipSlot widens legendW when focal pane has
//     h_classification; collapses to 0 otherwise
//   - Chip drawing branch: gated on (paneOffset === 0) AND
//     (state.candidate.h_classification.bands)
//   - Age overlay slot: looked up via candidate.age_my[bandIdx], drawn when
//     finite, skipped when null/missing
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
// 1. Source-level: helpers exist + window-exported
// ============================================================================
console.log('\n=== 1. Source-level definitions ===');

ok('_hlabelChipColor defined',     /function _hlabelChipColor\b/.test(html));
ok('_hlabelDrawChip defined',      /function _hlabelDrawChip\b/.test(html));
ok('_HLABEL_CHIP_COLOR constant',  /const _HLABEL_CHIP_COLOR\s*=\s*\{/.test(html));
ok('window._hlabelChipColor exported',
   /window\._hlabelChipColor\s*=\s*_hlabelChipColor/.test(html));
ok('window._hlabelDrawChip exported',
   /window\._hlabelDrawChip\s*=\s*_hlabelDrawChip/.test(html));
ok('window._HLABEL_CHIP_COLOR exported',
   /window\._HLABEL_CHIP_COLOR\s*=\s*_HLABEL_CHIP_COLOR/.test(html));

// ============================================================================
// 2. _HLABEL_CHIP_COLOR — palette content
// ============================================================================
console.log('\n=== 2. Chip color palette ===');

ok('palette defines HOM',        /HOM:\s*'#3a7dde'/.test(html));
ok('palette defines HET',        /HET:\s*'#d97842'/.test(html));
ok('palette defines AMBIGUOUS',  /AMBIGUOUS:\s*'#f5a524'/.test(html));
ok('palette defines NO_DOSAGE',  /NO_DOSAGE:\s*'#888888'/.test(html));

// ============================================================================
// 3. Behavioral: extract & exec _hlabelChipColor
// ============================================================================
console.log('\n=== 3. _hlabelChipColor behavior ===');

function extractFn(name) {
  const re = new RegExp('function ' + name + '\\([\\s\\S]*?\\n\\}', 'm');
  const m = html.match(re);
  return m ? m[0] : null;
}

const colorSrc = extractFn('_hlabelChipColor');
ok('source extractable', colorSrc !== null);

const ctx = vm.createContext({});
vm.runInContext(`
  const _HLABEL_CHIP_COLOR = {
    HOM:        '#3a7dde',
    HET:        '#d97842',
    AMBIGUOUS:  '#f5a524',
    NO_DOSAGE:  '#888888',
  };
  ${colorSrc}
`, ctx);

ok('HOM → blue',          ctx._hlabelChipColor('HOM') === '#3a7dde');
ok('HET → orange',        ctx._hlabelChipColor('HET') === '#d97842');
ok('AMBIGUOUS → yellow',  ctx._hlabelChipColor('AMBIGUOUS') === '#f5a524');
ok('NO_DOSAGE → grey',    ctx._hlabelChipColor('NO_DOSAGE') === '#888888');
ok('null → grey fallback', ctx._hlabelChipColor(null) === '#888888');
ok("'' → grey fallback (empty string ≈ falsy)",
   ctx._hlabelChipColor('') === '#888888');
ok('unknown value → grey fallback',
   ctx._hlabelChipColor('SOMETHING_ELSE') === '#888888');

// ============================================================================
// 4. _hlabelDrawChip — canvas calls
// ============================================================================
console.log('\n=== 4. _hlabelDrawChip canvas calls ===');

const drawSrc = extractFn('_hlabelDrawChip');
ok('source extractable', drawSrc !== null);

function makeMockCtx() {
  const calls = [];
  return {
    _calls: calls,
    fillStyle: null,
    strokeStyle: null,
    lineWidth: null,
    font: null,
    textAlign: null,
    textBaseline: null,
    beginPath() { calls.push(['beginPath']); },
    arc(x, y, r, a, b) { calls.push(['arc', x, y, r, a, b]); },
    fill() { calls.push(['fill', this.fillStyle]); },
    stroke() { calls.push(['stroke', this.strokeStyle, this.lineWidth]); },
    fillText(t, x, y) { calls.push(['fillText', t, x, y, this.fillStyle]); },
  };
}

const drawCtx = vm.createContext({
  Math, Number,
  _hlabelChipColor: ctx._hlabelChipColor,
});
vm.runInContext(`
  const _HLABEL_CHIP_COLOR = {
    HOM:        '#3a7dde',
    HET:        '#d97842',
    AMBIGUOUS:  '#f5a524',
    NO_DOSAGE:  '#888888',
  };
  ${drawSrc}
`, drawCtx);

// 4a. Chip without age — no fillText call
{
  const c = makeMockCtx();
  drawCtx._hlabelDrawChip(c, 100, 50, 4, 'HOM', null, '#aaa');
  const arcCalls = c._calls.filter(x => x[0] === 'arc');
  const fillCalls = c._calls.filter(x => x[0] === 'fill');
  const strokeCalls = c._calls.filter(x => x[0] === 'stroke');
  const textCalls = c._calls.filter(x => x[0] === 'fillText');
  ok('arc called twice (fill + stroke paths)', arcCalls.length === 2);
  ok('fill called once with HOM blue',
     fillCalls.length === 1 && fillCalls[0][1] === '#3a7dde');
  ok('stroke called once', strokeCalls.length === 1);
  ok('outline color = passed inkDim', strokeCalls[0][1] === '#aaa');
  ok('no fillText (no age provided)', textCalls.length === 0);
}

// 4b. Chip with age — fillText present
{
  const c = makeMockCtx();
  drawCtx._hlabelDrawChip(c, 100, 50, 4, 'HET', 2.3, '#aaa');
  const textCalls = c._calls.filter(x => x[0] === 'fillText');
  ok('fillText called when ageMy is finite', textCalls.length === 1);
  ok('text content is "2.3"', textCalls[0][1] === '2.3');
  ok('text uses white fill for legibility', textCalls[0][4] === '#fff');
  // The sequence should be: fill bg → arc again → stroke → font/textBaseline/etc → fillText
  // Just verify centered around (100, 50)
  ok('text x ≈ chip cx (100)', textCalls[0][2] === 100);
}

// 4c. Chip with NaN/undefined age — no fillText
{
  const c = makeMockCtx();
  drawCtx._hlabelDrawChip(c, 100, 50, 4, 'HOM', NaN, '#aaa');
  ok('NaN age skipped', c._calls.filter(x => x[0] === 'fillText').length === 0);

  const c2 = makeMockCtx();
  drawCtx._hlabelDrawChip(c2, 100, 50, 4, 'HOM', undefined, '#aaa');
  ok('undefined age skipped', c2._calls.filter(x => x[0] === 'fillText').length === 0);
}

// 4d. Defensive: ctx without arc — no-op rather than throw
{
  const stub = { fillStyle: null };   // no beginPath, no arc
  let threw = false;
  try { drawCtx._hlabelDrawChip(stub, 0, 0, 4, 'HOM', null, '#aaa'); }
  catch (_) { threw = true; }
  ok('no throw on minimal ctx', !threw);
}

// 4e. Unknown classification → grey chip (defensive fallback)
{
  const c = makeMockCtx();
  drawCtx._hlabelDrawChip(c, 0, 0, 4, 'WHO_KNOWS', null, '#aaa');
  const fills = c._calls.filter(x => x[0] === 'fill');
  ok('unknown classification falls back to grey',
     fills.length === 1 && fills[0][1] === '#888888');
}

// ============================================================================
// 5. Legend integration — chipSlot adapts to focal pane + h_classification
// ============================================================================
console.log('\n=== 5. Legend chip integration (source-level) ===');

ok('legend reads state.candidate for h_classification',
   /_stateCand[\s\S]{0,200}h_classification/.test(html));
ok('chip drawn only on paneOffset === 0',
   /paneOffset === 0[\s\S]{0,400}h_classification/.test(html));
ok('chipBands gated on Array.isArray(h_classification.bands)',
   /Array\.isArray\(_stateCand\.h_classification\.bands\)/.test(html));
ok('chipSlot is 0 when chipBands is null (collapses cleanly)',
   /chipSlot\s*=\s*chipBands\s*\?\s*\(?chipR\s*\*\s*2\s*\+\s*2\s*\)?\s*:\s*0/.test(html));
ok('legendW adds chipSlot',
   /legendW\s*=\s*swW\s*\+\s*4\s*\+\s*chipSlot\s*\+\s*widest/.test(html));
ok('label x-position offset by chipSlot',
   /fillText\(legendLabelFor\(k\),\s*lx\s*\+\s*swW\s*\+\s*4\s*\+\s*chipSlot/.test(html));
ok('chip cx anchored to swatch right edge',
   /cx\s*=\s*lx\s*\+\s*swW\s*\+\s*2\s*\+\s*chipR/.test(html));
ok('chip cy centered in swatch row',
   /cy\s*=\s*ry\s*\+\s*1\s*\+\s*swH\s*\/\s*2/.test(html));

// ============================================================================
// 6. Age overlay slot reserved (forward-compat with future age JSON layer)
// ============================================================================
console.log('\n=== 6. Age overlay slot ===');

ok('candidate.age_my looked up by band index',
   /_stateCand\.age_my\s*&&\s*_stateCand\.age_my\[k\]\s*!=\s*null/.test(html));
ok('age slot defaults to null when missing',
   /\?\s*_stateCand\.age_my\[k\]\s*:\s*null/.test(html));
// _hlabelDrawChip receives ageMy as 6th arg
ok('_hlabelDrawChip called with ageMy parameter',
   /_hlabelDrawChip\(ctx,\s*cx,\s*cy,\s*chipR,\s*bandClass,\s*ageMy,\s*themeColor/.test(html));

// ============================================================================
// 7. Defensive guards: try/catch around _hlabelDrawChip + typeof checks
// ============================================================================
console.log('\n=== 7. Defensive guards ===');

ok("typeof _hlabelDrawChip === 'function' guard",
   /typeof _hlabelDrawChip === 'function'/.test(html));
ok('try/catch around chip draw call',
   /try\s*\{[\s\S]{0,300}_hlabelDrawChip\(ctx[\s\S]{0,200}\}\s*catch\s*\(_\)/.test(html));

// ============================================================================
// 8. Summary
// ============================================================================
console.log('\n=== Summary ===');
console.log('PASS: ' + pass);
console.log('FAIL: ' + fail);
process.exit(fail === 0 ? 0 : 1);
