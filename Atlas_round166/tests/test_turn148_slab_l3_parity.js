// =============================================================================
// turn 148 — Slab-mode L3 panel parity with L2 mode (queued screenshot fix #3)
// =============================================================================
// Quentin (turn 144 screenshot 3): "1w/Nw L3 panel parity with L2 mode —
// renderL3PanelSlab needs the same rich pane structure as renderL3Panel.
// Currently the slab-mode (1w / 5w / 10w / Nw) renderer is a downgraded
// variant of the L2-mode renderer."
//
// What this turn ports from L2 mode into slab mode:
//   - Per-pane toolbar mirrors (K group, color group, recluster) via
//     _l3PaneHeaderToolsHtml() — clicks already wired by _wireL3PaneToolsDelegation
//   - Karyotype chips ABOVE the mini-PCA on focal panes via
//     _kSpecificMetaInlineHtml(cl, null) — l2idx=null skips the L2-only
//     sub-band-context chip cleanly
//   - Spotlight click setup on every mini-PCA via _setupL3MiniClick
//   - canvas.__l3_render populated by drawSlabMiniPCA (with isSlab:true and
//     slabRange) so the click handler can hit-test against slab-mean PCs
//   - Slab-aware hit-test branch in _setupL3MiniClick (uses aggregateSlab to
//     match what drawSlabMiniPCA actually plots, otherwise click would
//     mis-target at W>1)
//   - Band-selector strip on focal pane in candidate mode via
//     _l3BandSelectorHtml(K, bandCounts) — already K-and-counts-only so no
//     L2 dependency
//   - col.__l3_meta with isSlab:true + leftRange/rightRange + invPerm so
//     _applySpotlightHighlights can resolve cells
//   - Slab-aware cluster lookup in _applySpotlightHighlights (routes via
//     getSlabClusterAt(range[0], range[1], K) when meta.isSlab is set)
//   - _applySpotlightHighlights post-pass call from renderL3PanelSlab
//
// Scope discipline: this is a structural-parity port. The full focal-content
// body of L2 mode (heterozygosity ridgelines, family-purity diagnostics,
// invariant-meta hoisted block in K=both, etc.) remains slab's lighter
// slabFocalContentHtml — that's a much larger separate undertaking.
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
// 1. renderL3PanelSlab — per-pane toolbar parity
// ============================================================================
console.log('\n=== 1. Per-pane toolbar in slab mode ===');

// Extract the renderL3PanelSlab function body so we can pattern-check inside
// it without false matches from the L2 renderer just above.
const slabFnRe = /function\s+renderL3PanelSlab\s*\(\s*\)\s*\{([\s\S]*?)\n\}\s*\n\s*\/\/[^\n]*Compact focal-content HTML for a slab/;
const slabFnMatch = html.match(slabFnRe);
ok('renderL3PanelSlab body extracted from source',
   !!slabFnMatch, 'function boundary regex did not match');
const slabBody = slabFnMatch ? slabFnMatch[1] : '';

ok('slab renderer calls _l3PaneHeaderToolsHtml (per-pane toolbar)',
   /paneToolsHtml\s*=\s*\(typeof\s+_l3PaneHeaderToolsHtml\s*===\s*'function'\)\s*\?\s*_l3PaneHeaderToolsHtml\(\)/.test(slabBody));

ok('slab renderer interpolates paneToolsHtml into focal head innerHTML',
   /◆ slab focal[\s\S]{0,500}?(?:\$\{paneToolsHtml\}|\+\s*\n?\s*paneToolsHtml)/.test(slabBody));

ok('slab renderer interpolates paneToolsHtml into neighbor head innerHTML',
   /slab \$\{offset > 0[\s\S]{0,500}?(?:\$\{paneToolsHtml\}|\+\s*\n?\s*paneToolsHtml)/.test(slabBody));

ok('slab focal title wrapped in .l3-pane-title (matches L2 mode convention)',
   /class="l3-pane-title"[^>]*>◆ slab focal/.test(slabBody));

// The old "SIM —" placeholder was removed (was a meaningless L2-aping element
// that didn't reflect actual slab similarity). Confirm it's gone from the
// slab body.
ok('old "SIM —" placeholder span removed from slab focal head',
   !/SIM\s*—/.test(slabBody),
   'placeholder still present in slab head');

// ============================================================================
// 2. Karyotype chips ABOVE the mini-PCA (focal pane)
// ============================================================================
console.log('\n=== 2. Karyotype chips above mini-PCA ===');

ok('slab renderer calls _kSpecificMetaInlineHtml for focal',
   /_kSpecificMetaInlineHtml\(clFocal,\s*null\)/.test(slabBody));

ok('slab renderer appends focalChipsAbove BEFORE the mini-PCA wrapper',
   /focalChipsAbove[\s\S]{0,400}?col\.appendChild\(focalChipsAbove\)[\s\S]{0,800}?const miniWrap = document\.createElement\('div'\)/.test(slabBody));

ok('focalChipsAbove uses the small 9.5px font (parity with L2 mode)',
   /focalChipsAbove\.style\.cssText\s*=\s*'font-size: 9\.5px/.test(slabBody));

// Passing l2idx=null is critical — it suppresses the sub-band annotation
// chip which only makes sense for L2 envelopes in active draft context.
ok('null passed for l2idx (suppresses sub-band chip cleanly)',
   /_kSpecificMetaInlineHtml\(clFocal,\s*null\)/.test(slabBody));

// ============================================================================
// 3. Spotlight click setup wired
// ============================================================================
console.log('\n=== 3. Spotlight click on slab mini-PCAs ===');

ok('_setupL3MiniClick called on slab mini-PCA canvas',
   /if \(typeof _setupL3MiniClick === 'function'\) _setupL3MiniClick\(mini\)/.test(slabBody));

// ============================================================================
// 4. drawSlabMiniPCA populates __l3_render so click handler has hit-test data
// ============================================================================
console.log('\n=== 4. drawSlabMiniPCA __l3_render context ===');

const drawSlabRe = /function\s+drawSlabMiniPCA\s*\([^)]*\)\s*\{([\s\S]*?)\n\}\s*\n/;
const drawSlabMatch = html.match(drawSlabRe);
ok('drawSlabMiniPCA body extracted',
   !!drawSlabMatch, 'function boundary regex did not match');
const drawSlabBody = drawSlabMatch ? drawSlabMatch[1] : '';

ok('drawSlabMiniPCA assigns canvas.__l3_render',
   /canvas\.__l3_render\s*=\s*\{/.test(drawSlabBody));

ok('__l3_render carries isSlab:true marker',
   /canvas\.__l3_render\s*=\s*\{[\s\S]{0,400}?isSlab:\s*true/.test(drawSlabBody));

ok('__l3_render carries slabRange',
   /canvas\.__l3_render\s*=\s*\{[\s\S]{0,400}?slabRange:/.test(drawSlabBody));

ok('__l3_render carries pad / plotW / plotH / x-y bounds (parity with L2 mode)',
   /canvas\.__l3_render\s*=\s*\{[\s\S]{0,400}?pad,\s*plotW,\s*plotH,[\s\S]{0,200}?xMin,\s*xMax,\s*yMin,\s*yMax/.test(drawSlabBody));

ok('wMid is set to the slab center via bit-shift midpoint',
   /const slabMid\s*=\s*\(range\[0\]\s*\+\s*range\[1\]\)\s*>>\s*1/.test(drawSlabBody));

// ============================================================================
// 5. _setupL3MiniClick — slab-aware hit-test branch
// ============================================================================
console.log('\n=== 5. Slab-aware hit-test branch ===');

const clickFnRe = /function\s+_setupL3MiniClick\s*\(\s*canvas\s*\)\s*\{([\s\S]*?)\n\}\s*\n/;
const clickMatch = html.match(clickFnRe);
ok('_setupL3MiniClick body extracted',
   !!clickMatch);
const clickBody = clickMatch ? clickMatch[1] : '';

ok('click handler branches on ctx.isSlab',
   /if \(ctx\.isSlab\s*&&\s*ctx\.slabRange\)/.test(clickBody));

ok('slab branch calls aggregateSlab for slab-mean PCs',
   /aggregateSlab\(ctx\.slabRange\[0\],\s*ctx\.slabRange\[1\]\)/.test(clickBody));

ok('slab branch sets signMul = 1 (sign already applied in agg.xs)',
   /signMul\s*=\s*1/.test(clickBody));

ok('L2 branch retained: const pcRes = getPC(ctx.wMid)',
   /const pcRes\s*=\s*getPC\(ctx\.wMid\)/.test(clickBody));

ok('hit-test loop now uses pc1Plot / pc2Plot (not pc1 / pc2 directly)',
   /pc1Plot\[si\]\s*\*\s*signMul[\s\S]{0,80}?pc2Plot\[si\]/.test(clickBody));

// ============================================================================
// 6. Band-selector strip on focal pane in candidate mode
// ============================================================================
console.log('\n=== 6. Band-selector strip in candidate mode ===');

ok('slab renderer calls _l3BandSelectorHtml (focal in candidate mode)',
   /if \(isFocal && typeof _l3BandSelectorHtml === 'function'\)/.test(slabBody));

ok('slab band selector reads cluster labels for bandCounts',
   /bandCounts\s*=\s*new Array\(K\)\.fill\(0\)[\s\S]{0,300}?for \(let s = 0; s < labels\.length/.test(slabBody));

ok('slab band selector calls _wireBandSelectorClicks after insertion',
   /if \(typeof _wireBandSelectorClicks === 'function'\)\s*\{\s*_wireBandSelectorClicks\(\);\s*\}/.test(slabBody));

// ============================================================================
// 7. col.__l3_meta stash for slab columns
// ============================================================================
console.log('\n=== 7. col.__l3_meta on slab neighbor columns ===');

ok('slab renderer builds invPerm from cmp.perm',
   // turn 149: invPerm is now built against cmp's K (ctK) instead of the
   // section K, so kmeans-K6 reMode against l3KMode='k3' sections gets
   // the correct K=6 contingency dimensions.
   /const _invPerm\s*=\s*new Array\((?:K|ctK)\);[\s\S]{0,400}?_invPerm\[cmp\.perm\[r\]\]\s*=\s*r/.test(slabBody));

ok('slab __l3_meta carries isSlab:true',
   /col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?isSlab:\s*true/.test(slabBody));

ok('slab __l3_meta carries leftRange / rightRange (not leftIdx / rightIdx)',
   /col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?leftRange:[\s\S]{0,200}?rightRange:/.test(slabBody));

ok('slab __l3_meta carries invPerm',
   /col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?invPerm:\s*_invPerm/.test(slabBody));

// ============================================================================
// 8. _applySpotlightHighlights post-pass call from slab renderer
// ============================================================================
console.log('\n=== 8. _applySpotlightHighlights call from slab renderer ===');

ok('slab renderer calls _applySpotlightHighlights at end of layout loop',
   /_applySpotlightHighlights\(null\)/.test(slabBody));

// The call is wrapped in try/catch (parity with L2 mode line 44823)
ok('call wrapped in try/catch (fail-soft)',
   /try\s*\{\s*_applySpotlightHighlights\(null\);\s*\}\s*catch/.test(slabBody));

// ============================================================================
// 9. _applySpotlightHighlights — slab-aware cluster lookup
// ============================================================================
console.log('\n=== 9. _applySpotlightHighlights slab branch ===');

const spotFnRe = /function\s+_applySpotlightHighlights\s*\([^)]*\)\s*\{([\s\S]*?)\n\}\s*\n/;
const spotMatch = html.match(spotFnRe);
ok('_applySpotlightHighlights body extracted',
   !!spotMatch);
const spotBody = spotMatch ? spotMatch[1] : '';

ok('spotlight body branches on meta.isSlab',
   /if \(meta\.isSlab\)/.test(spotBody));

ok('slab branch uses getSlabClusterAt with meta.leftRange',
   /getSlabClusterAt\(meta\.leftRange\[0\],\s*meta\.leftRange\[1\],\s*meta\.K\)/.test(spotBody));

ok('slab branch uses getSlabClusterAt with meta.rightRange',
   /getSlabClusterAt\(meta\.rightRange\[0\],\s*meta\.rightRange\[1\],\s*meta\.K\)/.test(spotBody));

ok('L2 path retained: getL2Cluster(meta.leftIdx)',
   /getL2Cluster\(meta\.leftIdx\)/.test(spotBody));

ok('slab branch defends against missing getSlabClusterAt',
   /typeof getSlabClusterAt\s*!==\s*'function'/.test(spotBody));

// ============================================================================
// 10. Sandboxed run of _applySpotlightHighlights against a mock slab column
// ============================================================================
console.log('\n=== 10. Sandboxed slab-spotlight execution ===');

// Build a minimal sandbox that exercises the slab branch end-to-end. We
// extract _applySpotlightHighlights and stub the helpers it calls.
{
  const _calls = { getSlab: [], getL2: [] };

  // Mock cluster: 4 samples, K=3. Sample 0 is in band 0 in left slab and
  // band 1 in right slab. The aligned permutation maps right-band 1 → 0,
  // so sample 0 should land in cell (r=0, c=0) (the diagonal).
  const leftCluster  = { labels: [0, 1, 2, 0], fixedKLabels: null };
  const rightCluster = { labels: [1, 2, 0, 1], fixedKLabels: null };

  // Track which cells got tagged
  const tagged = [];
  const taggedAxes = [];

  // Mock TD / TH for the queried cells. classList carries no spotlight
  // class initially (so the upgrade path isn't blocked).
  function mkCell(r, c, kind) {
    const classes = new Set();
    return {
      kind, r, c,
      classList: {
        contains: (n) => classes.has(n),
        add: (n) => { classes.add(n); if (kind === 'td') tagged.push({r,c,n}); else taggedAxes.push({r,c,n,kind}); },
        remove: (n) => classes.delete(n),
      },
      _classes: classes,
    };
  }

  // Mock col with __l3_meta for a slab pane: K=3, offset=+1.
  // leftRange=[0,4] (focal), rightRange=[5,9] (slab+1).
  // invPerm: right-label 1 → row 0 (so sample 0 with leftL=0, rightL=1
  // → r=0, c=invPerm[1]=0)
  const cellTD = mkCell(0, 0, 'td');
  const cellTH_r = mkCell(0, null, 'th');
  const cellTH_c = mkCell(null, 0, 'th');

  const slabCol = {
    __l3_meta: {
      K: 3, offset: 1,
      isSlab: true,
      leftRange: [0, 4],
      rightRange: [5, 9],
      invPerm: [2, 0, 1],     // right band 1 → aligned 0 (diag with leftL=0)
    },
    querySelector: (sel) => {
      // strip the leading td/th and parse data attrs
      if (sel === 'td[data-r="0"][data-c="0"]') return cellTD;
      if (sel === 'th[data-r="0"]') return cellTH_r;
      if (sel === 'th[data-c="0"]') return cellTH_c;
      return null;
    },
    querySelectorAll: () => [],
  };

  const body = {
    querySelectorAll: (sel) => {
      if (sel.includes('.l3-col')) return [slabCol];
      if (sel.includes('spotlight-cell')) return [];     // no stale to clear
      return [];
    },
  };

  const sandbox = {
    document: {
      getElementById: (id) => (id === 'l3Body') ? body : null,
    },
    state: {
      data: { n_samples: 4 },
      spotlight: 0,                          // primary spotlight on sample 0
      spotlightTrackedAll: false,
      tracked: new Set(),
    },
    // Helper stubs
    getL2Cluster: (idx) => { _calls.getL2.push(idx); return null; },     // shouldn't be called for slab cols
    getSlabClusterAt: (s, e, K) => {
      _calls.getSlab.push([s, e, K]);
      if (s === 0 && e === 4) return leftCluster;
      if (s === 5 && e === 9) return rightCluster;
      return null;
    },
    console,
  };
  const ctx = vm.createContext(sandbox);

  // Extract the function source and run it in the sandbox so we exercise
  // the actual ported branch logic (not a re-implementation).
  const fnSourceMatch = html.match(/(function\s+_applySpotlightHighlights\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  ok('extracted _applySpotlightHighlights for sandbox run',
     !!fnSourceMatch);
  if (fnSourceMatch) {
    vm.runInContext(fnSourceMatch[1], ctx);
    // Now invoke. Pass null for curFocalL2 (slab mode doesn't have one).
    try {
      vm.runInContext('_applySpotlightHighlights(null)', ctx);

      ok('getSlabClusterAt called with leftRange',
         _calls.getSlab.some(c => c[0] === 0 && c[1] === 4 && c[2] === 3));
      ok('getSlabClusterAt called with rightRange',
         _calls.getSlab.some(c => c[0] === 5 && c[1] === 9 && c[2] === 3));
      ok('getL2Cluster NOT called for slab column',
         _calls.getL2.length === 0,
         'getL2Cluster was called: ' + JSON.stringify(_calls.getL2));

      // Sample 0: leftL=0, rightL=1, invPerm[1]=0 → (r=0, c=0). Diagonal cell
      // should have spotlight-cell-primary.
      ok('sample 0 cell (r=0, c=0) tagged spotlight-cell-primary',
         tagged.some(t => t.r === 0 && t.c === 0 && t.n === 'spotlight-cell-primary'));

      // Both axis headers should be tagged spotlight-axis (primary tagging).
      ok('row axis header tagged spotlight-axis',
         taggedAxes.some(t => t.kind === 'th' && t.n === 'spotlight-axis'));
      ok('col axis header tagged spotlight-axis',
         taggedAxes.length >= 2);
    } catch (e) {
      ok('sandbox invocation completed', false, e.message);
    }
  }
}

// ============================================================================
// 11. Sandboxed slab-spotlight: tracked-only mode (no primary)
// ============================================================================
console.log('\n=== 11. Sandboxed tracked-only spotlight ===');

{
  const tagged = [];
  function mkCell() {
    const classes = new Set();
    return {
      classList: {
        contains: (n) => classes.has(n),
        add: (n) => { classes.add(n); tagged.push(n); },
        remove: (n) => classes.delete(n),
      },
    };
  }
  const cell = mkCell();
  const slabCol = {
    __l3_meta: {
      K: 3, offset: 1, isSlab: true,
      leftRange: [0, 4], rightRange: [5, 9],
      invPerm: [2, 0, 1],
    },
    querySelector: (sel) => {
      if (sel === 'td[data-r="0"][data-c="0"]') return cell;
      return null;
    },
  };

  const body = {
    querySelectorAll: (sel) => {
      if (sel.includes('.l3-col')) return [slabCol];
      if (sel.includes('spotlight-cell')) return [];
      return [];
    },
  };

  const leftCluster  = { labels: [0, 1, 2, 0], fixedKLabels: null };
  const rightCluster = { labels: [1, 2, 0, 1], fixedKLabels: null };

  const sandbox = {
    document: { getElementById: (id) => (id === 'l3Body') ? body : null },
    state: {
      data: { n_samples: 4 },
      spotlight: null,
      spotlightTrackedAll: true,
      tracked: new Set([0]),                 // sample 0 is tracked
    },
    getL2Cluster: () => null,
    getSlabClusterAt: (s, e) => (s === 0 ? leftCluster : rightCluster),
    console,
  };
  const ctx = vm.createContext(sandbox);
  const fnSrc = html.match(/(function\s+_applySpotlightHighlights\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  if (fnSrc) {
    vm.runInContext(fnSrc[1], ctx);
    vm.runInContext('_applySpotlightHighlights(null)', ctx);
    ok('tracked-only mode tags cell with spotlight-cell-tracked (not -primary)',
       tagged.includes('spotlight-cell-tracked') && !tagged.includes('spotlight-cell-primary'));
  }
}

// ============================================================================
// 12. Defensive: missing __l3_meta returns early (no crash)
// ============================================================================
console.log('\n=== 12. Defensive paths ===');

{
  const colNoMeta = {
    __l3_meta: null,
    querySelector: () => null,
  };
  const body = {
    querySelectorAll: (sel) => {
      if (sel.includes('.l3-col')) return [colNoMeta];
      if (sel.includes('spotlight-cell')) return [];
      return [];
    },
  };
  const sandbox = {
    document: { getElementById: (id) => (id === 'l3Body') ? body : null },
    state: { data: { n_samples: 4 }, spotlight: 0, spotlightTrackedAll: false, tracked: new Set() },
    getL2Cluster: () => null,
    getSlabClusterAt: () => null,
    console,
  };
  const ctx = vm.createContext(sandbox);
  const fnSrc = html.match(/(function\s+_applySpotlightHighlights\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  if (fnSrc) {
    vm.runInContext(fnSrc[1], ctx);
    let crashed = false;
    try { vm.runInContext('_applySpotlightHighlights(null)', ctx); }
    catch (e) { crashed = true; }
    ok('column without __l3_meta does not crash', !crashed);
  }
}

// Slab column where getSlabClusterAt returns null (cluster failed)
{
  const slabCol = {
    __l3_meta: {
      K: 3, offset: 1, isSlab: true,
      leftRange: [0, 4], rightRange: [5, 9],
      invPerm: [0, 1, 2],
    },
    querySelector: () => null,
  };
  const body = {
    querySelectorAll: (sel) => sel.includes('.l3-col') ? [slabCol] : [],
  };
  const sandbox = {
    document: { getElementById: () => body },
    state: { data: { n_samples: 4 }, spotlight: 0, spotlightTrackedAll: false, tracked: new Set() },
    getL2Cluster: () => null,
    getSlabClusterAt: () => null,                     // simulate cluster failure
    console,
  };
  const ctx = vm.createContext(sandbox);
  const fnSrc = html.match(/(function\s+_applySpotlightHighlights\s*\([^)]*\)\s*\{[\s\S]*?\n\})\s*\n/);
  if (fnSrc) {
    vm.runInContext(fnSrc[1], ctx);
    let crashed = false;
    try { vm.runInContext('_applySpotlightHighlights(null)', ctx); }
    catch (e) { crashed = true; }
    ok('slab cluster lookup returning null does not crash', !crashed);
  }
}

// ============================================================================
// 13. Negative tests — confirm the L2 path was NOT broken by the slab edits
// ============================================================================
console.log('\n=== 13. L2 path preserved ===');

// L2 mode in renderL3Panel still does its own __l3_meta stash with leftIdx /
// rightIdx (NOT leftRange / rightRange).
const l2RendererRe = /function\s+renderL3Panel\s*\(\s*\)\s*\{([\s\S]*?)\n\}\s*\n\s*\/\/[^\n]*={3,}/;
const l2RendererMatch = html.match(l2RendererRe);
ok('renderL3Panel body extracted',
   !!l2RendererMatch);
const l2RendererBody = l2RendererMatch ? l2RendererMatch[1] : '';

ok('L2 mode still uses leftIdx in __l3_meta',
   /col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?leftIdx:/.test(l2RendererBody));

ok('L2 mode still uses rightIdx in __l3_meta',
   /col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?rightIdx:/.test(l2RendererBody));

ok('L2 mode does NOT set isSlab:true (only slab does)',
   !/col\.__l3_meta\s*=\s*\{[\s\S]{0,400}?isSlab:\s*true/.test(l2RendererBody));

ok('L2 mode still calls _applySpotlightHighlights(curL2)',
   /_applySpotlightHighlights\(curL2\)/.test(l2RendererBody));

// Slab branch of compareUnit short-circuit still routes to renderL3PanelSlab
ok('renderL3Panel still short-circuits to renderL3PanelSlab when compareUnit !== "L2"',
   /if \(state\.compareUnit && state\.compareUnit !== 'L2'\)\s*\{\s*return renderL3PanelSlab\(\)/.test(l2RendererBody));

// ============================================================================
// 14. Functionality preservation — pre-existing slab features still wired
// ============================================================================
console.log('\n=== 14. Pre-existing slab features still wired ===');

// Things that existed before this turn and must NOT have been removed:
ok('slab still calls slabFocalContentHtml',
   /content\.innerHTML\s*=\s*slabFocalContentHtml\(cl,\s*range,\s*K\)/.test(slabBody));
ok('slab still calls drawSlabMiniPCA on focal',
   /drawSlabMiniPCA\(mini,\s*range,/.test(slabBody));
ok('slab still calls drawSlabMiniPCA on neighbors',
   /drawSlabMiniPCA\(mini,\s*offsetSlab,/.test(slabBody));
ok('slab still calls compareSlabPair',
   /compareSlabPair\(/.test(slabBody));
ok('slab still calls alignedSlabLabelsTo',
   /alignedSlabLabelsTo\(range,\s*offsetSlab,\s*K\)/.test(slabBody));
ok('slab still calls ctHtml for neighbor contingency',
   /ctHtml\(cmp,\s*offset,\s*alignedLabels\)/.test(slabBody));

// l3KMode parity (was already implemented in turn 55, must remain)
ok('slab still respects state.l3KMode (k3 / k6 / both)',
   /const l3kMode\s*=\s*state\.l3KMode\s*\|\|\s*'k3'/.test(slabBody));
ok('slab still computes ksToRender from kMode',
   /ksToRender = \(l3kMode === 'k6'\) \? \[6\][\s\S]{0,200}?\(l3kMode === 'both'\) \? \[state\.k, 6\]/.test(slabBody));

// ============================================================================
// Summary
// ============================================================================
console.log(`\n=========================================`);
console.log(`  PASS: ${pass}   FAIL: ${fail}`);
console.log(`=========================================\n`);
process.exit(fail > 0 ? 1 : 0);
