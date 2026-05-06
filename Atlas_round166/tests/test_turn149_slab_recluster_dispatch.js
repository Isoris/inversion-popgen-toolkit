// =============================================================================
// turn 149 — compareSlabPair_byMode + slab renderer reMode dispatch
// =============================================================================
// Closes the gap turn 148 introduced: per-pane recluster dropdown was visible
// in slab mode but slab compare silently used kmeans-K3 regardless of
// state.l3ReclusterMode. This turn adds the dispatch + a visible fallback
// notice for U/V modes (which require a slab-aware rotation pipeline that's
// queued but not shipped).
//
// Coverage today:
//   - kmeans-K3 → existing compareSlabPair path, no change for default users
//   - kmeans-K6 → compareSlabPair_byMode → compareSlabPair(left, right, 6)
//   - U/V modes → compareSlabPair_byMode falls back to kmeans-K3 with
//     fellBack:true and requestedMode='uv-...'; renderer prepends a notice
//     div above the contingency table so the user knows the visible result
//     doesn't reflect the dropdown setting
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
// 1. compareSlabPair_byMode function exists with correct signature
// ============================================================================
console.log('\n=== 1. compareSlabPair_byMode signature ===');

ok('compareSlabPair_byMode function declared',
   /function\s+compareSlabPair_byMode\s*\(\s*leftRange\s*,\s*rightRange\s*,\s*mode\s*\)/.test(html));

const fnRe = /function\s+compareSlabPair_byMode\s*\([^)]*\)\s*\{([\s\S]*?)\n\}/;
const fnMatch = html.match(fnRe);
ok('compareSlabPair_byMode body extracted',
   !!fnMatch, 'function boundary regex did not match');
const fnBody = fnMatch ? fnMatch[1] : '';

// ============================================================================
// 2. Mode handling — kmeans-K3, kmeans-K6, U/V real-dispatch
// ============================================================================
console.log('\n=== 2. Mode dispatch within compareSlabPair_byMode ===');

// turn 150: K is now derived in one place at the top of the function
// (kmeans-K6 → K=6 ; everything else → state.k). The dispatch then
// branches: kmeans-K3/K6/null short-circuit through compareSlabPair;
// U/V goes through getSlabClusterByMode.
ok('K derivation: kmeans-K6 → K=6, else state.k (single ternary)',
   /const K\s*=\s*\(m === 'kmeans-K6'\)\s*\?\s*6\s*:\s*state\.k/.test(fnBody));

ok('kmeans-K3 / kmeans-K6 / null path uses compareSlabPair short-circuit',
   /m === 'kmeans-K3'\s*\|\|\s*!m\s*\|\|\s*m === 'kmeans-K6'[\s\S]{0,400}?compareSlabPair\(leftRange,\s*rightRange,\s*K\)/.test(fnBody));

ok('U/V branch fetches via getSlabClusterByMode for both ranges',
   /getSlabClusterByMode\(leftRange\[0\],\s+leftRange\[1\],\s+m\)[\s\S]{0,400}?getSlabClusterByMode\(rightRange\[0\],\s*rightRange\[1\],\s*m\)/.test(fnBody));

ok('U/V branch handles cluster failure → falls back to kmeans-K3 with fellBack:true',
   /if \(!cl\s*\|\|\s*!cr[\s\S]{0,400}?compareSlabPair\(leftRange,\s*rightRange,\s*state\.k\)[\s\S]{0,300}?fellBack\s*=\s*true/.test(fnBody));

ok('U/V branch builds contingency via alignLabels + chiSquare (no kmeans-K3 substitution)',
   /alignLabels\(llab,\s*rlab,\s*K\)[\s\S]{0,400}?chiSquare\(align\.table,\s*K\)/.test(fnBody));

// ============================================================================
// 3. Output shape annotations
// ============================================================================
console.log('\n=== 3. Output shape ===');

// turn 150: kmeans short-circuit always sets fellBack:false explicitly
ok('kmeans short-circuit sets fellBack:false explicitly',
   /reclusterMode\s*=\s*m\s*\|\|\s*'kmeans-K3'[\s\S]{0,200}?fellBack\s*=\s*false/.test(fnBody));

ok('U/V success path sets reclusterMode: m (preserved)',
   /reclusterMode:\s*m,/.test(fnBody));

ok('U/V failure-fallback sets reclusterMode = "kmeans-K3"',
   /reclusterMode\s*=\s*'kmeans-K3'[\s\S]{0,200}?fellBack\s*=\s*true/.test(fnBody));

ok('returns null when leftRange/rightRange missing',
   /if \(!leftRange\s*\|\|\s*!rightRange\) return null/.test(fnBody));

// ============================================================================
// 4. renderL3PanelSlab dispatches to compareSlabPair_byMode
// ============================================================================
console.log('\n=== 4. Slab renderer dispatch ===');

const slabFnRe = /function\s+renderL3PanelSlab\s*\(\s*\)\s*\{([\s\S]*?)\n\}\s*\n\s*\/\/[^\n]*Compact focal-content HTML for a slab/;
const slabFnMatch = html.match(slabFnRe);
ok('renderL3PanelSlab body extracted',
   !!slabFnMatch);
const slabBody = slabFnMatch ? slabFnMatch[1] : '';

ok('slab renderer reads state.l3ReclusterMode',
   /const reMode\s*=\s*state\.l3ReclusterMode\s*\|\|\s*'kmeans-K3'/.test(slabBody));

ok('slab renderer dispatches to compareSlabPair_byMode for non-default modes',
   /if \(reMode !== 'kmeans-K3'\)\s*\{[\s\S]{0,400}?compareSlabPair_byMode\(/.test(slabBody));

ok('slab renderer keeps default kmeans-K3 path on plain compareSlabPair',
   /\}\s*else\s*\{[\s\S]{0,400}?compareSlabPair\(offsetSlab,\s*range,\s*K\)/.test(slabBody));

// Both offset arms (negative and positive) covered
ok('byMode dispatch handles offset < 0 (offsetSlab on left, range on right)',
   /reMode !== 'kmeans-K3'[\s\S]{0,400}?\(offset\s*<\s*0\)[\s\S]{0,200}?compareSlabPair_byMode\(offsetSlab,\s*range,\s*reMode\)/.test(slabBody));

ok('byMode dispatch handles offset > 0 (range on left, offsetSlab on right)',
   /reMode !== 'kmeans-K3'[\s\S]{0,400}?compareSlabPair_byMode\(range,\s*offsetSlab,\s*reMode\)/.test(slabBody));

// ============================================================================
// 5. Fallback notice rendering
// ============================================================================
console.log('\n=== 5. Fallback notice in contingency content ===');

ok('renderer checks cmp.fellBack to surface notice',
   /if \(cmp\.fellBack/.test(slabBody));

ok('renderer reads cmp.requestedMode for the notice',
   /cmp\.requestedMode/.test(slabBody));

ok('notice mentions "slab fallback" so it is searchable',
   /slab fallback/.test(slabBody));

ok('notice is prepended to ctHtml output (not replacing it)',
   /content\.innerHTML\s*=\s*prefix\s*\+\s*ctHtml\(cmp,\s*offset,\s*alignedLabels\)/.test(slabBody));

ok('notice has descriptive title attribute (educates the user)',
   /title="The U\/V recluster modes/.test(slabBody));

ok('notice uses the accent border to flag fallback state',
   /border-left:\s*2px solid var\(--accent\)/.test(slabBody));

// ============================================================================
// 6. col.__l3_meta carries the contingency K (not section K) when reMode used
// ============================================================================
console.log('\n=== 6. __l3_meta uses contingency K (ctK) ===');

ok('renderer reads ctK from cmp.K (with fallback to section K)',
   /const ctK\s*=\s*\(cmp\s*&&\s*cmp\.K\)\s*\?\s*cmp\.K\s*:\s*K/.test(slabBody));

ok('invPerm sized by ctK',
   /const _invPerm\s*=\s*new Array\(ctK\)/.test(slabBody));

ok('__l3_meta.K set to ctK',
   /col\.__l3_meta\s*=\s*\{\s*K:\s*ctK/.test(slabBody));

// ============================================================================
// 7. Sandboxed run of compareSlabPair_byMode against mock environment
// ============================================================================
console.log('\n=== 7. Sandboxed compareSlabPair_byMode run ===');

{
  // Build a sandbox with the required helpers stubbed. compareSlabPair is
  // stubbed to return a recognizable shape for each (left, right, K) call so
  // we can verify which K the byMode dispatcher passed.
  const _calls = [];
  function makeCmpStub() {
    return (lr, rr, K) => {
      _calls.push({ lr: lr.slice(), rr: rr.slice(), K });
      return {
        leftRange: lr, rightRange: rr, K,
        table: Array.from({length: K}, () => new Array(K).fill(0)),
        perm: Array.from({length: K}, (_, i) => i),
        concord: 1.0,
        verdict: 'MERGE',
        cl_ok: true, cr_ok: true,
      };
    };
  }

  // Extract compareSlabPair_byMode source
  const fnSrcMatch = html.match(/(function\s+compareSlabPair_byMode\s*\([^)]*\)\s*\{[\s\S]*?\n\})/);
  ok('extracted compareSlabPair_byMode for sandbox',
     !!fnSrcMatch);

  if (fnSrcMatch) {
    // turn 150: compareSlabPair_byMode now calls getSlabClusterByMode for U/V
    // modes (no longer falls back). The sandbox needs to stub it. We also need
    // alignLabels, fisher2x2, chiSquare for the U/V code path that builds
    // the contingency directly.

    // Stubbed cluster output for U/V modes
    function makeSlabClusterStub(K) {
      const labels = new Int8Array(10);
      for (let i = 0; i < 10; i++) labels[i] = i % K;
      const npg = new Array(K).fill(0);
      for (let i = 0; i < 10; i++) npg[labels[i]]++;
      return {
        ok: true, reason: null, labels, fixedKLabels: labels,
        n_per_group: npg, usedK: K,
      };
    }

    const _modeCalls = [];
    const sandbox = {
      state: { k: 3, mergeThr: 0.85, alpha: 0.001 },
      compareSlabPair: makeCmpStub(),
      getSlabClusterByMode: (s, e, mode) => {
        _modeCalls.push({ s, e, mode });
        return makeSlabClusterStub(3);
      },
      alignLabels: (l, r, K) => {
        const table = Array.from({length: K}, () => new Array(K).fill(0));
        for (let i = 0; i < l.length; i++) table[l[i]][r[i]]++;
        return {
          table,
          perm: Array.from({length: K}, (_, i) => i),
          aligned: r,
          concord: 1.0,
        };
      },
      fisher2x2: () => 0.5,
      chiSquare: () => ({ p_approx: 0.5 }),
      console,
    };
    const ctx = vm.createContext(sandbox);
    vm.runInContext(fnSrcMatch[1], ctx);

    // Test 1: kmeans-K3 → K=3, fellBack=false
    _calls.length = 0; _modeCalls.length = 0;
    let r1 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], 'kmeans-K3')`, ctx);
    ok('kmeans-K3 → fellBack:false',
       r1.fellBack === false);
    ok('kmeans-K3 → routes through compareSlabPair (not getSlabClusterByMode)',
       _calls.length === 1 && _modeCalls.length === 0 && _calls[0].K === 3);
    ok('kmeans-K3 → reclusterMode preserved as kmeans-K3',
       r1.reclusterMode === 'kmeans-K3');
    ok('kmeans-K3 → requestedMode reflects input',
       r1.requestedMode === 'kmeans-K3');

    // Test 2: kmeans-K6 → K=6, fellBack=false
    _calls.length = 0; _modeCalls.length = 0;
    let r2 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], 'kmeans-K6')`, ctx);
    ok('kmeans-K6 → fellBack:false',
       r2.fellBack === false);
    ok('kmeans-K6 → routes through compareSlabPair at K=6',
       _calls.length === 1 && _calls[0].K === 6);
    ok('kmeans-K6 → reclusterMode = kmeans-K6 (preserved)',
       r2.reclusterMode === 'kmeans-K6');

    // Test 3: U/V modes — turn 150 changes from fellBack:true to fellBack:false.
    // U/V modes now go through getSlabClusterByMode and produce real
    // U/V-clustered contingency tables. Two assertions inverted from turn 149.
    for (const uvMode of ['uv-rotated', 'uv-denoise', 'uv-dbscan',
                          'uv-dist-rank', 'uv-dist-fuzzy', 'distance-uv']) {
      _calls.length = 0; _modeCalls.length = 0;
      let r = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], '${uvMode}')`, ctx);
      ok(`${uvMode} → fellBack:false (turn 150: actually runs the U/V cluster)`,
         r && r.fellBack === false);
      ok(`${uvMode} → routes through getSlabClusterByMode (not fallback)`,
         _modeCalls.length === 2 &&
         _modeCalls[0].mode === uvMode &&
         _modeCalls[1].mode === uvMode);
      ok(`${uvMode} → reclusterMode preserved (no longer collapsed to kmeans-K3)`,
         r && r.reclusterMode === uvMode);
      ok(`${uvMode} → requestedMode preserved as '${uvMode}'`,
         r && r.requestedMode === uvMode);
    }

    // Test 4: undefined / null mode → defaults to kmeans-K3
    _calls.length = 0; _modeCalls.length = 0;
    let r4 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], null)`, ctx);
    ok('null mode → K=3 default, fellBack:false',
       r4 && r4.fellBack === false && _calls[0].K === 3);

    // Test 5: missing ranges → null
    let r5 = vm.runInContext(`compareSlabPair_byMode(null, [5,9], 'kmeans-K3')`, ctx);
    ok('missing leftRange → returns null',
       r5 === null);

    // Test 6: compareSlabPair returns null on kmeans-K3 → byMode returns null
    sandbox.compareSlabPair = () => null;
    let r6 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], 'kmeans-K3')`, ctx);
    ok('inner compareSlabPair null → byMode returns null',
       r6 === null);

    // Test 7: U/V mode with cluster lookup failing → falls back to kmeans-K3
    // with fellBack:true. This tests the genuine-failure fallback path.
    sandbox.compareSlabPair = makeCmpStub();
    sandbox.getSlabClusterByMode = () => null;     // simulate rotation failure
    _calls.length = 0;
    let r7 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], 'uv-rotated')`, ctx);
    ok('U/V cluster failure → falls back to kmeans-K3 with fellBack:true',
       r7 && r7.fellBack === true && r7.reclusterMode === 'kmeans-K3');
    ok('U/V cluster failure → requestedMode still preserved',
       r7.requestedMode === 'uv-rotated');

    // Test 8: unknown future mode (not kmeans-K3, kmeans-K6, or known U/V) —
    // turn 150: routes through getSlabClusterByMode default branch (which
    // falls back to kmeans-K3 internally), so byMode returns the U/V-style
    // contingency, NOT the kmeans-K3 short-circuit. Whether fellBack is true
    // depends on whether getSlabClusterByMode returns success — it does in
    // our stub, so fellBack is false. The mode plumbing is now generic.
    sandbox.compareSlabPair = makeCmpStub();
    sandbox.getSlabClusterByMode = (s, e, mode) => {
      _modeCalls.push({ s, e, mode });
      return makeSlabClusterStub(3);
    };
    _calls.length = 0; _modeCalls.length = 0;
    let r8 = vm.runInContext(`compareSlabPair_byMode([0,4], [5,9], 'gmm-3')`, ctx);
    ok('unknown mode "gmm-3" → routes through getSlabClusterByMode (mode dispatch generic)',
       r8 && _modeCalls.length === 2);
    ok('unknown mode preserves requestedMode',
       r8 && r8.requestedMode === 'gmm-3');
  }
}

// ============================================================================
// 8. Integration — verify renderer wiring still extracts cleanly
// ============================================================================
console.log('\n=== 8. Renderer integration sanity ===');

// The renderer's reMode block should be a coherent if/else. Check that the
// branches don't get confused by our regex-based asserts.
ok('renderer has matched if/else branches (one for non-default, one for default)',
   /if \(reMode !== 'kmeans-K3'\)\s*\{[\s\S]{0,500}?\}\s*else\s*\{[\s\S]{0,500}?compareSlabPair\(offsetSlab,\s*range,\s*K\)[\s\S]{0,200}?compareSlabPair\(range,\s*offsetSlab,\s*K\)/.test(slabBody));

// The neighbor mini-PCA still uses alignedSlabLabelsTo with the section K
// (because the mini-PCA's palette is tied to the section, not the
// contingency). When reMode forces a different K, the mini-PCA stays at
// section K and only the contingency changes — same trade-off L2 mode accepts.
ok('neighbor mini-PCA still uses section-K aligned labels',
   /alignedSlabLabelsTo\(range,\s*offsetSlab,\s*K\)/.test(slabBody));

// ============================================================================
// 9. L2 path unchanged (turn 148 negative test still holds)
// ============================================================================
console.log('\n=== 9. L2 path preserved ===');

ok('compareL2Pair_byMode still exists (L2 mode unchanged)',
   /function\s+compareL2Pair_byMode\s*\(\s*leftIdx\s*,\s*rightIdx\s*,\s*mode\s*\)/.test(html));

ok('L2 mode renderer still uses compareL2Pair_byMode (not compareSlabPair_byMode)',
   /compareL2Pair_byMode\(l2idx,\s*curL2,\s*reMode\)|compareL2Pair_byMode\(curL2,\s*l2idx,\s*reMode\)/.test(html));

// ============================================================================
// 10. Source-level tests for the new function's location and structure
// ============================================================================
console.log('\n=== 10. Code organization ===');

// compareSlabPair_byMode should sit right after compareSlabPair so they're
// findable together. Pattern: end of compareSlabPair → blank/comment lines
// → compareSlabPair_byMode opener.
ok('compareSlabPair_byMode placed near compareSlabPair',
   /function\s+compareSlabPair\s*\([\s\S]{0,3000}?function\s+compareSlabPair_byMode/.test(html));

ok('compareSlabPair_byMode comment block explains the U/V fallback',
   /turn 149 — compareSlabPair_byMode[\s\S]{0,2000}?U\/V[\s\S]{0,500}?fall back to kmeans-K3/.test(html));

// ============================================================================
// Summary
// ============================================================================
console.log(`\n=========================================`);
console.log(`  PASS: ${pass}   FAIL: ${fail}`);
console.log(`=========================================\n`);
process.exit(fail > 0 ? 1 : 0);
