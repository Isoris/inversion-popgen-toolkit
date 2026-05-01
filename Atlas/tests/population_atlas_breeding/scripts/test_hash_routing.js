// Hash routing test (v4 turn 79).
// Verifies that:
//   1. Loading Inversion_atlas.html#page18 lands on the marker panel tab,
//      not whatever was last in localStorage.
//   2. Loading Population_atlas.html#page8 lands on the breeding tab.
//   3. Bogus hashes fall back to localStorage / first tab.
//   4. hashchange events fire correctly for in-page navigation.
//   5. The marker panel methods note now includes a link to the
//      Population Atlas breeding page.

const { JSDOM } = require('jsdom');
const fs = require('fs');

function fail(suite, msg) { console.error('[' + suite + '] FAIL:', msg); process.exit(1); }
function ok(suite, msg)   { console.log ('[' + suite + '] ok  :', msg); }

// =====================================================================
// Suite 1: Inversion Atlas hash routing
// =====================================================================
function testInversionAtlas() {
  const SUITE = 'inv-hash';
  const html = fs.readFileSync('/home/claude/work/Inversion_atlas.html', 'utf-8');

  // Test 1A: load with #page18 → marker panel tab activates
  const dom1 = new JSDOM(html, {
    runScripts: 'dangerously',
    pretendToBeVisual: true,
    url: 'http://localhost/Inversion_atlas.html#page18',
  });
  const w1 = dom1.window;
  // Wait for the microtask that runs hash routing to fire
  return new Promise(resolve => {
    setTimeout(() => {
      const activePage = w1.document.querySelector('.page.active');
      if (!activePage) fail(SUITE, 'no active page after hash routing');
      if (activePage.id !== 'page18') fail(SUITE, 'expected #page18 active, got #' + activePage.id);
      ok(SUITE, '#page18 hash → marker panel tab active');

      const activeTab = w1.document.querySelector('#tabBar button.active');
      if (!activeTab) fail(SUITE, 'no active tab button');
      if (activeTab.dataset.page !== 'page18') fail(SUITE, 'expected page18 tab active, got ' + activeTab.dataset.page);
      ok(SUITE, 'tab button matches: ' + (activeTab.textContent || '').trim().replace(/\s+/g,' '));

      // Test 1B: hashchange triggers re-routing
      w1.location.hash = '#page17';
      // dispatch hashchange manually since JSDOM doesn't auto-fire on .hash setter
      w1.dispatchEvent(new w1.Event('hashchange'));
      setTimeout(() => {
        const activePageAfter = w1.document.querySelector('.page.active');
        if (activePageAfter.id !== 'page17') fail(SUITE, 'hashchange to #page17 failed, got #' + activePageAfter.id);
        ok(SUITE, 'hashchange #page17 → stats profile tab activates');

        // Test 1C: bogus hash gracefully ignored (does not navigate away from current)
        const beforeBogus = w1.document.querySelector('.page.active').id;
        w1.location.hash = '#nonexistent_page';
        w1.dispatchEvent(new w1.Event('hashchange'));
        setTimeout(() => {
          const afterBogus = w1.document.querySelector('.page.active').id;
          if (afterBogus !== beforeBogus) fail(SUITE, 'bogus hash should not change tab, was #' + beforeBogus + ' now #' + afterBogus);
          ok(SUITE, 'bogus hash ignored, stays on #' + afterBogus);

          // Test 1D: marker panel methods note includes Population Atlas back-link
          // Activate page18 first by dispatching click
          const page18Btn = w1.document.querySelector('#tabBar button[data-page="page18"]');
          page18Btn.click();
          setTimeout(() => {
            const mpBody = w1.document.getElementById('mpBody');
            if (!mpBody) fail(SUITE, 'mpBody slot missing');
            if (!mpBody.innerHTML.includes('Population_atlas.html#page8')) {
              fail(SUITE, 'methods note missing Population Atlas back-link');
            }
            if (!mpBody.innerHTML.includes('Population Atlas \u2192 Breeding')) {
              fail(SUITE, 'methods note missing "Population Atlas → Breeding" text');
            }
            ok(SUITE, 'marker panel methods note includes Population Atlas → Breeding back-link');

            resolve();
          }, 30);
        }, 30);
      }, 30);
    }, 30);
  });
}

// =====================================================================
// Suite 2: Population Atlas hash routing
// =====================================================================
function testPopulationAtlas() {
  const SUITE = 'pop-hash';
  return new Promise(resolve => {
    const html = fs.readFileSync('/home/claude/work/Population_atlas.html', 'utf-8');

    // Test 2A: load with #page8 → breeding tab activates
    const dom = new JSDOM(html, {
      runScripts: 'dangerously',
      pretendToBeVisual: true,
      url: 'http://localhost/Population_atlas.html#page8',
    });
    const w = dom.window;
    setTimeout(() => {
      const activePage = w.document.querySelector('.page.active');
      if (!activePage) fail(SUITE, 'no active page');
      if (activePage.id !== 'page8') fail(SUITE, 'expected #page8 active, got #' + activePage.id);
      ok(SUITE, '#page8 hash → breeding tab active');

      // Test 2B: bp-empty visible (since no upstream data) — confirms render fired
      const slot = w.document.getElementById('bpTableSlot');
      if (!slot.innerHTML.includes('No rows yet')) fail(SUITE, 'breeding page did not render on hash-route');
      ok(SUITE, 'breeding page renders empty state correctly on hash-route entry');

      // Test 2C: hashchange to #page3 (families) works
      w.location.hash = '#page3';
      w.dispatchEvent(new w.Event('hashchange'));
      setTimeout(() => {
        const activeAfter = w.document.querySelector('.page.active');
        if (activeAfter.id !== 'page3') fail(SUITE, 'hashchange to #page3 failed');
        ok(SUITE, 'hashchange #page3 → families & clusters tab activates');

        // Test 2D: bogus hash ignored
        w.location.hash = '#bogus';
        w.dispatchEvent(new w.Event('hashchange'));
        setTimeout(() => {
          const after = w.document.querySelector('.page.active');
          if (after.id !== 'page3') fail(SUITE, 'bogus hash should not change tab');
          ok(SUITE, 'bogus hash ignored, stays on #page3');

          resolve();
        }, 30);
      }, 30);
    }, 30);
  });
}

// =====================================================================
// Suite 3: empty-hash fallback (no hash, fresh load → page1 visible)
// =====================================================================
function testEmptyHashFallback() {
  const SUITE = 'empty-hash';
  return new Promise(resolve => {
    // Fresh DOM with no hash and no localStorage
    const popHtml = fs.readFileSync('/home/claude/work/Population_atlas.html', 'utf-8');
    const dom = new JSDOM(popHtml, {
      runScripts: 'dangerously',
      pretendToBeVisual: true,
      url: 'http://localhost/Population_atlas.html',
    });
    const w = dom.window;
    setTimeout(() => {
      const active = w.document.querySelector('.page.active');
      if (!active) fail(SUITE, 'no active page on empty-hash load');
      // Should default to page1 (which has class="page active" baked into HTML)
      if (active.id !== 'page1') fail(SUITE, 'empty-hash should default to page1, got #' + active.id);
      ok(SUITE, 'empty hash + no localStorage → defaults to page1 (samples)');

      resolve();
    }, 30);
  });
}

// Run all suites sequentially
(async () => {
  await testInversionAtlas();
  await testPopulationAtlas();
  await testEmptyHashFallback();
  console.log('\n[hash-routing] ALL SUITES PASSED');
})().catch(err => {
  console.error('Unhandled:', err);
  process.exit(1);
});
