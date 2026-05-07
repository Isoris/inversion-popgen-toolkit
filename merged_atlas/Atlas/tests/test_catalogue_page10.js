// tests/test_catalogue_page10.js
// Smoke + behavioural tests for inversion_catalogue/page10.js (Marker panels).

import * as page10 from '../inversion_catalogue/page10.js';

let pass = 0, fail = 0;
function check(name, ok, detail) {
  if (ok) { console.log('  ✓', name); pass++; }
  else    { console.log('  ✗', name, detail ? '— ' + detail : ''); fail++; }
}

// 1. Module shape
check('exports wirePage10', typeof page10.wirePage10 === 'function');
check('exports default', typeof page10.default === 'function');

// 2. Wiring returns expected handles
const handle = page10.wirePage10({ data: {}, candidateList: [] });
check('wirePage10 returns object', handle && typeof handle === 'object');
check('handle has renderPage10', typeof handle.renderPage10 === 'function');
check('handle has renderMarkerPage alias', typeof handle.renderMarkerPage === 'function');
check('renderPage10 === renderMarkerPage (same fn)', handle.renderPage10 === handle.renderMarkerPage);

// 3. Empty-layers path renders empty state
{
  const slot = { innerHTML: '' };
  const subtitle = { textContent: '' };
  globalThis.document = {
    getElementById: id => {
      if (id === 'page10Content') return slot;
      if (id === 'page10Subtitle') return subtitle;
      return null;
    },
  };
  const state = { data: { _layers_present: [] }, candidateList: [] };
  const { renderPage10 } = page10.wirePage10(state);
  renderPage10();
  check('empty-layers subtitle set', subtitle.textContent === '(no marker layer loaded)');
  check('empty-layers shows guidance markup', slot.innerHTML.includes('No marker panels loaded'));
}

// 4. Missing slot DOM is tolerated (early return)
{
  globalThis.document = { getElementById: () => null };
  const { renderPage10 } = page10.wirePage10({ data: {}, candidateList: [] });
  let threw = false;
  try { renderPage10(); } catch (e) { threw = true; }
  check('renderPage10 with missing DOM does not throw', !threw);
}

// 5. Layer present but zero summaries — second empty-state branch
{
  const slot = { innerHTML: '' };
  const subtitle = { textContent: '' };
  globalThis.document = {
    getElementById: id => id === 'page10Content' ? slot : id === 'page10Subtitle' ? subtitle : null,
  };
  const state = {
    data: {
      _layers_present: ['marker_panel_summary'],
      marker_panel_summary: [],
      marker_catalogue: [],
      marker_primers: [],
    },
    candidateList: [],
  };
  const { renderPage10 } = page10.wirePage10(state);
  renderPage10();
  check('zero-summaries empty state rendered',
        slot.innerHTML.includes('Marker layer loaded but no panels emitted'));
}

console.log(`pass: ${pass}   fail: ${fail}`);
process.exit(fail === 0 ? 0 : 1);
