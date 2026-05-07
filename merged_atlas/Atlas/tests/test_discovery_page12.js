// Atlas/tests/test_discovery_page12.js
//
// Minimal smoke test — verify page12.js loads as ESM and exports the expected names.
// Per HANDOFF_BATCH_1: extracted bodies are not exercised here (they need DOM + state).
// The merge chat will add behavioural tests once the wiring lands.
//

import * as page12 from '../inversion_discovery/page12.js';

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log("--- page12 module loads ---");
check("module imports as object", typeof page12 === "object" && page12 !== null);

// 8 expected exports
check("_refreshThetaPiLayerStatus exported as function", typeof page12._refreshThetaPiLayerStatus === "function");
check("_refreshThetaPiPanelVisibility exported as function", typeof page12._refreshThetaPiPanelVisibility === "function");
check("_drawThCusumHero exported as function", typeof page12._drawThCusumHero === "function");
check("_drawThLinesPanel exported as function", typeof page12._drawThLinesPanel === "function");
check("_drawThSimMatPanel exported as function", typeof page12._drawThSimMatPanel === "function");
check("_drawThZPanel exported as function", typeof page12._drawThZPanel === "function");
check("_drawThAnchorStripPanel exported as function", typeof page12._drawThAnchorStripPanel === "function");
check("_drawThPcaPanel exported as function", typeof page12._drawThPcaPanel === "function");

console.log("");
console.log("=================");
console.log(`pass: ${pass}   fail: ${fail}`);
console.log("=================");
process.exit(fail === 0 ? 0 : 1);
