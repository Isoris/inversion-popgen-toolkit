// Atlas/tests/test_discovery_page15.js
//
// Minimal smoke test — verify page15.js loads as ESM and exports the expected names.
// Per HANDOFF_BATCH_1: extracted bodies are not exercised here (they need DOM + state).
// The merge chat will add behavioural tests once the wiring lands.
//

import * as page15 from '../inversion_discovery/page15.js';

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log("--- page15 module loads ---");
check("module imports as object", typeof page15 === "object" && page15 !== null);

// 1 expected exports
check("_refreshGhslLayerStatus exported as function", typeof page15._refreshGhslLayerStatus === "function");

console.log("");
console.log("=================");
console.log(`pass: ${pass}   fail: ${fail}`);
console.log("=================");
process.exit(fail === 0 ? 0 : 1);
