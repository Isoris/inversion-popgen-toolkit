// Atlas/tests/test_discovery_page8.js
//
// Minimal smoke test — verify page8.js loads as ESM and exports the expected names.
// Per HANDOFF_BATCH_1: extracted bodies are not exercised here (they need DOM + state).
// The merge chat will add behavioural tests once the wiring lands.
//

import * as page8 from '../inversion_discovery/page8.js';

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log("--- page8 module loads ---");
check("module imports as object", typeof page8 === "object" && page8 !== null);

// page has no JS handlers in legacy — pure HTML scaffold;
// merge chat / future batches will author them fresh.
check("no exports yet (scaffold only)", Object.keys(page8).length === 0);

console.log("");
console.log("=================");
console.log(`pass: ${pass}   fail: ${fail}`);
console.log("=================");
process.exit(fail === 0 ? 0 : 1);
