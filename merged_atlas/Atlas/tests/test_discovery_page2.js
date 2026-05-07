// Atlas/tests/test_discovery_page2.js
//
// Minimal smoke test — verify page2.js loads as ESM and exports the expected names.
// Per HANDOFF_BATCH_1: extracted bodies are not exercised here (they need DOM + state).
// The merge chat will add behavioural tests once the wiring lands.
//

import * as page2 from '../inversion_discovery/page2.js';

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log("--- page2 module loads ---");
check("module imports as object", typeof page2 === "object" && page2 !== null);

// 2 expected exports
check("renderCandidateMetadata exported as function", typeof page2.renderCandidateMetadata === "function");
check("wireCandidateNav exported as function", typeof page2.wireCandidateNav === "function");

console.log("");
console.log("=================");
console.log(`pass: ${pass}   fail: ${fail}`);
console.log("=================");
process.exit(fail === 0 ? 0 : 1);
