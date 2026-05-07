// Atlas/tests/test_discovery_page1.js
//
// Minimal smoke test — verify page1.js loads as ESM and exports the expected names.
// Per HANDOFF_BATCH_1: extracted bodies are not exercised here (they need DOM + state).
// The merge chat will add behavioural tests once the wiring lands.
//

import * as page1 from '../inversion_discovery/page1.js';

let pass = 0, fail = 0;
function check(name, cond, detail = "") {
  if (cond) { console.log(`  ✓ ${name}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log("--- page1 module loads ---");
check("module imports as object", typeof page1 === "object" && page1 !== null);

// 26 expected exports
check("drawSim exported as function", typeof page1.drawSim === "function");
check("drawSimMini exported as function", typeof page1.drawSimMini === "function");
check("drawZ exported as function", typeof page1.drawZ === "function");
check("drawLinesPanel exported as function", typeof page1.drawLinesPanel === "function");
check("drawPCA exported as function", typeof page1.drawPCA === "function");
check("drawAnchorStrip exported as function", typeof page1.drawAnchorStrip === "function");
check("renderL3Panel exported as function", typeof page1.renderL3Panel === "function");
check("renderL3PanelSlab exported as function", typeof page1.renderL3PanelSlab === "function");
check("renderL3PanelScaleStability exported as function", typeof page1.renderL3PanelScaleStability === "function");
check("updateWinLabel exported as function", typeof page1.updateWinLabel === "function");
check("setCur exported as function", typeof page1.setCur === "function");
check("autoPickRadial exported as function", typeof page1.autoPickRadial === "function");
check("applyData exported as function", typeof page1.applyData === "function");
check("onSimClick exported as function", typeof page1.onSimClick === "function");
check("onZClick exported as function", typeof page1.onZClick === "function");
check("onPCAClick exported as function", typeof page1.onPCAClick === "function");
check("togglePlay exported as function", typeof page1.togglePlay === "function");
check("cycleKAside exported as function", typeof page1.cycleKAside === "function");
check("renderTrackedList exported as function", typeof page1.renderTrackedList === "function");
check("renderManualGroupsList exported as function", typeof page1.renderManualGroupsList === "function");
check("buildLinesPanelCheckboxes exported as function", typeof page1.buildLinesPanelCheckboxes === "function");
check("buildLinesPanel exported as function", typeof page1.buildLinesPanel === "function");
check("buildTrackPanels exported as function", typeof page1.buildTrackPanels === "function");
check("drawTracks exported as function", typeof page1.drawTracks === "function");
check("refreshLinesColorMode exported as function", typeof page1.refreshLinesColorMode === "function");
check("setLinesPanelCandidateBands exported as function", typeof page1.setLinesPanelCandidateBands === "function");

console.log("");
console.log("=================");
console.log(`pass: ${pass}   fail: ${fail}`);
console.log("=================");
process.exit(fail === 0 ? 0 : 1);
