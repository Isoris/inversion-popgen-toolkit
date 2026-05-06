// =============================================================================
// turn 157 — dosage bridge auto-load + popstats/ancestry tab guard
// =============================================================================
// Quentin reported:
//   1. "dosage.js shim script is not found"
//   2. "click ancestry/popstats tab → can't click other tabs, doesn't refresh"
//   3. "popstats server connected but tracks not computed"
//
// Root cause for (1) + (3): js/atlas_dosage_bridge.js exists on disk but is
// NOT in the <script src> load list. The bridge's atlasInstallServerDosageBridge
// rewrites state.data.dosage_chunks to point at the live server's
// /api/dosage/chunk endpoint with templated URLs. Without it, dosage chunks
// remain whatever static index was loaded (typically nothing for live
// sessions), so the heatmap, stripe-quality, and the turn-156 V-shape
// diagnostic all stall.
//
// This turn ships:
//   1. <script src="js/atlas_dosage_bridge.js"> — the missing load tag.
//   2. _wireDosageBridgeAutoInstall — patches atlasServer.isAvailable so
//      every positive health probe triggers an idempotent
//      AtlasDosageBridge.install() with the active stateRef + atlasServer.
//      Failures are silent (server up but dosage subsystem disabled is a
//      legal partial state per turn 145).
//   3. try/catch around the page6/page7 RAF render entry-points so a
//      thrown error inside renderPopstatsPage / renderAncestryPage can't
//      starve the click event loop on subsequent tab clicks.
//
// Issues 4-8 from Quentin's message are real but each is its own turn.
// Captured in HANDOFF §8.
// =============================================================================

const fs = require('fs');
const path = require('path');

const ATLAS_PATH = path.resolve(__dirname, '..', 'Inversion_atlas.html');
const html = fs.readFileSync(ATLAS_PATH, 'utf8');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS ' + name); }
  else { fail++; console.log('  FAIL ' + name + (detail ? ' :: ' + detail : '')); }
}

// ============================================================================
// 1. Source-pattern: dosage bridge script loaded
// ============================================================================
console.log('\n=== 1. Dosage bridge script loaded ===');

ok('js/atlas_dosage_bridge.js script tag present',
   /<script\s+src="js\/atlas_dosage_bridge\.js"/.test(html));

// Ordering: dosage bridge must load AFTER atlas_request_layer (which
// it could in principle depend on), and the auto-install <script> must
// load AFTER the bridge UMD. Here we just verify it's after the other
// js/* loads — anywhere in the trailing block is fine.
const bridgeIdx   = html.indexOf('js/atlas_dosage_bridge.js');
const requestIdx  = html.indexOf('js/atlas_request_layer.js');
ok('dosage bridge loaded after request layer',
   bridgeIdx > requestIdx);

// ============================================================================
// 2. Source-pattern: auto-install hook
// ============================================================================
console.log('\n=== 2. Auto-install hook ===');

ok('_wireDosageBridgeAutoInstall IIFE present',
   /_wireDosageBridgeAutoInstall/.test(html));

ok('hook patches atlasServer.isAvailable',
   /window\.atlasServer\.isAvailable\s*=/.test(html));

ok('hook is idempotent via _dosageBridgeWired flag',
   /atlasServer\._dosageBridgeWired/.test(html));

ok('hook calls AtlasDosageBridge.install with stateRef + atlasServer',
   /AtlasDosageBridge\.install\s*\(\s*\{[\s\S]*?atlasServer:\s*window\.atlasServer[\s\S]*?stateRef:\s*window\.state[\s\S]*?\}\s*\)/.test(html));

ok('hook only installs on ok=true (server reachable)',
   /if\s*\(\s*ok\s*&&\s*window\.AtlasDosageBridge/.test(html));

ok('hook catches install promise rejections silently',
   /\.catch\s*\(\s*\(\s*e\s*\)\s*=>\s*\{[\s\S]{0,400}?console\.warn[\s\S]{0,200}?dosage bridge install failed/.test(html));

ok('hook catches synchronous install throws (defensive)',
   /try\s*\{[^}]*AtlasDosageBridge\.install/.test(html) &&
   /\}\s*catch\s*\(\s*e\s*\)\s*\{[\s\S]*?console\.warn\s*\(\s*'\[turn 157\]\s*dosage bridge install threw/.test(html));

// ============================================================================
// 3. Source-pattern: page6/page7 try/catch
// ============================================================================
console.log('\n=== 3. page6/page7 RAF render guarded ===');

// Find the click-handler section
const clickRe = /_tabBtns\.forEach\(\s*btn\s*=>\s*\{[\s\S]*?if\s*\(\s*target\s*===\s*'page6'/;
ok('click handler block extracted', clickRe.test(html));

// page6 wrapped in try/catch
ok('page6 RAF wrapped in try/catch with console.warn',
   /target\s*===\s*'page6'[\s\S]{0,400}?try\s*\{\s*renderPopstatsPage\s*\(\s*\)\s*;?\s*\}[\s\S]{0,200}?catch\s*\(\s*e\s*\)\s*\{[\s\S]{0,200}?renderPopstatsPage threw/.test(html));

ok('page7 RAF wrapped in try/catch with console.warn',
   /target\s*===\s*'page7'[\s\S]{0,400}?try\s*\{\s*renderAncestryPage\s*\(\s*\)\s*;?\s*\}[\s\S]{0,200}?catch\s*\(\s*e\s*\)\s*\{[\s\S]{0,200}?renderAncestryPage threw/.test(html));

// Note: the original RAF-scheduled call still happens (turn 157 only adds
// try/catch INSIDE the RAF callback; click handler itself is unchanged).
ok('page6 RAF schedule pattern preserved (still uses requestAnimationFrame)',
   /'page6'[\s\S]{0,400}?requestAnimationFrame/.test(html));

ok('page7 RAF schedule pattern preserved',
   /'page7'[\s\S]{0,400}?requestAnimationFrame/.test(html));

// ============================================================================
// 4. Existing flow preserved
// ============================================================================
console.log('\n=== 4. Existing flow preserved ===');

ok('atlasServer object still defined',
   /const\s+atlasServer\s*=\s*\{/.test(html));

ok('atlasServer.isAvailable function still defined',
   /async\s+isAvailable\s*\(\s*forceRefresh\s*\)/.test(html));

ok('renderPopstatsPage still defined',
   /^function\s+renderPopstatsPage\s*\(\s*\)/m.test(html));

ok('renderAncestryPage still defined',
   /^function\s+renderAncestryPage\s*\(\s*\)/m.test(html));

ok('_resolveChunkUrl (atlas-side dosage chunk URL substitution) still defined',
   /function\s+_resolveChunkUrl\s*\(\s*chunkRef\s*,\s*region\s*\)/.test(html));

// Turn 156 V-shape diagnostic must still be present
ok('turn 156 _buildVShapeData still defined',
   /function\s+_buildVShapeData\s*\(\s*candidate\s*,\s*chunk\s*,\s*sqRows\s*\)/.test(html));

// Turn 155 cache-key signature still extended
ok('turn 155 _inheritanceCacheKey signature still extended',
   /function\s+_inheritanceCacheKey\s*\(\s*items\s*,\s*mode\s*,\s*threshold\s*\)/.test(html));

// ============================================================================
// 5. Verify the bridge file exists on disk (the actual fix this enables)
// ============================================================================
console.log('\n=== 5. Bridge file exists on disk ===');

const bridgeFile = path.resolve(__dirname, '..', 'js', 'atlas_dosage_bridge.js');
ok('js/atlas_dosage_bridge.js exists',
   fs.existsSync(bridgeFile));

if (fs.existsSync(bridgeFile)) {
  const bridge = fs.readFileSync(bridgeFile, 'utf8');
  ok('bridge exports AtlasDosageBridge on root',
     /root\.AtlasDosageBridge\s*=\s*factory\(\)/.test(bridge));
  ok('bridge install() is async',
     /async\s+function\s+install\s*\(\s*opts\s*\)/.test(bridge));
  ok('bridge install accepts atlasServer + stateRef opts',
     /opts\.stateRef/.test(bridge) && /opts\.atlasServer/.test(bridge));
}

// ============================================================================
// SUMMARY
// ============================================================================
console.log('\n=== SUMMARY ===');
console.log(`PASS: ${pass}`);
console.log(`FAIL: ${fail}`);
process.exit(fail > 0 ? 1 : 0);
