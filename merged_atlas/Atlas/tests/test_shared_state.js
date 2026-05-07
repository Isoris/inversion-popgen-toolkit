// tests/test_shared_state.js

import {
  SLOT_REGISTRY, CROSS_ATLAS_SLOTS, PERSISTED_SLOTS,
  makeState, serializeState, deserializeState, mergeStateFromSession,
  readPersistedSlots, writePersistedSlot,
} from '../shared/state.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}

console.log('--- registry shape ---');
check('SLOT_REGISTRY is frozen',     Object.isFrozen(SLOT_REGISTRY));
check('has > 50 slots',              Object.keys(SLOT_REGISTRY).length > 50);
check('CROSS_ATLAS_SLOTS array',     Array.isArray(CROSS_ATLAS_SLOTS));
check('PERSISTED_SLOTS array',       Array.isArray(PERSISTED_SLOTS));
// Spot-check critical slots
const need = ['k', 'kMode', 'candidate', 'candidateList', 'lockedLabels', 'manualGroups',
              'tracked', 'data', 'cur', 'l3Layout', 'linesColorMode',
              'bandTraceFishSet', 'gPanelTab', 'l2SweepEnabled'];
for (const n of need) check(`SLOT_REGISTRY.${n}`, n in SLOT_REGISTRY);

console.log('\n--- tag distribution sanity ---');
const tags = {};
for (const [name, def] of Object.entries(SLOT_REGISTRY)) {
  tags[def.tag] = (tags[def.tag] || 0) + 1;
}
console.log('  tag distribution:', tags);
check('has cross_atlas slots',  (tags.cross_atlas || 0) > 5);
check('has persisted slots',    (tags.persisted   || 0) > 10);
check('has transient slots',    (tags.transient   || 0) > 10);

console.log('\n--- makeState defaults ---');
{
  const s = makeState();
  check('k defaults to 3',                  s.k === 3);
  check('kMode defaults to "fixed"',        s.kMode === 'fixed');
  check('candidate defaults to null',       s.candidate === null);
  check('candidateList defaults to []',     Array.isArray(s.candidateList) && s.candidateList.length === 0);
  check('tracked defaults to []',           Array.isArray(s.tracked));
  check('linesColorMode = "kmeans"',         s.linesColorMode === 'kmeans');
  check('viewControls is fresh object',      s.viewControls && Array.isArray(s.viewControls.pcaXY));
  check('viewControls.pcaXY = [pc1, pc2]',   s.viewControls.pcaXY[0] === 'pc1' && s.viewControls.pcaXY[1] === 'pc2');
}

console.log('\n--- makeState defaults are deep-cloned ---');
{
  const s1 = makeState();
  const s2 = makeState();
  s1.candidateList.push({ id: 'foo' });
  check('s1 push does not affect s2',       s2.candidateList.length === 0);
  s1.viewControls.pcaXY[0] = 'pc3';
  check('s1 viewControls mutation isolated', s2.viewControls.pcaXY[0] === 'pc1');
  s1.kRange[0] = 999;
  check('s1 kRange mutation isolated',       s2.kRange[0] === 2);
}

console.log('\n--- makeState curated subset ---');
{
  const s = makeState({ slots: ['k', 'candidate', 'candidateList'] });
  check('only requested slots',    Object.keys(s).length === 3);
  check('k present',                s.k === 3);
  check('candidate present',        s.candidate === null);
}

console.log('\n--- serialize/deserialize roundtrip (primitives) ---');
{
  const s = makeState();
  s.candidate = { id: 'cand_LG28_15Mb', chrom: 'LG28', start_bp: 15_000_000, end_bp: 16_500_000 };
  s.candidateList = [s.candidate, { id: 'cand_LG28_22Mb', chrom: 'LG28', start_bp: 22e6, end_bp: 23e6 }];
  s.tracked = [42, 17, 99];
  const snap = serializeState(s, { workflow: 'inversion' });
  const json = JSON.stringify(snap);
  const round = deserializeState(JSON.parse(json));
  check('candidate round-trip',         round.candidate.id === 'cand_LG28_15Mb');
  check('candidate.start_bp preserved', round.candidate.start_bp === 15_000_000);
  check('candidateList length',         round.candidateList.length === 2);
  check('candidateList[1].id',          round.candidateList[1].id === 'cand_LG28_22Mb');
  check('tracked round-trip',           Array.isArray(round.tracked) && round.tracked.length === 3 && round.tracked[1] === 17);
  check('snapshot.schema_version === 1', snap.schema_version === 1);
  check('snapshot.workflow === inversion', snap.workflow === 'inversion');
  check('snapshot.created_at is ISO string',
        typeof snap.created_at === 'string' && /^\d{4}-\d{2}-\d{2}T/.test(snap.created_at));
}

console.log('\n--- serialize/deserialize Set & Map ---');
{
  const s = makeState();
  s.activeSampleSet = new Set([1, 2, 3, 5, 8, 13]);
  s.activeSampleReasons = new Map([
    [1, 'low coverage'],
    [2, 'sex mismatch'],
    [3, 'duplicate sample'],
  ]);
  const round = deserializeState(JSON.parse(JSON.stringify(serializeState(s))));
  check('Set type preserved',           round.activeSampleSet instanceof Set);
  check('Set size correct',             round.activeSampleSet.size === 6);
  check('Set has(13)',                  round.activeSampleSet.has(13));
  check('Map type preserved',           round.activeSampleReasons instanceof Map);
  check('Map size correct',             round.activeSampleReasons.size === 3);
  check('Map.get(2) preserved',         round.activeSampleReasons.get(2) === 'sex mismatch');
}

console.log('\n--- serialize/deserialize typed arrays ---');
{
  const s = makeState();
  s.lockedLabels = Int8Array.from([0, 1, 2, 0, 1, 2, -1, 0, 1, 2]);
  const round = deserializeState(JSON.parse(JSON.stringify(serializeState(s))));
  check('lockedLabels is Int8Array',   round.lockedLabels instanceof Int8Array);
  check('lockedLabels length',          round.lockedLabels.length === 10);
  check('lockedLabels[6] = -1',         round.lockedLabels[6] === -1);
  check('lockedLabels[2] = 2',          round.lockedLabels[2] === 2);
}

console.log('\n--- serialize ignores transient slots ---');
{
  const s = makeState();
  s.candidate = { id: 'A' };          // cross_atlas — should be in snapshot
  s.cur = 1234;                        // transient — should NOT be
  s.k = 6;                             // persisted — should NOT be in snapshot (lives in localStorage)
  s.__hetRateCache = new Map([[1, 0.4]]);  // transient — should NOT be
  const snap = serializeState(s);
  check('candidate in snapshot',        'candidate' in snap.slots);
  check('cur NOT in snapshot',          !('cur' in snap.slots));
  check('k NOT in snapshot',            !('k' in snap.slots));
  check('__hetRateCache NOT in snapshot', !('__hetRateCache' in snap.slots));
}

console.log('\n--- mergeStateFromSession ---');
{
  // Sub-atlas A creates a snapshot with a candidate set.
  const sA = makeState();
  sA.candidate = { id: 'X', chrom: 'LG12' };
  sA.tracked = [1, 2, 3];
  const snap = serializeState(sA);

  // Sub-atlas B starts fresh, then loads the snapshot.
  const sB = makeState();
  sB.k = 6;                            // local persisted state
  sB.cur = 99;                          // local transient
  mergeStateFromSession(sB, snap);
  check('B picked up A.candidate',      sB.candidate.id === 'X');
  check('B picked up A.tracked',        sB.tracked.length === 3 && sB.tracked[2] === 3);
  check('B kept its local k',           sB.k === 6);
  check('B kept its local cur',         sB.cur === 99);
}

console.log('\n--- writePersistedSlot / readPersistedSlots ---');
{
  // Mock localStorage
  class MockStorage {
    constructor() { this._d = new Map(); }
    getItem(k)    { return this._d.has(k) ? this._d.get(k) : null; }
    setItem(k, v) { this._d.set(k, String(v)); }
    removeItem(k) { this._d.delete(k); }
  }
  const store = new MockStorage();

  check('writePersistedSlot k=6',           writePersistedSlot('k', 6, store) === true);
  check('writePersistedSlot l3HetColoring=true', writePersistedSlot('l3HetColoring', true, store) === true);
  check('writePersistedSlot mergeThr=0.9',  writePersistedSlot('mergeThr', 0.9, store) === true);
  check('writePersistedSlot viewControls (object)',
        writePersistedSlot('viewControls', { pcaXY: ['pc3','pc4'], linesYsources: ['pc3'], linked: false }, store) === true);
  check('non-persisted slot rejected',      writePersistedSlot('candidate', { id: 'X' }, store) === false);
  check('unknown slot rejected',            writePersistedSlot('zzznope', 'val', store) === false);

  // Read back
  const persisted = readPersistedSlots(store);
  check('k coerced to int',                 persisted.k === 6);
  check('l3HetColoring coerced to bool',    persisted.l3HetColoring === true);
  check('mergeThr coerced to float',        Math.abs(persisted.mergeThr - 0.9) < 1e-12);
  check('viewControls round-trips JSON',
        persisted.viewControls && persisted.viewControls.pcaXY[0] === 'pc3' && persisted.viewControls.linked === false);
}

console.log('\n--- snapshot version handling ---');
{
  // Future-version snapshot — should warn but accept
  const futureSnap = { schema_version: 99, workflow: 'inversion', slots: { tracked: [7,8] } };
  const orig = console.warn;
  let warned = false;
  console.warn = () => { warned = true; };
  try {
    const out = deserializeState(futureSnap);
    check('future-version still decodes',   Array.isArray(out.tracked) && out.tracked[0] === 7);
    check('future-version warned',          warned);
  } finally { console.warn = orig; }
}

console.log('\n--- empty / malformed input handling ---');
{
  check('deserialize null → {}',            Object.keys(deserializeState(null)).length === 0);
  check('deserialize undefined → {}',       Object.keys(deserializeState(undefined)).length === 0);
  check('deserialize {} → {}',              Object.keys(deserializeState({})).length === 0);
  check('deserialize {slots:{}} → {}',      Object.keys(deserializeState({ schema_version: 1, slots: {} })).length === 0);
  check('serialize empty state',            (function(){
    const snap = serializeState(makeState());
    return typeof snap === 'object' && snap.schema_version === 1;
  })());
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
