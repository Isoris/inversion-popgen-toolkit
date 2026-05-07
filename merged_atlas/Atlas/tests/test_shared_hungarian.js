// tests/test_shared_hungarian.js

import {
  alignLabels, permutations,
  hungarianChainProjection, concordanceMatrix,
  LINEAGE_CHAIN_BREAK_AGREEMENT,
} from '../shared/hungarian.js';

let pass = 0, fail = 0;
function check(name, cond, detail = '') {
  if (cond) { console.log(`  ✓ ${name}${detail ? '  (' + detail + ')' : ''}`); pass++; }
  else      { console.log(`  ✗ ${name}  ${detail}`); fail++; }
}
function approx(a, b, tol = 1e-9) { return Math.abs(a - b) <= tol; }

console.log('--- LINEAGE_CHAIN_BREAK_AGREEMENT ---');
check('exported as 0.50',  LINEAGE_CHAIN_BREAK_AGREEMENT === 0.50);

console.log('\n--- permutations ---');
{
  const p1 = permutations(1);
  check('K=1 → 1 perm',        p1.length === 1 && p1[0][0] === 0);
  const p2 = permutations(2);
  check('K=2 → 2 perms',       p2.length === 2);
  const p3 = permutations(3);
  check('K=3 → 6 perms',       p3.length === 6);
  // Check uniqueness
  const seen = new Set();
  for (const p of p3) seen.add(p.join(','));
  check('K=3 perms unique',    seen.size === 6);
  const p4 = permutations(4);
  check('K=4 → 24 perms',      p4.length === 24);
  const p6 = permutations(6);
  check('K=6 → 720 perms',     p6.length === 720);
}

console.log('\n--- alignLabels (perfect identity) ---');
{
  const g1 = [0, 0, 1, 1, 2, 2];
  const g2 = [0, 0, 1, 1, 2, 2];
  const r = alignLabels(g1, g2, 3);
  check('perm = identity',      r.perm[0] === 0 && r.perm[1] === 1 && r.perm[2] === 2);
  check('concord = 1',           r.concord === 1);
  check('aligned matches g2',    Array.from(r.aligned).every((v, i) => v === g2[i]));
}

console.log('\n--- alignLabels (rotated labels) ---');
{
  // Same partition but with labels rotated: 0→1, 1→2, 2→0
  const g1 = [0, 0, 1, 1, 2, 2];
  const g2 = [1, 1, 2, 2, 0, 0];
  const r = alignLabels(g1, g2, 3);
  check('concord = 1 after rotation',  r.concord === 1);
  // Aligned should equal g1 since the partition is the same
  check('aligned = g1',                Array.from(r.aligned).every((v, i) => v === g1[i]));
  // perm[0]=1 (g1 label 0 aligns with g2 label 1)
  check('perm[0] = 1',                 r.perm[0] === 1);
  check('perm[1] = 2',                 r.perm[1] === 2);
  check('perm[2] = 0',                 r.perm[2] === 0);
}

console.log('\n--- alignLabels (genuinely different partitions) ---');
{
  // g1 = [0,0,1,1] (split at i=2)
  // g2 = [0,1,0,1] (alternating)
  // Best perm yields concord = 0.5
  const g1 = [0, 0, 1, 1];
  const g2 = [0, 1, 0, 1];
  const r = alignLabels(g1, g2, 2);
  check('concord = 0.5 for crossed partition',  r.concord === 0.5);
}

console.log('\n--- hungarianChainProjection (clean chain) ---');
{
  // 4 L2s, 3 samples each, all label 0=A, 1=B (cohort stays put with rotated labels)
  const labelsByL2 = {
    10: new Int8Array([0, 0, 1]),    // canonical
    20: new Int8Array([1, 1, 0]),    // rotated; aligns to 10 via perm
    30: new Int8Array([0, 0, 1]),
    40: new Int8Array([1, 1, 0]),    // rotated again
  };
  const get = (idx) => labelsByL2[idx];
  const proj = hungarianChainProjection([10, 20, 30, 40], 2, get);
  check('one chain (all stable)',    proj.n_chains === 1);
  check('n_total_L2 = 4',             proj.n_total_L2 === 4);
  check('n_samples = 3',              proj.n_samples === 3);
  check('chain.l2_indices = [10,20,30,40]',
        proj.chains[0].l2_indices.join(',') === '10,20,30,40');
  // Aligned rows should all be [0,0,1] after projection (canonical)
  const flat = proj.chains[0].projected;
  for (let row = 0; row < 4; row++) {
    const r = [flat[row*3], flat[row*3+1], flat[row*3+2]];
    check(`row ${row} aligns to canonical [0,0,1]`,
          r[0] === 0 && r[1] === 0 && r[2] === 1, `got [${r}]`);
  }
}

console.log('\n--- hungarianChainProjection (chain break) ---');
{
  // First two L2s: same K=3 partition (3+3+3)
  // Third L2: collapsed to a single label — best Hungarian align gives
  //   concord = 3/9 = 0.333, well below LINEAGE_CHAIN_BREAK_AGREEMENT (0.5).
  const labelsByL2 = {
    1: new Int8Array([0, 0, 0, 1, 1, 1, 2, 2, 2]),
    2: new Int8Array([0, 0, 0, 1, 1, 1, 2, 2, 2]),
    3: new Int8Array([0, 0, 0, 0, 0, 0, 0, 0, 0]),  // K-collapse → forces break
  };
  const get = (idx) => labelsByL2[idx];
  const proj = hungarianChainProjection([1, 2, 3], 3, get);
  check('two chains after break',          proj.n_chains === 2);
  check('chain 0 has 2 L2s',                proj.chains[0] && proj.chains[0].l2_indices.length === 2);
  check('chain 1 has 1 L2',                 proj.chains[1] && proj.chains[1].l2_indices.length === 1);
  check('chain 1 starts at idx 3',          proj.chains[1] && proj.chains[1].l2_indices[0] === 3);
}

console.log('\n--- hungarianChainProjection (edge cases) ---');
{
  const get = () => new Int8Array([0, 0, 1]);
  check('throws on missing callback',
        (() => { try { hungarianChainProjection([1], 2); return false; } catch { return true; } })());
  check('empty l2_indices',                 hungarianChainProjection([], 2, get).n_chains === 0);
  // All L2s return null
  const proj = hungarianChainProjection([1, 2, 3], 2, () => null);
  check('all-null callback → empty',        proj.n_chains === 0 && proj.n_samples === 0);
}

console.log('\n--- concordanceMatrix ---');
{
  // 3 samples, 2 L2s, sample 0 and 1 always agree, sample 2 always disagrees
  const proj = {
    chains: [{
      projected: new Int8Array([
        0, 0, 1,    // L2 0
        0, 0, 1,    // L2 1
      ]),
      n_L2: 2,
    }],
    n_samples: 3,
    n_total_L2: 2,
  };
  const C = concordanceMatrix(proj);
  check('diagonal = 1',                  C[0] === 1 && C[4] === 1 && C[8] === 1);
  check('C[0][1] = 1 (always agree)',    C[0*3+1] === 1);
  check('C[1][0] = 1 (symmetric)',       C[1*3+0] === 1);
  check('C[0][2] = 0 (always disagree)', C[0*3+2] === 0);
  check('C[1][2] = 0',                   C[1*3+2] === 0);
}
{
  // -1 labels excluded from valid pair count
  const proj = {
    chains: [{
      projected: new Int8Array([
         0,  0, -1,    // sample 2 no-call at L2 0
         0,  0,  0,    // all called at L2 1
        -1,  0,  0,    // sample 0 no-call at L2 2
      ]),
      n_L2: 3,
    }],
    n_samples: 3,
    n_total_L2: 3,
  };
  const C = concordanceMatrix(proj);
  // sample 0 vs 1: valid at L2 0 (0,0=match), L2 1 (0,0=match) → 2/2 = 1
  // sample 0 vs 2: valid only at L2 1 (0,0=match) → 1/1 = 1
  // sample 1 vs 2: valid at L2 1 (0,0=match), L2 2 (0,0=match) → 2/2 = 1
  check('skips no-calls',  approx(C[0*3+1], 1) && approx(C[0*3+2], 1) && approx(C[1*3+2], 1));
}

console.log('\n--- chain projection + concordance integration ---');
{
  // Real-flow check: chain → concordance
  const labelsByL2 = {
    1: new Int8Array([0, 0, 1, 1]),
    2: new Int8Array([1, 1, 0, 0]),  // rotated; should align
    3: new Int8Array([0, 0, 1, 1]),
  };
  const proj = hungarianChainProjection([1, 2, 3], 2, idx => labelsByL2[idx]);
  const C = concordanceMatrix(proj);
  // Samples 0 and 1 always together; 2 and 3 always together; 0 vs 2 always apart
  check('integration: 0,1 fully concordant',  C[0*4+1] === 1);
  check('integration: 2,3 fully concordant',  C[2*4+3] === 1);
  check('integration: 0,2 fully discordant',  C[0*4+2] === 0);
}

console.log('\n=================');
console.log(`pass: ${pass}   fail: ${fail}`);
console.log('=================');
process.exit(fail === 0 ? 0 : 1);
