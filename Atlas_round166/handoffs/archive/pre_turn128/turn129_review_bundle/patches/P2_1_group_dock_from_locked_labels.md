# P2.1 — Group dock: build candidate dimensions from `locked_labels`

## Risk: medium
## Lines changed: ~80
## Depends on: P1.1 (popgen wiring needs collectPopstatsTracks intact)
## Verification: G key opens dock, candidate dims appear, Compute returns FST/dXY/theta from server

---

## What's wrong

You observed:
```
state.candidate exists
state.candidate.locked_labels.length = 226    // K=3 labels per sample
state.candidate.fish_calls = undefined        // not built
state.data.samples.length = 226
```

The group dock only shows global / lasso dims. Candidate dims like
"diploid_class@<id>" with H1/H1, H1/H2, H2/H2 chips don't appear.

You verified manually that conversion works:
- 0 → H1/H1 (61 samples)
- 1 → H1/H2 (103)
- 2 → H2/H2 (62)

So the data is correct. The group engine just doesn't synthesize
`fish_calls` from `locked_labels` when the former is missing.

## Where to fix

Two reasonable architectures. Pick one.

### Architecture A — group engine reads `locked_labels` directly (recommended)

The group engine builds candidate dimensions on demand. When it
encounters a candidate with `locked_labels` but no `fish_calls`, it
builds calls in-memory without mutating the candidate. This means:
- No risk of stale `fish_calls` getting persisted to session JSON.
- Group engine remains the single source of truth for "what dims
  exist for this candidate".

### Architecture B — synthesize `fish_calls` once, attach to candidate

When `locked_labels` is set on a candidate, also populate
`fish_calls`. The whole codebase keeps reading `fish_calls` and
nothing changes downstream.

Risk: anywhere else that builds `fish_calls` from real K=3
voting (the proper compute path) needs to win when both can
fire. Otherwise we end up overwriting good data with the
fallback shim.

**Pick Architecture A.** It's the smaller change and doesn't
introduce a write to the candidate object.

## The patch

Anchor strings (production atlas):
- `function _atlasGroupEngine_collectCandidateDims` or similar
- Or look for the function that returns `dimensions` from the
  group engine when `state.candidate` is set.
- Or grep for a string like `'per candidate ·'` or
  `'diploid_class@'` in the group dock render code.

If the helper doesn't exist yet, add it adjacent to the existing
candidate dim collector.

### New helper: `_atlasCandidateBuildSyntheticFishCalls`

Place this near other group-engine helpers. Pure function. Doesn't
mutate the candidate.

```js
// turn 129 P2.1: synthesize per-sample fish_calls from a candidate's
// locked_labels when fish_calls hasn't been built yet. This unblocks
// the group dock for candidates produced from the auto-pipeline (which
// stamps locked_labels but doesn't always run the full
// computeCandidateAssignments compute path that produces fish_calls).
//
// Conventions:
//   - Input: locked_labels: Int8Array | Array<int> length n_samples
//     where each entry is in [0, K-1] or -1 for ambiguous.
//   - Output: fish_calls: Array<{sample_idx, sample_id, regime,
//     confidence, n_supporting, n_intervals, ambiguous, votes}>
//     matching the existing fish_calls shape (see comment around
//     line 25971 in atlas v62k).
//
// K=3 mapping (the primary regime in this atlas):
//   0 -> H1/H1
//   1 -> H1/H2
//   2 -> H2/H2
//
// For other K values (e.g. K=6 substructure), labels stay generic
// 'g0' .. 'gK-1'. The dock label-formatting helper handles both.
function _atlasCandidateBuildSyntheticFishCalls(cand, samples) {
  if (!cand || !cand.locked_labels) return null;
  const labels = cand.locked_labels;
  const n = labels.length || 0;
  if (n === 0) return null;
  if (!samples || samples.length !== n) {
    // Defensive: sample list and label vector must align by index.
    console.warn('[group-engine] locked_labels length', n,
                 'does not match samples length',
                 samples ? samples.length : 'null');
    return null;
  }
  const out = new Array(n);
  for (let si = 0; si < n; si++) {
    const r = labels[si];
    out[si] = {
      sample_idx: si,
      sample_id:  samples[si] && (samples[si].cga || samples[si].ind || String(si)),
      regime:     (r == null || r < 0) ? -1 : Number(r),
      confidence: (r == null || r < 0) ? 0.0 : 1.0,
      // Synthetic — we don't have per-L2 votes when reconstructing
      // from locked_labels. Mark as 1/1 to indicate full agreement
      // with the locked label.
      n_supporting: 1,
      n_intervals:  1,
      ambiguous:    (r == null || r < 0),
      votes:        [r == null ? -1 : Number(r)],
      // Marker so downstream code can tell synthetic from real.
      _source: 'synthesized_from_locked_labels',
    };
  }
  return out;
}
```

### Helper: K-aware diploid_class label

```js
// turn 129 P2.1: derive a human-readable label for a regime index.
// K=3 -> H1/H1 H1/H2 H2/H2 (the canonical inversion karyotype).
// K=6 -> g0..g5 (substructure).
// Other K -> g0..gK-1.
function _atlasRegimeLabel(regimeIdx, K) {
  if (regimeIdx == null || regimeIdx < 0) return 'ambig';
  K = K || 3;
  if (K === 3) {
    if (regimeIdx === 0) return 'H1/H1';
    if (regimeIdx === 1) return 'H1/H2';
    if (regimeIdx === 2) return 'H2/H2';
  }
  return 'g' + regimeIdx;
}
```

### Hook the group engine

Anchor: find where the group engine builds the list of candidate
dimensions. Around the line that does
`if (cand.fish_calls) { ... }` (you saw this pattern in
`state.candidate.fish_calls = undefined` causing nothing to show).

Replacement pattern:

```js
// before
const fish_calls = cand && Array.isArray(cand.fish_calls)
                 ? cand.fish_calls : null;
if (!fish_calls) {
  // shows nothing
  return globalDimsOnly;
}

// after
let fish_calls = cand && Array.isArray(cand.fish_calls)
               ? cand.fish_calls : null;
if (!fish_calls && cand && cand.locked_labels && state.data && state.data.samples) {
  // turn 129 P2.1: candidates from auto-pipeline often have only
  // locked_labels, not fish_calls. Synthesize on-demand.
  fish_calls = _atlasCandidateBuildSyntheticFishCalls(cand, state.data.samples);
}
if (!fish_calls) {
  return globalDimsOnly;
}
```

Then the dim builder uses `fish_calls` as before. The K used in
labelling can come from `cand.K` (typically 3) or fall back to 3.

### Build the dimension entries

The group engine should produce one dim per (candidate, regime
class):

```js
function _atlasBuildCandidateDims(cand, fish_calls) {
  const K = cand.K || 3;
  const candId = cand.id || 'candidate';
  // Bucket samples by regime
  const buckets = new Map();   // regimeIdx -> Array<sample_idx>
  for (const fc of fish_calls) {
    if (!fc || fc.regime == null || fc.regime < 0) continue;
    if (!buckets.has(fc.regime)) buckets.set(fc.regime, []);
    buckets.get(fc.regime).push(fc.sample_idx);
  }
  const dims = [];
  for (const [regime, sampleIdxList] of buckets) {
    const label = _atlasRegimeLabel(regime, K);
    dims.push({
      // Group dock display
      group:     'per candidate · ' + candId,
      dim_id:    'diploid_class@' + candId,
      chip:      label,
      n_samples: sampleIdxList.length,
      // Compute payload — sample IDs the server expects
      sample_indices: sampleIdxList,
      sample_ids:     sampleIdxList.map(si =>
        state.data.samples[si] &&
        (state.data.samples[si].cga ||
         state.data.samples[si].ind ||
         String(si))).filter(Boolean),
      // Marker
      regime_idx: regime,
      candidate_id: candId,
    });
  }
  // Sort dims by regime_idx so chips render H1/H1 H1/H2 H2/H2 in order
  dims.sort((a, b) => a.regime_idx - b.regime_idx);
  return dims;
}
```

### Wire the Compute button to send `sample_ids`

Anchor: find the function that posts to `/api/popstats/groupwise`.
It should accept the dim payload and send explicit sample IDs:

```js
// Compute payload
{
  chrom: state.data.chrom,
  start_bp: cand.start_bp,
  end_bp:   cand.end_bp,
  groups: {
    HOM1: dim_for_H1H1.sample_ids,
    HET:  dim_for_H1H2.sample_ids,
    HOM2: dim_for_H2H2.sample_ids,
  },
  metrics: ['fst', 'dxy', 'theta_pi'],
}
```

Slot mapping (HOM1/HET/HOM2) comes from the dock's slot UI. You
drag a chip onto a slot, the slot remembers which dim is bound,
and Compute reads from that binding.

## Verification

1. Apply the patch.
2. Reload atlas.
3. Make a candidate active (one with `locked_labels`, e.g. on LG28).
4. Press G to open the group dock.
5. Expect: `per candidate · <id>` section with chips `H1/H1 (61)`,
   `H1/H2 (103)`, `H2/H2 (62)`.
6. Drag chips onto HOM1 / HET / HOM2 slots.
7. Click Compute.
8. Server returns FST, dXY, theta_pi rows; charts render.

## Test (source-level)

```js
// tests/test_p2_1_synthetic_fish_calls.js
const fs = require('fs');
const vm = require('vm');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

ok('helper _atlasCandidateBuildSyntheticFishCalls defined',
   /function _atlasCandidateBuildSyntheticFishCalls\(/.test(html));
ok('helper _atlasRegimeLabel defined',
   /function _atlasRegimeLabel\(/.test(html));
ok('label H1/H1 mapped from regime 0',
   /regimeIdx === 0[^\n]*'H1\/H1'/.test(html));
ok('label H1/H2 mapped from regime 1',
   /regimeIdx === 1[^\n]*'H1\/H2'/.test(html));
ok('label H2/H2 mapped from regime 2',
   /regimeIdx === 2[^\n]*'H2\/H2'/.test(html));
ok('K-fallback uses g<n>',
   /'g'\s*\+\s*regimeIdx/.test(html));
ok('synthetic marker present',
   /_source:\s*'synthesized_from_locked_labels'/.test(html));

// Behavioural: extract helpers and test conversion
function pullFn(src, name) {
  const re = new RegExp('^function\\s+' + name + '\\s*\\(', 'm');
  const m = src.match(re); if (!m) return null;
  let i = src.indexOf('{', m.index), depth = 1; i++;
  while (i < src.length && depth > 0) {
    const c = src[i];
    if (c === '{') depth++;
    else if (c === '}') depth--;
    i++;
  }
  return src.substring(m.index, i);
}
const fnSyn   = pullFn(html, '_atlasCandidateBuildSyntheticFishCalls');
const fnLabel = pullFn(html, '_atlasRegimeLabel');
const sandbox = `${fnSyn}\n${fnLabel}\nmodule.exports = { _atlasCandidateBuildSyntheticFishCalls, _atlasRegimeLabel };`;
const ctx = { module: { exports: {} }, console };
vm.createContext(ctx);
vm.runInContext(sandbox, ctx);
const m = ctx.module.exports;

const labels = new Int8Array([0, 0, 1, 1, 1, 2, -1]);
const samples = [{cga:'a'},{cga:'b'},{cga:'c'},{cga:'d'},{cga:'e'},{cga:'f'},{cga:'g'}];
const fc = m._atlasCandidateBuildSyntheticFishCalls({ locked_labels: labels }, samples);
ok('synth length matches input', fc && fc.length === 7);
ok('regime 0 sample has regime=0', fc[0].regime === 0);
ok('regime 1 sample has regime=1', fc[2].regime === 1);
ok('regime 2 sample has regime=2', fc[5].regime === 2);
ok('ambiguous sample marked',     fc[6].ambiguous === true);
ok('sample_id from cga',          fc[0].sample_id === 'a');
ok('K=3 label H1/H1',             m._atlasRegimeLabel(0, 3) === 'H1/H1');
ok('K=3 label H1/H2',             m._atlasRegimeLabel(1, 3) === 'H1/H2');
ok('K=3 label H2/H2',             m._atlasRegimeLabel(2, 3) === 'H2/H2');
ok('K=6 label g4',                m._atlasRegimeLabel(4, 6) === 'g4');
ok('ambig label',                 m._atlasRegimeLabel(-1, 3) === 'ambig');

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## Risk notes

- If the candidate already has `fish_calls` (real, computed from
  `computeCandidateAssignments`), the patch is a no-op — the
  fallback only fires when `fish_calls` is null/undefined.
- The synthesized calls have `n_supporting: 1, n_intervals: 1,
  votes: [r]` — a downstream caller that queries
  `fc.n_supporting / fc.n_intervals` for "agreement strength" will
  see 1.0 (full agreement with the locked label). This is honest:
  we don't have the per-L2 vote breakdown when reconstructing.
- The `_source: 'synthesized_from_locked_labels'` marker lets future
  diagnostics distinguish synthesized from real calls.
