// run_consensus_precomp.mjs — driver for the LANTA precomp JSON shape
// =====================================================================
// Reads the precomp JSON produced by the discovery pipeline (chrom,
// windows[].pc1[], l2_envelopes[]) and runs consensus_partition on a
// chosen locus.
//
// The locus is defined by one or more L2 envelopes. Selection options:
//   --envelope <id>          one specific L2 envelope by candidate_id
//   --envelopes <id,id,...>  comma-separated list of envelopes
//   --bp-range <start>-<end>  all L2 envelopes overlapping a bp interval
//   (default: every L2 envelope in the precomp, run as separate loci)
//
// Each (envelope × cluster) is a "band". Within an envelope, sample
// labels come from running 1-D K-means on the per-sample mean PC1
// across that envelope's windows. Across envelopes, source bands
// "vote" about target bands by Jaccard overlap on shared samples.
//
// USAGE:
//   node run_consensus_precomp.mjs --in <precomp.json> [options]
//
// OPTIONS:
//   --in <path>             precomp JSON (required)
//   --envelope <id>         single envelope by candidate_id
//   --envelopes <id,...>    multiple envelopes (forming one locus)
//   --bp-range <s>-<e>      bp interval; selects overlapping envelopes
//   --k <int>               clusters per envelope (default 3)
//   --out-dir <path>        write per-locus JSON reports here (default ./)
//   --k-cap <int>           bruteforce cap (default 10; locus ≥ 11 bands → skipped)
//   --partner-thr <f>       derive_partner_sets primaryThr (default 0.7)
//   --ambiguous-thr <f>     ambiguousBandThr (default 0.5)
//   --delta <f>             top-N adaptive δ (default 0.10)
//   --min-samples-per-cluster <int>  drop small clusters (default 5; matches minNGroup)
//
// OUTPUT:
//   <out-dir>/<locus_id>_consensus.json   — full per-locus report
//   stdout                                 — human summary table
// =====================================================================

import fs from 'node:fs';
import path from 'node:path';
import {
  consensus_partition,
  PATTERN_CLASS,
} from './Atlas/shared/band_tracking/index_min.js';
import { kmeans1D } from './Atlas/shared/band_tracking/_kmeans_imported.js';

// ---------- argv ----------
function parseArgs(argv) {
  const out = {
    in: null, outDir: '.',
    envelope: null, envelopes: null, bpRange: null,
    k: 3, kCap: 10, minSamplesPerCluster: 5,
    partnerThr: null, ambiguousThr: null, delta: null,
  };
  for (let i = 2; i < argv.length; i++) {
    const a = argv[i]; const next = () => argv[++i];
    if (a === '--in') out.in = next();
    else if (a === '--out-dir') out.outDir = next();
    else if (a === '--envelope') out.envelope = next();
    else if (a === '--envelopes') out.envelopes = next().split(',');
    else if (a === '--bp-range') {
      const m = next().match(/^(\d+)-(\d+)$/);
      if (!m) throw new Error('--bp-range expects format START-END');
      out.bpRange = [+m[1], +m[2]];
    }
    else if (a === '--k') out.k = +next();
    else if (a === '--k-cap') out.kCap = +next();
    else if (a === '--min-samples-per-cluster') out.minSamplesPerCluster = +next();
    else if (a === '--partner-thr') out.partnerThr = +next();
    else if (a === '--ambiguous-thr') out.ambiguousThr = +next();
    else if (a === '--delta') out.delta = +next();
    else if (a === '-h' || a === '--help') { printHelp(); process.exit(0); }
    else { console.error(`Unknown arg: ${a}`); printHelp(); process.exit(2); }
  }
  if (!out.in) { printHelp(); process.exit(2); }
  return out;
}

function printHelp() {
  console.error(
`run_consensus_precomp.mjs — feed precomp JSON into consensus_partition.

Required:
  --in <path>             precomp JSON (e.g. LG28.json)

Locus selection (pick one; default = all envelopes, each as own locus):
  --envelope <id>         single envelope by candidate_id (e.g. C_gar_LG28_d17L2_0010_05)
  --envelopes <id,id,...> multiple envelopes forming one locus
  --bp-range <s>-<e>      bp interval; locus = all overlapping envelopes

Tunables:
  --k <int>               clusters per envelope (default 3)
  --out-dir <path>        output dir (default ./)
  --k-cap <int>           bruteforce cap (default 10)
  --min-samples-per-cluster <int>  default 5
  --partner-thr <f>       default 0.7   --ambiguous-thr <f>  default 0.5
  --delta <f>             default 0.10
`);
}

// ---------- aggregation: per-sample mean PC1 across an envelope's windows ----------
//
// Reproduces the relevant slice of clusterL2 + aggregateL2 from
// Atlas/shared/per_l2_cluster.js, without dragging in the full state
// object that module expects. Aggregation mode = mean_pc1 (the default
// in the production atlas).
function aggregateEnvelope(precomp, env) {
  const N = precomp.n_samples;
  const sums = new Float64Array(N);
  const counts = new Int32Array(N);
  // Envelope is defined by start_w/end_w (inclusive window indices into
  // precomp.windows[]). The legacy code uses ._s0/._e0 alias fields,
  // but the static precomp uses start_w/end_w directly.
  const s0 = env.start_w;
  const e0 = env.end_w;
  for (let w = s0; w <= e0; w++) {
    const win = precomp.windows[w];
    if (!win || !win.pc1) continue;
    for (let s = 0; s < N; s++) {
      const v = win.pc1[s];
      if (v == null || Number.isNaN(v)) continue;
      sums[s] += v;
      counts[s]++;
    }
  }
  const xs = new Float64Array(N);
  for (let s = 0; s < N; s++) {
    xs[s] = counts[s] > 0 ? sums[s] / counts[s] : NaN;
  }
  return { xs, n_windows: e0 - s0 + 1 };
}

// ---------- cluster one envelope, returning per-sample labels ----------
function clusterEnvelope(precomp, env, k, minSamplesPerCluster) {
  const agg = aggregateEnvelope(precomp, env);
  // Drop NaNs by carrying a parallel index map
  const validIdx = [];
  const validXs = [];
  for (let s = 0; s < agg.xs.length; s++) {
    if (Number.isFinite(agg.xs[s])) {
      validIdx.push(s);
      validXs.push(agg.xs[s]);
    }
  }
  const n = validXs.length;
  if (n < k) return null;
  const r = kmeans1D(Float64Array.from(validXs), k);
  // Map labels back to original sample-index space
  const labelsBySample = new Int32Array(agg.xs.length).fill(-1);
  for (let i = 0; i < validIdx.length; i++) {
    labelsBySample[validIdx[i]] = r.labels[i];
  }
  // Reject clusters smaller than minSamplesPerCluster
  for (let kk = 0; kk < k; kk++) {
    if (r.n_per_group[kk] < minSamplesPerCluster) {
      return { labels: labelsBySample, centers: Array.from(r.centers),
               n_per_group: Array.from(r.n_per_group), n_windows: agg.n_windows,
               ok: false, reason: `cluster ${kk} too small (n=${r.n_per_group[kk]})` };
    }
  }
  return { labels: labelsBySample, centers: Array.from(r.centers),
           n_per_group: Array.from(r.n_per_group), n_windows: agg.n_windows,
           ok: true };
}

// ---------- locus → voteRecords ----------
//
// Each (envelope_idx_within_locus, cluster_id) is a band. A vote is
// emitted from each source (env_i, c_i) about every other (env_j, c_j):
// Jaccard ≥ 0.5 → visited, ≤ 0.2 → excluded.
function locusToVoteRecords(precomp, envelopes, k, minSamplesPerCluster) {
  const clusters = []; // [{ env, labels, ok, n_per_group, ... }]
  for (const env of envelopes) {
    const c = clusterEnvelope(precomp, env, k, minSamplesPerCluster);
    if (c === null) {
      clusters.push({ env, labels: null, ok: false, reason: 'too few valid samples for K-means' });
    } else {
      clusters.push({ env, ...c });
    }
  }

  // Build band id space: one band per (envelope, cluster_id) for ok envelopes
  const bandIdRev = []; // [{ env_idx_in_locus, env_id, cluster, n }]
  const bandIdBy = new Map(); // "envIdx:c" → bandId
  for (let ei = 0; ei < clusters.length; ei++) {
    const c = clusters[ei];
    if (!c.ok) continue;
    for (let cc = 0; cc < k; cc++) {
      const id = bandIdRev.length;
      bandIdBy.set(`${ei}:${cc}`, id);
      bandIdRev.push({
        env_idx_in_locus: ei,
        env_id: c.env.candidate_id,
        cluster: cc,
        n: c.n_per_group[cc],
      });
    }
  }

  // Per-band sample sets
  const bandSamples = bandIdRev.map(b => {
    const c = clusters[b.env_idx_in_locus];
    const set = new Set();
    for (let s = 0; s < c.labels.length; s++) {
      if (c.labels[s] === b.cluster) set.add(s);
    }
    return set;
  });

  const voteRecords = [];
  const J_VISITED = 0.5;
  const J_EXCLUDED = 0.2;
  for (let src = 0; src < bandIdRev.length; src++) {
    const srcSet = bandSamples[src];
    const srcEi = bandIdRev[src].env_idx_in_locus;
    const visited = [];
    const excluded = [];
    for (let tgt = 0; tgt < bandIdRev.length; tgt++) {
      if (tgt === src) continue;
      if (bandIdRev[tgt].env_idx_in_locus === srcEi) continue;  // same envelope
      const tgtSet = bandSamples[tgt];
      let inter = 0;
      for (const s of srcSet) if (tgtSet.has(s)) inter++;
      const uni = srcSet.size + tgtSet.size - inter;
      const J = uni > 0 ? inter / uni : 0;
      if (J >= J_VISITED) visited.push(tgt);
      else if (J <= J_EXCLUDED) excluded.push(tgt);
    }
    let pattern_class;
    const totalBands = bandIdRev.length;
    const mentioned = visited.length + excluded.length;
    if (visited.length === 0 && excluded.length === 0) pattern_class = PATTERN_CLASS.EMPTY;
    else if (visited.length === 1 && excluded.length >= 2) pattern_class = PATTERN_CLASS.SINGLE;
    else if (visited.length >= 2 && excluded.length >= 2) pattern_class = PATTERN_CLASS.SUBSET;
    else if (visited.length >= 2 && excluded.length === 0) pattern_class = PATTERN_CLASS.SUBSET_SPLIT;
    else if (mentioned > 0 && mentioned < totalBands * 0.4) pattern_class = PATTERN_CLASS.SPLIT_TWO;
    else pattern_class = PATTERN_CLASS.SCATTER;
    voteRecords.push({
      source_w: srcEi,
      source_k: src,
      visited_bands: visited,
      excluded_bands: excluded,
      pattern_class,
      n_source: srcSet.size,
    });
  }
  return { voteRecords, bandIdRev, clusters };
}

// ---------- locus selection ----------
function selectEnvelopes(precomp, args) {
  const all = precomp.l2_envelopes;
  if (args.envelope) {
    const e = all.find(x => x.candidate_id === args.envelope);
    if (!e) throw new Error(`Envelope not found: ${args.envelope}`);
    return [{ locus_id: args.envelope, envelopes: [e] }];
  }
  if (args.envelopes) {
    const set = new Set(args.envelopes);
    const sel = all.filter(x => set.has(x.candidate_id));
    if (sel.length !== args.envelopes.length) {
      const missing = args.envelopes.filter(id => !sel.find(e => e.candidate_id === id));
      throw new Error(`Envelopes not found: ${missing.join(', ')}`);
    }
    return [{ locus_id: args.envelopes.join('+'), envelopes: sel }];
  }
  if (args.bpRange) {
    const [s, e] = args.bpRange;
    const sel = all.filter(x => x.end_bp >= s && x.start_bp <= e);
    if (sel.length === 0) throw new Error(`No envelopes overlap ${s}-${e} bp`);
    const id = `${precomp.chrom}_${(s/1e6).toFixed(3)}-${(e/1e6).toFixed(3)}Mb`;
    return [{ locus_id: id, envelopes: sel }];
  }
  // Default: every envelope as its own locus
  return all.map(env => ({ locus_id: env.candidate_id, envelopes: [env] }));
}

// ---------- main ----------
const args = parseArgs(process.argv);
const precomp = JSON.parse(fs.readFileSync(args.in, 'utf8'));
console.log(`# precomp: ${args.in}  chrom=${precomp.chrom}  ` +
            `samples=${precomp.n_samples}  windows=${precomp.n_windows}  ` +
            `L2_envelopes=${precomp.l2_envelopes.length}`);

if (!fs.existsSync(args.outDir)) fs.mkdirSync(args.outDir, { recursive: true });

const loci = selectEnvelopes(precomp, args);
console.log(`# loci to process: ${loci.length}`);

function toJSONSafe(o, depth = 0) {
  if (depth > 8) return '[truncated]';
  if (o == null || typeof o !== 'object') return o;
  if (o instanceof Float64Array || o instanceof Float32Array || o instanceof Int32Array || o instanceof Int8Array) {
    return Array.from(o);
  }
  if (Array.isArray(o)) return o.map(x => toJSONSafe(x, depth + 1));
  if (o instanceof Map) {
    const m = {};
    for (const [k, v] of o.entries()) m[String(k)] = toJSONSafe(v, depth + 1);
    return m;
  }
  if (o instanceof Set) return Array.from(o);
  const out = {};
  for (const k of Object.keys(o)) out[k] = toJSONSafe(o[k], depth + 1);
  return out;
}

// Summary table header
console.log('');
console.log('locus_id'.padEnd(40) + '  ' +
            'n_env'.padStart(5) + '  ' +
            'K'.padStart(3) + '  ' +
            'class'.padEnd(22) + '  ' +
            'score'.padStart(7) + '  ' +
            'residual'.padStart(8) + '  ' +
            'top_blocks');
console.log('-'.repeat(140));

const opts = {};
if (args.partnerThr   != null) opts.partnerSetOpts = { primaryThr: args.partnerThr };
if (args.ambiguousThr != null) opts.classifyOpts   = { ambiguousBandThr: args.ambiguousThr };
if (args.kCap         != null) opts.enumOpts       = { K_cap: args.kCap };
if (args.delta        != null) opts.selectOpts     = { delta: args.delta };

const rows = [];
for (const { locus_id, envelopes } of loci) {
  const r = locusToVoteRecords(precomp, envelopes, args.k, args.minSamplesPerCluster);
  if (r.bandIdRev.length === 0) {
    console.log(`${locus_id}`.padEnd(40) + '  ' +
                String(envelopes.length).padStart(5) + '  ' +
                String(args.k).padStart(3) + '  ' +
                'NO_OK_ENVELOPES'.padEnd(22) + '  (all clusters too small)');
    continue;
  }
  // band_groups: bands from the same envelope can never be co-grouped
  // (they are alternative arrangements at one locus position).
  // bandIdRev[b].env_idx_in_locus is exactly the group id we want.
  const band_groups = r.bandIdRev.map(b => b.env_idx_in_locus);
  const result = consensus_partition(r.voteRecords, { ...opts, band_groups });
  const top = result.top_partitions[0];
  const topBlocks = top ? JSON.stringify(top.blocks) : '—';
  const score = top && Number.isFinite(top.score) ? top.score.toFixed(4) : 'NaN';
  console.log(`${locus_id}`.padEnd(40) + '  ' +
              String(envelopes.length).padStart(5) + '  ' +
              String(args.k).padStart(3) + '  ' +
              result.consensus_class.padEnd(22) + '  ' +
              score.padStart(7) + '  ' +
              result.pca_hidden_regime_residual.toFixed(4).padStart(8) + '  ' +
              topBlocks.slice(0, 80));
  rows.push({ locus_id, envelopes: envelopes.map(e => e.candidate_id),
              consensus_class: result.consensus_class,
              score: top ? top.score : null,
              residual: result.pca_hidden_regime_residual,
              n_bands: r.bandIdRev.length,
              top_blocks: top ? top.blocks : null,
              ambiguous_band_ids: result.ambiguous_band_ids });
  // Per-locus full report
  const outPath = path.join(args.outDir, `${locus_id}_consensus.json`);
  fs.writeFileSync(outPath, JSON.stringify({
    generated_at: new Date().toISOString(),
    locus_id, envelopes: envelopes.map(e => e.candidate_id),
    args,
    bandIdRev: r.bandIdRev,
    cluster_diagnostics: r.clusters.map(c => ({
      env_id: c.env.candidate_id,
      ok: c.ok, reason: c.reason || null,
      n_per_group: c.n_per_group || null,
      centers: c.centers || null,
      n_windows: c.n_windows || null,
    })),
    n_voteRecords: r.voteRecords.length,
    result: toJSONSafe(result),
  }, null, 2));
}

// Aggregate summary
const summaryPath = path.join(args.outDir, '_summary.json');
fs.writeFileSync(summaryPath, JSON.stringify({
  precomp: args.in,
  chrom: precomp.chrom,
  n_samples: precomp.n_samples,
  k: args.k,
  thresholds: {
    partnerThr: args.partnerThr ?? '(module default 0.7)',
    ambiguousThr: args.ambiguousThr ?? '(module default 0.5)',
    delta: args.delta ?? '(module default 0.10)',
    kCap: args.kCap,
  },
  rows,
}, null, 2));

console.log('');
console.log(`# wrote ${rows.length} per-locus reports to ${args.outDir}/`);
console.log(`# summary: ${summaryPath}`);

// Class histogram
const hist = {};
for (const row of rows) hist[row.consensus_class] = (hist[row.consensus_class] || 0) + 1;
console.log('');
console.log('# class histogram:');
for (const [cls, n] of Object.entries(hist).sort((a, b) => b[1] - a[1])) {
  console.log(`#   ${cls.padEnd(25)} ${n}`);
}
