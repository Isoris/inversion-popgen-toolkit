// run_consensus.mjs — generic consensus_partition driver
// =====================================================================
// Feeds any candidate's precomp JSON into consensus_partition and
// reports the result. Locus-agnostic by design: nothing about LG28,
// any chromosome, any inversion catalogue entry, any K_locus.
//
// Optional `candidate_id` field in the input JSON is echoed into the
// output filename and console header. If absent, the driver uses the
// input filename's stem.
//
// Three input modes — all end up as the same `voteRecords` array fed
// into consensus_partition.
//
// USAGE:
//   node run_consensus.mjs --mode <kmeans|projections|voterecords> --in <file>
//
// MODES:
//
// 1. --mode voterecords
//    Input is already a voteRecords array (or {voteRecords: [...]}):
//      [{source_w, source_k, visited_bands, excluded_bands,
//        pattern_class, n_source}, ...]
//    No reshaping. Use this when an upstream LANTA step emits Layer-2
//    projections in this form.
//
// 2. --mode projections
//    Input is the output of the per-source-band projection sweep:
//      { candidate_id?: string,
//        projections: [
//          { source_k, per_window: [
//              { w, ve: {visited_bands, excluded_bands},
//                cls: {pattern_class},
//                pv: {n_source, purity_vector?} },
//              ...
//          ]},
//          ...
//      ]}
//    Reshaped via voteRecords_from_projections from vote_evidence.js.
//
// 3. --mode kmeans
//    Raw per-window K-means labels:
//      { candidate_id?: string,
//        windows: [
//          { w, start_bp?, end_bp?,
//            band_labels: [int per sample],   // K-means cluster id
//            sample_ids?: [str]               // optional, for cross-window join
//          },
//          ...
//        ],
//        n_samples?: int }
//
//    Driver builds a Jaccard-based source-target projection per window
//    pair. Each (window, cluster) becomes a band id; voters at one
//    window project onto bands at OTHER windows.
//
//    *** kmeans mode is a fallback for data without an explicit
//    *** projection sweep. Read FINDING_2026-05-06_both_excluded_bug.md
//    *** for the K_arrangements > 2 caveat before relying on output.
//
// OUTPUT:
//   Writes a JSON report to ./<candidate_id>_consensus.json (or to
//   --out <path> if given). Human summary to stdout.
// =====================================================================

import fs from 'node:fs';
import path from 'node:path';
import {
  consensus_partition,
  voteRecords_from_projections,
  PATTERN_CLASS,
} from './Atlas/shared/band_tracking/index_min.js';

// ---------- argv parsing ----------
function parseArgs(argv) {
  const out = {
    mode: null, in: null, out: null,
    partnerThr: null, ambiguousThr: null, delta: null, kCap: null,
    candidateId: null,
  };
  for (let i = 2; i < argv.length; i++) {
    const a = argv[i];
    const next = () => argv[++i];
    if (a === '--mode') out.mode = next();
    else if (a === '--in') out.in = next();
    else if (a === '--out') out.out = next();
    else if (a === '--candidate-id') out.candidateId = next();
    else if (a === '--partner-thr') out.partnerThr = +next();
    else if (a === '--ambiguous-thr') out.ambiguousThr = +next();
    else if (a === '--delta') out.delta = +next();
    else if (a === '--k-cap') out.kCap = +next();
    else if (a === '-h' || a === '--help') { printHelp(); process.exit(0); }
    else { console.error(`Unknown arg: ${a}`); printHelp(); process.exit(2); }
  }
  if (!out.mode || !out.in) { printHelp(); process.exit(2); }
  return out;
}

function printHelp() {
  console.error(
`run_consensus.mjs — feed any candidate's precomp into consensus_partition.

Usage:
  node run_consensus.mjs --mode <kmeans|projections|voterecords> --in <file>

Options:
  --out <path>          output JSON path (default: <candidate_id>_consensus.json)
  --candidate-id <id>   override candidate id (otherwise read from input JSON
                        or derived from the input filename's stem)

Threshold overrides (optional — defaults inherited from the band_tracking
modules; passed through only when explicitly set so the modules' own
defaults remain the single source of truth):
  --partner-thr <f>      derive_partner_sets primaryThr (module default 0.7)
  --ambiguous-thr <f>    classify_consensus ambiguousBandThr (default 0.5)
  --delta <f>            select_top_partitions_adaptive delta (default 0.10)
  --k-cap <int>          enumerate_and_score_all_partitions K_cap (default 10)
`);
}

// ---------- mode 3: kmeans → voteRecords ----------
//
// See header comment + FINDING_2026-05-06_both_excluded_bug.md for caveats.
// Each (window, cluster) becomes a band id; sources project across
// windows via Jaccard on shared sample membership.
function kmeansToVoteRecords(j) {
  if (!j.windows || !Array.isArray(j.windows)) {
    throw new Error('kmeans input: missing top-level "windows" array');
  }
  const bandIdMap = new Map();
  const bandIdRev = [];
  function bandIdOf(w, c) {
    const k = `${w}:${c}`;
    if (bandIdMap.has(k)) return bandIdMap.get(k);
    const id = bandIdRev.length;
    bandIdMap.set(k, id);
    bandIdRev.push({ w, c });
    return id;
  }

  const setsByWindow = new Map();
  for (const win of j.windows) {
    const w = win.w;
    const m = new Map();
    setsByWindow.set(w, m);
    const labels = win.band_labels;
    if (!Array.isArray(labels)) {
      throw new Error(`kmeans input: window w=${w} has no band_labels array`);
    }
    for (let s = 0; s < labels.length; s++) {
      const c = labels[s];
      if (c == null || c < 0) continue;
      if (!m.has(c)) m.set(c, new Set());
      m.get(c).add(s);
      bandIdOf(w, c);
    }
  }

  const voteRecords = [];
  const J_VISITED = 0.5;
  const J_EXCLUDED = 0.2;
  for (const [src_w, srcMap] of setsByWindow.entries()) {
    for (const [src_c, srcSet] of srcMap.entries()) {
      const visited = [];
      const excluded = [];
      for (const [tgt_w, tgtMap] of setsByWindow.entries()) {
        if (tgt_w === src_w) continue;
        for (const [tgt_c, tgtSet] of tgtMap.entries()) {
          let inter = 0;
          for (const s of srcSet) if (tgtSet.has(s)) inter++;
          const uni = srcSet.size + tgtSet.size - inter;
          const J = uni > 0 ? inter / uni : 0;
          const tid = bandIdOf(tgt_w, tgt_c);
          if (J >= J_VISITED) visited.push(tid);
          else if (J <= J_EXCLUDED) excluded.push(tid);
        }
      }
      let pattern_class;
      const totalBands = bandIdRev.length;
      const mentioned = visited.length + excluded.length;
      if (visited.length === 0 && excluded.length === 0) {
        pattern_class = PATTERN_CLASS.EMPTY;
      } else if (visited.length === 1 && excluded.length >= 2) {
        pattern_class = PATTERN_CLASS.SINGLE;
      } else if (visited.length >= 2 && excluded.length >= 2) {
        pattern_class = PATTERN_CLASS.SUBSET;
      } else if (visited.length >= 2 && excluded.length === 0) {
        pattern_class = PATTERN_CLASS.SUBSET_SPLIT;
      } else if (mentioned > 0 && mentioned < totalBands * 0.4) {
        pattern_class = PATTERN_CLASS.SPLIT_TWO;
      } else {
        pattern_class = PATTERN_CLASS.SCATTER;
      }
      voteRecords.push({
        source_w: src_w, source_k: bandIdOf(src_w, src_c),
        visited_bands: visited,
        excluded_bands: excluded,
        pattern_class,
        n_source: srcSet.size,
      });
    }
  }
  return { voteRecords, bandIdRev, K_locus: bandIdRev.length };
}

// ---------- main ----------
const args = parseArgs(process.argv);
const raw = JSON.parse(fs.readFileSync(args.in, 'utf8'));

// Resolve candidate id
let candidateId = args.candidateId;
if (!candidateId && raw && typeof raw === 'object' && !Array.isArray(raw)) {
  candidateId = raw.candidate_id || null;
}
if (!candidateId) {
  candidateId = path.basename(args.in).replace(/\.[^.]+$/, '');
}

let voteRecords;
const derivation = { mode: args.mode, candidate_id: candidateId };

if (args.mode === 'voterecords') {
  if (Array.isArray(raw)) voteRecords = raw;
  else if (Array.isArray(raw.voteRecords)) voteRecords = raw.voteRecords;
  else throw new Error('voterecords mode: input must be an array or {voteRecords: [...]}');
  derivation.note = 'passthrough';
  derivation.n_records = voteRecords.length;
} else if (args.mode === 'projections') {
  const projs = raw.projections || raw;
  if (!Array.isArray(projs)) {
    throw new Error('projections mode: input must have {projections: [...]} or be an array');
  }
  voteRecords = voteRecords_from_projections(projs);
  derivation.note = 'reshaped via voteRecords_from_projections';
  derivation.n_source_bands = projs.length;
  derivation.n_records = voteRecords.length;
} else if (args.mode === 'kmeans') {
  const r = kmeansToVoteRecords(raw);
  voteRecords = r.voteRecords;
  derivation.note = 'Jaccard-based source-target projection; band_groups=window (see FINDING)';
  derivation.K_locus_inferred = r.K_locus;
  derivation.bandIdMap = r.bandIdRev;
  derivation.n_records = voteRecords.length;
  // Same-window pairs cannot be co-grouped (one cluster per arrangement
  // per window). Use window index as the group id.
  derivation.band_groups = r.bandIdRev.map(b => b.w);
} else {
  throw new Error(`Unknown mode: ${args.mode}`);
}

// Build option dicts only with explicitly-set values, so the band_tracking
// modules' own defaults remain the single source of truth.
const opts = {};
if (args.partnerThr   != null) opts.partnerSetOpts = { primaryThr: args.partnerThr };
if (args.ambiguousThr != null) opts.classifyOpts   = { ambiguousBandThr: args.ambiguousThr };
if (args.kCap         != null) opts.enumOpts       = { K_cap: args.kCap };
if (args.delta        != null) opts.selectOpts     = { delta: args.delta };
if (derivation.band_groups)    opts.band_groups    = derivation.band_groups;

const result = consensus_partition(voteRecords, opts);

// Strip typed arrays for JSON serialization
function toJSONSafe(o, depth = 0) {
  if (depth > 8) return '[truncated]';
  if (o == null || typeof o !== 'object') return o;
  if (o instanceof Float64Array || o instanceof Float32Array || o instanceof Int32Array) {
    return Array.from(o);
  }
  if (Array.isArray(o)) return o.map(x => toJSONSafe(x, depth + 1));
  if (o instanceof Map) {
    const m = {};
    for (const [k, v] of o.entries()) m[String(k)] = toJSONSafe(v, depth + 1);
    return m;
  }
  const out = {};
  for (const k of Object.keys(o)) out[k] = toJSONSafe(o[k], depth + 1);
  return out;
}

const outPath = args.out || `${candidateId}_consensus.json`;
const safeResult = toJSONSafe(result);
fs.writeFileSync(outPath, JSON.stringify({
  generated_at: new Date().toISOString(),
  args, derivation, result: safeResult,
}, null, 2));

// Human summary
console.log(`\n=== consensus_partition: ${candidateId} (mode=${args.mode}) ===`);
console.log(`  voteRecords:                 ${voteRecords.length}`);
console.log(`  K_locus:                     ${result.K_locus}`);
console.log(`  bruteforce_used:             ${result.bruteforce_used}`);
console.log(`  n_partitions_total:          ${result.n_partitions_total}`);
console.log(`  consensus_class:             ${result.consensus_class}`);
console.log(`  pca_vote_consensus_score:    ${result.pca_vote_consensus_score.toFixed(4)}`);
console.log(`  pca_hidden_regime_residual:  ${result.pca_hidden_regime_residual.toFixed(4)}`);
console.log(`  pca_partition_entropy:       ${result.pca_partition_entropy.toFixed(4)}`);
console.log(`  pca_overlap_conflict_score:  ${result.pca_overlap_conflict_score.toFixed(4)}`);
console.log(`  pca_resolving_power_class:   ${result.pca_resolving_power_class}`);
if (result.ambiguous_band_ids && result.ambiguous_band_ids.length) {
  console.log(`  ambiguous_band_ids:          [${result.ambiguous_band_ids.join(', ')}]`);
}
const topN = Math.min(5, result.top_partitions.length);
if (topN > 0) {
  console.log(`\nTop ${topN} partition${topN > 1 ? 's' : ''}:`);
  for (let i = 0; i < topN; i++) {
    const p = result.top_partitions[i];
    console.log(`  #${i+1}  score=${p.score.toFixed(4)}  blocks=${JSON.stringify(p.blocks)}`);
  }
}
console.log(`\nPer-band summary (top-10 by hidden_regime_residual desc):`);
const sorted = result.per_band.slice().sort((a, b) =>
  b.hidden_regime_residual - a.hidden_regime_residual);
for (const b of sorted.slice(0, 10)) {
  console.log(
    `  b${b.band}  residual=${b.hidden_regime_residual.toFixed(3)}  ` +
    `consensus=${b.voter_consensus.toFixed(3)}  ` +
    `evidence=${b.evidence_base.toFixed(2)}  ` +
    `primary=[${b.primary_partners.join(',')}]  ` +
    `secondary=[${b.secondary_partners.join(',')}]  ` +
    `${b.has_conflict ? 'CONFLICT' : ''}`
  );
}
console.log(`\nFull report: ${outPath}`);
console.log(`\nReasoning trace:`);
for (const r of result.reasoning) console.log(`  ${r}`);
