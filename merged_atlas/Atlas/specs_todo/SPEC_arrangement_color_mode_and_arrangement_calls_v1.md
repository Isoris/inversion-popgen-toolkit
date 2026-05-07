# SPEC — arrangement-color mode for per-sample-lines panel + arrangement_calls_v1 JSON

**Status**: drafted 2026-05-06
**Depends on**: shipped `band_tracking` modules (`vote_evidence.js`,
`band_voters.js`, `partition_enumerate.js`, `partition_consensus.js`,
plus the 2026-05-06 patches: `band_groups` correction + best-anchored
multi-layer detection)
**Consumes**: the existing `_resolveSampleColorByMode` /
`_resolveSampleScopeColor` dispatch infrastructure already wired into
`Inversion_atlas.html` (turn 130 Slice 2 era).
**Producer**: NEW — a small atlas-side runner OR a LANTA-side step
emitting `arrangement_calls_v1.json`. See §6.

---

## 0. Naming decisions still open

These three terms appear throughout this spec; before implementation,
pick a single canonical word for each. The spec body uses provisional
names tagged `[NAME-N]` so global rename is one find/replace.

| Tag | Provisional | What it means | Notes |
|---|---|---|---|
| `[NAME-1]` | **arr-cluster** | One K-means cluster within one probe (= one row of K_arrangement labels). Also called a "band" in the band_tracking math. | "band" already used in the math for the same thing — consider keeping "band" if no collision elsewhere |
| `[NAME-2]` | **candidate-region** | The genomic span being analyzed by `consensus_partition`, made up of ≥2 probes. | Existing pipeline already says "candidate"; adding "-region" disambiguates from "candidate as registry entry" |
| `[NAME-3]` | **probe** | The unit of clustering: a contiguous window-range on which K-means is run. Could be 1 window, 5 windows, an L2 envelope, half an L2 envelope, a custom bp range. | NEW concept — the existing pipeline has L2 envelopes but those are one specific kind of probe |

The collision Quentin flagged (2026-05-06): the current
`run_consensus_precomp.mjs` driver and the FINDING document use
"envelope" for both L1/L2 envelopes (genomic regions defined by
boundary-scan / Z-score insulation) AND for the analytic clustering
unit. They are different things. The L2 envelope is *one possible
choice* of `[NAME-3]`, not the only one.

## 1. Why a color mode (not a new panel)

The per-sample-lines panel (`drawLinesPanel()`, line ~34894 of
`Inversion_atlas.html`) is already the right place to display
arrangement calls. It has:

- One horizontal line per sample (226 lines, one row each)
- x-axis is window index along the chromosome
- already supports multiple coloring schemes via the
  `_resolveSampleColorByMode(si, gi, mode)` dispatch hook (lines 33256+)
- already has a UI mode picker (`linesColorModeSelect`, eight modes
  registered at line 33139)
- the v3.99 / turn 14e infrastructure left dispatch stubs and a
  `state.linesColorMode` slot ready for new modes

So: **arrangement coloring slots in as the ninth mode**, not a new
panel. Same canvas, same gridlines, same hover, same lasso, same
collapse-button — just one new color resolver function.

This is the lowest-LOC, highest-coherence way to surface the
band-tracking output. ~80–150 LOC into `Inversion_atlas.html`,
zero new panel scaffolding.

## 2. The atlas-side surface

### 2.1 New entry in `_LINES_COLOR_MODES`

```js
{ id: 'arrangement', layer: 'arrangement_calls', label: 'arrangement' },
```

Tooltip text (matching the style of existing modes):
> Arrangement calls from the band-tracking layer. Recolors each sample
> by the arrangement (= cross-probe sample partition) recovered by
> `consensus_partition` for the candidate-region under the cursor.
> Empty state until the `arrangement_calls_v1.json` layer is loaded.
> See SCHEMA §16.

### 2.2 New dispatch in `_resolveSampleColorByMode`

```js
if (mode === 'arrangement') {
  return _arrangementColor(si, gi);   // window-aware: gi == window index
}
```

And in `_resolveSampleScopeColor` (per-sample, no window axis):

```js
if (mode === 'arrangement') {
  return _arrangementScopeColor(si);  // uses currently-selected candidate
}
```

### 2.3 The two color functions

```js
function _arrangementColor(si, gi) {
  // gi is the window index along the chromosome.
  // 1. Find which candidate-region (if any) this window belongs to.
  //    state.candidateAtWindow[gi] → candidate_id, or null if outside any candidate.
  const cid = state.candidateAtWindow ? state.candidateAtWindow[gi] : null;
  if (cid == null) return null;     // outside any candidate → fall back to kmeans
  return _arrangementColorForCandidate(si, cid);
}

function _arrangementScopeColor(si) {
  // Per-sample, no window axis: use the currently-focused candidate
  // (state.activeCandidateId, set by the page-3 cursor / candidate selector)
  const cid = state.activeCandidateId;
  if (cid == null) return null;
  return _arrangementColorForCandidate(si, cid);
}

function _arrangementColorForCandidate(si, cid) {
  const calls = state.arrangementCalls;
  if (!calls || !calls[cid]) return null;
  const arr = calls[cid];
  const aid = arr.arrangement_per_sample[si];
  if (aid == null || aid < 0) return ARR_COLOR_UNCALLED;   // grey
  // n_arrangements is small (typically 2-6); use a fixed palette
  // in canonical order (deepest negative PC1 = arrangement 0).
  return ARR_PALETTE[aid % ARR_PALETTE.length];
}
```

### 2.4 The arrangement palette

Six deterministic colors, matching the existing K-means palette (so
the K=3 case looks identical to the kmeans mode at first glance,
which is the truth: when consensus is clean, arrangement IS K-means
post-Hungarian-aligned). Diverges only when band_tracking detects
something K-means alone couldn't see (e.g. cluster-id permutation
across probes).

```js
const ARR_PALETTE = [
  '#3074C8',   // arrangement 0 — deep blue (= existing K-means cluster 0)
  '#2BAA50',   // arrangement 1 — green
  '#D04545',   // arrangement 2 — red
  '#A060B8',   // arrangement 3 — purple
  '#D8A030',   // arrangement 4 — gold
  '#3DB5C0',   // arrangement 5 — teal
];
const ARR_COLOR_UNCALLED = 'rgba(140, 140, 140, 0.35)';
```

Beyond 6 arrangements: golden-angle rotation, same as `_lineageColor`.

### 2.5 New state slots

```js
state.arrangementCalls = null;     // keyed by candidate_id
state.candidateAtWindow = null;    // Int32Array(n_windows), set on candidate load
```

Persisted via the existing `state_io.js` machinery — both slots
round-trip localStorage on the same key the rest of `state` uses.

### 2.6 New JSON layer registration in `detectSchemaAndLayers`

Add `'arrangement_calls'` to the layers presence set; the mode-picker
auto-enables when the layer is loaded, auto-greys with tooltip
"Needs `arrangement_calls_v1.json` JSON layer" when missing — same
pattern as every other gated mode.

## 3. The JSON contract: `arrangement_calls_v1.json`

One file per chromosome (like `LG28.json` in the precomp). Loaded
into `state.arrangementCalls` keyed by candidate_id.

```jsonc
{
  "tool": "arrangement_calls_v1",
  "schema_version": 1,
  "chrom": "C_gar_LG28",
  "n_samples": 226,
  "generated_at": "2026-MM-DDTHH:MM:SSZ",
  "candidates": {
    "C_gar_LG28_cand_15Mb": {
      // The candidate-region [NAME-2] this entry describes.
      "candidate_id": "C_gar_LG28_cand_15Mb",
      "start_bp": 15115000,
      "end_bp": 18005000,
      "start_w": 3253,
      "end_w": 3970,

      // The probes [NAME-3] that voted on this candidate. Each probe
      // is the unit of K-means clustering. Probes can be:
      //   - L2 envelopes (most common — { kind: "l2_envelope", l2_id: "..." })
      //   - L1 envelopes (rare; for very wide regions)
      //   - sliding window groups ({ kind: "windows_n5" | "windows_n10" }, etc.)
      //   - explicit bp ranges ({ kind: "custom_range", start_bp, end_bp })
      // Order: chromosome order (start_bp ascending).
      "probes": [
        {
          "probe_idx": 0,
          "kind": "l2_envelope",
          "l2_id": "C_gar_LG28_d17L2_0010_01",
          "start_w": 3224, "end_w": 3252,
          "start_bp": 15087364, "end_bp": 15208345,
          "k_used": 3,
          "centers": [-0.082, 0.006, 0.072],     // sorted ascending (PC1 K-means convention)
          "n_per_arr_cluster": [62, 102, 62],
          "ok": true, "reason": null
        },
        // ... one entry per probe
      ],

      // Cross-probe arrangement assignment per sample.
      // arrangement_per_sample[si] ∈ [0, n_arrangements) ∪ { -1: uncalled }
      // -1 means: sample's arr-cluster identity didn't track cleanly
      // across probes (e.g. it switched arrangement between probes —
      // possible recombinant, or coverage gap, or noise).
      "n_arrangements": 3,
      "arrangement_per_sample": [0, 0, 0, /* ... */ 1, 1, 1, /* ... */ 2, 2, /* ... */],

      // Per-arrangement size + the order convention. By default,
      // arrangement 0 is the one whose PC1-mean across the candidate
      // is most negative (= "lowest" K-means cluster). This makes
      // arrangement labels stable across reruns.
      "arrangement_sizes": [60, 106, 60],

      // Per-(probe, arr_cluster) → arrangement mapping.
      // arr_cluster_to_arrangement[probe_idx][k_within_probe] = arrangement_id
      // This is the "decoded permutation" — what consensus_partition
      // recovered for which K-means label at which probe maps to
      // which arrangement.
      // Same shape regardless of n_probes; arrangement_id ∈ [0, n_arrangements).
      "arr_cluster_to_arrangement": [
        [0, 1, 2],   // probe 0: cluster 0 → arrangement 0, cluster 1 → arrangement 1, cluster 2 → arrangement 2
        [2, 0, 1],   // probe 1: cluster 0 → arrangement 2 (permuted!), cluster 1 → arrangement 0, cluster 2 → arrangement 1
        [1, 2, 0],   // probe 2: ...
        // ...
      ],

      // QC + classification fields from consensus_partition.
      // These get displayed inline (a small badge under the panel header).
      "consensus_class": "CLEAN_PARTITION",
      "pca_vote_consensus_score": 0.8667,
      "pca_hidden_regime_residual": 0.2667,
      "pca_partition_entropy": 0.0,
      "pca_overlap_conflict_score": 0.0,
      "pca_resolving_power_class": "COMPLEX_BUT_RESOLVABLE",
      "ambiguous_band_ids": null,    // or [int, ...] for AMBIGUOUS_BAND class
      "bruteforce_used": true,
      "n_partitions_total": 21147,
      "reasoning": [/* full trace from consensus_partition */]
    }
    // ... one entry per candidate-region
  }
}
```

Schema-version 1. If the band_tracking algorithm changes in a way
that affects field semantics, bump the version.

## 4. The producer side: how `arrangement_calls_v1.json` gets written

Two production paths, equivalent output:

### 4.1 Path A — atlas-side runner (browser)

A new module `arrangement_call_runner.js` that, given the precomp
JSON (already loaded as `state.data`) and a candidate registry entry,
runs:

1. Pick probes for the candidate (default: L2 envelopes inside it;
   user can override via UI).
2. For each probe, run `clusterL2`-equivalent K-means on per-sample
   mean PC1 (the existing `per_l2_cluster.js` already does this).
3. Build voteRecords + band_groups from the per-probe cluster labels
   (Jaccard-based projection — the same logic in
   `run_consensus_precomp.mjs:locusToVoteRecords`).
4. Call `consensus_partition(voteRecords, { band_groups })`.
5. Decode the result's top partition into per-sample arrangement
   ids (the `arr_cluster_to_arrangement` table). The decoder follows
   the convention: arrangement 0 = the block whose mean PC1 across
   the candidate is most negative.
6. Stash into `state.arrangementCalls[candidate_id]`.

Triggered by:
- A button on the candidate-page (page 3) cursor: "Run band-tracking
  on this candidate".
- Auto-run on first focus of any candidate in the registry, if not
  already cached.

Cost: ~5–50 ms per candidate (bruteforce dominates; K_locus ≤ 10
caps it at 115k partition evaluations, all cheap).

### 4.2 Path B — LANTA-side batch step

A standalone Node/Python script that reads the precomp JSON +
candidate registry + a probe-strategy spec, writes the JSON. Pros:
runs once, no browser cost; fits the existing
`STEP_C01f_*` pipeline shape. Cons: requires re-running when
candidates change.

Suggested step name: `STEP_C01f_g_emit_arrangement_calls.py` (same
naming convention as the BUSCO 4D step from the parallel age-spec
work). ~150 lines of Python that calls a small Node CLI wrapping
`consensus_partition` (or a Python port of the math — but the JS
implementation is the reference).

Recommendation: ship Path A first (browser-side, reactive, no LANTA
dependency), then add Path B for batch / reproducibility / manuscript
bundle exports.

### 4.3 Probe-strategy choices

The producer decides how to slice the candidate-region into probes.
The spec is agnostic; default strategy:

1. If the candidate is contained in a single L2 envelope: use that
   envelope as the only probe → consensus_partition has only 1 probe
   to vote with → no cross-probe voting → falls through to the
   per-band view only (which is informative).

   *Special handling*: when `n_probes == 1`, the producer can OPT to
   subdivide the L2 into 2–3 sub-probes (window groups of size n=ceil(N/2)
   or n=ceil(N/3)) so band_tracking has cross-probe evidence.
   Configurable.

2. If the candidate spans 2+ L2 envelopes: each L2 is a probe. This
   is the "natural" multi-probe case.

3. User override: explicit list of probes via UI (page-3 panel
   "probes for this candidate" sub-control), or via a JSON sidecar
   alongside the registry entry.

The probe-strategy choice IS recorded in the output JSON
(`probes[i].kind` field), so downstream consumers know what the unit
was.

## 5. Probes that subdivide an L2 envelope

Quentin's flag (2026-05-06): "in our logic we could also split L2 if
necessary so it should be more like using our old system per window
1, 5, 10 and also per L2 idk but its really not only per L2 only
because it may be a bit off."

This is correct. The `[NAME-3]` abstraction handles it:

- A probe of `kind: "windows_n5"` is a contiguous group of 5 windows.
  K-means on per-sample mean PC1 across those 5 windows.
- A probe of `kind: "windows_n10"` is 10 windows.
- A probe of `kind: "l2_envelope"` is whatever the L2 envelope spans
  (variable, often 100s of windows).
- A probe of `kind: "custom_range"` is an explicit `[start_bp,
  end_bp]` chosen by the user.

All four kinds produce the same shape of output: K cluster labels
per sample, mean PC1 per cluster (centers), n_per_cluster. Therefore
all four are interchangeable inputs to band_tracking. The producer
mixes-and-matches as appropriate for the candidate.

For the LG28 prototype: the right default is probably "L2 envelopes
01–06 of L1 region 0010" (six probes, K_locus=18 — exceeds bruteforce
cap, so until greedy fallback is wired in, restrict to 3 probes at a
time). With sub-probes, the prototype could be sliced into 6 × 3 =
18 short probes at K_locus=18 still, OR into a uniform "every 100
windows" tiling of ~7 probes at K_locus=21 — same cap problem,
better spatial resolution.

The takeaway: probe strategy is **a tunable parameter of the
producer**, not a fixed pipeline. The spec defines the wire format;
the producer's UX defines what defaults work for which discovery
scenarios.

## 6. UI surface (page 3 candidate panel)

Two small controls added to the candidate panel header (no new
panels):

### 6.1 "arrangement" mode in the existing color-mode picker

Already specified in §2 — auto-greyed when the layer is missing,
auto-enabled when loaded. The mode-picker is the user's "show me
arrangements" toggle.

### 6.2 Probe-strategy mini-control (expander)

A small expander labeled **probes** under the candidate panel header,
collapsed by default. When expanded, shows:

```
probes for this candidate:
  ◉ L2 envelopes (default — 6 probes)
  ○ sliding windows  [n=  5 ▾ ]   →  preview: 12 probes
  ○ sliding windows  [n= 10 ▾ ]   →  preview:  6 probes
  ○ custom — opens probe-edit dialog
  [recompute] [save as default for this candidate]
```

Selecting a strategy + recompute re-runs the producer on the active
candidate, updates `state.arrangementCalls[cid]`, redraws the
panel. The choice is persisted into the candidate registry so it
sticks.

### 6.3 QC badges under the panel header

When the arrangement mode is active for the focused candidate, a
small status strip displays the consensus_partition QC fields,
inline:

```
arrangement: CLEAN_PARTITION  ·  vote_consensus 0.87  ·  residual 0.27
              [details ▾]
```

Click the "details" disclosure to expand the full reasoning trace +
top-N partition list. This is the diagnostic view for when the
class isn't CLEAN.

The QC strip is small (~24 px tall when collapsed). Not a new
panel, just header chrome on the existing one.

## 7. Decoding `consensus_partition` output into arrangement IDs

The band_tracking core returns top partitions like
`[[0,3,6],[1,4,7],[2,5,8]]` — block ids are arbitrary (the order
they were first encountered during enumeration). The atlas needs a
**deterministic, manuscript-stable** mapping from block id to
arrangement id.

Convention (from the spec's invariant "arrangement 0 has lowest
mean PC1"):

1. Run consensus_partition. Pick the top partition.
2. For each block, compute the mean PC1 of its arr-clusters (using
   the K-means centers from each probe).
3. Sort blocks by ascending mean PC1.
4. Assign arrangement_id 0 to the most-negative-mean block, 1 to
   the next, etc.
5. Build `arr_cluster_to_arrangement[probe][k_within_probe]` from
   the (probe, k) → block_id mapping in the chosen partition.
6. Build `arrangement_per_sample[si]`: for each sample, look up its
   K-means cluster label at every probe; the sample's arrangement is
   the one that the majority of its (probe, cluster) lookups map to.
   If no majority (sample's clusters cross arrangements between
   probes), label as `-1` (uncalled — likely recombinant or noisy).

A sample that lands in arrangement A at every probe gets `arrangement_per_sample[si] = 0` (clean). A sample that lands in arrangement A at probe 0 but arrangement B at probe 1 gets `-1` (uncalled, flagged for review).

This decoder is ~40 LOC and can live in a new
`arrangement_decode.js` module in `Atlas/shared/band_tracking/`.

## 8. The label-ambiguity case (F12) — UX

When n_probes = 2 and K_arrangement = 3 (or any case where multiple
top partitions tie at score 1.0 due to label ambiguity), the
classifier returns MULTI_LAYER_STRUCTURE (per fixture F12).

UX response: when consensus_class is MULTI_LAYER_STRUCTURE AND
all top partitions tie at score = 1.0 AND residual = 0, the QC strip
shows:

```
arrangement: LABEL_AMBIGUOUS  ·  3 valid partitions tie at 1.0  ·  pick one ▾
              [needs 3rd anchor probe to disambiguate]
```

The user can pick which of the tied partitions to display (the
"pick one" dropdown lists them). The atlas does NOT silently pick
one and pretend the ambiguity is gone.

This is a distinct visual state, not a 7th classifier class — the
algorithmic class stays MULTI_LAYER_STRUCTURE; only the UI label
softens it. (The decision to add a 7th class
`LABEL_AMBIGUOUS_CLEAN_PARTITION` to the algorithm itself is parked
— see open decision §11.1.)

## 9. The "uncalled" sample marker

Samples with `arrangement_per_sample[si] = -1` (didn't track cleanly
across probes) are colored `ARR_COLOR_UNCALLED` (grey, alpha 0.35)
in the per-sample-lines panel. This is intentionally low-contrast
so they don't dominate the visual.

A small counter appears in the QC strip:
```
arrangement: CLEAN_PARTITION  ·  vote_consensus 0.87  ·  3 of 226 uncalled
```

Three uncalled samples out of 226 is fine. Thirty would be a flag
that something's off — possible recombinants, possible probe
boundary crossing the breakpoint, etc. The user can click the
counter to highlight just those samples in the panel.

## 10. Implementation order (next chat)

| # | Step | LOC est | Dependencies |
|---|---|---|---|
| 1 | Naming decisions (§0) → global rename across spec, FINDING, drivers, tests | ~50 | nothing |
| 2 | `arrangement_decode.js` module: top-partition → arrangement_per_sample + arr_cluster_to_arrangement | ~80 | shipped band_tracking |
| 3 | `arrangement_call_runner.js` (Path A): wraps clusterL2 + voteRecords assembly + consensus_partition + decode | ~150 | shipped band_tracking, existing per_l2_cluster.js |
| 4 | `_LINES_COLOR_MODES` registration + `_arrangementColor` + `_arrangementScopeColor` + state slots | ~80 | none |
| 5 | Probe-strategy mini-control (page 3 expander) | ~100 | runner from step 3 |
| 6 | QC badges under candidate panel header | ~60 | runner from step 3 |
| 7 | Label-ambiguity UX (F12 case) | ~40 | step 6 |
| 8 | Tests for the decoder + runner (synthetic + real-LG28-precomp fixtures) | ~200 | steps 2-3 |
| 9 | LANTA-side batch step (Path B) | ~150 Python | optional, not blocking |

**Total atlas-side**: ~510 LOC + ~200 LOC tests. About 60% the size
of the age-layer patch (which was ~880 LOC).

**Critical path**: 1 → 2 → 3 → 4. Steps 5–9 are independent
follow-ups. The minimum viable visualization is steps 1+2+3+4:
a working "arrangement" mode showing the band_tracking output for
any candidate that has L2-envelope-default probes.

## 11. Open decisions

### 11.1 LABEL_AMBIGUOUS as a 7th classifier class?

When `pca_vote_consensus_score == 1.0`, `pca_hidden_regime_residual
== 0`, and ≥2 picked partitions tie at score 1.0, the data is
mathematically saturated but cross-probe label assignment is
ambiguous (e.g. F12 case: 2 probes, K=3 arrangements → 6 valid
permutations). Currently classifies MULTI_LAYER_STRUCTURE. UX
softens this (§8). Algorithm-level: should this become its own
class?

Pros: cleaner semantics; downstream consumers know the difference.
Cons: 7-class taxonomy harder to remember; the spec's UX-level
softening already handles the user-facing case.

Recommendation: defer until a non-test case actually fires
LABEL_AMBIGUOUS and the UX softening proves insufficient.

### 11.2 Auto-subdivide single-L2 candidates?

If a candidate is contained in 1 L2 envelope, default behavior is
"can't run band_tracking, fall through to K-means". Optional:
producer auto-subdivides into 2–3 sub-probes. Convenient but adds
producer-side complexity. UI flag for whether to auto-subdivide?

Recommendation: ship without auto-subdivide; let user pick via the
probe-strategy expander (§6.2) when they need it.

### 11.3 Probe-strategy persistence

The probe-strategy choice (§6.2) is per-candidate. Where does it
persist? Three options:
(a) candidate registry sidecar JSON
(b) per-user localStorage
(c) part of `arrangement_calls_v1.json` itself (the producer
    records what strategy it used)

Recommendation: (c) — the producer records its choice in the output;
the UI defaults to whatever the last successful run used. No new
sidecar files.

### 11.4 Color palette match with existing K-means palette

§2.4 proposes the arrangement palette use the same first 3 colors
as the existing K-means palette so the K=3 clean-consensus case
looks visually identical to the kmeans mode. Confirm this is
desired (vs. a deliberately distinct palette so the user can
visually tell the modes apart at a glance).

Recommendation: match the K-means palette for arrangement IDs 0–5,
diverge only via golden-angle rotation past 6.

### 11.5 What does "arrangement" mean when consensus_class != CLEAN?

If the class is OVERLAPPING_VOTES, AMBIGUOUS_BAND, or
NO_CLEAN_CONSENSUS, do we still emit `arrangement_per_sample`? Two
options:
(a) Always emit using the top partition; let the QC class warn the
    user that the call may be unreliable.
(b) Only emit when class ∈ {CLEAN_PARTITION, SOFT_PARTITION}; for
    other classes, set `arrangement_per_sample` to all -1.

Recommendation: (a) — always emit. Greying out the panel hides
information that may still be useful for review. The QC class +
residual numbers tell the user how much to trust it.

## 12. Cross-references

- `band_voters.js`, `partition_consensus.js`, `vote_evidence.js`,
  `partition_enumerate.js` — the band_tracking math (shipped, with
  2026-05-06 patches).
- `Atlas/shared/per_l2_cluster.js` — provides `clusterL2` /
  `clusterL2AtK`, the atlas-side K-means per L2 envelope. The
  arrangement runner reuses this for L2-kind probes.
- `Atlas/shared/kmeans.js` — `kmeans1D`, used directly for
  windows_n5 / windows_n10 / custom_range probes that don't have an
  existing per-L2 cache.
- `_resolveSampleColorByMode`, `_resolveSampleScopeColor`,
  `_LINES_COLOR_MODES`, `refreshLinesColorMode` — existing
  per-sample-lines color-mode infrastructure (lines ~33139–33330 of
  `Inversion_atlas.html`).
- `FINDING_2026-05-06_both_excluded_bug.md` — the bug fix that
  enables this spec to deliver correct K≥3 classifications on real
  data. Without that fix, every K=3 candidate would render with
  MULTI_LAYER_STRUCTURE QC despite being clean.

## 13. What this is NOT

- **Not a new panel**. Reuses `drawLinesPanel()` end to end. Only
  adds: 1 color mode entry, 2 dispatch hook bodies, 1 mini-control
  expander, 1 QC strip.
- **Not a replacement for the K-means mode**. `kmeans` and
  `arrangement` modes coexist. K-means shows raw labels (cluster id
  permutes between probes); arrangement shows recovered identity
  (constant across probes).
- **Not a per-probe re-clustering UI**. The atlas does not let the
  user tune K-means parameters from the probe-strategy expander —
  that already lives elsewhere (sidebar K selector, etc.). The
  expander only chooses **which probes** to feed to band_tracking,
  not how K-means is run inside each probe.
- **Not L2-locked**. Probes can be L2 envelopes, but they can also
  be sliding window groups, half-L2 splits, or custom bp ranges.
  The producer is flexible by design (see §5).
- **Not a fix for the K_locus > 10 bruteforce limit**. That's a
  separate parked item — needs the greedy fallback
  (`kt_infer_macro_band_groups` from band_layers tarball #3) wired
  in. The arrangement spec works at K_locus ≤ 10; for larger
  candidates the runner will refuse and the QC strip will say so.
- **Not a fix for the F12 label-ambiguity case at the algorithmic
  level**. The UX softening (§8) handles user-facing display; the
  algorithm still classifies MULTI_LAYER_STRUCTURE.

End of spec.
