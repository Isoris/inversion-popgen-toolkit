# OPTION B finding — `both-excluded → together` is unsafe when K_arrangements > 2

**Date**: 2026-05-06
**Confirmed on real data**: 2026-05-06 (same day, follow-up run).

**Synthetic run**: `synth_LG28_kmeans.json` — 226 samples (60/106/60
LG28 prototype shape), 3 windows, K=3 clusters per window, with
window-permuted cluster ids.
**Real-data run**: `LG28.json` precomp from
`Atlas_turn166_round2_2026-05-05_tar-1.gz` — full LG28 chromosome,
4302 windows, 226 samples, 48 L2 envelopes; tested on the prototype
interval 15.115–18.005 Mb (envelopes _0010_01 through _0010_06).
**Drivers**: `run_consensus.mjs --mode kmeans` (synthetic),
`run_consensus_precomp.mjs` (real precomp).

## Headline

The K=3 arrangement structure is **mathematically recovered** on both
synthetic and real data — the top partition correctly groups the
HOM_REF / HET / HOM_INV cluster ids across windows / envelopes despite
K-means cluster-id permutations. The Jaccard-based collapse (synthetic
driver) and the cross-envelope projection (real-data driver) both
work.

But the **classification ceilings at `MULTI_LAYER_STRUCTURE`** with
`pca_vote_consensus_score ≈ 0.67–0.73` (not 1.0) and
`pca_hidden_regime_residual ≈ 0.53–0.67` (not 0.0). The class is wrong:
this is a CLEAN three-arrangement structure, not a multi-layer locus.

## Real-data confirmation (added 2026-05-06)

Running `run_consensus_precomp.mjs` on the LG28 precomp JSON, with
locus = pairs/triples of L2 envelopes around the prototype interval:

| Locus | n_env | K_locus | top blocks | score | class |
|---|---|---|---|---|---|
| envs 01+02 | 2 | 6 | `[[0,3],[1,4],[2,5]]` | 0.6667 | MULTI_LAYER_STRUCTURE |
| envs 02+03+04 | 3 | 9 | `[[0,3,6],[1,4,7],[2,5,8]]` | 0.7333 | MULTI_LAYER_STRUCTURE |
| envs 04+05 | 2 | 6 | `[[0,3],[1,4],[2,5]]` | 0.6667 | MULTI_LAYER_STRUCTURE |

Per-envelope K-means cluster sizes (registered prototype is 60/106/60):

- _0010_01: 62/102/62 ✓
- _0010_02: 61/103/62 ✓
- _0010_03: 62/104/60 ✓
- _0010_04: 71/95/60  ≈
- _0010_05: 91/67/68  ✗ (anomalous — possible breakpoint/regime change at 17.46–17.94 Mb)
- _0010_06: 76/96/54  ≈

The band-tracking layer's bruteforce **still ranks the correct K=3
partition first** on every pair/triple, with the score
ceiling (0.6667 / 0.7333) identical to the synthetic K=3 fixture's.

The bug fires identically on real fish data.

## Operational note: K_locus exceeds bruteforce K_cap on the full prototype

The full prototype (envelopes 01–06, six envelopes × K=3 clusters =
**K_locus=18**) exceeds the bruteforce hard_fail at K=13 and the
default K_cap at K=10. The handoff's queued fallback is
`kt_infer_macro_band_groups` (greedy, in the not-yet-shipped
band_layers tarball #3). Until that's wired into the driver, the
locus-centric classification only runs on smaller groupings (single
envelope = no consensus to compute since voting is cross-envelope;
2–3 envelopes = bruteforce works).

For the manuscript: the per-band view (Layer A) runs at any K_locus
and shows the K=3 structure cleanly across all 18 bands of the
6-envelope prototype (verified 2026-05-06). The locus-centric Layer C
classification needs either greedy fallback OR same-arrangement-merging
across envelopes to push K_locus down before the score is meaningful.

## Anomaly worth flagging: envelope `_0010_05` (17.46–17.94 Mb)

This envelope's K=3 split is 91/67/68, not the expected ~60/106/60.
Either there's a real signal change in this sub-region (possibly the
breakpoint vicinity — the registered prototype ends at 18.005 Mb),
or PC1 has a different shape there that pushes the K-means clustering
into a different local optimum. The Jaccard-based cross-envelope
projection still groups its bands with the correct counterparts in
other envelopes (verified — band 12 from this envelope is in the
[[0,3,6,9,12,15],...] group with the other "first cluster" bands), so
sample identities are preserved across the size-shape change. But
the cluster-size mismatch is real and worth investigating in the
manuscript.

## Root cause

`build_coassociation_matrix` in `vote_evidence.js` uses the rule:

```js
if ((iVis && jVis) || (iExc && jExc)) together[idx] += w;
else                                    apart[idx]    += w;
```

(lines 199–203). The "both excluded → together" arm assumes the vote
treats `excluded_bands` as a single coherent "other side." That holds
in K=2 arrangements: a HOM_REF voter's excluded set IS uniformly
HOM_INV (the only other side), so any two bands in the excluded set
genuinely are alike.

In K=3 arrangements (HOM_REF / HET / HOM_INV — i.e., the LG28
prototype with 60/106/60 karyotypes), a HOM_REF voter excludes both
HET clusters AND HOM_INV clusters in the other windows. The rule
counts them as co-grouped, but they're actually different
arrangements that the voter just doesn't distinguish.

## Coassociation matrix evidence

For the synthetic K=3 fixture, the matrix shows:

```
        0       1       2       3       4       5       6       7       8
  0   1.00/ 0  0.33/ 6  0.33/ 6  1.00/ 3  0.33/ 3  0.33/ 3  1.00/ 3  0.33/ 3  0.33/ 3
  1   0.33/ 6  1.00/ 0  0.33/ 6  0.33/ 3  1.00/ 3  0.33/ 3  0.33/ 3  1.00/ 3  0.33/ 3
  ...
```

(format: `coassoc / total_weight`)

- Same-arrangement pairs (e.g. 0–3, 0–6, 3–6 = HOM_REF across windows):
  `coassoc = 1.00` ✓ — correct.
- Cross-arrangement pairs (e.g. 0–1 = HOM_REF / HET in window 0):
  `coassoc = 0.33` ✗ — should be 0 (or unmentioned). The 0.33 = 2 / 6:
  two of the six votes mentioning both bands place them on the same
  excluded side; four place them apart.

The "two together votes" come from the third-arrangement source
(HOM_INV from the other windows) which excludes both 0 (HOM_REF) and
1 (HET) jointly. That's the bug surface.

## Why the existing 991-test suite didn't catch this

Every fixture in `test_band_consensus.js` uses K=2 arrangements
(visited = {0,1,2}, excluded = {3,4,5}). The synthetic test data
never includes a vote whose excluded set spans two distinct
arrangements. So the bug never fires in tests, and the score ceiling
of ~0.73 in K=3 cases is invisible.

## Implications

1. **The algorithm core is sound.** The bruteforce enumeration still
   ranks the correct partition at the top — the K=3 answer
   `[[0,3,6],[1,4,7],[2,5,8]]` outscores all 21,146 alternatives.
   What's wrong is the absolute score level and the downstream
   classification.

2. **Calibration thresholds are not the right knob.** Lowering
   `cleanScoreThr` from 0.85 to 0.70 would re-class this run as
   CLEAN_PARTITION but it'd also re-class genuinely soft cases as
   clean. The fix has to be upstream in coassociation, not downstream
   in thresholds.

3. **All synthetic test fixtures need a K=3 sibling.** Once the fix
   is in, the test suite should pin "K=3 arrangements with permuted
   cluster ids → CLEAN_PARTITION at score ≥ 0.95" as a regression test.

4. **The Option A age-layer patch can proceed independently.** It
   touches none of the band-tracking code; the bug is at the K_locus
   layer, not in age estimation.

## Proposed fix (sketch — for next session, not this turn)

Two options:

### Fix option 1 — partition the excluded set by source

Each vote already knows its source band (`source_k`). When projecting,
the `excluded_bands` set could be subdivided by which arrangement they
belong to (recoverable from purity_vector or from the upstream
projection sweep's per-target-band assignment). Then the rule becomes:

```js
if ((iVis && jVis) || (iExc && jExc && i and j in same excluded subgroup))
   together[idx] += w
else
   apart[idx] += w
```

Requires the per-target-band assignment to be carried through
voteRecords. Cleanest semantics; needs upstream pipeline to provide
the subgrouping.

### Fix option 2 — drop "both excluded → together" entirely

Treat `excluded` as evidence of "I think these are not me," not "I
think these are each other." Counts every excluded-side pair as
`apart`.

**Tested on the same K=3 fixture.** Result: **fix option 2
over-corrects.** Without the both-excluded reinforcement:

- coassoc[0,1] (cross-arrangement, was 0.33) → 0.00 ✓ correct
- coassoc[0,3] (same-arrangement HOM_REF, was 1.00) → 0.33 ✗ broken

Why: same-arrangement evidence for (3, 6) comes from two sources:
the HOM_REF voter's "both visited" (1 contribution) AND the 6 other
voters' "both excluded" (6 contributions). Stripping the second arm
removes 6/7 of the same-arrangement signal. Bruteforce now ranks the
fully-singleton partition slightly above the K=3 truth (best score
0.9333 vs 0.9259 for any 2-merger). The correct K=3 partition is no
longer at the top.

So **fix option 2 is wrong**.

### Fix option 3 — discount weight on both-excluded contributions

Count `together_excluded` with weight w' < 1. The K=2 case ceiling
falls modestly (probably 1.0 → ~0.95 at w'=0.5); the K=3 case rises
modestly (0.73 → ?). This tunes between under- and over-counting but
adds a parameter and doesn't fix the underlying ambiguity.

### Fix option 4 — same-window subgrouping

Within a single voter, partition `excluded_bands` into per-target-window
groups. Bands from the same target window are by construction in
different arrangements (a window has one cluster per arrangement), so
two excluded bands from the same target window should be `apart`.
Bands from different target windows in the same voter's excluded set
... could still be the same arrangement (HET in window 1 and HET in
window 2 both excluded by HOM_REF source). So same-window-pair → apart,
cross-window-pair → uncertain.

This requires the voter to know the target's window structure, which
is fine — it's already implicit in the per-band → per-window mapping
that the projection sweep produces.

## Recommendation

**Fix option 1 (per-source subgrouping) is the right fix** if the
upstream pipeline can carry per-target-band arrangement assignment.
Cleanest semantics; no information loss.

**Fix option 4 (same-window subgrouping)** is the right fix if it
can't, because it uses information already present in the voter's
view (target windows are known) without requiring additional pipeline
output.

Either way, **fix option 2 is ruled out empirically** by this run.
Both options 1 and 4 preserve the "same-arrangement reinforcement"
evidence (bands 3 and 6 still co-associate at 1.0 because HOM_REF
voters visit both, HET voters exclude both from the same arrangement,
HOM_INV voters exclude both from the same arrangement — all three
correctly count toward together).

For now, **the finding is the deliverable**: bug localized to one
rule in `build_coassociation_matrix`; naive fix tested and ruled out;
two viable fix paths identified.

## Caveats

This run used the **kmeans-mode driver**, which itself uses a
Jaccard-based source-target collapse (see comments in `run_LG28.mjs`).
The collapse may interact with the underlying bug in ways that
amplify or mute it on real data. The same-vote-shape problem
(`excluded` lumps unlike arrangements) would also fire in
projections-mode if the upstream `bp_project_bandset_to_target_bands`
emits multi-arrangement excluded sets, which it almost certainly does
for K=3.

To be sure, the next step on real data should run **both** modes:
1. The actual production projections JSON via `--mode projections`.
2. The raw K-means labels via `--mode kmeans`.

If both classify the LG28 prototype as something other than
CLEAN_PARTITION when the locus is a clean K=3 inversion, the fix is
needed regardless of input mode.
