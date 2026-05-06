# SPEC — Single-band contingency rows in L3 panel + automated band-track extraction

**Status**: forward-looking spec. Not implemented. Atlas-side feature
extending the existing L3 contingency tables panel. Implementation
estimate: 2.5–3.5 turns split across slices.

**Extended this session** with §1.4–1.5 (het band as principal
track + lines-panel as the band-track view) and §4A (het-anchored
mode with carrier-continuity skeleton extraction, internal-dropout
tolerance, and the §4A.5 explanation of why this lives inside the
symmetric-mode spec rather than a separate one). §6 gains TSVs E
(het skeletons) and F (het carrier edges); §7 gains §7.4
(lines-panel `band_track_*` color modes) and §7.5 (skeleton
overlay); slices and tests extended to match. §3.2 adds the
HPC-conceptual ↔ atlas-primitive input-column mapping so a future
R-module mirror has a stable input contract. §4A.6 records the
seven possible biological + technical causes of het-signal
dropout (recombination, gene conversion, marker dropout, local
low diversity, shared ancestral haplotype, ancestry mosaicism,
noisy K-means, lack of informative SNPs); §4A.7 records the
canonical manuscript slogans **"dropout does not necessarily mean
boundary"** and **"different downstream carriers = different
system."** The spec now covers everything from per-band rows in
the L3 K×K view to the full het-skeleton-driven candidate
skeleton + lines coloring + manuscript-grade vocabulary.

**Reading order**:
- this spec → `SPEC_distant_band_concordance_fish_trajectory.md`
  (the dual representation — fish paths instead of cohort paths;
  Slices 1+2 already shipped in atlas turn 130)
- `SPEC_l2_sweep_inheritance.md` (system-level merger this layer
  sits inside)
- band-trace work shipped turns 160-164 (single-cohort tracker;
  the seed primitive for this spec's interactive mode)

**Trigger** (Quentin, this session):
> *"In complex regions, a structural system may span a broad interval,
> but individual sample bands may start/end at different positions.
> Right now I manually pick one K-means band and visually track whether
> the same samples stay together across windows. This works, but it is
> manual and not reproducible. I want this automated."*
>
> *"Be careful with double crossovers or A-B-A patterns. A band can
> disappear in the middle of a genomic region and then reappear later.
> This does not necessarily mean it is a completely unrelated signal.
> It may mean the same cohort returns after an internal block."*
>
> *"Its really to upgrade our L3 contingency table panels of atlas.
> To also have single bands contingencies smth like that. Idk how to
> explain. Its complementary."*
>
> *"I am not sure but its like the per sample lines but with lines
> colored by their 'intervals'."*
>
> *"Add this to our single band. Its basically we track the het band
> from dosage bc its the skeleton for our candidate inversion
> interval … focus specifically on the heterozygous/divergent band
> and track sample identity across windows. The inversion 'skeleton'
> or 'backbone' is the chain of windows/blocks where the same
> heterozygous carrier set persists. Do not define the backbone only
> by visual band shape. Define it by sample identity continuity."*

---

## 1. The problem this solves

The existing L3 contingency panel renders, for each adjacent
window/L2 pair, a full K×K table with row/column sums plus the
secondary metric (Cramér's V / NMI / AMI / ARI) summarising the
table as one number. This is **system-level**: it answers *"are
these two windows similar enough to merge?"*.

What it does NOT show: **per-band persistence**. Inside a single
K×K table, each row is a cohort whose fate to the next window
might be very different from the cohort one row over:

```
                     window w+1
                  c0    c1    c2    c3    c4    c5
window w  f0     56    0     0     0     0     0     ← cohort stays put: g0→c0
          f1     2     38    0     0     1     0     ← cohort stays put: g1→c1
          f2     0     0     12    18    0     0     ← cohort SPLITS: g2 → {c2, c3}
          f3     0     0     0     2     0     33    ← cohort jumps label: g3→c5
          f4     ...
          f5     ...
```

The K×K table shows this if you read each row carefully, but the
single-number metric (Cramér's V = 0.78) collapses it into "yes,
similar." The user has to eyeball every row to see that g2 is
splitting while g0/g1/g3/g5 are stable.

**Quentin's manual workflow** today is to pick one K-means band
(say g2), follow the same fish-set forward by clicking through
windows, and watch where they stay together vs where they split.
That's exactly what `state.bandTraceFishSet` (turn 161) automates
for one cohort at a time, but it operates at L2-envelope
granularity and has no UI in the L3 panel itself.

This spec proposes **per-band rows** in the L3 contingency panel:
for each of the K bands at each L3 pane, a one-row summary of how
that single cohort behaves into the next pane. The K×K table
remains; the per-band rows are an additional surface alongside it.

The same compute, run across all adjacent panes in the L3 strip,
yields **band tracks**: a band track is a cohort that stays
together across consecutive panes, even when its raw K-means
label changes window-to-window. That's the automated version of
"pick a band, follow the fish."

### 1.1 Why this isn't already covered

- **Existing K×K tables** show all cohorts at once, summarized to
  one Cramér's V. The per-band rows decompose that summary.
- **`state.bandTraceFishSet` (turn 161)** tracks one user-picked
  cohort at L2 granularity. It doesn't enumerate bands; it
  doesn't run automatically across all bands at once; and it's
  invisible in the L3 panel itself (lives in the lines panel).
- **`SPEC_distant_band_concordance_fish_trajectory.md`** operates
  on **fish paths** (one row per fish, one column per L2). This
  spec operates on **cohort paths** (one row per band cohort, one
  column per pane). They're the dual representations of the same
  partition matrix; both are useful for different questions.

The **prior framing** that makes this clean: a band cohort is a
distribution prior over the fish-set ("samples currently in this
band"). At each adjacent window/L2, that prior gets pushed
through the K×K transition matrix. The single-band contingency
row is the prior's posterior — what fraction of the cohort
landed in each of the K destination bands.

| Mode of the posterior | What it means |
|---|---|
| Sharply peaked on one band | cohort persists |
| Spread across two bands | cohort splits / unstable assignment |
| Diffuse | cohort dissolves / loses coherence |
| Re-concentrates downstream after a diffuse window | A-B-A / double-crossover candidate |

### 1.2 Relationship to the existing band-trace work

Turns 160–164 shipped a **single-cohort** band-trace pipeline at
**L2-envelope** granularity:
- `state.bandTraceFishSet` — the cohort being tracked
- `_bandTraceForFishSet` — projects that cohort through every L2
  on the chromosome
- `_bandTraceRegimeRuns` — finds runs where the cohort stays
  cleanly in one band

This spec is the **multi-cohort** version at **L3-pane** granularity:
- Every band of every L3 pane gets a row
- Adjacent panes' bands get connected when the prior persists
- Connected components of those edges are the band tracks
- A-B-A detection looks for cohorts that re-concentrate after a
  diffuse gap

The two views compose. The band-trace pipeline is the
single-cohort interactive ("show me this one cohort across the
whole chromosome"). This spec is the multi-cohort enumerative
("automatically find every coherent cohort within this candidate's
L3 strip").

### 1.3 Where this lives

**Atlas-side**, not HPC. All inputs already exist in the browser:
- `state.candidate.l2_indices` — the panes in the L3 strip
- `getL2Cluster(l2idx).fixedKLabels` — per-sample K-means
  assignments per pane (the K=3 default, or K=6 if the user has
  switched K)
- `_buildContingency(labelsA, labelsB, KA, KB)` — already exists
  (line 12250)
- `_ssContingencyTableHtml(ct)` — already exists (line 12654);
  the new per-band rows render alongside or below this
- `state._btraceHits`, `_bandTraceForFishSet` — primitives we can
  reuse for the cohort projection

The output JSON layer (band_tracks.json, see §6) is optional —
the in-session compute is enough for the L3 panel UI. JSON
emission becomes useful only when Quentin wants to ship band
tracks as supplementary data alongside the manuscript.

### 1.4 The het band is the principal track (manuscript backbone)

Of all the band tracks the algorithm extracts, **the het-carrier
track is biologically privileged**. It is the inversion's
**skeleton**: the chain of windows where the same heterozygous
carrier set persists. Heterozygous carriers — fish carrying one
ancestral and one derived arrangement — are the operational
signature of the inversion. Where they remain coherent, the
inversion's signal is detectable; where they dissolve, the
detectable signal is gone. The het track is what the manuscript
will cite as "this candidate spans Mb X to Mb Y."

This means the spec runs in **two modes**:

1. **Symmetric mode** — every band gets a track (the §5 path
   construction operates on all K bands equally). Useful for
   complex regions, nested patterns, and the per-band L3 rows
   (§7.1).

2. **Het-anchored mode** — the het band is identified per pane
   first; its carrier track is extracted as the **principal
   skeleton**; other bands' tracks are computed but reported as
   **subsidiary** to the skeleton. This is the manuscript-grade
   view.

Het-anchored mode answers: *"Where does the heterozygous carrier
set first form a coherent cohort, where does it stay coherent,
where does it dissolve, and (when the same carriers reappear
downstream) where is the same skeleton continuing across an
internal interruption?"*

**Critical wording discipline.** A loss of het signal in a
window does **not** mean the inversion ends there. Possible
reasons for a het-signal dropout that don't imply an interval
boundary:
- recombination or double crossover
- gene conversion
- marker dropout / locally low informative-SNP density
- shared ancestral haplotype in this window
- ancestry mosaicism
- noisy K-means band assignment
- local low diversity

The skeleton-extraction algorithm therefore distinguishes
**dropout** (signal stops, then resumes with the same carriers)
from **boundary** (signal stops, and downstream carriers are
genuinely different). The §5.3 split-track logic generalises to
this case — a skeleton with an internal dropout is one
multi-segment skeleton, classified `het_skeleton_with_dropout`,
not two separate skeletons. Operational rule: the algorithm
flags candidate dropouts; the user (or downstream test) decides.
Vocabulary used in outputs:
- *"detectable divergent haplotype signal stops here"* (correct)
- *"inversion ends here"* (NOT used by the algorithm; only Quentin
  can write that in the manuscript after reviewing)

### 1.5 The lines panel as the band-track view

The atlas's per-sample-lines panel currently colours each fish's
line by its raw K-means group at every pane. K-means group ids
rotate window-to-window (Hungarian alignment fixes this for
display in the L3 K×K table, but the lines panel does its own
coloring), so a fish that stays in "the same band" biologically
gets a different color in every pane. The lines panel is hard
to read for exactly this reason.

**Band-track membership gives every fish a stable color.** A
fish in the het skeleton from pane 0 to pane 8 paints one solid
color across the whole interval. A fish that drops out of the
skeleton at pane 4 paints the skeleton color in panes 0–3, then
something else (its actual band-track) in panes 4–8. A fish in
a split-track that disappears at panes 4–5 and reappears at
pane 6 paints the same color in panes 0–3 and panes 6–8 with a
neutral grey at 4–5.

The lines panel becomes the band-track display:
- A coherent inversion = a thick band of synchronized colors
  spanning the candidate's interval
- A complex region = many short stripes, fish lines that change
  color rapidly
- A double-crossover-like pattern = a fish line that returns to
  the skeleton color after a brief detour
- An internal dropout = a fish line that goes neutral grey
  briefly, then returns to the skeleton color

Implementation: this is a new value for the existing
`state.linesColorMode` slot — `'band_track'` (or, more
specifically, `'band_track_het'` when in het-anchored mode).
The resolver `_resolveSampleColorByMode(sample_idx, window)`
looks up the sample's track membership at the window from
`state.bandTracks[cand.id]` and returns the track's color. When
the sample isn't in any track at the window, returns the
"neutral grey" color (`var(--ink-dimmer)` at low alpha).

The existing per-sample-lines color-mode picker scaffolding from
v3.99 turn 14e+ already exists; only the resolver body needs
filling in.

---

## 2. Vocabulary discipline

This spec touches manuscript-grade claims, so vocabulary matters:

- **band cohort** — the set of samples assigned to one K-means
  group at one window / L2 pane. Per-window object. NEVER use
  this to mean a multi-window structural object.
- **band track** — a cohort that stays together across consecutive
  panes, possibly with raw K-means labels changing pane-to-pane.
  Multi-pane object. The output of this spec.
- **track segment** — a contiguous run of panes a band track
  occupies. A track may have one segment (continuous) or multiple
  (split / A-B-A).
- **persists** — the operational definition: cohort A at pane i
  and cohort B at pane i+1 are connected if the connection rule
  (§4) is satisfied. Implies neither biology nor inversion
  ancestry — purely an empirical sample-set overlap.
- **A-B-A** / **split track** / **reappearing track** — three
  names for the same thing: cohort C appears in panes i..j, is
  absent in panes j+1..k-1, and reappears with high core-sample
  overlap in panes k..end. Use **split track** in code, **A-B-A**
  in figure captions / docs.
- **het band** — at one pane, the K-means band whose members
  carry the heterozygous / divergent haplotype signal. Identified
  per §4A.1. May be one band, two bands, or zero bands per pane;
  never forced to exactly one when the data is ambiguous.
- **het skeleton** / **het-carrier skeleton** / **same-carrier
  backbone** — the band track formed by chaining het bands across
  panes via the connection rule. The principal track for the
  candidate. May be continuous, may have internal dropouts, may
  be split / reappearing.
- **het carrier** — a sample that is currently a member of a het
  band. Per-pane membership; a sample can be a carrier in some
  panes and not others.
- **carrier set** — the sample-set membership of a het band at
  one pane. The §4A continuity check operates on carrier sets.
- **dropout** — a pane where the het signal is absent / weak but
  flanking panes have detectable het carriers with high overlap.
  Operationally: a gap inside an otherwise-coherent skeleton.
  Distinct from a **boundary**, where downstream carriers are
  genuinely different.
- **detectable divergent haplotype signal** — preferred phrasing
  for what the het signal measures. Avoids implying that signal
  loss = inversion absence.
- **"het signal stops"** vs **"inversion stops"** — the algorithm
  detects the former (an operational signal property); the latter
  is a biological claim. The spec, the UI tooltips, and the TSV
  outputs all use the operational phrasing only. The biological
  interpretation is reserved for the manuscript text.
- **"dropout does not necessarily mean boundary"** — canonical
  slogan for the §4A.3 internal-dropout-tolerance behaviour.
  Use in tooltips when rendering a `skeleton_with_dropout`, and
  in the manuscript methods section. Complementary slogan:
  **"different downstream carriers = different system."**
- **carrier stability score** — fraction of non-dropout panes in
  which every consensus carrier was a member of the het band.
  Stored on the het skeleton (§4A.4); reported in
  `system_band_summary.tsv`. Range [0, 1]; 1.0 = the consensus
  carriers are het in every non-dropout pane.

Avoid:
- "this band is HOM_REF / HET / HOM_INV" — band tracks have no
  karyotype labels in this spec; that's the L-label vocabulary's
  job in the H-label classification spec.
- "this is a double crossover" — the algorithm can flag cases
  that *look like* double crossovers (A-B-A pattern with high
  core-sample overlap), but biological double crossover is a
  claim that needs separate validation. Use `class:
  split_track_reappearing` in outputs and **possible_interpretation:
  double_crossover_like** as a separate field.
- "track A is the same as track B" — even with high core-sample
  overlap across a gap, the algorithm reports a split-track
  candidate; only Quentin's review (or a downstream test) decides
  whether they're the same biological cohort.
- **"the inversion ends here"** — operational outputs only say
  "het signal stops here" or "carrier-set continuity breaks
  here." The biological interpretation is reserved for the
  manuscript text, written by the author.

---

## 3. The data flow

```
state.candidate.l2_indices            // panes in the L3 strip
        ↓
getL2Cluster(l2idx).fixedKLabels      // per-pane sample → group
        ↓
S[pane_i, g] = Set<sample_idx>         // band cohorts (§3.1)
        ↓
edges[(i, g_i, i+1, g_j)] = {          // adjacent-pane connections (§4)
  overlap, jaccard, retention_forward, retention_reverse,
  connected: bool
}
        ↓
band_tracks = build_paths_from_edges()  // connected components (§5)
        ↓
classify_tracks() + detect_split_tracks() // multi-segment merger (§5.3)
        ↓
classify_system()                       // overall L3-strip class (§5.4)
        ↓
UI: per-band rows in the L3 panel       // primary deliverable (§7)
optional: band_tracks.json / TSV exports // supplementary data (§6)
```

### 3.1 Band cohorts

For each L3 pane (= L2 in `state.candidate.l2_indices`):

```js
// At L2 = panes[i], with K = state.k (3 default, 6 optional)
const cluster = getL2Cluster(panes[i]);
const labels = cluster.fixedKLabels;   // Int8Array[n_samples], values 0..K-1 or -1
// Build sample-set per group
const S_i = new Array(K);
for (let g = 0; g < K; g++) S_i[g] = new Set();
for (let s = 0; s < labels.length; s++) {
  const g = labels[s];
  if (g >= 0 && g < K) S_i[g].add(s);
}
```

Edge cases:
- `label === -1` → sample is uncalled at this pane → not in any
  cohort, doesn't participate in any edge for this pane.
- Cohort with `|S| < min_band_size` (default 5) → still a cohort,
  but flagged with `unstable: true` because Jaccard ratios get
  noisy on tiny groups. Edges from/to unstable cohorts are
  computed but the connection rule's threshold is relaxed.

Reuses the same min-band-size threshold as `_AUTO_PROMOTE_MIN_BAND_SIZE`
from the L2-sweep spec (5 fish) for consistency.

### 3.2 Input columns (HPC-conceptual ↔ atlas-primitive)

Quentin's request enumerated the per-window inputs for the R-module
form of this algorithm:

> *"chromosome, start, end, sample_id, band_label, optional PC1/PC2,
> optional dosage summary, optional heterozygosity count/fraction,
> optional missingness, optional local ancestry Q metrics"*

The atlas-side equivalents already exist as in-browser primitives —
the table below is the canonical mapping. New compute is needed only
for the local-ancestry Q column, which is gated on a JSON layer
that hasn't shipped yet.

| HPC column | Atlas primitive | Status |
|---|---|---|
| `chromosome` | `state.data.chrom` | shipped (every JSON load) |
| `start`, `end` | `state.data.windows[w].start_bp`, `.end_bp` | shipped |
| `sample_id` | `state.sample_index[s]` | shipped |
| `band_label` | `getL2Cluster(l2idx).fixedKLabels[s]` | shipped |
| `PC1`, `PC2` | `getL2Cluster(l2idx).pca.U[s][0]`, `.U[s][1]` | shipped |
| `dosage summary` | `window.popgenDosage.getCachedChunk(...)` | shipped (turn 14e) |
| `heterozygosity count/fraction` | `_computeHetRateForL2(l2idx)[s]` | shipped (turn 129) |
| `missingness` | `fixedKLabels[s] === -1` (uncalled at this pane) | shipped |
| `local ancestry Q` | new JSON layer `local_ancestry_q` | **not shipped**, optional input |

The spec does not require the local-ancestry Q layer for the
het-skeleton compute. When present, it serves as an additional
het-enrichment-score signal (§4A.1) — bands with admixed-Q
profiles can be flagged as `het_candidate` even when raw
het-rate is moderate. When absent, the algorithm falls back to
the het-rate / dosage / PC1-middle hierarchy without warning.

The HPC column is enumerated here so a future R-module mirror
of this spec (per the trigger's *"design this as a robust R
module that can plug into my existing inversion pipeline"*)
has a stable input contract. The atlas implementation and the
HPC implementation should produce comparable outputs from the
same input set; differences should only arise from threshold
choices documented in the TSV header rows.

---

## 4. The connection rule (per Quentin's design)

For every adjacent pane pair (i, i+1) and every (g_i, g_j) pair:

```js
const overlap = sizeIntersection(S[i][g_i], S[i+1][g_j]);
const union   = sizeUnion(S[i][g_i], S[i+1][g_j]);
const jaccard = overlap / union;
const rf      = overlap / S[i][g_i].size;     // retention forward
const rr      = overlap / S[i+1][g_j].size;   // retention reverse
```

**Connection predicate** (Quentin's recommendation; thresholds
configurable):

```js
const connected =
  (jaccard >= JACCARD_MIN) ||
  (rf >= RETENTION_MIN && rr >= RETENTION_MIN);
```

Defaults:
- `JACCARD_MIN = 0.65`
- `RETENTION_MIN = 0.75`

Why both predicates: Jaccard penalises group-size changes
symmetrically (a cohort growing from 30 to 50 fish drops Jaccard
even when retention is 100%). Retention pair handles asymmetric
changes correctly: cohort A (n=30) entirely contained in cohort B
(n=50) gets `rf=1.0, rr=0.6, jaccard=0.6` — Jaccard says no, but
retention-pair says yes (rf > 0.75 but rr < 0.75 → no), so the
**OR** with Jaccard catches the case where Jaccard is acceptable
even though retention isn't perfectly symmetric.

Edge case: empty cohort. If `|S[i][g_i]| === 0`, all edges from
that node have `rf` undefined. Treat as `connected = false` and
`jaccard = NaN`. The cohort node still exists (it's just empty)
but has no outgoing edges.

### 4.1 Why not Hungarian alone

Hungarian projection (already shipped in the atlas via
`_hungarianChainProjection`) gives one global label-permutation
per pane pair. That's the right answer for *system-level* merging
("what's the best way to relabel pane i+1 so its labels match pane
i?"). It picks one permutation, locks it in, and computes Cramér's V
on the diagonal-aligned table.

For per-band tracking, Hungarian is the wrong tool because:
- It assumes a 1-to-1 cohort mapping pane-to-pane. Real data has
  splits (cohort A at pane i → two cohorts at pane i+1) and
  merges (two at pane i → one at pane i+1) that Hungarian can't
  represent.
- It commits to the best total assignment, not the best per-band
  assignment. A row's "true" successor might be the second-best
  Hungarian match.

Connection edges are 1-to-many and many-to-1 by construction —
exactly what a split or merge needs.

The two are complementary: Hungarian for system-level pane
alignment, connection edges for per-band tracking. Outputs
should mention which Hungarian permutation was applied at each
pane boundary so the user can cross-reference.

---

## 4A. Het-anchored mode (the principal-track variant)

The §4 connection rule is symmetric — every band gets equal
treatment. **Het-anchored mode** breaks that symmetry by
identifying the het band at each pane *first*, then extracting
its carrier track as the principal skeleton. This is the mode
that produces the manuscript-grade interval claim.

### 4A.1 Identifying the het band(s) per pane

For each pane i with K bands, compute a **het-enrichment score**
per band, then classify each band as:

- `het` — high het-enrichment, large enough to be a real cohort
- `het_candidate` — moderate het-enrichment, secondary to the
  primary het band but tracked as a parallel skeleton
- `non_het` — low het-enrichment (HOM_REF or HOM_INV)

**Het-enrichment score per band**, computed in priority order
of available signals:

1. **Mean per-sample het rate** (primary, when `dosage_chunks`
   is present): leverages the existing `_computeHetRateForL2(l2idx)`
   primitive (turn 129) which returns `Float32Array[n_samples]`
   of het rates ∈ [0, 1]. Per band, take the mean over members:
   ```js
   const hetRate = _computeHetRateForL2(l2idx);
   const bandMembers = S[i][g];
   const meanHet = mean(hetRate[s] for s in bandMembers);
   ```
   Bands with `meanHet >= 0.35` are het-enriched (the
   diallelic-locus expected-HET-rate is 0.5, but real bands
   with mixed populations come in around 0.35–0.45).

2. **Mean dosage** (fallback when `dosage_chunks` is present
   but het rates are degenerate): per band, mean dosage value
   ∈ [0, 2]. Bands with `meanDosage` in `[0.7, 1.3]` are
   het-enriched.

3. **Middle-PC1 heuristic** (fallback when no dosage layer): at
   K=3, the middle PC1 band is canonically the het band. At
   K=6 with no dosage, no per-band het score is available — the
   algorithm falls back to symmetric mode and emits a warning
   in `state.bandTracks[cand.id].warnings`.

The het-enrichment score is normalised to [0, 1] regardless of
which signal source produced it (the middle-PC1 heuristic emits
1.0 for the middle band, 0.0 otherwise). Stored per band:

```js
{
  pane_idx, group_id,
  het_enrichment: 0.42,         // [0, 1]
  het_class: 'het' | 'het_candidate' | 'non_het',
  het_signal_source: 'het_rate' | 'dosage' | 'pc1_middle',
}
```

**Multiple het bands per pane.** At K=6 a window can have two
plausibly-het bands (for example, two slightly different het
arrangements with different ancestral allele frequencies).
Don't force a unique het band: classify as `het` any band with
`het_enrichment >= 0.55`, and `het_candidate` any band with
`0.35 <= het_enrichment < 0.55`. Each one seeds a parallel
skeleton; the persistence step (§4A.2) decides which seeds
survive into named skeletons by checking carrier-continuity.

**No het band identified.** If no band passes the het threshold
at a pane, mark the pane as a **dropout candidate** for the
skeleton: it's a pane where the het signal is absent, which is
the §1.4 ambiguity ("dropout vs boundary"). The skeleton-build
step (§4A.3) decides whether the dropout is internal (resumes
later with same carriers) or terminal (genuinely a boundary).

### 4A.2 Carrier-continuity edges

Once het bands are identified per pane, build edges using the
same connection rule from §4 — but only between het bands at
adjacent panes (or panes within `MAX_DROPOUT_GAP`, see §4A.3).

The metrics are exactly Quentin's wording:

```
J(A, B) = |A ∩ B| / |A ∪ B|                    // Jaccard
O(A, B) = |A ∩ B| / min(|A|, |B|)              // overlap coefficient
R(A→B) = |A ∩ B| / |A|                          // forward retention
R(B→A) = |A ∩ B| / |B|                          // reverse retention
```

**Connection predicate** for het skeletons:

```js
const connected =
  (jaccard >= HET_JACCARD_MIN) ||
  (overlap_coef >= HET_OVERLAP_MIN) ||
  (rf >= HET_RETENTION_MIN && rr >= HET_RETENTION_MIN);
```

The overlap coefficient is added to the §4 predicate because it
handles the **subset/dropout** case better than Jaccard alone:
a het band of 30 carriers shrinking to a het band of 18 carriers
that are all subsets of the previous 30 produces
`jaccard = 18/30 = 0.6, overlap = 18/18 = 1.0, rf = 0.6, rr = 1.0`.
Jaccard says marginal; overlap says clean subset.

Defaults (Quentin's recommended values):
- `HET_JACCARD_MIN = 0.5` (looser than the symmetric §4 default
  of 0.65 — het skeletons can grow/shrink across panes due to
  marker dropout, so we accept lower Jaccard when overlap is
  high)
- `HET_OVERLAP_MIN = 0.75`
- `HET_RETENTION_MIN = 0.75`

### 4A.3 Skeleton extraction with internal-dropout tolerance

The §5.1 path-construction algorithm runs over the het edges,
but with one modification: the forward walk can **skip up to
MAX_DROPOUT_GAP panes** when the next het band is missing or
disconnected, *as long as* a downstream het band re-establishes
high carrier-overlap with the current band.

```js
const MAX_DROPOUT_GAP = 3;       // panes (configurable)
const DROPOUT_OVERLAP_MIN = 0.7; // jaccard between A and C across the gap

// At pane i with het band g_i = current skeleton tip:
//   If pane i+1 has a connected het band → extend normally.
//   Else, look ahead up to MAX_DROPOUT_GAP panes for a het band
//   whose carrier set has Jaccard >= DROPOUT_OVERLAP_MIN with g_i.
//   If found at pane i+k → skip panes i+1 .. i+k-1, mark them
//   as internal dropouts of the skeleton, extend to i+k.
//   Else → terminate the skeleton at pane i.
```

This is the §1.4 dropout-vs-boundary distinction made
operational. Internal dropouts get tagged on the skeleton:

```js
{
  skeleton_id: 'HET_SK_<chrom>_<system>_01',
  segments: [
    { start_pane: 0, end_pane: 4 },
    { start_pane: 6, end_pane: 9 },          // resumed after dropout
  ],
  dropout_panes: [5],                         // the gap
  dropout_diagnosis: {
    n_dropout_panes: 1,
    pane_5_n_carriers: 0,                    // no het band identified
    cross_gap_jaccard: 0.86,                 // carriers at pane 4 vs pane 6
    classification: 'internal_dropout',     // or 'boundary' if jaccard fails
  },
}
```

When the cross-gap Jaccard falls below `DROPOUT_OVERLAP_MIN`,
the algorithm does NOT bridge the gap — the two segments stay
as separate skeletons (and may form an A-B-A pair via §5.3 if
the §5.3 looser thresholds catch them).

### 4A.4 Distinguished outputs for het skeletons

Each candidate's `state.bandTracks[cand.id]` gains a
`het_skeleton` field alongside the symmetric `tracks` array:

```js
state.bandTracks[cand.id] = {
  // ... existing fields from §6.1 ...
  tracks: BandTrack[],            // symmetric mode (every band)
  het_skeleton: HetSkeleton | null,  // principal track in het-anchored mode
  het_skeleton_secondary: HetSkeleton[],  // additional het_candidate skeletons
  het_signal_source: 'het_rate' | 'dosage' | 'pc1_middle',
  warnings: string[],             // e.g. 'no_dosage_layer_falling_back_to_pc1'
}
```

`HetSkeleton` shape extends `BandTrack` (§5.2) with:
- `principal: true` (vs `principal: false` for secondary skeletons)
- `dropout_panes: number[]`
- `dropout_diagnosis: {...}`
- `n_carriers_per_pane: number[]` (length = n_panes including dropouts; 0 at dropout panes)
- `consensus_carriers: Set<sample_idx>` (intersection across all non-dropout panes)
- `consensus_carrier_count: number`
- `carrier_stability_score: number` ∈ [0, 1] — fraction of panes
  where every consensus_carrier was a member of the het band
- `het_enrichment_per_pane: number[]`
- `mean_jaccard`, `mean_overlap`, `mean_retention_forward`, `mean_retention_reverse`
- `class: 'clean_skeleton' | 'skeleton_with_dropout' | 'reappearing_skeleton' | 'mosaic' | 'fragmented'`
- `class_reason: string`

Skeleton-level classification rules:

```js
if (n_dropout_panes === 0 && n_segments === 1 && stability >= 0.9)
  class = 'clean_skeleton';
else if (n_dropout_panes >= 1 && n_segments === 1 && stability >= 0.8)
  class = 'skeleton_with_dropout';
else if (n_segments >= 2 && cross_segment_jaccard >= 0.7)
  class = 'reappearing_skeleton';
else if (stability < 0.6)
  class = 'mosaic';
else
  class = 'fragmented';
```

### 4A.5 Why this lives inside the symmetric-mode spec, not separately

The het-skeleton extraction reuses every primitive from §3–§5:
- the cohort-set construction (§3.1) is identical
- the connection rule (§4) generalises with one extra metric
  (overlap coefficient) and looser thresholds
- the path construction (§5.1) generalises with one extra step
  (dropout-skip lookahead)
- the §5.3 split-track detection runs as-is on the resulting
  segments
- the §5.4 system classification gets the het skeleton's class
  as its dominant input when in het-anchored mode

So het-anchored mode is a **specialisation**, not a separate
algorithm. The implementation is one extra `mode: 'symmetric' |
'het_anchored'` parameter on `_buildBandTracksForCandidate`,
plus the het-band identifier (§4A.1) and the skeleton-output
post-processing (§4A.4).

### 4A.6 Dropout causes (why "het signal stops" ≠ "inversion stops")

The §4A.3 dropout-tolerance logic exists because **het signal
loss has many possible biological and technical causes**, only
one of which is "the inversion ended." Operationally, the
algorithm cannot distinguish these — but the spec records the
causes here so the reported `class: skeleton_with_dropout` is
read as "the het signal weakened but the carriers re-establish
downstream", not as a positive claim about a specific mechanism.

Possible causes of an internal dropout (Quentin's enumeration,
verbatim from the trigger):

- **recombination or double crossover** — a recombination event
  inside the inversion (rare for true inversions but possible
  near boundaries) breaks the linkage that produced the het
  signal in this window
- **gene conversion** — a small-scale gene-conversion tract
  homogenises a stretch of the divergent haplotype, locally
  erasing the heterozygosity that the K-means classifier reads
- **marker dropout** — the SNP set in this window is
  unfortunate (mostly invariant within the carriers), so the
  PCA / clustering can't separate carriers from non-carriers
  even though the underlying haplotype structure is intact
- **local low diversity** — the window happens to be in a
  conserved region, so all samples look similar regardless of
  inversion state
- **shared ancestral haplotype** — the inversion's two
  arrangements share a common ancestor stretch in this window;
  carriers and non-carriers genuinely have similar dosage here
- **ancestry mosaicism** — admixture or introgression has
  imported a tract from a different population, locally
  scrambling the cluster identities
- **noisy K-means band assignment** — at K=6, K-means
  occasionally lands in a local optimum that splits the het
  band incorrectly; reseeding would recover the signal but the
  algorithm doesn't get to retry
- **local lack of informative SNPs** — too few segregating sites
  in this window to support the clustering at the requested K

The first three are biological, the last four are technical or
statistical. The algorithm cannot tell them apart. The spec's
operational rule:

- **"het signal stops"** = the algorithm cannot detect divergent
  haplotype signal in this window
- **"inversion stops"** = a *biological* claim about where the
  inverted segment ends in the genome

The skeleton is a het-signal object; it answers "how far does
the detectable divergent haplotype signal extend?" The
inversion's biological boundary may extend further (when the
algorithm misses signal due to one of the technical causes)
or shorter (when carriers genuinely share haplotypes inside the
inversion through gene conversion or shared ancestry). The
manuscript text should report the het-skeleton interval as the
**operational interval** and discuss the biological-interval
caveat separately.

This is why §6.4's `sample_band_paths.tsv` consensus_band rule
requires both `dominant_fraction >= 0.85` AND `num_switches <=
1` — without that guard, samples in the dropout panes would get
force-assigned to the skeleton's principal band even though
their actual K-means assignment in the dropout pane was
ambiguous. The guard preserves the distinction between "this
sample is consistently in the het skeleton" (clean) and "this
sample is mostly in the het skeleton but had ambiguous panes"
(uncertain — needs case-by-case review).

### 4A.7 The dropout-vs-boundary slogan

For figure captions and manuscript text, the canonical phrasing
is:

> "Loss of detectable divergent haplotype signal does not
> automatically indicate a structural boundary. The het-carrier
> skeleton bridges internal dropouts when carrier identity is
> recovered downstream; only when downstream carriers are
> genuinely different does the algorithm report a boundary."

Short slogan: **"dropout does not necessarily mean boundary."**
Use this phrasing in tooltips, in the L3 panel UI when a
skeleton-with-dropout is rendered, and in the manuscript's
methods section describing how the inversion intervals were
called.

The complementary slogan for the boundary case:
**"different downstream carriers = different system."**
Use when the carrier set after a gap has low Jaccard with the
carrier set before the gap (and the two segments stay separate
skeletons rather than merging).

---

## 5. Building band tracks from edges

### 5.1 Path construction

A band track is a connected component of the (band-cohort node)
× (connection edge) graph, restricted to consecutive panes.

```
nodes:    (pane_idx, group_id) for every (pane, K) cell
edges:    connection edges from §4 between adjacent panes only
```

Naive connected-components on this graph would lump all of pane i
and pane i+1 together if any one cohort straddles them. That's
wrong — band tracks are paths, not arbitrary CC's. Algorithm:

```
For each starting pane (left-to-right):
  For each cohort node not yet assigned to a track:
    Greedy forward walk:
      From current node, find the best outgoing edge
      (= max-overlap connected edge to the next pane).
      If found, append the destination node to the track.
      If multiple equally-good edges, mark the cohort as splitting
      and start two tracks (one continuing on each branch).
      Stop when the next pane has no connected destination.
    Greedy backward walk (symmetric) from the start node.
    The result is one band track.
```

Implementation note: implement the forward-walk as the canonical
direction; ensure determinism by tie-breaking on (a) destination
group_id ascending, (b) destination cohort size descending.

Splits naturally produce two tracks sharing a prefix. To avoid
double-counting cohorts:
- Track parent_track_id when a split occurs.
- A cohort's "membership" is the latest-walked track that
  includes it; earlier tracks see the cohort as a branch point
  with `splits_into: [track_id1, track_id2]`.

Merges (two tracks → one cohort at pane k) are recorded similarly:
the merged track stores `merged_from: [track_id1, track_id2]`.

### 5.2 Continuous track summary

For each track segment (contiguous run of panes the track
occupies):

```js
{
  track_id: 'BT_<chrom>_<system>_<seq>',
  segment_id: 0,                          // 0 if single-segment
  start_pane: 0, end_pane: 4,
  start_l2: 142, end_l2: 156,
  start_bp: 15_100_000, end_bp: 15_800_000,
  n_panes: 5,
  // edge stats (averaged over the n_panes - 1 connections)
  mean_jaccard, mean_retention_forward, mean_retention_reverse,
  // cohort-size stats (averaged over the n_panes cohorts)
  mean_group_size, median_group_size,
  // sample composition
  core_samples: Set<sample_idx>,            // intersection across all panes
  variable_samples: Set<sample_idx>,        // (union − intersection)
  n_core, n_variable,
  // per-pane raw K-means group ids it actually used
  pane_group_path: [g_0, g_1, g_2, ...],    // length = n_panes
  // bookkeeping
  splits_into: [track_id, ...] | null,
  merged_from: [track_id, ...] | null,
}
```

`core_samples` is the foundational object: it's the set of fish
that are members of every cohort in the track. `variable_samples`
is the set that's a member of at least one but not all. Quentin's
manuscript-grade output should report `n_core` and `n_variable`
together so readers can see how stable the cohort is.

### 5.3 Split-track / A-B-A detection

After all continuous tracks are built, scan for pairs where:

```
track_T1: panes [a..b], core_samples C1
track_T2: panes [c..d], core_samples C2

c > b + 1                                     // gap of ≥ 1 pane
genomic_distance(b, c) <= MAX_GAP_BP          // within reasonable distance
|C1 ∩ C2| / |C1 ∪ C2| >= GAP_OVERLAP_MIN      // cores overlap heavily
|C1 ∩ C2| >= MIN_CORE_OVERLAP_N               // absolute floor
```

When matched, emit a multi-segment band track:

```js
{
  track_id: 'BT_<chrom>_<system>_<seq>',
  class: 'split_track_reappearing',
  possible_interpretation: 'double_crossover_like',
  segments: [
    { segment_id: 0, ...summary of T1... },
    { segment_id: 1, ...summary of T2... },
  ],
  gap_bp: c.start_bp - b.end_bp,
  core_sample_overlap_between_segments: jaccard(C1, C2),
  // What was happening in the gap?
  gap_panes: [b+1 .. c-1],
  gap_diagnosis: {
    // Did the prior go diffuse, or did it concentrate on a
    // different cohort?
    concentrated_on_other: bool,            // the cores show up in another track
    other_track_id: track_id | null,
    // Or did the prior just dissolve (no high-mass cohort in the gap)?
    diffuse: bool,
  },
}
```

Defaults:
- `MAX_GAP_BP = 2_000_000` (2 Mb — within a typical inversion span)
- `GAP_OVERLAP_MIN = 0.7` (Jaccard between cores)
- `MIN_CORE_OVERLAP_N = 8` (don't merge tiny tracks just because
  3 fish happen to recur)

Critical guard: **never auto-merge two tracks across a gap if the
gap_diagnosis shows the cohort was concentrated on a *different*
cohort in the middle**. That's the "the same fish briefly joined
another band, then came back" pattern, which is exactly the
double-crossover signal — but it's the user's call whether to
treat that as A-B-A or as two separate tracks. The algorithm
flags it; the UI presents it; the user decides.

### 5.4 System classification

After all tracks (continuous + multi-segment) are catalogued for
the L3 strip, classify the whole strip:

```js
const class_buckets = {
  full_span:        tracks where (n_panes covered by track) >= 0.85 * total_panes,
  internal:         tracks fully inside the strip, not at edges, n_panes >= 3,
  short:            tracks with n_panes < 3,
  split:            multi-segment tracks (class === 'split_track_reappearing'),
};

let classification, reason;
if (class_buckets.full_span >= 2 && class_buckets.split === 0
    && class_buckets.internal === 0) {
  classification = 'simple_single_inversion';
  reason = `${class_buckets.full_span} full-span tracks, no internal/split`;
} else if (class_buckets.full_span >= 1 && class_buckets.internal >= 1) {
  classification = 'nested_or_internal_component';
  reason = `${class_buckets.full_span} full-span + ${class_buckets.internal} internal`;
} else if (class_buckets.split >= 1) {
  classification = 'double_crossover_like';
  reason = `${class_buckets.split} split track(s) detected`;
} else if (class_buckets.short >= total_tracks * 0.5) {
  classification = 'fragmented_or_noisy';
  reason = `${class_buckets.short} of ${total_tracks} tracks shorter than 3 panes`;
} else {
  classification = 'complex';
  reason = 'mixed full-span, internal, and short tracks';
}
```

`adjacent_systems` from Quentin's request requires multi-system
context (tracks side-by-side spanning *different* L3 strips); it's
out of scope for this single-strip classifier and would be added
when the spec extends to whole-chromosome compute.

---

## 6. Output objects (in-session + optional JSON/TSV)

### 6.1 In-session state

```js
state.bandTracks: {
  [systemKey]: {                             // systemKey = candidate.id
    chrom: string,
    panes: number[],                         // L2 indices in this strip
    mode: 'symmetric' | 'het_anchored',      // which compute mode produced this
    tracks: BandTrack[],                     // §5.2 / §5.3 shape (always present)
    // Het-anchored mode adds the following fields (§4A):
    het_skeleton: HetSkeleton | null,        // principal track when mode='het_anchored'
    het_skeleton_secondary: HetSkeleton[],   // additional `het_candidate` skeletons
    het_signal_source: 'het_rate' | 'dosage' | 'pc1_middle' | null,
    warnings: string[],                       // e.g. 'no_dosage_layer_falling_back_to_pc1'
    // System-level
    classification: string,
    classification_reason: string,
    computed_at_K: number,                   // K=3 default, K=6 if Quentin switched
    cache_key: string,                       // invalidates on candidate / K / mode change
  }
}
```

Compute is per-candidate (one L3 strip at a time). Cache key
includes the candidate id, panes, K, the connection-rule
thresholds, and the mode (symmetric vs het_anchored). Switching
mode invalidates the cached result.

Het-anchored mode is the default when `dosage_chunks` is
present; symmetric mode is the fallback when no het-enrichment
score can be computed (and is also user-toggleable via the
mode chip in §7.3).

### 6.2 Output files (optional, supplementary-table use)

Per Quentin's request — four TSVs. Each row identified by
`(system_id, band_track_id, segment_id)` triple where applicable.

**A. `band_tracks.tsv`** — one row per (track, segment):

| column | type | source |
|---|---|---|
| `system_id` | str | `state.candidate.id` |
| `band_track_id` | str | `BT_<chrom>_<system>_<seq>` |
| `chrom` | str | `state.candidate.chrom` |
| `track_class` | str | `simple` / `internal` / `short` / `split_track_reappearing` |
| `segment_id` | int | 0 for single-segment; 0,1,... for multi-segment |
| `start_pane` | int | local pane index in `state.candidate.l2_indices` |
| `end_pane` | int | inclusive |
| `start_l2` | int | absolute L2 index |
| `end_l2` | int | inclusive |
| `start_bp` | int | from L2 envelope |
| `end_bp` | int | from L2 envelope |
| `n_panes` | int | end_pane - start_pane + 1 |
| `mean_jaccard` | float | averaged over edges in this segment |
| `mean_retention_forward` | float | |
| `mean_retention_reverse` | float | |
| `mean_group_size` | float | averaged over panes in segment |
| `median_group_size` | int | |
| `n_core_samples` | int | |
| `core_sample_ids` | str | comma-joined sample ids; empty if `n_core_samples > 50` (use export TSV instead) |
| `pane_group_path` | str | comma-joined raw K-means group ids the track used |
| `splits_into` | str | comma-joined other track_ids when this track branches |
| `merged_from` | str | comma-joined other track_ids when this track is a merge |
| `gap_bp_to_next_segment` | int | for multi-segment tracks; NA otherwise |
| `core_sample_overlap_between_segments` | float | for multi-segment; NA otherwise |
| `notes` | str | free-text flags (e.g., `unstable_band` for tiny cohorts) |

**B. `band_track_edges.tsv`** — one row per adjacent-pane edge:

| column | type |
|---|---|
| `system_id` | str |
| `chrom` | str |
| `pane_i` | int |
| `pane_j` | int (= pane_i + 1) |
| `l2_i`, `l2_j` | int |
| `group_i`, `group_j` | int |
| `n_i`, `n_j` | int (cohort sizes) |
| `overlap` | int |
| `union_size` | int |
| `jaccard` | float |
| `retention_forward` | float |
| `retention_reverse` | float |
| `connected` | bool |
| `connection_reason` | str (`jaccard` / `retention` / `both` / `none`) |
| `track_id` | str | NA if edge not in any track |

This is the largest of the four files (K² × n_panes - 1 rows per
system) but is the most diagnostically useful — every connection
decision is auditable.

**C. `system_band_summary.tsv`** — one row per system (candidate):

| column | type |
|---|---|
| `system_id` | str |
| `chrom` | str |
| `system_start_bp`, `system_end_bp` | int |
| `n_panes` | int |
| `K` | int (the K used for compute) |
| `n_band_tracks` | int |
| `n_full_span_tracks` | int |
| `n_internal_tracks` | int |
| `n_split_tracks` | int |
| `n_short_tracks` | int |
| `classification` | str |
| `classification_reason` | str |
| `connection_thresholds` | str (e.g., `jaccard=0.65;retention=0.75`) |

**D. `sample_band_paths.tsv`** — one row per (system, sample):

| column | type |
|---|---|
| `system_id` | str |
| `sample_id` | str |
| `band_path` | str | comma-joined track_id per pane (`BT_001,BT_001,—,BT_002,BT_002`) |
| `dominant_track` | str | majority track over all panes the sample appears in |
| `dominant_fraction` | float | fraction of panes assigned to dominant_track |
| `num_switches` | int | number of pane boundaries where sample changed track |
| `path_class` | str | `clean` / `clean_het` / `double_xo` / `nested` / `uncertain` |
| `consensus_band` | str | dominant_track if `dominant_fraction >= 0.85 AND num_switches <= 1`, else NA |

The **consensus_band rule** (Quentin's exact wording, kept verbatim
in code comments): `if dominant_fraction >= 0.85 and num_switches
<= 1: assign consensus_band else: consensus_band = NA`. Without
this guard, complex/A-B-A regions get force-assigned to one band
and the manuscript ends up reporting genuinely ambiguous fish as
clean — exactly the failure mode this spec exists to fix.

**E. `het_skeletons.tsv`** — one row per (system, het skeleton):

| column | type |
|---|---|
| `system_id` | str |
| `skeleton_id` | str | `HET_SK_<chrom>_<system>_<seq>` |
| `chrom` | str |
| `principal` | bool | true for the dominant skeleton, false for `het_candidate` skeletons |
| `het_signal_source` | str | `het_rate` / `dosage` / `pc1_middle` |
| `start_pane`, `end_pane` | int |
| `start_l2`, `end_l2` | int |
| `start_bp`, `end_bp` | int |
| `n_panes` | int (segments + dropouts) |
| `n_segments` | int |
| `n_dropout_panes` | int |
| `dropout_panes` | str (comma-joined) |
| `n_consensus_carriers` | int |
| `consensus_carrier_ids` | str (comma-joined; empty if `> 50` — separate carriers TSV otherwise) |
| `carrier_stability_score` | float ∈ [0, 1] |
| `mean_jaccard`, `mean_overlap` | float |
| `mean_retention_forward`, `mean_retention_reverse` | float |
| `mean_het_enrichment` | float |
| `class` | str | `clean_skeleton` / `skeleton_with_dropout` / `reappearing_skeleton` / `mosaic` / `fragmented` |
| `class_reason` | str |
| `cross_gap_jaccard` | float (NA if single-segment) |
| `notes` | str | free-text flags |

**F. `het_carrier_edges.tsv`** — one row per adjacent het-band
edge (subset of `band_track_edges.tsv` filtered to het bands):

| column | type |
|---|---|
| `system_id`, `chrom` | str |
| `pane_i`, `pane_j` | int |
| `l2_i`, `l2_j` | int |
| `het_class_i`, `het_class_j` | str | `het` / `het_candidate` |
| `n_carriers_i`, `n_carriers_j` | int |
| `overlap` | int |
| `union_size` | int |
| `jaccard` | float |
| `overlap_coef` | float |
| `retention_forward`, `retention_reverse` | float |
| `connected` | bool |
| `connection_reason` | str | `jaccard` / `overlap` / `retention` / `dropout_skip` / `none` |
| `is_dropout_skip` | bool | true if this edge bridges a dropout gap |
| `n_dropout_panes_skipped` | int (0 for direct edges, ≥ 1 for skip edges) |
| `skeleton_id` | str (NA if edge not in any skeleton) |

Note that **F is not a duplicate of B**: B's `band_track_edges`
includes ALL edges (every K-band pair, all panes). F is the
het-anchored projection — only the edges between het / het_candidate
bands plus any dropout-skip edges. Smaller, more readable for
manuscript Q&A.

### 6.3 TSV writer naming + location

Following the band-trace TSV convention (turns 162-164):

```
band_tracks_<chrom>_<system_id>_K<K>.tsv                  (A)
band_track_edges_<chrom>_<system_id>_K<K>.tsv             (B)
system_band_summary_<chrom>_<system_id>_K<K>.tsv          (C)
sample_band_paths_<chrom>_<system_id>_K<K>.tsv            (D)
het_skeletons_<chrom>_<system_id>_K<K>.tsv                (E)
het_carrier_edges_<chrom>_<system_id>_K<K>.tsv            (F)
```

Each writer is a pure serializer (`_bandTracksToTSV(state)` →
string) plus a download orchestrator (`_bandTracksDownloadTSV()`
→ filename, headless-safe per the established pattern). All
serializers share a `# `-prefixed comment header documenting
the parameters (K, JACCARD_MIN, RETENTION_MIN, MAX_GAP_BP,
HET_JACCARD_MIN, HET_OVERLAP_MIN, MAX_DROPOUT_GAP, etc.) so the
supplementary files are self-describing.

---

## 7. UI surfaces in the L3 panel

Three surfaces, in priority order:

### 7.1 Single-band rows under each K×K table (slice 1, primary deliverable)

For each pane boundary in the L3 strip, *below* the existing K×K
contingency table rendered by `_ssContingencyTableHtml`, render
**K extra single-row mini-tables** — one per band cohort at the
left pane:

```
[K×K table for panes (i, i+1)]

g0 (n=56) → 56 0  0  0  0  0    [c0]   ← persists; rf=1.00 jaccard=0.96
g1 (n=39) → 2 38  0  0  0  0    [c1]   ← persists; rf=0.97 jaccard=0.86
g2 (n=30) → 0  0 12 18  0  0    [c2/c3] ← splits; rf=1.00 (12+18=30) jaccard=0.40+0.50
g3 (n=35) → 0  0  0  2  0 33    [c5]   ← jumps label; rf=0.94 jaccard=0.62
g4 (n=8)  → tiny; flagged unstable
g5 (n=12) → ...
```

Each row:
- Labels each destination cell with the K-means group id
- Highlights the destination(s) with the most mass
- Shows `connected=true/false` for each destination via a small
  badge (✓ green for connected, ✗ red for not)
- Tooltip on each cell: full edge stats from §4
- Tooltip on the row: the connection summary (`persists` /
  `splits` / `dissolves`) and which destination(s) the rule fired on

Rendering helper: `_ssBandRowHtml(M, rowIdx, KA, KB, edgeStats)`
mirroring the `_ssContingencyTableHtml` shape; the edgeStats map
is computed once per pane pair and reused across rows.

Default ON for `state.l3SingleBandRows` toggle (persisted to
localStorage), with a checkbox in the L3 toolbar next to the
existing metric chip:

```
[ ] single-band rows
```

### 7.2 Track strip below the L3 strip (slice 2)

A horizontal strip below the bottom of the L3 panel, one row per
detected band track. Each row is a series of colored segments
matching the panes the track occupies; gaps render as dim
hatched zones. Click a track row → highlights the corresponding
cohort's fish in the per-sample-lines panel (fish-set source,
reuses `setBandTraceFishSet`); double-click → opens a popover
with the track's full summary (§5.2 fields).

Visual treatment:
- Continuous tracks: solid fill, color hashed from track_id
- Split tracks: two segments connected by a dashed line at the
  gap, with the gap diagnosis as a small inline icon
- Unstable tracks (small cohorts): reduced opacity
- Hover on a track row: shows core_samples count and class
- Track rows are sorted by start_pane ascending, then by
  n_panes descending (longest first within a start)

### 7.3 System classification chip in the L3 toolbar (slice 3)

Next to the existing Cramér's V chip, a new chip showing the
system classification:

```
[Cramér's V: 0.72 medium ▸]  [class: simple_single_inversion ▸]
```

Click cycles through the alternative classifications the
algorithm considered (when ambiguous between e.g.
`simple_single_inversion` and `complex`, both are evaluated with
their reasons); the chip shows all candidates with their scores
on hover. This surfaces the classifier's uncertainty, which is
genuine and should not be hidden behind a single label.

A second chip next to it shows the active mode and het-signal
source:

```
[mode: het_anchored ▸ source: het_rate]
```

Click cycles between `het_anchored` and `symmetric`. The
`het_signal_source` field is read-only — it reports which
fallback path the algorithm took (`het_rate` when dosage_chunks
present, `dosage` when only mean dosage available, `pc1_middle`
as the heuristic-only fallback). Mode change invalidates the
cache and re-renders all surfaces.

### 7.4 Lines-panel coloring by track membership (slice 2)

The atlas's per-sample-lines panel (§1.5) gets a new
`state.linesColorMode` value: `'band_track'`. When active, each
fish's line is painted at every pane by its band-track
membership color rather than its raw K-means group. Adjacent
panes where the fish stays in the same track render as one
visual stripe; pane boundaries where the fish jumps tracks
render as a color change.

Four sub-modes (cycle via the existing color-mode picker):

| sub-mode | color rule |
|---|---|
| `band_track_all` | every track gets a unique color from a hashed palette; track-membership-per-fish at every pane resolves to that track's color |
| `band_track_het` | the het skeleton paints in a distinguished saturated color (the existing `_HET_RAMP` warm end), all other tracks dim grey, dropouts neutral |
| `band_track_principal_only` | only the principal het skeleton paints; everything else is dim grey |
| `band_track_with_het_overlay` | `band_track_all` but the het skeleton's track is forced to the saturated het color, overriding its hashed color |

Resolver (`_resolveSampleColorByMode('band_track', sampleIdx, paneIdx)`):
- Look up `state.bandTracks[cand.id]` for the active candidate
- Find the track containing `sampleIdx` at `paneIdx`; null if dropout
- Map track to color via the active sub-mode's color rule
- Return RGB string consumable by the lines-panel painter

When no candidate is active (chromosome-wide view), the resolver
returns the neutral grey for all fish — band tracks are
per-candidate.

The lines-panel rendering pipeline is already in place from v3.99
turn 14e+; the resolver body and the four sub-mode entries in
the color-mode picker are the missing pieces.

### 7.5 Het skeleton overlay on the candidate strip + lines panel (slice 2)

Visual overlays for the **principal het skeleton** of the active
candidate, distinct from the §7.4 line coloring:

- **On the candidate strip** (page 1 top): paint a thin
  saturated-warm strip *underneath* the candidate's existing bp
  range marker, spanning the het skeleton's segments only. Gaps
  in the skeleton (dropouts) appear as visible breaks. The
  reader sees: "this candidate's full bp range is X–Y, but the
  het signal is only solid in [start–dropout, dropout–end]
  segments."
- **On the lines panel** (when `band_track_het` is active): the
  het-skeleton's bp range gets a faint warm-tinted background
  band (alpha ≈ 0.08), drawn beneath the line traces. Reader
  sees the het-anchored interval as a visual emphasis, with
  individual fish lines coloured by their actual track
  membership on top.
- **On the L3 strip**: each pane of the L3 strip that's part of
  the het skeleton gets a small warm-tinted dot in its header
  (over the L2 index). Dropout panes get a hollow / outlined dot.
  Skeleton-end panes get a vertical edge marker.

These overlays are **observation-only**, mirroring the
manuscript-grade vocabulary (§2): the warm tint marks "het signal
is detectable here" not "this is the inversion." Tooltips on the
overlay surfaces explicitly say "het skeleton segment" rather
than "inversion interval."

### 7.6 What this is NOT

- **Not a heatmap visualization of the K×K table.** The existing
  table cells are already heat-shaded; adding another heatmap is
  redundant. Single-band rows give the per-row decomposition
  that the K×K view collapses.
- **Not a replacement for the K×K table.** Both render together;
  the K×K view answers "are these panes similar overall?" and
  the per-band rows answer "what happens to each cohort?".
  Different questions, both useful.
- **Not a karyotype assignment.** No band track gets called
  HOM_REF / HET / HOM_INV in this view. Karyotype labels are the
  H-label classifier's job; band tracks are the upstream object
  the H-label classifier could consume.
- **Not multi-system.** Tracks are computed per-L3-strip (per
  candidate). Cross-system tracks ("the same cohort spans I3 and
  I7 on LG28") are the linkage-table job (turn 164) — a different
  primitive operating on confirmed candidates' locked_labels,
  not on raw K-means cohorts.

---

## 8. Implementation slices

### Slice 1 — single-band rows + het-band identifier (~1 turn)

- `_buildBandEdgesForPanes(panes, K)` — computes the edge
  bundle for every adjacent pane pair in the L3 strip, returns
  `{[pane_pair_key]: {[g_i_g_j_key]: edgeStats}}`. Pure compute.
  Edge stats include `jaccard`, `overlap_coef`, `retention_forward`,
  `retention_reverse`, `connected`, `connection_reason`.
- `_classifyHetBandsForPane(l2idx, K, opts)` — per-pane het-band
  identifier per §4A.1. Returns
  `[{group_id, het_enrichment, het_class, het_signal_source}]`
  for each band at the pane. Reuses `_computeHetRateForL2` (turn 129
  primitive) when `dosage_chunks` is loaded; falls back to mean
  dosage; falls back to PC1-middle heuristic. The function is
  pure (takes l2idx + K, reads state for the dosage layer).
- `_ssBandRowHtml(M, rowIdx, KA, KB, edgeStats, hetClass)` —
  renderer helper. Mirrors `_ssContingencyTableHtml`'s style;
  the `hetClass` parameter tints the row's left-margin label
  (warm tint for `het`, dim tint for `het_candidate`).
- `_ssBandRowsForPaneHtml(ct, edgesForPair, leftK, hetClasses)` —
  wraps K band rows for one pane boundary into one container.
- `state.l3SingleBandRows` toggle + localStorage persistence,
  default ON (persisted at `pca_scrubber_v3.l3SingleBandRows`).
- Toolbar checkbox in the L3 head, between the metric chip and
  the recluster-mode picker.
- `renderL3Panel()` calls `_ssBandRowsForPaneHtml` after
  `_ssContingencyTableHtml` when `state.l3SingleBandRows` is on.
- Cache `state.__bandEdgesCache: { [candidate_id + K + mode]: edges }`
  with invalidation on candidate change, K change, or mode change.
- Tests: edge compute on synthetic 6-sample-K=3 fixture;
  het-band identifier on (a) `dosage_chunks`-present case, (b)
  middle-PC1 fallback case; row HTML asserts (including het tint);
  toggle persistence; cache invalidation.

### Slice 2 — band-track extraction + skeleton + lines coloring + track strip (~1.5 turns)

- `_buildBandTracksForCandidate(cand, opts)` — runs the path
  construction (§5.1) symmetric across all bands, returns
  `BandTrack[]` per §5.2 shape. Cached at
  `state.bandTracks[cand.id]`. Accepts `opts.mode: 'symmetric'
  | 'het_anchored'`; when `het_anchored`, calls
  `_classifyHetBandsForPane` per pane and runs the §4A.3
  skeleton extraction with dropout-skip. Default mode is
  `'het_anchored'` when `dosage_chunks` is loaded; falls back
  to `'symmetric'` otherwise. Mode is recorded on the result.
- `_extractHetSkeleton(panes, K, hetClassesByPane)` — the
  §4A.3 algorithm. Returns `{principal: HetSkeleton,
  secondary: HetSkeleton[], warnings: string[]}`. Pure compute.
- `_detectSplitTracks(tracks, opts)` — A-B-A detection per §5.3.
  Returns the merged-tracks list. Same algorithm operates on
  het skeletons (treated as one-track input).
- `_classifyBandSystem(tracks, hetSkeleton, panes)` — system-level
  classifier per §5.4. When `hetSkeleton` is present, it
  dominates the classification (skeleton class → system class
  mapping per §4A.4 + §5.4).
- `_resolveSampleColorByMode('band_track', sampleIdx, paneIdx)`
  — fills in the v3.99 turn 14e+ stub. Reads
  `state.bandTracks[cand.id]`, finds the track containing
  `sampleIdx` at `paneIdx`, returns the color per the active
  sub-mode (`band_track_all` / `band_track_het` /
  `band_track_principal_only` / `band_track_with_het_overlay`).
  Returns null (→ neutral grey) for dropout panes.
- Lines-panel color-mode picker gains the four `band_track_*`
  sub-modes.
- `_drawBandTrackStrip(canvas, opts)` — track strip renderer
  below the L3 panel. Same coordinate-pixel pattern as the
  band-trace strip from turn 161 (CSS pixels throughout,
  fitCanvas-aware DPR scaling). Renders the principal het
  skeleton first (top row, saturated), then secondary
  skeletons, then symmetric tracks below.
- Het skeleton overlays per §7.5: candidate strip warm-tint
  segment paint, lines-panel warm-tint background band,
  L3-strip pane-header dots.
- Track rows clickable → `setBandTraceFishSet(track.core_samples)`
  populates the existing band-trace fish-set; the lines panel
  and linkage table immediately reflect the cohort.
- `_invalidateBandTracksCache()` hooked into existing candidate
  switch path + K change + mode change.
- Tests: synthetic A-B-A fixture; classification on simple /
  nested / split / fragmented fixtures; het-band identifier
  threshold gating (het / het_candidate / non_het); skeleton
  extraction with dropout-skip on synthetic fixture (carriers
  vanish at panes 4-5, reappear at 6 with high overlap → one
  multi-segment skeleton); skeleton with non-bridgeable gap
  (carriers at 6 are different → two separate skeletons);
  lines-panel resolver returns track color for matched
  fish-pane, neutral grey for dropout; click-to-fish-set
  source-pattern check.

### Slice 3 — TSV exports + classification chip + mode chip (~0.5 turn)

- Six serializers (`_bandTracksToTSV`, `_bandTrackEdgesToTSV`,
  `_systemBandSummaryToTSV`, `_sampleBandPathsToTSV`,
  `_hetSkeletonsToTSV`, `_hetCarrierEdgesToTSV`).
- One download orchestrator (`_bandTracksDownloadAllTSV()`)
  emitting all six TSVs with a common timestamp prefix into a
  zipped bundle (or six separate downloads — depending on
  whether headless Blob can zip; the band-trace turns 162-164
  shipped four-file orchestrators successfully).
- Classification chip in the L3 toolbar; cycle on click between
  the candidates with their scores on hover.
- Mode chip in the L3 toolbar (per §7.3): cycles between
  `het_anchored` and `symmetric`. Mode change invalidates the
  band-tracks cache and re-renders.
- TSV-export button in the L3 toolbar (`📊 tracks` next to the
  existing exports).
- Tests: serializer round-trip on synthetic fixtures (all six
  TSVs); consensus_band rule honoured; chip displays correct
  class; mode-chip cycle invalidates cache; het skeleton TSV
  contains skeleton id, segment list, dropout panes,
  classification.

---

## 9. Open questions

1. **K=3 vs K=6.** Quentin works at K=6 in the trigger
   description but the atlas defaults to K=3 (set by
   `state.k`). Should this spec compute at the user's current K,
   or always at K=6? Recommendation: compute at `state.k`
   (whatever the user has set); document in the chip that K=6
   produces denser tracks. K can be added to the cache key so
   switching K bumps the result.

2. **Pane source.** L3 panes currently mean the candidate's
   `l2_indices` array. Should the spec also support arbitrary
   pane lists (e.g. a sliding window of K=6 panes Quentin picks
   manually)? Recommendation: defer. The L3 strip is the natural
   home; a free-form pane picker is a separate feature.

3. **Cross-system band tracks.** A future extension would chain
   tracks across multiple L3 strips on the same chromosome (or
   even across chromosomes). That's `SPEC_distant_band_concordance_fish_trajectory.md`'s
   territory operating at L2 envelopes. Defer this spec to per-strip;
   note the integration point.

4. **Auto-label proposals.** With band tracks plus their core
   samples, the H-label classifier could pre-fill karyotype
   proposals (low-PC1 track → HOM_REF candidate, mid-PC1 → HET
   candidate, etc.). Out of scope here but worth a follow-up
   spec hook in the karyotype tab of the G-panel.

5. **Persistence.** Should `state.bandTracks` survive session
   save/load? Recommendation: no — it's a derived object and
   the in-session compute is fast (~30 panes × K² edges =
   ~1k operations). Persisting derived state introduces
   cache-invalidation surface area for no savings.

6. **Render cost.** Slice 1 adds K extra rows per pane boundary
   (typically 6 panes × 5 boundaries × 6 rows = 180 rows), each
   a small HTML table fragment. Negligible compared to the
   existing K×K rendering. Slice 2's track strip is one canvas;
   ~10 tracks × ~30 panes ≈ 300 rects — also negligible. No
   performance worry expected.

---

## 10. Tests

Source-pattern + behavioural for each slice:

### Slice 1 tests (`test_turnN_band_rows_in_l3.js`):

- `_buildBandEdgesForPanes` — synthetic 6-sample-K=3 fixture.
  Verify edge stats for: clean continuation (jaccard=1, rf=rr=1);
  pure split (jaccard=0.5, rf=1, rr=1); jump-label (high rf+rr,
  Jaccard depends on size).
- Connection rule: jaccard above and below threshold; retention
  pair above and below; OR predicate's two branches both fire.
- Overlap coefficient: subset case (cohort 30→18 all-subset)
  produces overlap=1.0, jaccard=0.6 → connected via overlap branch.
- `_classifyHetBandsForPane` — with `dosage_chunks` mock present
  and synthetic per-sample het rates: middle-PC1 band gets the
  highest mean-het score; gets `het_class === 'het'`. With no
  dosage layer, falls back to PC1-middle heuristic and emits the
  `'no_dosage_layer_falling_back_to_pc1'` warning.
- `_classifyHetBandsForPane` — at K=6 with two plausibly-het
  bands (mean-het 0.45 and 0.40), one is `het`, the other is
  `het_candidate`.
- `_ssBandRowHtml` — output contains correct cell highlights for
  the connection's destination(s); shows `n=N` on the row label;
  `connected` badge present; `hetClass` parameter tints the row
  label.
- Toggle: `state.l3SingleBandRows` reads localStorage on init;
  flipping toggles the row visibility.
- Cache: same `cand.id + K + mode` → same result identity;
  switching K → new compute; switching mode → new compute.

### Slice 2 tests (`test_turnN+1_band_tracks_extraction.js`):

Symmetric mode (existing scope):
- `_buildBandTracksForCandidate` — synthetic 5-pane fixture with
  one full-span track + one internal track produces 2 tracks
  with correct start/end_pane.
- Path construction: split — fixture where g0→{g0,g1} produces
  two tracks sharing the prefix.
- Path construction: merge — fixture where {g0,g1}→g0 records
  `merged_from`.
- A-B-A: synthetic fixture with cohort vanishing at panes 3-4
  and reappearing at pane 5 produces one multi-segment track
  with `class: split_track_reappearing` and the right
  `gap_bp_to_next_segment`.
- A-B-A guard: fixture where the cohort *concentrated on a
  different track* in the gap correctly does NOT auto-merge
  (or auto-merges but flags `gap_diagnosis.concentrated_on_other`).
- Classification: `simple_single_inversion` / `nested` /
  `split_track_reappearing` / `fragmented` / `complex`
  fixtures all return the right class + reason.

Het-anchored mode (new scope):
- `_extractHetSkeleton` on a synthetic 9-pane fixture: het
  carriers persist {1,2,3,4,5} across all panes →
  `class: 'clean_skeleton'`, `n_segments: 1`, `dropout_panes: []`.
- Skeleton with internal dropout: het carriers {1,2,3,4,5} at
  panes 0–3, no het band at panes 4–5, het carriers {1,2,3,4,5}
  at panes 6–8 → one skeleton, `class: 'skeleton_with_dropout'`,
  `dropout_panes: [4, 5]`, `n_segments: 1` (the dropout-skip
  bridges them within MAX_DROPOUT_GAP).
- Skeleton with non-bridgeable dropout: het carriers {1,2,3,4,5}
  at panes 0–3, no het signal panes 4–8 → one skeleton ending
  at pane 3, no skip.
- Two separate skeletons: carriers {1,2,3,4,5} at panes 0–3,
  carriers {20,21,22,23,24} at panes 5–8, gap at pane 4 →
  cross-gap Jaccard = 0 → two separate skeletons (NOT bridged).
- Reappearing skeleton: carriers {1,2,3,4,5} at panes 0–2, gap
  at panes 3–5 (longer than MAX_DROPOUT_GAP), carriers
  {1,2,3,4,5} at panes 6–8 → split into two skeletons by §4A.3,
  then §5.3 split-track detection re-merges them into one
  `class: 'reappearing_skeleton'` with high cross-gap Jaccard.
- Multiple het bands per pane (K=6): synthetic fixture with two
  het-class bands at every pane → produces one principal
  skeleton + one secondary skeleton (`het_skeleton_secondary`
  array length 1).
- Het signal source fallback: `dosage_chunks` absent → uses
  PC1-middle heuristic, sets `het_signal_source: 'pc1_middle'`,
  warnings array contains the fallback message.
- `consensus_carriers` correctness: intersection across all
  non-dropout panes; `consensus_carrier_count` matches.
- `carrier_stability_score` ∈ [0, 1] and equals 1.0 on a
  perfectly-clean skeleton.

Lines coloring:
- `_resolveSampleColorByMode('band_track', sampleIdx, paneIdx)`
  — returns the track color when sampleIdx ∈ track at paneIdx;
  returns null (→ neutral grey) when in a dropout pane.
- Sub-mode `band_track_het`: returns the warm `_HET_RAMP` color
  for samples in the principal het skeleton, dim grey otherwise.
- Sub-mode `band_track_principal_only`: returns null for samples
  not in the principal skeleton.
- Track strip click → `setBandTraceFishSet` called with
  `core_samples` (source-pattern + behavioural via stub).

### Slice 3 tests (`test_turnN+2_band_tracks_export_and_chip.js`):

- Six TSV serializers — round-trip a synthetic compute result
  through each serializer, verify column counts, row counts
  (n_tracks for A; n_edges for B; 1 for C; n_samples for D;
  n_skeletons for E; n_het_edges for F), and the consensus_band
  rule's threshold gating.
- TSV files self-describing: comment header includes K,
  JACCARD_MIN, RETENTION_MIN, MAX_GAP_BP, HET_JACCARD_MIN,
  HET_OVERLAP_MIN, MAX_DROPOUT_GAP.
- Het-skeleton TSV (E): row contains skeleton_id, segment list
  (joined), dropout_panes (joined), classification, cross_gap_jaccard.
- Het-edges TSV (F): one row per het-band edge; `is_dropout_skip`
  set correctly for skip edges; `n_dropout_panes_skipped > 0`
  for skip edges.
- Classification chip: source-pattern check that it reads
  `state.bandTracks[cand.id].classification`; tooltip lists
  alternatives.
- Mode chip: cycle invalidates band-tracks cache; the new
  compute reflects the new mode.

---

## 11. Cross-references

- Existing primitives this spec consumes:
  - `getL2Cluster(l2idx).fixedKLabels` — per-pane labels
  - `_buildContingency(labelsA, labelsB, KA, KB)` — K×K table
  - `_ssContingencyTableHtml(ct)` — HTML rendering hook (line
    12654)
  - `renderL3Panel()` — top-level render (line 48685)
  - `_hungarianChainProjection`, `_concordanceMatrix` (turn 130)
  - `setBandTraceFishSet` (turn 161) — the click-target setter
  - `_bandTraceForFishSet`, `_bandTraceRegimeRuns` (turn 160) —
    single-cohort analogs at L2 granularity, complementary view
  - `_computeHetRateForL2(l2idx)` (turn 129) — per-sample het rate
    `Float32Array[n_samples]` per L2; the primary input to the
    het-band identifier (§4A.1)
  - `_HET_RAMP`, `_hetRateColor` (turn 129) — the warm/cold
    divergent ramp; reused for the het skeleton's saturated
    color (§7.4 `band_track_het` sub-mode + §7.5 overlays)
  - `state.linesColorMode` + `_resolveSampleColorByMode` stub
    (v3.99 turn 14e+) — the lines-panel color resolver hook
    that this spec fills in for `band_track_*` modes

- Specs this spec interacts with:
  - `SPEC_distant_band_concordance_fish_trajectory.md` — dual
    representation (fish paths instead of cohort paths). Both
    operate on the same per-pane partition matrices.
  - `SPEC_l2_sweep_inheritance.md` — system-level merger
    (already shipped). Band-track extraction sits inside an
    L3 strip; the strip itself is a system from the L2-sweep.
  - `SPEC_review_surfaces_auto_and_lineages.md` — classifications
    surface here would extend the dashed-outline "auto" tab with
    a track-class column once both ship.
  - `SPEC_lines_panel_candidate_bands.md` — Slice 2 of that spec
    paints alpha intervals for confirmed candidates; band-track
    output gives a finer-grained intra-candidate version of the
    same idea.
  - `SPEC_g_panel_unified_groups.md` (karyotype tab) — band
    tracks' core_samples could pre-fill manual-group proposals
    in a future hook.

- Cited atlas turns:
  - Turns 130 / 152 / 156: Hungarian + Cramér's V machinery
  - Turns 160-164: band-trace (single-cohort) pipeline
  - Turns 133-134: L2-sweep auto-promote (the system-level
    producer)
  - Turn 165: review-surfaces auto tab

---

## 12. Manuscript figure / supplementary integration

For the Nature Communications submission draft (manuscript
v19→v20 architecture per the project memory):

- **Result 3** (inversion discovery atlas) gets a track-strip
  thumbnail next to each candidate's mini-PCA, showing the
  band-track decomposition. One row per track; colors hashed
  from track_id; multi-segment tracks rendered with the dashed
  gap.
- **Result 4** (validation) consumes `system_band_summary.tsv`
  for the per-candidate classification table (column
  `classification`).
- **Result 5** (diversity / integration for breeding) consumes
  `sample_band_paths.tsv`'s `consensus_band` column where it's
  defined (and `path_class` where it isn't); this is the safe
  per-fish karyotype-call source for downstream breeding-utility
  tables.
- **Supplementary tables** ship all four TSVs as `Table_S<N>`
  attachments.

The advanced panel for the manuscript figure (per Quentin's
"Publication/Atlas" notes in the trigger):

```
Candidate    clean_A    clean_het    clean_B    double_xo    nested    uncertain
INV_LG28     60         106          60         0            0         0
INV_LG12     45         98           52         12           8         11
```

is a direct group-by on `sample_band_paths.tsv`'s `path_class`
column. The simple-band-counts panel:

```
Candidate    band_1    band_2    band_3    mixed
INV_LG28     60        106       60        0
```

is the same group-by on `consensus_band` (with NA → mixed). Both
panels are derivable from the one TSV — the consensus_band rule
in §6.2 ensures the simple panel only counts samples where the
assignment is safe.

---

End of spec.
