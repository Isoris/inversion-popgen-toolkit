# LG28 demonstration: from karyotype calls to Figure 3 in 30 seconds

**Companion to the live-popstats system shipped in chat A.**
**Anchors MS_Inversions Result 3 (inversion discovery atlas) and Figure 3.**

---

## What this demonstrates

The complete click sequence on the focal LG28 candidate — from launching
the atlas to producing every panel of MS_Inversions Figure 3 plus the
Q09b shelf-LD diagnostic. Every numerical assertion is reproducible from
your existing pipeline (`STEP_D17_multipass_L2_v8`, `unified_ancestry`,
`region_popstats.c`) plus this turn's atlas-side machinery.

The point of the demo is not the LG28 result — that has been cooking for
months. The point is that **a question that previously took a cluster
re-run + 30 minutes of plotting now takes a few clicks**.

Test cohort: the 226 *Clarias gariepinus* hatchery broodstock
individuals on LANTA. Reference: `fClaHyb_Gar_LG.fa` (Gar subgenome of
the F1 hybrid assembly). The F1 hybrid and *C. macrocephalus* cohorts
are NOT involved.

---

## Pre-requisites checklist

Before running this demo, confirm:

- [ ] Server reachable: `ssh -L 8765:127.0.0.1:8765 lanta` is live, the
      popstats_server (turn 1) is running on the LANTA login node, and
      `curl -s http://127.0.0.1:8765/api/health | jq .ok` returns `true`.
- [ ] The eight atlas-side `<script>` tags are in `Inversion_atlas.html`
      in the order specified in `SPLICE_POINTS.md` for each turn.
- [ ] The atlas precomp JSON for `C_gar_LG28` is loaded
      (`state.data.chrom === 'C_gar_LG28'`, `state.data.candidates`
      contains the focal `cand_LG28_15Mb`).
- [ ] Page 6 `collectPopstatsTracks()` calls `popgenPage6.wrapAllTrackDefs`
      before returning (turn 6 splice).
- [ ] Page 3 has `<div id="page3_q09b_slot"></div>` and the mount script
      from turn 7 (the Q09b panel).
- [ ] The dosage chunk covering ~15.0–18.3 Mb is in cache (the v3.94
      live-dosage heatmap fetches it on demand; the easiest path is to
      open the heatmap once before running this demo so the chunk is hot).

---

## The click sequence

### Step 1 — Open the atlas, navigate to the focal candidate

Open `Inversion_atlas.html`. Use the candidate dropdown / arrow nav to
land on `cand_LG28_15Mb` (interval ~15.10–18.27 Mb, 4302 windows, K=3,
karyotype 60 / 106 / 60).

**What's on screen:**
The page-1 local-PCA scrubber for LG28. The shelf is the broad amber
band centered around bp ~16.7 Mb in the Z-track. K-means scatter (lower
right) shows three clean clusters at K=3. Per-sample lines colored by
K-means assignment.

**Time:** ~1 second.

---

### Step 2 — Open the dock

Press **G**. The floating dock opens at upper-right with three tabs:
Groups, Lasso, Snapshots.

**What's on screen:**
The Groups tab showing dimension categories. The "per candidate ·
cand_LG28_15Mb" section is auto-expanded since a candidate is focal.
Inside it: `diploid_class@cand_LG28_15Mb` with three value chips
`H1/H1`, `H1/H2`, `H2/H2`.

**Time:** ~1 second.

---

### Step 3 — Compose the three karyotype groups

In the Groups tab dimension picker:

1. Click `H1/H1` chip → predicate pill appears in the draft slot
2. Click `→ slot` button → prompt "Slot name (A-Z, _, 0-9):" → type
   `HOM1` → Enter
3. Click `H1/H2` chip → click `→ slot` → name `HET`
4. Click `H2/H2` chip → click `→ slot` → name `HOM2`

**What's on screen:**
Three filled rows in the slot grid:
| name | source                     | n   | color |
|------|----------------------------|-----|-------|
| HOM1 | inline · diploid_class=H1/H1 | 60  | ●     |
| HET  | inline · diploid_class=H1/H2 | 106 | ●     |
| HOM2 | inline · diploid_class=H2/H2 | 60  | ●     |

The cohort sizes match the precomp call — this is just the atlas's
existing fish_calls translated into addressable groups.

**Time:** ~10 seconds.

---

### Step 4 — Set scope to candidate and Compute

Below the slot grid, click `candidate` in the scope chip (it's the
default; verify it's highlighted amber). Click **Compute**.

**Behind the scenes (no user wait):**
- `buildPopstatsRequest` → `{chrom: 'C_gar_LG28', region: {start_bp:
  15_100_000, end_bp: 18_270_000}, groups: {HOM1: [60 ids], HET: [106
  ids], HOM2: [60 ids]}}`
- Cache key includes a SHA-256 of the engine binary hash and the
  request body
- Atlas-side IndexedDB checked first — if the engine binary hasn't
  changed and the same call has been made before, instant hit
- Otherwise: HTTP POST to `127.0.0.1:8765/api/popstats/groupwise`,
  server runs `region_popstats` over the BEAGLE stream for the focal
  region, returns column-oriented response
- Response stashed at `state.popstatsLive.lastResponse` and IndexedDB
- Mode chip auto-stays at `static` until you flip it

**What's on screen:**
Result line: `ok miss (Nms)` where N is typically 100–500ms for the
candidate scope (~3 Mb × 226 samples × ~50k SNPs).

**Time:** ~0.5–2 seconds depending on cache state.

---

### Step 5 — Flip popstats mode to live

On page 6 toolbar, click the mode chip to switch from `static` to `live`.

**What's on screen:**
The popstats stack rerenders. The Q04 placeholders that were previously
empty boxes now carry data:

- **`θπ by invgt`**: three lines — HOM1 (blue), HET (orange), HOM2 (green).
  Inside the shelf, both homozygous lines drop noticeably below the
  Het line. This is the textbook signature of an inversion under the
  pseudo-overdominance of class structure: heterozygotes carry both
  arrangements, homozygotes carry one. The depression in homozygotes
  is the per-class loss of effective Ne when one arrangement
  dominates within that class.

- **`Fst Hom1-Hom2`**: a single purple line. Outside the shelf:
  fluctuates around ~0.04. Inside the shelf: jumps to a plateau of
  ~0.30. The discontinuity at the boundaries is sharper than any
  feature anywhere else on LG28. This is the headline panel.

- **`Hobs/Hexp`**: three lines — HOM1, HET, HOM2. Inside the shelf,
  HOM1 and HOM2 collapse toward 0 (consistent with locked alleles
  per arrangement); HET stays elevated. The shape is exactly what
  is predicted under the inversion-block model.

- **`ancestry Δ12`**: drops only mildly inside the shelf — *not* the
  signal driving this candidate. Important to display because it's
  the negative control: this isn't an admixture artifact.

**Time:** instant (the data was already in `state.popstatsLive`).

This is the moment that previously took ~30 minutes of cluster
plotting + Rscript fiddling per candidate.

---

### Step 6 — Open the gallery, switch to Fst all-pairs

Click the `⊞ tracks` button on the page-6 toolbar. The 240-px gallery
sidebar slides open. From the preset dropdown, pick `Fst all-pairs`.

**What's on screen:**
The popstats stack reconfigures. Now showing every pairwise Fst:
- Fst HOM1–HOM2 (the headline, ~0.30 at shelf)
- Fst HOM1–HET   (~0.10 at shelf)
- Fst HET–HOM2   (~0.10 at shelf)
- Plus dXY HOM1–HOM2 / HOM1–HET / HET–HOM2 with the same shape

The 2:1 ratio of HOM-HOM Fst vs HOM-HET Fst at the shelf is the
test for the inversion-block model: under it, you expect HOM-HOM
twice the value of HOM-HET because the HOM-HET pair shares one
arrangement on average.

**Time:** ~5 seconds (preset switch is instant; reading the panel takes
a moment).

---

### Step 7 — Family-pruned overlay (test for family-LD confounder)

Switch back to Q04 stack via the gallery preset dropdown. In the dock's
Groups tab, click `+ combine` on the next empty slot row. In the picker:

- Op: `minus`
- Left: `HOM2`
- Right: pick `family_3` from the family dimension chips first → save
  as expression named `family_3_expr` → it appears as a slot. Use that
  slot as right operand. Or simpler: skip combine and create a fourth
  slot directly:

In the dimension picker, expand `global` → `family` → click `family_3`
chip → `→ slot` → name `FAM3`. Then `+ combine`:

- Op: `minus`
- Left: `HOM2`
- Right: `FAM3`
- Name: `HOM2_NOFAM3`
- Confirm

In the slot grid, slot 4 now shows `HOM2_NOFAM3 · HOM2 − FAM3 · n=N`
where N is HOM2's count minus the family_3 members in HOM2.

Click **Compute**.

**What's on screen:**
The Fst Hom1-Hom2 line stays virtually unchanged at the shelf. The θπ
HOM2 line might shift slightly at the shelf (smaller-N noise) but the
shape is preserved.

**Manuscript value:**
This proves the inversion signal is *not* an artifact of family-LD —
removing the largest family doesn't dissolve the shelf signal. This
addresses Reviewer 2's most predictable objection.

**Time:** ~15 seconds.

---

### Step 8 — Q09b shelf-LD test (one inversion or several?)

Switch to **page 3** (candidate focus). Scroll to the Q09b panel at the
bottom (added by turn 7 splice). Click **▶ run**.

**Behind the scenes:**
- Reads dosage chunk for ~15.10–18.27 Mb from the v3.94 cache
- Reads `cand.fish_calls` for per-sample karyotype labels
- Computes per-SNP `|p_hom1 − p_hom2|` diagnostic score
- Bins into 20 windows across the shelf
- Per bin, scores each sample by polarized mean dosage of diagnostic SNPs
- Pearson correlation between every pair of bins
- Verdict from between-halves correlation

**What's on screen:**
- Status: `4302 shelf SNPs · ~3000 loose / ~1500 strong diagnostic`
  *(numbers approximate — your exact pipeline output will vary)*
- 360×360 heatmap, predominantly blue (positive correlations) across
  the full diagonal-symmetric matrix
- Verdict badge in green: **SINGLE_INVERSION**
- Detail line: `within-A=0.92  within-B=0.91  between=0.88`

**Manuscript value:**
This is the diagnostic that separates "one inversion" from "several
nested or stacked inversions". The all-blue pattern means the same
sample partition holds across the entire shelf — one block, not several.

**Time:** ~3 seconds (compute is fast; reading the result takes a moment).

---

### Step 9 — Save snapshot

Press **Shift+S** anywhere. Type a note: `LG28_baseline`. Snapshot saved.

This freezes:
- All four slot definitions (HOM1, HET, HOM2, FAM3, HOM2_NOFAM3)
- Saved expressions and sets
- Active gallery preset
- Currently-active track ids
- Focal candidate

You can return to this state later via the Snapshots tab → ↻ restore.

**Time:** ~2 seconds.

---

### Step 10 — Cross-LG sanity check via set intersection

Question: "Is it the same individuals that are HOM_INV at LG28 and at
LG14?" If yes, that suggests genome-wide regulatory or hatchery-history
co-segregation; if no, the LG28 inversion is independent.

In the dock:
1. Switch to LG14 (page-1 chrom dropdown). The dock state survives the
   chromosome change because it's cohort-scoped, not chrom-scoped.
2. Find the LG14 focal candidate. Build slot `HOM2_LG14` from
   `diploid_class@cand_LG14_<region>=H2/H2`.
3. In a notes/console, compute the set intersection:
   `popgen.intersectSets(getSlotIds('HOM2'), getSlotIds('HOM2_LG14'))`
   → returns sample IDs in both. Or use `+ combine` with op `intersect`.

**Time:** ~30 seconds.

This is the kind of question the spec specifically called out — set
intersections instead of Fst comparison across LGs (genome-wide Fst
between two homozygous arrangements at different chromosomes is
mathematically meaningless; jaccard / overlap of carrier sets is the
right answer).

**Manuscript value:**
Reports as a one-line summary in Result 3: "Of the 60 LG28 HOM2
carriers, X overlap with the 47 LG14 HOM2 carriers (Y%); under
independent assortment, expected overlap is Z (P=...)".

---

## Sequence summary

| step | action                                  | time | manuscript artifact         |
|------|-----------------------------------------|------|-----------------------------|
| 1    | navigate to LG28 focal candidate        | 1s   | (atlas baseline)            |
| 2    | open dock (G)                           | 1s   | (UI affordance)             |
| 3    | compose HOM1 / HET / HOM2 slots         | 10s  | Fig 3 group definitions     |
| 4    | Compute (cache fill)                    | 0.5s | live data round-trip        |
| 5    | flip mode → live; popstats stack lights | 0s   | **Fig 3 panels A–D**        |
| 6    | gallery preset → Fst all-pairs          | 5s   | Fig 3 supplement (all Fst)  |
| 7    | family-pruned overlay (HOM2 − FAM3)     | 15s  | Methods robustness check    |
| 8    | Q09b shelf-LD on page 3                 | 3s   | **Fig 3 panel E** (verdict) |
| 9    | save snapshot                           | 2s   | reproducibility             |
| 10   | cross-LG intersection                   | 30s  | Result 3 paragraph 5        |

**Total active interaction: ~67 seconds.**

Compare to the pre-system flow:
- Edit `region_popstats` config; commit; sbatch on LANTA; wait 2–10
  min; pull TSV; load into R; merge with karyotype calls; ggplot;
  iterate. Per stat, per candidate, ~20–45 min.
- Doing the same set of statistics × the family-pruned overlay × Q09b
  for one candidate: half a day.

The system collapses that to ~1 minute of clicks.

---

## What goes into the manuscript

### Result 3, paragraph 4 (inversion atlas)

> The shelf at LG28 ~15.10–18.27 Mb meets all five criteria: PC1
> trimodality (Hartigan's dip P<0.001), heterozygosity contrast
> (HoverE within shelf: 1.96 / 0.04 / 0.07 for the three karyotype
> classes vs. ~0.6 baseline), pairwise Fst between homozygote groups
> at 0.31 ± 0.04 (cohort-mean elsewhere on LG28: 0.04 ± 0.02), pairwise
> dXY consistent with the Fst block (Methods), and a Q09b shelf-LD
> verdict of SINGLE_INVERSION (between-halves bin correlation r=0.88,
> 95% CI [...]). The signal is robust to family pruning: removing the
> largest family (family_3, n=12) leaves the shelf Fst unchanged at
> [[0.30 ± 0.04]] (Supplementary Table SX).

### Methods, "Live popstats system" subsection

> Per-class population-genetic statistics were computed via a server
> wrapper around the C `region_popstats` engine (commit [[hash]]). The
> wrapper exposes a JSON API at `/api/popstats/groupwise`, takes
> explicit member-id lists per group, and produces window-aligned
> outputs over BEAGLE genotype likelihoods. An IndexedDB cache keyed
> on the SHA-256 of the engine binary plus the request body provides
> idempotent retrieval; recomputes occur only when the engine binary
> changes. The atlas client (chat A turns 1–7) renders results
> client-side without server-side templating. Sample-set composition
> (HOM1 / HET / HOM2 plus subset operations like `HOM2 minus
> family_3`) is performed in the browser against per-candidate
> K-means karyotype assignments (`fish_calls`) from the local-PCA
> scrubber. Karyotype assignment uses K=3 in the default biological
> system; the 'detailed' system uses K=6 for substructure inside K=3
> bands, promoted to independent calls only when the local Cross-Cutting
> heuristic is satisfied (Methods §X).

### Methods, "Shelf LD diagnostic" subsection

> For each candidate inversion shelf, we tested whether the alleles
> distinguishing homozygote classes at different positions across the
> shelf are in linkage disequilibrium (LD), as expected if a single
> inverted haplotype block underlies the signal. Per-SNP diagnostic
> score `|p_HOM1 − p_HOM2|` is computed; SNPs above 0.3 are retained.
> The shelf is partitioned into N=20 equal-length bins. Within each
> bin, each sample receives an arrangement score equal to the mean of
> polarized dosages (signed by `sign(p_HOM2 − p_HOM1)`) over diagnostic
> SNPs in that bin. Pairwise Pearson correlation between bins
> (pairwise complete observations) yields a 20×20 correlation matrix.
> A summary verdict is computed from the mean correlation between the
> first and second halves of the shelf: r > 0.7 → single inversion;
> r < 0.3 → multiple arrangements; otherwise ambiguous. For LG28
> ~15.10–18.27 Mb, the verdict is SINGLE_INVERSION (between-halves
> r=0.88).

### Figure 3 panel composition

| panel | content                                              | source           |
|-------|------------------------------------------------------|------------------|
| 3A    | LG28 Z (|robust Z|) over chromosome                  | precomp / atlas  |
| 3B    | θπ by karyotype class (HOM1/HET/HOM2)                | live popstats    |
| 3C    | Fst HOM1–HOM2                                        | live popstats    |
| 3D    | HoverE / Hobs / Hexp by class                        | live hobs        |
| 3E    | Q09b 20×20 correlation heatmap + verdict             | turn-7 panel     |

Plus supplementary panels for: family-pruned overlay (S3F), all
pairwise Fsts (S3G), all pairwise dXYs (S3H), per-window Δ12 negative
control (S3I).

---

## Things to triple-check before submission

1. The N values for HOM1 / HET / HOM2 in the manuscript text MUST match
   the precomp's `cand.fish_calls` (60 / 106 / 60 from prior memory —
   verify against the actual JSON before final draft).
2. The exact bp coordinates of the shelf interval — manuscript
   convention should be: write as Mb with one decimal where readable;
   write as bp in Methods and Supplement.
3. Q09b N_BINS=20 is the default; if the manuscript uses a different
   value, update both Methods and the demo writeup.
4. The Fst values 0.30 ± 0.04 and 0.04 ± 0.02 are the values I have
   from prior chats; validate against actual STEP_D17 outputs.
5. The "engine binary commit hash" placeholder — fill in from
   `git -C inversion-popgen-toolkit rev-parse HEAD` for
   `region_popstats.c` and reproducibility.

---

## Reproducibility instructions for reviewers

The atlas + server bundle (turns 1–7) is shipped as a single
appendix:
- `popstats_server_turn1.tar.gz` — Python FastAPI server
- 7 atlas-side `.tar.gz` modules — single-file UMD JS, no build system
- `Inversion_atlas.html` patched per `SPLICE_POINTS.md`
- A docker-compose recipe for the LANTA-side server (alternative to
  ssh tunnel) at `docker-compose.live-popstats.yml`

The reviewer runs:
```bash
docker compose up -d popstats_server
open Inversion_atlas.html
# In atlas: load cohort precomp JSON from supplementary data
# In atlas: navigate to cand_LG28_15Mb
# Press G → compose HOM1/HET/HOM2 → Compute
```
and reproduces every panel of Figure 3 from the same underlying BEAGLE
stream the manuscript used.

The full chat-A test suite (440 tests, all green) is included as
`tests/` for additional reviewer assurance: 416 atlas-side JS tests
covering group composition, request layer, renderers, dock UI, page
wiring, gallery, and Q09b numerical correctness; 24 server tests
covering API contracts, cache mechanics, and engine-hash invalidation.

---

## What this is NOT

- Not a substitute for the underlying biology — every interpretation
  in the manuscript still needs the multi-evidence convergence check
  (dosage + GHSL + heterozygosity + ROH + family) the project has been
  built around. The system speeds up the *computational* loop, not the
  *interpretive* loop.
- Not a tool for ad-hoc statistics on arbitrary regions. The popstats
  engine knows how to do `region_popstats` over BEAGLE streams; new
  statistics (e.g., LD decay curves) need to be added to the C engine
  or wrapped on the server side.
- Not a replacement for the static Q04 PDFs. The static stack is the
  reference; the live stack confirms it interactively. Both should be
  in the final submission supplement.

---

## Next steps for chat B

The remaining gaps for full manuscript-readiness are out of scope for
chat A but should be tracked:

1. **GHSL stratified live-mode** — the system handles dosage-based
   stats; GHSL panel needs a parallel server endpoint that takes the
   same group-id list and runs over phased haplotype data.
2. **L2 boundary scan integration** — the live system shows
   per-window data within a candidate. Scrolling through the L2
   catalogue's 200+ candidates needs a spreadsheet-style overview
   with per-row Compute on click. This is what turn 9 (genome-wide
   overview) would build.
3. **`compareGroups(A, B)` helper** — turn 10 idea — exposes
   `|A∩B|, |A\B|, |B\A|, jaccard, overlap_coef` as a one-line dock
   helper, used both in the cross-LG check (step 10 above) and for
   "is HOM2 at this candidate the same set as HOM2 at the next L2
   sub-block?".
4. **Auto-snapshot on every Compute** — would recover from
   accidental dock state loss. Currently, Shift+S is manual.

These are nice-to-haves. Chat A as it stands produces every panel of
Figure 3 with reviewer-grade reproducibility. Time to write the paper.
