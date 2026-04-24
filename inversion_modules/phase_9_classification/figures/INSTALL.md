# Chat 15 — figures 1a, 2, 3a (partial) + test harness

## Files

All four drop into:

```
inversion_modules/phase_9_classification/figures/
```

| File | Status | Purpose |
|---|---|---|
| `fig4e_gene_count_by_biotype.R` | static-reviewed | Figure 2: biotype × region gene-count table-heatmap |
| `fig4e_class_along_inversion.R` | static-reviewed | Figure 1a (pragmatic): karyotype-class proportion along an inversion |
| `fig4e_cumulative_burden_per_group.R` | static-reviewed | Figure 3a: cumulative del-gene discovery, stacked by group |
| `test_figures.R` | new | Smoke test that runs all three on synthetic data |

## Critical context

**Figure 1a is the pragmatic analog of what you originally asked for.**
The "ancestry proportion along the inversion" figure as specced would
require a per-sample × per-window local-ancestry Q matrix. Module 2B
in your repo runs genome-wide NGSadmix only — no per-window local
ancestry is computed. Building figure 1a against data that doesn't
exist is aspirational-ship, which I've been doing too much of.

Instead, `fig4e_class_along_inversion.R` uses the **GHSL v6 per-window
karyotype calls** that do exist (`snake3v6_karyotype_calls.tsv.gz`)
and produces the same visual grammar — stacked area of class
proportions along genomic position. The scientific question is
related: "does the cohort partition shift across the inversion," i.e.
is this composite?

If you want the true ancestry-along-inversion figure, the upstream
step is a local-ancestry compute (ELAI / RFMix / windowed NGSadmix
over the inversion ± flanks). That script doesn't exist in the repo
yet. Adding it is its own chat.

**Figure 3a is wired for parameterization into 3b and 3c.**
Currently `--entity gene` is the only implemented mode.
`--entity variant` and `--entity gene_inv` will error with a message
pointing at the additional inputs they'd need. This is deliberate:
I don't want to ship 3b and 3c unverified when the data contracts
aren't confirmed.

## First thing to do on LANTA

Before you run any of this on real data:

```bash
cd inversion_modules/phase_9_classification/figures/
Rscript test_figures.R
```

This is the **first test harness that actually exercises these
scripts end-to-end.** The last three chats accumulated unverified
scripts without a test cycle. If `test_figures.R` fails, feed the
error back and I'll fix in-place before we add anything new.

## Schema assumptions (pinned here in one place)

These are the column assumptions across the three scripts. If any
differ from your real files, tell me and I'll fix in-place.

### `gene_map.tsv` (fig 2)
From MODULE_CONSERVATION STEP 00:
```
gene_id  chr  gene_start  gene_end  biotype  [strand cds_length is_canonical]
```
Extra columns ignored. Multi-transcript duplicates deduped on gene_id.

### `regions.tsv` (fig 2)
User-supplied:
```
region_id  chrom  start  end  [category]  [label]
```

### `candidate_bed` (fig 1a)
```
chrom  start  end
```
First row used; others ignored.

### `karyo_calls` (fig 1a, from STEP_C04b_snake3_ghsl_classify.R)
```
sample_id  window_start  window_end  n_windows  call  mean_rank
```
`call` ∈ {INV_INV, INV_nonINV}. Window indices are per-chrom (NOT
per-region slice).

### `window_coords` (fig 1a)
Per-chrom, NOT pre-subset to a region:
```
chrom  global_window_id  start_bp  end_bp
```

### `burden_roh` (fig 3a, from MODULE_CONSERVATION STEP 15)
```
sample_id  chr  roh_start  roh_end  roh_length  n_del_variants
  n_del_hom  sum_priority_score  n_classA  n_classB  genes_hit
```
`genes_hit` is `;`-joined; "none" for empty.

### `sample_groups` (fig 3a)
Two columns, no header:
```
sample_id  group_id
```

## Honest-status notes

Status ladder: *not run* → *smoke test passed* → *run on real data
OK*. Everything in this tarball is at level 1 (not run). The test
harness is the path to level 2. I can't take anything to level 3
without your feedback.

Bugs I caught and fixed in static review this turn:

- Fig 2: `colorRamp2` degenerate case (max_val < 2). Fixed with
  2-point ramp fallback.
- Fig 2: `setNames(...[seq_along(unique(cat_vec))], ...)` silently
  gave NA colors for 6+ categories. Fixed with a recycle-safe
  palette.
- Fig 3a: `by = .(row_idx = .I)` trick for exploding `genes_hit`
  depended on subtle data.table semantics I wasn't sure about. Fixed
  with an explicit `rep() + strsplit()` that has obvious semantics.
- Fig 3a: first-hit attribution tiebreak was row-order-dependent and
  therefore not reproducible across runs. Fixed with
  sample_groups-file-order tiebreak.
- Fig 2: zero-overlap case (no genes in any region) now writes a
  placeholder PDF with an informative message instead of erroring.

Things I didn't do because they'd be aspirational:

- `--entity variant` (3b) and `--entity gene_inv` (3c) are stubbed
  with error messages pointing at the additional inputs. Will be
  implemented once 3a is verified on real data.
- Figure 1b (group-composition-per-sub-inversion) not built this
  turn. You said sequential, and 1a needs to be verified first.
- The chat-14 tubemap preview (`fig4e_recombinant_tubemap.R`) is
  still unverified. Schema confirmations needed.

## Where each figure came from

- **Figure 2 (gene biotype table):** Nature paper's ROH Chr 6 vs
  Chr 12 gene-count side-table, generalized to N user-supplied
  regions.
- **Figure 1a (class-along-inversion):** inforiver "5 types of
  stacked charts" visual grammar; specifically the US music sales
  stacked-area chart. Adapted to show karyotype class proportions
  along an inversion interval.
- **Figure 3a (cumulative burden by group):** Nature paper's
  Afghani/Jordanian/… cumulative gene count panel. Direct adaptation
  with `--attribution first` being the default that matches the paper's
  stack semantics.
