# Step 6 design note — dosage × karyotype, shape, and Procrustes

**Date:** 2026-05-04
**Authored:** end of step-7 session, in response to Quentin's three questions.
**Read alongside:** `HANDOFF_2026-05-04_FINAL.md` (status) and `from_turn129/S6_dosage_heatmap_streaming_viewer.md` (existing renderer contract).

This note answers three architectural questions before step 6 (sample × SV heatmap) gets built. Settle these first; the implementation is ~30 minutes once the design is locked.

---

## Q1 — "Procrustes / PCA-on-the-side: nonsense or useful?"

**Verdict: nonsense for this page.** Skip it.

### Why ChatGPT's suggestion doesn't fit

Procrustes asks "do these two configurations of the same samples agree after rotation/scaling/translation?" It's the right tool when you have **two independent ordinations** and want to test whether they tell the same story:

- SNP PCA × SV PCA — useful, both are independent reductions.
- Genotype PCA × phenotype PCA — useful, asks whether genetics predicts phenotype.
- Population A's PCA × Population B's PCA — useful for cross-cohort comparison.

For the SV evidence page, the karyotype labels (H1/H1, H1/H2, H2/H2) come from page 21's locked NGSadmix K=3 — they were locked using the same SV-supporting markers that Fisher tests now run against. A Procrustes on PCA(dosage) vs PCA(karyotype-encoding) would be asking *"do the SV calls recover the karyotype that defined them?"* That's circular. You'd be cross-validating with the same labels twice.

### What's actually useful

The page already shows the structure-vs-dosage relationship more directly through the **karyotype-banded heatmap** (step 6's main view). When samples are sorted by karyotype group:

- Real candidate → three near-uniform horizontal bands.
- Noisy candidate / artefact → static-like.

A reviewer can read this in two seconds. No statistic needed for the figure itself.

If a reviewer demands a number, the right one is **silhouette score of karyotype-group assignment in dosage PCA space** — measures whether the three groups separate well in PCA coordinates derived from the SV-dosage matrix. One scalar, defensible, easy to report. Specifically:

```python
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

pcs = PCA(n_components=2).fit_transform(dosage_matrix)  # n_samples × 2
sil = silhouette_score(pcs, karyotype_labels)            # one float in [-1, 1]
```

Range interpretation: ≥ 0.5 = strong separation; 0.25-0.5 = moderate; < 0.25 = weak. Add as a small text annotation on the heatmap header if at all.

But **no Procrustes**, no separate PCA panel. The heatmap is the figure.

---

## Q2 — "Do SV/breakpoint coordinates have the same shape as karyotype coordinates?"

**No.** They share the sample axis only. This is the core architectural question and the answer drives the heatmap design.

### Shapes

| layer | row index | col index | values | shape |
|-------|-----------|-----------|--------|-------|
| Karyotype (page 21 lock) | sample_id | — | label ∈ {H1/H1, H1/H2, H2/H2} | `(n_samples,)` |
| SV dosage / support | sample_id | sv_id | int ∈ {0=AA, 1=AB, 2=BB, -1=miss} | `(n_samples, n_svs)` |
| SV genotype counts (`sv_genotype_counts_v1`) | sv_id | group × {AA,AB,BB,miss} | int count | `(n_svs, 3, 4)` aggregate |
| Combinations (`sv_evidence_combinations_v1`) | combination | members + samples + size | mixed | per-combination |

So:

- **Karyotype is 1D**, indexed by sample.
- **SV support is 2D**, indexed by (sample, SV).
- **Genotype counts is the 2D matrix collapsed across the sample axis by group**, losing per-sample info.
- **Combinations is the same 2D matrix collapsed across the SV axis by evidence pattern**, losing per-SV info.

Karyotype and SV support share **only** the `sample_id` axis. That's the join.

### The relationship that matters

**SV support is the high-resolution view; karyotype is the low-resolution summary.**

You can compute karyotype from SV-support (it's the dominant haplotype across boundary-marker SVs, basically what NGSadmix did) but you can't reconstruct SV-support from karyotype (information is lost). One direction only.

This means in the heatmap:
- **Rows = samples**, sortable by karyotype group → three banded blocks.
- **Columns = SVs**, sortable by position / FDR / pattern label.
- **Karyotype shows up as a left-side row-annotation strip** — three coloured swatches running down the left edge, one cell per row, encoding which group each sample is in.

That ordering is the trick that makes the figure readable.

### Why both views (UpSet + heatmap) belong on this page

Both read **SV support** but collapse it differently:

- **UpSet** (step 5, shipped): collapse by *evidence-type pattern* across samples → "which evidence patterns exist, and how many fish in each".
- **Heatmap** (step 6, TODO): no collapse — show the full 2D matrix → "which sample carries which SV at what dosage".

They answer different questions about the same underlying matrix. UpSet gives a per-candidate summary ("is this candidate well-supported by which evidence layer?"); heatmap gives the per-sample × per-SV detail ("which fish exactly is in which substructure?").

---

## Q3 — Layer for step 6: schema decision

The third per-candidate layer needs to carry the per-sample × per-SV matrix. The S6 spec already specs this for the **chromosome-wide** dosage viewer (chunked, indexed, streaming for performance). For the SV-evidence page, we need the **candidate-scoped** version.

### Two options

**Option A — Reuse S6's dosage chunk format.** The atlas page reads from `state.data.dosage_chunks` directly, scoping to the candidate window. Pros: zero new schema, the bridge backend in S6 spec gives us the data. Cons: chunk format is per-chromosome, not per-candidate; we'd be reading more than we need; the per-candidate folder structure breaks.

**Option B — New per-candidate layer `sv_support_by_sample_v1`.** Producer-emitted JSON, sibling of `sv_genotype_counts.json` and `sv_evidence_combinations.json` in the per-candidate folder. Pros: per-candidate scoping, drag-drop one folder gets all three layers, clean local-first workflow. Cons: new schema, new producer, slight duplication of the dosage chunk's data.

**Recommendation: B.**

The per-candidate folder is Quentin's preferred mental model and the existing UpSet bar already reads from it. The SV-evidence page is candidate-scoped throughout — adding a chromosome-wide reader would break that consistency. The slight duplication with the dosage chunks is fine because the candidate scope is small (typically 200-2000 SVs × 226 samples = 45k-450k cells, well under 1 MB JSON).

There's already a fixture file in the test tree (`fixture_sv_support_by_sample_v1.json`) from a prior session, suggesting this was the consensus direction. Use that as the schema seed.

### Proposed schema

There's already a thought-out fixture at `tests/sv_evidence/fixture_sv_support_by_sample_v1.json` (from a prior session — clearly someone hashed this out before). It uses a smarter compact form than I initially proposed:

```json
{
  "format_version": "sv_support_by_sample_v1",
  "candidate_id":   "INV_LG28_003",
  "encoding":       "0=AA, 1=AB, 2=BB, -1=miss",
  "samples":        ["CGA_001", "CGA_002", "..."],   // length = n_samples; canonical row order
  "sv_ids":         ["SV001", "SV002", "..."],         // length = n_svs;     canonical col order
  "row_groups":     {                                  // karyotype banding for fast row sorting
    "H1/H1": [0,    60],                               //   [start_idx, end_idx] inclusive
    "H1/H2": [61,  163],
    "H2/H2": [164, 225]
  },
  "dosage_compact": [                                  // ROW-AS-STRING form
    "01210...",                                        //   '0'=AA, '1'=AB, '2'=BB, '.'=miss
    "11201...",                                        //   one char per SV; row length = n_svs
    ...                                                //   length = n_samples
  ]
}
```

**Why row-as-string is the right call here:**

For 226 × 8 SVs:
- Nested int array `[[0,1,2,...], ...]`: ~5-6 KB JSON
- Compact string `["01210...", ...]`: ~2 KB JSON

For 226 × 1000 SVs:
- Nested int array: ~1 MB
- Compact string: ~226 KB

The compact form is 4-5× smaller because there's no per-cell `,` overhead. Browser-side parsing is also faster — one `JSON.parse` then per-row `.charAt(j)` indexing, no nested array traversal.

**When to switch to nested arrays:** when n_svs > ~50 *and* the producer wants per-cell missing-ness encoding richer than a single char. The fixture's own comment says this. For now, compact strings handle 99% of candidates.

The `row_groups` key is a per-group `[start_idx, end_idx]` index into `samples`. Pre-sorted by the producer so the heatmap can render bands with three slice operations — no per-row group lookups in the hot path.

### Producer

`STEP_SV_SUPPORT_emit_support_by_sample.py`. Inputs: same VCF as `STEP_SV_GT_AGG`, candidate JSON, karyotype TSV (just for sample-order canonicalisation). Output: the JSON above.

This producer is mechanically simple — just transcribe per-sample × per-SV genotypes. Most of the work is in `STEP_SV_GT_AGG` already (the GT parser); we'd factor out a shared `_parse_calls_and_samples()` helper.

---

## Q4 (implicit) — How does the heatmap unlock the deferred features?

Once `_state.supportLayer` is loaded:

1. **Recount table genotype columns within `selectedSamples`.** When the user clicks an UpSet bar, walk the matrix rows for samples in selectedSamples, count AA/AB/BB/miss per SV per karyotype group, replace the table's `genotype_counts` columns. ~20 lines.
2. **Dim non-selected SV glyphs in locus.** For each SV, compute "does this SV have any carriers in selectedSamples?" by reading the matrix column. If no carriers in selection, dim the glyph. ~10 lines.
3. **Per-cell hover tooltip on heatmap.** Just read `matrix[i][j]` on hover. Native to the heatmap's own renderer.

All three become trivial. This is why step 6 is the unlock.

---

## What step 6 should do, concretely

In order:

1. **Validate** the `sv_support_by_sample_v1` schema with `_validateSupportLayer` (mirrors `_validateCombinationsLayer`).
2. **Load** via the same folder-walk drag-drop: `_ingestJsonText` adds a third branch.
3. **Render** the heatmap in two views:
   - **Right-rail compact** (300px-wide thumbnail below UpSet): aggregated, e.g. mean dosage per group × SV → shows three rows × n_svs heat.
   - **Main-area large**: full samples × SVs, sortable rows/cols, hover tooltips, click-to-highlight, karyotype-banded sample ordering by default.
4. **Wire UpSet selection → heatmap dim**: when `selectedSamples` is set, dim non-selected sample rows.
5. **Wire heatmap cell click → locus glyph highlight + table row highlight**: the cell's SV becomes `_state.highlightedSvId`.
6. **Unlock the three deferred features** (table recount, locus glyph dim, hover tooltip).
7. **Producer:** `STEP_SV_SUPPORT_emit_support_by_sample.py` + extend `run_sv_evidence_pipeline.slurm` to run all three.
8. **Test:** extend `tests/producers/test_step7_producers.py` to also produce + validate the support layer.

That's step 6 + the related step-7 extension in one bundle, matching the same "per-step bundle" rhythm we've used.

---

## Quick reference: what *not* to do

- **No Procrustes panel.** The heatmap with karyotype-banded ordering is the figure.
- **No PCA-on-the-side panel.** Same reason — circular with the karyotype lock.
- **No "ChatGPT said X" without checking against the spec.** ChatGPT doesn't have the manuscript context or the existing locked-labels framework.
- **Don't merge SV-support with the dosage chunks.** Keep candidate-scoped per-candidate-folder layers.
- **Don't try to reconstruct karyotype from SV-support in the page.** Read the lock from `state.candidateState[cid].locked_labels`, same as the rest of the atlas does.

---

## When you start step 6, the first thing to do is

Decide between Option A (reuse dosage chunks) and Option B (new per-candidate layer). I'm recommending B above, but Quentin should sign off before any code is written. If B, the existing `fixture_sv_support_by_sample_v1.json` in `tests/sv_evidence/` is the schema seed — open it, check it matches the schema in §Q3 above, then proceed.
