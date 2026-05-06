# SPEC — Dosage-based recombinant detector (per-sample changepoint scan)

**Status**: drafted turn 130 final session. Not yet implemented in
the atlas. **Quentin already has a LANTA-side script** (`STEP_C01h_recombinant_scanner.R`
from chat `b7bc2608`). This spec covers atlas-side **visualization +
review** of recombinant calls, plus the architectural decision Quentin
made to **replace** the existing PCA-based recombinant signal in
phase 4b.

**Trigger** (Quentin, chat `7d167416`):
> *"Its nonsense to use PCA for that. Its ok but like ok so if its a
> recombinant normally it moves from like INV INV to INV non INV right?
> Or from INV INV to ref ref its basically a sample that was in a band
> and close to other samples but moved for some intervals right."*
>
> *"The decompose module should be upgraded to use dosage-based
> changepoint detection as the primary recombinant signal, with GHSL
> as confirmation and PCA demoted to a diagnostic-only role, and the
> unsupervised k-means fallback should be removed entirely."*

A `DECOMPOSE_UPGRADE_DESIGN_chat9.md` exists from that chat with the
full design. This spec brings the design into the atlas spec folder
so future sessions don't lose it.

---

## 1. Why PCA-based recombinant detection is wrong

PCA finds axes of maximum variance across samples per window. It's
designed to separate **groups** (REF vs HET vs INV), not to detect
**within-sample transitions**.

The current pipeline:
1. Compute per-window PC1 score per sample.
2. Collapse across windows (`colMeans(pc1_mat)`) to get one PC1 per
   sample → assigns to k-means band.
3. Re-assign each window to nearest of the three global cluster
   centers → categorical label per sample per window.
4. `find_switch_segment` looks for runs of label disagreement.

This converts a real-valued changepoint problem (does the dosage track
have a step?) into a categorical run-length problem on a noisy
re-labeled sequence. At 9× coverage PC1 noise easily flips windows
between HET and HOM centers when the sample's dosage didn't actually
change. **False positives are systematic.**

## 2. The right algorithm — dosage-based changepoint

For each sample, you have a dosage vector along the candidate interval
— one value per SNP, `d_i(pos)` for sample i. A **non-recombinant**
sample's dosage is approximately constant across the interval. A
**recombinant** sample has a step in `d_i(pos)`.

### 2.1 Per-sample changepoint statistic

At each candidate split position p:

```
left_mean  = mean(d_i(pos)) for pos < p
right_mean = mean(d_i(pos)) for pos >= p
delta(p)   = |left_mean - right_mean|
```

Scan p across all positions in the candidate interval; the maximum
delta is the sample's strongest changepoint.

### 2.2 Significance via permutation null

For each sample, generate B=1000 permutations of its dosage vector
(shuffle SNP positions). Compute `max(delta)` for each permutation.
The permutation null gives a sample-specific p-value:

```
p_value = (1 + n_perm_with_delta_higher) / (1 + B)
```

Bonferroni-correct across n_samples. Threshold p_adj < 0.05 →
recombinant call.

### 2.3 Recombinant classes

After thresholding, classify the changepoint pattern:

- **Simple recombinant** (1 changepoint): `INV INV ... INV REF REF`
- **Double crossover** (2 changepoints): `INV INV REF REF INV INV`
- **Complex mosaic** (≥3 changepoints): mixed pattern

Use BIC-based segmentation (PELT or binary segmentation) to count.

## 3. Atlas-side surface

The atlas already has the **dosage heatmap** (turn 94, FIG_C08-style).
This spec adds:

### 3.1 Recombinant call overlay on the heatmap

When `state.recombinantCalls` is loaded (from a LANTA-side TSV), each
sample's row in the heatmap gets:

- Vertical tick marks at detected changepoint positions
- Side annotation: REF→INV, INV→REF, double-crossover, etc.
- Filter: show only recombinants (hide non-recombinant rows)

### 3.2 Recombinant tab in G-panel

Tab `recombinants` in the G-panel showing:

- Per-sample changepoint position(s) for the active candidate
- Recombinant class
- Confidence (p_adj)
- Action: confirm / dismiss / inspect (jump to that sample's row)

### 3.3 Per-candidate breeding implication

In the breeding-readiness card (Atlas 5), a column flags recombinants:

> *"Of the {n_het} HET samples, {n_recom} are recombinant (simple:
> {n_simple}, double-crossover: {n_dco}, complex: {n_complex}). True
> HET count for breeding purposes: {n_het - n_recom}."*

## 4. Implementation slices

### Slice 1 — JSON layer schema (~0.3 turn, atlas-side)
- Define the `recombinant_calls_v1` JSON layer schema
- Per-candidate, per-sample: `{ sample, changepoints: [pos], class, p_adj }`
- Detect on JSON load via `detectSchemaAndLayers`

### Slice 2 — heatmap overlay (~0.5 turn)
- Extend `drawDosageHeatmap` to render changepoint ticks
- Per-sample side annotation column

### Slice 3 — G-panel recombinants tab (~0.5 turn)
- Tab visibility gated on layer presence
- Per-row actions

### Slice 4 — LANTA-side producer (~outside atlas)
- Run `STEP_C01h_recombinant_scanner.R` on Clair3 dosage matrices per
  candidate → emit `recombinant_calls_v1.json`

## 5. Open questions

1. **Tier interaction**: a candidate where most "HET" samples are
   recombinants — does that downgrade the tier? Probably yes, because
   true HET count is what matters for breeding usability.
2. **GHSL as confirmation**: the LANTA spec uses GHSL karyotype runs
   as a tier-3 confirmation layer. Atlas-side, we can show both
   dosage-changepoint AND GHSL-run-split in the same view as
   independent evidence.
3. **PCA demotion**: agreed in the LANTA spec. Atlas should NOT
   surface PCA-based recombinant calls anymore. The existing per-window
   class track can stay as a diagnostic but never feed the recombinant
   classification.

## 6. What this is NOT

- Not a replacement for the karyotype assignment itself. Recombinants
  are a sub-class of HET that needs separate handling.
- Not a way to discover *new* inversions. It refines the genotype
  call within an already-discovered candidate.
- Not implementable atlas-side without a LANTA-side producer. The
  changepoint scan needs raw dosage matrices that are too big for the
  browser.

## 7. Tests

- Synthetic candidate with 30 HET samples, 3 of which have planted
  changepoints at fixed positions → producer detects all 3, p_adj < 0.05.
- Dosage heatmap overlay renders ticks at correct positions.
- G-panel recombinants tab populates from synthetic JSON.

## 8. Dependencies

- LANTA-side `STEP_C01h_recombinant_scanner.R` (already exists per
  chat `b7bc2608`).
- Atlas dosage heatmap (shipped turn 94).
- G-panel scaffold (`SPEC_g_panel_unified_groups.md` Slice 1).

## 9. Cross-references

- `DECOMPOSE_UPGRADE_DESIGN_chat9.md` (LANTA-side, may exist as a
  standalone file in Quentin's project tree).
- `SPEC_per_candidate_breeding_readiness_card.md` for the breeding
  flag column.
- `SPEC_review_surfaces_auto_and_lineages.md` for the G-panel pattern.
