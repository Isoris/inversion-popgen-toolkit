# S2 — Indel slope / burden layer

**Status:** spec only. No code this turn.
**Used by:** P5.1 Section 4 (edge tracks) — as one of the layers,
labelled "indel-burden support" not "breakpoint proof".
**Source scripts:** `cheat9_indel_slope.R`,
`STEP_P04_export_weak_indel_candidates.py` (you uploaded these earlier).

---

## What this layer is

Per (candidate × sample) burden of indels, broken into:

- heterozygous indels per kb
- homozygous-alt indels per kb
- het:hom ratio

Computed from the cohort's existing indel call set (BCF/VCF) without
needing BAMs. The signal is closer to **het dosage / GHSL** than to
PCA/MDS — it tells you *which samples carry structural variation*,
not *where the breakpoint is*.

## Interpretation rules (your discipline notes)

| H_class | expected indel-burden pattern |
|---|---|
| HOM_REF (H1/H1) | low heterozygous, lower hom-alt burden |
| HET (H1/H2) | many heterozygous indels (the structural carriers) |
| HOM_ALT/HOM_INV (H2/H2) | more homozygous-alt indels |

If a candidate has 226 samples and the indel-burden pattern across
H1/H1 / H1/H2 / H2/H2 follows the expectation above, that's
**indel-burden support** for the karyotype call — NOT proof of
where the boundary is.

## Atlas placement

Boundaries + SV evidence page (P5.1), Section 4 (edge-evidence
tracks). The track is labelled:

> **Indel-burden support** *(not breakpoint proof — see spec S2)*

…with a hover tooltip explaining the discipline.

A second placement: Section 5 (dosage heatmap) can sort/colour
samples by indel burden as an alternative to dosage. This is opt-in
via a heatmap-mode toggle.

## Data layer schema

Two TSVs under `${base}/indel_burden/`:

### `indel_burden_per_sample.tsv.gz`

One row per (candidate_id, sample_id). The "burden per sample".

| column | type | description |
|---|---|---|
| `candidate_id` | string | candidate this row belongs to |
| `sample_id` | string | sample (matches `samples.ind`) |
| `n_het_indels` | int | heterozygous indels in candidate region |
| `n_homalt_indels` | int | homozygous-alt indels in candidate region |
| `n_homref_indels` | int | homozygous-ref indels in candidate region |
| `n_missing` | int | indels with missing genotype |
| `region_kb` | float | candidate region width in kb |
| `het_per_kb` | float | `n_het_indels / region_kb` |
| `homalt_per_kb` | float | `n_homalt_indels / region_kb` |
| `het_to_hom_ratio` | float | `n_het_indels / max(1, n_homalt_indels)` |
| `H_class` | enum | inherited from candidate_group_membership for convenience |

### `indel_burden_sliding.tsv.gz` (optional, for boundary support)

If the upstream pipeline emits sliding-window per-sample indel
burden, the atlas can use it for **boundary transition support**
in P5.1 Section 4.

| column | type | description |
|---|---|---|
| `candidate_id` | string | |
| `chrom` | string | |
| `window_start` | int | bp |
| `window_end` | int | bp |
| `sample_id` | string | |
| `het_per_kb` | float | window-local |
| `homalt_per_kb` | float | window-local |

Reading: per-sample sliding tracks aligned to candidate coords.
Useful for showing the boundary transition where
heterozygous-indel slope changes sharply.

## Layer JSON shape (for atlas-side ingestion)

When loaded as a JSON enrichment, the atlas should expect:

```json
{
  "indel_burden": {
    "per_sample": [
      {
        "candidate_id": "lg28-shelf",
        "sample_id": "Cgar_001",
        "n_het_indels": 47,
        "n_homalt_indels": 3,
        "het_per_kb": 0.018,
        "homalt_per_kb": 0.0012,
        "het_to_hom_ratio": 15.7,
        "H_class": "H1/H2"
      }
    ],
    "sliding": null,
    "schema_version": 1
  }
}
```

`sliding` is null when the per-window track isn't available.

## Computed views

### View 1 — burden distribution by H_class

For an active candidate, plot `het_per_kb` and `homalt_per_kb`
distributions, three boxplots side by side (H1/H1, H1/H2, H2/H2).
Expected pattern: H1/H2 has highest het, H2/H2 has highest hom-alt,
H1/H1 has lowest both.

If pattern matches: chip "indel-burden support: ✓".
If pattern violates: chip "indel-burden support: ✗ — burden does not
track the karyotype call".

### View 2 — within-H_class outliers

Samples whose burden departs from their H_class's expected range —
candidates for re-genotyping or contamination check. Render as a
sortable table or a scatter (sample_idx vs burden, coloured by
H_class).

### View 3 — boundary transition (sliding only)

If `indel_burden.sliding` is loaded, the layer can render an edge
slope: how sharply does het_per_kb jump at the candidate boundaries?

## Server endpoint (future)

```
GET /api/indel_burden/per_sample?candidate_id=<id>
GET /api/indel_burden/sliding?chrom=&start=&end=&cap=200
```

These don't exist yet. The atlas can also accept the layer as a
JSON enrichment (drag-drop) without server.

## Spec is closed

This layer is **deliberately not implemented this turn** because:

1. The two TSV files don't exist yet — they require running
   `cheat9_indel_slope.R` against the cohort BCF/VCF.
2. The cluster R-script needs to emit those TSVs for the canonical
   set of candidates first.
3. The atlas-side renderer for "burden boxplots by H_class" is a
   straightforward view but should land alongside Section 4 of the
   Boundaries + SV page (P5.1) — not before.

When the upstream R script ships the TSVs, write a follow-up patch
that:

- Adds `_isIndelBurdenJSON()` schema detector.
- Adds the boxplot renderer to P5.1 Section 4.
- Optionally adds a server endpoint for live serving.

## Risk notes

- **Don't claim breakpoint proof.** Every label, tooltip, and
  manuscript reference for this layer must say "indel-burden
  support" or "indel-burden does/doesn't track the karyotype". Never
  "indel marks the breakpoint at X".
- **Family confounding.** Indel burden can be inflated by hatchery
  family effects (high relatedness within a family → shared
  artifact indels). The boxplot view should be split by family or
  use family-corrected residuals before claiming H_class
  enrichment. Spec note for the renderer.
- **Region size effect.** Don't compare raw `n_het_indels` between
  candidates of different widths. Always per-kb.
