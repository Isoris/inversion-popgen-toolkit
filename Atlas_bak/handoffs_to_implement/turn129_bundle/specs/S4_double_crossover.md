# S4 — Double crossover / recombinant extension

**Status:** spec only. Future module — explicitly do not implement
this turn.
**Used by:** an extension of S3's read-evidence framework, applied
to internal transitions instead of just left/right boundaries.

---

## What this module is

Most of the atlas's machinery treats inversions as
**simple bipartite events**: one left boundary + one right boundary,
samples are H1/H1 or H1/H2 or H2/H2 across the whole interval.

But real cohorts contain **recombinants**: samples whose haplotype
state changes mid-interval. A classic case is double crossover:

```
position →
    L_boundary           internal_1   internal_2          R_boundary
        |                    |             |                  |
        v                    v             v                  v
H1 ─────[───────H2 segment──[─H1 segment──]─H2 segment───────]─────── H1
```

The sample carries a sub-segment that "reverted" to H1 inside an
otherwise-H2 inversion. This is a **recombinant structural haplotype**.

## Why it matters for the manuscript

Recombinants are evidence that:
1. The inversion is old enough to have accumulated meiotic
   recombinants between H1 and H2 (constrains age / Ne).
2. There's at least some gene flow between H_classes (otherwise
   recombinants couldn't form).
3. The break point itself is well-defined — recombinants have
   sharp transitions inside the inversion body.

Detecting them requires the same machinery as S3, but applied to
**paired internal transitions** rather than just edges.

## The discipline

> Need paired internal transitions:
>   entry into alternate segment
>   exit from alternate segment
>
> Do not call double crossover from one edge alone.

That last line is the rule. A single transition inside the body is
**not** a double crossover. It's either:
- A single crossover (recombinant where the inversion ends mid-body)
- A genotyping error
- A real but unpaired structural feature
- Noise

You only get to call "double crossover" when you observe **both** an
entry and an exit within the body of the inversion, in the same
sample, with consistent H_class flips.

## Algorithm shape

Per candidate, per sample, scan the candidate body in sliding
windows for **H_class transitions**:

1. Compute per-window H_class evidence (from PCA, dosage, GHSL,
   indel burden, or read evidence — whatever's available).
2. Detect transitions: positions where the per-window H_class call
   flips from one state to another at a sharp threshold.
3. For each pair of detected transitions in the same sample within
   the candidate body:
   - If the pair forms entry-then-exit pattern (e.g. H2 → H1 →
     H2), label `paired_internal_transition`.
   - Apply S3's beta-binomial scoring to each transition
     individually for confidence.
4. Roll up to candidate level.

## Pattern labels

The vocabulary you defined:

| label | when it fires |
|---|---|
| `double_crossover_supported` | Two paired internal transitions in same sample with consistent H_class flips and read-evidence support at both transitions |
| `recombinant_structural_haplotype` | Same as above, but observed in multiple samples with consistent breakpoints — suggests a stable recombinant haplotype circulating in the population |
| `paired_internal_transition` | Two paired transitions detected, but read-evidence support is below the threshold — interesting candidate, not yet supported |
| `single_internal_edge_only` | One internal transition only — could be partial inversion, genotyping error, or noise. Do not call as double crossover. |
| `mosaic_without_read_support` | Multiple transitions with no read-evidence support across any of them — likely genotyping noise or somatic mosaicism (rare) |
| `nested_or_overlapping_system` | Transitions look paired but H_class flips are inconsistent (e.g. H1→H2→H3) — suggests a second inversion overlapping this one |
| `mapq0_repeat_ambiguity` | Transitions colocate with MAPQ0 / repeat — exclude |
| `family_artifact` | Transitions correlate with hatchery family rather than being independent observations |
| `unresolved_complex` | Transitions exist but don't fit any of the above — fallback bucket |

## Inputs (shared with S3)

- Per-sample, per-window H_class evidence (the same matrix used in
  Sections 4-5 of P5.1).
- Read-evidence clusters within candidate bodies (BAM-derived,
  per-sample).
- `candidate_group_membership.tsv` for nominal H_class assignments.

## Outputs

`recombinant_events.tsv`:

| column | type | description |
|---|---|---|
| `candidate_id` | string | parent candidate |
| `sample_id` | string | sample carrying the event |
| `transition_1_bp` | int | position of first transition |
| `transition_2_bp` | int | position of second transition |
| `width_bp` | int | distance between transitions |
| `H_class_outer` | enum | H_class outside the recombinant segment |
| `H_class_inner` | enum | H_class inside the recombinant segment |
| `evidence_score_1` | float | S3 posterior at transition 1 |
| `evidence_score_2` | float | S3 posterior at transition 2 |
| `pattern_label` | enum | from vocabulary above |

Plus a roll-up `recombinant_summary_per_candidate.tsv` for the
catalogue (P6.1):

| column | description |
|---|---|
| `candidate_id` | |
| `n_paired_internal_transitions` | how many `paired_internal_transition` events found |
| `n_double_crossover_supported` | how many fully supported |
| `recombinant_haplotype_count` | how many distinct recombinant H states observed |
| `recombination_class` | summary label for the candidate as a whole: `simple_inversion` / `recombinant_haplotype_population` / `nested_system` / `unresolved_complex` |

## Atlas placement

P5.1 Section 8 (read-evidence clusters) becomes the home for both
S3 (boundary breakpoint scoring) and S4 (internal transition
scoring). They share the underlying read-evidence framework.

The candidate catalogue (P6.1) gets two new columns:
- `recombinant_summary` — the roll-up label
- `n_recombinants` — count

These appear in the Catalogue table once this layer ships.

## Why this is spec-only

Same reasons as S3 plus:

1. Detection of paired internal transitions requires confidence in
   the per-window H_class calls — which currently come from PCA /
   dosage / GHSL, all of which are noisy mid-interval.
2. The discipline rule ("two paired transitions, not one edge")
   needs validation against a known recombinant. Without ground
   truth, false-positive rate is unknown.
3. The classification of `recombinant_haplotype_population` (a
   stable circulating recombinant) requires multiple samples with
   matching breakpoints — needs cohort-level analysis above the
   per-candidate level.

When this ships, the renderer in P5.1 Section 8 should:

- Show per-sample transition maps (small genome browser strip per
  sample within the candidate)
- Highlight paired transitions in green, single in grey
- Display the `recombination_class` chip on the candidate header
- Link to a recombinant-detail panel showing the underlying read
  evidence

## Risk notes

- **Don't claim recombinants from PCA alone.** Mid-interval H_class
  flips in PCA can come from local LD breakdown, not from real
  recombinants. Require S3-level read evidence for the
  `_supported` labels.
- **Beware confounded signals.** A second overlapping inversion
  (nested system) produces the same H_class flip pattern as a
  recombinant. The `nested_or_overlapping_system` label exists
  precisely for this case — be willing to use it.
- **Population frequency matters.** A single sample with paired
  internal transitions could be one recombinant individual. The
  same transitions in many samples suggest a stable recombinant
  haplotype circulating in the cohort — the `recombinant_haplotype_population`
  label captures this distinction.
