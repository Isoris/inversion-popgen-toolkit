# S1 — SV evidence tables (5-table schema)

**Status:** spec only. No code this turn.
**Used by:** P5.1 Section 6 (SV evidence) and Section 9 (group association).

---

## Why tables first

> Use tables first, visualization second.

Visualizations come and go. Tables are the canonical artifact —
they're what gets exported, reviewed, and ends up in Supplementary
Tables of the manuscript. Build the data layer right, then layer
visualizations on top.

## The 5 tables

### Table 1 — `sv_variant_catalog.tsv`

One row per SV record observed in the cohort. The "what" of each SV.

| column | type | description |
|---|---|---|
| `sv_id` | string | unique ID, e.g. `lg28-inv-0007` |
| `chrom` | string | chromosome name (matches reference) |
| `pos_start` | int | 1-based start coordinate |
| `pos_end` | int | end coordinate; for BND, may be the mate position |
| `sv_type` | enum | `INV`, `BND`, `DEL`, `DUP`, `INS`, `CNV`, `OTHER` |
| `caller` | string | `manta`, `delly`, `lumpy`, `nanovar`, etc. |
| `quality` | float | caller-provided QUAL or a quality proxy |
| `candidate_id` | string \| null | ID of the cohort candidate this SV is associated with (nullable for SVs not yet linked) |
| `zone` | enum | `left_boundary`, `right_boundary`, `body`, `left_flank`, `right_flank` (relative to candidate) |
| `distance_to_edge` | int | bp distance to the nearest candidate boundary (signed, negative = inside body) |
| `event_id` | string \| null | for BND mate-pair grouping (Manta `EVENT`, MATEID, PARID) |
| `mate_sv_id` | string \| null | if this is a BND with a paired mate, the mate's `sv_id` |

### Table 2 — `sv_sample_genotypes.tsv`

One row per (sv_id, sample_id). The "who has this SV".

| column | type | description |
|---|---|---|
| `sv_id` | string | foreign key → Table 1 |
| `sample_id` | string | foreign key → Table 3 (`candidate_group_membership.tsv`) |
| `sv_gt` | enum | `0/0`, `0/1`, `1/1`, `./.`, or specific BND genotype |
| `sv_dosage` | int \| null | 0/1/2 alt-allele dosage where defined |
| `present` | bool | true if sample carries any alt evidence (covers BND presence/absence too) |

### Table 3 — `candidate_group_membership.tsv`

One row per (candidate_id, sample_id). Defines what "group" each
sample belongs to for a given candidate. This drives the H_class /
family / ancestry stratification.

| column | type | description |
|---|---|---|
| `candidate_id` | string | foreign key |
| `sample_id` | string | foreign key |
| `H_system` | string | which haplotype system this candidate is on, e.g. `LG28-shelf` |
| `H_class` | enum | `H1/H1`, `H1/H2`, `H2/H2`, `ambig`, `excluded` |
| `family` | string | hatchery family / cluster ID (from NAToRA / family JSON) |
| `ancestry_group` | string | NGSadmix K cluster label (e.g. `K8-c3`) |

### Table 4 — `sv_group_genotype_counts.tsv`

One row per (candidate_id, sv_id, H_class). Aggregate counts for
each SV stratified by structural haplotype class.

| column | type | description |
|---|---|---|
| `candidate_id` | string | |
| `sv_id` | string | |
| `H_class` | enum | `H1/H1`, `H1/H2`, `H2/H2` |
| `n_00` | int | count of samples with sv_gt 0/0 in this H_class |
| `n_01` | int | count of 0/1 |
| `n_11` | int | count of 1/1 |
| `n_missing` | int | count of ./. or unobserved |
| `freq_alt` | float | observed alt allele frequency in this H_class |
| `present_rate` | float | fraction of samples in this H_class with `present=true` |

This table is **derived** from Tables 2 + 3 — recomputable any time
either source changes. Pipeline: join sv_sample_genotypes with
candidate_group_membership on sample_id, group by (candidate_id,
sv_id, H_class), aggregate counts.

### Table 5 — `sv_group_enrichment_results.tsv`

One row per (candidate_id, sv_id, test). The statistical question:
does this SV's genotype distribution differ by H_class?

| column | type | description |
|---|---|---|
| `candidate_id` | string | |
| `sv_id` | string | |
| `zone` | enum | inherited from Table 1 |
| `sv_type` | enum | inherited from Table 1 |
| `test` | enum | `fisher_exact`, `chisq`, `cochran_armitage`, `binom_glm` |
| `group_a` | string | e.g. `H1/H1` |
| `group_b` | string | e.g. `H2/H2` |
| `odds_ratio` | float | from 2×2 contingency on (group_a, group_b) × (present, absent) |
| `p_value` | float | raw p |
| `fdr` | float | Benjamini-Hochberg adjusted across all SVs in this candidate |
| `pattern_label` | enum | see vocabulary below |

## Pattern vocabulary

The `pattern_label` column is the human-readable interpretation of
the group enrichment result.

| label | when it fires | what it means |
|---|---|---|
| `canonical_breakpoint_marker` | SV in left_boundary or right_boundary; H1/H1 ≈ 0/0, H2/H2 ≈ 1/1, H1/H2 ≈ 0/1 | the SV alt allele tracks the structural haplotype dose 1:1 — could be the breakpoint itself or a perfectly linked variant |
| `dominant_or_presence_marker` | SV present (any alt) in H1/H2 + H2/H2, absent in H1/H1 | dominant (or near-dominant) inheritance pattern; useful as PCR presence/absence assay |
| `subhaplotype_marker` | SV present only in some subset of H2/H2 carriers | identifies a subhaplotype within H2 — can split H2 into H2.1, H2.2 |
| `het_specific_marker` | SV strongly enriched in H1/H2 but rare in both homozygous classes | rare; possibly a recombinant-junction marker |
| `internal_linked_marker` | SV in body zone (not boundary), still associated with H_class | linked cargo within the inversion; not breakpoint proof |
| `uninformative` | enrichment p > 0.05 or odds ratio ~1 | SV doesn't differentiate the haplotype classes |
| `family_artifact` | SV strongly correlated with a family ID, weakly with H_class | likely PCR/library artifact within a hatchery family |
| `repeat_noise` | SV in repeat-flagged region or with high MAPQ0 nearby | likely false positive from repeat-mediated alignment errors |

The atlas should never display "this is a breakpoint" without a
pattern label that supports it. `canonical_breakpoint_marker` only
fires on boundary-zone SVs with the strict 0/1/2 dose pattern.

## Important reminder (from your message)

> H classes are structural-haplotype states, not guaranteed AA/AB/BB
> genotypes. Do not force H1=A and H2=B. Ask whether SV genotype
> distribution differs by H class and what pattern it follows.

So the tests in Table 5 ask "is genotype distribution different
between H_class groups?" — they don't presume "H1 carries the alt
allele". The pattern label is what tells you which way the
asymmetry goes.

## File layout on disk

```
${base}/sv_evidence/
├── sv_variant_catalog.tsv.gz
├── sv_sample_genotypes.tsv.gz
├── candidate_group_membership.tsv.gz   # written by atlas, not from upstream
├── sv_group_genotype_counts.tsv.gz     # derived; can be re-derived live
└── sv_group_enrichment_results.tsv.gz  # derived; can be re-derived live
```

`candidate_group_membership` is special: it's written by the atlas
when the user finalizes a candidate's `locked_labels` (i.e. it's
*output* of the atlas workflow, not input). The atlas server can
read this back when ranking samples for the dosage heatmap row
order.

## Atlas-side ingestion plan

The atlas needs an ingestion path that:

1. Reads `sv_variant_catalog.tsv.gz` once on data load.
2. Reads `sv_sample_genotypes.tsv.gz` lazily, by sv_id range or by
   active candidate.
3. Reads or computes `sv_group_genotype_counts` per active
   candidate.
4. Reads or computes `sv_group_enrichment_results` per active
   candidate (FDR scoped to within-candidate, not global).

This belongs in a future patch (call it P5.2 when ready).

## Server-side endpoints (future)

```
GET /api/sv/catalog?chrom=&start=&end=&candidate_id=
GET /api/sv/genotypes?sv_id=
GET /api/sv/enrichment?candidate_id=
```

These don't exist yet. Spec only.

## Source-of-truth question

If both the cluster pipeline AND the atlas can compute Tables 4-5,
which wins?

**Answer:** the cluster pipeline writes them once in a snapshot.
The atlas can recompute on the fly when:
- the snapshot is missing
- the user manually changes group memberships (i.e. retriggers
  candidate_group_membership)

When both exist, the atlas reads the cluster's version by default
and shows a "stale" warning if the candidate_group_membership has
been edited since.
