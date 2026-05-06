# S7 — Karyotype-stratified breakpoint & internal-transition evidence scoring

**Status:** spec only. No code this turn.
**Consolidates:** two user prompts ("boundary breakpoint scoring" + "internal transition / double crossover scoring").
**Relationship to existing specs:**
- **Extends** S3 (Bayesian breakpoint scoring) — S3 defines the model; S7 defines how it composes with the read-extraction layer that already exists.
- **Extends** S4 (double crossover) — S4 defines the discipline; S7 defines paired-edge scoring as a single entry/exit object, not two correlated singletons.
- **Builds on** S1 (5-table SV schema) for `candidate_group_membership.tsv`.
- **Honors** S5 (interpretation rules) — every `_supported` label requires Tier-1 evidence per S5 Rule 1/2.
- **Fed by** existing scripts in `inversion_modules/phase_8_evidence_biology/` — see "Reuse, don't reinvent" below.

**Used by:** P5.1 Section 8 (read-evidence clusters) when this module ships.

---

## TL;DR

The user's prompt asked for a new `karyotype_breakpoint_evidence.py` script that scans markdup BAMs, scores breakpoint signals by H_class, and writes atlas-ready JSON. Plus a sibling `internal_transition_evidence.py`.

Reading the project shows this would **duplicate ~80%** of `STEP_03_per_read_evidence.py` (which already does pysam-based BAM scanning, persists per-read evidence with read names, and writes 21 of the 25 Q7B keys). The Bayesian scoring layer, the MAPQ0 evidence stream, the cluster object, the karyotype-conditional posteriors, and the internal-transition pairing **are genuinely new** and worth building.

Proposal: ship two **new** scripts that consume STEP_03's output, plus targeted upgrades to STEP_03 to fill gaps. Boundary and internal transitions share the scoring core because they're the same statistical question at different positions.

```
existing                                  new (this spec)
────────                                  ────────────────
phase_6/run_pipeline.sh
  → consensus boundary TSVs
  → STEP_03B bridge (writes registry)
                                          STEP_E01 cluster builder
phase_8/STEP_03_per_read_evidence.py  ─→  STEP_E02 coverage normalizer
  → evidence_reads_bam.tsv.gz             STEP_E03 karyotype-stratified
  → evidence_reads_vcf.tsv.gz                      Bayesian scorer
                                          STEP_E04 internal-transition
  + (UPGRADE) MAPQ0 emitter                         pairing scorer
  + (UPGRADE) cov_local emitter           STEP_E05 atlas JSON emitter
                                                   for P5.1 §8
```

Five new scripts in a new directory `inversion_modules/phase_8_evidence_biology/q7b_karyotype_evidence/`. Two upgrades to STEP_03. One atlas JSON layer. Total: a few hundred lines of new code, not a few thousand.

---

## Reuse, don't reinvent

The user prompt proposed `karyotype_breakpoint_evidence.py` with these responsibilities:

| Prompt asked for | Where it already lives | What we actually need to add |
|---|---|---|
| Load candidate / group / BAM / SV inputs | `STEP_03B_bp_pipeline_bridge.py` (registry block) + `breakpoint_evidence_audit.py` (VCF counts) | nothing — read the registry block + ledger |
| Define boundary zones | `phase_6_breakpoint_refinement` (KDE + bootstrap CIs) | nothing — read `q3_final_left_bp` / `q3_final_right_bp` from registry |
| pysam scan: split / clipped / discordant mate | `STEP_03_per_read_evidence.py` (`evidence_reads_bam.tsv.gz`) | nothing for these three |
| pysam scan: MAPQ0 enrichment | **gap** — STEP_03 only emits MAPQ for kept reads, not MAPQ0 density | **upgrade to STEP_03**: emit per-window MAPQ0 read count + MAPQ0 fraction |
| Local coverage normalization | **gap** — STEP_03 emits row-level reads, not denominator | **upgrade to STEP_03**: emit `cov_local` per (sample × bp_side × window) |
| Cluster co-located evidence | **gap** — STEP_03 is row-per-read, no clustering | **new STEP_E01** |
| Karyotype-stratified Bayesian score | **gap** — neither STEP_03 nor `breakpoint_evidence_audit` does posterior modeling | **new STEP_E03** |
| Pattern-label classification | partial — `STEP_04_assign_structural_class.py` does global tier; not cluster-level | **new in STEP_E03** |
| Internal-transition pairing | **gap** — entirely new | **new STEP_E04** |
| Atlas JSON emit | partial — registry blocks exist; atlas doesn't read them yet | **new STEP_E05** |

Net new code is ~5 small scripts. The pysam-heavy pieces are already production-tested; we don't rewrite them.

### Why two scripts (not one merged)

Boundary scoring and internal-transition scoring share the **scoring core** but differ in:

- **Where they look** — boundary windows (one per side, ~Δ kb) vs internal transition windows (one per detected segment edge, often dozens per candidate).
- **How they pair** — boundary scoring is per-cluster; internal-transition scoring is per-pair (entry + exit).
- **What they classify** — boundary labels vs recombinant labels.

Sharing the scoring core via a `q7b_score.py` helper module + two thin driver scripts (`STEP_E03` and `STEP_E04`) is cleaner than one monolith.

---

## Critical groundwork: vocabulary alignment

The project has **three different naming conventions for the same thing**:

| Where | Names |
|---|---|
| Cluster pipeline (`STEP_03`, registry `q7b_*`) | `INV / HET / REF` |
| Atlas runtime (`_buildSampleLookups`) | `HOMO_1 / HET / HOMO_2` |
| Spec S1 (`candidate_group_membership.tsv`) | `H1/H1, H1/H2, H2/H2` |

This module's TSV outputs **must pick one**. Proposal: emit **both** the cluster-side (`INV/HET/REF`) and the spec-S1 (`H1/H1, H1/H2, H2/H2`) names in adjacent columns, and document the mapping in a `q7b_group_naming.md` glossary. The atlas reads whichever it wants. Neither side has to refactor.

The mapping is not symmetric across all candidates — `INV` could mean "carries the inverted arrangement" (cluster convention) or "H2/H2" (spec S1) depending on which arrangement is the reference state. The candidate registry block must record `H1_polarity` (which of `INV`/`REF` corresponds to `H1/H1`) at write time so the mapping is unambiguous downstream.

This is a small but real source of bugs if left implicit. Solve once, here.

---

## Inputs (read, don't redefine)

### From the registry (already exists)

Per candidate, via `read_block(cand_id, "per_read_evidence_summary")`:

- `q3_chrom`, `q3_final_left_bp`, `q3_final_right_bp` — boundaries (from STEP_03B bridge)
- `q7b_bam_ledger_path` — path to `evidence_reads_bam.tsv.gz`
- `q7b_vcf_sidecar_path` — path to `evidence_reads_vcf.tsv.gz`
- `q7b_window_bp`, `q7b_min_mapq` — extraction parameters

### From `evidence_reads_bam.tsv.gz` (already exists, 29 columns)

Most relevant columns for this module (from STEP_03 docstring lines 74–116):

```
candidate_id, chrom, sample_id, group ∈ {INV, HET, REF, unknown},
bp_side ∈ {left, right}, bp_pos,
read_name, evidence_type ∈ {discordant_FF, discordant_RR, split_read, soft_clip, mate_at_partner},
clip_end_template ∈ {5prime, 3prime, both, none},
clip_end_mapped ∈ {left, right, both, none},
clip_len_left, clip_len_right, read_pos, read_end, read_strand,
cigar, mate_chrom, mate_pos, mate_strand, mapq,
nm, as_score, xs_score, sa, mc, xa,
is_supplementary, is_primary_of_pair, flag
```

This file is the canonical row-level audit trail. Don't re-extract; read it.

### From `evidence_reads_vcf.tsv.gz` (already exists)

Per (candidate, sample, side, caller, support_type): aggregate counts (DV / RV / PR_alt / SR_alt) and genotype. Used to flag `caller_concordant` clusters.

### From `candidate_group_membership.tsv` (S1 schema; aspirational)

Per (candidate, sample): `H_system`, `H_class`, `family`, `ancestry_group`. Until S1 lands, use the cluster-side `--samples_tsv` (which already has `sample_id, group` per candidate from Q1 decomposition output) as a fallback, and **fail loudly if `H_class` and `group` disagree** — that's a flag, not a silent override.

### NEW upgrades to STEP_03 (this spec requires)

#### Upgrade 1 — `cov_local` per sample per side

Add an emitter that does a separate pysam pileup over `[bp_pos − Δ, bp_pos + Δ]` (Δ = `q7b_window_bp`, currently 500 by default in STEP_03), counts reads passing `--min-mapq` divided by window width, and writes:

```
cov_local.tsv.gz
   columns: candidate_id, sample_id, bp_side, window_bp, n_reads_pass, depth_mean, depth_capped
   one row per (candidate × sample × bp_side)
```

`depth_capped` = `min(depth_mean, 50)` — the cap is the key guardrail (S3 Risk Note 1). Document the cap value in a registry key `q7b_depth_cap`.

#### Upgrade 2 — MAPQ0 stream

STEP_03 currently filters reads at `--min-mapq` (default 20). MAPQ0 reads never reach the ledger. Add a parallel emitter with no MAPQ filter that counts MAPQ0 reads in the boundary windows:

```
mapq0_density.tsv.gz
   columns: candidate_id, sample_id, bp_side, window_bp, n_mapq0, n_total, mapq0_fraction
```

`mapq0_fraction` = `n_mapq0 / n_total`. Used by the `mapq0_repeat_ambiguity` label (S3, S4).

Both upgrades are ~30 lines each in STEP_03. Net additional pysam passes: 1 (the no-filter MAPQ0 pass; cov_local can ride the existing ledger if MAPQ0s are off-by-one — but a separate pass is cleaner).

---

## STEP_E01 — cluster builder

Read `evidence_reads_bam.tsv.gz` and group reads into clusters.

### Clustering rule

Two reads belong in the same cluster if **any** of:

1. `|read_pos_a − read_pos_b| < ε_pos` (default ε_pos = 200 bp)
2. They share a `read_name` (a split read paired with its supplementary alignment)
3. They share an SA-tag partner coordinate within ε_pos
4. They share a `mate_pos` within ε_pos (discordant pair → same target)

Then add the position constraint: clusters must lie within ±Δ_cluster of the candidate's boundaries (Δ_cluster = `q7b_window_bp`, default 500). Reads outside this window are flagged `out_of_zone` and dropped before clustering.

### Output: `boundary_evidence_clusters.tsv`

| column | type | description |
|---|---|---|
| `candidate_id` | string | |
| `cluster_id` | string | e.g. `LG28_cand_1__bp_left__cl_003` |
| `chrom` | string | |
| `bp_side` | enum | `left`, `right`, `internal_<idx>` (internal = STEP_E04) |
| `pos_anchor_bp` | int | median read_pos in the cluster |
| `partner_chrom` | string \| null | most-common SA / mate chrom (or null) |
| `partner_pos_bp` | int \| null | median SA / mate pos (if consistent across reads) |
| `partner_orient` | enum | `compatible_INV`, `compatible_DEL`, `compatible_DUP`, `compatible_TRA`, `incompatible`, `unknown` |
| `n_reads_total` | int | reads in cluster |
| `n_split` | int | of which `evidence_type == split_read` |
| `n_discord_FF_RR` | int | of which `discordant_FF` or `discordant_RR` |
| `n_softclip` | int | `soft_clip` only (no SA partner) |
| `n_mapq0` | int | from MAPQ0 stream, joined by position |
| `n_samples` | int | distinct sample_ids contributing |
| `dominant_evidence_type` | enum | most-common evidence_type |
| `caller_concordant` | enum | `manta`, `delly`, `both`, `neither` (joined from VCF sidecar by ±100bp) |

### Output: `boundary_sample_support.tsv`

Per (cluster × sample). The user's prompt called for this; it's the audit row.

| column | type | description |
|---|---|---|
| `cluster_id` | string | |
| `sample_id` | string | |
| `group` | enum | from cluster pipeline (INV/HET/REF) |
| `H_class` | enum | from S1 (H1/H1, H1/H2, H2/H2) — null if S1 not yet wired |
| `family` | string | from family JSON, null if not provided |
| `n_split_i` | int | reads in this cluster from this sample, evidence_type==split_read |
| `n_discord_i` | int | discordant_FF or discordant_RR |
| `n_softclip_i` | int | soft_clip with no SA |
| `n_mapq0_i` | int | from MAPQ0 stream |
| `cov_local_capped` | float | from cov_local upgrade |

---

## STEP_E02 — coverage normalizer

Joins the per-sample support table with the cov_local table. Writes `boundary_sample_support_normalized.tsv` with two derived columns:

- `support_total_i` = `n_split_i + n_discord_i + n_softclip_i` (MAPQ0 deliberately separate — it's not "support")
- `cov_denominator_i` = `cov_local_capped × window_bp` (= total expected reads in the window if the cap is hit)

These are the per-sample numerator + denominator that feed the beta-binomial.

---

## STEP_E03 — karyotype-stratified Bayesian scorer

Implements S3's beta-binomial model, applied to clusters from STEP_E01.

### Model (per cluster k, per H_class g)

```
Y_{k,g} = Σ_{i in g}  support_total_i        (numerator)
C_{k,g} = Σ_{i in g}  cov_denominator_i      (denominator)
p_{k,g} | data ~ Beta(α₀ + Y_{k,g}, β₀ + C_{k,g} − Y_{k,g})
```

Default prior: `Beta(α₀=1, β₀=99)`. Configurable via `--prior-alpha` / `--prior-beta`. Stored in registry block as `q7b_score_prior_alpha`, `q7b_score_prior_beta` so reproducibility doesn't depend on default values.

### Why beta-binomial and not negative-binomial / Poisson

- Negative binomial is appropriate for raw count overdispersion when the denominator is implicit. Here the denominator is explicit (`cov_local`).
- Poisson rate models give a rate but don't naturally give a credible interval on a probability bounded in [0, 1].
- Beta-binomial gives a posterior on the **rate of breakpoint-supporting reads per unit coverage**, which is what we actually care about.
- Conjugate update means the posterior is closed-form — no MCMC, no PyMC dependency. One line in scipy.

The cohort isn't large enough (226 samples, often <20 per H_class for rare arrangements) to justify hierarchical modeling. A flat beta prior with weak informativeness is honest about what we know.

### Reported summaries (per cluster × H_class)

```
posterior_mean      = (α₀ + Y) / (α₀ + β₀ + C)
posterior_var       = α·β / ((α+β)² · (α+β+1))      # standard Beta variance
cred_interval_95    = [scipy.stats.beta.ppf(0.025, α', β'),
                       scipy.stats.beta.ppf(0.975, α', β')]
```

### Cross-class summaries (per cluster)

For 3-class (`H1/H1, H1/H2, H2/H2`):

- `P(top > middle > bottom)` from a 10,000-draw Monte Carlo against the three Beta posteriors
- `enrichment_ratio` = `posterior_mean(top) / posterior_mean(bottom)`
- `P(intermediate)` = `P(p_HET ∈ (min(p_HOMO_1, p_HOMO_2), max(p_HOMO_1, p_HOMO_2)))` from the same draws — tests the "het carriers are intermediate" expectation. **Reported, not labeled** — a high `P(intermediate)` is consistent with a real breakpoint marker but doesn't prove one.

For ≠ 3-class candidates, only `enrichment_ratio` and pairwise `P(g_a > g_b)` are reported. Document the limitation.

### Family-aware diagnostics (NOT correction)

Computing posteriors after dropping the largest family is a robustness check, not a primary report. Add columns:

- `posterior_mean_drop_top_family` — recompute after dropping the family with the most cluster reads
- `top_family_share` — fraction of cluster reads from the largest family
- Flag `family_dominant` if `top_family_share > 0.5` AND only one family contributes — fires the `family_artifact` label later

This sidesteps "should we correct for family?" by being honest: we report both numbers and let the label discipline (next section) handle the interpretation. The user already does T9 family jackknife inside C01f for ancestry stats; we mirror that here.

### Pattern label assignment

Following S3 vocabulary (verbatim — keep stable):

| label | rule |
|---|---|
| `breakpoint_resolved_candidate` | One H_class strongly carries cluster (`P(top > bottom) > 0.99`), het is between (`P(intermediate) > 0.95` if 3-class), AND `n_split ≥ 1` AND `caller_concordant ∈ {manta, delly, both}` |
| `breakpoint_supported_interval` | Same posteriors, but `n_split == 0` OR `caller_concordant == neither` (Tier-2/3 evidence per S5) |
| `boundary_support_signal` | H_class enrichment present (`P(top > bottom) > 0.95`) but no single-class dominance |
| `linked_internal_marker` | H_class enrichment present but `pos_anchor_bp` is `> 5 kb` inside the candidate body — cargo, not breakpoint |
| `mapq0_repeat_ambiguity` | `n_mapq0 / (n_reads_total + n_mapq0) > 0.5` — drop the cluster from breakpoint claims; still reported in audit |
| `family_artifact` | `family_dominant == true` AND `posterior_mean_drop_top_family / posterior_mean < 0.3` |
| `unresolved_noise` | none of the above; default fallback |

The label is **a property of the cluster**, not the candidate. A candidate can have one resolved cluster on the left and a `mapq0_repeat_ambiguity` on the right.

### Output: `boundary_group_scores.tsv` and `boundary_cluster_classification.tsv`

`boundary_group_scores.tsv` — per (cluster × H_class):
```
cluster_id, H_class, n_samples_in_class,
Y, C, posterior_mean, posterior_var, cred95_lo, cred95_hi,
posterior_mean_drop_top_family
```

`boundary_cluster_classification.tsv` — per cluster:
```
cluster_id, candidate_id, bp_side, dominant_class, P_top_gt_bottom, enrichment_ratio,
P_intermediate, top_family_share, family_dominant, mapq0_fraction, pattern_label
```

Plus a candidate-level rollup (best label across clusters per side):
`candidate_boundary_summary.tsv`:
```
candidate_id, side, n_clusters, best_label, best_cluster_id, best_posterior_top
```

---

## STEP_E04 — internal-transition pairing scorer

Same scoring core, applied to **paired** internal transitions.

### Inputs

The user's prompt described an `internal_transition_table` with `transition_id`, `transition_bp`, `evidence_source`. **Where does this come from?** Two options:

**Option A — derive transitions from existing per-window H_class evidence** (recommended). The atlas already has per-window H_class evidence inside the candidate body from PCA / dosage / GHSL (the "lines" in the atlas with K=3 colours). Detect transitions as positions where the per-sample H_class call flips between adjacent windows. Persist as a new TSV `internal_transitions.tsv` written by a new step `STEP_E04a_detect_internal_transitions.R` (consumes existing C01f outputs).

**Option B — require the user to provide `internal_transition_table.tsv`**. Unclear where it comes from in their pipeline; punts the detection problem.

Pick A. Wire detection from C01f. Schema:

```
internal_transitions.tsv:
  candidate_id, sample_id, transition_id, transition_bp,
  state_before ∈ {H1_like, H2_like, HET_like, unknown},
  state_after  ∈ {same vocab},
  evidence_source ∈ {dosage, localPCA, GHSL, thetaPi_CUSUM},
  confidence  float [0,1] from the per-window method
```

### Pairing rule (the discipline from S4)

> Need paired internal transitions: entry into alternate segment + exit from alternate segment. Do not call double crossover from one edge alone.

For each sample within each candidate:

1. Sort the sample's transitions by `transition_bp`.
2. Build pairs `(t_a, t_b)` where:
   - both lie inside `[start_bp + ε_inner, end_bp − ε_inner]` (truly internal; ε_inner ~ 50 kb so we don't double-count boundary noise)
   - `state_before(t_a) == state_after(t_b)` (return to same outer state)
   - `state_after(t_a) == state_before(t_b)` (consistent inner state)
   - `state_after(t_a) != state_before(t_a)` (it's actually a flip)
   - `transition_bp(t_b) − transition_bp(t_a) >= ε_pair_min` (default 1 kb — avoid noise)
3. The pair is a **putative recombinant segment** for that sample.

### Scoring

Apply STEP_E03 to each transition independently, treating the transition window as a "boundary" (`bp_side = internal_<idx>`). This gives:

- `posterior(t_a) per recomb_class`
- `posterior(t_b) per recomb_class`

Then aggregate per pair:

```
paired_score = min(P(supported at t_a), P(supported at t_b))   # weakest link
shared_samples = samples with both transitions detected
overlap_with_recomb_set = |shared_samples ∩ recombinant_class_g| / |recombinant_class_g|
```

The "recombinant class" is a sample-level switch class:

| class | meaning |
|---|---|
| `no_switch_H1` | sample stays H1-like across the body |
| `no_switch_H2` | sample stays H2-like across the body |
| `switch_H1_to_H2` | one transition only (single crossover or partial) |
| `switch_H2_to_H1` | same, opposite direction |
| `recombinant_H1_H2_H1` | paired transitions, H1 → H2 → H1 |
| `recombinant_H2_H1_H2` | paired transitions, H2 → H1 → H2 |
| `unknown_mixed` | inconsistent state calls |

Classification operates **on samples**, not on candidates. A candidate has a distribution of samples across these classes.

### Pattern labels (per pair, per candidate)

Following S4 vocabulary verbatim:

| label | rule |
|---|---|
| `double_crossover_supported` | Pair has `paired_score > 0.95` AND `n_shared_samples_with_split_evidence ≥ 2` |
| `recombinant_structural_haplotype` | Same as above, but `n_shared_samples ≥ 5` AND breakpoints align across samples (within ε_pos × 5 = 1 kb) |
| `paired_internal_transition` | Pair detected, `paired_score ∈ [0.5, 0.95]`, no Tier-1 (split-read) at one or both transitions |
| `single_internal_edge_only` | Sample has one transition only |
| `mosaic_without_read_support` | Sample has paired transitions in dosage/PCA but `n_split + n_discord == 0` at both windows |
| `nested_or_overlapping_system` | Pair detected but `state_after(t_a) ∉ {H1_like, H2_like}` (e.g. `H3_like`) — signals a second inversion |
| `mapq0_repeat_ambiguity` | Either transition window has `mapq0_fraction > 0.5` |
| `family_artifact` | Recombinant samples are concentrated in one family (`top_family_share > 0.5`) |
| `unresolved_complex` | Pair detected but doesn't fit any category |

### Outputs

`internal_transition_clusters.tsv`, `internal_transition_sample_support.tsv`, `internal_transition_group_scores.tsv` — same shapes as boundary versions, with `bp_side = internal_<idx>`.

`paired_transition_scores.tsv`:
```
candidate_id, pair_id, transition_id_1, transition_id_2, t1_bp, t2_bp, width_bp,
shared_supporting_samples (int), recomb_class, posterior_t1, posterior_t2,
paired_score, pattern_label
```

`recombinant_summary_per_candidate.tsv` — one row per candidate (S4 schema, kept verbatim):
```
candidate_id, n_paired_internal_transitions, n_double_crossover_supported,
recombinant_haplotype_count, recombination_class
```

where `recombination_class ∈ {simple_inversion, recombinant_haplotype_population, nested_system, unresolved_complex}`.

---

## STEP_E05 — atlas JSON emitter

Reads the TSVs above; writes one JSON layer for the atlas.

### Single combined JSON, NOT two

The user's prompts proposed `candidate_boundary_evidence.json` and `recombinant_transition_evidence.json` as two separate outputs. Reading the existing atlas pages (P5.1 §8): both feed into Section 8 (read-evidence clusters). Splitting into two JSONs creates two index entries, two schema-detection paths, and two atlas-side fetchers for what is conceptually one layer.

Single layer. One JSON per candidate (or one combined per chrom). Schema:

```jsonc
{
  "schema_version": "1.0",
  "candidate_id": "LG28_cand_1",
  "chrom": "C_gar_LG28",
  "boundaries": {
    "left":  { "bp": 15115243, "n_clusters": 3, "best_label": "breakpoint_supported_interval", ... },
    "right": { "bp": 18006104, "n_clusters": 1, "best_label": "boundary_support_signal", ... }
  },
  "clusters": [                                      // boundary + internal in one array
    {
      "cluster_id": "LG28_cand_1__bp_left__cl_001",
      "kind": "boundary",                            // "boundary" | "internal_transition"
      "bp_side": "left",
      "pos_anchor_bp": 15115201,
      "partner_orient": "compatible_INV",
      "n_reads": { "split": 4, "discord": 12, "softclip": 7, "mapq0": 2 },
      "n_samples": 18,
      "caller_concordant": "both",
      "scores_by_H_class": [
        { "H_class": "H1/H1", "n_samples_in_class": 60,
          "posterior_mean": 0.0021, "cred95": [0.0008, 0.0048] },
        { "H_class": "H1/H2", "n_samples_in_class": 106,
          "posterior_mean": 0.0140, "cred95": [0.0098, 0.0193] },
        { "H_class": "H2/H2", "n_samples_in_class": 60,
          "posterior_mean": 0.0287, "cred95": [0.0223, 0.0361] }
      ],
      "P_top_gt_bottom": 0.999,
      "P_intermediate": 0.97,
      "enrichment_ratio": 13.7,
      "family_dominant": false,
      "mapq0_fraction": 0.06,
      "pattern_label": "breakpoint_supported_interval"
    },
    {
      "cluster_id": "LG28_cand_1__internal_002__cl_001",
      "kind": "internal_transition",
      "bp_side": "internal_2",
      ...
    }
  ],
  "paired_transitions": [
    { "pair_id": "LG28_cand_1__pair_01",
      "t1_bp": 16203100, "t2_bp": 16804220,
      "shared_samples": ["S017", "S031", "S082"],
      "recomb_class": "recombinant_H1_H2_H1",
      "paired_score": 0.92,
      "pattern_label": "paired_internal_transition" }
  ],
  "candidate_summary": {
    "best_left_label": "breakpoint_supported_interval",
    "best_right_label": "boundary_support_signal",
    "recombination_class": "simple_inversion",
    "n_paired_internal_transitions": 0,
    "n_double_crossover_supported": 0
  },
  "params": {
    "prior_alpha": 1, "prior_beta": 99, "depth_cap": 50,
    "window_bp": 500, "epsilon_pos_bp": 200,
    "epsilon_pair_min_bp": 1000, "epsilon_inner_bp": 50000
  }
}
```

The atlas reads this as `state.data.boundary_evidence` (new layer in `detectSchemaAndLayers`).

### Atlas rendering (P5.1 §8)

The atlas **does not redesign the page**. It populates the existing Section 8 scaffolding from the JSON:

- Boundary evidence track: clusters as ticks at `pos_anchor_bp`, sized by `n_reads_total`, coloured by `pattern_label`
- MAPQ0 density track: `mapq0_fraction` across the candidate window
- SV-overlap chips: `caller_concordant ∈ {manta, delly, both}` per cluster
- Karyotype-group posterior table: `scores_by_H_class` rendered as a small bar chart with credible-interval whiskers
- Sample × cluster heatmap: rows = samples (ordered by `H_class` then `family`), cols = clusters; cell = `support_total_i / cov_denominator_i`
- Paired transition strip: paired bars connecting `t1_bp` and `t2_bp`, coloured by `pattern_label`

No new page. Section 8 already exists as scaffold (per P5.1).

---

## Decision points (need user confirmation)

1. **Where does this live in the codebase?** Proposal: `inversion_modules/phase_8_evidence_biology/q7b_karyotype_evidence/` (new sibling to `q7_existence_audit/`). Rationale: this is downstream of Q7B BAM extraction and conceptually about the *karyotype-conditional* version of the same evidence — `q7b` is the right neighborhood. Alternative: stick it in `q7_existence_audit/` to keep all Q7B logic together (cluttered but simpler).

2. **MAPQ0 — separate stream or first-class evidence_type?** Proposal: separate stream (`mapq0_density.tsv.gz`), join by position. **Don't add MAPQ0 to STEP_03's `evidence_type` enum** — that would mean MAPQ0 reads enter the per-read ledger and downstream consumers would either have to filter it back out or accidentally count it as breakpoint support (S3 Risk Note 3, S5 Rule 4).

3. **Vocabulary alignment — emit both names or pick one?** Proposal: emit both (`group` for cluster convention, `H_class` for spec S1) in the TSVs, document the polarity at write time. The atlas can pick whichever it wants. The downside: every TSV has redundant columns. Worth it for the cross-pipeline clarity.

4. **Internal transition detection source?** Proposal: derive from existing C01f per-window evidence (Option A above), don't require a new user-provided table. Means writing a small R helper `STEP_E04a_detect_internal_transitions.R` that reads existing dosage / PCA / GHSL outputs and emits `internal_transitions.tsv`.

5. **Coverage cap value — 30× or 50×?** S3 says "30× or 50×". Proposal: **50×** as default (less aggressive, keeps signal from 60× outliers without nuking them). Configurable via `--depth-cap`.

6. **Standalone-first, atlas-later?** User said "standalone first, before atlas integration." Agreed. STEP_E01–E04 are all command-line scripts that produce TSVs. STEP_E05 (atlas JSON) is the final, optional step. Ship E01–E04 with TSV outputs first, validate on LG28 + 2 other resolved candidates, then E05.

---

## Test plan

### STEP_03 upgrades

- `test_cov_local_emits_capped_value` — synthetic high-depth sample (80×) → `depth_capped == 50`.
- `test_mapq0_density_excludes_passing_reads` — reads with MAPQ ≥ 20 not in MAPQ0 stream.

### STEP_E01

- `test_split_with_sa_clusters_with_partner_position` — split read at pos 100 with SA pointing to pos 5000 clusters with a discordant pair landing at pos 5050.
- `test_clip_orientation_drives_partner_compatibility` — soft-clip at CIGAR[0] + mate at later pos → `partner_orient == compatible_INV`.
- `test_out_of_zone_reads_dropped` — read at `bp_pos + 5000` with `q7b_window_bp == 500` → not in cluster.

### STEP_E03

- `test_beta_posterior_matches_scipy` — analytical test: known α₀, β₀, Y, C → known mean, cred95.
- `test_p_top_gt_bottom_monte_carlo` — two non-overlapping posteriors → P > 0.99 in 10k draws.
- `test_p_intermediate_three_class` — synthesize 3-class posteriors with HET in middle → `P_intermediate > 0.95`.
- `test_family_drop_recovers_signal_when_clean` — same data after dropping a family with no contribution → posterior unchanged.
- `test_pattern_label_decision_tree` — known synthetic clusters → expected labels for each rule branch.

### STEP_E04

- `test_pair_detection_requires_consistent_state_flip` — `H1 → H2 → H1` paired; `H1 → H2 → H3` not paired (nested signal).
- `test_pair_min_distance_filter` — two transitions 500 bp apart with `ε_pair_min == 1000` → not paired.
- `test_recombinant_class_assignment` — sample with paired H1→H2→H1 → class `recombinant_H1_H2_H1`.
- `test_double_crossover_supported_requires_split_at_both` — only one transition has Tier-1 → label drops to `paired_internal_transition`.

### STEP_E05

- `test_json_round_trip` — deserialize → assert clusters / paired_transitions arrays present.
- `test_schema_detector_recognizes_boundary_evidence` — atlas-side detector loads JSON.

---

## Assumptions, caveats, and what we deliberately don't do

### Assumptions

- The boundary registry block (`q3_final_left_bp`, `q3_final_right_bp`) is correct. STEP_03B already trusts it; we inherit that trust.
- `evidence_reads_bam.tsv.gz` is the canonical row-level BAM extraction. We don't re-pileup.
- Per-sample H_class assignments exist (from cluster pipeline `samples_tsv` or S1 `candidate_group_membership.tsv`). The module fails loudly if missing, doesn't silently default.
- The 3-class case (H1/H1, H1/H2, H2/H2) is the common one. The code degrades gracefully for K-class candidates but `P_intermediate` is only meaningful for K=3.

### Caveats baked into outputs

- `pattern_label` is per-cluster, not per-candidate. The candidate-level rollup is a max over clusters.
- `breakpoint_supported_interval` does **not** mean "the breakpoint is at `pos_anchor_bp`". The molecular breakpoint still requires assembly junction (see Q4 mechanism, `cheat29_junction_forensics.R`). S5 Rule 1 / 2.
- `recombinant_haplotype_population` is a strong claim. Demanding `n_shared_samples ≥ 5` AND breakpoint alignment within 1 kb is intentionally conservative. Over-claiming this label compromises the manuscript's recombination-population genetics section.
- `family_dominant` is a flag, not a correction. Drop-top-family posterior is reported; the user/atlas decides what to do with it.

### Out of scope

- Phasing reads with PEPPER / WhatsHap. Long-read or trio data not assumed.
- BND mate-pair graph reconstruction beyond `partner_pos_bp` lookup. `STEP_02_bnd_sided_support.py` already does single-sided rescue; we don't redo it.
- Ground-truth calibration of α₀, β₀ priors. Need a known-truth candidate before tuning. LG28 ~17 Mb shelf is the obvious calibration target.
- A new atlas page. P5.1 §8 is the home; this spec only fills its data.
- Replacing `STEP_03_per_read_evidence.py`. We add to it; we don't refactor it.

### Why this is spec-only

- The MAPQ0 + cov_local upgrades to STEP_03 need testing on real LANTA data, not synthetic.
- Internal-transition detection (STEP_E04a) depends on per-window H_class outputs from C01f that may need a small TSV emit step before this module can read them.
- The β prior calibration needs at least one fully resolved candidate validated against assembly junctions (Q4 mechanism work). Without that, "supported" labels are claims about a model, not the genome.

When this lands, P5.1 §8 stops being empty-state.

---

## Implementation order (when this leaves spec-mode)

1. STEP_03 upgrades — `cov_local.tsv.gz`, `mapq0_density.tsv.gz`. Re-run on LANTA for the 226-cohort subset of resolved candidates. (1 session.)
2. `STEP_E04a_detect_internal_transitions.R` reading C01f per-window outputs. (1 session.)
3. STEP_E01 cluster builder. Tests against synthetic + LG28 reads. (1 session.)
4. STEP_E02 + STEP_E03 scoring core, with full test suite. (1–2 sessions.)
5. STEP_E04 pairing scorer. (1 session.)
6. Validate end-to-end on LG28 ~17 Mb shelf, plus 2 other resolved candidates. Calibrate priors. (manual review + 1 session of tuning.)
7. STEP_E05 atlas JSON emit. P5.1 §8 wires the new layer. (1 session.)

Net: ~7 sessions, with assembly-junction validation in the middle as the gate. Don't ship "supported" labels into the manuscript before that gate.
