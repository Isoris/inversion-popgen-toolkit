# BK_KEYS_EXPLAINED.md

Per-key reference for the three canonical phase-4b block schemas
(chat-15 BK). For each key: **where the value comes from** (which
script, which computation on what input), **what it measures about
the genome** (the biological quantity), and **how it is used
downstream** (which consumer reads it, what decision it informs).

Keys are grouped by their emitting block. The block is the Tier-2
JSON evidence object written via `reg$evidence$write_block()`; the
schema's `keys_extracted` directive extracts a flat scalar from the
block and writes it to `<candidate_outdir>/keys.tsv` for aggregate
consumption (C01d scoring, compute_candidate_status, downstream
tables).

---

## Table of contents

1. [Reading this document](#reading-this-document)
2. [Block 1: `internal_dynamics` — PCA-based preliminary class assignment](#block-1-internal_dynamics)
3. [Block 2: `recombinant_map` — per-sample R + G gate recombinant calls](#block-2-recombinant_map)
4. [Block 3: `internal_ancestry_composition` — nested ancestry structure per candidate](#block-3-internal_ancestry_composition)
5. [Pre-rename / post-rename key map](#pre-rename--post-rename-key-map)

---

## Reading this document

**What is a "candidate"?** A candidate is a genomic interval on one
chromosome (chrom, start_bp, end_bp) that the phase-2 inv_detect
pipeline flagged as potentially harbouring an inversion polymorphism.
Every key is scoped per-candidate.

**What is a "sample"?** One of the 226 *Clarias gariepinus* broodstock
individuals in the cohort (verified as pure *C. gariepinus* by Mash
screening against both parental subgenome references). Each sample has
one diploid genotype per candidate.

**What do Q1 / Q2 / Q6 mean?** The key prefix identifies which of
the seven characterisation questions the key informs:
- **Q1** — What is the shape of the signal? Does it look like an inversion?
- **Q2** — What's happening inside the interval? Internal dynamics,
  recombinants, nested ancestry.
- **Q6** — What's the frequency and class composition?

**Per-sample vs cohort-level.** The schemas carry both per-sample
arrays (rich; variable length; not flattened) and cohort-level scalars
(flat; one value per candidate; these are what `keys_extracted`
produces). Downstream scoring reads the scalars; downstream plotting
and the manuscript text consume the arrays directly from the block
JSONs.

---

## Block 1: `internal_dynamics`

**Writer.** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R`

**Inputs.**
- PCA local-window embeddings (per-candidate, ~20 kb windows)
  produced by upstream phase-2 steps (sim_mat / local-PCA).
- Seed sets: *either* flashlight v2 seeds (per-sample anchor classes
  HOM_REF / HET / HOM_INV from independent SV breakpoint evidence)
  OR the step03 cohort-wide deconvolution seed table, combined via
  `lib_step03_seed_loader.R` with drop-on-conflict reconciliation.
- cheat2 het-DEL-at-breakpoint constraint (samples with a het
  deletion signature spanning the inversion breakpoint are forced
  HOM_INV → HET).

**What it does.** For a candidate interval:
1. Loads the per-window × per-sample PC1 coordinate.
2. Seeds a k-means (k=2 or k=3) using the loaded seed sets; chat-9
   design REQUIRES seeding — if seeds are sparse and `seeded_kmeans`
   falls back to unsupervised, the block is stamped `no_seeding`
   and C01i_d_seal skips the candidate (we do NOT ship unsupervised).
3. Computes a preliminary class label per sample from the mean PC1
   across the interval's windows.
4. Quantifies cluster quality (silhouette, BIC gap k=3 vs k=2) and
   per-sample consistency (fraction of windows whose local class
   matches the sample's global class).
5. Integrates HET-support from WhatsHap phase blocks and flags
   samples where the SV anchor disagrees with the PCA cluster.

**What it is NOT.** This block does not call recombinants (that's
`recombinant_map`), does not characterise internal ancestry structure
(that's `internal_ancestry_composition`), and does not write the
final group assignment (that's `STEP_C01i_d_seal`, which synthesises
the three blocks plus GHSL into the `sample_registry` groups table).

### Keys

#### `q2_pca_decomp_status`
- **Comes from.** `block_data$decomp_status` (string). One of
  `seeded_ok` or `no_seeding`.
- **Measures.** Whether the seeded k-means ran to completion with
  enough seeds in every class. This is an execution state, but it
  has biological meaning: `no_seeding` flags candidates where
  flashlight + step03 could not agree on enough anchor samples,
  which itself diagnoses either (a) sparse SV evidence for this
  interval or (b) cohort-wide disagreement between SV-based and
  deconvolution-based class assignments.
- **Used by.** C01i_d_seal gates synthesis on
  `q2_pca_decomp_status == "seeded_ok"`; candidates with
  `no_seeding` are dropped to the synthesis-skipped set rather than
  being synthesised with unsupervised labels (chat-9 design).

#### `q2_pca_decomp_skip_reason`
- **Comes from.** `block_data$decomp_reason`; populated only when
  `decomp_status == "no_seeding"`. String constants include the
  seed_loader verdict ("all samples dropped by conflict filter",
  "insufficient flashlight coverage and no step03 table") or
  `seeded_kmeans fell back to unsupervised — not shipping`.
- **Measures.** Root cause of the skip. Diagnostic only.
- **Used by.** Downstream completion-accounting tables (shows why
  a candidate did not receive a preliminary class assignment).

#### `q2_pca_seed_source`
- **Comes from.** `block_data$seed_source`. String values:
  `flashlight+step03`, `step03`, `flashlight`, or `flashlight_legacy`.
- **Measures.** Which anchor-evidence source produced the seeds used
  to initialise the k-means. `flashlight+step03` means both agreed
  on the same samples; single-source values mean only one was
  informative or available.
- **Used by.** Quality stratification — candidates seeded from both
  sources carry more confidence than single-source ones.

#### `q2_pca_seed_mode`
- **Comes from.** `block_data$flashlight_mode`. Historical name
  preserved inside the block; the flat key uses the scientific name.
  Value mirrors `seed_source` with finer-grained provenance (e.g.
  distinguishes `flashlight` from `flashlight_legacy` by whether
  the v2 SV-breakpoint-proximity rule was applied).
- **Measures.** Same axis as `seed_source` but tracks the seeding
  algorithm's mode flag. Redundant with seed_source in most cases.
- **Used by.** Reproducibility tracking; no scoring consumer.

#### `q2_n_seed_class_conflicts`
- **Comes from.** `block_data$n_seed_conflict`. Count of samples
  where flashlight v2 and step03 disagreed on the class and the
  drop-on-conflict rule removed them from the seed set.
- **Measures.** Per-candidate cohort-level disagreement between the
  two independent class-assignment systems (SV-based and PCA-based).
  A high value flags intervals where the two evidence streams
  contradict each other — often a signal of complex internal
  dynamics (nested inversions, recombinants, composite structure).
- **Used by.** Diagnostic. Not currently a scoring gate; chat-9
  Finding AY notes that the drop-on-conflict implementation is
  stricter than the original chat-9 design's "priority-pick" spec.

#### `q2_n_seed_HOM_REF` / `q2_n_seed_HET` / `q2_n_seed_HOM_INV`
- **Comes from.** `block_data$n_seed_HOM_REF` etc. Class-stratified
  counts of surviving seeds.
- **Measures.** How balanced the seed set is across the three
  genotype classes. `seeded_kmeans` requires a minimum per-class
  count or it falls back to unsupervised.
- **Used by.** Diagnostic; ablation checks on seeding robustness.

#### `q2_pca_cluster_silhouette`
- **Comes from.** `block_data$silhouette_score`. Mean silhouette
  coefficient across samples, in [-1, +1], rounded to 4 dp.
- **Measures.** How well-separated the preliminary HOM_REF / HET /
  HOM_INV clusters are on the per-sample mean PC1 axis. Values near
  +1 indicate clean three-cluster structure (the inversion cleanly
  partitions samples); values near 0 indicate the clusters overlap
  substantially. Values near -1 indicate misclassification, but are
  rare given the seeded k-means constraint.
- **Used by.** Derives `q2_pca_cluster_separation_flag` (the clean /
  noisy binary). Direct scoring by C01d is intentionally NOT
  enabled for the first wiring pass (chat-13 Finding AV); the raw
  scalar is archived for later promotion decisions.

#### `q2_pca_cluster_separation_flag`
- **Comes from.** `block_data$decomp_quality`. Categorical:
  `clean` (silhouette ≥ 0.40), `noisy` (silhouette < 0.40), or
  `null` when silhouette itself is NA.
- **Measures.** Binarised version of the silhouette score for
  quick gating. The 0.40 cutoff is a v5-inherited default; per user
  direction it stands unless observed distributions break something.
- **Used by.** Annotation-only at the first wiring pass (chat-13 AV).
  Promotion to C01d / C01f scoring is a chat-16+ decision.

#### `q2_bic_gap_k3_vs_k2`
- **Comes from.** `block_data$bic_gap_k3_vs_k2`. Signed scalar:
  `BIC(k=2) - BIC(k=3)` under the PC1 mixture-model fit.
- **Measures.** Whether three clusters (HOM_REF / HET / HOM_INV) fit
  better than two (HOM_REF / HOM_INV collapsed). Positive values
  indicate k=3 wins — the three-cluster structure is real. Negative
  or near-zero values indicate the inversion may be effectively
  biallelic from the PCA perspective at this coverage (HETs either
  too few or too similar to one homozygous class to be resolvable).
- **Used by.** Diagnostic on the appropriateness of k=3 seeding at
  this candidate. Complements the silhouette — the two together
  discriminate "well-separated three clusters" from "well-separated
  two clusters, HET class is a label artefact".

#### `q2_het_phase_support_fraction`
- **Comes from.** `block_data$phase_concordance`. Number in [0, 1].
  Fraction of samples labelled HET by the PCA whose Clair3 WhatsHap
  phase blocks span the inversion interval with alternating
  haplotype identity across the breakpoints.
- **Measures.** Phase-level confirmation of heterozygous state.
  A true HET carries one inverted and one non-inverted haplotype
  and should show a phase-block switch at or near the breakpoints;
  a false-HET (e.g. a sample wrongly assigned HET by PC1 noise)
  will not. High values (>0.7) are strong concurrent evidence for
  true biallelic inversion polymorphism.
- **Used by.** Annotation; informs manuscript discussion of HET
  validity per candidate.

#### `q2_n_het_phase_supported`
- **Comes from.** `block_data$n_phase_supported`. Integer.
- **Measures.** Numerator of `q2_het_phase_support_fraction`.
  Count of HET samples with phase-block support.
- **Used by.** Sanity-check the fraction against the denominator
  (`q6_n_HET_prelim`).

#### `q2_n_samples_sv_pca_discordant`
- **Comes from.** `block_data$n_discordant`. Integer. Count of
  samples where the SV breakpoint evidence (flashlight anchor)
  places them in a class different from the PCA-cluster call.
- **Measures.** Per-candidate conflict between two orthogonal
  signals: SV-level (read alignment, split-read, discordant-pair
  signatures at the inversion breakpoints) and SNP-level (PC1 from
  LD / allele-frequency patterns). Low values are the norm; high
  values flag candidates where either the SV calls are noisy (e.g.
  BND/Manta false positives) or the PCA structure at the interval
  reflects non-inversion sources (family LD, residual admixture).
- **Used by.** Downstream audit table (see chat-11 carrier
  reconciliation). Not a scoring input for C01d at this pass.

#### `q2_n_hetDEL_breakpoint_reclass`
- **Comes from.** `block_data$n_cheat2_reclassified`. Integer.
- **Measures.** Number of samples the cheat2 filter moved from
  HOM_INV to HET because they carry a het-DEL signature at the
  inversion breakpoint region. Biologically, a het-DEL at the
  breakpoint suggests the sample is actually heterozygous for the
  inversion (one inverted + one non-inverted allele), and the
  DEL artefact reflects the reference-disagreement of the inverted
  allele's breakpoint in a short-read alignment.
- **Used by.** Provenance of per-sample class assignments;
  flagged in C01i_d_seal logs for reproducibility of the
  reclassification decisions.

#### `q2_pca_k_used`
- **Comes from.** `block_data$k_used`. Enum {2, 3}.
- **Measures.** Which k the seeded k-means ran at. k=3 by default
  (for candidates with enough HET seeds); k=2 as a fallback when
  the HET class has too few seeds to support a stable third cluster.
- **Used by.** Diagnostic on the class-resolution level for this
  candidate. q6_n_HET_prelim is NA (or implicitly zero) when k=2.

#### `q6_n_HOM_REF_prelim` / `q6_n_HET_prelim` / `q6_n_HOM_INV_prelim`
- **Comes from.** `block_data$n_HOM_REF_prelim` etc. Count of
  samples assigned to each class by the preliminary seeded k-means,
  before C01i_d_seal applies the RECOMBINANT labels.
- **Measures.** Raw per-class sample counts at the PCA preliminary
  stage. The `_prelim` suffix distinguishes them from the FINAL
  counts (`q6_n_HOM_REF` etc.) which come from the seal's
  synthesised `sample_registry` groups table and include the
  recombinant-subgroup carve-outs.
- **Used by.** Frequency calculation (`q6_freq_inv`) uses final
  counts; the `_prelim` versions are tracked separately so the
  decompose-to-seal pipeline stage is auditable per-candidate.

---

## Block 2: `recombinant_map`

**Writer.** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R`

**Inputs.**
- The DAG-derived per-sample regime-deviation signal R, produced by
  `lib_recomb_combination.R::derive_R_from_regime` against the C01j
  `regime_memberships.tsv.gz` sidecar. Each sample's per-window
  regime labels across the candidate interval are collapsed via
  `rle()`; each run is a graph node with a window-count weight, and
  transitions are edges. A sample fires R if its
  `deviation_fraction ≥ 0.15` AND `longest_deviation_bp ≥ 50 kb`.
- The per-sample GHSL call G from
  `lib_ghsl_confirmation.R::ghsl_per_sample_in_interval` —
  categorical per-sample {SPLIT, UNIFORM, UNINFORMATIVE, ...}
  derived from within-sample haplotype divergence at the candidate
  interval.
- The per-SNP arrangement-discordant tract detector (per-sample
  gene-conversion-tract count + tract widths) from
  `gene_conversion_detector.R`.
- (Optional) cheat24 recombinant prior
  (`cheat24_recombinant_prior.R::classify_recombinant_event`) for
  per-sample posterior probability; if unreachable, posterior is
  NA cohort-wide (chat-12 Finding V — no inline fallback, to avoid
  two candidates processed on the same HPC run getting different
  posteriors depending purely on cheat24 file reachability).

**What it does.**
1. For every sample, combines R and G via the gate rule:
   - R ∧ G → `RECOMBINANT`, confidence HIGH
   - R ∧ ¬G, GHSL resolution insufficient → `RECOMBINANT`, MEDIUM
   - R ∧ ¬G, GHSL resolution sufficient → `recomb_disputed`
   - ¬R ∧ G → `recomb_ghsl_only`
2. Within RECOMBINANT, classifies event_class by
   `longest_deviation_bp` relative to `gate_params.min_dco_bp`
   (default 200 kb) and by presence of gene-conversion tracts:
   - `longest_deviation_bp ≥ 200 kb` → `double_crossover`
   - else if GC tracts present → `gene_conversion_embedded`
   - else → `single_crossover_or_ambiguous`
3. Records per-sample s1 (per-window class-switch segment from the
   decompose per_window_class RDS), s2 (WhatsHap phase switches),
   s3 (flashlight hemizygous segments) as three independent
   confirming signals.

**What it is NOT.** This block does not write group membership
(that is seal's job via `sample_registry`). It characterises
per-sample recombinant STATE, including the disputed / ghsl_only
cohorts that do NOT become RECOMBINANT in the final group registration.

### Keys

#### `q2_n_recombinant_samples`
- **Comes from.** `block_data$n_recombinants`. Integer. Count of
  samples with `recomb_status == "RECOMBINANT"` (HIGH + MEDIUM
  confidence). Does NOT count `recomb_disputed` or `recomb_ghsl_only`.
- **Measures.** Cohort-level prevalence of crossover events inside
  the candidate interval. A non-zero value is the signature that
  recombination is occurring across the inversion breakpoints in
  at least some samples, which itself constrains the inversion
  mechanism interpretation (recent-origin inversions often show
  more recombinants than ancient ones).
- **Used by.** C01d Layer C scoring, seal's group-registration
  decision tree, Q2 reporting in the final candidate status.

#### `q2_n_recombinant_gc_samples`
- **Comes from.** `block_data$n_recombinant_gc`. Count of RECOMBINANTs
  with `event_class == "gene_conversion_embedded"` AND
  `recomb_subgroup == "RECOMBINANT_GC"`.
- **Measures.** How many samples have short arrangement-discordant
  IBS tracts consistent with gene conversion events between inverted
  and non-inverted haplotypes. Gene conversion at inversion
  boundaries is a known mechanism of boundary erosion.
- **Used by.** Seal registers an `inv_<cid>_RECOMBINANT_GC` group
  for these samples. Q2 reporting; manuscript mechanism discussion.

#### `q2_n_recombinant_dco_samples`
- **Comes from.** `block_data$n_recombinant_dco`. Count of RECOMBINANTs
  with `event_class == "double_crossover"` AND
  `recomb_subgroup == "RECOMBINANT_DCO"`.
- **Measures.** How many samples show evidence of double-crossover
  events inside the inversion — longer (>min_dco_bp) mosaic tracts
  consistent with two independent crossovers flanking an internal
  retained segment.
- **Used by.** Seal registers an `inv_<cid>_RECOMBINANT_DCO` group.

#### `q2_n_regime_ghsl_disputed_samples`
- **Comes from.** `block_data$n_disputed`. Count of samples where R
  fired (regime DAG flagged a deviation) but GHSL reported a
  non-SPLIT call with sufficient resolution to argue against
  recombination.
- **Measures.** Per-candidate disagreement between the two
  orthogonal recombinant-detection signals. A high value flags
  candidates where either the regime-DAG signal is sensitive to
  family LD structure (false positives on R) or the GHSL divergence
  scheme is insensitive to the recombinant flavour at this interval.
- **Used by.** Diagnostic; informs the threshold-review discussion
  for AW (min_dco_bp) and the R-gate parameters.

#### `q2_n_ghsl_split_only_samples`
- **Comes from.** `block_data$n_ghsl_only`. Count of samples where
  GHSL called SPLIT but the regime DAG did NOT fire R.
- **Measures.** The complement of the disputed set — cases where
  haplotype-level divergence suggests within-sample recombination
  but the cohort-level regime membership does not change across the
  interval. Possible biological reads: small-scale events below
  the regime-window resolution, or short GC tracts that show at
  the SNP level without altering the regime label.
- **Used by.** Diagnostic; not registered as a formal group.

#### `q2_recomb_posterior_source`
- **Comes from.** `block_data$cheat24_version`. String.
  `cheat24_recombinant_prior.R` when the external cheat24 was
  reachable and the posterior was computed; `cheat24_unavailable_posterior_na`
  when it was not.
- **Measures.** Provenance marker for the per-sample posterior
  probability column in the `recombinants` array. When
  `_unavailable_posterior_na`, every posterior is NA and no
  cohort-wide aggregate of posteriors is meaningful.
- **Used by.** Audit; flags runs where posterior-dependent keys
  (none in the current canonical schema) would need to be
  quarantined.

#### `q2_recomb_gate_rule_version`
- **Comes from.** `block_data$combination_rule`. Constant string
  `lib_recomb_combination.R (DAG R + G gate, chat-12)` for this
  schema version.
- **Measures.** Provenance; records which combination-rule
  implementation produced the recombinant calls.
- **Used by.** Reproducibility tracking across chats/versions.

#### `q2_recomb_min_regime_dev_fraction`
- **Comes from.** `block_data$gate_params.min_deviation_fraction`.
  Default 0.15.
- **Measures.** The R-fire threshold on the per-sample
  `deviation_fraction` (= 1 − dominant_run_windows / total_windows).
  A sample fires R if AT LEAST this fraction of its interval is
  in a non-dominant regime.
- **Used by.** Reproducibility; enables per-candidate sensitivity
  analysis if the threshold is tuned.

#### `q2_recomb_min_regime_dev_bp`
- **Comes from.** `block_data$gate_params.min_deviation_bp`.
  Default 50 kb.
- **Measures.** The R-fire minimum bp length of the longest
  contiguous non-dominant regime run. Together with the fraction
  threshold, rules out fraction-only calls on very short deviations
  that are likely to be noise.
- **Used by.** Reproducibility.

#### `q2_recomb_dco_threshold_bp`
- **Comes from.** `block_data$gate_params.min_dco_bp`.
  Default 200 kb.
- **Measures.** The threshold above which a RECOMBINANT's
  `longest_deviation_bp` is labelled `double_crossover` rather than
  `gene_conversion_embedded` or `single_crossover_or_ambiguous`.
  A "double crossover" interpretation requires the mosaic tract
  to be long enough that two independent crossover events are a
  more parsimonious explanation than one crossover plus a gene
  conversion tract.
- **Used by.** Reproducibility; the 200 kb default is inherited
  from v5 as an educated guess (chat-12 notes).

---

## Block 3: `internal_ancestry_composition`

**Writer.** `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py`

**Inputs.**
- Engine B per-window per-sample ancestry labels from the local-Q
  cache (`<q_cache_dir>/<chrom>.local_Q_samples.tsv.gz`). Each
  label is an integer population index produced by the Engine B
  K-mode ancestry deconvolution.
- Candidate intervals (chrom, start, end) from the C01d candidate
  table.

**What it does.** For each parent candidate:
1. Slices the per-window ancestry labels inside the interval for
   every sample.
2. Classifies each sample's within-interval pattern into one of:
   - `homogeneous` — one label dominates
   - `dominant_plus_secondary` — one label dominates but a
     smaller secondary is persistent
   - `two_block_composite` — the interval is cleanly split into
     two ancestry blocks (suggests a real internal structural
     break within the candidate)
   - `continuous_gradient` — ancestry labels shift gradually
     along the interval (likely recombination-driven, not a
     structural inversion boundary)
   - `multi_block_fragmented` — three or more alternating blocks
   - `diffuse_mixed` — no clear pattern
3. Rolls cohort-level percentages up and derives a
   `composite_flag`:
   - `clean` — homogeneous-like patterns dominate (≥80%)
   - `maybe_composite` — some two_block / multi_block / gradient
     evidence (5-20%)
   - `likely_composite` — >20% two_block or multi_block
   - `unknown_no_engine_b` — Engine B cache not available; stub
   - `unknown_no_data` — no per-window data in interval

**Why this matters.** "Composite intervals" are candidates where
the inversion signal may be an artefact of ancestry substructure
rather than a real inversion. Two inversions within different
ancestry backgrounds can mimic a single inversion at a PCA window;
the composite flag catches these and caps them at UNCERTAIN so
C01f cannot promote them to SUPPORTED / VALIDATED.

### Keys

#### `q1_composite_flag`
- **Comes from.** `block_data$composite_flag`. Categorical:
  `clean`, `maybe_composite`, `likely_composite`,
  `unknown_no_engine_b`, `unknown_no_data`.
- **Measures.** Per-candidate verdict on whether the "inversion"
  signal is confounded by ancestry substructure. `likely_composite`
  is the hard stop.
- **Used by.** `STEP_C01f` promotion gate — `likely_composite`
  caps the candidate at UNCERTAIN regardless of other evidence.
  `maybe_composite` allows promotion but logs it. This is the
  single most consequential flat key for the group-validation
  decision tree. Schema key name is `q1_composite_flag` (v10.1
  canonical rename; chat-9 Finding X removed the duplicate
  `q1_ancestry_composite_flag` alias).

#### `q1_ancestry_dominant_pattern`
- **Comes from.** `block_data$dominant_structure_type`. The most
  common per-sample structure label across the cohort.
- **Measures.** The modal within-sample ancestry pattern for this
  candidate. `homogeneous` modes are strong evidence for a single-
  population signal; `two_block_composite` modes flag candidates
  where most samples individually show two-block structure (likely
  a real internal structural break).
- **Used by.** Manuscript Q1 reporting; informs mechanism discussion.

#### `q2_pct_samples_two_block_ancestry`
- **Comes from.** `block_data$pct_two_block_composite`. Fraction
  in [0, 1].
- **Measures.** Cohort-level prevalence of the two-block ancestry
  pattern at this candidate. If >20% of samples show two-block
  composition, the cohort-level signal is likely a composite of two
  distinct structural variants at the same interval.
- **Used by.** Direct input to `composite_flag` derivation.

#### `q2_pct_samples_multi_block_ancestry`
- **Comes from.** `block_data$pct_multi_block_fragmented`.
- **Measures.** Cohort-level prevalence of the multi-block
  fragmented pattern (three or more alternating ancestry blocks
  within the interval). High values indicate either highly complex
  internal structure or strong within-sample recombination / noise.
- **Used by.** Direct input to `composite_flag`.

#### `q2_pct_samples_homogeneous_ancestry`
- **Comes from.** `block_data$pct_homogeneous_like`. Sum of
  `homogeneous` and `dominant_plus_secondary`.
- **Measures.** Prevalence of "clean" per-sample patterns —
  samples whose within-interval ancestry is dominated by one
  label throughout. High values (≥80%) are one of the gating
  conditions for `composite_flag == "clean"`.
- **Used by.** Direct input to `composite_flag`.

#### `q2_pct_samples_gradient_ancestry`
- **Comes from.** `block_data$pct_continuous_gradient`.
- **Measures.** Prevalence of gradient patterns — samples whose
  ancestry labels shift smoothly along the interval rather than
  at a sharp boundary. Gradients are the canonical recombination
  signature; high values argue against a single sharp inversion
  boundary and for recombinant cohorts dominating the interval.
- **Used by.** `composite_flag`; Q2 reporting.

#### `q2_pct_samples_diffuse_ancestry`
- **Comes from.** `block_data$pct_diffuse_mixed`.
- **Measures.** Prevalence of diffuse / mixed patterns — samples
  with no discernible structure. Elevated values often flag
  low-coverage samples, Engine B cache gaps, or interval mismatches.
- **Used by.** Diagnostic; informs quality review.

#### `q2_mean_ancestry_fragmentation`
- **Comes from.** `block_data$mean_fragmentation`. Cohort mean
  of per-sample `fragmentation_score`.
- **Measures.** Average number of distinct ancestry blocks per
  sample within the interval, normalised. Higher values indicate
  more fragmentation (more boundaries within each individual
  sample). Complements the categorical patterns with a continuous
  quantification.
- **Used by.** Q2 scoring (if promoted by a future chat); currently
  annotation-only.

#### `q2_mean_ancestry_entropy_within_sample`
- **Comes from.** `block_data$mean_internal_entropy`. Cohort mean
  of per-sample Shannon entropy over the ancestry-label distribution
  within the interval.
- **Measures.** Per-sample within-interval label diversity.
  Low entropy ≈ one-label dominance; high entropy ≈ even mix across
  labels. Averaged across samples.
- **Used by.** Q2 reporting.

#### `q2_mean_ancestry_switches_per_sample`
- **Comes from.** `block_data$mean_switches`. Cohort mean of the
  per-sample count of adjacent-window ancestry-label changes.
- **Measures.** Average number of ancestry-label switches per
  sample within the interval. More direct measure of fragmentation
  than the normalised fragmentation_score.
- **Used by.** Q2 reporting.

#### `q2_ancestry_K_used`
- **Comes from.** `block_data$K_used`. Engine B K. Default 8.
  Value 0 in stub blocks.
- **Measures.** Which K-mode Engine B deconvolution was used for
  the ancestry labels — determines the granularity of the ancestry
  call. Higher K means finer subdivision but more noise.
- **Used by.** Reproducibility; comparability across candidates.

#### `q2_ancestry_n_samples_analyzed`
- **Comes from.** `block_data$n_samples_analyzed`. Integer.
- **Measures.** Number of samples that had enough per-window
  Engine B coverage inside the interval to be classified. Can be
  less than 226 if some samples have missing Q cache or too few
  windows.
- **Used by.** Denominator sanity-check for the `pct_*` keys.

---

## Pre-rename / post-rename key map

Chat-15 BK applied a systematic rename from procedural / mechanical
labels to scientific names that say what each key measures. The full
map:

### internal_dynamics (12 renamed / 7 preserved)

| Old | New | Reason |
|-----|-----|--------|
| `q2_decomp_status` | `q2_pca_decomp_status` | qualifies which decomposition |
| `q2_decomp_reason` | `q2_pca_decomp_skip_reason` | fires only on skip |
| `q2_seed_source` | `q2_pca_seed_source` | clarifies what was seeded |
| `q2_n_seed_conflict` | `q2_n_seed_class_conflicts` | says "class" conflicts |
| `q2_silhouette_score` | `q2_pca_cluster_silhouette` | names what is measured |
| `q2_decomp_quality` | `q2_pca_cluster_separation_flag` | names what is flagged |
| `q2_phase_concordance` | `q2_het_phase_support_fraction` | fraction of what |
| `q2_n_phase_supported` | `q2_n_het_phase_supported` | specifies samples |
| `q2_flashlight_mode` | `q2_pca_seed_mode` | mode of what |
| `q2_n_discordant` | `q2_n_samples_sv_pca_discordant` | which discordance |
| `q2_n_cheat2_reclass` | `q2_n_hetDEL_breakpoint_reclass` | biological trigger |
| `q2_k_used` | `q2_pca_k_used` | qualifies which k |

Preserved: `q2_bic_gap_k3_vs_k2`, `q2_n_seed_HOM_REF`, `q2_n_seed_HET`,
`q2_n_seed_HOM_INV`, `q6_n_HOM_REF_prelim`, `q6_n_HET_prelim`,
`q6_n_HOM_INV_prelim` (the names already say what they measure).

### recombinant_map (10 renamed / 0 preserved)

| Old | New | Reason |
|-----|-----|--------|
| `q2_n_recombinants` | `q2_n_recombinant_samples` | count of what |
| `q2_n_recombinant_gc` | `q2_n_recombinant_gc_samples` | parallel form |
| `q2_n_recombinant_dco` | `q2_n_recombinant_dco_samples` | parallel form |
| `q2_n_recomb_disputed` | `q2_n_regime_ghsl_disputed_samples` | names the conflict |
| `q2_n_recomb_ghsl_only` | `q2_n_ghsl_split_only_samples` | what only GHSL saw |
| `q2_cheat24_version` | `q2_recomb_posterior_source` | what the version tracks |
| `q2_combination_rule` | `q2_recomb_gate_rule_version` | what rule and why |
| `q2_recomb_min_dev_fraction` | `q2_recomb_min_regime_dev_fraction` | fraction of what |
| `q2_recomb_min_dev_bp` | `q2_recomb_min_regime_dev_bp` | parallel |
| `q2_recomb_min_dco_bp` | `q2_recomb_dco_threshold_bp` | what the threshold is for |

### internal_ancestry_composition (9 renamed / 3 preserved)

| Old | New | Reason |
|-----|-----|--------|
| `q1_ancestry_dominant_structure` | `q1_ancestry_dominant_pattern` | "pattern" used in paper |
| `q2_pct_samples_two_block` | `q2_pct_samples_two_block_ancestry` | what kind of two-block |
| `q2_pct_samples_multi_block` | `q2_pct_samples_multi_block_ancestry` | parallel |
| `q2_pct_samples_homogeneous` | `q2_pct_samples_homogeneous_ancestry` | parallel |
| `q2_pct_samples_gradient` | `q2_pct_samples_gradient_ancestry` | parallel |
| `q2_pct_samples_diffuse` | `q2_pct_samples_diffuse_ancestry` | parallel |
| `q2_mean_fragmentation` | `q2_mean_ancestry_fragmentation` | fragmentation of what |
| `q2_mean_internal_entropy` | `q2_mean_ancestry_entropy_within_sample` | entropy of what, where |
| `q2_mean_switches` | `q2_mean_ancestry_switches_per_sample` | what is counted |

Preserved: `q1_composite_flag` (v10.1 canonical, chat-9 Finding X already
did this rename the other direction), `q2_ancestry_K_used`,
`q2_ancestry_n_samples_analyzed` (already carry the `ancestry_` prefix).
