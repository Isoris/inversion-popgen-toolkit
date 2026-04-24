# =============================================================================
# DESIGN RATIONALE — gene_conversion_detector and GC/recombinant separation
# =============================================================================
# Date:   2026-04-24 (consolidated from chat-10, chat-12, chat-17, chat-D
#         discussions)
# Scope:  Why the gene_conversion_detector looks the way it does. What
#         choices are literature-based, what choices are operational, what
#         choices are admitted unknowns. Written so sections below can be
#         lifted into the manuscript Methods with minimal rewriting.
# Files:  inversion_modules/phase_7_karyotype_groups/proposal/
#           gene_conversion_detector.R
#           RUN_gene_conversion_detector.R
#           cheats/cheat24_recombinant_prior.R
#         registries/schemas/structured_block_schemas/
#           gene_conversion_tracts.schema.json
# =============================================================================


# ─────────────────────────────────────────────────────────────────────────────
# 1. What the detector actually measures (honest naming)
# ─────────────────────────────────────────────────────────────────────────────

The detector is named `gene_conversion_detector.R` and writes a block
called `gene_conversion_tracts`. The name is **shorthand** kept for
pipeline back-compatibility. What the detector actually surfaces is:

  SHORT ARRANGEMENT-DISCORDANT TRACTS within a candidate inversion:
  stretches of ~2–10 consecutive diagnostic SNPs where a single sample's
  dosage matches the OPPOSITE arrangement's baseline, then reverts.

What we can say from the sequence data:
  "Sample X, baseline HOM_REF, looks INV-like across SNPs at positions
   p1..pk inside this interval."

What we CANNOT say without additional experiments (trios, long-read
phased haplotypes, mutation-accumulation panel):
  "Sample X underwent a gene-conversion event at this locus, with the
   INV arrangement serving as the donor template."

This matters for the MANUSCRIPT. Candidate causes of an
arrangement-discordant IBS tract include:
  - Gene conversion (classical interpretation; Harris 1993 GC tracts
    18–774 bp; Korunes & Noor 2017 Drosophila/mouse literature)
  - Short double crossover
  - Paralog mismapping at the reference
  - Coincidental IBS with the other arrangement at a handful of
    diagnostic sites

**Manuscript wording rule.** Use "arrangement-discordant IBS tracts
consistent with gene conversion" or "GC-like tracts" in figure legends
and Methods/Results. Reserve "gene conversion events" for cases
supported by independent evidence.

Ref for this decision: 2026-04-17 audit finding BF. Block/file names
kept for pipeline back-compat; terminology choice localized to
manuscript text.


# ─────────────────────────────────────────────────────────────────────────────
# 2. Why per-SNP run-length, not windowed binning
# ─────────────────────────────────────────────────────────────────────────────

The detector was REWRITTEN in chat-12 (2026-04-18). The prior
implementation was a windowed-binning detector:
  - 40-SNP windows, 10-SNP step
  - CUSUM-style excursion sweep on window-mean dosage
  - `max_tract_bins` gate on the excursion length

Why this was replaced. At our SNP density (~1.6 kb/SNP after MAF
filtering at 9× coverage), a 40-SNP window spans ~60 kb. That is
LARGER than real GC biology:
  - Harris 1993: GC tracts 18–774 bp in yeast S. cerevisiae; upper
    bound ~5 kb is atypical.
  - Korunes & Noor 2017: Drosophila mean tract length ~380 bp;
    distribution geometric with p ≈ 0.0025.
  - Hellenthal & Stephens 2006: human PRDM9-associated tracts
    ~200-300 bp.

A 60-kb detection floor therefore made the real signal INVISIBLE.
Worse, the signals it DID surface were often SNP-column artefacts
(vertical stripes in the heatmap) mislabelled as tracts.

The rewrite uses per-SNP run-length with tolerance:
  - `min_run = 2` SNPs (a minimum of 2 consecutive flagged SNPs)
  - `max_tolerance = 1` (allow 1 intervening non-flagged SNP inside
    a tract, to accommodate missing calls and occasional genotyping
    noise)
  - Confidence scored by run length under a geometric prior:
    P(k consecutive genotyping errors at diagnostic SNPs) ≈ p^k at
    per-SNP error rate p ≈ 0.01 at 9× coverage. run_length=2 gets
    LOW confidence, 3 MEDIUM, ≥4 HIGH.

This matches the heatmap-visual signal: short horizontal streaks in
the sample × SNP dosage plot, one sample wide, 2-3 SNPs long.


# ─────────────────────────────────────────────────────────────────────────────
# 3. Two separate threshold logics — lower vs upper bound
# ─────────────────────────────────────────────────────────────────────────────

The detector has length gates at both ends. They are set by DIFFERENT
logic and should not be conflated.

## 3a. Lower bound — literature-justified (real biology)

`min_run = 2` SNPs. `min_delta_af = 0.5` for the diagnostic mask.

Justification: we care to CATCH real GC biology. Harris 1993 and
follow-up literature establish typical GC tract lengths at 50 bp – ~800
bp with an atypical tail up to ~5 kb. At our SNP density (~1.6 kb/SNP),
a 500 bp tract often covers 0 or 1 diagnostic SNPs and is unrecoverable.
A 1-2 kb tract covers 1-2 SNPs. A 3-5 kb tract covers 2-3 SNPs.

`min_run = 2` therefore sits right at the detection floor our marker
density permits. Lower than 2 would be pure single-SNP genotyping
noise; higher than 2 would miss the biologically typical short end.

This is a CARE-ABOUT-BIOLOGY setting.

## 3b. Upper bound — operational routing (NOT biology)

`max_span_bp = 20000`. `max_flagged_snps = 10`.

Justification: this is a **handoff boundary** between the GC detector
and C01j's recombinant machinery. Tracts shorter than this are routed
to the GC detector; tracts longer than this are treated as recombinant
events (short double crossover or single crossover with breakpoint
rescue) and handled by C01j.

From the detector's own header: "Longer events are recombinants,
handled by C01j — discard here."

Nothing about the 20 kb number is a biological ceiling. We do not
know, from sequence data alone, whether a 15 kb arrangement-discordant
tract in our cohort is:
  - A genuine gene conversion (atypically long)
  - A short double crossover
  - Paralog mismapping over a longer region
  - Coincidental IBS at many diagnostic sites

The 20 kb line is drawn for one purpose: **do not double-count the
same event in both pipelines**. Without a cutoff, a 15 kb
arrangement-discordant tract would appear as both a GC tract (here)
and a recombinant event (C01j), inflating both counts.

## 3c. How 20 kb was picked — heatmap-visual separability

The 20 kb number came from the DOSAGE HEATMAP, not the literature.
The chat-12 rewrite was motivated by wanting to match what is visually
separable in a sample × SNP dosage heatmap at our SNP density:

  Visual signal in heatmap     →  What it is              →  Routed to
  ───────────────────────────────────────────────────────────────────
  Short horizontal stripe         Single sample, 1-3         GC detector
  (one row, 1-3 columns)          diagnostic-SNP run
  ───────────────────────────────────────────────────────────────────
  Vertical stripe                 Non-diagnostic SNP or     Filtered by
  (many rows, 1 column)           paralog mismapping         diagnostic_snp_mask
  ───────────────────────────────────────────────────────────────────
  Wide horizontal flip            One sample, many SNP      C01j
  (one row, many columns)         columns = recombinant
  ───────────────────────────────────────────────────────────────────
  Color block covering many       Real karyotype class      C01i_decompose
  rows (cross-sample signal)      signal (the inversion)

At ~1.6 kb/SNP, 20 kb ≈ 12 SNPs. Anything wider than that is clearly
a recombinant-shaped stripe flip in the heatmap, not a fleck. This is
how the number was set. It is a heatmap-visual threshold
operationalized in bp.


# ─────────────────────────────────────────────────────────────────────────────
# 4. Why per-sample run-length and NOT correlation-core-block
# ─────────────────────────────────────────────────────────────────────────────

The phase 6 breakpoint refinement module (`01_dosage_signal.R`) uses a
DIFFERENT detection algorithm: correlation-core-block + extend. Both
modules run on the same dosage matrix; both use a diagnostic-SNP
notion. One could ask: should the GC detector use correlation too, or
phase 6 use run-length?

No. The two algorithms are paired with their tasks by algorithmic
structure, not by accident.

Phase 6 is looking for a CROSS-SAMPLE signal: every HOM_INV sample
deviates from every HOM_REF sample at the same markers, across the
whole inversion. Correlation-core-block is the right tool because
marker-to-marker correlation across all samples is exactly what that
signal looks like. A per-SNP run-length scan would scan each sample
independently, then have to aggregate, which throws away the
cross-sample correlation structure that separates signal from noise
at Mb scale.

The GC detector is looking for a SINGLE-SAMPLE event: one carrier's
short tract of opposite-arrangement-looking genotypes. Per-SNP
run-length scans each sample independently, which is the right
shape. Correlation-core-block averages across samples, which is
precisely what erases single-sample deviations — one GC tract in one
sample out of 226 is a 1/226 perturbation at that marker, invisible
to cross-sample correlation.

**Structural summary:**

  Task                          Signal shape        Right algorithm
  ─────────────────────────────────────────────────────────────────
  Locate the inversion          Cross-sample,       correlation-
  (phase 6)                     Mb-scale, shared    core-block + extend
  ─────────────────────────────────────────────────────────────────
  Locate short discordant       Single-sample,      per-SNP
  tracts inside it              kb-scale, local     run-length
  (GC detector)

Swapping would put each algorithm on the task it is worst at. The
algorithms are correctly sited; only the *diagnostic-SNP mask
computation* is genuinely shareable between the two modules (same
math, currently duplicated — see CODE_REVIEW_TO_FIGURE_OUT.md §5).


# ─────────────────────────────────────────────────────────────────────────────
# 5. Dosage vs karyotype class vs arrangement — what the detector is reading
# ─────────────────────────────────────────────────────────────────────────────

This section disentangles three concepts that are easy to confuse when
reading the detector code.

## 5a. Dosage is per-SNP, per-sample. It is NOT the karyotype class.

At one SNP, for one sample:
  dosage = P(AB) + 2·P(BB)
where A is the major allele and B is the minor allele at that SNP (per
ANGSD convention in STEP_A01_beagle_to_dosage_by_chr.py). So:
  dosage ≈ 0  → confident AA (homozygous major)
  dosage ≈ 1  → confident AB (heterozygous)
  dosage ≈ 2  → confident BB (homozygous minor)
  intermediate values reflect genotype-likelihood uncertainty.

Dosage is a per-allele measurement at a per-SNP position. It says
nothing about arrangement on its own.

## 5b. Karyotype class (HOM_REF/HET/HOM_INV) is a derived region-level label.

HOM_REF / HET / HOM_INV are labels that describe which inversion
arrangement a sample carries across the WHOLE candidate interval.
They are derived by clustering samples on a region-level signal (PC1
loadings from the interval's dosage matrix, k-means on the loading
distribution, etc. — see STEP_C01i_decompose.R). One karyotype label
per sample per candidate, not per SNP.

## 5c. Arrangement predicts dosage AT DIAGNOSTIC SNPs only.

The relationship between karyotype class and per-SNP dosage is
probabilistic and SNP-specific:

  At a DIAGNOSTIC SNP (|AF_HOM_REF − AF_HOM_INV| ≥ 0.5):
    - A HOM_REF sample is expected to have dosage near 2·AF_HOM_REF[i]
      (which is close to 0 or close to 2, depending on which allele
      happens to be major on the REF background at SNP i).
    - A HOM_INV sample is expected to have dosage near 2·AF_HOM_INV[i]
      (close to the opposite end).
    - HET samples fall in the middle.

  At a NON-DIAGNOSTIC SNP (allele frequencies similar across
  arrangements):
    - A sample's dosage is not predicted by its arrangement. The
      sample carries whatever alleles it inherited regardless of
      arrangement.

So at diagnostic SNPs, the arrangement → dosage relationship is a
*prediction*, not a deterministic map. Most of the time the prediction
holds. When it fails at a localized stretch inside a sample, that is
the GC-like signal.

## 5d. Why a single sample has different dosages across a range of markers

Each SNP has its own pair of alleles and its own per-arrangement
frequencies. Even inside an inversion, consecutive diagnostic SNPs can
point in opposite directions:

  SNP 1   AF_REF = 0.05, AF_INV = 0.95  → REF predicts dosage ≈ 0.1,
                                          INV predicts dosage ≈ 1.9
  SNP 2   AF_REF = 0.90, AF_INV = 0.10  → REF predicts dosage ≈ 1.8,
                                          INV predicts dosage ≈ 0.2
  SNP 3   AF_REF = 0.50, AF_INV = 0.50  → NON-DIAGNOSTIC

The detector treats per-SNP group means as the prediction and flags
samples by distance to the opposite prediction. A HOM_REF sample is
"INV-consistent at SNP i" iff its dosage at SNP i is closer to
AF_HOM_INV[i]·2 than to AF_HOM_REF[i]·2. Diagnostic SNPs are exactly
those where the two predictions are far apart, so the comparison is
crisp.

**Implementation note — this was fixed 2026-04-24 (chat D).** The
chat-12 rewrite used fixed dosage thresholds (≤ 0.5 = consistent, ≥
1.5 = flag for HOM_REF; symmetric for HOM_INV). That worked when
ANGSD's major/minor call aligned with arrangement polarity — which
holds when the REF arrangement is much more common — but failed at
reverse-polarity diagnostic SNPs (where a HOM_REF sample correctly
carries dosage ≈ 2). The fix reads af_ref / af_inv off the diagnostic
mask (attached as attributes in diagnostic_snp_mask) and does the
polarity-aware comparison at each SNP. Fixed-threshold fallback path
retained for callers that don't pass af_ref / af_inv, with a warning.
Test_gc_detector.R tests 11, 12, 13 exercise the reverse-polarity and
mixed-polarity cases.

## 5e. Why a single sample can carry a GC signal

GC events happen in a single meiosis in a single individual's
ancestry. They do not propagate across the cohort. So the signal is
by construction SINGLE-SAMPLE: one carrier shows a localized stretch
where its dosages match the opposite arrangement's predictions, then
reverts to its own arrangement's predictions.

Concretely, for a HOM_REF sample X at a HOM_REF-predicted-value-near-0
stretch of diagnostic SNPs:

  SNP 39: dosage 0.05  ✓ HOM_REF-like (matches AF_HOM_REF × 2)
  SNP 40: dosage 0.10  ✓ HOM_REF-like
  SNP 41: dosage 1.85  ✗ HOM_INV-like (flag)   ┐
  SNP 42: dosage 1.90  ✗ HOM_INV-like (flag)   │ the tract
  SNP 43: dosage 1.75  ✗ HOM_INV-like (flag)   ┘
  SNP 44: dosage 0.12  ✓ HOM_REF-like (revert)
  SNP 45: dosage 0.08  ✓ HOM_REF-like

Three consecutive flagged SNPs, spanning ~5 kb at our marker density,
all pointing in the same direction (sample X looks INV-typical when
it should look REF-typical). That is a GC-like tract.

Biologically: sample X carries the HOM_REF arrangement on both
chromatids. At a tiny stretch inside the inversion, one chromatid
acquired a short IBS match to the HOM_INV arrangement's typical
alleles. Candidate mechanisms are gene conversion, short double
crossover, paralog mismapping, or coincidental IBS — see §1. The
detector surfaces the pattern; cause is not recoverable from
sequence alone.

## 5f. Why HET samples are excluded from the scan

A HET sample has one REF chromatid and one INV chromatid. At a
diagnostic SNP, its dosage is expected to be intermediate
(~AF_HOM_REF + AF_HOM_INV). A GC tract on one of its two chromatids
would shift the dosage by only half a haplotype — from ~1.0 toward
~0.5 or ~1.5. At 9× coverage with unphased dosages, intermediate
dosages (0.5–1.5) are within the noise envelope of a normal HET
call. We cannot distinguish "HET sample with GC tract on one
chromatid" from "HET sample with noisy genotype calls" without
phased haplotypes.

This is a coverage limitation, not a biological one. At higher
coverage with phased data the scan extends to HET samples naturally
— the detector's per-sample structure already supports it.

## 5g. Why cross-sample detectors cannot see this

Per §4 above. GC is a 1/N perturbation at a given marker where N is
the cohort size (226 samples). Marker-to-marker correlation across
all samples is dominated by the (N-1) samples whose dosages match
their baseline prediction. The single carrier's deviation is well
below the correlation noise floor. Per-sample run-length scans can
see it because each sample is scanned independently against its own
baseline prediction. This is why the detection algorithm is
structurally per-sample, not why it happens to be implemented that
way.


# ─────────────────────────────────────────────────────────────────────────────
# 6. Why diagnostic SNPs only
# ─────────────────────────────────────────────────────────────────────────────

Only diagnostic SNPs — those with `|AF(HOM_REF) − AF(HOM_INV)| ≥ 0.5`
— carry information about arrangement. Non-diagnostic SNPs (shared
polymorphisms between arrangements) are skipped.

At a non-diagnostic SNP, both arrangements have similar allele
frequencies (e.g., REF AF = 0.45, INV AF = 0.48). A sample's dosage
at that SNP does not predict its arrangement, and no arrangement
predicts the sample's dosage. Including non-diagnostic SNPs in the
scan would flag random dosage variation as "arrangement-discordant,"
which is meaningless — the arrangement doesn't disagree with the
sample because the arrangement doesn't make a prediction here.

The vertical stripes often visible in raw dosage heatmaps are
non-diagnostic SNPs (and, separately, paralog-mismapping artefacts).
They are NOT GC events. The diagnostic_snp_mask plus excess-het
filtering removes both before the run-length scan.

**Implementation note.** The diagnostic mask is recomputed per
candidate, using the actual HOM_REF and HOM_INV assignments inside
that candidate's interval. It is NOT a precomputed genome-wide
"diagnostic SNP list" — diagnostic-ness is candidate-specific.


# ─────────────────────────────────────────────────────────────────────────────
# 7. Why HET samples are skipped from the scan
# ─────────────────────────────────────────────────────────────────────────────

Expanded above in §5f. Short version: HET samples' per-haplotype
signal is below the noise envelope at 9× with unphased dosages. A GC
tract on one of two chromatids shifts the dosage by half a
haplotype — from ~1.0 toward ~0.5 or ~1.5, which overlaps the normal
HET noise distribution. Without phased data we cannot distinguish
the two. This is a coverage-limited decision; the per-sample
detector structure readily extends to HET samples if/when phased
higher-coverage data becomes available.


# ─────────────────────────────────────────────────────────────────────────────
# 8. Confidence levels, not hard cutoffs
# ─────────────────────────────────────────────────────────────────────────────

Rather than a single hard threshold for "real tract vs noise," the
detector emits `confidence` ∈ {LOW, MEDIUM, HIGH} per tract based on
run length:

  run_length = 2     LOW    P(spurious) ~ p^2 = 1e-4 per position
  run_length = 3     MEDIUM P(spurious) ~ p^3 = 1e-6 per position
  run_length ≥ 4     HIGH   P(spurious) ~ p^4 = 1e-8 per position

With per-SNP genotyping error p ≈ 0.01 at 9× coverage. Note these are
PER-POSITION error probabilities; expected false-positive count across
the genome is `p^k × (n_diagnostic_SNPs × n_samples)`, which at
thousands of diagnostic SNPs × 226 samples puts a dozen or so
`run_length=2` false positives per candidate if no other filtering
applies.

Downstream consumers (C01i_b, C01f) can choose their own confidence
threshold depending on context. Figures that want to report "high
confidence GC tracts" should filter to `confidence = HIGH`. Sensitivity
analyses that want broader coverage can include MEDIUM or LOW.

This replaces the earlier binary keep/discard logic. It's more honest
about the per-sample-per-tract noise floor and lets the manuscript
tables report tract counts at stated confidence levels.


# ─────────────────────────────────────────────────────────────────────────────
# 9. Relationship to C01j (recombinant) and cheat24 (classification)
# ─────────────────────────────────────────────────────────────────────────────

Three modules touch the recombinant / GC space. They do different jobs.

  gene_conversion_detector   DETECTS short tracts (2–10 SNPs, ≤20 kb).
                              Surfaces arrangement-discordant runs.
                              Does NOT classify event cause.

  C01j compatibility engine   DETECTS long recombinants (Mb-scale
                              regime shifts). Surfaces recombinant
                              samples via per-window compatibility
                              clustering.

  cheat24_recombinant_prior   CLASSIFIES observed events by tract
                              length + position relative to
                              breakpoints:
                                short + near breakpoint → gene
                                  conversion prior (~1e-3, constant)
                                long + center             → double
                                  crossover prior (~1e-2, logistic
                                  increase to center)
                                long + near breakpoint    → suspicious,
                                  flag
                                short + center            → ambiguous

The detector hands its output to cheat24 via the recombinant_map
pipeline; cheat24 assigns the classification. The detector's length
gate (max_span_bp = 20000) and cheat24's classification thresholds
(MOSAIC_SHORT_BP = 50000, MOSAIC_LONG_BP = 200000) are NOT the same
numbers because they're answering different questions:

  Detector threshold   "Is this short enough for me to scan, or do I
                       hand it to C01j?"
  cheat24 thresholds   "Given the event exists at this length and
                       position, how should I classify its mechanism?"

The detector doesn't need to match cheat24 because the detector's job
ends at "here's a tract, here's how long it is" — cheat24 takes over
from there with its own thresholds applied to events the detector
surfaced.


# ─────────────────────────────────────────────────────────────────────────────
# 10. Draft Methods paragraph (lift and adapt)
# ─────────────────────────────────────────────────────────────────────────────

Below is a draft paragraph that can be lifted into Methods. Citations
need expansion; adjust parameter values if these change in the final
run.

-----------------------------------------------------------------------
Arrangement-discordant IBS tracts (hereafter "GC-like tracts") were
detected within each candidate inversion using a per-SNP run-length
scan on expected-dosage data. For each candidate, we first identified
diagnostic SNPs as those with |AF(HOM_REF) − AF(HOM_INV)| ≥ 0.5 (AF
computed within the samples assigned to each homozygous class by
C01i_decompose; diagnostic-ness is candidate-specific). At each
diagnostic SNP a sample was flagged when its expected dosage
(P(AB) + 2·P(BB), from ANGSD genotype likelihoods) was closer to the
OPPOSITE arrangement's per-SNP group mean than to its OWN. This
polarity-aware comparison is necessary because at a balanced cohort
(e.g., 60/106/60 at LG28) a non-trivial fraction of diagnostic SNPs
have reverse polarity (the REF-arrangement-typical allele is called
MINOR by ANGSD). HET samples were excluded because at 9× coverage
per-haplotype ambiguity precludes tract resolution. Consecutive runs
of ≥ 2 flagged SNPs with ≤ 1 intervening non-flagged SNP and total
span ≤ 20 kb were called as GC-like tracts. Longer discordant regions
were routed to the recombinant classification pipeline (C01j). Tract
confidence was scored by run length under a geometric per-SNP error
prior (p ≈ 0.01): run_length = 2 LOW, 3 MEDIUM, ≥ 4 HIGH. Tract
lengths at our marker density (~1.6 kb/SNP) resolve the typical
biological range reported for gene conversion in the literature
(Harris 1993: 18–774 bp; Korunes & Noor 2017: Drosophila mean
~380 bp). We report tracts as "arrangement-discordant IBS tracts
consistent with gene conversion" rather than "gene conversion
events" because the underlying cause cannot be distinguished from
short double crossovers, paralog mismapping, or coincidental IBS
from sequence data alone.
-----------------------------------------------------------------------


# ─────────────────────────────────────────────────────────────────────────────
# 11. Open questions — see CODE_REVIEW_TO_FIGURE_OUT.md
# ─────────────────────────────────────────────────────────────────────────────

The design choices above are the ones we're committed to. The doc
CODE_REVIEW_TO_FIGURE_OUT.md tracks open concerns, including:
  §5  Duplicated diagnostic-SNP mask computation between phase 6
      and the GC detector
  §6  Why breakpoint refinement runs before final decompose
  (others cross-reference the GC detector tangentially)

This file is the positive counterpart — decisions made, why we made
them, and what stands behind them. If anything here stops being true
in a future pass, update this file; do not let Methods text lose the
reasoning.
