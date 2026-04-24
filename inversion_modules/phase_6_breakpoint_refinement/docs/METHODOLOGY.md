# Methodology — Dosage-Based Inversion Breakpoint Refinement

**Project**: MS_Inversions_North_african_catfish
**Target resolution**: SNP-level with per-carrier uncertainty quantification
**Replaces**: 50 kb Fst-binned breakpoint estimation
**Reference implementation**: STEP_BP_01 through STEP_BP_04 in the followup pipeline

---

## 1. Motivation

Structural-variant callers (DELLY, Manta) and Fst-based shelf estimation each give coarse or unreliable breakpoint estimates for polymorphic inversions in mid-coverage resequencing cohorts. In our 226-sample pure *Clarias gariepinus* hatchery cohort at ~9× coverage:

- DELLY called an INV at 14.87–14.94 Mb on LG28, 200 kb proximal to the true breakpoint that we independently localized to 15.115 Mb from Fst shelf edge sharpening.
- Manta produced no call for the same inversion.
- Fst at 50 kb bins localizes breakpoints only to the resolution of one bin — we could say "between 15.10 and 15.15 Mb," no better.

We therefore developed a dosage-based breakpoint refinement procedure that uses the per-sample, per-SNP inheritance signal directly — requiring neither split-read evidence nor Fst binning — and that carries native per-carrier uncertainty.

---

## 2. Conceptual basis

An inversion arrangement, once it fixed in one lineage, propagates as an intact haplotype through subsequent generations except where recombination has cut across it. In a segregating inversion at frequency ~0.5, recombination within the inversion is suppressed in heterozygotes but can occur in homozygotes and, rarely, via double-crossover or gene conversion in heterozygotes.

Each inversion-carrying sample alive today therefore harbors a **segment of continuous inversion-ancestry** that may be shorter than the full inversion if any of that sample's ancestors underwent recombination that disrupted the inversion block. The individual extent of this segment is what we call the sample's **ancestral fragment**.

Three consequences follow:

**C1.** The union of all per-sample ancestral fragment boundaries on one side of the inversion clusters around the true breakpoint on that side. Samples whose ancestors never recombined vote for the true breakpoint position. Samples whose ancestors recombined inside the inversion vote for a position interior to the true breakpoint.

**C2.** The distribution of per-sample boundaries on each side is therefore right-skewed on the left side and left-skewed on the right side (toward the interior in both cases). The **modal** position of this distribution estimates the true breakpoint. The **spread** estimates uncertainty. The **tail** quantifies recombinant history.

**C3.** Ancestral fragment boundaries computed from dosage data have SNP-level resolution, bounded only by SNP density. At our ~8.6×10⁴ SNPs across LG28 (~20 Mb), mean SNP spacing is ~230 bp. Fragment-boundary precision of a few kb is achievable when sufficient carriers contribute.

This consequence-chain generalizes a seed observation already present in our inversion decomposition pipeline: **STEP_C01i::compute_ancestral_fragments()** performs the per-sample scan. What was missing was (a) aggregating the resulting fragment distributions into a population-level breakpoint estimate with CI, (b) fusing this estimate with independent dosage-based signals, and (c) replacing all Fst-based breakpoints in downstream outputs with these refined estimates.

---

## 3. Algorithm

### 3.1 Input definitions

Let $X \in \mathbb{R}^{M \times N}$ be the dosage matrix for a candidate region: $M$ SNPs in rows (sorted by genomic position $p_1 < p_2 < \ldots < p_M$), $N$ samples in columns. Entries $X_{ij} \in [0, 2]$ are BEAGLE posterior expected allele counts.

Let $G \in \{\text{HOMO\_1}, \text{HET}, \text{HOMO\_2}\}^N$ be the group assignment from STEP21's rotated PCA (`candidate_pca_rotated.tsv`, column `coarse_group_refined`).

Denote by $S^{\text{H1}}, S^{\text{HET}}, S^{\text{H2}}$ the index sets of samples in each group.

### 3.2 Informative markers (STEP_BP_01)

For each SNP $i$, compute the per-group mean dosage:

$$
\bar X^{\text{H1}}_i = \frac{1}{|S^{\text{H1}}|} \sum_{j \in S^{\text{H1}}} X_{ij}
\qquad
\bar X^{\text{H2}}_i = \frac{1}{|S^{\text{H2}}|} \sum_{j \in S^{\text{H2}}} X_{ij}
$$

(means taken over non-missing entries). The signed informativeness is:

$$
\delta_i = \bar X^{\text{H2}}_i - \bar X^{\text{H1}}_i
$$

A SNP is **informative** if $|\delta_i| \geq \delta_{\min}$. Default $\delta_{\min} = 1.0$ (i.e., the two homokaryotypes disagree by at least one expected dose). This threshold is strict — a fully segregating diallelic SNP tagging the inversion perfectly would produce $|\delta_i| \approx 2$.

### 3.3 Core block (STEP_BP_01, adapted from C01i)

Let $\mathcal{I}$ be the set of informative SNP indices. Compute the marker–marker correlation matrix $R$ restricted to $\mathcal{I}$:

$$
R_{ii'} = \text{cor}(X_{i,\cdot}, X_{i',\cdot}) \qquad i, i' \in \mathcal{I}
$$

Build an adjacency graph with edges where $|R_{ii'}| \geq \rho_{\text{block}}$ (default $\rho_{\text{block}} = 0.7$), then extract connected components. The **core block** is the largest connected component with at least $M_{\min} = 20$ markers. Denote its extent $[i_L, i_R]$ in SNP indices.

### 3.4 Block boundary extension (STEP_BP_01, adapted from C01i)

Compute the core consensus vector — the mean dosage of each sample over core markers:

$$
c_j = \frac{1}{|\mathcal{C}|} \sum_{i \in \mathcal{C}} X_{ij}
$$

where $\mathcal{C} \subseteq [i_L, i_R]$ is the top-tier core (markers in the 95th percentile of in-block correlation with the block consensus).

For each candidate extension marker $i \notin [i_L, i_R]$, compute:

$$
r_i = \left| \text{cor}(X_{i, \cdot}, c) \right|
$$

Starting from $i_L$, scan leftward ($i_L - 1, i_L - 2, \ldots$) and set $i_L^{\text{ext}} = i$ as long as $r_i \geq \rho_{\text{ext}}$ (default $0.6$). Stop at the first $i$ where $r_i < \rho_{\text{ext}}$; that $i$ is outside the inversion. Symmetrically extend rightward to $i_R^{\text{ext}}$.

The population-level block breakpoints are $p_{i_L^{\text{ext}}}$ (left) and $p_{i_R^{\text{ext}}}$ (right).

### 3.5 Per-carrier ancestral fragment scan (STEP_BP_02, from C01i)

Let $\mathcal{K} = S^{\text{HET}} \cup S^{\text{H2}}$ be the set of inversion carriers (assuming HOMO_2 is the minor arrangement — swap otherwise based on $|S^{\text{H1}}|$ vs $|S^{\text{H2}}|$).

For each carrier $j \in \mathcal{K}$ and each side:

**Left-side scan**: starting at $i = i_L^{\text{ext}}$ and decrementing, at each position $i$ test whether:

$$
\left| \text{cor}(X_{i, \cdot}, c) \right| \geq \rho_{\text{frag}}
$$

using only samples with non-missing genotype at SNP $i$ (default $\rho_{\text{frag}} = 0.5$, more permissive than $\rho_{\text{ext}}$ because we're making a per-sample call against a population consensus). The scan stops at the first $i$ where the correlation drops. That $i$ is the **left boundary for sample $j$**, $\ell_j$.

Analogously for the right boundary $r_j$.

The per-sample fragment is $[p_{\ell_j}, p_{r_j}]$.

### 3.6 Population-level breakpoint from fragment distribution (STEP_BP_02)

Define:

$$
B^{\text{left}} = \{ p_{\ell_j} : j \in \mathcal{K} \}
\qquad
B^{\text{right}} = \{ p_{r_j} : j \in \mathcal{K} \}
$$

Because recombinant samples have their boundaries **interior** to the true breakpoint (shorter fragments), $B^{\text{left}}$ is right-skewed (mass piles up at the true left breakpoint, with a tail extending into the inversion interior), and $B^{\text{right}}$ is left-skewed.

The **modal** breakpoint on each side is estimated by kernel density estimation on $B^{\text{left}}$ and $B^{\text{right}}$:

$$
\hat p^{\text{left}} = \arg\max_p \hat f_{B^{\text{left}}}(p)
$$

with a Gaussian kernel of bandwidth chosen by Silverman's rule on the non-tail portion (points within 1 MAD of the median, avoiding the recombinant tail biasing the bandwidth).

The **CI** on each breakpoint is bootstrap-derived. Resample carriers with replacement 1000 times, recompute the mode each time, report the 2.5th and 97.5th percentiles of the bootstrap distribution.

### 3.7 Consensus merge (STEP_BP_03)

Multiple independent sources produce breakpoint estimates (see `WORKFLOW_DIAGRAM.md`). Each source $s$ contributes a left-breakpoint estimate $\hat p_s^{\text{left}}$ and a right-breakpoint estimate $\hat p_s^{\text{right}}$ with associated weight $w_s$:

| Source | Weight $w_s$ |
|---|---|
| Ancestral-fragment mode (STEP_BP_02) | 3.0 |
| Block-boundary extension (STEP_BP_01) | 2.0 |
| C01j regime transitions | 2.0 |
| STEP40 coherence breaks | 1.0 |
| STEP41 switch pile-ups | 1.0 |
| C01l Delta_12 flanking drops | 1.0 |
| STEP37 SV cluster medians | 0.5 |

These weights reflect:
- per-carrier-voted methods at SNP resolution are most informative (weight 3)
- population-level dosage signals at SNP resolution are next (weight 2)
- window-level signals are less precise (weight 1)
- SV callers are known to miss inversions or place them tens of kb off (weight 0.5, user instruction)

The consensus left breakpoint is the weighted median:

$$
\hat p^{\text{left}}_{\text{cons}} = \text{wmedian}\left( \hat p_s^{\text{left}}, \, w_s \right)
$$

The CI is derived from the weighted median absolute deviation (wMAD) scaled by 1.4826. Methods within $2 \times \text{wMAD}$ of $\hat p^{\text{left}}_{\text{cons}}$ count as "agreeing"; agreement count is reported alongside the consensus.

### 3.8 Fallback behaviour when optional sources are missing

If only STEP_BP_01 and STEP_BP_02 outputs are available (Clair3 VCF–dependent scripts haven't run for this chromosome), the consensus degenerates to a two-method weighted average with a CI dominated by the STEP_BP_02 bootstrap CI. This is noted in the output under column `sources_available`.

### 3.9 Four-encoding robustness check (STEP_BP_05, diagnostic layer)

This step does NOT produce breakpoint estimates. It is an independent diagnostic that tests whether the sample clustering within a candidate region is robust to arbitrary allele-polarity conventions.

For each candidate region, construct four encoded versions of the dosage matrix:

1. **Minor**: raw dosage $X$ (minor-allele count convention)
2. **Major**: $2 - X$ (major-allele count convention)
3. **012**: hard discretization via threshold ($x < 0.5 \to 0$; $0.5 \le x < 1.5 \to 1$; $x \ge 1.5 \to 2$)
4. **Polarized**: apply the STEP29 L1+L2 per-marker sign to flip markers whose polarity disagrees with the homozygote-anchored consensus

For each encoding, compute pairwise sample Manhattan distance (normalized by marker count), then apply `hclust(method = "ward.D2")` and `cutree(k = 3)`. This yields four cluster assignments per sample.

**Inter-encoding agreement is measured by Adjusted Rand Index (ARI).** For every pair of encodings, compute ARI between the two cluster assignments. Report the full $4 \times 4$ agreement matrix.

Interpretation:
- **Min ARI ≥ 0.9 across all pairs**: sample clustering is robust to encoding choice. The biological structure is strong enough to survive arbitrary polarity conventions. Commit to any one encoding for downstream work.
- **Min ARI between 0.5 and 0.9**: some encodings diverge. Inspect the MDS panels. Usually either (a) the 012 discretization disagrees with the continuous encodings due to threshold effects at borderline sites, or (b) the polarized encoding disagrees with the raw encodings, indicating STEP29 is capturing or correcting a real polarity signal the raw encodings miss.
- **Min ARI < 0.5**: the structure depends strongly on encoding choice. The candidate should be treated as ambiguous or composite. No single clustering is trustworthy without further information.

**Caveat: Manhattan distance is partially polarity-invariant.** A flipped marker contributes $|x_i - x_j|$ to sample-sample distance regardless of which allele is labeled "minor." For strong 3-band signals, all four encodings give ARI ≈ 1 by construction. The method is therefore most informative for subtle structure (sub-variants within a karyotype, noisy sites, recombinant boundaries) rather than for broad 3-band separation.

This step is an optional diagnostic (enabled by setting `RUN_BP05=1` in the wrapper environment). It does not enter the consensus merge in STEP_BP_03.

---

## 4. Output semantics

### 4.1 `candidate_breakpoints_consensus.tsv`

One row per candidate. Columns:

| Column | Meaning |
|---|---|
| `candidate_id`, `chrom` | identifier |
| `final_left_bp`, `final_right_bp` | consensus breakpoint estimates (bp) |
| `left_ci_low`, `left_ci_high` | 95% CI bounds on left breakpoint |
| `right_ci_low`, `right_ci_high` | 95% CI bounds on right breakpoint |
| `left_ci_width_kb`, `right_ci_width_kb` | CI width summary |
| `n_methods_agreeing_left` | count of methods within 2×wMAD of consensus left |
| `n_methods_agreeing_right` | same for right |
| `sources_available` | comma-separated list of which upstream signals were available |
| `primary_source` | the source with the highest single-method weight that was available |
| `inversion_size_bp` | `final_right_bp - final_left_bp` |
| `refined_vs_input_shift_left_kb` | signed shift from the input CANDIDATE_TABLE `start_bp` |
| `refined_vs_input_shift_right_kb` | signed shift from the input CANDIDATE_TABLE `end_bp` |

### 4.2 `candidate_ancestral_fragments.tsv`

One row per (candidate, sample) for each inversion carrier. Columns:

| Column | Meaning |
|---|---|
| `candidate_id`, `sample` | identifier |
| `coarse_group` | HET or HOMO_2 (or HOMO_1 if HOMO_1 is minor) |
| `frag_start_bp`, `frag_end_bp` | per-sample fragment boundary positions |
| `frag_length_bp` | fragment length |
| `extended_left_snps`, `extended_right_snps` | number of SNPs scanned past block edge before correlation dropped |
| `left_cor_at_boundary`, `right_cor_at_boundary` | correlation at the stopping position |

This table is the raw per-sample vote. The distribution of `frag_start_bp` across carriers IS the breakpoint estimate's empirical distribution.

### 4.3 `candidate_breakpoints_per_method.tsv`

One row per (candidate, method, side). Columns:

| Column | Meaning |
|---|---|
| `candidate_id`, `chrom` | identifier |
| `method` | source name (e.g., `ancestral_fragments`, `block_extension`, `c01j_transitions`, etc.) |
| `side` | `left` or `right` |
| `breakpoint_bp` | that method's estimate for that side |
| `weight` | weight applied in consensus |
| `ci_low`, `ci_high` | method-specific CI if the method natively produces one; else NA |
| `distance_to_consensus_kb` | signed distance from the final consensus |
| `within_agreement_band` | boolean: is this method within 2×wMAD of consensus? |

---

## 5. Validation for LG28

Known characteristics of the LG28 inversion (from prior session work):
- Shelf interval 15.115–18.005 Mb (2.89 Mb span)
- Intermediate frequency (p ≈ 0.5), HWE-consistent three-band genotype structure
- 60/106/60 karyotype counts
- DELLY INV at 14.87–14.94 Mb (misplaced by ~200 kb — use as negative control for the SV-weight setting)
- 360 kb assembly gap at 15.33–15.70 Mb internal to the inversion

**Expected outcomes when the workflow is run on LG28**:

1. `final_left_bp` within 30 kb of 15.115 Mb
2. `final_right_bp` within 30 kb of 18.005 Mb
3. Left CI width < 100 kb; right CI width < 100 kb
4. `n_methods_agreeing_left ≥ 3` (ancestral fragments + block extension + at least one window-level method)
5. The DELLY call at 14.87–14.94 Mb should appear in `candidate_breakpoints_per_method.tsv` with `method = step37_sv`, `weight = 0.5`, and `within_agreement_band = FALSE` — demonstrating that the consensus is robust to SV-caller errors.

If outcome 5 fails (DELLY dragged the consensus left of 15.115), reduce the STEP37 weight to 0.25 in `STEP_BP_03` config.

---

## 6. Pipeline placement

The four scripts live in the **followup pipeline**, not in `phase_qc_shelf`. Rationale:

- `phase_qc_shelf` operates on raw BEAGLE per chromosome, discovering candidates.
- The followup pipeline operates on pre-computed dosage per candidate, characterizing discovered candidates.
- Breakpoint refinement is characterization, not discovery.

A shell wrapper `STEP_Q10_breakpoint_refinement.sh` chains the followup pipeline's breakpoint scripts for a single chromosome, providing phase_qc_shelf users with a one-command entry point.

---

## 7. Reference implementations and provenance

Every core algorithmic step in this methodology already had a working implementation in the user's codebase prior to this session. The new scripts are a composition and extension, not an invention:

| Step | Source function | Location |
|---|---|---|
| §3.2 informative markers | `find_informative_markers` | `STEP_C01h_recombinant_scanner.R:142–184` |
| §3.3 core block via marker correlation | `find_coseg_blocks` | `STEP_C01i_multi_inversion_decomposition.R:168–270` |
| §3.4 block boundary extension | `extend_block_boundaries` | `STEP_C01i_multi_inversion_decomposition.R:398–453` |
| §3.5 per-carrier ancestral fragment scan | `compute_ancestral_fragments` | `STEP_C01i_multi_inversion_decomposition.R:501–577` |
| §3.6 modal breakpoint from distribution | new this session | `STEP_BP_02_ancestral_fragments.R` |
| §3.7 consensus merger | new this session | `STEP_BP_03_consensus_merge.R` |
| §3.8 fallback handling | new this session | in STEP_BP_03 |
| §3.9 four-encoding comparison | original design from April sample-belonging engine, ported this session | `STEP_BP_05_four_encoding_comparison.R` |

Items §3.2–§3.5 were already operating inside the C01 multi-inversion decomposition pipeline but were *not* being used to produce population-level breakpoints with CI — they were being used for system decomposition and per-sample classification only. This session made the population-breakpoint repurposing explicit.

---

## 8. Limitations

1. **Minor-arrangement carrier count**: fragment distribution has $|\mathcal{K}|$ = n_HET + n_HOMO_2 ≈ 166 samples for LG28, which is ample for a smooth KDE. For inversions at frequency ≪ 0.5, the fragment distribution has fewer samples and the CI widens accordingly.

2. **Imputation softness**: BEAGLE dosages at low-coverage sites have inherent uncertainty that propagates to the per-sample boundary calls. A boundary that depends on a handful of low-confidence sites near the edge may be unreliable. STEP_BP_02 optionally restricts to sites with BEAGLE posterior max > 0.8 as a sensitivity check.

3. **Recombination asymmetry**: if one side of the inversion has experienced systematically more historical recombination than the other, the fragment distributions on the two sides have different shapes. The modal-estimation approach handles this naturally, but tail statistics differ between sides and should not be compared directly without normalization.

4. **Assembly-gap adjacency**: the LG28 internal gap at 15.33–15.70 Mb creates a SNP-dropout zone. Fragments from carriers with recent recombination inside the gap cannot be precisely localized and may pile up at the gap edges rather than at the recombination position. Output flags any carrier whose fragment boundary falls within 50 kb of an assembly gap for manual review.

5. **Does not detect breakpoints of inversions not already in CANDIDATE_TABLE**: this is a refinement pipeline, not a discovery pipeline. New inversion candidates must come from the upstream phase_qc_shelf or followup detection flow.

---

## 9. Future extensions (not in scope this session)

- Multi-system handling: when C01i detects ≥2 overlapping inversions in the same candidate region, produce one consensus breakpoint set per system rather than one set per region. STEP_BP_01 is compatible with this but STEP_BP_03 currently merges to one consensus; needs a system-aware pass.
- Cross-cohort validation: run the same workflow on the C. macrocephalus wild cohort (planned paper) and compare breakpoint positions at homologous inversions across species.
- Registry push: once the ancestry registry schema is confirmed, add a thin adapter that loads the three output TSVs into registry tables.

---

*This document is a living reference. Update it when thresholds change or when new evidence sources are added to the consensus merge.*
