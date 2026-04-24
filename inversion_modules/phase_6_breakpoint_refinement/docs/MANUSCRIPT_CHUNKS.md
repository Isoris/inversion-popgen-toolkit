# Manuscript Chunks — Detection Floors & Deer Mice Comparison

**Status**: draft snippets discovered in the 2026-04-20 session. Not final prose — raw material for the manuscript. Numbers are estimates based on cohort properties (226 samples, 9× coverage, p=0.5 for LG28); replace with empirical values after running the detection-floor experiment on LANTA.

---

## Where each chunk goes

| Chunk | Section | Length | Status |
|---|---|---|---|
| Detection-floor one-liner | Methods | 1–2 sentences | ready (placeholder numbers) |
| Detection-floor experiment | Results subsection | 1 figure + 3–4 sentences | needs empirical run |
| Repeat-context enrichment | Results subsection | 1 figure + 1 paragraph | needs empirical run |
| Deer mice comparison | Discussion | 1 paragraph | ready (needs final numbers) |
| Limitations paragraph | Limitations | 1 paragraph | ready |

---

## Methods chunk (short, insert into methods section on breakpoint detection)

> Inversion detection sensitivity was bounded by two factors: informative-marker density within the candidate region and the count of minor-homozygote carriers available for the per-marker signal (mean(HOMO_2) − mean(HOMO_1)). Based on the cohort properties of 226 individuals at ~9× coverage and the informative-marker rate observed in the 226-sample data, the method reliably detects inversions ≥50 kb in length with minor-allele frequency ≥15%. Below these thresholds, the core-block correlation structure is unstable, and breakpoint confidence intervals become meaningless. Empirical sensitivity was evaluated by controlled subsampling of the LG28 inversion (see Results § [sensitivity]).

---

## Results chunk — detection-floor experiment (needs empirical run)

**Proposed figure**: two panels.
- **Panel A**: detection probability vs candidate size. X-axis = simulated inversion size (50 kb, 100 kb, 200 kb, 500 kb, 1 Mb, 2 Mb). Y-axis = fraction of runs (out of 20 replicates) where breakpoint was recovered within 30 kb.
- **Panel B**: detection probability vs MAF. X-axis = MAF (0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30). Y-axis = fraction recovered.

**Proposed prose** (placeholder; replace with empirical numbers):

> We assessed detection sensitivity by progressively subsampling the LG28 inversion. Size sensitivity was evaluated by truncating the candidate region from its full 2.89 Mb to 1 Mb, 500 kb, 200 kb, 100 kb, and 50 kb. MAF sensitivity was evaluated by randomly subsetting carriers from the full 226-sample cohort down to 60, 40, 30, 20, and 10 minor-homozygote carriers plus proportional heterozygotes. For each condition, 20 replicate subsamples were drawn. Detection was scored as recovery of the known breakpoint within 30 kb. The method achieved >90% detection probability down to 50 kb in size and down to ~15% MAF; below these floors, detection degraded sharply (Fig. [SX]).

---

## Results chunk — repeat-context enrichment (needs empirical run)

**Proposed figure**: one panel.
- X-axis: distance from breakpoint (kb), from −500 to +500.
- Y-axis: density of each repeat class (segmental duplications, LINEs, SINEs, simple repeats).
- Overlay: genome-wide random-region null density as a grey ribbon.
- Kolmogorov–Smirnov test for SD enrichment against null, p-value in figure legend.

**Proposed prose** (placeholder; replace with empirical numbers):

> Across the [N] inversion breakpoints refined in this study, breakpoint regions (±500 bp) were enriched for segmental duplications compared to size-matched random intervals (Kolmogorov–Smirnov test: P < [X]). A subset of breakpoints harboured tandemly structured SDs (Fig. [SY]), consistent with a role for ectopic recombination between paralogous sequences. The inversion breakpoint landscape in *C. gariepinus* thus parallels observations in other vertebrate systems, where segmental duplications have been repeatedly implicated as fragile sites (e.g., [deer mice citation]).

---

## Discussion chunk — comparison to deer mice and other studies

> Our detection floor of ~50 kb compares favourably to approaches based solely on read-pair or split-read evidence, which typically require inversions ≥1 Mb to accumulate sufficient mapping anomalies ([deer mice citation: Peichel et al.]). This advantage derives from the use of linkage-disequilibrium signal encoded in BEAGLE dosage rather than direct read evidence: an inversion as small as 50 kb still contains ~100 informative markers under typical segregating variation in the 226-sample cohort, which is sufficient for stable correlation-based core-block detection even when per-sample read-level evidence is too sparse to call. The MAF floor of ~15% is set not by the method but by the cohort size: at MAF 15%, the expected count of minor homozygotes under HWE is 226 × 0.15² ≈ 5, which is at the floor for reliable per-marker homozygote-mean estimation. Larger cohorts would push this floor lower without requiring methodological changes.

---

## Limitations paragraph

> Several classes of events are undetectable by the present method. (i) Inversions with MAF < 15% produce fewer than five minor homozygotes in this cohort, insufficient for stable core-block detection. (ii) Inversions shorter than ~50 kb fall below the informative-marker count threshold for stable correlation structure. (iii) Gene conversion tracts and double crossover events within the inversion interior are collapsed into the same signal as the primary breakpoints; distinguishing them would require phased long-read data. (iv) Segmental-duplication-rich regions produce mapping ambiguity that can inflate the breakpoint confidence interval beyond what the dosage signal predicts; such regions are flagged in the output but not corrected. (v) We did not compute a genome-wide split-read or discordant-pair null distribution, because baseline rates in repetitive and low-complexity regions are high enough that any per-position enrichment test would have a fat null tail dominated by mapping artifacts rather than real rearrangements; our breakpoint calls instead rely on orthogonal evidence (dosage signal, DELLY/Manta concordance, regime-compatibility continuity).

---

## Unresolved numbers (fill in from LANTA runs)

- [ ] **Total N inversions called** in the 28-chromosome cohort (deer mice reports 21; comparable number for *C. gariepinus* pending)
- [ ] **Size distribution** of called inversions (mean, median, range)
- [ ] **MAF distribution** of called inversions
- [ ] **Empirical detection floor** from sensitivity experiment (expected: ~50 kb, ~15% MAF — confirm)
- [ ] **SD enrichment p-value** (Kolmogorov–Smirnov)
- [ ] **SD-harbouring breakpoint fraction** (deer mice: 50% in top 90th percentile)
- [ ] **Proportion of breakpoints with tandem-structured SDs** vs interspersed vs neither

---

## Future script ideas (not for next session — wiring first)

Two future scripts that would produce the Results figures above. NEITHER is to be built in the wiring chat. Recorded here so they aren't forgotten.

### Future script — sensitivity sweep
- Input: full LG28 dosage + sites, list of carriers
- Procedure: loop over (size truncation × MAF subsample × replicate), invoke scripts 01–03 on each subsample, record whether `final_left_bp`/`final_right_bp` fall within 30 kb of ground truth
- Output: `sensitivity_detection_floor.tsv` with columns (size_kb, maf, replicate_id, detected, left_error_kb, right_error_kb) + a 2-panel figure
- Estimated size: ~150 lines R

### Future script — breakpoint repeat enrichment
- Input: final breakpoint list from `breakpoints_catalog.tsv`, RepeatMasker BED, assembly fasta index
- Procedure: for each breakpoint, compute repeat-class density in ±500 bp, compare to N random matched-size intervals from the same chromosome. Kolmogorov–Smirnov test per class.
- Output: `breakpoint_repeat_enrichment.tsv` + ridge density figure
- Estimated size: ~100 lines R

Both fit naturally as Module 2 post-processing (they consume Module 2 outputs). Alternatively they could be a new Module 4 (manuscript statistics), separate from the core pipeline.
