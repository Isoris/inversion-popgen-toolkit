# Session Audit — 2026-04-20

Full record of what was decided, when, and why. Written for the next Claude session so it can understand the state of the codebase without re-asking the user.

---

## Session timeline (condensed)

### Part 1 — Heatmap classifier rabbit hole

User asked for a classifier that decomposes composite inversion heatmaps into architecture classes (SIMPLE_STRONG, COMPOSITE_INTERNAL, etc.) and implements "separate first, correct second" polarity harmonization.

**Built**:
- `STEP40_architecture_classifier.R` — 10-feature classifier with 7 architecture classes + 4 tiers
- `STEP41_separate_first_heatmaps.R` — per-subgroup locally-harmonized heatmaps
- `q08b_shelf_heatmap_architecture.R` — upgraded shelf heatmap with architecture overlays + D1-D3 separate-first panels
- `STEP_Q08b_shelf_heatmap_architecture.sh` — bash driver

**Self-test**: Python port of classifier logic. 3/5 scenarios passed on first run (SIMPLE_STRONG, COMPOSITE_INTERNAL, COMPOSITE_OVERLAP). Two failures traced to test-harness RNG reuse issues, not classifier bugs.

**User feedback that redirected the work**: User pointed out "family_ld_correlation" was a mislabel — I had named a sample-partition similarity check "LD" which was wrong. Renamed in final docs (the R code itself still has the old name; this is a follow-up for next session).

### Part 2 — User challenges on heatmap interpretation

User pushed back several times on my reasoning about composite heatmaps:

1. "Why family LD? Heatmaps have nothing to do with LD." → I corrected: the feature was a sample-space similarity, poorly named.

2. "Where do polarization errors come from?" → I framed as markers disagreeing with karyotype expectation. User pushed: "how can samples disagree with themselves?" → reframed correctly: a cell (sample i, marker k) deviates from the group-level expectation at marker k; the unit of disagreement is cell-level, not sample-self-level.

3. Looking at images: User observed "some polarization-error columns stop at karyotype boundaries, others cross them, not constant, not periodic." I proposed Case 1/2/3 classification framework. User rejected: "its strange, its mixed" — suggested maybe many crossovers shared across samples, which would require deep ancestry. Reasonable skepticism.

4. User proposed: "maybe the problem is we didn't order the y axis top-down." → I checked. User was right. STEP28/STEP29/STEP35/q08 all use `setorder(ri, coarse_group_refined, u)` without within-karyotype clustering. `cluster_rows=FALSE` in all ComplexHeatmap calls.

### Part 3 — The u/v rotation "rediscovery"

User: "I coded the PCA reprojection... the u v transform putting Hom1 as 0... why I'm so genius... search for the u v transform."

I found it in their own codebase:
- `STEP20b::compute_anchor_features` — Het-centered, Hom1→Hom2 axis rotation
- `STEP21::rotate_to_stripe_axis` — applies the rotation per candidate, writes `candidate_pca_rotated.tsv`

My earlier claim that "u is saturated within each karyotype" was wrong. u IS a meaningful within-karyotype ordering. I apologized and corrected.

**Real conclusion**: u is correct for simple candidates. Composite candidates (LG28) have substructure orthogonal to u that u-sorting cannot expose. Solution: u primary + secondary tiebreaker within each karyotype (hclust, PC2, residual PC1, or error_profile).

**Built**:
- `within_group_ordering.R` — helper with 4 secondary-axis options, u always primary
- `STEP36_ordering_comparison_heatmaps.R` — 3-version diagnostic (V1 raw PC1 / V2 u-only / V3 u+hclust)
- Patches for STEP28, STEP35 to use the helper
- Upgraded STEP41 and q08b to use the helper

User asked for 3 layouts: individual full-size PDFs + 3-page combined + side-by-side single-page. All three provided.

### Part 4 — Breakpoint question

User asked: "How to define breakpoints for inversions? Fst too coarse. Dosage-based? Marker-level?"

I read 10 additional scripts user provided:
- STEP32 (multiscale windows), STEP33 (within-stripe analysis), STEP34 (anchor resolution), STEP35 (heatmaps — already seen), STEP36 (Clair3 local signatures — different from my STEP36 above; name collision), STEP37 (SV breakpoints), STEP38 (composite figure), STEP39 (SV visualization), STEP40 (internal coherence — different from my STEP40 above; name collision), STEP41 (membership trajectories — different from my STEP41 above; name collision)

Then user provided 5 more: STEP_C01h (recombinant scanner), STEP_C01i (multi-inversion decomposition), STEP_C01j (regime compatibility), STEP_C01l (local structure segments), STEP_C01m (distance concordance).

Found the key functions: `find_informative_markers` (C01h), `find_coseg_blocks` (C01i), `extend_block_boundaries` (C01i), `compute_ancestral_fragments` (C01i), `segment_regimes` (C01j), `define_segments` (C01l), `compute_concordance_at_distance` (C01m), plus `detect_changepoints` (C01h, binary segmentation with double-crossover extension).

**User confirmed**: "use mostly the dosage stuff. Put sv on top if you have energy."

**User confirmed**: the complete workflow was already pre-shot in their codebase — this session's job is wiring, not invention.

User mentioned: "they plug into the ancestry registry... there was something about ancestry but I forgot." I searched the scripts and found `compute_ancestral_fragments` in C01i — this is the ancestry hook. Per-carrier fragment lengths ARE the ancestry signal.

### Part 5 — Final workflow design

User instruction: "Create the full handoff and the tallbar and the session summary and the diagram that u did and the detailed document md file for the methodology. And also create the audit of this session so the next chat code it all. But if you have energy do it all in this chat."

Decided to do both documentation AND code. Started with docs (HANDOFF, WORKFLOW_DIAGRAM, METHODOLOGY, SESSION_AUDIT) because those are the highest-priority artifacts if I run out of time/tokens.

---

## Key decisions made and their rationale

### D1. No Fst-based breakpoint step

User was explicit: Fst at 50 kb is too coarse. Dosage is the signal. Do not add Fst refinement.

### D2. SV calls weighted 0.5, used as bottom-layer only

User: "37 39 the SVs we see it doesn't really work." DELLY called LG28 INV at 14.87–14.94 Mb, ~200 kb off. Weight 0.5 so it can't drag consensus when other methods disagree. Plotted with low alpha in diagnostic figure.

### D3. Per-carrier ancestral fragments as the primary signal (weight 3.0)

User's own code in `C01i::compute_ancestral_fragments` already does the per-sample scan. It had not been lifted to population-level breakpoints with CI; this session does that.

Rationale for weight 3.0: only per-carrier-voted method with native CI. Other methods give point estimates.

### D4. Four scripts + one wrapper, not more

User said: "Yes thats exactly it" to the four-script proposal. Do not add scripts. Keeps the pipeline graph manageable.

### D5. Scripts live in followup pipeline, not phase_qc_shelf

phase_qc_shelf discovers candidates. Followup characterizes them. Breakpoint refinement is characterization. The shell wrapper `STEP_Q10_breakpoint_refinement.sh` bridges phase_qc_shelf users into the followup pipeline for one chromosome.

### D6. Three separate output TSVs, not a registry schema

User could not remember the ancestry registry schema. Safe choice: three separate TSVs (per-candidate, per-sample-per-candidate, per-method-per-candidate). Whatever the registry wants, it can pull from these.

### D7. KDE modal + bootstrap CI for fragment-distribution breakpoint

Not mean (recombinant tails bias it). Not median (modal is more appropriate for a right-skewed distribution with a sharp peak). KDE mode + bootstrap is the right combination for robust point estimate with uncertainty.

### D8. Weighted median, not weighted mean, for consensus across methods

Median is robust to outliers (e.g., a single bad SV call 200 kb off). The user specifically wanted SV calls not to be able to drag the consensus.

### D9. N_methods_agreeing as a quality flag

Rather than producing a "confidence score" that collapses information, expose the raw count of agreeing methods. User can filter downstream.

### D10. Do not implement multi-system breakpoint calling yet

If C01i detects 2+ overlapping inversions in one candidate, current scripts merge to one consensus. Multi-system handling is deferred — noted in METHODOLOGY §9.

### D11. Do not attempt runtime testing in Claude's sandbox

No R available in this sandbox. All scripts pass static brace-paren balance checks. Runtime testing deferred to user's LANTA cluster run.

### D12. Write R scripts before completing every optional feature

Completeness of the core path > perfection of each script. STEP_BP_01–04 will do the minimum correctly; edge cases flagged in comments for follow-up.

### D13. Add BP_05 four-encoding diagnostic as an OPTIONAL layer (not a consensus source)

Late in the session, the user recalled a previous architectural choice from the April sample-belonging engine: computing sample × sample distance under four encodings (minor dosage, major dosage, 012 discrete, and optionally STEP29-polarized) to test whether sample clustering is robust to arbitrary allele-polarity conventions.

This was ported as `STEP_BP_05_four_encoding_comparison.R`. Key choices:

- **NOT added to the BP_03 consensus**. BP_05 produces cluster assignments, not breakpoints. It answers a different question (is the structure robust to encoding choice?) than the breakpoint workflow (where are the edges?).
- **Optional by default**. Controlled by `RUN_BP05=1` environment variable in the wrapper. Skipped by default to avoid multiplying output files.
- **Dual-mode**: supports standalone invocation with explicit file paths (Mode A) OR config-driven invocation that iterates candidates (Mode B), matching the BP_01–04 call pattern.
- **Synthetic validation during development revealed an important caveat**: Manhattan sample-sample distance is partially polarity-invariant by construction. A flipped marker contributes `|x_i - x_j|` to sample distance regardless of which allele is labeled "minor." For strong 3-band signals, all four encodings give ARI ≈ 1 by construction. BP_05 is therefore most informative for subtle structure (sub-variants within a karyotype) rather than for broad 3-band separation. This is documented in both the methodology and the workflow diagram.

### D14. Acknowledge what BP_05 doesn't solve

The user originally reached for the 4-encoding approach to handle the "how do we compare genotype vectors across samples without arbitrary allele choice?" question. The honest finding, after implementation and synthetic testing: for the primary 3-band structure detection, the question is largely moot because Manhattan distance on dosage is already polarity-robust. For sub-karyotype structure (BB_1 vs BB_2), the 4-encoding comparison IS useful — which aligns with the user's original April motivation, where the engine lived inside a "within-stripe sub-cluster detection" workflow (not breakpoint refinement).

This means:
- For breakpoint refinement: BP_05 doesn't add much; the primary signal from BP_01/BP_02 is already polarity-robust.
- For composite-candidate sub-structure (e.g., LG28): BP_05 is a real diagnostic; it can distinguish real sub-structure from encoding artifacts.

Next session should decide whether to promote BP_05 into a separate diagnostic module for composite candidates specifically, or leave it as an optional layer in the breakpoint workflow.

### D15. Integration with the unified_ancestry module is NOT done this session

The user mentioned wanting to wire the regime / ancestry analysis into their existing `unified_ancestry/` module (with Engine A = NGSadmix full survey, Engine B = `instant_q.cpp` fixed-F EM, and `region_popstats.c` for Fst/dXY/theta). This was discussed but not implemented. The unified_ancestry module was described in past conversations but the actual current `region_stats_dispatcher.R`, `instant_q.R`, and the latest C01j script were not uploaded in this session.

Recommendation for next session: upload those three files, then decide whether C01j (or BP_05) should consume `get_region_stats()` via the dispatcher rather than computing its own distances. This is the proper ancestry-engine integration.

### D16. Rolling seed-and-extend regime analysis — IDEA RECORDED, NOT IMPLEMENTED

At the very end of the session the user raised an idea that is architecturally more important than the breakpoint workflow: instead of asking "where are the edges" (BP_01–04), slide along the chromosome tracking compatibility-group composition, and record regime transitions along with their *shape* (sharp flip vs gradual vs bubble).

This is strictly more informative than the breakpoint workflow: it contains breakpoints as a derived output AND captures recombination events, family LD transitions, and gene conversion bubbles that the breakpoint workflow collapses to a single position.

**Decision**: do NOT code this tonight. Record it in HANDOFF.md §"Priority next direction" with the five questions that must be decided before coding (state representation, block definition, transition-shape characterization, output schema, computational budget). The next session should address those questions before writing any code.

**Why not tonight**: the user raised this while tired; coding it without rested decisions on the five questions would produce a module that needs a rewrite. Estimated real size: 500–1000 lines R or 300–500 lines C + R wrapper. Multi-day task.

---

## End-of-session state (post-BP_05)

Files produced or updated this session:

```
R/
  STEP_BP_01_dosage_signal.R                (new)
  STEP_BP_02_ancestral_fragments.R          (new)
  STEP_BP_03_consensus_merge.R              (new)
  STEP_BP_04_diagnostic_figure.R            (new)
  STEP_BP_05_four_encoding_comparison.R     (new, dual-mode)
  STEP_C01j_stream_graph.R                  (new, standalone regime visualization)
  STEP_Q10_breakpoint_refinement.sh         (new, chains BP_01-04 + optional BP_05)

docs/
  HANDOFF.md                                (updated — reflects BP_05 + open items)
  METHODOLOGY.md                            (updated — §3.9 added)
  WORKFLOW_DIAGRAM.md                       (breakpoint workflow, unchanged)
  WORKFLOW_DIAGRAM_BP05.md                  (new — 4-encoding module tall-bar)
  SESSION_AUDIT.md                          (this file)
```

All R scripts pass brace/paren balance checks. Python synthetic tests validate:
- Cluster matching across windows (stream graph): 2<->1 local relabeling correctly collapsed to same global ID
- 4-encoding equivalence (BP_05): all four encodings give ARI = 1.0 on clean 3-band signal even with 50% of markers polarity-flipped

No R runtime tests performed (no R in sandbox).

---

## Things I got wrong this session and had to correct

1. Named a feature "family_ld_correlation" when it was actually a kinship-partition similarity check. Not LD. User caught this.

2. Said "u is saturated within karyotype" before checking. It isn't. User-corrected by pointing to their own STEP20b/STEP21 code.

3. Kept inventing classifiers and features instead of reading the user's existing code first. User corrected: "read first, propose second."

4. Proposed Case 1/2/3 for polarization error columns based on pattern-matching an image, without a testable definition. User correctly pushed back with "it's mixed, it's strange."

5. When user asked "why some stripes stop at boundaries, some cross" I explained but my initial framing was still not backed by the data — I should have proposed a diagnostic (pull numeric values, not interpret an image) earlier.

---

## Things the next session should double-check

1. Confirm the ancestry registry schema with the user before writing registry adapters.
2. Runtime-test `STEP_BP_01` through `STEP_BP_04` on LG28 — expected outputs documented in METHODOLOGY.md §5.
3. If `n_methods_agreeing_left` < 3 for most candidates, lower the `2×wMAD` agreement band threshold or raise the weight of weaker methods.
4. The `cor_min` threshold (0.5) in the per-carrier fragment scan is conservative. If CIs come out wider than 100 kb on LG28, consider lowering to 0.45.
5. STEP28 and STEP35 patches from the earlier session are documented in `STEP28_patch.txt` and `STEP35_patch.txt` — these were never applied to the user's actual STEP28/STEP35 files. The within_group_ordering.R helper is stand-alone.
6. There are three naming collisions in the user's codebase with my work from earlier turns:
   - User has `STEP36_candidate_clair3_local_signatures.R`; my new `STEP36_ordering_comparison_heatmaps.R` uses the same number. **Rename mine to STEP42 or similar to avoid collision.**
   - User has `STEP40_candidate_internal_coherence.R`; my `STEP40_architecture_classifier.R` collides. **Rename mine to STEP43 or similar.**
   - User has `STEP41_candidate_membership_trajectories.R`; my `STEP41_separate_first_heatmaps.R` collides. **Rename mine to STEP44 or similar.**
7. The new breakpoint scripts `STEP_BP_01` through `STEP_BP_04` use `BP_` prefix to avoid further numeric collisions. Good.

---

## Late-session additions (beyond D16)

### D17. Panel D replica (script 07) added

User uploaded an image showing a publication-style multi-panel inversion figure (Chr12 example with genes, regional PCA, sample genotype heatmap, breakpoint evidence panel D, gene/repeat context). Asked to reproduce panel D for his own data.

Built `07_breakpoint_evidence_pileup.py` — per-sample stacked read evidence at both breakpoints. Ports the MODULE_5A2 / C01f pysam extraction logic. Dual-mode: (A) reads pre-extracted evidence TSV; (B) extracts fresh from BAMs via pysam. Samples ordered INV → HET → REF, labels optionally anonymized as sample1/sample2/sample3.

**Runtime-tested in the sandbox**: Python synthetic evidence generation for 30 samples → 5319 evidence rows → script renders PDF + PNG successfully. Output saved as `docs/example_panel_D_output.png` for reference.

**Known limitation flagged at D17**: the assembly-alignment track is a schematic placeholder. User has only 9× short-read Illumina coverage per-sample, so real assembly at per-sample level isn't possible. Either drop this track or leave schematic — don't pretend to have data that doesn't exist.

**Open schema question flagged at D17**: script 07 expects long-format per-read-event TSV (`sample, bp_side, read_type, read_pos, read_end, ...`) but user's C01f STEP02 produces wide-format per-sample summary (`sample, group, support, pe, sr, n_discordant_FF, ...`). Mode A won't work against current C01f output. Solutions: (A) just use Mode B = fresh pysam; (B) modify C01f STEP02 to ALSO write a long-format TSV. Deferred to next session.

### D18. Three-module framing (user explicit)

User explicitly reframed the session's output as three modules at the end:
1. **Module 1** = marker dosage heatmaps (the heatmap ordering upgrade in `_archive/`)
2. **Module 2** = breakpoint finding + visualization (scripts 01–04 + 05 + 07)
3. **Module 3** = regime / ancestry analysis (script 06 + C01j + rolling regime idea + unified_ancestry integration)

Modules share an L1/L2 polarity layer (STEP29 machinery). This framing is correct and should be preserved in the next session's reorganization. The existing PIPELINE_DIAGRAM.md mostly covers Modules 2 + 3's visualization; Module 1 is scaffolded in the archive but not polished.

### D19. Next-chat mandate: WIRING, not new code

User's explicit direction for next chat: *"wire everything to the registries and api so it collects all genome locations and help each script to work together so our results is all automatic bc otherwise its too tired."*

Translation: no more module-building. Connect existing modules to:
- `results_registry` (read inputs by key, write outputs with provenance + sha256)
- `region_stats_dispatcher.R` (Module 3 consumes `get_region_stats()` instead of computing own distances)
- A top-level genome-location collector that walks all candidates × all chromosomes

**The trap to avoid**: "let's improve X while wiring." No. Every script stays as-is. Wiring is mechanical and must finish in one session. Improvements go in future sessions.

Files the user must upload at start of next chat (listed in HANDOFF.md): `DATABASE_DESIGN.md`, `region_stats_dispatcher.R`, `instant_q.R`, current C01j, one example C01f STEP02 output TSV, one example config file.

### Additional mistakes corrected in late session

6. In the 4-encoding script design I initially claimed STEP29 polarization was "the solution we used before" when the user had actually designed a different multi-encoding comparison in April. User pushed back; I searched past conversations and corrected — user's original solution was the 3-encoding Layer B2 (minor / major / 012), not STEP29. The 4th polarized encoding is NEW, not a re-derivation.

7. Suggested "code it all" routes multiple times when user was tired. Should have pushed back harder earlier. Eventually did push back and wrote a self-contained reference implementation instead of guessing at pipeline internals. User flagged this too: "but this idea it's been months i have it in my head."

8. Stopped mid-finalization on one turn when user asked for tarball; produced individual files only. User had to ask "where is the gz." Next session: always generate the tarball if final-bundle is requested.

### D20. Manuscript chunks document created for deer mice comparison

User showed a screenshot from the deer mice inversions paper (Peichel et al.) — the breakpoint SD-enrichment result and the size/MAF floors ("conservative estimate because our approach was limited to identifying inversions >1 Mb in length with a minimum allele frequency of ~10%"). User asked: "for us how to know what is our size in precision and what is our MAF?"

Worked through the theoretical detection floors for the 226-sample cohort at ~9× coverage:

- **Size floor ≈ 50 kb** (limited by informative-marker density; needs ~15–20 informative markers inside the inversion; informative markers occur at ~300–600 bp spacing)
- **MAF floor ≈ 15%** (limited by minor-homozygote count; at MAF 15%, N·p² ≈ 5 which is the floor for stable per-marker homozygote-mean estimation)

Compared to deer mice floors (~1 Mb / ~10%): dosage method reaches smaller inversions than their read-pair-only approach; MAF floor is set by cohort size not methodology.

User then asked about genome-wide null distribution for split-read / discordant-pair evidence at inversion positions. Honest response: probably not worth it for LG28 specifically (dosage is already stronger evidence) and noise-dominated in general (repeat regions inflate the fat tail). Recommended instead: (1) sensitivity sweep via subsampling LG28, (2) breakpoint repeat-context enrichment vs random intervals. Skipped genome-wide null distribution.

User asked where this material goes in the manuscript. Structured the split:

- **Methods**: short paragraph stating the detection floor as a number
- **Results**: sensitivity figure (size × MAF) and repeat-enrichment figure, each as a short subsection
- **Discussion**: comparison paragraph to deer mice
- **Limitations**: paragraph listing undetectable event classes

Created `docs/MANUSCRIPT_CHUNKS.md` containing:
- Draft paragraph text for each of the four manuscript locations above
- Proposed figure layouts (sensitivity sweep = 2-panel; repeat enrichment = ridge density + KS test)
- List of unresolved numbers to fill in from LANTA runs
- Future-script sketches for the two missing analyses (sensitivity sweep, repeat enrichment) — NOT for next chat, recorded for a later session

**Decision**: these are manuscript-augmentation scripts, not pipeline-integration. They fit as Module 2 post-processing OR as a possible new Module 4 (manuscript statistics). Next chat must NOT start building them — wiring comes first.

---

## File naming convention used

- User's existing scripts: `STEP<NN>_<n>.R` and `STEP_C01<letter>_<n>.R`
- New breakpoint scripts (this session): originally `STEP_BP_<NN>_<n>.R`, renamed to `<NN>_<n>.R` for the final bundle (01, 02, 03, 04, 05, 06, 07)
- New heatmap scripts (earlier this session): originally named `STEP36`, `STEP40`, `STEP41` — **collide with user scripts, must be renamed next session** (kept in `_archive/heatmap_ordering_upgrade/`)

---

## Files produced this session (current bundle layout)

```
breakpoint_pipeline/
  README.md
  PIPELINE_DIAGRAM.md
  run_pipeline.sh
  01_dosage_signal.R                     Module 2 core
  02_ancestral_fragments.R               Module 2 core
  03_consensus_merge.R                   Module 2 core
  04_diagnostic_figure.R                 Module 2 core
  05_four_encoding_diagnostic.R          Module 2 optional
  06_regime_stream_graph.R               Module 3 visualization
  07_breakpoint_evidence_pileup.py       Module 2 panel D (runtime-tested)
  test_07_synthetic.py                   smoke test for 07
  docs/
    HANDOFF.md
    METHODOLOGY.md
    SESSION_AUDIT.md                     (this file)
    MANUSCRIPT_CHUNKS.md                 draft Methods/Results/Discussion/Limitations paragraphs
    example_panel_D_output.png
  _archive/
    README.md
    heatmap_ordering_upgrade/            Module 1 scaffolding
    architecture_classifier/             abandoned
    original_workflow_diagrams/          superseded
```

All R scripts: balance depth=0. Python scripts: `ast.parse` clean. Wrapper: `bash -n` clean.

Python synthetic tests validated this session:
- Cluster matching across windows (script 06): 2↔1 local relabeling correctly collapsed to same global ID
- 4-encoding equivalence (script 05): all four encodings give ARI = 1.0 on clean 3-band signal even with 50% markers polarity-flipped
- Panel D pileup (script 07): 30-sample render succeeded; INV→HET→REF visual ordering correct; breakpoint red-dashed-line visible; coverage dip visible; sample labels anonymized as sample1..sample30

No R runtime tests. Runtime validation deferred to LANTA.

---

## End of audit
