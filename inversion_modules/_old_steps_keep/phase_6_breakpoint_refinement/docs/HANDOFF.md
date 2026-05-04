# Session Handoff — Inversion Analysis Workflow (three modules)

**Date**: 2026-04-20
**Project**: `MS_Inversions_North_african_catfish`
**Cluster**: LANTA
**Cohort**: 226 samples, pure *C. gariepinus* hatchery broodstock, ~9× coverage
**Target inversion**: LG28 at 15.115–18.005 Mb (≈2.89 Mb, p≈0.5, 60/106/60 karyotype counts)

---

## Purpose of this handoff

This document tells the next Claude session exactly where work stopped, what decisions were made, what was built vs designed, and how to continue without re-deriving context. If you are that next session: **read this file first**, then `docs/METHODOLOGY.md`, then `docs/SESSION_AUDIT.md`. Do not re-explain the architecture to the user — he knows it, he designed most of it.

---

## One-line summary

Built three connected modules for inversion analysis on the 226-sample catfish cohort. Modules have been designed and scaffolded; next chat wires them into the results registry and ancestry engine for automatic end-to-end runs.

---

## The three modules

The user explicitly framed the work as three modules at the end of the session. Respect this framing.

### Module 1 — Marker dosage heatmaps
Renders sample × marker dosage matrices with principled row ordering and L2-block polarity harmonization. Lives in `_archive/heatmap_ordering_upgrade/`:
- `within_group_ordering.R` — shared helper for row ordering (primary = u coordinate from STEP21 rotated PCA, secondary = configurable: hclust, v, residual_pc1, error_profile, or none)
- `STEP36_ordering_comparison_heatmaps.R` — V1/V2/V3 diagnostic comparison
- `STEP41_separate_first_heatmaps.R` — upgraded to use the helper
- `q08b_shelf_heatmap_architecture.R` — upgraded to use the helper
- `STEP_Q08b_shelf_heatmap_architecture.sh` — existing wrapper
- `STEP28_patch.txt` / `STEP35_patch.txt` — notes-only patches to apply to user's existing scripts

Not yet packaged as a clean standalone module. The archive is scaffolding, not polished.

**Known issue**: the scripts are named STEP36 / STEP40 / STEP41 but those numbers already exist in the user's codebase doing different things. Must be renumbered (likely STEP42/43/44) before committing.

### Module 2 — Breakpoint finding + visualization
The main output of this session. Lives in the top-level of `breakpoint_pipeline/`.

Core (always run in order):
- `01_dosage_signal.R` — informative-marker detection (|delta| ≥ 1.0), core block via marker correlation, extension via correlation walk
- `02_ancestral_fragments.R` — per-carrier fragment scan, modal position, bootstrap CI. This is the primary breakpoint signal.
- `03_consensus_merge.R` — weighted merger across sources; produces `candidate_breakpoints_consensus.tsv` (the main answer)
- `04_diagnostic_figure.R` — multi-track PDF per candidate

Optional diagnostics:
- `05_four_encoding_diagnostic.R` — four-encoding robustness (minor / major / 012 / STEP29-polarized) with ARI agreement matrix. Dual-mode (standalone + config-driven). NOT in consensus. Enable via `RUN_05=1`.
- `07_breakpoint_evidence_pileup.py` — Python replica of panel D (per-sample stacked read evidence at both breakpoints). Runtime-tested in sandbox with synthetic data; output example at `docs/example_panel_D_output.png`.
- `test_07_synthetic.py` — smoke test for script 07 without needing real BAMs

Wrapper: `run_pipeline.sh` chains 01→02→03→04, adds 05 if `RUN_05=1`.

### Module 3 — Regime / ancestry analysis (partially designed, not yet wired)
Chromosome-wide rather than candidate-centric. Captures regime shape, not just edges.

What exists:
- `06_regime_stream_graph.R` — standalone stacked-area visualization consuming C01j output (regime_memberships.tsv.gz + regime_windows.tsv.gz). Verified with Python synthetic test: cluster-matching across windows correctly collapses local 2↔1 relabeling to same global ID.
- User's existing C01j regime compatibility engine (not in this bundle — lives in user's codebase)

What is designed but not coded:
- Rolling seed-and-extend regime analysis (see §"Deferred design idea" below)
- Integration with user's `unified_ancestry/` module:
  - Engine A = NGSadmix wrapper
  - Engine B = `instant_q.cpp` fixed-F EM (clean-room Skotte 2013)
  - `region_popstats.c` for Hudson Fst / dXY / theta / MI
  - `region_stats_dispatcher.R` gateway that Module 3 should consume via `get_region_stats(what=c("Q"))` instead of computing its own distances

---

## What the shared polarity layer does

Both Module 1 and Module 2 use the same L1/L2 polarity machinery (originally STEP29):

1. Per-marker L1 sign: `sign(mean(HOMO_2) - mean(HOMO_1))` anchored on the homozygotes only (anti-circular)
2. Block-level L2 consensus: smooth L1 signs into blocks of ≥5 consecutive same-sign markers
3. Apply the L2 sign uniformly to all samples at each marker: `dosage_polarized = (dosage - 1) × L2_sign + 1`

After polarization, the canonical dosage scale is: HOMO_1 ≈ 0, HET ≈ 1, HOMO_2 ≈ 2 at every marker.

**Module 1 uses this to build canonically oriented heatmaps.**
**Module 2 uses this:**
- As the 4th encoding in script 05 (where raw minor / raw major / 012 are the other three)
- Potentially as the preferred input to script 02's per-carrier walk (not wired yet — flagged in METHODOLOGY)

The user's insight at the end of the session, captured for next time: "the L1/L2 boundaries work for coordinates too" — meaning the same block-consensus machinery that orders heatmaps can inform breakpoint refinement, because a genuine biological boundary should show up as both a polarity transition AND a coherence break.

---

## Priority next direction — WIRING, not new code

The user's explicit direction for the next chat: **"wire everything to the registries and api so it collects all genome locations and help each script to work together so our results is all automatic bc otherwise its too tired."**

Translation: stop building modules. Connect them to the existing infrastructure. The modules are designed enough.

### Specifically the next chat should

1. **Registry adapter**. Each of the 7 scripts gets a thin registry layer: read inputs by key from `results_registry`, write outputs with provenance + sha256. Stop using loose TSVs in candidate folders as the primary interface.

2. **Genome-location collector / driver**. A top-level invoker that walks all candidates across all 28 chromosomes, runs the modules in order, aggregates results to one genome-wide table. One row per candidate with every evidence column populated.

3. **Cross-module API (remove hardcoded paths)**. Script 02 needs STEP29 polarity → asks the registry. Script 05 needs the same polarity → asks the registry. Script 07 needs STEP21 karyotype groups → asks the registry. Stop duplicating path config across scripts.

4. **Module 3 plugs into unified_ancestry**. C01j (and by extension the rolling regime analysis if implemented) consumes `get_region_stats(what=c("Q"))` via the dispatcher rather than computing its own Q-matrices. Engine B becomes the Q source.

### Files the user must upload at the start of the next chat

Without these, the next chat will guess and waste time. With them, wiring is mechanical.

- `DATABASE_DESIGN.md` — current registry schema (source of truth)
- Actual `results_registry` schema files (whatever SQL/TSV format)
- `region_stats_dispatcher.R` — current version
- `instant_q.R` (or the cpp + R wrapper pair)
- Current latest `C01j` — the regime compatibility engine
- **One example C01f STEP02 output TSV** — to resolve the per-sample evidence schema question (see below)
- One example `config_inversion_followup.R` — for path/variable naming conventions

### The trap to avoid in the next chat

"Let's also improve X while we're at it." No. The next chat's job is wiring, not improving. If something is wrong with a script, note it for a later chat. Every script stays as-is; we hook it to the registry and the dispatcher. Otherwise wiring never finishes.

---

## Deferred design idea — rolling seed-and-extend regime analysis

Recorded here so the next chat doesn't forget, but explicitly NOT for the next chat (which is wiring, not new code).

**The idea**. Instead of asking "where are the edges of this inversion," slide along the chromosome tracking compatibility-group composition. Extend a window while composition is stable. When composition changes, record the position AND the shape of the change:
- Sharp flip (1 window) → inversion breakpoint
- Gradual (several windows) → recombinant edge
- Bubble (change then return) → gene conversion or double crossover
- Unstable composition consolidating → family LD transition

This is strictly more informative than the current breakpoint workflow: it contains breakpoints as a derived output AND captures the three other phenomena the breakpoint workflow collapses to a single position.

**Relationship to existing Module 3**. C01j already computes compatibility groups per window. What's missing: explicit seed-and-extend block protocol, transition-shape characterization, per-sample trajectory within blocks, chromosome-wide deployment.

**Five decisions required before coding** (do not make these on the user's behalf):
1. State representation per window (single encoding / C01j compatibility groups / parallel 4-encoding)
2. Regime-block definition (persistence drop threshold? compositional entropy?)
3. Transition-shape characterization (block width? per-sample trajectory? both?)
4. Output schema (per-block + per-window + per-sample rows, or single denormalized?)
5. Computational budget (500–1000 lines R, or 300–500 lines C + R wrapper)

**This comes AFTER the wiring session, not during it.**

---

## Manuscript chunks — detection floors, deer mice comparison

Documented in `docs/MANUSCRIPT_CHUNKS.md`. Contains:
- Methods paragraph: detection floor ≈ 50 kb, MAF ≈ 15% (placeholder, needs empirical confirmation)
- Results subsection drafts: sensitivity sweep figure + breakpoint repeat-context enrichment figure
- Discussion paragraph: comparison to deer mice (Peichel et al.) — dosage method reaches smaller inversions than read-pair-only approaches
- Limitations paragraph: what the method cannot detect and why

Also contains a shortlist of **unresolved numbers** the user needs to fill in from LANTA runs (total N inversions, size/MAF distributions, SD enrichment p-value, etc.).

**Two future script ideas** are documented in that same file:
- Sensitivity sweep script (size × MAF subsampling on LG28)
- Breakpoint repeat-context enrichment script

Both are manuscript-augmentation scripts, not pipeline-integration. They fit naturally as Module 2 post-processing OR as a new Module 4 (manuscript statistics). **NOT for the next chat** — wiring comes first. Recorded here so the user and next-chat-Claude don't forget them.

---

## Open question: C01f STEP02 schema vs script 07 input format

Script 07 (breakpoint evidence pileup) expects a long-format per-read-event TSV with columns: `sample, bp_side, read_type, read_pos, read_end, mate_pos, strand, clip_len, side, orient, depth`.

The user's existing C01f STEP02 output is a wide-format per-sample summary: `sample, group, support, pe, sr, n_discordant_FF, n_discordant_RR, n_split_reads, n_soft_clip, n_mate_at_partner, supporting_read_names`.

**These don't align.** Wide summary has COUNTS, not POSITIONS — but 07 needs per-read positions to draw tick marks and arcs.

**Two options for next chat:**

- **Option A (simple)**: script 07 runs in Mode B only (fresh pysam extraction from BAMs). Works always. ~2 minutes per candidate on 30 carrier BAMs — fine for hero figures, not for bulk.
- **Option B (structural fix)**: modify C01f STEP02 to ALSO write a long-format per-read TSV alongside its per-sample summary. Small change to the existing pysam loop. Then 07 Mode A works directly.

**Recommendation**: start with Option A. If panel D works visually on LG28 via Mode B, leave it. Only do Option B if bulk rendering (all 28 chromosomes × N candidates) becomes a need.

---

## What validation looks like for LG28

If the pipeline works, these are the checks that should pass on LG28:

- `final_left_bp` within 30 kb of 15.115 Mb
- `final_right_bp` within 30 kb of 18.005 Mb
- `left_ci_width_kb` < 100, `right_ci_width_kb` < 100
- `n_methods_agreeing_left` ≥ 3 and `n_methods_agreeing_right` ≥ 3
- DELLY 14.87 Mb call flagged with `within_agreement_band = FALSE`

Defined in full in METHODOLOGY §5.

---

## What is NOT done (explicit non-goals for next session)

Do not spend time on:
- Runtime testing scripts in the Claude sandbox (no R interpreter; just run on LANTA)
- Rebuilding the architecture classifier (`_archive/architecture_classifier/`). Already explicitly archived; user told me mid-session to stop building more classifiers.
- Adding Fst-based breakpoint estimation. User wants dosage-based.
- Using STEP37/STEP39 SV data as primary signal. User: "doesn't work for LG28 (DELLY called 14.87–14.94 Mb when real breakpoint is 15.115 Mb — 200 kb off)." Weight stays at 0.5.
- Inventing a registry schema. Ask the user; he has one.
- Coding the rolling regime analysis. Five design questions must be decided first.
- Modifying the assembly-alignment track in script 07. User doesn't have assembly data at 9× per-sample; the track is placeholder. Either drop it or leave the schematic — don't try to fill it with real data that doesn't exist.

---

## Files produced this session (complete list)

```
breakpoint_pipeline/
  README.md                                 plain-English entry point
  PIPELINE_DIAGRAM.md                       unified tall-bar (Modules 2 + 3 visualization)
  run_pipeline.sh                           wrapper: 01-04 + optional 05
  01_dosage_signal.R                        Module 2
  02_ancestral_fragments.R                  Module 2
  03_consensus_merge.R                      Module 2
  04_diagnostic_figure.R                    Module 2
  05_four_encoding_diagnostic.R             Module 2 optional
  06_regime_stream_graph.R                  Module 3 visualization
  07_breakpoint_evidence_pileup.py          Module 2 visualization (panel D)
  test_07_synthetic.py                      smoke test for 07

  docs/
    HANDOFF.md (this file)                  where you stopped
    METHODOLOGY.md                          math + provenance (§3.9 covers 05)
    SESSION_AUDIT.md                        decisions D1-D19 + mistakes corrected
    MANUSCRIPT_CHUNKS.md                    draft paragraphs + future script ideas
    example_panel_D_output.png              reference render of script 07 output

  _archive/
    README.md                               what each archive subfolder is
    heatmap_ordering_upgrade/               Module 1 scaffolding
      within_group_ordering.R
      STEP36_ordering_comparison_heatmaps.R
      STEP41_separate_first_heatmaps.R
      q08b_shelf_heatmap_architecture.R
      STEP_Q08b_shelf_heatmap_architecture.sh
    architecture_classifier/                abandoned mid-session
      STEP40_architecture_classifier.R
      test_classifier_selftest.py
      STEP28_patch.txt
      STEP35_patch.txt
    original_workflow_diagrams/             superseded by PIPELINE_DIAGRAM.md
      WORKFLOW_DIAGRAM.md
      WORKFLOW_DIAGRAM_BP05.md
```

**Load order for the next chat**: this file → METHODOLOGY.md → SESSION_AUDIT.md → PIPELINE_DIAGRAM.md if visual layout needed.
