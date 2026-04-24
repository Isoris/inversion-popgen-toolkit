# `phase_6_breakpoint_refinement/` — bp-resolution breakpoint refinement from dosage

(Formerly the top-level `breakpoint_pipeline/` module. Moved into
phase_4 as sub-block 4c in pass 12, 2026-04-24. See
`docs/PHASE4_RENUMBER_PROPOSAL.md`.)

Find the exact edges of a chromosomal inversion from dosage data, at
SNP resolution, with confidence intervals.

---

## TL;DR — what to run

```bash
./run_pipeline.sh C_gar_LG28 config_inversion_followup.R 1
```

One chromosome, one config, one candidate ID. Done.

Your answer is in:
```
FOLLOWUP_DIR/C_gar_LG28.candidate_1/candidate_breakpoints_consensus.tsv
```

That file has `final_left_bp`, `final_right_bp`, and confidence intervals.

---

## What's in this folder

**Core pipeline** — these are what you need. Run in order (the wrapper does this for you).

| # | Script | What it does |
|---|---|---|
| 01 | `01_dosage_signal.R` | Finds informative markers (where HOMO_1 and HOMO_2 disagree), builds the core block, extends outward until correlation drops |
| 02 | `02_ancestral_fragments.R` | For each inversion carrier, scans outward and records where their personal inversion ancestry ends — the core of the method |
| 03 | `03_consensus_merge.R` | Combines 01 and 02 outputs (plus optional upstream signals) into the final breakpoint call with CI |
| 04 | `04_diagnostic_figure.R` | Multi-track PDF showing the signal, the per-carrier votes, and the final breakpoint call |

**Optional diagnostics** — run only if you need them.

| # | Script | When to run it |
|---|---|---|
| 05 | `05_four_encoding_diagnostic.R` | When you have a composite candidate (like LG28) and want to check whether claimed sub-structure holds under different allele-coding conventions |
| 06 | `06_regime_stream_graph.R` | When you have `regime_memberships.tsv.gz` from C01j and want to visualize the regime landscape across a whole chromosome — separate tool, not part of the breakpoint call |
| 07 | `07_breakpoint_evidence_pileup.py` | Per-sample stacked read evidence at the breakpoints (panel D replica). Python, reads C01f output OR extracts fresh from BAMs. Manuscript-grade figure. |

**Wrapper**

`run_pipeline.sh` — chains 01→02→03→04 for you. Adds 05 if you set `RUN_05=1`.

**Archive**

`_archive/` — other files produced during this session that aren't part of the core pipeline: the heatmap row-ordering fix, an experimental architecture classifier (probably don't use), and earlier workflow diagrams that have been superseded. See `_archive/README.md` for details on each.

---

## Minimal required inputs

Your config file must define these paths:

- `FOLLOWUP_DIR` — where per-candidate outputs go
- `DOSAGE_DIR` — contains `<chr>.dosage.tsv.gz` + `<chr>.sites.tsv.gz`
- `CANDIDATE_TABLE` — list of candidates with `chrom, start_bp, end_bp`
- `PLOTS_DIR` — where the diagnostic PDFs go

Each candidate must already have `candidate_pca_rotated.tsv` from STEP21 (gives the three-group karyotype assignments: HOMO_1, HET, HOMO_2).

Everything else is optional and picked up automatically if present.

---

## Optional upstream inputs (improves the consensus)

If these directories are set in config, the consensus merger will read them:

- `C01I_DIR` — `ancestral_fragments.tsv.gz` from Clair3-based decomposition
- `C01J_DIR` — `regime_transitions.tsv` from C01j regime engine
- `C01L_DIR` — `segment_summary.tsv.gz` from C01l local structure
- `STEP37_DIR` — SV caller breakpoints (used as bottom-layer sanity check, weight 0.5)

The pipeline works fine without any of these. The more you have, the tighter the CI.

---

## Run modes

**Core breakpoint call only** (default):
```bash
./run_pipeline.sh C_gar_LG28 config_inversion_followup.R 1
```

**With the robustness diagnostic** (when you want to check sub-structure claims):
```bash
RUN_05=1 ./run_pipeline.sh C_gar_LG28 config_inversion_followup.R 1
```

**All candidates on one chromosome** (instead of just `cid=1`):
```bash
./run_pipeline.sh C_gar_LG28 config_inversion_followup.R all
```

**Regime stream graph** (separate, not part of the pipeline):
```bash
Rscript 06_regime_stream_graph.R \
  --memberships results/c01j/regime_memberships.tsv.gz \
  --windows     results/c01j/regime_windows.tsv.gz \
  --chrom       C_gar_LG28 \
  --out         results/c01j/stream_LG28
```

---

## What the main output looks like

`candidate_breakpoints_consensus.tsv` has one row per candidate with these columns:

```
candidate_id           1
chrom                  C_gar_LG28
final_left_bp          15115243           <-- THE ANSWER
final_right_bp         18005891           <-- THE ANSWER
left_ci_low            15092301
left_ci_high           15138102
right_ci_low           17981544
right_ci_high           18030118
left_ci_width_kb       45.8
right_ci_width_kb      48.6
n_methods_agreeing_left    4
n_methods_agreeing_right   3
primary_source_left    ancestral_fragments
primary_source_right   ancestral_fragments
```

The two numbers you want are `final_left_bp` and `final_right_bp`. The CIs tell you how confident to be.

---

## How it works (one paragraph)

Each inversion carrier's ancestors eventually recombined somewhere. The undisturbed stretch of inversion ancestry each carrier inherited ends at that ancestral recombination (or at the true breakpoint, whichever is more interior). So each carrier gives you a personal breakpoint estimate on each side. The **modal position** across carriers = the true breakpoint. The **spread** = the CI. Samples with short fragments = old recombinant ancestry; samples with long fragments = recent uninterrupted inheritance.

More detail in `docs/METHODOLOGY.md` and the tall-bar diagram in `PIPELINE_DIAGRAM.md`. Draft manuscript paragraphs on detection floors, the deer-mice comparison, and limitations are in `docs/MANUSCRIPT_CHUNKS.md`.

---

## If something goes wrong

The wrapper logs to `logs/breakpoint_pipeline/`. Each step's output ends up there.

Common issues:
- **"no rotated PCA"** → You need to run STEP21 first to produce `candidate_pca_rotated.tsv`.
- **"insufficient informative markers"** → The candidate region might not be a real inversion, or your delta_min threshold is too strict. Defaults at the top of `01_dosage_signal.R`.
- **"too few carriers"** → Inversion is rare (< ~5 minor homozygotes + heterozygotes). Fragment distribution isn't meaningful below this count.
- **Wide CIs (> 100 kb each side)** → Either the inversion has lots of old recombinants (real biology), or SNP density near the breakpoint is low (tool limitation). Check the diagnostic PDF.

---

## Known limitations

Read `docs/METHODOLOGY.md` §8 for the full list. Short version:

- Rare inversions (< 30 carriers) have unreliable CIs
- Gene conversion looks like recombination
- Recent admixture can make the fragment distribution bimodal (method assumes unimodal)
- Dosage softness at low-coverage sites propagates to wider boundaries
- Assembly gaps at breakpoints will confound the scan — flagged in output for manual review

The method is honest about these. If your candidate hits one, the diagnostic figure will show it.
