# `phase_8_evidence_biology/q4_mechanism/` — formation mechanism (NAHR / NHEJ / MMBIR)

## What this sub-directory answers

**Q4 of the 7-question framework** — *how did this inversion form?*
For each confirmed candidate, classify the breakpoint junction into one
of the known formation-mechanism classes: NAHR (non-allelic homologous
recombination), NHEJ (non-homologous end joining), MMEJ/alt-NHEJ
(microhomology-mediated), or MMBIR/FoSTeS (replication-based template
switching).

## Script name lookup

The filenames use a `cheatN` numbering inherited from the
`flashlight_cheats` development phase (preserved for provenance — the
numbers link to chat history and `_archive_superseded/` entries like
`cheat17_fossil_detection/`). The table below maps them to descriptive
short names so the files are self-documenting without renaming:

| Short name | File | Role |
|---|---|---|
| **sd_substrate_minimap2** | `sd_substrate/sd_substrate_minimap2.R` | **Angle A — trusted NAHR classifier.** minimap2 self-alignment on ±50 kb flanks (`-ax asm10 -X --eqx --secondary=yes -N 50 -p 0.1`). Outputs sorted/indexed BAM + strand-colored BED for IGV inspection. Writes `sd_substrate_minimap2` block. |
| **sd_substrate_biser2** | `sd_substrate/sd_substrate_biser2.R` | **Angle B — BISER2 catalog lookup.** Interval overlap on pre-computed SD catalog. Writes `sd_substrate_biser2` block. Fast (~0.1 s/cand). |
| **sd_substrate_concordance** | `sd_substrate/sd_substrate_concordance.R` | **Joint verdict.** Reads A + B blocks, emits concordance + confidence per side and overall. Writes `sd_substrate_concordance` block (mirrors into `q4_sd_concordance` / `q4_mechanism_confidence` for backward compat). |
| **mmbir_substrate** | `cheat28_tandem_repeat_context.R` | TRF + GMATA tandem-repeat context for MMBIR/FoSTeS |
| **junction_ref** | `cheat29_junction_forensics.R` | Reference-genome junction classifier (BLUNT / MH / INSERTION / DUP / DEL) |
| **junction_asm** | `STEP_01_assembled_junction.py` | Assembled-junction classifier from DELLY CONSENSUS / Manta PRECISE. Listed as STEP 1 in the cohort-level `run_evidence_biology.sh` chain. |

When you see `cheat28`, `cheat29`, or `cheat29b` in a registry block or
log line, use this table to map to its role. The `cheatN` numbering is
historical provenance — the filenames stay because they link to chat
history and `_archive_superseded/` entries.

`STEP_01_assembled_junction.py` uses the STEP_NN naming because it's
invoked directly by the cohort-level `run_evidence_biology.sh` driver
as step 1 of 4. The three sd_substrate scripts and the three `cheatN`
scripts are per-candidate, dispatched by `LAUNCH_group_cheats.sh`.

## Validation gate

**NONE** — every script here runs regardless of `q6_group_validation`.
Mechanism is a property of the junction sequence, not of the carrier
groups, so it can be called for every candidate that has a
reference-genome coordinate, including SUSPECT ones.

This is why the LAUNCH orchestrator (`../orchestrator/LAUNCH_group_cheats.sh`)
dispatches these four scripts unconditionally.

## Contents (detailed)

| File | Language | Detailed role | Primary inputs |
|---|---|---|---|
| `sd_substrate/sd_substrate_minimap2.R` *(Angle A)* | R | **Trusted NAHR classifier.** minimap2 self-alignment on ±50 kb flanks. Flags: `-ax asm10 -X --eqx --secondary=yes -N 50 -p 0.1`. Emits sorted/indexed BAM (`self_align.bam[.bai]`), exact command log (`self_align.cmd.txt`), parsed hits TSV, strand-colored BED9 (`inverted_hits.bed` red=inv, green=direct) and breakpoint region BED (blue). ~30 s/candidate. | Reference FASTA, candidate chrom + coords |
| `sd_substrate/sd_substrate_biser2.R` *(Angle B)* | R | **BISER2 catalog cross-check.** Interval overlap on pre-computed SD catalog. Finds SDs straddling both breakpoints. Emits `flanking_sds.tsv` + `flanking_sds.bed` (red=inv, green=direct). Pure TSV — no external binary. <0.1 s/candidate. | BISER2 catalog TSV, candidate coords |
| `sd_substrate/sd_substrate_concordance.R` *(joint)* | R | **Joint verdict.** Reads Angle A + Angle B blocks, produces per-side concordance (agree_nahr/agree_nhej/disagree_*/*_only/complex_flagged/no_data) + confidence (high/medium/low). Candidate-level confidence is worst of {left, right}. Mirrors to `q4_sd_concordance` + `q4_mechanism_confidence` for STEP_04 consumer. | Angles A and B blocks in `<cid>/structured/` |
| `cheat28_tandem_repeat_context.R` *(mmbir_substrate)* | R | **TR/VNTR context** — dual role: (1) flag MMBIR/FoSTeS candidates via tandem-repeat expansions at junctions; (2) mark windows with high STR density as theta noise for Cheat 12. | TRF BED, GMATA SSR BED, boundary catalog |
| `cheat29_junction_forensics.R` *(junction_ref)* | R | **Reference-junction forensics** — extract junction sequence from the reference at the estimated boundary, classify as BLUNT / MH / INSERTION / DUPLICATION / DELETION per Porubsky 2022 + Sultana 2017. | Reference FASTA, breakpoint pair table |
| `STEP_01_assembled_junction.py` *(junction_asm)* | Python | **Assembled-junction forensics** — read DELLY CONSENSUS + Manta PRECISE ALT to get the *actual* assembled junction sequence from split reads. Strictly better than `junction_ref` when a PRECISE record is available. Runs as STEP 1 in `run_evidence_biology.sh`. | DELLY INV/BND VCFs, Manta INV VCF |

## How the scripts fit together

```
SD substrate:       sd_substrate_minimap2 (Angle A, trusted) ─┐
                                                              ├── concordance ──▶ q4_sd_* keys
                    sd_substrate_biser2   (Angle B, catalog) ─┘

Junction:           junction_ref (ref) ── agree/disagree ─────┐
                    junction_asm (asm) ── primary if PRECISE ─┤
                                                              ▼
TR context:         mmbir_substrate (TRF/GMATA) ── 3rd class ─┤
                                                              ▼
                                    mechanism.schema.json block
                                    q4_* flat keys
```

Precedence rules:
1. **sd_substrate_minimap2 is trusted over sd_substrate_biser2.** When
   they disagree, minimap2 wins the `q4_sd_primary_call` but the
   disagreement is surfaced in the concordance block and the overall
   confidence drops accordingly. BISER2 is never allowed to override
   a positive minimap2 NAHR call to "no NAHR"; it can only reinforce
   or flag it.
2. **junction_asm** (assembled junction) overrides junction_ref whenever
   a DELLY CONSENSUS or Manta PRECISE record is available —
   junction_ref describes reference context, which can be off by
   kilobases at NAHR-like events.
3. **mmbir_substrate** (TR context) is a *third mechanism class*
   (MMBIR/FoSTeS), orthogonal to the NAHR-vs-NHEJ axis — it feeds a
   separate flag rather than voting in the NAHR/NHEJ call.
4. `mechanism_class` (the top-level Q4 verdict) is assigned by junction
   analysis (junction_ref + junction_asm) + sd_substrate evidence. No
   single script writes it alone.

## Writes

Each script writes its own Tier-2 block:

- **sd_substrate_minimap2** (Angle A) writes
  `sd_substrate_minimap2.schema.json`. Flat keys: `q4a_mm2_*`.
- **sd_substrate_biser2** (Angle B) writes
  `sd_substrate_biser2.schema.json`. Flat keys: `q4b_biser2_*`.
- **sd_substrate_concordance** (joint) writes
  `sd_substrate_concordance.schema.json`. Flat keys:
  `q4_sd_concordance_{left,right,overall}`,
  `q4_sd_confidence_{left,right,overall}`,
  `q4_sd_primary_call`, `q4_sd_primary_call_source`. Also mirrors
  `q4_sd_concordance` and `q4_mechanism_confidence` into
  `mechanism.schema.json` for backward compat with
  `STEP_04_assign_structural_class.py`.
- **junction_ref wired** — junction class, microhomology length,
  inserted sequence. Block: `mechanism` (junction portion).
- **junction_asm wired** — `q4b_asm_*` sub-block of
  `mechanism_assembled.schema.json`.
- **mmbir_substrate aspirational** — 10 keys (`q4_n_tr_at_bp_{left,right}`,
  `q4_dominant_tr_class_*`, `q4_mmbir_signal_*`). Computations exist;
  flat-key extraction not yet wired.

## Running

Dispatched per-candidate by `../orchestrator/LAUNCH_group_cheats.sh`
with env vars:

```bash
# Both angles + concordance
REF_FASTA=/path/to/ref.fasta \
BISER2_TSV=/path/to/biser2.tsv \
sbatch ../orchestrator/LAUNCH_group_cheats.sh --chroms LG01,LG28

# Fast BISER2-only sweep (skip minimap2)
NO_MINIMAP2=1 \
BISER2_TSV=/path/to/biser2.tsv \
sbatch ../orchestrator/LAUNCH_group_cheats.sh --chroms all

# Targeted minimap2 with stricter preset
SD_MM2_PRESET=asm5 \
REF_FASTA=/path/to/ref.fasta \
BISER2_TSV=/path/to/biser2.tsv \
sbatch ../orchestrator/LAUNCH_group_cheats.sh --chroms LG28
```

Standalone per script:

```bash
# Angle A — minimap2
Rscript sd_substrate/sd_substrate_minimap2.R \
    --candidate LG28_cand_1 --chrom C_gar_LG28 \
    --left_bp 15115243 --right_bp 18005891 \
    --ref_fasta /path/to/ref.fasta \
    --outdir <per_candidate_output_root>

# Angle B — BISER2
Rscript sd_substrate/sd_substrate_biser2.R \
    --candidate LG28_cand_1 --chrom C_gar_LG28 \
    --left_bp 15115243 --right_bp 18005891 \
    --biser2_tsv /path/to/biser2_results.tsv \
    --outdir <per_candidate_output_root>

# Concordance (reads both blocks from structured/)
Rscript sd_substrate/sd_substrate_concordance.R \
    --candidate LG28_cand_1 \
    --outdir <per_candidate_output_root>

# junction_ref
Rscript cheat29_junction_forensics.R \
    --candidate <cid> \
    --outdir <per_candidate_raw_dir>
```

Intermediate files (BAM, BED, TSV) persist under
`<outdir>/<candidate_id>/sd_substrate/{minimap2,biser2}/`. Nothing
is written to tmp/. Re-running with the same args is idempotent —
the BAM is reused unless `--overwrite` is passed.

See `sd_substrate/README.md` for IGV loading instructions and
concordance matrix details.

## Upstream dependencies

- Phase 4 catalog: `candidate_id`, breakpoint coordinates
- Phase 6 breakpoint refinement: `candidate_breakpoints_consensus.tsv`
  (final_left_bp, final_right_bp with CIs), wired via
  `../bp_bridge/STEP_03B_bp_pipeline_bridge.py`
- Reference assembly FASTA + index
- External annotation: BISER2 SD pairs, TRF BED, GMATA SSR BED,
  RepeatMasker (for mmbir_substrate noise filter)
- SV caller VCFs (for junction_asm): DELLY INV/BND, Manta INV

## Downstream consumers

- `../../phase_9_classification/characterize_candidate.R::characterize_q4()`
  reads `q4_*` keys to answer "what is the formation mechanism?"
- `../../phase_9_classification/STEP_04_assign_structural_class.py`
  reads the `mechanism` + `mechanism_assembled` blocks to pick the
  `_NAHR_like_*` / `_NHEJ_like_*` / `_MMEJ_like_*` suffix on supported
  balanced inversions, and to decide `_supported_by_assembly` vs
  `_hypothesis` on that suffix.
- Axis 11 (`mechanism_class`: NAHR / NHEJ / MMBIR / unknown) in the
  14-axis final classification.

## Known issues / TODO

- **mmbir_substrate flat-key extraction** is the single biggest gap in
  q4 (10 aspirational keys). Wiring target: schema `keys_extracted`
  extension or direct `add_evidence()` call.
- **Reference-based vs assembly-based disagreement logging.** When
  junction_ref and junction_asm give different junction classes, the
  current code picks junction_asm silently. A disagreement flag (e.g.
  `q4_junction_ref_vs_asm_disagree`) would help manuscript-time audit.
- **sd_substrate standalone validation** — the merged implementation
  hasn't been smoke-tested end-to-end on LANTA yet. Candidate for
  first-use verification: LG28 with its 226-sample cohort's BISER2
  catalog + reference FASTA; should produce `agree_nahr` or
  `agree_nhej` at high confidence.

## See also

- `../orchestrator/LAUNCH_group_cheats.sh` — per-candidate dispatcher
- `../q5_age_and_origin/` — age + population-genetic origin (sister Q)
- `../cross_species/` — phylogenetic conservation (different axis)
- `../q6_group_robustness/` — group-trust diagnostic (different axis)
- `../q7_existence_audit/` — existence evidence (different axis)
- `../bp_bridge/` — breakpoint-pipeline registry wiring
- `../../phase_9_classification/SPEC_VS_REALITY.md` — aspirational key
  tracker, see §"cheat28 tandem-repeat context"
