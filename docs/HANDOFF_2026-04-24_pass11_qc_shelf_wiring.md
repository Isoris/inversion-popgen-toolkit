# HANDOFF — phase_qc_shelf scale-out wiring (partial pass 11)

**Context chat:** inversion-toolkit reorg — this chat covered pass 9 (phase_3 Layer rename), pass 10 (phase_6 MODULE_5E archive), and started pass 11 (phase_qc_shelf scale-out).

**Status:** 1 of 4 tasks shipped. Context budget ran out mid-pass.

---

## The architectural question you raised

> "phase_qc_shelf is just to try on 1 inversion. But once it's ready we need it to work for every inversion in the full genome. IDK where it fits."

We agreed on this answer:

**phase_qc_shelf is a SIBLING of phase_4_postprocessing, not a subfolder of it.**

The flow:

```
phase_2 (discover) ─→ phase_3 (SV evidence) ─→ phase_4a (catalog-birth: candidate_summary.tsv)
                                                          │
                                              ┌───────────┴───────────┐
                                              ▼                       ▼
                                      phase_4b/c/d/e          phase_qc_shelf
                                      (scoring, group         (per-candidate
                                       validation,             data-quality
                                       characterization)       diagnostic)
                                              │                       │
                                              └───────────┬───────────┘
                                                          ▼
                                               phase_4e/compute_candidate_status
                                               reads BOTH registries, assigns
                                               tier + q_qc_shelf_flag
                                               (messy inversions still get tiered,
                                                just flagged "messy")
```

Key design principle: the integration is **SOFT**. A messy inversion
still gets characterized; phase_qc_shelf just adds a flag column. No
hard gate. This matches Quentin's explicit guidance: "a messy
inversion is still worth characterizing, just flagged as messy."

---

## Task breakdown (4 subtasks from pass 11 plan)

| # | Task | Status | Where |
|---|---|---|---|
| 1 | **Bridge script**: phase_4a `candidate_summary.tsv` → phase_qc_shelf `shelf_coords.tsv` | ✅ **DONE + smoke-tested** | `4b_qc_triage/scripts/bridge_from_phase4.sh` |
| 2 | Rewrite `phase_qc_shelf/README.md` documenting both modes + bridge | ⏳ **PENDING** | `phase_qc_shelf/README.md` |
| 3 | Update `docs/MODULE_MAP.md` to place phase_qc_shelf as phase_4 sibling | ⏳ **PENDING** | `docs/MODULE_MAP.md` |
| 4 | 4e reader: new `q_qc_shelf_*` keys + `characterize_q_qc_shelf()` helper | ⏳ **PENDING** | `phase_4_postprocessing/4g_final_classification/compute_candidate_status.R` + `characterize_candidate.R` |

---

## Task 1 — DONE (in pass 11 tarball)

`4b_qc_triage/scripts/bridge_from_phase4.sh` (~90 lines bash+awk).

**What it does:** Reads `candidate_summary.tsv` (output of
`STEP_C01d_candidate_scoring_wired_25_v934_registry.R`, columns:
`chrom, interval_id, start_mb, end_mb, span_mb, final_score,
dim_positive, tier, pattern, shape_class`) and emits
`shelf_coords.tsv` in the format expected by
`slurm/array_28chrom.sh` (columns: `chrom, shelf_start_mb,
shelf_end_mb, bp1_mb, bp2_mb`).

**Filters:**
- `--min-tier N` — keep only tier ≤ N (default: no filter)
- `--min-span-kb N` — keep only span ≥ N kb (default: 50 kb)

**Safety:** Backs up existing `shelf_coords.tsv` to
`shelf_coords.tsv.bk_YYYYMMDD_HHMMSS` before overwriting.

**Smoke-tested in sandbox** with a 4-row synthetic input:
- Tier-1 candidates (LG28, LG12) kept
- Tier-3 50kb candidate (LG01) kept with default filters, filtered at
  `--min-span-kb 51`
- Tier-4 10kb candidate (LG05) correctly excluded
- Backup file created on rerun
- Output has correct header comments + tab-separated rows
- Per-chromosome count summary printed

**Known limitation** (documented in script header):
`run_all_28chrom.sh` dispatches one SLURM array task per chromosome,
so if `shelf_coords.tsv` has three rows for one chromosome, only the
first is QC'd. Full per-candidate dispatch is a follow-up item.

**Important nuance on bp1/bp2:** At phase_4a catalog-birth time,
breakpoints are NOT bp-resolution refined (that's phase_3's job via
STEP_D03 Layer D). The bridge uses `start_mb`/`end_mb` as both shelf
coords AND as initial bp1/bp2. If phase_3 has already produced
refined breakpoints via the evidence registry, a `--refined-bp` flag
could pull those instead — flagged as a TODO in the script header,
not implemented yet.

---

## Task 2 — PENDING: rewrite phase_qc_shelf README

**Current state:** 112 lines, documents only the single-candidate
LG28 Quick Start. No mention of `run_all_28chrom.sh`, the bridge, or
the 2-mode architecture.

**Rewrite structure needed:**

```
# phase_qc_shelf — data-quality QC for inversion candidates

[1-paragraph intro: what the module does, why it exists]

## Two modes

### Mode A: single-candidate prototype (quick diagnosis of one shelf)

    SHELF_START_MB=15 SHELF_END_MB=18 \
      BP1_MB=15.115 BP2_MB=18.005 \
      bash run_chrom.sh C_gar_LG28

Use when you have one specific candidate and want the full
Q01→Q09 diagnostic suite on it.

### Mode B: genome-wide production run (all candidates in the catalog)

    # Step 1: bridge phase_4 candidate list into shelf_coords.tsv
    bash scripts/bridge_from_phase4.sh \
      /path/to/phase_4a/candidate_summary.tsv \
      --min-span-kb 50

    # Step 2: review shelf_coords.tsv, then run the full sweep
    bash run_all_28chrom.sh --resume

The pipeline dispatches one SLURM array task per chromosome,
populates the 4 registries (interval/sample/results/evidence) via
Q10, and produces per-chromosome diagnostic PDFs.

## Layout (preserve current tree)

[preserve existing tree diagram + Q-step table]

## Integration with phase_4

Candidates that pass through phase_qc_shelf get a `q_qc_shelf_*` flag
set read by `phase_4e/compute_candidate_status.R`. The integration is
soft — a "messy" candidate still gets a tier, just flagged as messy.
See `docs/MODULE_MAP.md` for the sibling relationship with
phase_4_postprocessing.

## Engine B (Q06) and Engine F (Q07) — preserve existing section

## Scrubber — preserve

## Skip patterns — preserve

## Dependencies — preserve
```

---

## Task 3 — PENDING: MODULE_MAP update

Current line in `docs/MODULE_MAP.md` phase table:

> `phase_qc_shelf/` | MODULE_QC_ShelfDiagnosis | Systematic QC tracks to diagnose Z-plateau authenticity | §3.3 supplement

Needs updating to show:

1. Its role as a phase_4 sibling (not a supplement)
2. The bridge arrow: `phase_4a/candidate_summary.tsv` → `bridge_from_phase4.sh` → `shelf_coords.tsv` → `run_all_28chrom.sh`
3. The soft integration at phase_4e via `q_qc_shelf_*` keys

Add a new subsection after the phase_3/phase_4 sub-block tables:

```
### phase_qc_shelf/ mode matrix

| Mode | When to use | Entry point | Output |
|---|---|---|---|
| **A: single-candidate prototype** | One specific candidate, quick diagnosis | `run_chrom.sh <CHR>` with SHELF_* env | per-chrom Q04/Q08 PDFs |
| **B: genome-wide production** | Full cohort QC sweep | `bridge_from_phase4.sh` then `run_all_28chrom.sh --resume` | populated Q10 registry for all candidates |

Integration with phase_4: the Q10 registry entries written in mode B
are consumed by `phase_4e/compute_candidate_status.R` as
`q_qc_shelf_*` flag keys (added via pass 11 task 4 — pending).
```

---

## Task 4 — PENDING: 4e reader (the most invasive task)

**Goal:** Make `phase_4e/compute_candidate_status.R` read the
phase_qc_shelf Q10 registry output and set flag columns on each
candidate.

**What Q10 writes** (from `STEP_Q10_register.sh`):
1. `interval_registry` — `chrom, start_bp, end_bp, method="phase_qc_shelf", confidence`
2. `sample_registry` — three groups per candidate: `inv_<cid>_HOM1 / HET / HOM2` with `dimension="karyotype"`
3. `results_registry` — `kind=pairwise, stat=fst, who=hom1_vs_hom2, where=<cid>` + summary stats
4. `evidence_registry` — diagnostic PDFs + `gap_features.tsv` + `popstats_invgt.tsv` + `invgt_assignments.tsv` + `summary.json`

**New keys to materialize** (from `summary.json`):
- `q_qc_shelf_flag` — one of: clean / low_snp / high_uncertain / coverage_artifact / messy
- `q_qc_shelf_snp_density` — median SNPs per 10kb window in shelf
- `q_qc_shelf_uncertain_frac` — mean BEAGLE uncertainty in shelf
- `q_qc_shelf_cov_cv` — coefficient of variation of coverage across samples in shelf
- `q_qc_shelf_fst_hom1_hom2` — Hudson Fst between Hom1 and Hom2 groups in shelf

**Files to touch** (both already edited by pass 8 — read CAREFULLY
before editing to avoid regressing axis 5 wiring):

1. `phase_4_postprocessing/4g_final_classification/compute_candidate_status.R`:
   - Add new Q-block to `build_key_spec()` function
   - Add helper `read_qc_shelf_registry()` that reads method=phase_qc_shelf entries
   - Wire the flag columns into the main candidate loop (ADDITIVE — do not change axis 1-3 tier logic)

2. `phase_4_postprocessing/4g_final_classification/characterize_candidate.R`:
   - Add `characterize_q_qc_shelf()` that returns status/label/evidence
   - Include in the characterize_all() dispatcher

3. `phase_4_postprocessing/4g_final_classification/LAUNCH_characterize_classify.sh`:
   - Pass `QC_SHELF_DIR` env var through (new — similar pattern to the `V7_FINAL_DIR` env from pass 8)

**Test before shipping:** Write a minimal unit test in
`phase_4_postprocessing/tests/test_qc_shelf_reader.py` that mocks a
`summary.json` and verifies the keys end up in the output table.

**Critical**: Axis 5 wiring from pass 8 (env-gated by `V7_FINAL_DIR`)
must NOT be regressed. Verify its tests still pass after task 4 edits.

---

## Where to start next chat

**Suggested opening prompt to Claude:**

> Continue pass 11. Bridge script done + smoke-tested. I have the
> handoff doc (paste contents). Please do tasks 2 + 3 (README +
> MODULE_MAP doc pass, low risk). Skip task 4 for a later chat after
> I've re-read compute_candidate_status.R to confirm the right
> integration point.

Or if you want task 4 first (high value, moderate risk):

> Continue pass 11 task 4: add q_qc_shelf_* reader to
> compute_candidate_status.R. Read the existing axis 5 wiring block
> first (V7_FINAL_DIR, pass 8) to make sure my additions don't
> regress it. Write a minimal unit test.

---

## What's in the pass 11 partial tarball

- `4b_qc_triage/scripts/bridge_from_phase4.sh` (new, chmod +x)
- `DROP_README.md` (pass 11 partial)
- `docs/HANDOFF_2026-04-24_pass11_qc_shelf_wiring.md` (this file)
- `docs/SESSION_AUDIT_2026-04-24_chat3.md` (session audit, separate file)

No edits to existing files. Zero regression risk.

---

## Pass 9 / pass 10 / pass 11 progress

- ✅ **Pass 9**: phase_3_refine Layer-tagged rename (STEP_{A,B,D}NN_) — SHIPPED
- ✅ **Pass 10**: phase_6 MODULE_5E archive — SHIPPED
- 🟡 **Pass 11 (partial)**: phase_qc_shelf scale-out bridge — SHIPPED; tasks 2/3/4 pending

## Remaining phase 3-6 reorg work (post pass 11)

- phase_4b: 10 flat files, needs sub-grouping
- phase_4d: 14 flat files, needs sub-grouping
- phase_5: multi-chat reorg

## OPEN_TODOS referenced

- T1 (axis 5 wiring) — CLOSED in pass 8
- T4 (Gap 2 GDS) — CLOSED in pass 8
- T2 / T3 — BREEDING + Q7B wiring (needs prior-chat spec)
- T7 — RENAMING Cat 2 sed pass
- T12 — manuscript MAPQ inconsistency
- T5 / T6 — HPC-blocked (full v7 run + Engine H validation)

## Cohort identity reminder (critical — never conflate)

1. F1 hybrid (*C. gariepinus* × *C. macrocephalus*) — assembly paper ONLY
2. **226 pure *C. gariepinus* hatchery broodstock on LANTA — THIS paper, K=8 = broodline structure not hybrid**
3. Pure *C. macrocephalus* wild — future paper

User's full name: **Quentin Andres** (NOT Pasquet).
