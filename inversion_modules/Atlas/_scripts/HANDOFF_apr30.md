# HANDOFF — next session entry point

**Date written:** April 30, 2026
**Status:** Atlas architecture fully specified. Zero code tested against
real data. Cluster work pending.

If you (next-Claude) only read one file: read this, then read
`SESSION_AUDIT_apr30.md` if you want detail.

---

## What Quentin should do FIRST (in order)

1. **LG28 dry-run on cluster** (~1 hour). Follow
   `RUNBOOK_produce_phase2_jsons.md` Section 1 step-by-step. Produces
   `LG28_ghsl.json` and `LG28_theta.json`.

2. **Browser smoke test**. Drag-drop the two JSONs into `Atlas/Inversion_atlas.html`,
   open browser console, verify:
   - No JS errors
   - Page 3 panels render (zone bar, color stripe, dosage L1/L2 bars, GHSL Z, GHSL lines, GHSL PCA, GHSL MDS)
   - Page 12 renders (currently text-only header, will say "θπ tracks ready" if JSON loaded)

3. **Send back evidence.** Screenshot the atlas, paste any console
   errors, and a structural dump:
   ```bash
   jq 'paths(scalars) | . as $p | {path: $p, type: (getpath($p) | type)}' \
       LG28_ghsl.json | head -100
   ```

4. **Tell me what broke.** Almost certainly something will. The first-run
   bugs are usually field-name mismatches between the R exporter and the
   atlas renderer. Fix is fast once we see the actual error.

---

## Decisions NOT to redo

These are locked. Do not re-litigate without strong reason:

- **4 discovery streams**: dosage + GHSL + θπ + ancestry-Q.
  NOT 5 (rare-MAF demoted to diversity atlas). NOT 3.
- **Q page + Ancestry page split** (ADR-13). Asymmetric naming is
  intentional. NOT a tab toggle.
- **CUSUM as per-sample resolution layer**. Three CUSUMs (R02, T05,
  SC01 cell cliff-walker), shared kernel. CUSUM compute is cluster-side;
  JS just renders. **No 4th CUSUM for Q-regime.**
- **NGSadmix stays genome-wide K=8** on NAToRA-pruned 81 unrelated.
  Don't re-run per-chromosome to "fix" LG28.
- **GHSL grid 5kb fixed**, θπ at win10000.step2000 (~16,500 windows
  on LG28).
- **226 pure C. gariepinus only.** NOT F1 hybrid. NOT C. macrocephalus.

If a future session proposes any of these reversals: refer them to
`ARCHITECTURE_DECISIONS.md`, `ADR-13_q_and_ancestry_pages.md`,
`specs/CUSUM_SPEC.md`, or `specs/QR03_SPEC_v2.md` as appropriate.

---

## Specs ready to implement (gated on cluster output)

| Spec | What to build | Estimate |
|---|---|---|
| `JSON_CONTRACT.md` | Replace 4 stub renderers on page 3 | ~2 days |
| `JSON_CONTRACT.md` | STEP_TR_B v5 retrofit (sim_mat + MDS for θπ) | ~2 days |
| `JSON_CONTRACT.md` | Page 12 panel renderers | ~2 days |
| `specs/CUSUM_SPEC.md` | `cusum_core.R` + R02 + T05 + DC06 | ~3 days |
| `ADR-13` | Ancestry page (11-track shelf + concordance grid) | ~3-5 days |
| `specs/QR03_SPEC_v2.md` §5 | LG28 validation script (one afternoon) | gates QR02 |
| `specs/QR03_SPEC_v2.md` | STEP_QR02 + QR03 + Q page | ~7-8 days |

---

## What NOT to start

- Manuscript methods sections — wait for real results
- Synthetic test JSON — Quentin was not interested
- Diversity atlas — separate product, not blocking inversion work
- Re-running NGSadmix per-chromosome — see `specs/QR03_SPEC_v2.md` §6
- Adding rare-MAF back as 5th stream — see `specs/QR03_SPEC_v2.md` §8

---

## Critical biological reminders

- 226 pure C. gariepinus hatchery cohort
- 9× coverage, hatchery F1 violates HWE
- K = broodlines, not geography
- Phase blocks WhatsHap on short reads: 44–819 bp
- BAM manifest col3 = filtered BAM (`pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv`)
- Rscript: `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`, env `assembly`
- Working dir: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/`

---

## Manuscript angle locked

**First published inversion pipeline to properly handle hatchery /
breeding-program structure** (per `specs/QR03_SPEC_v2.md` §7). When
writing methods, lead with: "we developed this pipeline specifically
for aquaculture genomics in cohorts with known pedigree structure."
The Q-regime + dual-clustering + 11-track shelf-QC are the
methodological contributions reviewers care about.

---

## When in doubt

Read `SESSION_AUDIT_apr30.md` for the comprehensive trace of what was
decided this session, including the failure-mode patterns to avoid.
