# `phase_8_evidence_biology/q7_existence_audit/` — existence evidence audit

## What this sub-directory answers

**Q7 of the 7-question framework** — *does this inversion really exist?*
More specifically: for a candidate that already passed the 4-layer
existence test (Layer A dosage / Layer B SV concordance / Layer C GHSL /
Layer D Fisher OR), what's the *detailed evidence audit* behind the
verdict? This is where the Q7B breakpoint-audit keys are written.

Distinct from:

- **Phase_7 validation** (writes `q6_group_validation`) — that's about
  whether the *carrier groups* are trustworthy.
- **`../q6_group_robustness/`** (leave-one-family-out Fst) — that's a
  *robustness diagnostic* on the groups, feeding `q6_*` keys.

This folder is purely about whether the **inversion call itself** is
well-supported at read and population level.

## Validation gate

Both remaining scripts run at **NONE** — they count read evidence
per carrier and don't need group statistics.

## Contents

| Script | Language | Role |
|---|---|---|
| `STEP_02_bnd_sided_support.py` | Python | **Single-sided BND rescue scoring** (cohort). For every Phase 2/3 candidate, checks whether there's a CT=3to3 BND near the left boundary and/or CT=5to5 BND near the right boundary, regardless of whether pairing (STEP06) was possible. Half-sided support is strictly stronger than no support at low coverage. Runs as STEP 2 in `run_evidence_biology.sh`. |
| `STEP_03_per_read_evidence.py` | Python | **Per-read BAM evidence ledger + VCF sidecar** (per-candidate, added pass 21). Extracts FF/RR discordants, split reads, soft-clips directly from BAM via pysam — persists read names + tags (NM/AS/SA/MC/XA) + 5'/3' clip orientation. Writes row-level `evidence_reads_bam.tsv.gz`, aggregate-only `evidence_reads_vcf.tsv.gz`, and a `per_read_evidence_summary` registry block with `q7b_bam_n_{split,discord,soft}_{left,right}_{INV,HET,REF}` counts. Coords read from registry. Closes the Q7B flat-key extraction gap. |
| `breakpoint_evidence_audit.py` | Python | **Population-level SV-caller audit.** Reads DELLY/Manta VCF FORMAT fields (DV, RV, PR, SR). Counts PCA carriers with strong / weak / no SV support, plus non-carriers with unexpected support. Complementary to STEP_03: STEP_03 reads BAM directly (row-level evidence with read names); this reads the SV caller's aggregate summaries. 25 Q7B keys still aspirational here — STEP_03 fills part of the Q7B gap directly (bam-side counts). |

## Evidence source split

The two row-level/aggregate questions are handled by different scripts:

| Source | Script | What it extracts | Has read names? |
|---|---|---|---|
| **BAM pileup (pysam)** | `STEP_03_per_read_evidence` | Row-level FF/RR/split/clip reads with coords, strand, CIGAR, tags, 5'/3' clip end | Yes — the ledger |
| **VCF FORMAT fields** | `STEP_03_per_read_evidence` (sidecar) + `breakpoint_evidence_audit` | Aggregate per-sample counts (DV, RV, PR, SR) + genotypes | No — DELLY/Manta don't export names |

STEP_03 emits both in one invocation but stores them in separate TSVs with different schemas, joinable by `(candidate_id, sample_id, bp_side)`.

## What was here, what's gone

- **`cheat6_ancestry_jackknife_v934.R`** moved out in pass 17
  (2026-04-24) → `../q6_group_robustness/`. Reason: cheat6 asks "is the
  Fst signal robust to dropping one ancestry family at a time?" That's
  a *group-trust robustness diagnostic*, not an existence test. The
  production jackknife (T9) is implemented inline inside C01f at
  `../../phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R`;
  the standalone at `../q6_group_robustness/` is kept for per-candidate
  manuscript plots and audit of the family-dropout pattern.

## Writes

| Key family | Writer | Status |
|---|---|---|
| `q7b_bam_n_split_{left,right}_{INV,HET,REF}` (6 keys) | STEP_03_per_read_evidence | **wired (pass 21)** |
| `q7b_bam_n_discord_{left,right}_{INV,HET,REF}` (6 keys) | STEP_03_per_read_evidence | **wired (pass 21)** |
| `q7b_bam_n_soft_{left,right}_INV` (2 keys) | STEP_03_per_read_evidence | **wired (pass 21)** |
| `q7b_bam_ledger_path`, `q7b_vcf_sidecar_path`, `q7b_window_bp`, `q7b_min_mapq`, `q7b_per_read_source` (5 keys) | STEP_03_per_read_evidence | **wired (pass 21)** |
| `q7b_delly_*_carriers`, `q7b_manta_*_carriers`, `q7b_*_carrier_loss` (dropout) | breakpoint_evidence_audit | aspirational (7 keys) |
| `q7b_delly_bnd_{3to3,5to5}`, `q7b_manta_bnd_inv{3,5}`, `q7b_bnd_rescued`, `q7b_bnd_rescue_concordance` (fragmentation) | breakpoint_evidence_audit + bnd_sided_support | aspirational (8 keys) |
| `q7b_delly_site_*`, `q7b_manta_site_*` (filter-pass) | breakpoint_evidence_audit | aspirational (7 keys) |
| `q7b_pop_prior_*` (3 keys) | breakpoint_evidence_audit | aspirational (3 keys) |
| `q7b_bnd_{left,right}_support` | bnd_sided_support | wired |

Q7B status after pass 21: 19 keys wired from STEP_03 (BAM row-level
evidence) + 2 wired from bnd_sided_support = 21 wired. 25 VCF-aggregate
keys still aspirational inside breakpoint_evidence_audit — those cover
carrier-dropout / filter-pass / fragmentation patterns that need the
DELLY/Manta VCF-side logic, not the raw BAM side.

## Running

Both scripts are cohort-level and run via `../run_evidence_biology.sh`:

```bash
DELLY_INV_VCF=... DELLY_BND_VCF=... MANTA_INV_VCF=... MANTA_RAW_VCF=... \
    ../run_evidence_biology.sh
```

## Upstream dependencies

- DELLY and Manta VCFs (INV, BND, raw)
- Phase 4 candidate catalog (coordinates)
- Phase 7 C01i decomposition output (for PCA carrier lists)

## Downstream consumers

- `../../phase_9_classification/characterize_candidate.R::characterize_q7()`
  — reads `q7_*` + `q7b_*` keys for the existence verdict.
- `STEP_04_assign_structural_class.py` (phase_9) — reads `q7_verdict` +
  `q7_tier`.
- Axis 1–4 in the 14-axis final classification (four existence layers).

## Known issues / TODO

- **Q7B flat-key extraction — partial (pass 21).** STEP_03 now wires
  19 of the 25 Q7B keys directly from BAM read evidence. The 6 remaining
  aspirational keys cover DELLY/Manta VCF-side behaviors (carrier
  dropout, filter-pass, fragmentation) that belong inside
  `breakpoint_evidence_audit.py` and still need `add_evidence()` wiring.
- **bnd_sided_support Q7B extension not fully wired.** `q7b_bnd_rescued`
  and `q7b_bnd_rescue_concordance` depend on coordinating with phase_3
  STEP06's rescue pipeline; the coordination is aspirational.
- **Integration of STEP_03 output into `07_breakpoint_evidence_pileup.py`
  (phase_6 figure generator).** The two scripts currently do independent
  pysam extractions of the same data. A future pass should refactor 07
  to read STEP_03's `evidence_reads_bam.tsv.gz` (mode A) and drop its
  own extractor. Non-urgent — both work correctly in isolation.

## See also

- `../q6_group_robustness/` — jackknife robustness diagnostic (where
  cheat6 lives now)
- `../orchestrator/LAUNCH_group_cheats.sh` — per-candidate dispatcher
- `../run_evidence_biology.sh` — cohort-level driver (this folder's
  scripts live here)
- `../../phase_3_refine/STEP_D03_statistical_tests_and_seeds.py` —
  writes Layer D Fisher OR (feeds `q7_layer_d_*`)
- `../../phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R` —
  runs the production T9 jackknife internally (feeds `q6_family_linkage`)
- `../../phase_9_classification/SPEC_VS_REALITY.md` §"breakpoint_evidence_audit.py
  (Q7B audit)" — the 25 aspirational keys list
