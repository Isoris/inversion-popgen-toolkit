# AUDIT_LOG_chat12_2026-04-18.md

Session log for chat 12 (phase-4 coherence audit + DAG rewrite).

**Date:** 2026-04-18
**Preceded by:** chat 11 (registry library buildout) + chat 11.5
(decompose redesign via C01j/C01l/C01m).
**Follows into:** chat 13 (registry wiring), chat 14 (first HPC run).

---

## Primary deliverables

1. `AUDIT_PHASE4_COHERENCE_2026-04-18.md` — the coherence audit.
   Per-script section, flow-of-evidence trace, pattern-enum coverage,
   threshold drift inventory, findings AO–BE, chat-13 wiring
   recommendations.

2. **Code:**
   - `lib_recomb_combination.R::derive_R_from_regime` rewritten with
     DAG formulation. 50/50 tests passing.
   - `gene_conversion_detector.R` rewritten as per-SNP run-length
     detector. 33/33 tests passing.
   - `STEP_C01i_decompose.R` patched (seed loader wired, unsupervised
     fallback removed).
   - `STEP_C01i_b_multi_recomb.R` full rewrite (R+G gate via
     lib_recomb_combination, 100 kb threshold removed, inline cheat24
     deleted, Signal 1 demoted to Tier-4 diagnostic).
   - `STEP_C01f_hypothesis_tests.R::comp_from_registry` collapsed
     to one `get_groups_for_candidate` call.
   - `plot_sample_regime_dag.R` (new).
   - `tests/test_dag_derive_R.R` (new, 50/50 PASS).
   - `tests/test_gc_detector.R` (new, 33/33 PASS).

3. **Schemas:**
   - `regime_sample_dag.schema.json` (new). 5 keys_extracted.
   - `gene_conversion_tracts.schema.json` updated with
     `run_length_flagged_snps`, `tolerance_snps`, `confidence`,
     `span_bp`, `snp_qc` accounting block, `params` block.
     Legacy `tract_width_bp` preserved as alias of `span_bp`.

---

## Findings closed this session

- **V** (cheat24 inline fallback divergence): inline fallback deleted
  from multi_recomb. If real cheat24 unavailable, posterior=NA.
- **AO** (C01d doc drift: 10 vs 12 dims): mis-doc identified, trivial
  fix deferred as text-only; code is correct.
- **AU** (decompose internal fallback path): closed — chat-12 patch
  now also catches `km$seeded = FALSE` and emits `no_seeding`.
- **AX** (nested_composition vs C01l scope overlap): resolved — no
  conflict, complementary, wire both.
- **Z'** (GC detector wrong scale): closed by rewrite.

---

## Findings opened this session (AO+)

| # | Severity | Owner | Summary |
|---|----------|-------|---------|
| AO | doc | chat 15 | C01d header says "10 dimensions" in one place |
| AP | feature | chat 14 | Pattern enum missing `partial_overlap_neighbor` |
| AQ | bug risk | chat 13 | C01e Panel H path after C01j directory move |
| AR | feature | chat 14 | C01g+C01l integration needs boundary_v2 schema |
| AS | calibration | chat 14 | C01j structure_score cutoffs on Quentin's cohort |
| AT | thresh | chat 13 | C01l flank_bp should scale with candidate span |
| AU | bug | chat 12 | CLOSED |
| AV | thresh | chat 13 | Silhouette cutoff 0.25 → 0.40 for 9x |
| AW | calibration | chat 14 | min_dco_bp = 200 kb guess |
| AX | scope | chat 12 | CLOSED |
| AY | doc | chat 13 | seed_loader drop-conflict is stricter than doc |
| AZ | thresh | low | STEP03 seed FDR |
| BA | thresh | low | STEP03 seed support scale by span |
| BB | bug risk | chat 13 | ghsl_confirmation annot_dt uniqueness assertion |
| BC | metric | chat 13 | compute key_spec wired-fraction |
| BD | bug risk | chat 13 | C01k path verification |
| BE | feature | chat 14 | Add 3 pattern-enum buckets |
| BF | terminology/doc | manuscript | GC shorthand → "arrangement-discordant IBS tracts consistent with GC". Schema desc + detector header updated chat 12; fig legends & ms text to follow |

---

## Test results

**test_dag_derive_R.R** — 50/50 PASS.

Key patterns verified from the handoff spec:

- `AAABBA`: dominant=A, n_distinct=2, n_edges=2,
  deviation_fraction=1/3, longest_dev_windows=2, terminal_matches_start
  TRUE, R_fired TRUE under default gates.
- `AAABABBA`: deviation_fraction=3/8, longest_dev_windows=2, n_edges=4.
  Bp gate protects against firing R when window size small (250 kb
  threshold rejects; 150 kb accepts).
- `AAAAAA`: all zeros, R_fired FALSE.
- `AAABBB`: dominant=A (first-occurrence tie-break),
  terminal_matches_start FALSE (weak recombinant signal).

Cohort aggregates, interval filtering, empty/NULL inputs, and both
back-compat paths (`as_regime_dt` unwrap + `combine_cohort_recomb`
list acceptance) all covered.

**test_gc_detector.R** — 33/33 PASS.

Key behaviours verified:

- 3-SNP tract on HOM_REF sample → detected, MEDIUM confidence,
  direction INV_in_REF_context.
- 1 interior HET-looking SNP tolerated inside a 4-SNP run.
- 2 consecutive skips break the tract.
- max_flagged_snps gate rejects 8-flag runs (recombinant territory).
- max_span_bp gate rejects 6 kb spans at 5 kb threshold.
- Paralog artefact SNP (all-het column) dropped by QC.
- Non-diagnostic SNPs excluded.
- HET sample produces 0 tracts (baseline class skip).
- HOM_INV baseline with REF-tract → direction REF_in_INV_context.

**plot_sample_regime_dag.R** — smoke-tested on synthetic 5-sample
cohort; 4.8 KB PDF produced successfully.

---

## Parse-check status

All modified R files pass `parse()`:

- STEP_C01f_hypothesis_tests.R: 80 exprs OK
- STEP_C01i_b_multi_recomb.R: 34 exprs OK
- STEP_C01i_decompose.R: 33 exprs OK
- lib_recomb_combination.R: 7 exprs OK
- gene_conversion_detector.R: 7 exprs OK
- plot_sample_regime_dag.R: sources cleanly

---

## Posture for chat 13

Chat 12 produced a plan that chat 13 can execute mechanically without
second-guessing the architecture. The audit is the architectural
blessing; chat 13 does registry-block routing for C01j/C01l/C01m,
seal updates for the new subgroup names, and characterize Q2
enrichment.

Chat 14 is the first HPC run on LANTA. Calibrate structure_score
cutoffs (AS) and min_dco_bp (AW) from observed distributions.
