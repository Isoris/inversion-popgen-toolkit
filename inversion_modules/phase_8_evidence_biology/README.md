# `phase_8_evidence_biology/` — group-dependent evidence (requires SUPPORTED+)

This sub-phase is **not yet populated**. Upload your existing
cheat / Q5 / Q6 scripts here when ready.

## Gate

Everything in this sub-phase is guarded by:

```r
if (q6_group_validation %in% c("SUPPORTED", "VALIDATED")) {
  # run the test
}
```

Candidates at UNCERTAIN or SUSPECT skip this sub-phase entirely. This
is the mechanism that keeps composite intervals (capped at UNCERTAIN
by 4b.4) from producing misleading Fst / age / burden estimates.

## What runs here

| Script | Purpose |
|---|---|
| `cheat6_ancestry_jackknife.R` | already ran once in 4c as Layer D test — may re-run here for final reporting |
| `cheat20_junction_classifier.R` | breakpoint microhomology / TE / inverted repeat junction classification |
| `cheat30_GDS.R` | genotype disequilibrium statistic (age proxy) |
| `Q5_age_fst.R` | age estimation via between-arrangement Fst |
| `Q6_burden.R` | per-sample mutational burden regression on arrangement state |

Other `cheatNN` scripts from `phase_4_catalog/cheats/` that read
validated groups belong here too.

## Writes

Each script writes its corresponding Tier-2 block:

- `age_evidence.json` (from Q5)
- `burden.json` (from Q6)
- `mechanism.json` (from cheat20)
- Additional block types per script

## Upload checklist

When you upload, also include:

- The current source of each cheat script you want to integrate
- The Q5 / Q6 drivers if separate from the cheats
- One example JSON output from each, if available

Then I can migrate each to use the registry library instead of
ad-hoc TSV reads.
