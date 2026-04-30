# Diversity Atlas v2.3 — supplementary table coverage closed

## You were right — 6 source tables were missing

Audit of `MODULE3_Supplementary_Tables_7_.xlsx` (32 sheets) vs the v2.2 atlas
revealed that 6 supplementary tables had no visible representation:

| Sheet | Real data? | Where it lives now |
|---|---|---|
| **M3_S4b** per-chr × per-sample H | `[fill]` placeholders only | Tab 7 placeholder card with schema |
| **M3_S8c** all 681,286 ROH tracts | Real data (28 MB JSON if all loaded) | Tab 5 long-ROH atlas — ≥1 Mb subset (5,850 tracts, 0.30 MB) |
| **M3_S12** het in/out ROH | `[fill]` in xlsx — **but derivable from S1!** | Tab 5 — fully populated with 226 rows from S1 θ_in / θ_out / callable_bp |
| **M3_ST4** per-window θπ pointer | Pointer to external 30 MB TSV | Tab 7 pointer card with full schema |
| **M3_ST5** multiscale θπ summary | `[fill]` placeholders only | Tab 7 placeholder card with schema |
| **M3_REF** inversion-support pointers | `[fill]` placeholders only | Tab 7 placeholder card with schema |

## What this turn shipped (Diversity Atlas v2.3, 1.7 MB)

**S12 — Heterozygosity inside vs outside ROH** (fully populated, derived from S1)
- 226-row sortable + filterable table on Tab 5
- Cohort median ratio out/in = **4.27**, range 1.01 – 8.30 (Q05 = 2.85)
- 12 outliers flagged (≤ Q05) — candidates for ngsF-HMM mis-calls or coverage artefacts
- Filter pills: all 226 / outliers only
- Computed entirely from S1's `th_in`, `th_out`, `roh_total_bp`, and `callable_bp` fields
  — `het_in_ROH` = th_in, `het_out_ROH` = th_out, `callable_bp_in_ROH` = roh_total_bp,
  `callable_bp_out_ROH` = callable_bp − roh_total_bp, ratio = th_out / th_in

**S8c — Long-ROH tract atlas (≥ 1 Mb subset)** (real data)
- 5,850 individual ROH tracts, sortable + filterable on Tab 5
- Filters: sample name, chromosome (all 28 LGs), length threshold (≥1 / ≥2 / ≥4 / ≥8 Mb)
- Display capped at 1,000 rows for DOM responsiveness; filter to drill down
- Why ≥1 Mb subset and not all 681,286: the full table would be 28 MB of JSON. The ≥1 Mb subset is the *manuscript-relevant* slice (matches the F_ROH ≥1 Mb threshold in the pipeline).
- Spot-check: LG28 has 133 long tracts; only 1 of those is ≥ 8 Mb (the very-recent inbreeding signal)

**Supplementary table coverage map** (Tab 7 About)
- 24-row table mapping every source xlsx sheet to its location in the atlas
- Status badges: ✓ shipped (green) / ◐ partial (accent) / ⚪ pending (dim) / 📦 external (accent)
- Reviewer-facing transparency about what's surfaced where

**Placeholder/pointer cards** for S4b, ST4, ST5, REF (Tab 7 About)
- Each card shows the table's title, description, expected row count, source, and column schema
- Status tag (`⚪ pending` for [fill] tables, `📦 external` for ST4)
- Renders the schema as a typed-column reference table so the data structure is clear even without the values

## Updated manuscript paragraphs

The `[BRACES 0.20–0.30]` placeholder for the cohort-mean θ in/out ratio
is now filled with real numbers from S12:

> θ inside ROH was depressed to roughly **1/4.3** of θ outside ROH per
> individual (cohort median ratio out/in = 4.27, range 1.01–8.30; 12/226
> samples flagged as low-ratio outliers below the 5th percentile of 2.85,
> suggesting possible ngsF-HMM mis-calls or coverage artefacts in those
> individuals; Supplementary Table S12).

## Validated

JSDOM smoke test: 0 errors, 0 warnings.
- 25 tables on screen, 2,448 rendered table rows total
- All 8 tabs cycle cleanly
- All 7 new UI elements (S12 / S8c / coverage / 4 placeholder cards) rendering
- S12 outlier filter: shows exactly 12 ✓
- S8c LG28 + ≥8 Mb filter: shows exactly 1 ✓
- All cross-link hashes (`#page1` through `#page8`) navigate correctly

## What's still genuinely missing (upstream computation needed)

These will not be filled until cluster runs complete:

1. **S4b** per-chr × per-sample H — needs per-sample × per-chrom thetaStat batch (same SAFs already exist; one batch run away)
2. **ST5** multiscale θπ summary — multiscale tracks computed but cohort summary not yet emitted (R aggregation script needed)
3. **REF** inversion-support θπ tracks — pipeline integration with MODULE_5E pending; θ tracks exist but candidate-relative reformatting not yet done

If Quentin wants, I can write the build scripts for any of these. The atlas
will pick up the data automatically once the JSON is in place — the
loaders use the same `loadJSON('dt_X')` pattern as everything else,
and the placeholder cards already document the expected schema.

## Where this leaves things

The Diversity Atlas now surfaces **every supplementary table from the source xlsx** —
real data where available, placeholder schema where pending. A reviewer
opening the atlas can verify supplementary-data coverage by reading the
coverage map on Tab 7, then drill into any specific table from there.

For Quentin's manuscript: numerical sources are now traceable end-to-end
from claim → atlas tab → source supplementary table → S1 / S2 / etc raw
data. The earlier `[BRACES]` placeholder for θ in/out is now filled.
