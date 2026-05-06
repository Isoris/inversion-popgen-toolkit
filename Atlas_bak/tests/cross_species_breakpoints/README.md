# Cross-species breakpoints — LANTA pipeline

Two scripts that take the Cgar and Cmac haplotype FASTAs (already on disk
from the F₁ hybrid assembly), align them with wfmash, and emit a
`cs_breakpoints_v1` JSON ready to drop into the Inversion Atlas
"cross-species" page.

## Files

```
scripts/
├── run_wfmash_gar_vs_mac.sh         bash driver for wfmash 1-to-1 alignment
├── STEP_CS01_extract_breakpoints.py PAF -> breakpoints JSON, with TEfull
│                                    flanking density annotation
├── test_STEP_CS01.py                synthetic end-to-end test (run before
│                                    trusting the real output)
└── test_atlas_page16.js             JSDOM smoke test for the atlas-side
                                     cross-species page renderer
```

The test scripts run independently without LANTA-side data; they validate the
pipeline logic and the atlas-side renderer respectively. Both should pass
before the pipeline is trusted on real data.

## Inputs you need on LANTA

```
/project/lt200308-agbsci/01-catfish_assembly/02-annot/
├── fClaHyb_Gar_LG.fa     28 LGs, ~960 Mb (the QUERY)
└── fClaHyb_Mac_LG.fa     27 LGs, ~960 Mb (the TARGET)
```

Plus, optionally (for the flanking-repeat annotation):
- `Cgar_TE_density/LG{01..28}_repeat_density_TEfull.json` — the per-chrom
  TEfull JSONs you've already generated for the boundaries page
- `Cmac_TE_density/LG{01..27}_repeat_density_TEfull.json` — same shape but
  for the Mac haplotype (these may not exist yet — if not, the script
  emits breakpoints without Mac-side flanking annotation)

Plus, optionally (for the polymorphic-candidate overlap):
- `inversion_localpca.candidate_regions.tsv` — the 226-sample candidate
  intervals (cols: candidate_id, chrom, start_bp, end_bp; same format
  the atlas already loads)

## Step 1 — wfmash alignment

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04
conda activate assembly   # (or whichever env has wfmash)

# Test with a single chrom pair first (5-10 min):
mkdir -p test_wfmash_LG28
samtools faidx /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Gar_LG.fa LG28 \
  > test_wfmash_LG28/gar_LG28.fa
samtools faidx /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Mac_LG.fa LG28 \
  > test_wfmash_LG28/mac_LG28.fa
samtools faidx test_wfmash_LG28/gar_LG28.fa
samtools faidx test_wfmash_LG28/mac_LG28.fa
wfmash test_wfmash_LG28/mac_LG28.fa test_wfmash_LG28/gar_LG28.fa \
  -X -p 90 -s 50000 -n 1 -t 4 > test_wfmash_LG28/test.paf 2> test_wfmash_LG28/test.log
wc -l test_wfmash_LG28/test.paf
# Expected: a few hundred PAF rows for one chromosome pair with the LG28 sub-telomeric
# inversion. Open test.paf and grep -c '^LG28' to sanity-check coverage.

# Full run (30 min — both haplotypes are ~1 GB):
bash scripts/run_wfmash_gar_vs_mac.sh
# Output: /scratch/.../cross_species_breakpoints/01_wfmash/{gar_vs_mac.paf, wfmash.log, run_metadata.json}

# Tuning notes:
#   P=85  bash run_wfmash_gar_vs_mac.sh    # if -p 90 produces sparse coverage at LG ends
#   S=25000  bash run_wfmash_gar_vs_mac.sh # for finer breakpoint resolution
#   N=5  bash run_wfmash_gar_vs_mac.sh     # if 1-to-1 misses inversions in repetitive regions
```

If wfmash isn't installed:
```bash
mamba install -c bioconda wfmash
```

## Step 2 — breakpoint extraction

```bash
mkdir -p /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/cross_species_breakpoints/02_breakpoints

python3 scripts/STEP_CS01_extract_breakpoints.py \
  --paf        /scratch/.../cross_species_breakpoints/01_wfmash/gar_vs_mac.paf \
  --te-dir     /scratch/.../het_roh/15_repeat_density/Cgar_TE_density \
  --mac-te-dir /scratch/.../het_roh/15_repeat_density/Cmac_TE_density \
  --candidates /scratch/.../inversion_localpca_v7/06_mds_candidates/inversion_localpca.candidate_regions.tsv \
  --out        /scratch/.../cross_species_breakpoints/02_breakpoints/cs_breakpoints_v1.json
```

(Adjust the `--te-dir` / `--mac-te-dir` paths to wherever your TEfull
JSONs live; if the Mac dir doesn't exist yet, omit `--mac-te-dir`.)

Stderr will show:
- PAF row count and post-filter survivors
- Raw breakpoint count
- Post-clustering count (typically 30-50% reduction depending on radius)
- TE density chrom counts loaded
- Candidate overlap count
- Final n_breakpoints + event-type histogram

## Step 3 — verify

Before trusting the output:

```bash
# Run the synthetic LANTA-side test (no real data required):
cd scripts
python3 test_STEP_CS01.py
# Expected: "[test] PASS — n_breakpoints=3, event_types={...}"
```

The atlas-side renderer has its own JSDOM smoke test
(`scripts/test_atlas_page16.js`) — useful if you've modified the atlas HTML.
Run with:
```bash
# requires jsdom (npm install jsdom — anywhere on disk)
NODE_PATH=/path/to/node_modules node scripts/test_atlas_page16.js
# Expected: "[smoke] ALL CHECKS PASSED" (35 checks across breakpoint + synteny layers)
```

Then sanity-check the real output:
```bash
python3 -c '
import json
data = json.load(open("/scratch/.../02_breakpoints/cs_breakpoints_v1.json"))
print("n_breakpoints:", data["n_breakpoints"])
print("by event type:", data["n_by_event_type"])
print("first 3 breakpoints:")
for bp in data["breakpoints"][:3]:
    print(f"  {bp[\"id\"]:<25} {bp[\"gar_chr\"]:<6} @ {bp[\"gar_pos_mb\"]:>7.2f} Mb  "
          f"{bp.get(\"event_type_refined\", bp[\"event_type\"]):<25}")
'
```

Pick 5 random breakpoints, open the Cgar BAM at those positions in IGV,
and confirm:
- discordant read pairs are present near the breakpoint
- coverage doesn't drop to zero (would indicate a wfmash artifact)
- the breakpoint position roughly matches what IGV shows for split reads

## Step 4 — load into the atlas

Drag `cs_breakpoints_v1.json` into the Inversion Atlas. The "cross-species"
tab will populate. Each breakpoint becomes a row in the catalogue list;
clicking a row loads it into the focus panel with both ideogram strips
and the flanking-repeat density panels.

## Tuning advice (from the handoff)

- Default `-p 90` is appropriate for Cgar↔Cmac (~10-15 Mya divergence,
  expected 92-95% ANI in syntenic regions). Drop to `-p 85` only if
  LG ends look sparse.
- `-s 50000` matches the manuscript scale (most polymorphic candidates
  span 100 kb–10 Mb). Drop to `-s 25000` for finer resolution at the
  cost of more compute and a noisier PAF.
- `-n 1` keeps things clean for the catalogue. If you suspect inversions
  in repeat-rich regions are being missed, rerun with `-n 5` and compare.
- `--cluster-radius-bp 50000` is the default; raise it if you see many
  closely-spaced breakpoints that should be one event, lower it if real
  events are getting merged.

## What the output JSON looks like

See the handoff's `cs_breakpoints_v1` schema. Each breakpoint has:
- `id`, `event_type` (raw), `event_type_refined` (post-context)
- Gar position (start, end, midpoint Mb)
- prev_block + next_block (full Mac coordinates with strand and mapq)
- flanking_repeat_density_gar (mean/max/n_windows per requested class)
- flanking_repeat_density_mac (per side: prev block end, next block start)
- candidate_overlap (list of candidate_id strings, when --candidates given)
- manuscript_note (Spalax-style flag, when local all_TE >= 1.5x chrom mean)

## Schema v2 — macro-synteny extension

`STEP_CS01_extract_breakpoints.py` emits `schema_version: 2`, which adds
four fields to the JSON for the atlas-side macro-synteny analysis:

```jsonc
{
  "n_synteny_blocks": 412,
  "synteny_blocks": [
    { "gar_chr": "LG28", "gar_start": 0, "gar_end": 3000000,
      "mac_chr": "LG28", "mac_start": 0, "mac_end": 3000000,
      "strand": "+", "block_size_bp": 3000000, "mapping_quality": 60 },
    // ... one row per post-filter PAF block, sorted by (gar_chr, gar_start)
  ],
  "chrom_lengths_query":  { "LG01": ..., "LG28": ..., ... },
  "chrom_lengths_target": { "LG01": ..., "LG27": ..., ... }
}
```

These fields drive the **macro-synteny section** below the per-breakpoint
focus panel on page 16:

1. **3-tier Sankey diagram** — Cgar LGs (left) ↔ classification labels
   (middle: `1:1 conserved`, `fusion in target`, `fission in target`,
   `many-to-many`) ↔ Cmac LGs (right). Flow widths are weighted by aligned
   bp. Hover any ribbon for per-flow Mb / % stats.

2. **Per-chromosome relationships table** — for every Cgar chromosome,
   lists its dominant Cmac partners (≥5% aligned-bp threshold) and the
   classification label.

3. **Inversion candidate context table** — for each promoted candidate in
   `state.candidateList`, computes:
   - distance from each breakpoint to the nearest synteny edge
     (left/right shown separately, smaller value bolded)
   - whether the candidate sits inside a single synteny block or spans
     multiple blocks (and how many)
   - the chromosome's macro-synteny class
   - a short evolutionary interpretation phrase (e.g. "inside conserved
     1:1 syntenic block — likely recent within-lineage polymorphism" vs
     "breakpoints near synteny edges on a fusion/fission chrom — possible
     reuse of ancient rearrangement boundary").

4. **Permutation enrichment test** (button at the bottom of the inversion
   table). Test statistic per candidate:
   `dist_min_edge_bp = min(dist from start to nearest edge, dist from end
   to nearest edge)`. The null distribution is built per candidate by
   drawing 10,000 random START positions uniformly on the same chromosome
   (preserving the inversion length), recomputing dist_min_edge for each.
   Per candidate: percentile + smoothed p-value + `near_edge` indicator
   (observed ≤ 25th percentile of null). Globally: one-sided binomial
   test for excess "near edge" hits over the q=0.25 expected rate.
   Interpretation:
   - p<0.05 + enrichment_ratio>1 → "ancient/fragile architecture reuse"
   - p<0.05 + enrichment_ratio<1 → "recent within-lineage polymorphism"
   - else → "mixed signal — inspect per-candidate"

The classification is gar-centric. "fission_in_target" means **1 Cgar →
multiple Cmac partners** as observed; without an outgroup we can't
polarise this as a Cmac fission vs a Cgar fusion. The "in target" suffix
acknowledges this: it's the shape of the relationship as seen with Cmac
as the target reference.

All synteny analysis is derived on-the-fly in the browser — no new file
format. The atlas caches the derived structures (`state._csSyntenyCache`,
`_csSyntenyEdgesCache`, `_csInversionContextCache`) and invalidates them
on JSON reload or `_clearCrossSpecies()`.

A v1 JSON (without `synteny_blocks`) still loads cleanly — the macro-synteny
section is hidden, the per-breakpoint catalogue and focus panels work as
before.

## Tests

```bash
# LANTA-side: synthetic PAF → JSON pipeline
cd scripts && python3 test_STEP_CS01.py
# Expected: "[test] v2 schema check OK" + "[test] PASS"

# Atlas-side: cross-species breakpoints + macro-synteny (35 checks)
NODE_PATH=/path/to/node_modules node scripts/test_atlas_page16.js
# Expected: "[smoke] ALL CHECKS PASSED"

# Atlas-side: stats profile synthesis page (24 checks)
NODE_PATH=/path/to/node_modules node scripts/test_atlas_page17.js
# Expected: "[sp-smoke] ALL CHECKS PASSED"

# Atlas-side: marker readiness panel (26 checks, private-indel architecture)
NODE_PATH=/path/to/node_modules node scripts/test_atlas_page18.js
# Expected: "[mp-smoke] ALL CHECKS PASSED"
```

## Marker readiness panel (page 18)

A new tab in the atlas dedicated to the **per-inversion marker confidence
tier table**, rebuilt around **private indels with clean dosage signatures**
\u2014 not breakpoint PCR. Designed so a 500-THB PCR test risk is small: users
start with a Tier 1 pilot panel of validated private markers, with positive
and negative controls in every PCR run.

### Why private indels, not breakpoint PCR

Breakpoint PCR markers are not as reliable as previously assumed because
breakpoint precision itself is uncertain \u2014 even the cs_breakpoint catalogue
can mis-localise breakpoints by 10s of kb. A primer pair designed against a
mis-localised breakpoint will fail. Private indels with a clean
allele-frequency dosage signature (AF_STD\u00A0\u22480, AF_HET\u00A0\u22480.5,
AF_INV\u00A0\u22481) are far more robust: they tag the structural haplotype
itself, are gel-visible if the size difference is in [20,300] bp, and don't
depend on knowing exactly where the inversion ends.

### Four tiers

- **Tier 1 \u2014 private indel/SNP tag (highest confidence)**. Clean dosage:
  AF_STD\u00A0\u22640.02, AF_HET\u00A0\u2208\u00A0[0.25,\u00A00.75],
  AF_INV\u00A0\u22650.80. Plus positive controls in INV/INV and HET groups,
  negative controls in STD/STD, and bighead specificity confirmed (the marker
  must NOT amplify the bighead orthologous region in the same way \u2014 critical
  because the breeding system uses interspecific F\u2081 hybrids).
- **Tier 2 \u2014 multi-marker panel OR strong tag (medium)**. Either \u22653 linked
  private markers across the inversion (left/middle/right \u2014 if all agree,
  the haplotype call is robust) OR a single strong tag
  (AF_STD\u00A0\u22640.05, AF_INV\u00A0\u22650.70) with controls but imperfect
  het.
- **Tier 3 \u2014 breakpoint PCR candidate (DEMOTED)**. Breakpoint within
  \u00B150 kb of a cs_breakpoint, controls identified, wet-lab validation
  pending. Lower confidence than private indels because breakpoint precision
  is uncertain.
- **Tier 4 \u2014 exploratory (lowest)**. Association exists but controls
  incomplete or weak. For research-only follow-up; do not use for breeding.

### How tiers get assigned

The atlas computes tiers in two stages:

1. **Atlas auto-tier** (no AF data): every promoted candidate starts at
   Tier 4, except those within \u00B150 kb of a `cs_breakpoint` which start at
   Tier 3 (breakpoint PCR candidate).
2. **AF promotion** (after `variant_afs.json` is loaded): each variant is
   scored \u2014 `private_score = AF_INV - AF_STD`,
   `dosage_score = 1 - |AF_HET - 0.5|*2`, `final_score` is their sum
   (range 0\u20132, \u22651.7 = clean). The variant's AF-derived tier promotes the
   candidate's row tier when better. \u22653 Tier-1 variants on one inversion
   triggers a "multi-marker panel feasible" evidence note.
3. **User override** (after `marker_panel.json` is loaded): primer
   sequences, amplicon sizes, validation status, cross-species control
   results. The override can also explicitly set the tier; the original
   auto/AF tier is preserved in the evidence trail for transparency.

### Auto-suggested controls from karyotype

When a candidate carries `assignments[sample] = 0|1|2` (the K-means PCA
karyotype assignments from page 1), the panel auto-suggests up to 3
positive-INV, 3 heterozygote, and 3 negative-STD control samples per
marker. Cross-species control fields (bighead negative status, F\u2081 hybrid
expected pattern, no-template control requirement) are populated from a
default template and overridable via `marker_panel.json`.

### What the page shows

1. **Tier definition cards** \u2014 four side-by-side cards (T1, T2, T3, T4)
   with their AF/control criteria and expected use.
2. **Summary cards** \u2014 count per tier + count with controls identified +
   count validated/pilot-tested.
3. **Filter toolbar** \u2014 by tier, validation status, chromosome; load
   `variant_afs.json`, load `marker_panel.json`, export CSV, reset.
4. **Main table** \u2014 one row per inversion candidate. Columns: inversion_id,
   chrom \u00b7 position, tier, marker type, **top scored variants (AF)** with
   STD/HET/INV breakdown and gel-visibility badges, **controls** column with
   auto-suggested INV/HET/STD samples and cross-species note, validation
   status pill, tier reason.
5. **Pilot validation plan** \u2014 5-step checklist updated for cross-species
   controls:
   1. Run each marker on positive controls (INV/INV and HET samples).
   2. Run on negative controls (STD/STD) and a no-template control (water).
   3. Run on bighead catfish DNA to confirm cross-species specificity.
   4. Test on 20\u201330 unrelated broodstock; keep markers with clean separation.
   5. Expand to larger populations.
6. **Methods note** \u2014 includes both protective sentences:
   *"Because the target breeding system involves interspecific F\u2081 hybrids,
   candidate markers were additionally screened against the bighead catfish
   orthologous genome to distinguish African inversion tags from
   maternal-species sequence differences."* AND *"Candidate markers are
   provided as a starting resource and should be locally validated before
   routine breeding deployment."*

### File formats

**`variant_afs.json`** (drives the AF tier promotion):
```json
{
  "metadata": {...},
  "variants_by_inversion": {
    "CAND_LG28_subtelo": [
      { "chr": "LG28", "pos": 9234567, "ref": "A", "alt": "ATTGGCC...",
        "type": "indel", "indel_size_bp": 60,
        "af_std": 0.01, "af_het": 0.48, "af_inv": 0.96 }
    ]
  }
}
```

**`marker_panel.json`** (overlay with primers + controls):
```json
{
  "metadata": {...},
  "markers": [
    {
      "inversion_id": "CAND_LG28_subtelo",
      "tier": 1,
      "validation_status": "pilot_tested",
      "primer_F": "CTGAACGGAATGCAACTGGT",
      "primer_R": "TGGCAACTGAATGCGAACGT",
      "amplicon_bp_state_A": 240,
      "amplicon_bp_state_B": 180,
      "n_carriers_tested": 24,
      "controls": {
        "positive_controls_INV": ["CGA021", "CGA044"],
        "heterozygote_controls": ["CGA005"],
        "negative_controls_STD": ["CGA009", "CGA010"],
        "bighead_negative_status": "tested_no_amplification",
        "bighead_orthologous_sequence_status": "absent",
        "f1_hybrid_expected_pattern": "..."
      }
    }
  ]
}
```

Demo files at `variant_afs.demo.json` and `marker_panel.demo.json` in this
directory.

### Stats profile integration

Page 17's `markerability` row reads the live promoted panel from
`state._markerPanel.panel` (so it reflects AF promotion + user overrides)
and reports the tier breakdown directly:
*"Tier 1: 5 \u00B7 Tier 2: 12 \u00B7 Tier 3: 8 \u00B7 Tier 4: 10 (out of 35 markers,
live AF-aware promotion; 22 confirmed candidates) \u2014 see page 14 marker
panel"*. The headline metric is the count of Tier 1 markers.

## Stats profile synthesis page (page 17)

A new tab in the atlas, sitting beside page 16 cross-species. It answers
**"what is statistically special about inversion regions, and are the patterns
shared between Cgar and Cmac?"** in a single comparative table.

The 12 default rows mirror the comparative-inversion-phenotype mockup:

1. **Breakpoint architecture** — distance to nearest synteny edge; fusion/fission proximity
2. **Genomic composition** — repeat density flanking breakpoints; telomere/centromere/gap proximity; gene density inside inversions
3. **Functional cargo** — GO/KEGG/Pfam enrichment; breeding-candidate genes
4. **Population variation** — heterozygosity, ROH overlap, deleterious burden, FST/dXY between haplotypes
5. **Breeding utility** — markerability (count of confirmed inversions ready for marker design)

Each row has a six-column layout (detailed) or four-column layout (compact):
*category · statistic/test · null comparison · Cgar result · Cmac result · biological interpretation.*

### How rows get populated

The page operates in three states per cell:

- **`atlas`** — the row is computed live from atlas state. Currently auto-derived: breakpoints near synteny edges (runs the page 16 length-matched permutation test), fusion/fission proximity (counts cs_breakpoints by event type), repeat density flanks (counts Spalax-flagged breakpoints + reports Mac-side mean), markerability (counts confirmed candidates).
- **`loaded`** — the row was filled from a stats_profile JSON or TSV upload. Use this for rows the atlas can't compute itself: GO/KEGG enrichment, ROH overlap, deleterious burden, FST between haplotypes, gene density (requires GFF), telomere/centromere proximity (requires annotation track).
- **`add data`** — placeholder. The cell shows a hint pointing to the analysis or pipeline that would populate it (e.g. *"requires SIFT4G / VESM_650M output (MODULE_CONSERVATION)"*).

### File formats

The page accepts two input formats via the *load JSON / TSV* button:

**JSON** (canonical):
```json
{
  "metadata": { "species_1": "C. gariepinus", "species_2": "C. macrocephalus", "date": "..." },
  "summary_rows": [
    {
      "id": "gene_density_inversions",
      "african_result": { "effect": "lower_than_background", "observed": 14.2, "expected": 18.7,
                          "p_value": 0.012, "n": 8, "label": "..." },
      "bighead_result": { "effect": "not_significant", ... },
      "interpretation": "..."
    }
  ]
}
```

**TSV** (for spreadsheet workflows):
```
category	statistic	test	null_comparison	african_label	african_effect	african_observed	african_expected	african_p	african_q	bighead_label	bighead_effect	bighead_observed	bighead_expected	bighead_p	bighead_q	interpretation
```

Use empty string or `NA` for missing values. The TSV parser autodetects which
columns are present from the header row.

A complete demo file is at `stats_profile.demo.json` in this directory \u2014 drag
it into the atlas to see all 7 sample rows render with realistic placeholder
values (gene density, KEGG enrichment, breeding-candidate genes, het, ROH,
deleterious burden, FST). Replace each result block with values from your own
analyses.

### Effect tags

Subtle palette tied to the rest of the atlas:
- `enriched` / `closer_than_expected` / `higher_than_background` / `common` \u2014 green (#3cc08a)
- `depleted` / `farther_than_expected` / `lower_than_background` / `rare` \u2014 amber (#f0a35e)
- `not_significant` \u2014 muted grey
- `not_tested` \u2014 darker grey
- `african_only` / `cross_species` \u2014 cyan / purple scope tags

### Filters and exports

- View mode: detailed (6 cols) / compact (4 cols)
- Filter by category, scope (all / African-only / cross-species), significant-only
- Export current view as CSV (full 17 columns regardless of view mode)
- Export rendered table as SVG (single-file, manuscript-quality, foreignObject-wrapped)

### Methods note shown in-page

> Statistics are compared against chromosome-, size-, or feature-matched genomic backgrounds where possible. Population-based tests are performed in the African catfish broodstock panel; bighead catfish values refer to orthologous genomic context unless population resequencing data are available. Atlas-derived rows recompute live from `state.crossSpecies` and `state.candidateList`; rows requiring annotation tracks (gene density, GO/KEGG, ROH, deleterious burden, FST) accept a stats_profile.json or .tsv overlay.
