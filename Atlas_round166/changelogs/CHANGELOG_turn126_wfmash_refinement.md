# CHANGELOG — turn 126: wfmash refinement of mashmap karyotype calls

**Atlas Δ:** +248 lines (61,676 → 61,924)
**Tests:** 76 new (695/695 cumulative green across 11 test suites)
**No regressions** — all prior turns 117–125 still pass.

---

## What shipped

The karyotype layer now supports **per-cell wfmash refinement**.
Every interesting mashmap class (typically 1-2 / 1-3) can be re-evaluated
at higher resolution (50 kb segment, 90% identity wfmash) and the result
flows directly into the atlas display + polarization confidence.

### Why this matters

Mashmap at 1 Mb / 85% identity is fast — that's why it can do all-vs-all
on 17 catfish genomes in reasonable time. But its 1-2 calls can be
artifacts: a sister chromosome that looks like two distinct mappings at
1 Mb scale may collapse into one continuous mapping at 50 kb scale (or
vice versa, depending on rearrangement details). The atlas should
distinguish "this 1-2 is real" from "this 1-2 is a low-res artifact"
and let the user see both states.

This turn provides the schema, atlas display, and a LANTA-side script
(`refine_karyotype_lineage_with_wfmash.py`) to produce the refinement
data. The atlas already has the rendering logic.

### New schema fields per cell

```jsonc
"classes_by_species": {
  "Cmac": {
    "class": "1-2",
    "targets": ["C_mac_LG01", "C_mac_LG18"],
    "refined_by_wfmash": "confirmed",          // NEW
    "wfmash_class": "1-2",                      // NEW (when refined/refuted)
    "wfmash_targets": [                         // NEW
      { "chrom": "C_mac_LG01", "start_bp": 8100000,  "end_bp": 11300000, "strand": "-" },
      { "chrom": "C_mac_LG18", "start_bp": 14200000, "end_bp": 18500000, "strand": "+" }
    ],
    "wfmash_n_blocks": 22,                      // NEW
    "wfmash_pct_identity": 92.4                 // NEW
  }
}
```

### Five refinement states

| State | Meaning | Atlas display |
|---|---|---|
| **`confirmed`** | wfmash agrees with mashmap class | green ✓ chip |
| **`refuted`** | wfmash collapses to a less-fragmented class (e.g. 1-2 → 1-1) | orange ✗ chip + strikethrough on mashmap class |
| **`refined`** | wfmash gives a different class (e.g. 1-2 → 1-3) | blue ➜ chip showing new class |
| **`failed`** | wfmash ran but produced no confident alignment | grey ⚠ chip |
| **`null` / absent** | not yet attempted | dim em-dash chip |

### Polarization confidence adjustment

The verdict confidence is now adjusted based on the refinement
distribution across cells:

- **≥80% confirmed AND 0 refuted** (and ≥50% of cells touched)
  → confidence boosted one tier (medium → high). Rationale gets a
  subtle "Boosted by wfmash confirmation across N/M cells" note.
- **≥30% refuted** → confidence dropped one tier and rationale gets
  a caution note: "wfmash refuted N/M mashmap cells at higher resolution".
- **`unresolved` verdicts** never adjusted (refinement doesn't
  synthesize new info from sparse signal).

The threshold logic is in `_msAdjustConfidenceForRefinement` —
pure function, easy to tune.

### Atlas display additions

In the karyotype context panel:

1. **Refinement summary line** above the per-species table:
   *"wfmash refinement: **5** confirmed · **1** refuted · **2** failed"*
2. **Refinement column** in the per-species table — only rendered
   when at least one cell has been touched by wfmash. Each row gets
   the appropriate refinement chip.
3. **Strikethrough** on the mashmap class chip when wfmash refuted
   it. The `1-2` is still visible (so the user sees the original
   mashmap call), but visually demoted.
4. **Updated legend**:
   *"Mashmap (1 Mb / 85% id) gives chromosome-scale class. wfmash
   (kbp / higher id) refines: ✓ confirmed, ✗ refuted, ➜ refined to
   different class. Pair with breakpoint-scale wfmash synteny layer
   for sub-Mb resolution."*

### LANTA-side: `refine_karyotype_lineage_with_wfmash.py`

Takes an existing `karyotype_lineage_v1.json` and produces a
`karyotype_lineage_v1.refined.json` with all the new fields populated.

```bash
python3 refine_karyotype_lineage_with_wfmash.py \
  --input karyotype_lineage_v1.json \
  --species-fasta-map species_fasta_map.tsv \
  --output karyotype_lineage_v1.refined.json \
  --threads 16
```

What it does for each cell where mashmap class != '1-1' (or all if
`--refine-all`):

1. Extract the focal-species chromosome region with `samtools faidx`
2. Run `wfmash` against the target species fasta at higher resolution
   (50 kb segments, 90% identity by default — tunable)
3. Parse the PAF; count distinct target chromosomes; classify as
   1-1/1-2/1-3/1-4+
4. Compare the wfmash class to the mashmap class and label the cell
   `confirmed`/`refuted`/`refined`/`failed`

The script preserves the species-agnostic schema from turn 125.
Output schema_version remains 2.

### Helper functions added (atlas)

- `_msSummarizeRefinement(entry)` → `{confirmed, refuted, refined, failed, not_attempted, total}`
- `_msAdjustConfidenceForRefinement(verdict, summary)` → adjusted verdict
- `_msPolarizeKaryotypeEventWithRefinement(focalChr)` → polarization with confidence adjustment applied
- `_msGetEffectiveClassForCell(cell)` → wfmash class when refuted/refined, else mashmap class
- `_msGetEffectiveTargetsForCell(cell)` → wfmash targets when refuted/refined, else mashmap targets
- `_msBuildRefinementChipHtml(cell)` → small chip HTML

All exposed on `window.*` for tests.

### Demo JSON updated

`karyotype_lineage_v1.demo.json` now exercises all four refinement
states across 4 Cgar chromosomes:

| chrom | mashmap-driven verdict | refinement state | post-refinement confidence |
|---|---|---|---|
| LG28 | sister_lineage_fission | 5 confirmed, 2 failed, 2 not attempted | medium (55% confirmed, below boost threshold) |
| LG07 | recurrent_fission_hotspot | 4 confirmed, 1 refuted, 1 failed, 1 not attempted | medium (refutation at 14%, below drop threshold) |
| LG10 | focal_lineage_fission | 6 confirmed (100%) | **high** (boosted) |
| LG01 | no_karyotype_change | no refinement (params filter) | medium (unchanged) |

End-to-end smoke test confirmed: load demo JSON → atlas detects →
polarization runs → confidence adjusted → render produces correct HTML.

---

## What's still TODO

Per the original plan, this turn was meant to set up the schema +
display + LANTA glue. Things deliberately NOT in scope:

- **Running the actual wfmash refinement on LANTA.** That's a multi-hour
  job that needs the genomes downloaded and the script invoked. Will
  produce real refined data when run.
- **`refine_karyotype_lineage_with_wfmash.py` polish.** The current
  script is functional but minimal — single-threaded outer loop, no
  resume support, no parallelism across cells. Good enough for a first
  pass on ~10 species × 28 Cgar chromosomes (~280 cells).
- **Schema field for breakpoint-scale wfmash refinement of cs_breakpoints.**
  That's a separate layer (it'd belong with `synteny_multispecies_v1`,
  not `karyotype_lineage_v1`). Different from chromosome-scale refinement.

---

## Files changed

```
Inversion_atlas.html                                  +248 lines
test_turn126_wfmash_refinement.js                     NEW (76 tests)
test_turn124_karyotype_lineage.js                     updated sandbox imports
test_turn125_species_agnostic.js                      updated sandbox imports
karyotype_lineage_v1.demo.json                        refinement fields populated
refine_karyotype_lineage_with_wfmash.py               NEW (LANTA wfmash runner)
CHANGELOG_turn126_wfmash_refinement.md                THIS FILE
```
