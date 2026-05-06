# Inversion Atlas — handoff document (after turn 127)

**Last commit point:** turn 127 — header reorganization with three colour-coded folders
**Atlas size:** 62,023 lines (`Inversion_atlas.html`)
**Test status:** 762/762 green across 12 test suites

---

## You are

Quentin Andres, PhD researcher in computational genomics and aquaculture
genetics at Kasetsart University, Bangkok. Working on LANTA HPC (account
`lt200308`). Primary manuscript:
**MS_Inversions_North_african_catfish**, target: Nature Communications.

You are direct and terse. You work across multiple parallel chat sessions
simultaneously and have strong scientific intuition. French native, English
fluent. You want "read first, propose second" before coding, no code
written for methods not yet implemented.

## Three catfish cohorts (NEVER conflate)

1. **F1 hybrid** (C. gariepinus × C. macrocephalus) — ONLY for the genome
   assembly paper (a separate manuscript).
2. **226-sample pure C. gariepinus hatchery cohort on LANTA** — current
   work, the inversion manuscript. K clusters reflect hatchery broodline
   structure, NOT species admixture.
3. **Pure C. macrocephalus wild cohort** — future paper.

The atlas's main page works on cohort #2 (226 samples on LG28).

---

## Atlas architecture (high level)

`/home/claude/work/build/Inversion_atlas.html` is a single 62k-line
self-contained HTML application. Pages 1–17 are the manuscript-driving
views; the multi-species cockpit at page 16b is the latest synthesis page.

Key infrastructure constants:
- LANTA path: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/`
- Reference: `fClaHyb_Gar_LG.fa` (28 pseudochromosomes)
- NAToRA-pruned 81 unrelated samples
- NGSadmix K=8
- KING/Manichaikul thresholds; `mergeThr 0.85`; `purity_threshold 0.80`;
  `state.k 3` (primary karyotype regime)
- Scrubber JSON layers in `scrubber/data/<chrom>/`

State management uses a global `state` object + localStorage namespace
`inversion_atlas.<key>` for persistence.

---

## Recent turns 121–127 (what shipped session by session)

### Turn 121 — page 16b multi-species cockpit (+1,068 lines, 110 tests)
Tab `13b multi-species`, three-column page (clickable phylo tree left /
classification chips center / dotplot+detail right). State slots
`state.syntenyMultispecies`, `state.phyloTree`,
`state._multiSpeciesUI.active_species`. Layers:
`synteny_multispecies_v1`, `phylo_tree_v1`. Architecture-class
auto-suggest A/B/C/D/E/F. `_MS_DEFAULT_SPECIES` 10-species fallback (Tros,
Smer, Tfulv, Ipun, Hwyc, Phyp, Capus, Cfus, Cmac, Cgar). Min Newick
parser, custom SVG tree, no library. Active breakpoint from
`state._crossSpeciesUI.active_id`.

### Turn 122 — age model + override + TSV export
`dxy_per_inversion_v1` layer, age-model auto-suggest
(MULTI-AGE-HOTSPOT > LINEAGE-KARYO > OLD-BP-YOUNG-INV > OLD-POLY > YOUNG-POP),
manual override editor, `_msBuildClassificationTSV`, TE-fragility
focal-vs-bg slot in detail column.

### Turn 123 — comparative TE fragility (+218 lines, 73 tests)
`comparative_te_breakpoint_fragility_v1` layer. Per-species TE density at
homologous region in focal radius vs chrom-wide median. Chip thresholds:
red ≥ 2.5×, orange 1.5–2.5×, neutral ~1×, depleted < 1×. On-screen
disclaimer: "Fragility proxy only — not evidence of polymorphism in
non-focal species."

### Turn 124 — karyotype lineage layer (+427 lines, 89 tests)
`karyotype_lineage_v1` JSON layer from mashmap one-to-one (1 Mb / 85%
identity) all-vs-all. Six polarization verdicts. Drives RULE 0 in
age-model auto-suggest at high confidence. Lineage table merges
karyotype-class column. Honest disclaimer about resolution. LANTA pipeline
at `/project/lt200308-agbsci/01-catfish_assembly/05_ancestral_karyotype/mashmap_fusion_scan_lines/`
(`mashmap_precurl_ancestral_events.slurm` + `get_best_refs_and_1to2.sh`).

### Turn 125 — species-agnostic refactor (+106 lines, 50 tests)
Hardcoded Cgar/Cmac removed from karyotype polarization. Schema v2:
`focal_species`, `sister_species[]`, `outgroup_species[]`,
`per_focal_chr[]`. Verdict keys: `focal_lineage_fission`,
`sister_lineage_fission` (legacy `cgar_/cmac_` accepted as backward-compat
aliases). Verdict label substitutes actual species name at render time
("Cgar lineage fission"). Polarization rules read roles from JSON.
focal_species inferred from data when not declared. Atlas now reusable
for Cmac paper without code changes.

### Turn 126 — wfmash refinement (+248 lines, 76 tests)
Per-cell refinement fields: `refined_by_wfmash`
(`confirmed`|`refuted`|`refined`|`failed`|null), `wfmash_class`,
`wfmash_targets[{chrom,start_bp,end_bp,strand}]`, `wfmash_n_blocks`,
`wfmash_pct_identity`. Display: refinement summary line above table,
refinement column with chips, strikethrough on refuted classes.
`_msAdjustConfidenceForRefinement`: ≥ 80% confirmed AND 0 refuted →
boost (medium → high); ≥ 30% refuted → drop; < 50% touched → unchanged;
unresolved verdicts never adjusted. New helpers: `_msSummarizeRefinement`,
`_msPolarizeKaryotypeEventWithRefinement`,
`_msGetEffectiveClassForCell/Targets`, `_msBuildRefinementChipHtml`.

### Turn 127 — header folders (+99 lines, 67 tests) — JUST SHIPPED
Top-row buttons grouped into three colour-coded folders:
- 🔵 **session** (blue): save · load · reset layout · fixed/free
- 🟢 **mode** (green): candidate mode · auto-fill · I·g labels
- 🟡 **data** (yellow): export · matrix · active · server

Right-side row of buttons (`atlasToolsGroup`) moved out of the tab bar
into the yellow folder. Metadata strip trimmed (no L1/L2/peaks).
Hover-reveal pattern mirrors `#atlasModeIndicator`.

---

## File system at handoff

In `/home/claude/work/build/`:

```
Inversion_atlas.html                                 62,023 lines
test_turn117_integration.js                          15 tests
test_turn118_enriched_bundle.js                      82 tests
test_turn119_age_origin.js                           74 tests
test_turn120_scree_inset.js                          33 tests
test_turn121_multispecies.js                        110 tests
test_turn122_age_override_export.js                  79 tests
test_turn123_te_fragility.js                         73 tests
test_turn124_karyotype_lineage.js                    89 tests
test_turn125_species_agnostic.js                     50 tests
test_turn126_wfmash_refinement.js                    76 tests
test_turn127_header_folders.js                       67 tests
test_dotplot_hover_only.js                           14 tests
                                                    ───
                                                    762 tests
```

Plus LANTA-side glue:
```
build_karyotype_lineage_v1_json.py    converter for karyotype lineage layer (turn 124-125)
refine_karyotype_lineage_with_wfmash.py    wfmash refinement runner (turn 126)
species_map.tsv                       genome basename → atlas species id mapping
```

Demo JSONs:
```
synteny_multispecies_v1.demo.json
phylo_tree_v1.demo.json
dxy_per_inversion_v1.demo.json
comparative_te_breakpoint_fragility_v1.demo.json
karyotype_lineage_v1.demo.json    (v2 schema with refinement fields populated)
```

---

## How to continue in a new chat

### To run the test suite:
```
cd /home/claude/work/build
for t in test_turn117_integration.js test_turn118_enriched_bundle.js \
         test_turn119_age_origin.js test_turn120_scree_inset.js \
         test_turn121_multispecies.js test_turn122_age_override_export.js \
         test_turn123_te_fragility.js test_turn124_karyotype_lineage.js \
         test_turn125_species_agnostic.js test_turn126_wfmash_refinement.js \
         test_turn127_header_folders.js test_dotplot_hover_only.js; do
  result=$(node "$t" 2>&1 | grep -E 'PASSED:|passed:.*failed:' | tail -1)
  echo "$t  ::  $result"
done
```

Expected: 762/762 green.

### Key idioms used by tests

Source-level checks use regex against the raw HTML:
```js
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');
ok('feature X wired',
   /someUniqueAnchor/.test(html));
```

Behavioural tests use a vm sandbox with extracted functions:
```js
function pullFunction(src, fnName) { /* brace-balanced extractor with
  comment + string awareness. Copy from test_turn126_wfmash_refinement.js */ }
```

The sandbox needs a fake `state` object, fake `_esc`, and any constants the
extracted functions reference. See `test_turn126_wfmash_refinement.js` lines
77-180 for the canonical recipe.

### Where things live in the atlas

- **Header bar markup:** lines 4509–4690 (`<header>...</header>`)
- **Tab bar:** lines 4691–4830 (`<nav id="tabBar">`)
- **CSS for header folders:** lines 195–290 (the `header .header-folder...`
  rules)
- **Multi-species karyotype layer:** functions live around lines
  23062–23615 (`_isKaryotypeLineageJSON` through
  `_msBuildKaryotypeContextHtml`)
- **Age-model auto-suggest:** function `_msAutoSuggestAgeModel` around
  line 21750. RULE 0 (LINEAGE-KARYO override) starts at line 21895.
- **Server status badge:** dynamically appended by
  `_atlasServerInitBadge()` at line 48136. Append target:
  `document.getElementById('atlasToolsGroup')`.
- **Metadata strip build:** line 41635
  (`document.getElementById('headerMeta').innerHTML = ...`).

### Things the atlas marks as TODO

1. **TSV export columns** for the page 16b classification TSV
   (`_msBuildClassificationTSV`) — should add columns for the karyotype
   verdict and wfmash refinement state. Columns to add:
   `karyo_verdict`, `karyo_confidence`, `wfmash_n_confirmed`,
   `wfmash_n_refuted`, `wfmash_n_failed`. The function is around line
   22700 area; I haven't shipped this yet but it's small.
2. **Active-folder-glow JS hookup** (turn 127 follow-up): when
   `candidateModeBtn[data-active="1"]` flips, the parent green folder
   should get `data-has-active="1"` so the CSS ring rule fires. The CSS
   rule already exists. Find the `setCandidateMode` function (or its
   equivalent toggle handler) and add a sibling line. Trivial.
3. **Page 16c chromosome-karyotype overview**: a 28-LG grid showing
   verdict chips for every Cgar chromosome. Separate scope from page 16b
   (which is per-breakpoint). Mentioned in earlier sessions.
4. **Run actual chromosome-scale wfmash refinement on LANTA**: the
   `refine_karyotype_lineage_with_wfmash.py` script is functional but
   nobody has run it on real data yet. Multi-hour job.
5. **BUSCO integration**: you said you have BUSCO data to show but it's
   "quite challenging — first you need 2-3 more turns right to solve
   that. then only then we can try BUSCO." Now turns 125, 126, 127 have
   shipped, so you said the next session you'll show me BUSCO data.

### What you'll likely ask next

Given the sequence of recent turns, the most likely next request is one of:

1. **BUSCO integration** — the big one you've been holding back.
2. **TSV export columns** — small finishing touch on the manuscript export.
3. **Page 16c overview** — the 28-LG karyotype matrix view.
4. **A bug or polish issue you spot in the new header folders** — possible
   if the popovers don't behave on real data.

In all cases: ask to see the data first, agree on schema, then ship the
layer + tests.

---

## Critical rules (DO NOT BREAK)

1. **Never conflate the three cohorts.** F1 hybrid is the genome assembly
   paper, NOT the inversion paper. The 226-sample cohort is pure
   C. gariepinus hatchery, K-clusters reflect broodline structure, not
   species admixture. C. macrocephalus is a future paper.
2. **Never invent a surname for Quentin.** The full name is exactly
   "Quentin Andres".
3. **Vocabulary discipline (chat 1b4d8e12):**
   - "polymorphic in Gar" allowed only with 226-sample carrier evidence
   - "candidate fragile region in species X" / "predicted polymorphic
     hotspot" allowed for non-resequenced species
   - NEVER say "structural polymorphism in species X" when X is not the
     focal cohort
4. **Don't write code for methods not yet implemented.** Read first,
   propose second.
5. **Tests must use source-level regex + vm-sandbox extraction**, not
   require + jsdom. The atlas is a single 62k-line HTML file with a giant
   IIFE — too big to load as a module.
6. **Modular factory pattern.** Reuse `popgenDotplot` and
   `popgenFocalVsBg` for new ribbon/PCA panels rather than reinventing.
7. **Honest claims.** Polarization is at 1 Mb resolution — say so. wfmash
   refinement is at 50 kb — say so. No magic.

---

## Manuscript framing carried forward

Two-layer classification (page 5 help):

**Architecture class** — what the breakpoint looks like:
- A = simple inversion
- B = synteny-boundary
- C = fusion/fission-associated
- D = terminal translocation
- E = recurrent rearrangement hotspot ("strongest Nature-style mechanism class")
- F = ambiguous

**Age model** — when the breakpoint formed:
- YOUNG-POP — recent population polymorphism
- OLD-POLY — ancient polymorphism, multi-haplotype
- OLD-BP-YOUNG-INV — old breakpoint reused by young inversion
- LINEAGE-KARYO — fixed in lineage, karyotype-level event
- MULTI-AGE-HOTSPOT — recurrent reuse, multiple ages

The five Result-section architecture for the manuscript:
1. Find inversions
2. Define structural haplotypes
3. Design PCR / marker assays
4. Enable breeding
5. (Comparative / lineage context — the work shipped in turns 121-126)

---

## Schemas (the JSON layers the atlas reads)

### `karyotype_lineage_v1` (current schema_version: 2)

```jsonc
{
  "tool": "karyotype_lineage_v1",
  "schema_version": 2,
  "params": {
    "method": "mashmap_one_to_one",
    "segment_bp": 1000000,
    "identity_pct": 85,
    "wfmash_refinement": {  // optional, populated by refine script
      "method": "wfmash",
      "segment_bp": 50000,
      "identity_pct": 90
    }
  },
  "focal_species": "Cgar",
  "sister_species": ["Cmac"],
  "outgroup_species": ["Tros"],
  "per_focal_chr": [
    {
      "focal_chr": "C_gar_LG28",
      "classes_by_species": {
        "Cgar": {
          "class": "1-1",
          "targets": ["C_gar_LG28"],
          "refined_by_wfmash": "confirmed",       // optional
          "wfmash_class": "1-1",                  // optional
          "wfmash_targets": [...],                // optional
          "wfmash_n_blocks": 28,                  // optional
          "wfmash_pct_identity": 100.0            // optional
        },
        ...
      }
    }
  ]
}
```

Backward compat: `per_cgar_chr` accepted as alias for `per_focal_chr`;
`cgar_chr` accepted as alias for `focal_chr`. focal_species inferred
from data when not declared.

### `synteny_multispecies_v1` (turn 121)

Per-breakpoint, per-species coordinate-resolved alignments. Ribbon view
on the dotplot. See `synteny_multispecies_v1.demo.json` for shape.

### `phylo_tree_v1` (turn 121)

Newick tree string + optional per-node metadata. Parsed by minimal
in-atlas Newick parser (no library).

### `dxy_per_inversion_v1` (turn 122)

Per-(candidate, species-pair) dXY estimates. Drives age-model
auto-suggest.

### `comparative_te_breakpoint_fragility_v1` (turn 123)

Per-species, per-breakpoint TE density at the homologous region vs
chromosome-wide median. Strict on vocabulary: "fragility proxy" only.

---

## Communication style notes

- Direct and terse responses. Don't over-explain.
- Acknowledge mistakes briefly, fix them, move on.
- "Read first, propose second" — for any non-trivial change, look at the
  affected source before patching.
- Tarballs for multi-file deliverables. GitHub Desktop drag-drop on
  Windows; HPC clone on LANTA.
- Handoff documents inside each tarball.
- `STEP_` naming convention with central config files for R/python
  pipeline scripts.

---

## TL;DR for the next chat

```
Atlas is at /home/claude/work/build/Inversion_atlas.html — 62,023 lines.
762 tests across 12 suites all green.
Last turn (127) reorganized the header into three colour-coded folders.
Next likely request: BUSCO integration, TSV export columns, or page 16c.
Run tests with the loop in this doc. Add new turn = new test_turnNNN.js.
Modular factory pattern. Source-level regex + vm-sandbox tests. Vocabulary
discipline. Three cohorts never conflated.
```
