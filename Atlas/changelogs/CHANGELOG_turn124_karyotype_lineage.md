# CHANGELOG — turn 124: karyotype-context lineage layer (page16b)

**Atlas Δ:** +427 lines (61,068 → 61,495)
**Tests:** 88 new (568/568 cumulative green across 9 test suites)
**No regressions** — all prior turns 117–123 still pass.

---

## What shipped

A new JSON layer **`karyotype_lineage_v1`** plus a karyotype-context
panel at the **top** of page 16b's center column. This wires up the
chromosome-scale evidence from your existing
`mashmap_fusion_scan_lines/` pipeline into the atlas.

### The motivation, in one paragraph

Your already-running mashmap one-to-one all-vs-all scan
(`mashmap_precurl_ancestral_events.slurm` →
`get_best_refs_and_1to2.sh`) produces per-(query, ref) summaries
classifying each chromosome's mapping as **1-1**, **1-2**, **1-3**,
or **1-4+**. This is *chromosome-painting resolution* (1 Mb / 85% identity)
— it cannot resolve breakpoint coordinates, but it **does** answer
the question "did Cgar LG28 fission/fuse on the Cgar lineage,
the Cmac lineage, ancestrally, or recurrently?" That polarization
verdict was previously not in the atlas; now it drives the
LINEAGE-KARYO and MULTI-AGE-HOTSPOT age-model classifications.

### Two scales, side by side

The center column of page 16b now shows evidence at two scales:

| Scale | Layer | Resolution | Question |
|---|---|---|---|
| **Chromosome (top)** | `karyotype_lineage_v1` | 1 Mb / 85% identity | Did this Cgar LG fission/fuse on a specific lineage? |
| **Breakpoint (middle)** | `synteny_multispecies_v1` | wfmash kbp | Is the boundary present at this exact position in each species? |
| **TE fragility (bottom)** | `comparative_te_breakpoint_fragility_v1` | per-species, focal radius | Is the homologous region TE-rich (fragility proxy)? |

The karyotype panel sits **first** because it's the chromosome-scale
context that frames everything else. If the verdict is
"no_karyotype_change" then the breakpoint is intra-chromosomal — an
inversion. If the verdict is "cgar_lineage_fission" then it's a
karyotype-level event and the age-model classification jumps to
**LINEAGE-KARYO** with high confidence.

### Polarization rules (six verdicts)

| Verdict | Trigger | Manuscript implication |
|---|---|---|
| **`cgar_lineage_fission`** | Cgar 1-2, Cmac 1-1, ≥1 outgroup 1-1, ≥2 species 1-1 | Cgar-lineage karyotype event |
| **`cmac_lineage_fission`** | Cmac 1-2, Cgar 1-1, ≥1 outgroup 1-1 | Cmac-lineage karyotype event |
| **`ancestral_split_both_retained`** | ≥3 species 1-2 with concordant target sets | Ancestral fission preserved across most lineages |
| **`recurrent_fission_hotspot`** | ≥3 species 1-2 with **discordant** target sets | Independent reuse of same fragile region — strongest mechanism |
| **`no_karyotype_change`** | All species 1-1 (≥3 species) | No karyotype event at this resolution; intra-chromosomal only |
| **`unresolved`** | Sparse / mixed / 1-3 / 1-4+ dominant | Insufficient data at this resolution |

### Hooks into age-model auto-suggest

`_msAutoSuggestAgeModel` now takes a 4th argument `karyoVerdict`. New
**RULE 0** (highest priority) consumes the verdict directly:

- `cgar_lineage_fission` → **LINEAGE-KARYO** at high confidence
- `cmac_lineage_fission` → **LINEAGE-KARYO** at high confidence
- `ancestral_split_both_retained` → **LINEAGE-KARYO** at high confidence
- `recurrent_fission_hotspot` → **MULTI-AGE-HOTSPOT** at high confidence

Both call sites (center column rendering + TSV export) are updated to
pass the verdict.

### Lineage table now has karyo class column

When **both** layers (`synteny_multispecies_v1` AND
`karyotype_lineage_v1`) are loaded, the lineage distribution table
gains an extra column showing each species' karyotype class
(`1-1` / `1-2` / `1-3` / `1-4+`) next to the boundary status.
Without the karyotype layer, the column is hidden — clean graceful
degradation.

### New JSON layer schema

```jsonc
{
  "tool": "karyotype_lineage_v1",
  "schema_version": 1,
  "params": {
    "method": "mashmap_one_to_one",
    "segment_bp": 1000000,
    "identity_pct": 85,
    "filter": "one-to-one",
    "source_pipeline": "/.../mashmap_fusion_scan_lines"
  },
  "outgroup_species": ["Tros"],
  "per_cgar_chr": [
    {
      "cgar_chr": "C_gar_LG28",
      "classes_by_species": {
        "Cgar":  { "class": "1-1", "targets": ["C_gar_LG28"] },
        "Cmac":  { "class": "1-2", "targets": ["C_mac_LG01", "C_mac_LG18"] },
        "Cfus":  { "class": "1-1", "targets": ["C_fus_LG07"] },
        "Tros":  { "class": "1-1", "targets": ["Tros_chr08"] }
      }
    }
  ]
}
```

Tool name accepts both `karyotype_lineage_v1` and
`mashmap_karyotype_lineage_v1` so pipeline scripts can use either.

### LANTA-side converter

`build_karyotype_lineage_v1_json.py` reads your existing
`mashmap_fusion_scan_lines/summaries/*.summary.tsv` files directly and
emits the atlas JSON. Round-trip tested: synthetic summary files →
converter → atlas detector returns true. After you re-run the
all-vs-all scan with `fClaHyb_Gar_LG.fa` as a query (currently it's
been run with A_melas etc.), invoke as:

```bash
python3 build_karyotype_lineage_v1_json.py \
  --summaries-dir /project/lt200308-agbsci/01-catfish_assembly/05_ancestral_karyotype/mashmap_fusion_scan_lines/summaries \
  --query-genome fClaHyb_Gar_LG.fa \
  --outgroup T_rosablanca.fa \
  --species-map species_map.tsv \
  --output karyotype_lineage_v1.json
```

A `species_map.tsv` template matching your real `genome_list.txt` is
included.

### Honest limitation, on-screen

The karyotype panel's legend says explicitly:

> *Resolution: 1 Mb segments at 85% identity (mashmap one-to-one).
> Class is chromosome-scale only — cannot resolve breakpoint coordinates.
> Pair with the breakpoint-scale wfmash synteny layer for kbp-resolution.*

This is the vocabulary discipline reviewer 2 will check.

---

## State + persistence

| Slot | Purpose | Persisted |
|---|---|---|
| `state.karyotypeLineage` | parsed `karyotype_lineage_v1.json` | localStorage `inversion_atlas.karyotypeLineage.v1` |

Wired into `_classifyJSONKind`, `loadMultipleJSONs` dispatch chain,
and the restore-on-load chain alongside other multi-species layers.

---

## Demo file

`karyotype_lineage_v1.demo.json` exercises all four "interesting"
verdicts:

- **`C_gar_LG28`**: Cmac 1-2, everyone else 1-1 → `cmac_lineage_fission`
- **`C_gar_LG07`**: 3 Clariidae species 1-2 with discordant targets → `recurrent_fission_hotspot`
- **`C_gar_LG10`**: Cgar 1-2, others 1-1 → `cgar_lineage_fission`
- **`C_gar_LG01`**: All 1-1 → `no_karyotype_change`

Drop into the atlas + a `cs_breakpoints_v1.json` covering those LGs;
each breakpoint will get a different karyotype verdict and the age
model will react accordingly.

---

## What's left for next turn

- **Run the actual Cgar-as-query mashmap scan on LANTA.** The current
  `mashmap_fusion_scan_lines/` data only has A_melas-as-query rows
  (and similar non-Cgar queries). The atlas needs Cgar-as-query
  summaries to populate `karyotype_lineage_v1.json` for real.
- **Targeted wfmash refinement around breakpoints.** Per your note,
  this is the breakpoint-scale wfmash run (small genome regions
  around each Cgar↔Cmac breakpoint, expanded to all comparative
  genomes). That feeds the existing `synteny_multispecies_v1` layer,
  not the new karyotype layer.
- **Page 16c (chromosome-karyotype atlas page).** A dedicated page
  showing all 28 Cgar LGs with their verdict chips on a single
  overview. Different scope from page 16b's per-breakpoint cockpit.
- **TSV export columns for karyo verdict** — `_msBuildClassificationTSV`
  doesn't yet have a `karyo_verdict` column. Half an hour of work.
- **busco_anchors_v1** layer for deep-divergence dotplot backbone.

---

## Files changed

```
Inversion_atlas.html                      +427 lines
test_turn124_karyotype_lineage.js         NEW (88 tests)
karyotype_lineage_v1.demo.json            NEW (demo)
build_karyotype_lineage_v1_json.py        NEW (LANTA converter)
species_map.tsv                           NEW (LANTA template)
CHANGELOG_turn124_karyotype_lineage.md    THIS FILE
```
