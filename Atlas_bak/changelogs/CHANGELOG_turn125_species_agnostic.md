# CHANGELOG — turn 125: species-agnostic karyotype layer

**Atlas Δ:** +106 lines (61,570 → 61,676)
**Tests:** 50 new (619/619 cumulative green across 10 test suites)
**No regressions** — all prior turns 117–124 still pass.

---

## What shipped

The karyotype layer (`karyotype_lineage_v1`) is now **species-agnostic**.
The atlas no longer hardcodes `Cgar` and `Cmac` — instead, the JSON declares
`focal_species`, `sister_species`, and `outgroup_species`, and the
polarization rules read those roles from the JSON itself.

This means: the same atlas page 16b code will work for the **Cmac
manuscript** (when you pivot the focal to Cmac) without code changes.
Just point the converter at fClaHyb_Mac_LG.fa as the query genome,
declare Cgar as sister, and the polarization re-derives correctly.

### Schema changes

**v2 schema (preferred)** uses the species-agnostic field names:

```jsonc
{
  "tool": "karyotype_lineage_v1",
  "schema_version": 2,
  "params": { ... },
  "focal_species": "Cgar",
  "sister_species": ["Cmac"],
  "outgroup_species": ["Tros"],
  "per_focal_chr": [
    {
      "focal_chr": "C_gar_LG28",
      "classes_by_species": {
        "Cgar": { "class": "1-1", "targets": ["C_gar_LG28"] },
        "Cmac": { "class": "1-2", "targets": ["C_mac_LG01", "C_mac_LG18"] },
        ...
      }
    }
  ]
}
```

**v1 schema is still accepted** as backward-compat. `per_cgar_chr` is
auto-normalized to `per_focal_chr`, and `focal_species` is inferred
from data (the species id whose 1-1 target equals the focal_chr).

### Verdict naming generalized

| Old verdict (still accepted as alias) | New verdict (species-agnostic) |
|---|---|
| `cgar_lineage_fission` | `focal_lineage_fission` |
| `cmac_lineage_fission` | `sister_lineage_fission` |
| `ancestral_split_both_retained` | unchanged |
| `recurrent_fission_hotspot` | unchanged |
| `no_karyotype_change` | unchanged |
| `unresolved` | unchanged |

The age-model auto-suggest **RULE 0** accepts all 6 keys (new + legacy
aliases), so old demo JSONs keep working without modification.

### Verdict label substitution

When the panel renders, the label substitutes the actual species name
into the verdict. The user sees:

- "Cgar lineage fission" (not literal "focal lineage fission")
- "Cmac lineage fission" (not literal "sister lineage fission")

If `sister_species` lists multiple species, the label uses whichever
sister actually carries the 1-2 — e.g. "Cmac/Cfus lineage fission".

### Signals object refactored

Old (hardcoded):
```js
signals.cgar_class  // string
signals.cmac_class  // string
```

New (species-agnostic):
```js
signals.focal_species          // string from JSON
signals.focal_class            // string
signals.sister_species         // array from JSON
signals.sister_classes         // dict species_id → class
signals.sister_one_to_one      // count
signals.sister_one_to_two      // count
signals.outgroup_classes       // dict species_id → class (already existed)
signals.outgroup_one_to_one    // count
signals.outgroup_one_to_two    // count
```

### Multi-outgroup support

The converter now accepts comma-separated outgroups:
```bash
--outgroup T_rosablanca.fa,N_graeffei.fa
```

The atlas was already storing outgroups as an array — this just makes
the LANTA-side script align with that.

### Sister-group support in the converter

New `--sister` flag for sister-species declaration:
```bash
--sister fClaHyb_Mac_LG.fa
```

Drives the `sister_lineage_fission` verdict.

### Honest limitation, on-screen (unchanged)

The legend still says:

> *Resolution: 1 Mb segments at 85% identity (mashmap one-to-one).
> Class is chromosome-scale only — cannot resolve breakpoint coordinates.
> Pair with the breakpoint-scale wfmash synteny layer for kbp-resolution.*

## Files changed

```
Inversion_atlas.html                          +106 lines
test_turn125_species_agnostic.js              NEW (50 tests)
test_turn124_karyotype_lineage.js             updated for new schema
build_karyotype_lineage_v1_json.py            v2 emitter, --sister flag, multi-outgroup
karyotype_lineage_v1.demo.json                v2 schema with focal/sister/outgroup
CHANGELOG_turn125_species_agnostic.md         THIS FILE
```

## What's left

This is the species-agnostic refactor. **Next turn (126)** is the
**wfmash refinement** layer — every interesting mashmap cell (1-2,
1-3) deserves a confirmation pass at higher resolution. Schema
addition: each `classes_by_species[sp]` entry gains an optional
`refined_by_wfmash: true|false|null` flag plus `refined_targets` /
`refined_coords` fields. Atlas displays the refinement state alongside
the mashmap class.

After that, **turn 127 cleanup** + then **BUSCO integration**.
