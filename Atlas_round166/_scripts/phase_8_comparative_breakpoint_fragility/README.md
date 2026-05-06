# phase_8_comparative_breakpoint_fragility

**Comparative repeat-density / breakpoint-fragility layer for the Inversion Atlas.**

Quentin Andres · MS_Inversions_North_african_catfish · LANTA `lt200308`
Codebase v8.5 · Phase 8 · v0.1 (initial scaffold)

---

## 1. What this module does

For every Gar inversion candidate, this module asks two distinct questions:

1. **Focal-species (Gar) question.** Is the breakpoint TE-enriched relative to chromosome and local background? Are *both* breakpoints enriched, or only one?
2. **Comparative question.** When the homologous region is located in another catfish assembly (Mac, Pangasius, Silurus, Ictalurus, …) using a synteny mapping, does that *homologous boundary* also sit in TE-rich, repeat-rich architecture? Is the boundary conserved across lineages? Is there orientation discordance? Chromosome-context change? Terminal context?

It produces a **comparative repeat-density evidence layer** — one JSON per candidate, plus chromosome-level rollups and a flat summary TSV — that plugs into the existing Atlas (`pca_scrubber_v3` / `repeat_density_v2` schema family) and supports downstream hypothesis testing.

The module does not call inversions, does not call SVs, does not run synteny — it is **strictly a repeat-density evidence layer** built on top of inputs that already exist.

---

## 2. What this module does NOT prove

This is the most important section. Read it before writing any text using the outputs.

- **TE density is not polymorphism.** A TE-rich region is fragile, not segregating.
- **A homologous boundary in Pangasius is not a Pangasius inversion.** Without population resequencing, heterozygous-assembly evidence, or read-level SV support, we cannot claim polymorphism in any non-resequenced species.
- **Conserved synteny edges are not shared inversions.** A boundary that exists in both Gar and Silurus tells us the *boundary* is old, not the *inversion*.
- **Inversion age ≠ breakpoint age.** A modern Gar inversion can reuse an ancient TE-rich boundary. The reverse is also possible.
- **Architecture class E (recurrent rearrangement hotspot) is a prediction, not a finding.** It says: "this region looks like the kind of region where inversions recur," not "this region is polymorphic in N species."

The output JSON is deliberately structured to make these distinctions hard to lose. Two fields per candidate are mandatory:

```json
"polymorphism_confirmed_in_species": ["C_gariepinus"],
"polymorphism_unknown_in_species": ["C_macrocephalus", "Pangasius_hypophthalmus", ...]
```

Anything written about a species in `polymorphism_unknown_in_species` must use hedged language (see §6).

---

## 3. Input requirements

Inputs live in `config/`. None are produced by this module — they are produced upstream and registered here.

### 3.1 Gar inversion candidates — `config/candidate_breakpoints.tsv`

Tab-delimited, header required. Columns:

| column                | required | type   | notes                                                  |
|-----------------------|----------|--------|--------------------------------------------------------|
| `candidate_id`        | yes      | str    | unique, e.g. `LG28_INV_001`                            |
| `chrom`               | yes      | str    | Gar chromosome, e.g. `C_gar_LG28`                      |
| `start`               | yes      | int    | 1-based inclusive interval start (Gar reference)       |
| `end`                 | yes      | int    | 1-based inclusive interval end                         |
| `left_breakpoint`     | yes      | int    | 1-based left breakpoint coordinate                     |
| `right_breakpoint`    | yes      | int    | 1-based right breakpoint coordinate                    |
| `confidence_tier`     | no       | str    | `tier1` / `tier2` / `tier3` / free text                |
| `evidence_layers`     | no       | str    | comma-list, e.g. `local_pca,mds,sv,fst`                |
| `arrangement_label`   | no       | str    | optional, e.g. `REF/HET/INV` band labels               |
| `gar_polymorphic`     | no       | bool   | default `true`; set `false` for non-polymorphic tests  |

Coordinates are converted to 0-based half-open internally for BED/bedtools/pyranges.

### 3.2 Genome FASTA index per species — registered in `config/species_manifest.tsv`

| column          | required | notes                                                      |
|-----------------|----------|------------------------------------------------------------|
| `species_id`    | yes      | short canonical id, e.g. `C_gariepinus`, `C_macrocephalus` |
| `assembly_name` | yes      | e.g. `fClaHyb_Gar_LG`                                      |
| `fai_path`      | yes      | absolute path to `genome.fa.fai`                           |
| `is_focal`      | yes      | `true` only for Gar                                        |
| `notes`         | no       | free text                                                  |

The `.fai` is read with the standard layout: `chrom \t length \t offset \t linebases \t linewidth`.

### 3.3 TE / repeat annotation files — `config/te_file_manifest.tsv`

Either built **automatically** by `00_manifest/build_te_manifest.py` over a messy root directory, or curated by hand. Columns:

| column             | required | notes                                                          |
|--------------------|----------|----------------------------------------------------------------|
| `species_id`       | yes      | must match `species_manifest.tsv`                              |
| `assembly_name`    | yes      | matches `species_manifest.tsv`                                 |
| `file_path`        | yes      | absolute path                                                  |
| `file_format`      | yes      | one of: `edta_intact`, `edta_anno`, `gff3`, `repeatmasker_out`, `bed`, `unknown` |
| `priority`         | no       | int; lower = preferred when multiple files cover same assembly |
| `notes`            | no       | free text                                                      |

If `file_format = unknown`, the file is logged and skipped at normalization. Manifest is the single source of truth — anything not listed is invisible to the pipeline.

### 3.4 Cross-species synteny mapping (optional) — `config/synteny_mapping.tsv`

If absent, the pipeline still runs the **focal-only** layer (Gar TE density at Gar breakpoints). When present:

| column              | required | notes                                                       |
|---------------------|----------|-------------------------------------------------------------|
| `candidate_id`      | yes      | FK to `candidate_breakpoints.tsv`                           |
| `focal_species`     | yes      | typically `C_gariepinus`                                    |
| `focal_chrom`       | yes      |                                                             |
| `focal_start`       | yes      | 1-based                                                     |
| `focal_end`         | yes      | 1-based                                                     |
| `target_species`    | yes      | must match `species_manifest.tsv`                           |
| `target_chrom`      | yes      |                                                             |
| `target_start`      | yes      | 1-based                                                     |
| `target_end`        | yes      | 1-based                                                     |
| `orientation`       | yes      | `+` or `-` relative to focal                                |
| `synteny_block_id`  | no       |                                                             |
| `boundary_type`     | no       | `internal` / `terminal` / `block_edge` / `chr_context_change` |
| `mapping_method`    | no       | e.g. `wfmash`, `minimap2`, `nucmer`, `liftover`             |

Multiple rows per `(candidate_id, target_species)` are allowed (e.g. fragmented synteny across multiple target scaffolds).

---

## 4. Output JSON schema (per candidate)

One file per candidate at `output/per_candidate_json/<candidate_id>.json`. Schema mirrors the `repeat_density_v2`-family conventions used by the Atlas.

```json
{
  "schema_version": "comparative_breakpoint_fragility_v0.1",
  "candidate_id": "LG28_INV_001",
  "focal_species": "C_gariepinus",
  "focal_assembly": "fClaHyb_Gar_LG",
  "focal_interval": {
    "chrom": "C_gar_LG28",
    "start": 123456,
    "end": 456789,
    "left_breakpoint": 123456,
    "right_breakpoint": 456789,
    "confidence_tier": "tier1",
    "evidence_layers": ["local_pca", "mds", "sv", "fst"]
  },
  "gar_population_status": {
    "polymorphic_in_gar": true,
    "evidence_layers": ["local_pca", "mds", "sv", "fst"],
    "confidence": "high"
  },
  "polymorphism_confirmed_in_species": ["C_gariepinus"],
  "polymorphism_unknown_in_species": ["C_macrocephalus", "Pangasius_hypophthalmus"],

  "repeat_density": {
    "left_breakpoint":  { "window_50kb":{...}, "window_100kb":{...}, "window_500kb":{...} },
    "right_breakpoint": { "window_50kb":{...}, "window_100kb":{...}, "window_500kb":{...} },
    "inside_interval": {...},
    "left_flank_outside":  {...},
    "right_flank_outside": {...},
    "chromosome_background": {...},
    "local_background_2mb": {...}
  },

  "comparative_species": [
    {
      "species_id": "C_macrocephalus",
      "assembly_name": "Mac_assembly_v1",
      "homologous_region": {
        "chrom": "Mac_LG01",
        "start": 111111,
        "end": 999999,
        "orientation_relative_to_focal": "-",
        "synteny_block_id": "block_abc",
        "mapping_method": "wfmash"
      },
      "boundary_status": {
        "left_boundary_present": true,
        "right_boundary_present": true,
        "orientation_state": "alternative",
        "chromosome_context": "different_chromosome_context",
        "terminal_or_internal": "internal"
      },
      "repeat_density": {
        "left_boundary":  { "window_50kb":{...}, "window_100kb":{...}, "window_500kb":{...} },
        "right_boundary": { "window_50kb":{...}, "window_100kb":{...}, "window_500kb":{...} },
        "inside_interval": {...},
        "chromosome_background": {...},
        "local_background_2mb": {...}
      },
      "interpretation": "TE-rich homologous breakpoint in C_macrocephalus assembly; supports candidate recurrent hotspot, does not demonstrate polymorphism."
    }
  ],

  "classification": {
    "architecture_class": "E",
    "architecture_label": "recurrent_rearrangement_hotspot",
    "age_model": "OLD_BP_YOUNG_INV",
    "prediction_label": "predicted_polymorphic_hotspot",
    "confidence": "medium",
    "modifiers": ["conserved_boundary", "orientation_discordance"]
  },

  "safe_conclusion": "Confirmed polymorphic inversion in C. gariepinus; homologous TE-rich breakpoint architecture in comparative assemblies suggests a recurrent rearrangement hotspot but does not demonstrate polymorphism in non-resequenced species."
}
```

Each window block has the same shape:

```json
{
  "window_bp": 100000,
  "callable_bp": 99812,
  "te_bp_overlap": 32104,
  "te_density": 0.3217,
  "te_density_by_class": {
    "DNA": 0.18, "LTR": 0.07, "LINE": 0.02, "SINE": 0.01,
    "Helitron": 0.02, "TIR": 0.01, "MITE": 0.00,
    "satellite": 0.00, "low_complexity": 0.00, "unknown": 0.01
  },
  "fold_vs_chromosome": 1.79,
  "fold_vs_local_2mb": 1.53,
  "percentile_chr": 0.91,
  "breakpoint_TE_hotspot_score": 1.66
}
```

`breakpoint_TE_hotspot_score` is the simple first-pass score:
`mean(fold_vs_chromosome, fold_vs_local_2mb)` weighted by `percentile_chr`.

---

## 5. Classification labels

### 5.1 Architecture class

| code | label                                  | meaning                                                                                                 |
|------|----------------------------------------|---------------------------------------------------------------------------------------------------------|
| A    | simple_population_inversion            | clean Gar polymorphic inversion; comparative architecture not strongly indicative                       |
| B    | synteny_boundary_inversion             | breakpoints overlap conserved synteny-block edges                                                       |
| C    | fusion_fission_associated              | breakpoints near lineage-specific fusion/fission scars                                                  |
| D    | terminal_translocation_boundary        | breakpoint near chromosome/scaffold tip in ≥1 species                                                   |
| E    | recurrent_rearrangement_hotspot        | TE-rich architecture in multiple species + boundary modifiers; the strongest "predicted hotspot" call   |
| F    | ambiguous_assembly_risk                | breakpoint overlaps gaps / N-rich / collapsed-repeat / low-mappability — interpret with caution         |
| G    | ancient_conserved_inversion_candidate  | shared orientation switch across deep lineages; possible OLD_POLY                                       |

A candidate can have **one** primary class. Modifiers (below) refine it.

### 5.2 Age / model tag

| tag                              | meaning                                                                              |
|----------------------------------|--------------------------------------------------------------------------------------|
| `YOUNG_POP`                      | confirmed recent / population-level inversion in Gar                                 |
| `OLD_POLY`                       | inversion itself may be old (deep arrangement divergence / shared orientation switch)|
| `OLD_BP_YOUNG_INV`               | old breakpoint substrate reused by younger Gar inversion                             |
| `LINEAGE_KARYO`                  | species-lineage karyotype breakpoint                                                 |
| `MULTI_AGE_HOTSPOT`              | several different-age rearrangement signals reuse the same region                    |
| `ANCIENT_BOUNDARY`               | conserved synteny / block boundary, but inversion not shown in old species           |
| `PREDICTED_POLYMORPHIC_HOTSPOT`  | comparative architecture predicts a fragile region; polymorphism unknown in non-Gar  |

### 5.3 Modifiers (orthogonal flags)

`conserved_boundary` · `orientation_discordance` · `chromosome_context_change` · `terminal_context` · `assembly_risk`

### 5.4 Hotspot tier (derived from the score)

| tier      | rule                                                                                          |
|-----------|-----------------------------------------------------------------------------------------------|
| `high`    | both breakpoints enriched ≥2× chromosome background **or** ≥90th percentile chromosome-wide   |
| `medium`  | one breakpoint enriched, or both mildly enriched (1.5× – 2× background)                       |
| `low`     | no enrichment                                                                                  |

---

## 6. Safe interpretation language

Use these phrasings. Don't drift.

**Confirmed (only with population data, currently Gar):**
- "polymorphic inversion in *C. gariepinus* broodstock"
- "segregating in the resequenced cohort"

**Allowed comparative inference (assembly + TE evidence only):**
- "homologous breakpoint is TE-rich"
- "breakpoint boundary appears conserved across catfish assemblies"
- "candidate recurrent rearrangement hotspot"
- "possible ancient breakpoint substrate"
- "predicted polymorphic hotspot"
- "the homologous *Pangasius* region carries TE-rich architecture consistent with…"

**Forbidden:**
- "*Pangasius* is polymorphic for this inversion." ❌
- "shared inversion across catfish." ❌  (unless population SV evidence in each species)
- "*Silurus* carries this inversion in the heterokaryotypic state." ❌

**Manuscript paragraph (default boilerplate):**

> Although structural polymorphism could only be directly tested in *C. gariepinus*, we used TE density as a first-pass proxy for breakpoint fragility across homologous regions in comparative catfish assemblies. Several Gar polymorphic inversion boundaries overlapped TE-enriched synteny-block edges that were conserved or rearranged across deeper catfish lineages, suggesting that present-day structural haplotypes can reuse ancient, repeat-rich chromosome architecture.

**Short version:**

> We confirm polymorphism in Gar, and we predict which homologous breakpoint regions in other catfish are likely unstable based on TE-rich conserved synteny boundaries.

---

## 7. How this fits the Nature Communications five-Result architecture

| layer          | role                                                                                                |
|----------------|-----------------------------------------------------------------------------------------------------|
| **applied**    | Result 5 (diversity / breeding integration) — fragile breakpoint atlas as a target list for marker design and broodstock-decision support |
| **mechanistic**| Result 4 (validation / evolution) — TE-rich substrate explains *why* these inversions exist where they do |
| **evolutionary**| Result 4 (validation / evolution) — ancient chromosome architecture seeds modern structural haplotypes; *C. gariepinus* inversions sit on top of catfish-wide fragile boundaries |

The module does not change Result 1–3 (genome coordinate system; hatchery population context; inversion discovery atlas). It feeds Result 4 directly and Result 5 as a candidate-prioritisation layer.

---

## 8. Pipeline DAG

```
00_manifest/build_te_manifest.py
        │
        ▼
config/te_file_manifest.tsv ──┐
                              │
01_normalize_te/normalize_te_annotations.py
        │
        ▼
output/normalized_te/<species>.te.bed.gz  (sorted, indexed)
                              │
config/candidate_breakpoints.tsv ──┐
config/species_manifest.tsv (.fai) ─┤
                              │     │
                              ▼     ▼
02_breakpoint_windows/make_breakpoint_windows.py
        │
        ▼
output/breakpoint_windows/<species>.windows.bed
                              │
        ┌─────────────────────┴─────────────────────┐
        ▼                                           ▼
03_density/compute_te_density.py    04_comparative/compute_comparative_te_density.py
        │                                           │  (consumes synteny_mapping.tsv via
        │                                           │   04_comparative/map_homologous_breakpoints.py)
        ▼                                           ▼
output/density/focal/<candidate>.json  output/density/comparative/<candidate>__<target_species>.json
                              │                     │
                              └─────────┬───────────┘
                                        ▼
                       05_classification/classify_breakpoint_fragility.py
                                        │
                                        ▼
                       06_json_export/export_breakpoint_fragility_json.py
                                        │
                                        ▼
        output/per_candidate_json/<candidate_id>.json
        output/per_chromosome_json/<chrom>.json
        output/summary_tables/breakpoint_fragility_summary.tsv
                                        │
                                        ▼
                          07_plots/plot_breakpoint_te_density.R
```

Logs land in `output/logs/<step>__<timestamp>.log`. Every step writes a sidecar `.run.json` with input hashes, parameters, and counts.

---

## 9. Hypothesis tests supported by the output

The summary TSV (`output/summary_tables/breakpoint_fragility_summary.tsv`) is laid out so each of these tests is a one-liner in R:

| H | Hypothesis                                                                              | Test                                                             |
|---|------------------------------------------------------------------------------------------|------------------------------------------------------------------|
| 1 | Gar breakpoints are TE-enriched vs. chromosome background                                | `left_TE_density_100kb` & `right_TE_density_100kb` vs. chrom bg  |
| 2 | High-confidence Gar candidates are more TE-enriched than low-confidence                  | `confidence_tier` × `breakpoint_TE_hotspot_score`                |
| 3 | Cross-species chromosome-context change ⇒ higher Gar breakpoint TE density                | `mac_chromosome_context == "different"` × `*_fold_vs_chr`        |
| 4 | Some Gar inversions reuse ancient synteny boundaries                                     | `num_species_with_boundary` ≥ 2 × architecture class B/G          |
| 5 | Some homologous old-species regions are predicted polymorphic hotspots                   | `num_species_TE_enriched` ≥ 2 + `prediction_label == "PPH"`       |

---

## 10. Reproducibility / HPC notes

- All scripts are pure Python ≥ 3.9 + R ≥ 4.1, no exotic dependencies. `pyranges` preferred; falls back to `bedtools` shell calls if `pyranges` is missing.
- One LANTA SLURM submission script per step in `submit/` (not yet shipped — to be filled in as inputs stabilise).
- Working directory on LANTA: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_8_comparative_breakpoint_fragility/`.
- File manifests + run manifests are JSON; deletion-safe; nothing in `output/` is required for re-running upstream.

---

## 11. Roadmap (post-v0.1)

- **v0.2**: callable-mask integration (gaps, low-mappability) — currently `callable_bp ≈ window_bp`.
- **v0.3**: per-class hotspot tests (LTR-only, TIR-only) — schema already supports it.
- **v0.4**: bootstrap empirical p-values for percentile_chr (matched-GC / matched-mappability windows).
- **v0.5**: TE age (Kimura div) layer, restricted to candidates flagged hotspot tier ≥ medium.
- **v0.6**: integration with `pca_scrubber_v3` candidate page (compact track block under FIG_C08).

---

## 12. License / attribution

Internal to MS_Inversions_North_african_catfish. PI: Prapansak Srisapoome.
Module author: Quentin Andres. Pipeline scaffolding generated with assistance from Claude.
