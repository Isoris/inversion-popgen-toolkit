# SPEC — Comparative TE-density / breakpoint fragility layer

**Status**: forward-looking spec. None implemented. Created end of turn 113-tail.

**Reading order**: this spec → `SPEC_OVERVIEW_multispecies_architecture.md` (sees how it fits) → existing boundaries-page TE density code in atlas (`_renderRepeatDensityPanel` line ~15770) for the Gar implementation this generalizes.

---

## Purpose

The atlas already does **per-Gar-chromosome TE density** on the boundaries page. This spec extends that capability to **other catfish species' assemblies**, where we have:
- genome assemblies + chromosome lengths
- TE / repeat annotations (EDTA, GFF3, RepeatMasker — messy and scattered)
- NO population resequencing for most species

The output is a JSON layer the atlas consumes, exactly the same way it consumes Gar's TE density today — just with multiple species available.

## Why this matters (the motivation, not buried in the technical detail)

Many catfish species we care about (Pangasius, Silurus, Ictalurus, Heterobranchus, Trichomycterus) live in places it would take years and substantial money to resequence properly — Vietnam, Peru, the Amazon basin, central Africa. Population genomics for those species is not happening on this manuscript's timeline, possibly never on this PI's career timeline.

But assemblies + TE annotations DO exist for all of them already. So instead of "we don't know," we can say "we don't know about polymorphism, but we know whether the homologous region is repeat-rich and structurally fragile-looking." That's a real comparative inference.

## What this layer DOES claim

- **Confirmed polymorphic inversion in Gar** — we have 226-sample resequencing
- **Homologous region in species X is TE-rich** — we have the assembly
- **The breakpoint substrate is conserved or rearranged across catfish lineages** — we have multiple assemblies + synteny
- **The homologous region is a candidate recurrent rearrangement hotspot** — TE enrichment is a fragility proxy

## What this layer does NOT claim

- "Pangasius is polymorphic for this inversion" — unsupported without resequencing
- "This inversion exists in Silurus" — unsupported without per-individual data
- TE density alone proves anything beyond "this region has a lot of repeats"

The polymorphism claim requires per-individual data. TE density supports a **fragility hypothesis**, which is weaker but legitimate.

## Vocabulary discipline

| Confirmed | Allowed comparative inference | NEVER |
|---|---|---|
| polymorphic inversion in Gar broodstock | homologous breakpoint is TE-rich | Pangasius is polymorphic |
| 226-sample carrier evidence | breakpoint boundary conserved across catfish assemblies | this inversion exists in species X |
| fixed in Gar lineage | candidate recurrent rearrangement hotspot | structural polymorphism in non-resequenced species |
| | possible ancient breakpoint substrate | |
| | predicted polymorphic hotspot | |

Use "may", "suggests", "candidate", "predicted" for non-resequenced species. Reviewer 2 will catch overclaiming if it's there.

---

## Module structure

Suggested layout (illustrative — Quentin may reorganize):

```
phase_8_comparative_breakpoint_fragility/
  README.md
  00_manifest/                 build_te_manifest.py
  01_normalize_te/             normalize_te_annotations.py
  02_breakpoint_windows/       make_breakpoint_windows.py
  03_density/                  compute_te_density.py
  04_comparative/              map_homologous_breakpoints.py
                               compute_comparative_te_density.py
  05_classification/           classify_breakpoint_fragility.py
  06_json_export/              export_breakpoint_fragility_json.py
  07_plots/                    plot_breakpoint_te_density.R
  config/                      species_manifest.tsv
                               te_file_manifest.tsv
                               candidate_breakpoints.tsv
                               synteny_mapping.tsv
  output/                      per_candidate_json/
                               per_chromosome_json/
                               summary_tables/
                               logs/
```

The pipeline is **modular intentionally**. Each step writes its outputs to disk, so steps can be re-run independently when bugs are found or new species added.

---

## Inputs

### 1. Gar inversion candidate table

Already exists in atlas state (`state.candidates`). Pipeline reads either an exported TSV or the atlas's JSON dump. Required columns:
- `candidate_id`, `chrom`, `start_bp`, `end_bp`, `left_breakpoint`, `right_breakpoint`
- Optional: `confidence_tier`, evidence layer flags, arrangement labels

### 2. Per-species genome FASTA index / scaffold sizes

`samtools faidx` output (`.fai`) per species. Provides `chrom`, `length`, `offset`. The pipeline only needs `chrom` and `length` columns.

### 3. TE / repeat annotation files

The messy part. EDTA produces `.TEanno.gff3` and `.intact.gff3`. RepeatMasker produces `.out`. Custom libraries produce variable-format GFF3. The manifest builder (`build_te_manifest.py`) does recursive search + format detection + species inference.

### 4. Cross-species synteny / homology table (optional)

If exists: maps Gar candidate intervals to homologous coords in other species. Schema:
```
candidate_id, focal_species, focal_chrom, focal_start, focal_end,
target_species, target_chrom, target_start, target_end,
orientation, synteny_block_id, boundary_type
```

If absent: pipeline still runs Gar-only TE density. Cross-species enrichment fields are populated later when synteny lands.

---

## Pipeline steps

### Step 0 — Manifest builder (`build_te_manifest.py`)

- Recursively walk a root TE-results directory
- Identify files by extension + content sniff:
  - `*.gff3`, `*.gff` (with EDTA-style attributes)
  - `*.bed`
  - `*.out` (RepeatMasker)
  - `EDTA.TEanno.gff3`, `EDTA.intact.gff3`
- Infer species from path components (regex against known species names)
- Allow manual override via `species_manifest.tsv` for ambiguous cases
- Output `te_file_manifest.tsv`: `species`, `file_path`, `format`, `inferred_or_manual`, `notes`
- Log skipped files + unrecognized formats to `logs/manifest_unrecognized.tsv`

**Tolerance principle**: never crash on a malformed file. Log it, skip it, continue.

### Step 1 — TE normalizer (`normalize_te_annotations.py`)

- Convert each format to standard BED:
  ```
  chrom    start0    end    repeat_id    score    strand    repeat_class    repeat_family    source_file    species
  ```
- 1-based GFF3 → 0-based BED conversion done correctly (off-by-one is the most common bug here)
- Handle EDTA's `Classification=` tag for class/family
- Handle RepeatMasker's column 11 / 10 mapping to class/family
- Drop malformed lines but log them per-file
- Write per-species: `01_normalized/<species>.te.bed.gz`

### Step 2 — Breakpoint window builder (`make_breakpoint_windows.py`)

- Read candidate breakpoints
- Build window BEDs for each candidate × each window radius
- Window radii: **don't be religious about specific kb values**. Use a configurable list (e.g. 50, 100, 500 kb) but the framework tolerates any list; downstream consumers parameterize. The numbers are illustrative until empirically tuned on real data.
- Output: `02_windows/<candidate_id>_breakpoint_windows.bed`

### Step 3 — Density calculator (`compute_te_density.py`)

For each candidate × each window × each species:
- Compute `te_bp_overlap` via `bedtools intersect` or `pyranges`
- Compute `te_density = te_bp / window_size`
- Compute by repeat class if class info is available (DNA, LTR, LINE, SINE, Helitron, TIR, MITE, satellite, simple, low-complexity, unknown)
- Compute chromosome background (full-chromosome density) and local background (e.g. ±2 Mb excluding the breakpoint windows)
- Compute fold-enrichment vs both backgrounds
- Compute percentile rank within chromosome (sliding-window distribution)
- Output: `03_density/<candidate_id>_te_density.json`

### Step 4 — Comparative mapping (`map_homologous_breakpoints.py` + `compute_comparative_te_density.py`)

If synteny mapping exists:
- For each candidate, find homologous interval in each target species
- Run Step 3 logic on the homologous region using that species' TE BED
- Capture orientation + synteny boundary status
- Output per-species comparative density into the per-candidate JSON

### Step 5 — Classification (`classify_breakpoint_fragility.py`)

Assign architecture class A–G + age model + hotspot score per candidate. Logic detailed in next section.

### Step 6 — JSON export (`export_breakpoint_fragility_json.py`)

- One JSON per candidate (atlas-compatible)
- One combined JSON per chromosome
- One summary JSON across all candidates
- One TSV summary for plotting/statistics
- Atlas reads these as a new layer (key: `comparative_te_breakpoint_fragility`)

### Step 7 — Plots (`plot_breakpoint_te_density.R`)

Optional. Standalone diagnostic plots for QC; atlas page can be the primary visualization once layer is wired.

---

## Classification

### Architecture classes (Layer 1, extends the existing A–F schema with G)

| Class | Label | When |
|---|---|---|
| A | Simple inversion | Two clean inversion edges, no synteny boundary, no other rearrangement signal |
| B | Synteny-boundary inversion | Inversion edge sits at edge of conserved syntenic block |
| C | Fusion/fission-associated | One breakpoint overlaps cross-species chromosome split/join |
| D | Terminal translocation / reversed attachment | Region near chromosome end + homologous segment on different chromosome in another species, possibly reversed |
| E | Recurrent rearrangement hotspot | Same region carries multiple rearrangement signals + repeat-rich + multiple species show different rearrangements |
| F | Ambiguous / assembly-risk | Gap, low mappability, collapsed repeat, scaffold join, weak SV support, only one breakpoint visible |
| **G** | **Ancient conserved inversion** | **Conserved orientation switch across multiple species; inversion itself may be old** |

Class G is new in this spec and adds to the schema in the existing
"Breakpoint classification framework" on atlas page 5. Update the help-page
table to include it when this spec is implemented.

### Age model tags (Layer 2, extended with two new tags)

| Tag | Meaning |
|---|---|
| YOUNG_POP | Confirmed recent / population-level inversion in Gar |
| OLD_POLY | Inversion itself may be old (deep arrangement divergence or shared orientation switch) |
| OLD_BP_YOUNG_INV | Old breakpoint substrate reused by younger Gar inversion |
| LINEAGE_KARYO | Species-lineage karyotype breakpoint |
| MULTI_AGE_HOTSPOT | Several different-age rearrangement signals reuse same region |
| **ANCIENT_BOUNDARY** | **Conserved synteny / block boundary, but inversion not shown in old species** |
| **PREDICTED_POLYMORPHIC_HOTSPOT** | **Candidate region where other species may harbor undetected polymorphism — but not confirmed** |

Two new tags: ANCIENT_BOUNDARY (for cases where the boundary is shared but no inversion signal in older species) and PREDICTED_POLYMORPHIC_HOTSPOT (for explicit recording of "we predict, we did not test"). The PREDICTED tag is critical for vocabulary discipline — it's the wording reviewer 2 will look for.

### Hotspot score (illustrative thresholds)

For each breakpoint window:
- High-confidence: both breakpoints TE-enriched (>2× chromosome background or >90th percentile chromosome-wide)
- Medium-confidence: one breakpoint TE-enriched, OR both 1.5–2× background
- Low-confidence: no TE enrichment

Plus modifiers:
- `+conserved_boundary` if homologous block edge exists in older species
- `+orientation_discordance` if homologous block orientation differs across species
- `+chromosome_context_change` if homologous block on different chromosome / scaffold across species
- `+terminal_context` if breakpoint near chromosome / scaffold tip
- `+assembly_risk` if breakpoint overlaps gap / N-rich / low mappability

The thresholds (2×, 1.5×, 90th percentile) are starting points. Tune on real data.

---

## Output JSON schema

Per-candidate JSON. Schema is the contract; specific window sizes are illustrative:

```jsonc
{
  "candidate_id": "LG28_INV_001",
  "focal_species": "C_gariepinus",
  "focal_interval": {
    "chrom": "Gar_LG28",
    "start": 15115000,
    "end": 18005000,
    "left_breakpoint": 15115000,
    "right_breakpoint": 18005000,
    "orientation": "inverted_in_some_Gar_haplotypes"
  },
  "gar_population_status": {
    "polymorphic_in_gar": true,
    "evidence_layers": ["local_pca", "mds", "sv", "fst", "dxy"],
    "confidence": "high"
  },
  "repeat_density": {
    "left_breakpoint": {
      "window_50kb":  { "te_bp": 32000, "window_bp": 100000, "te_density": 0.32,
                         "fold_vs_chr": 1.8, "fold_vs_local": 1.5, "percentile_chr": 0.91 },
      "window_100kb": { "...": "..." },
      "window_500kb": { "...": "..." }
    },
    "right_breakpoint": { "...": "..." },
    "inside_interval":  { "te_density": 0.24, "fold_vs_chr": 1.2 },
    "chromosome_background": { "te_density": 0.18 },
    "local_background":      { "te_density": 0.21 }
  },
  "comparative_species": [
    {
      "species": "C_macrocephalus",
      "assembly": "fClaHyb_Mac_LG",
      "homologous_region": {
        "chrom": "C_mac_LG01", "start": 8101200, "end": 11332100,
        "orientation_relative_to_gar": "-",
        "synteny_block_id": "blk_0001"
      },
      "boundary_status": {
        "left_boundary_present": true,
        "right_boundary_present": true,
        "orientation_state": "alternative",
        "chromosome_context": "different_chromosome_context",
        "terminal_or_internal": "internal"
      },
      "repeat_density": {
        "left_boundary_100kb":  { "te_density": 0.41, "fold_vs_chr": 2.3, "percentile_chr": 0.96 },
        "right_boundary_100kb": { "te_density": 0.35, "fold_vs_chr": 1.9, "percentile_chr": 0.90 }
      },
      "polymorphism_confirmed_in_species": false,
      "polymorphism_unknown_in_species": true,
      "interpretation": "TE-rich homologous breakpoint; comparative architecture supports recurrent hotspot, not confirmed polymorphism"
    },
    {
      "species": "Pangasius_hypophthalmus",
      "assembly": "Phyp_assembly_v1",
      "homologous_region": { "...": "..." },
      "boundary_status": { "...": "..." },
      "repeat_density": { "...": "..." },
      "polymorphism_confirmed_in_species": false,
      "polymorphism_unknown_in_species": true,
      "interpretation": "candidate ancient boundary; polymorphism unknown"
    }
  ],
  "classification": {
    "architecture_class": "E",
    "architecture_label": "recurrent_rearrangement_hotspot",
    "age_model": "OLD_BP_YOUNG_INV",
    "prediction_label": "predicted_polymorphic_hotspot",
    "confidence": "medium"
  },
  "safe_conclusion": "Confirmed polymorphic inversion in Gar; homologous TE-rich breakpoint architecture in comparative assemblies suggests a recurrent rearrangement hotspot but does not demonstrate polymorphism in non-resequenced species."
}
```

The `polymorphism_confirmed_in_species` and `polymorphism_unknown_in_species` fields are explicit anti-overclaiming markers. They're in the JSON for every species so the atlas can color-code accordingly and so manuscript text generators can't accidentally claim otherwise.

---

## TSV summary columns

For quick R/Python plotting + reviewer responses:

```
candidate_id chrom start end left_bp right_bp gar_polymorphic
left_TE_density_100kb right_TE_density_100kb
left_fold_vs_chr right_fold_vs_chr
left_percentile_chr right_percentile_chr
both_breakpoints_TE_enriched
inside_TE_density chrom_TE_density local_TE_density
mac_boundary_present mac_orientation
mac_left_TE_density_100kb mac_right_TE_density_100kb
pangasius_boundary_present pangasius_orientation
pangasius_left_TE_density_100kb pangasius_right_TE_density_100kb
num_species_with_boundary num_species_TE_enriched
architecture_class age_model prediction_label confidence
safe_conclusion
```

---

## Hypothesis tests this layer enables

H1 — **Gar inversion breakpoints are TE-enriched vs chromosome background.**
Test: breakpoint ±100 kb TE density vs chromosome-wide sliding-window distribution. One-sided test. Bonferroni or BH correction across candidates.

H2 — **High-confidence inversion candidates are more TE-enriched than low-confidence ones.**
Test: candidate confidence tier vs breakpoint TE percentile. Rank correlation.

H3 — **Candidates with cross-species chromosome-context changes have higher breakpoint TE density.**
Test: context_change yes/no vs breakpoint TE enrichment. Two-group comparison.

H4 — **Some Gar polymorphic inversions reuse ancient synteny boundaries.**
Test: candidate breakpoint overlaps conserved synteny-block edge in Mac / Pangasius / Ictalurus / Silurus assemblies. Expected proportion vs random chromosomal positions.

H5 — **Some homologous old-species regions are predicted polymorphic hotspots.**
Test: homologous breakpoint TE-rich AND conserved boundary AND orientation-or-context discordance — but mark polymorphism status explicitly as unknown in those species. This is the hypothesis the layer enables but does not directly test.

---

## Atlas integration

When this layer's JSON is loaded:
- Atlas `_layers_present` includes `comparative_te_breakpoint_fragility`
- Boundaries page (`_renderRepeatDensityPanel`) gets a new mode toggle: "Compare to other species" — overlays homologous-breakpoint TE density from each loaded comparative species
- Page 16 (cross-species breakpoints) catalogue rows get a new chip showing the comparative architecture class
- Page 14 hypothesis registry (see `SPEC_hypothesis_registry_and_multispecies.md`) consumes this layer as evidence for hypotheses with `predictions[].layer === "comparative_te_breakpoint_fragility"`

---

## How this fits the manuscript story

Three narrative levels:

1. **Applied** (Results 1-5): Gar inversion atlas for breeding. Same as existing manuscript outline.
2. **Mechanistic** (Results 4 or new section): TE-rich breakpoint reuse — these breakpoints are in repeat-rich regions, consistent with recombination-prone fragile sites.
3. **Evolutionary** (Discussion): Ancient chromosome-architecture boundaries may seed modern structural haplotypes. Comparative TE evidence supports this without overclaiming polymorphism in unsequenced species.

Suggested manuscript paragraph (cautious version):
> *"Although structural polymorphism could only be directly tested in C. gariepinus, we used TE density as a first-pass proxy for breakpoint fragility across homologous regions in comparative catfish assemblies. Several Gar polymorphic inversion boundaries overlapped TE-enriched synteny-block edges that were conserved or rearranged across deeper catfish lineages, suggesting that present-day structural haplotypes can reuse ancient, repeat-rich chromosome architecture."*

Even more cautious:
> *"We confirm polymorphism in C. gariepinus, and we predict which homologous breakpoint regions in other catfish are likely structurally unstable based on TE-rich conserved synteny boundaries."*

---

## Connection to other specs

- `SPEC_OVERVIEW_multispecies_architecture.md` — top-level integration showing how this layer composes with miniprot anchors, BUSCO anchors, phylogenetic tree, and synteny graph
- `SPEC_phylogenetic_tree_integration.md` — uses the tree to place each rearrangement event on a branch (lineage-specific, ancestral, etc.)
- `SPEC_busco_anchors.md` — alternative cross-species anchor signal; sparser than miniprot at breakpoints but works at deeper divergences
- `SPEC_hypothesis_registry_and_multispecies.md` — hypothesis-and-evidence object model that consumes this layer

---

## Effort estimate

Per module (rough):
- Manifest builder: 1 day
- TE normalizer: 1-2 days (the formats are messy)
- Breakpoint window builder: half a day
- Density calculator: 1 day
- Comparative mapping + density: 2-3 days (depends on whether synteny is ready)
- Classification: 1 day
- JSON export: half a day
- Plots: half a day

Total: ~7-10 days of focused work, NOT one chat session. Best done across several chats, with intermediate checkpoints (e.g., manifest + normalizer working before density).

The atlas-side wiring (consume the JSON in the boundaries page) is much smaller — probably 1 chat once the JSON exists.
