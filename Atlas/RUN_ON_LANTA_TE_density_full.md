# Runbook — Comprehensive TE density layer

**Layer**: 78 per-class densities + 5 aggregates + 2 TSD layers per chromosome
**Source**: EDTA `TEanno.gff3` (homology + structural) + `intact.gff3` (TSD coords)
**Output**: per-species dirs at
  `/project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results/{Cgar|Cmac}_TE_density/`
**Schema**: repeat_density_v2 (loadable by Inversion Atlas drag-drop)
**Compute**: login-node only, ~3-5 min per species
**Status**: ⚪ ready to run · login node accessible during HPC maintenance

---

## What this layer is — and how it differs from the TSD-only layer

The TSD-only layer (turn 1, already shipped) gave you ~4k Gar / ~1.8k Mac
TSD events at 5-kb resolution. Sparse but breakpoint-diagnostic.

This comprehensive layer is the FULL TE density picture, derived from
EDTA's master `TEanno.gff3` (homology + structural, all 26 classes) plus
the `intact.gff3` for the explicit-coord TSDs. Per chromosome JSON:

| Layer group | Count | What it shows |
|---|---|---|
| Per-class density | 26 | bp-coverage fraction per 5-kb window, per TE class |
| Per-class young | up to 26 | same but identity ≥ 0.95 (recently inserted) |
| Per-class old | up to 26 | same but identity < 0.95 (ancient remnants) |
| `all_TE` | 1 | sum of all class densities (saturates → 1) |
| `young_TE_all` | 1 | sum of young across classes |
| `old_TE_all` | 1 | sum of old across classes |
| `insertion_count` | 1 | raw count of TE records per window |
| `intact_element_count` | 1 | count of intact full-length elements (from intact.gff3) |
| `target_site_duplication` | 1 | TSDs with explicit coords (= turn 1 layer) |
| `target_site_duplication_all` | 1 | + TIR-attribute-derived TSD positions |

Total ~85 layers per chromosome JSON. File size ~6-8 MB per chrom.
Total across both species: ~330-450 MB on disk.

## When to use which layer in the boundaries panel

- **General overview** (is there ANY TE here?) → `all_TE`
- **Substrate for ectopic recombination** (where could TE-mediated
  rearrangements form?) → `young_TE_all` (recently active TEs are the
  template for ectopic recombination)
- **Class-specific hypothesis** (this candidate looks Gypsy-mediated) →
  `Gypsy_LTR_retrotransposon` or its `__young` variant
- **Insertion event density** (is this a TE-rich integration zone?) →
  `insertion_count`
- **Intact full-length elements** (any complete TEs that could provide
  homologous templates?) → `intact_element_count`
- **Breakpoint footprint detection** → `target_site_duplication_all`
  (the denser TSD layer)

## Step-by-step

### 1. SSH and place the scripts

```bash
ssh qandres@lanta-login01
cd /project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results

# Drop these two files into this directory (or any folder):
#   STEP_RD_TE_density_full.R
#   RUN_RD_TE_FULL_BOTH_SPECIES.sh
# The wrapper finds the R script via its own dirname.

mamba activate assembly
which Rscript    # /lustrefs/disk/.../assembly/bin/Rscript
```

### 2. Verify input files

```bash
# Both species' TEanno.gff3 (~270-290 MB each)
ls -lh fClaHyb_Gar_LG.fa.mod.EDTA.TEanno.gff3
ls -lh fClaHyb_Mac_LG.fa.mod.EDTA.TEanno.gff3

# Intact.gff3 already verified during turn 1 — but recheck:
ls -lh fClaHyb_Gar_LG.fa.mod.EDTA.intact.gff3
ls -lh fClaHyb_Mac_LG.fa.mod.EDTA.intact.gff3

# .fai files
ls -lh fClaHyb_Gar_LG.fa.mod.fai
ls -lh fClaHyb_Mac_LG.fa.fai
```

### 3. Run

```bash
bash RUN_RD_TE_FULL_BOTH_SPECIES.sh
```

This is a longer-running job than turn 1. Expected output per species:

```
═══════════════════════════════════════════════════════════════
  Clarias gariepinus
═══════════════════════════════════════════════════════════════
  teanno_gff3 : .../fClaHyb_Gar_LG.fa.mod.EDTA.TEanno.gff3  (276M)
  intact_gff3 : .../fClaHyb_Gar_LG.fa.mod.EDTA.intact.gff3  (4.3M)
  fai         : .../fClaHyb_Gar_LG.fa.mod.fai
  out_dir     : .../Cgar_TE_density
[STEP_RD_TE_FULL] starting
[STEP_RD_TE_FULL] reading TEanno.gff3 ...
[STEP_RD_TE_FULL] TEanno records: ~470000
[STEP_RD_TE_FULL] records with identity field: NN/470000
[STEP_RD_TE_FULL] reading intact.gff3 ...
[STEP_RD_TE_FULL] TE classes detected: 26
  - Bel_Pao_LTR_retrotransposon
  - CACTA_TIR_transposon
  ...
  wrote C_gar_LG01_repeat_density_TEfull.json · n_windows=10225 · classes=83 · size=6500KB
  wrote C_gar_LG02_repeat_density_TEfull.json · ...
  ...
  ✓ Clarias gariepinus: 28 JSONs (~180MB) in 230s
```

If the script appears to hang for > 10 min on "reading TEanno.gff3 ...",
that's normal — fread on a 280 MB GFF takes 30-60 sec, and the per-class
foverlaps loops follow. Use `top` in another shell to confirm Rscript is
churning CPU.

### 4. Spot-check one JSON per species

```bash
JSON_DIR_GAR=/project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results/Cgar_TE_density
JSON_DIR_MAC=/project/lt200308-agbsci/02-TE_catfish/EDTA2.2_results/Cmac_TE_density

# How many classes ended up in LG01?
jq '.chromosomes[0] | {chrom, n_windows, n_classes,
     classes_emitted: (.by_class | keys | length),
     all_TE_max: .by_class.all_TE.max_density,
     intact_count: .by_class.intact_element_count.n_records,
     tsd_intact_count: .by_class.target_site_duplication.n_records,
     tsd_all_count: .by_class.target_site_duplication_all.n_records}' \
   $JSON_DIR_GAR/C_gar_LG01_repeat_density_TEfull.json

# What does Gypsy young density look like at the busiest window?
jq '[.chromosomes[0].by_class.Gypsy_LTR_retrotransposon__young.densities[]
     | select(. > 0)] | sort | reverse | .[0:10]' \
   $JSON_DIR_GAR/C_gar_LG01_repeat_density_TEfull.json

# Quick top-class scan: which classes have the highest max_density?
jq -r '.chromosomes[0].by_class | to_entries
       | map({key: .key, max: .value.max_density})
       | sort_by(.max) | reverse | .[0:15]
       | .[] | "\(.max | tostring)  \(.key)"' \
   $JSON_DIR_GAR/C_gar_LG01_repeat_density_TEfull.json
```

Sanity expectations:
- **n_classes** in JSON = 26 base + ~26 young + ~26 old + 5 aggregates + 2 TSD ≈ 80-85.
  The exact count depends on which classes have records on each chromosome.
- **all_TE max_density** should be 0.5-0.95 in TE-rich regions
  (saturating well below 1.0 because not all classes co-occur).
- **n_classes_emitted varies per chrom** (e.g. some rare classes may be
  absent on small chromosomes). The atlas's class picker discovers all
  available classes per JSON automatically.

### 5. File size check

```bash
du -sh $JSON_DIR_GAR $JSON_DIR_MAC
ls -lhS $JSON_DIR_GAR | head -5    # largest first
ls -lhS $JSON_DIR_MAC | head -5
```

If any chromosome's JSON is > 12 MB, the script over-emitted (some bug)
or that chrom has unusually many records. If a JSON is < 1 MB, it under-
emitted (likely because TEanno had no records for that chrom — verify).

### 6. Load into the Inversion Atlas

The 6-8 MB JSON per chrom is bigger than any prior layer. Two approaches:

**(a) One chromosome at a time** (recommended): drag the LG of interest's
TEfull JSON onto the boundaries page. The atlas merges it into state
within ~1 second.

**(b) Bulk-load via IndexedDB**: drag all 28 (or 55 across species) onto
the load zone. IndexedDB write takes ~10-20 sec for 200 MB total. Auto-
restore on next page load works fine but the initial cache populate is
visible.

If atlas performance suffers (slow class picker, slow chrom switches),
the layer is too large for IDB on your browser. Mitigation: load only
the chroms you're actively analyzing.

### 7. Use on the boundaries panel

1. Tab → 3 (boundaries)
2. Active candidate selected
3. Repeat density panel
4. Class picker — now showing 80+ classes alphabetically
5. Suggested test sequence:
   - Start with `all_TE` to see overall TE topography
   - Switch to `young_TE_all` to see recently-active TE distribution
   - Click the class summary auto-scan flagged classes (top 5-10 will
     surface the candidate-relevant ones)
   - For TE-mediated breakpoint hypothesis: pick the flagged class's
     `__young` variant — that's the recombination substrate

---

## Troubleshooting

**"Object 'foverlaps' not found"** — data.table is too old. Update with
`install.packages('data.table')` in the assembly env or use a newer
mamba env.

**Rscript runs out of memory** — TEanno.gff3 read is the peak (~2-3 GB).
If your login-node session has a memory limit, check with `ulimit -v`
and bump if needed. On LANTA login nodes 4 GB is typically available.

**JSON file too big (> 12 MB)** — open one and check what's eating space.
Most likely culprit: a class with millions of records on a long chrom.
Possible mitigation: drop the `__old` stratification (most records are
old). Tell me and I'll add a `--skip_old` flag.

**Atlas drag-drop fails** — check browser console. The atlas's existing
loader expects `version=2` and `binning_source=scrubber_windows` — both
are emitted by this script. If the loader rejects, paste the console
error and I'll diagnose.

**Class picker overwhelmed (80 classes)** — this is real UX friction.
Mitigation in a follow-up turn: filter the picker to show only classes
with `max_density > 0.05` on the current chrom by default, with a
"show all" toggle.

---

## What next

After this layer ships, planned follow-ups:

- **TE-class auto-flag in the candidate-summary panel**: extend the
  existing class-summary auto-scan to flag young/old strata separately,
  and call out classes with the highest `delta_vs_chrom` in the candidate.
- **Saturation curve for TSD coverage**: sanity that the
  `target_site_duplication_all` layer doesn't miss obvious breakpoint
  flanks vs the literature.
- **Cross-species coverage diff**: load both Gar and Mac TEfull JSONs for
  the same syntenic region; compute a per-window coverage difference;
  potential signal of lineage-specific TE expansion.

Each is a follow-up turn — not part of this run.
