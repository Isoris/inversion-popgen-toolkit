# SV Evidence Producers

**Inputs:** HPC-side BAMs, VCFs, BAM-evidence track JSON, locked karyotype labels
**Outputs:** Per-candidate folder layers consumable by the atlas SV evidence page (drag-drop)
**Cohort:** 226-sample pure *C. gariepinus* hatchery only (NOT F1 hybrid; NOT *C. macrocephalus* wild)
**Run on:** LANTA HPC, account `lt200308`, working dir `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/`

## Files

```
producers/sv_evidence/
├── write_candidate_folder.py            ← shared library: per-candidate folder writer + manifest
├── STEP_SV_GT_AGG_aggregate_genotype_counts.py   ← emits sv_genotype_counts.json
├── STEP_SV_EVID_COMB_emit_combinations.py        ← emits sv_evidence_combinations.json
└── run_sv_evidence_pipeline.slurm       ← LANTA SLURM wrapper that runs both
```

## Output layout

Per-candidate folder under `--out-root`:

```
<out-root>/
└── C_gar_LG28/
    ├── manifest.json                              ← chrom-level index of candidates + layers
    └── candidates/
        └── INV_LG28_003/
            ├── sv_genotype_counts.json            ← steps 1-4 atlas page
            ├── sv_evidence_combinations.json      ← step 5 UpSet panel
            └── sv_support_by_sample.json          ← future (step 6 heatmap)
```

The atlas page accepts **either** drag-dropping the candidate folder directly (folder-walk via `webkitGetAsEntry()` picks up every JSON inside) **or** dropping a single layer at a time.

## Step-1 producer: `STEP_SV_GT_AGG`

Aggregates per-SV genotype counts within a candidate's flanking + boundary + body zones. For each SV:

1. Bin samples into H1/H1, H1/H2, H2/H2 from the locked karyotype labels.
2. Compute the AA / AB / BB / miss counts per group.
3. Run Fisher's exact test (H1/H1 vs H2/H2, log-space hypergeometric) to get OR + raw p-value.
4. Apply Benjamini-Hochberg FDR correction across all SVs in the candidate window.
5. Classify zone (left_flank / left_boundary / inversion_body / right_boundary / right_flank).
6. Apply the pattern-label decision rule to assign one of:
   - `canonical_breakpoint_marker`
   - `dominant_presence_marker`
   - `het_specific_marker` (NB: gates *before* FDR — see explanation below)
   - `sub_haplotype_marker`
   - `internal_linked_marker`
   - `uninformative`
7. Build `boundary_summary.{left,right}.by_sv_type` aggregate counts for the right-rail tables.

### Important: het-specific markers and FDR

The pattern classifier's `het_specific_marker` rule fires *before* the FDR gate because by construction a het-specific marker is invisible to the H1/H1-vs-H2/H2 Fisher test (both homozygote groups carry at the same rate, near zero), so its FDR will be ~1 even though the marker is real and important. Recognized by group-fraction structure (h12 carriers ≥ 50%, both homs ≤ 10%) instead of the OR.

If you want het-specific markers gated by a real test, you'd need a separate Fisher contrast (het vs homs combined) — not in v1 of the schema, but a clean addition for v2.

### Inputs

| flag | required | description |
|------|----------|-------------|
| `--vcf` or `--tsv` | one | merged DELLY+Manta VCF, OR tabular TSV with `sv_id chrom position_bp end_bp sv_type sample_id GT [quality] [callers]` |
| `--candidate` | yes | candidate JSON (must include `candidate_id`, `chrom`, `boundary_left_bp`, `boundary_right_bp`, optional `zone_definitions_bp`) |
| `--karyotype` | yes | TSV: `sample_id <tab> label` where label ∈ {HOMO_1, HET, HOMO_2} |
| `--out-root` | yes | output root directory |
| `--fdr-cutoff` | no | FDR threshold for "associated" classification (default 0.05) |
| `--indent` | no | JSON indent (default: compact one-line) |

### Example

```bash
python3 producers/sv_evidence/STEP_SV_GT_AGG_aggregate_genotype_counts.py \
  --vcf /scratch/lt200308-agbsci/.../sv_calls/INV_LG28_003.merged.vcf \
  --candidate /scratch/lt200308-agbsci/.../candidates/INV_LG28_003.json \
  --karyotype /scratch/lt200308-agbsci/.../karyotype/INV_LG28_003.locked_labels.tsv \
  --out-root  /scratch/lt200308-agbsci/.../atlas_outputs/sv_evidence
```

## Step-2 producer: `STEP_SV_EVID_COMB`

Builds the UpSet panel layer. Reads a per-sample × per-evidence-type 0/1 matrix (TSV) and groups samples by their evidence-set; the top-N intersection sizes become the bars in the UpSet.

### Inputs

| flag | required | description |
|------|----------|-------------|
| `--evidence` | yes | per-sample × per-evidence-type 0/1 matrix TSV |
| `--evidence-types` | yes | JSON list of `{id, label, side, kind, tier}` — canonical row order |
| `--candidate-id` | yes | embedded in output JSON |
| `--chrom` | yes | output folder path |
| `--out-root` | yes | output root |
| `--top-n` | no | keep top-N combinations (default 20) |
| `--indent` | no | JSON indent |

### Producing the evidence matrix

The 0/1 matrix is upstream of this producer — it comes from MODULE_5A2 / S7 (your `phase_8_comparative_breakpoint_fragility` work). For each sample × candidate, you need 0/1 flags for each of:

- `left_SA` — split-read clipping at left boundary (samtools `SA` tag clusters)
- `right_SA` — same for right boundary
- `left_PE` — paired-end discordant orientation, left
- `right_PE` — same for right
- `Manta_INV_GT` — Manta INV genotype call passes
- `DELLY_INV_GT` — DELLY2 INV genotype call passes
- `MAPQ0_left` — low-MAPQ region present at left boundary
- `MAPQ0_right` — same for right

Quentin's own evidence-types vocabulary may differ — the producer doesn't hard-code these. Whatever you put in `--evidence-types` JSON, the producer projects the matrix columns onto that order.

### Example

```bash
python3 producers/sv_evidence/STEP_SV_EVID_COMB_emit_combinations.py \
  --evidence       /scratch/lt200308-agbsci/.../bam_evidence/INV_LG28_003.evidence_matrix.tsv \
  --evidence-types /scratch/lt200308-agbsci/.../bam_evidence/evidence_types.json \
  --candidate-id   INV_LG28_003 \
  --chrom          C_gar_LG28 \
  --out-root       /scratch/lt200308-agbsci/.../atlas_outputs/sv_evidence
```

## SLURM wrapper

`run_sv_evidence_pipeline.slurm` wires both producers up for one candidate:

```bash
sbatch producers/sv_evidence/run_sv_evidence_pipeline.slurm INV_LG28_003
```

Or batch-submit a whole chromosome:

```bash
for cid in INV_LG28_001 INV_LG28_002 INV_LG28_003 ...; do
  sbatch producers/sv_evidence/run_sv_evidence_pipeline.slurm "$cid"
done
```

Edit the `PROJECT_ROOT` / `INPUTS` / `OUT_ROOT` paths at the top of the SLURM script to match your tree.

## Atomic writes

Both producers write through `write_candidate_folder.atomic_write_json`, which:
1. Validates the payload via `json.dumps` *before* touching disk.
2. Writes to a sibling `.tmp` file.
3. `os.replace` to the final name (atomic on POSIX).

This means the atlas page can never read a half-written JSON — interrupted producer jobs leave `.tmp` files (which the atlas ignores) but never partial layers.

## Manifest

After each write the producers regenerate a chrom-level `manifest.json` listing every candidate folder and which layers each has. The atlas can use this for a future "load the whole chromosome" mode (out of scope for now).

## End-to-end test

```bash
python3 tests/producers/test_step7_producers.py
```

Synthesizes a small but realistic test set (50 fish, 4 SVs, 8 evidence types, 4 distinct evidence patterns), runs both producers, and validates the output JSONs against the same atlas-side validators that run in the browser. **36/36 passing.**

## Round-trip with atlas page

After running the producers, drag-drop the per-candidate folder onto the atlas SV evidence tab:

```
data/C_gar_LG28/candidates/INV_LG28_003/    ← drag this folder
```

The atlas walks the folder, dispatches each `*.json` by `format_version`, and populates both layers in one drop.
