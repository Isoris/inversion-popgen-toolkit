# sd_substrate/ — NAHR substrate detection (two angles + concordance)

## Scope

Detect inverted segmental duplications flanking a candidate's two
breakpoints. SD pairs in inverted orientation are the NAHR substrate;
their presence → recurrent inversion via non-allelic homologous
recombination. Absence → NHEJ mechanism (single-origin).

## Architecture (pass 19, 2026-04-24)

Three independent scripts. Each writes its own structured block
(separate provenance, separate registry entry). Disk layout per
candidate:

```
<outdir>/<candidate_id>/
├── structured/
│   ├── sd_substrate_minimap2.json     ← Angle A block
│   ├── sd_substrate_biser2.json       ← Angle B block
│   └── sd_substrate_concordance.json  ← derived verdict
└── sd_substrate/
    ├── minimap2/
    │   ├── flanks.fa                  ±50 kb extracted flanks
    │   ├── flanks.fa.fai
    │   ├── self_align.bam             sorted indexed BAM (load in IGV)
    │   ├── self_align.bam.bai
    │   ├── self_align.cmd.txt         exact invocation for audit
    │   ├── hits.tsv                   parsed alignment table
    │   ├── inverted_hits.bed          BED9 itemRgb: red=-, green=+
    │   └── region.bed                 candidate breakpoints in blue
    └── biser2/
        ├── flanking_sds.tsv           BISER2 rows straddling both sides
        └── flanking_sds.bed           BED9: red=inverted, green=direct
```

No tmp/. Everything persists to the candidate's output directory.

## Scripts

| Script | Role | Cost | Output block |
|---|---|---|---|
| `sd_substrate_minimap2.R` | Angle A: de novo minimap2 self-alignment on ±50 kb flanks | ~30 s/candidate | `sd_substrate_minimap2` |
| `sd_substrate_biser2.R` | Angle B: BISER2 catalog lookup for straddling SDs | <0.1 s/candidate | `sd_substrate_biser2` |
| `sd_substrate_concordance.R` | Reads A + B, emits joint verdict | <0.1 s/candidate | `sd_substrate_concordance` |

## minimap2 invocation (Angle A)

```
minimap2 -ax asm10 -X --eqx --secondary=yes -N 50 -p 0.1 \
    flanks.fa flanks.fa | samtools sort -O bam -o self_align.bam
samtools index self_align.bam
```

Flag rationale (audit trail preserved in `self_align.cmd.txt`):

- `-ax asm10` — assembly-to-assembly preset, ≤10% divergence. Catches
  SDs at ≥90% identity, the biologically meaningful NAHR substrate
  range. Override with `--mm2_preset asm5` (stricter, >95% identity)
  or `asm20` (broadest).
- `-X` — skip self/dual mappings. Essential for self-alignment,
  otherwise every sequence produces a trivial 100% hit to itself.
- `--eqx` — use = / X CIGAR ops so identity is exactly computable.
- `--secondary=yes` — keep secondary alignments. **Necessary** for SD
  detection: an SD pair produces one primary + one secondary alignment
  (A→B primary means B→A comes back as secondary). Without this, half
  of every SD pair is lost.
- `-N 50` — report up to 50 secondaries per query (default 5 is too
  few for repeat-rich regions).
- `-p 0.1` — secondary-to-primary score ratio threshold (default 0.8
  drops partial-identity SD copies that have score < 80% of primary;
  post-filter handles quality at higher resolution).

## BISER2 invocation (Angle B)

No external command — pure TSV parsing + interval overlap. Reads the
BISER2 output catalog (9 cols: `chr1 start1 end1 chr2 start2 end2
orientation identity sd_type`; 8-col variant also handled) and finds
SD pairs where one copy lies in `[left_bp - flank, left_bp + flank]`
and the other in `[right_bp - flank, right_bp + flank]` (either order).

## Concordance verdicts

Per-side matrix (left / right treated independently, then combined):

| Angle A says | Angle B says | Concordance | Confidence |
|---|---|---|---|
| NAHR_CANDIDATE | strong/weak | `agree_nahr` | high |
| NAHR_CANDIDATE | subthreshold | `minimap2_stronger_nahr` | medium |
| NAHR_CANDIDATE | direct_only / none | `disagree_biser2_missing` | high (A wins) |
| NHEJ_CANDIDATE | none / direct_only | `agree_nhej` | high |
| NHEJ_CANDIDATE | strong / weak | `disagree_minimap2_missing` | medium (review) |
| NAHR_POSSIBLE | strong / weak | `agree_nahr_weak` | medium |
| NAHR_POSSIBLE | none | `minimap2_only_weak` | low |
| COMPLEX | any | `complex_flagged` | low |
| (A missing) | strong/weak | `biser2_only_nahr` | medium |
| (A missing) | none/direct_only | `biser2_only_nhej` | medium |
| (B missing) | NAHR_* | `minimap2_only_nahr` | medium |
| (B missing) | NHEJ_* | `minimap2_only_nhej` | medium |
| (both missing) | | `no_data` | low |

Candidate-level confidence is **the worst of left and right**. Policy:
minimap2 is the trusted angle — when both disagree, minimap2 wins for
the primary call, but the disagreement is preserved in the block for
manuscript-level review.

## Running

### Both angles + concordance (typical)

```bash
# Angle A: minimap2
Rscript sd_substrate_minimap2.R \
    --candidate LG28_cand_1 --chrom C_gar_LG28 \
    --left_bp 15115243 --right_bp 18005891 \
    --ref_fasta /path/to/ref.fasta \
    --outdir /path/to/phase8_blocks

# Angle B: BISER2
Rscript sd_substrate_biser2.R \
    --candidate LG28_cand_1 --chrom C_gar_LG28 \
    --left_bp 15115243 --right_bp 18005891 \
    --biser2_tsv /path/to/biser2_results.tsv \
    --outdir /path/to/phase8_blocks

# Concordance
Rscript sd_substrate_concordance.R \
    --candidate LG28_cand_1 \
    --outdir /path/to/phase8_blocks
```

### Fast BISER2-only sweep (genome-wide first pass)

Skip Angle A entirely:

```bash
Rscript sd_substrate_biser2.R --candidate ... --biser2_tsv ... --outdir ...
Rscript sd_substrate_concordance.R --candidate ... --outdir ...
# concordance says biser2_only_nahr / biser2_only_nhej
```

### Targeted minimap2 confirmation

After identifying top candidates from the BISER2 pass, run Angle A
only on them and re-run concordance to promote to `agree_nahr` /
high confidence:

```bash
for cid in LG28_cand_1 LG28_cand_2 LG01_cand_5; do
    Rscript sd_substrate_minimap2.R --candidate $cid ... --ref_fasta ...
    Rscript sd_substrate_concordance.R --candidate $cid --outdir ...
done
```

### Orchestrated via LAUNCH_group_cheats.sh

The phase_8 per-candidate dispatcher `../../orchestrator/LAUNCH_group_cheats.sh`
calls all three scripts in order (A → B → concordance). Set `NO_MINIMAP2=1`
to skip Angle A for the fast sweep mode.

## Inspecting results in IGV / JBrowse

`self_align.bam` is sorted and indexed, ready to load. Load
`inverted_hits.bed` and `region.bed` as annotation tracks — they use
BED9 itemRgb so red/green/blue colouring applies automatically in
IGV's "feature colour from file" mode.

Red bars on the BED track are the NAHR-relevant inverted hits. Green
bars are direct repeats (not NAHR substrate but useful to see
architecture). Blue bars are the candidate breakpoint positions
themselves — navigate to them first, then zoom in to see the
surrounding SD architecture.

## Why three scripts instead of one merged

Earlier merged design (pass 18's first attempt) baked both angles into
one script. Two problems:

1. **Independent provenance.** If one angle fails or is skipped,
   the other's result should still be visible in the registry as its
   own block with its own `source_script` tag, not silently hidden
   inside a merged block with some fields NA.
2. **Inspection workflow.** Running BISER2-only as a fast sweep, then
   cherry-picking candidates for minimap2, is the realistic usage
   pattern. Separate scripts make that natural; a merged script would
   require conditional flags and an --only_angle_b kind of option
   that still hides its results under a shared block type.

The concordance script is where the two angles are combined — but
it's a *reader* of two separate blocks, not a writer that merges
them.

## Registry contracts

Each block schema has `keys_extracted` pulling the flat keys out. The
concordance block additionally mirrors its verdict into the existing
`mechanism.schema.json` keys (`q4_sd_concordance`,
`q4_mechanism_confidence`) for backward compatibility with
`STEP_04_assign_structural_class.py`.

## Archive provenance

The legacy pre-pass-18 scripts `cheat14_self_align.sh`,
`cheat14_repeat_architecture.R`, and `cheat27_sd_nahr_substrate.R`
are preserved under `_archive_superseded/cheat14_cheat27_merged_into_sd_substrate/`.

The pass-18 single-file merged `sd_substrate.R` was a transitional
design that pass-19 split back into three (per user request for
independent blocks and inspectable intermediates). It has been
removed from the live tree; the split scripts in this directory
supersede it.
