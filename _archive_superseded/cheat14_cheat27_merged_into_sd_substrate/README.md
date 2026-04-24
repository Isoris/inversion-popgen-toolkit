# cheat14 + cheat27 → merged into `sd_substrate.R` (pass 18, 2026-04-24)

## What this archive contains

Three files preserved for provenance after being merged into one live
script:

- **`cheat27_sd_nahr_substrate.R`** — moved from
  `inversion_modules/phase_8_evidence_biology/q4_mechanism/`. Was a
  function-library-only file (no CLI, no driver) that exported
  `load_biser2()`, `find_flanking_sds()`, `classify_sd_mechanism()`,
  `run_cheat27()`. The launcher in `LAUNCH_group_cheats.sh` called it
  with `Rscript ... --candidate ...`, but because the file had no
  CLI parser or main block, that call was a silent no-op (script
  parsed function definitions and exited). Has been a latent bug since
  pass 15.

- **`cheat14_self_align.sh`** — bash wrapper that extracts ±50 kb
  flanks with `samtools faidx` and runs `minimap2 -X -c --eqx` for
  self-alignment. Was only ever in `_archive/` (never in the live
  tree after the pass-12 restructuring). Copy preserved here for
  comparison with `sd_substrate.R`.

- **`cheat14_repeat_architecture.R`** — R library with
  `extract_breakpoint_flanks()`, `self_align_flanks()`,
  `classify_mechanism()` (NAHR_CANDIDATE / NHEJ_CANDIDATE / ...),
  `annotate_te_context()`, and a `run_cheat14()` wrapper. Was only
  ever in `_archive/` in the live tree. Copy preserved here.

## Why they were merged

`docs/toolkit_audit.md` §D.1 (pre-existing recommendation) noted that
cheat14 and cheat27 answer the same biological question — "Is there
an inverted SD pair flanking this candidate's breakpoints?" — from
two independent angles:

- **Angle A (cheat14):** *de novo* minimap2 self-alignment of ±50 kb
  flanks. Finds inverted repeats directly from the reference FASTA.
- **Angle B (cheat27):** lookup into the pre-computed BISER2 SD catalog.

Their **whole design** was to be combined into a `concordance` verdict:
agree_nahr / agree_nhej / disagree_biser2_nahr / disagree_minimap2_nahr.
cheat27 had a `cheat14_mechanisms` argument specifically for that
cross-check, but there was no plumbing to actually feed cheat14's
output into cheat27. Keeping them separate meant two drivers + external
coordination for one logical output.

The merge (`q4_mechanism/sd_substrate.R`) runs both angles and emits
the concordance verdict in one pass. Angle A can be skipped with
`--no_minimap2` when LANTA is busy and only the fast BISER2 lookup is
needed; the concordance column will then report `biser2_only_nahr` /
`biser2_only_nhej` to make the partial-coverage state visible.

## How to resurrect if needed

If the merged `sd_substrate.R` ever turns out to be wrong and you need
the original separated logic back:

```bash
# Restore cheat27 as function library
cp _archive_superseded/cheat14_cheat27_merged_into_sd_substrate/cheat27_sd_nahr_substrate.R \
   inversion_modules/phase_8_evidence_biology/q4_mechanism/

# Restore cheat14 into a q4_mechanism sibling location
cp _archive_superseded/cheat14_cheat27_merged_into_sd_substrate/cheat14_repeat_architecture.R \
   inversion_modules/phase_8_evidence_biology/q4_mechanism/
```

Then re-wire LAUNCH_group_cheats.sh back to the two separate calls
and remove `sd_substrate.R`.

## See also

- `inversion_modules/phase_8_evidence_biology/q4_mechanism/sd_substrate/`
  — the live implementation. Pass 18's single-file merged script was
  split into three (minimap2, biser2, concordance) in pass 19.
- `inversion_modules/phase_8_evidence_biology/q4_mechanism/README.md` —
  updated to reference the three sd_substrate/ scripts
- `docs/toolkit_audit.md` §D.1 — the pre-existing merge proposal that
  pass 18 acted on
