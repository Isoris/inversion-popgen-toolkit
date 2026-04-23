# SUPERSEDED — MODULE_5E_Inversion_HOBS

**Archived:** 2026-04-24
**Superseded by:** `inversion_modules/phase_qc_shelf/STEP_Q07b_hobs_per_group.sh` + `STEP_Q07c_hobs_windower.sh` (Engine H, shipped in phase_qc_shelf v3.1)

## Why archived

MODULE_5E was a simple, single-subset wrapper around Claire Mérot's
Hobs/HWE sliding-window logic. It took an already-computed ANGSD
`-doHWE` output and a `subset_label` and wrote windowed Hobs/Hexp/F
summaries with outlier burden.

The Q07b + Q07c pair in `phase_qc_shelf/` supersedes it scientifically:

| Dimension | MODULE_5E (archived) | Q07b + Q07c (canonical) |
|---|---|---|
| **Group awareness** | One subset at a time, caller decides pooling | Runs ANGSD `-doHWE` **per karyotype group** (Hom1 / Het / Hom2) separately using `invgt_assignments` — this is the scientifically correct Mérot approach because `Hexp` at an inversion SNP should be computed from *within-group* allele frequencies |
| **Engine** | Pure R, single pass | Patched ANGSD (`Isoris/angsd_fixed_HWE`) + the `hobs_windower` C binary |
| **Scales** | One windowed summary | Seven multi-scale sliding windows |
| **Signature tested** | Site-level F outlier burden | `HoverE_HET → 2.0`, `HoverE_HOM → 0` at arrangement-differentiating SNPs (the crisp diagnostic for a real balanced inversion) |
| **Downstream integration** | Nothing reads its outputs | Writes `hobs_sites.*.tsv.gz` + `hobs_win.*.<SCALE>.tsv.gz` consumed by `phase_5_followup/analysis/04_candidate_overlay.py` and `phase_5_followup/figures/05_plot_hobs_hwe.R` |
| **Active callers in repo** | Zero | Q07b/Q07c are wired in the phase_qc_shelf runner chain |

Grep on 2026-04-24 confirmed MODULE_5E had **zero callers** anywhere in
the live tree — no `LAUNCH_*.slurm`, no orchestrator, no README
referencing it as part of a chain. It wasn't broken; it was simply
replaced by a more correct implementation and never deleted.

## Preserved rather than deleted

- The R logic in `run_hobs_confirmation_module.R` is a compact
  reference for the basic single-subset Hobs windowing approach, in
  case a future use case ever needs a lightweight Hobs check without
  the per-group refinement.
- The filename / module identity is a recognisable pointer from older
  chat logs ("the 5E module", "HOBSDIR") — preserving it by name
  helps when grepping historical design docs.

## If you need per-group Hobs

Use `phase_qc_shelf/STEP_Q07b_hobs_per_group.sh` followed by
`phase_qc_shelf/STEP_Q07c_hobs_windower.sh`. See the phase_qc_shelf
README for the full chain.
