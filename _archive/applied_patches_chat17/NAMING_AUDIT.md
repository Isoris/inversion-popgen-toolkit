# Naming audit — non-scientific / confusing names

Systematic audit of the active tree for names that trip up readers. Based
on RENAMING.md's stated goal ("The manuscript and the codebase use
scientific vocabulary; the rename makes the codebase match.") and the
specific case the user flagged in chat 17 (`inv_like_dt` / `window_dt`).

**RENAMING.md already exists** under `phase_2_discovery/2c_precomp/` and
lists the master rename plan. Most filename renames are done; the
**internal identifiers and CLI flags were not applied beyond three core
scripts**. This audit confirms what's still pending and groups by risk.

## Tier A — Ship with chat 17 (safe, mechanical, high signal/noise)

All already planned in RENAMING.md §2.3 and §2.4.

### 1. `flashlight` → `sv_prior` everywhere (17 active files)

**Status:** filenames renamed (`STEP_C00_build_sv_prior.R`), but internal
identifiers and CLI flags still use `flashlight`. Notable offenders:

- `STEP_C01d_candidate_scoring_*.R` — `--flashlight_dir` CLI arg (accepted-but-ignored per its own bugfix comment)
- `STEP_C01g_boundary_catalog_*.R` — `--flashlight` CLI, `flashlight_dir` R var, reads `sv_flashlight_<chr>.rds` as fallback
- `phase_2_discovery/2c_precomp/patches/patch_C01*_flashlight.R` × 6 — all 6 flashlight patch files
- `4b_group_proposal/lib_step03_seed_loader.R` — "Tier 1 seeding from STEP03 + flashlight"
- `4e_final_classification/compute_candidate_status.R` — legacy BK-key name mapping notes

**Patch approach.** Literal string replace: `flashlight` → `sv_prior`,
`FLASHLIGHT_LOADER` → `SV_PRIOR_LOADER`, CLI flag
`--flashlight_dir` → `--sv_prior_dir`, `--flashlight` → `--sv_prior`.
Keep back-compat: the chat-17 C01g launcher still passes `--flashlight`;
will need to update that too.

**Expected substitutions:** ~50–80 total across 17 files.

### 2. `FLASHLIGHT_LOADER` env var → `SV_PRIOR_LOADER` (6 files)

Pure env var rename. Touches `00_inversion_config.sh`, `pipeline_bridge.sh`,
and 4 consumer scripts.

### 3. `snake_id` / `snake_phase` → `region_id` / `extension_phase` (5 files)

Already partially applied. Remaining active consumers are mostly in
phase 4 / phase 5 figures that read the regions table.

### 4. `core_family` → `scale_tier` (4 files)

Per RENAMING.md §2.3: "these are parameter scales, not biological
families — 'family' was misleading." Mechanical sed.

## Tier B — Ship with chat 17 but confirm once (small blast radius)

### 5. `snake_regions_multiscale/` directory name (8 files reference it)

RENAMING.md §2.1 suggests `regions/` but marks it "NOT YET APPLIED" with
a deliberate note that changing the dir name touches file paths in many
launchers. The chat-17 launchers all use this path. Recommendation:
**defer to chat 18 or 19.** If renamed now, every launcher written in
chat 17 needs its default path updated. The 8 references are all
fixable but the cost/benefit is weaker than the identifier renames.

### 6. `cheat*` files in `4d_group_dependent/` (5 files)

The 5 `cheat27–30` files + `cheat6_ancestry_jackknife_v934.R`. Per
RENAMING.md §3, these become `test27`, `test28`, etc. **Filename rename
is a separate question from identifier rename** — if you rename the
files, every launcher calling them needs updating too. The chat-17
`LAUNCH_group_cheats.sh` (existing) already calls them by their cheat*
names.

**Recommendation:** apply the in-file identifier rename (`cheat` → `test_NN`
in variable/column names where non-filename) and **keep the filenames
as-is** until a coordinated filename pass.

### 7. `bloc_` / `blocs_dt` / `bloc_ribbon` in fst_dxy plot scripts (3 files)

RENAMING.md §12 lists `BLOC_NEAR_FRAC` → `BLOCK_NEAR_FRAC` etc. —
plain English "block" is scientific, "bloc" was a typo/shorthand that
stuck. 3 files in `fst_dxy/` and `phase_5_followup/figures/`.

## Tier C — Leave alone (accurately named or too risky)

### `inv_likeness` (column, not table)

The column itself IS the inv-likeness score. Don't rename the column.
(The variable `inv_like_dt` → `window_dt` rename already queued.)

### `peel` / `peeling`

Documented scientifically in `STEP_D09n_peeling_diagnostic.R` header —
"sample peeling" means zeroing out specific samples to test block
robustness. Idiosyncratic but internally documented.

### `snake3` in GHSL scripts (`STEP_C04_snake3_ghsl_v6.R`)

GHSL = "Group-level Haplotype Similarity by Locus" is the scientific
name. `snake3` refers to the historical Snake 1/2/3 architecture. The
filename still carries it for continuity with the v5 → v6 rewrite
history. Could rename `STEP_C04_ghsl_v6_compute.R` but the chat-17
LAUNCH launchers reference the current name; rename is a coordinated
move with the launchers.

### `core` in other modules (MODULE_4A score_breeding_concern, etc.)

Different meaning outside the inversion pipeline — "conservation core"
is the Fish10K 464-species core, "breeding concern" is a scientific
metric. Unrelated to the inversion `core` → `seeded_region` rename.
Not part of this audit.

### `hatchery` / `wild` (`--mode` flag)

Scientific terms for population type. Not renaming.

### `STAIR_*`, `NN_*`, `CONSENSUS_*`, `LANDSCAPE_*`, `PEEL_*`

Per RENAMING.md §12: "these describe the parameter's function, not a
historical nickname." Keep.

## Recommended chat-17 scope

Ship one patch covering Tier A (#1–#4):

- `flashlight` → `sv_prior` (identifier + CLI flag rename; back-compat
  on both the `--flashlight_dir` and the reader fallback)
- `FLASHLIGHT_LOADER` → `SV_PRIOR_LOADER`
- `snake_id` / `snake_phase` → `region_id` / `extension_phase`
- `core_family` → `scale_tier`

Defer Tier B and skip Tier C.

Expected blast radius: ~25 files, ~100–150 substitutions. Same pattern
as `patch_rename_window_dt.sh`: Python literal replace, idempotent,
`.bk_chat17` backup, `_rcheck.py` parse-check all touched files before
and after.
