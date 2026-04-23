# Local Install — phase_qc_shelf v2.1 + v3.0 + v3.1

This bundle syncs your LOCAL clone of `inversion-popgen-toolkit` to match what
got installed on HPC on 2026-04-20, plus the v3.0 (registry) and v3.1 (Hobs /
Engine H) layers that shipped in the same tarball but weren't tested on LANTA
before the maintenance window.

## What you're adding

- **v2.1**: No-data strips in Q04, two-layer ideogram stipple, Q09 gap
  characterization, Q09b LD-structure check, `run_chrom.sh` per-chrom driver,
  standardized logging helpers in `00_config.sh`, standardized Q01.
- **v3.0**: Q10 registry bridge (populates interval/sample/results/evidence
  registries), 28-chrom SLURM array driver, `shelf_coords.tsv` seed file.
- **v3.1**: Q07b + Q07c Hobs per-group (Engine H / Mérot), merged into a new
  Q04 panel p10 showing three-line HoverE (Hom1/Het/Hom2).

## Bug fixes included in this bundle (vs the HPC-applied version)

- **Bug 3 (QC_OUT default path)** — `00_config.sh` default patched to
  `${BASE}/inversion_localpca_v7/phase_qc_shelf_results` (was
  `${BASE}/inversion_modules/phase_qc_shelf/results` which caused outputs
  to land outside the toolkit tree).
- **Bug 1 (install_update_2026-04-20.sh LIVE_MODULE path)** — patched to
  `${PWD}/inversion-popgen-toolkit/inversion_modules/phase_qc_shelf` so the
  installer can find the live module when run from the project base dir.

The `config.local.sh` that was in the original bundle has been renamed to
`config.local.sh.example` — your real `config.local.sh` stays on HPC and in
`.gitignore`, so your personal paths don't go to GitHub.

## Staging layout

```
staging/
├── LOCAL_INSTALL_README.md              (this file)
├── inversion_modules/
│   └── phase_qc_shelf/                  (45 files — drop into your local clone)
├── tools/
│   ├── build_inventory.py
│   ├── standardize_batch.py
│   └── install_scripts/                 (for future HPC installs)
│       ├── install_update_2026-04-20.sh (bug 1 patched)
│       └── fix_output_paths_2026-04-20.sh
└── docs/
    ├── ENGINES.md
    ├── STANDARDIZE_FOLDER.md
    ├── STANDARDIZE_SCRIPT.md
    └── phase_qc_shelf_CHANGELOG.md
```

## Installation steps (your local machine)

### 1. Back up your current local clone (belt + braces)

```bash
cd /path/to/parent/of/your/clone
cp -r inversion-popgen-toolkit inversion-popgen-toolkit.backup_$(date +%Y%m%d)
```

### 2. Drop the staging contents into your clone

Unpack the tarball, then rsync the contents of `staging/` into your local
`inversion-popgen-toolkit/` directory. Use `--ignore-existing=false` so it
overwrites existing files (we want that — the whole point is to update them).

```bash
# Unpack the bundle somewhere temporary
tar xzf phase_qc_shelf_v3.1_local.tar.gz -C /tmp/

# Drop into your clone (rsync preserves structure, copies everything)
rsync -av --progress /tmp/staging/ /path/to/inversion-popgen-toolkit/

# OR with plain cp:
cp -rv /tmp/staging/inversion_modules/phase_qc_shelf/* \
       /path/to/inversion-popgen-toolkit/inversion_modules/phase_qc_shelf/
cp -rv /tmp/staging/tools/* \
       /path/to/inversion-popgen-toolkit/tools/   # mkdir first if absent
cp -rv /tmp/staging/docs/* \
       /path/to/inversion-popgen-toolkit/docs/
```

### 3. Check what git sees

```bash
cd /path/to/inversion-popgen-toolkit
git status
```

Expected output: a long list of `modified:` (existing files updated) and
`new file:` (Q09, Q09b, Q07b, Q07c, Q10, run_chrom.sh, etc.). Nothing
deleted. Nothing unexpected.

### 4. Sanity-check syntax on the shell scripts

```bash
cd /path/to/inversion-popgen-toolkit/inversion_modules/phase_qc_shelf
for f in *.sh; do bash -n "$f" && echo "OK  $f" || echo "FAIL $f"; done
```

All should print `OK`. I've already run this on the staged files; this is
just a paranoia check after the copy.

### 5. Commit

```bash
cd /path/to/inversion-popgen-toolkit
git add -A
git status   # take one last look
git commit -m "phase_qc_shelf v3.1: Hobs (Engine H) + registry wiring + v2.1 features

- v2.1: Q09 gap characterization, no-data strips in Q04, two-layer ideogram
  stipple, Q09b LD-structure check, run_chrom.sh per-chrom driver,
  standardized logging helpers in 00_config.sh
- v3.0: Q10 registry bridge, 28-chrom SLURM array, shelf_coords.tsv seed
  (LG28 15-18 Mb)
- v3.1: Q07b/c Hobs per karyotype group (Mérot), Q04 panel p10 (HoverE lines)

Fixes QC_OUT default path (bug 3) to inversion_localpca_v7/phase_qc_shelf_results.
config.local.sh kept gitignored; renamed bundle copy to .example template.

Tested on LANTA for v2.1 (LG28 inversion confirmed, figures produced).
v3.0 + v3.1 shipped untested from this session — will validate when HPC returns."
```

### 6. Push

```bash
git push origin main   # or whatever branch you use
```

### 7. On HPC (when it comes back ~2026-04-27)

```bash
cd /scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion-popgen-toolkit
git pull

# Your config.local.sh on HPC is untouched (it's gitignored). Good.
# Verify the QC_OUT path is correct now:
grep QC_OUT inversion_modules/phase_qc_shelf/00_config.sh
# Should show: : "${QC_OUT:=${BASE}/inversion_localpca_v7/phase_qc_shelf_results}"

# Smoke test on LG28:
cd inversion_modules/phase_qc_shelf
bash STEP_Q01_snp_density.sh C_gar_LG28 2>&1 | head -40

# Full pipeline on LG28 with new Q07b/c + Q10 registry:
SHELF_START_MB=15 SHELF_END_MB=18 BP1_MB=15.115 BP2_MB=18.005 \
  bash run_chrom.sh C_gar_LG28

# Check registry population:
bash scripts/registry_query.sh describe C_gar_LG28_15000000_18000000
```

## If anything goes wrong

- Restore from backup: `rm -rf inversion-popgen-toolkit && mv inversion-popgen-toolkit.backup_YYYYMMDD inversion-popgen-toolkit`
- Restore from git: `git reset --hard HEAD~1` (before push), or `git revert` (after push)
- Individual file rollback: `git checkout HEAD~1 -- path/to/file`

## Known caveats

- **v3.1 (Hobs / Engine H) is untested end-to-end on HPC.** It needs:
  - The patched ANGSD fork built at `$BASE/angsd_fixed_HWE/angsd`
  - The `hobs_windower` C binary built at
    `$BASE/unified_ancestry/engines/hobs_hwe/scripts/hobs_windower`
  If these aren't present, Q07b/c will fail with a clear error message and the
  pipeline's other steps still work.
- **v3.0 registry** requires the registry APIs to exist in
  `inversion-popgen-toolkit/utils/` (interval_registry, sample_registry,
  results_registry, evidence_registry). Those should be in your existing
  toolkit; Q10 will fail gracefully with a path-not-found message if not.
- **The HoverE panel p10 in Q04** only renders when the Q07b/c output file
  `hobs_merged.<CHR>.tsv.gz` exists. Q04 falls back to the 9-panel layout
  gracefully when absent.
