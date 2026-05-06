# Atlas/_scripts/

Infrastructure for producing and maintaining the JSON files in `../json/`.
You should not need to touch anything here for normal atlas use.

## What's in here

### Cluster-side R scripts (run on LANTA, produce JSONs)

```
STEP_C04c_ghsl_local_pca.R       local PCA on STEP_C04 v6 GHSL matrices
                                  produces *.ghsl_v6_localpca.rds per chrom
STEP_C04d_ghsl_d17_wrapper.R     wraps STEP_D17 cross-block detector for GHSL
                                  produces 4 TSVs per chrom (L1/L2 envelopes + boundaries)
export_ghsl_to_json_v3.R          consolidates STEP_C04+C04c+C04d outputs
                                  produces <chr>_phase2_ghsl.json
STEP_TR_C_theta_d17_wrapper.R    wraps STEP_D17 for θπ
                                  produces 4 TSVs per chrom
STEP_TR_D_augment_theta_json.R    adds theta_d17_envelopes to existing θπ JSON
```

### SLURM launchers (one task per chromosome, 28-task arrays)

```
run_ghsl_array.sh                 chains C04c + C04d + export_ghsl_to_json_v3
run_theta_array.sh                chains TR_A + TR_B + TR_C + TR_D
```

### File-naming bridge (cluster names ↔ atlas names)

```
sync_jsons.sh                     copy from cluster paths to ../json/
                                  with the short LG##_<stream>.json convention
```

### Documentation

```
SESSION_MANIFEST.md               READ THIS FIRST — what's done, what's pending
ARCHITECTURE_DECISIONS.md         12 ADRs — why decisions were made
JSON_CONTRACT.md                  per-layer field specs for all precomp JSONs
RUNBOOK_produce_phase2_jsons.md   end-to-end walkthrough (LG28 dry-run + full array)
```

### Archive

```
_archive/                         atlas.html patchers, kept for reference.
                                   These are one-shot mutators: once applied to
                                   atlas.html, they're history. Re-apply only
                                   if you need to rebuild atlas.html from a
                                   fresh pca_scrubber_v4 base.
```

## Daily workflow

You will not run anything in here daily. The cluster pipeline runs once
per chromosome batch. The atlas patchers run once per atlas rebuild.

## When something here matters

- New chromosome batch produced on cluster → run `sync_jsons.sh`
- Need to add a new chromosome to `../json/` → cluster pipeline (see runbook)
- atlas.html broken or you want to rebuild from upstream → re-apply patchers
  from `_archive/` in order: scaffold → color_mode → panels → merge_fix

## What does NOT go here

- The atlas itself (`atlas.html` lives one level up)
- Your work products (candidate TSVs, screenshots — those go in `../saved/`)
- The JSON files (those go in `../json/`)
- Test scripts — testing happens during atlas development, not in production use.
  Tests that come with this session are in `_archive/tests/`.
