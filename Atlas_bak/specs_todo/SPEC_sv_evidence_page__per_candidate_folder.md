# SPEC patch — per-candidate folder layout (from Quentin's step-4 chat turn)

**Date:** 2026-05-04
**Status:** documenting Quentin's preferred local-first folder structure
**Patches:** `specs_todo/SPEC_sv_evidence_page.md` §2 (data layout)

## Decision

Per-candidate folders, one folder per candidate, with fixed filenames inside.

```
data/
└── LG28/                                      ← one folder per chromosome (optional grouping)
    └── candidates/
        ├── INV_LG28_001/                       ← one folder per candidate
        │   ├── sv_genotype_counts.json         ← step 1-2-3 main layer
        │   ├── sv_evidence_combinations.json   ← step 5 layer (Quentin's UpSet redirect)
        │   ├── sv_support_by_sample.json       ← step 5+ S7 read-evidence (per-sample × per-evidence-type 0/1 matrix)
        │   └── (future) bam_evidence.json      ← S7 BAM-level cluster evidence
        ├── INV_LG28_002/
        │   └── ...
        └── INV_LG28_003/
            └── ...
```

## Why this is right for the local-first workflow

- **One mental unit per candidate.** The user thinks "I want to look at INV_LG28_003" — drag-drop the folder, atlas walks it, populates whatever layers it finds.
- **Schema versioning per file.** Each JSON carries its own `format_version`; the atlas just reads. No coordination across files.
- **Partial-state friendly.** If only `sv_genotype_counts.json` exists, the atlas renders steps 1-4. If `sv_evidence_combinations.json` shows up later, the UpSet panel populates.
- **Candidate scope clean.** All artefacts for one candidate live together. Easy to share, archive, version, or move.
- **Supports Quentin's "download from HPC, run locally" model.** A reviewer can be sent one folder for one candidate; the atlas does the rest.

## Atlas-side behaviour

### Reading (already partly works)

Step 1's drag-drop already accepts any single JSON with `format_version: sv_genotype_counts_v1`. Step 5 will need:

1. **Folder drop** — when the user drops a folder (not a file), iterate `dataTransfer.items`, walk the folder via `webkitGetAsEntry()`, identify each `*.json` by reading the first ~8 KB and checking `format_version`.
2. **Layer routing** — based on `format_version`, populate the matching `_state.*` slot:
   - `sv_genotype_counts_v1` → `_state.layer`
   - `sv_evidence_combinations_v1` → `_state.combinationsLayer` (new)
   - `sv_support_by_sample_v1` → `_state.supportLayer` (new, future S7)
3. **Candidate id resolution** — if any layer carries `candidate_id`, all layers in the folder must agree (or the folder name overrides). Mismatched ids = error.
4. **Backwards-compatible single-file path** — current behaviour preserved. Drag-drop one JSON, get one layer populated.

### Producer side (HPC)

A writer that builds the per-candidate folder. For each candidate:

```python
import json, os
def write_candidate_dir(cand_id, out_root, **layers):
    """layers maps short-name → dict; writes each as <short_name>.json."""
    cand_dir = os.path.join(out_root, cand_id)
    os.makedirs(cand_dir, exist_ok=True)
    for name, payload in layers.items():
        with open(os.path.join(cand_dir, f"{name}.json"), "w") as f:
            json.dump(payload, f, separators=(",", ":"))
```

`STEP_SV_GT_AGG` (step 7) and `STEP_SV_EVID_COMB` (the future step-5 producer) both write into the same per-candidate folder. They're independent — either can ship without the other.

## Implementation timing

- **Step 5** is the right place to add the folder-walk drag-drop, because that's when there's a second layer to coordinate.
- **Step 1-4** keep the current single-file `fetch('json/sv_genotype_counts/<cid>.json')` path (still works against `atlas_server.py` if a user wants a server-backed setup).
- **Step 7** writes both layers into one per-candidate folder.

## Backwards compatibility

The current `json/sv_genotype_counts/<cid>.json` flat layout keeps working. Step 5 adds folder-walk as an additional path, not a replacement. Users who already have flat layouts just get fewer features (UpSet panel stays empty until they have the second layer too).

## Open question

Should the per-candidate folder also carry a small `candidate.json` manifest with metadata (chrom, boundaries, span, tier, validation status)? It would be redundant with the embedded fields in each layer, but it answers "what is this folder?" without parsing every JSON. **Recommendation: yes, lightweight ~20-line manifest.** Not blocking for step 5; can add in step 7's writer.
