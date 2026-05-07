# Atlas — canonical layout

This doc is the single source of truth for **where files go** in this
repository. Two locked pieces:

1. **Naming**: Atlases are `<workflow>_<phase>.html`
2. **Data folder**: five subfolders, separated by **who writes them**

If you're adding a new file and you're not sure where it goes, this
document tells you. If it doesn't, the doc is wrong — fix it.

---

## 1. Naming convention for Atlases

`<workflow>_<phase>.html`

| Workflow | Discovery | Review | Catalogue / Synthesis | Compare |
|---|---|---|---|---|
| **inversion** | `inversion_discovery.html` | `inversion_review.html` | `inversion_catalogue.html` | `inversion_comparative.html` |
| **diversity** | `diversity_discovery.html` | `diversity_review.html` | `diversity_catalogue.html` | `diversity_comparative.html` |
| **population** | `population_discovery.html` | `population_review.html` | `population_catalogue.html` | (n/a) |
| **assembly** | `assembly_qc.html` | `assembly_review.html` | `assembly_catalogue.html` | `assembly_synteny.html` |

Each Atlas has a sibling folder with the same name for its JS modules:

```
inversion_discovery.html
inversion_discovery/
    main.js
    panels/
    pages/
```

`shared/` is at the repo root — modules in there (Hungarian solver,
contingency-table builder, color ramps, JSON I/O) are cross-workflow.

---

## 2. The data folder layout — locked

```
data/
├── precomp/        # R pipeline output, per-chromosome.    READ-ONLY for atlases.
├── cohort/         # Cohort-level, chromosome-independent. READ-ONLY for atlases.
├── candidates/     # Per-candidate evidence stacks.        READ-ONLY for atlases.
├── comparative/    # Cross-species / cross-genome.         READ-ONLY for atlases.
└── review/         # Atlas-side curation state.            ATLASES WRITE HERE.
```

**The rule that keeps this sane**: you can `rm -rf data/precomp/<chrom>/`
and re-run the upstream R pipeline without losing any curation work.
Anything an atlas writes goes under `review/`, full stop.

### `precomp/` — per-chromosome R pipeline output

One subfolder per chromosome. Files inside named
`<chrom>.<layer>.<ext>` so a flat search across chromosomes works:

```
find data/precomp -name '*.band_tracks.*'
```

```
data/precomp/
├── LG28/
│   ├── LG28.json                              # main scrubber data: windows, envelopes,
│   │                                          # tracks, samples, cluster labels, ancestry,
│   │                                          # ghsl, theta_pi, dosage_chunks index
│   ├── LG28.repeat_density.scrubber_windows.json
│   │
│   │   # SPEC BLOCK 1 R-module outputs (per-chromosome band tracks):
│   ├── LG28.band_nodes.tsv
│   ├── LG28.band_edges.tsv
│   ├── LG28.transition_events.tsv
│   ├── LG28.band_trajectories.tsv
│   ├── LG28.het_band_backbones.tsv
│   ├── LG28.candidate_track_proposals.tsv
│   ├── LG28.candidate_tracks.json
│   └── LG28.manual_review_queue.tsv
│
├── LG12/
│   └── ...
└── LG07/
    └── ...
```

An atlas detects "chromosome X exists" by the presence of
`precomp/<chrom>/<chrom>.json`. All other files are optional layers; the
atlas degrades gracefully when a layer is missing (per-layer empty-state
UI; nothing else breaks).

### `cohort/` — cohort-level metadata

The 226-fish cohort doesn't change per chromosome. Loaded once at
atlas startup, reused everywhere.

```
data/cohort/
├── relatedness.json                # KING/Manichaikul kinship; hub_id_1st/2nd
├── cohort_diversity_v1.json        # GenomeScope baseline H, Hobs/Hexp
├── sample_froh.json                # per-sample F_ROH, ROH intervals
├── sample_manifest.tsv             # 226 rows: sample_id, family, sex, hatchery
├── sample_groups.tsv               # named groups (e.g. NAToRA-pruned 81 unrelated)
└── natora_pruned.tsv               # the unrelated sample set
```

### `candidates/` — per-candidate evidence

One subfolder per **promoted** candidate. Naming convention:
`cand_<chrom>_<approx_Mb>` (matches the pre-existing
`sv_genotype_counts/cand_LG28_15Mb.json` style).

```
data/candidates/
├── cand_LG28_15Mb/
│   ├── sv_genotype_counts.json
│   ├── boundaries_refined.json
│   ├── gene_cargo.json
│   ├── marker_primers.json
│   ├── breeding_readiness_card.json
│   └── final_classification.json   # 14-axis tier classification
└── cand_LG12_8Mb/
    └── ...
```

Atlases read from here but never write here. Candidates are produced by
the upstream R pipeline and committed to this folder. The user's
accept/reject/lock/override decisions go to `review/` instead.

### `comparative/` — cross-species

Different cohort, different question (cross-species karyotype comparison
for the manuscript). Held separate from `precomp/` because mixing them
breaks scripts that loop over chromosomes-of-the-target-cohort.

```
data/comparative/
├── cs_breakpoints_v1.json
├── synteny_multispecies_v1.json
├── phylo_tree_v1.json
├── comparative_breakpoint_fragility/
│   └── (per-comparison files)
└── te_fragility_v1.json
```

### `review/` — atlas-side curation, the only writable folder

Subdivided by workflow so different Atlases don't collide on filenames.

```
data/review/
├── inversion/
│   ├── manual_overrides.json           # force-merge / force-cut / lock_candidate / allow_cross_L2
│   ├── candidate_review_decisions.json # accept / reject / split / merge / mark_complex
│   ├── locked_karyotype_groups.json    # H1/H1, H1/H2, H2/H2 lockings
│   ├── confirmed_candidates.json       # the curated final list
│   └── sessions/
│       └── session_2026-05-05_HHMM.json    # full-state autosaves
├── diversity/
│   └── ... (future)
├── population/
│   └── ... (future)
└── assembly/
    └── ... (future)
```

Files under `review/<workflow>/` are JSON. They're the only
non-deterministic data in the repo (every other file is reproducible
from the upstream pipeline). Treat them as you would lab notebook
entries — back them up before any destructive operation.

---

## 3. Layer detection (atlas-side)

Atlases use `shared/state_io.js`'s `detectLayers(precompDir)` to know
which layers exist for a chromosome. The function reads the directory
listing, matches against the canonical filename pattern, and returns
`Set<string>` of layer names. UIs key off that set to enable/disable
panels. **No layer name is required**; the atlas always degrades
gracefully when a layer is absent.

The canonical layer-name list is in
`shared/state_io.js` (constant `KNOWN_LAYERS`). Adding a new layer = add
its name + filename pattern there + a one-line entry in this README.

---

## 4. Adding a new workflow (e.g. diversity)

1. Pick a workflow name (single word, lowercase, no underscores).
2. Add new HTML files: `<workflow>_discovery.html`,
   `<workflow>_review.html`, `<workflow>_catalogue.html`,
   `<workflow>_comparative.html` (skip phases you don't need).
3. Create matching `<workflow>_<phase>/` JS module folder for each.
4. Inside `data/precomp/<chrom>/`, your workflow will have its own
   layer files (e.g. `LG28.diversity_panels.json`). Add their detection
   patterns to `shared/state_io.js`.
5. Inside `data/review/`, create `data/review/<workflow>/`.
6. Update this README's section 1 table with the new row.

That's the whole protocol. Cohort, candidates, comparative, and the
precomp folder structure are reused as-is.

---

## 5. What does NOT live here

- **R scripts / R packages** — those live in
  `inversion_codebase_v8.5/MODULE_*` and on LANTA, not in the Atlas
  repo. The Atlas consumes their outputs via the layout above.
- **Manuscript drafts** — separate repo / Drive folder.
- **Slurm submission scripts** — in `inversion_codebase_v8.5/RUN_*.sh`.
- **HPC-side scratch / intermediate files** — never committed; only
  the final precomp/cohort/candidates outputs cross over.

The Atlas repo is the **interactive viewer + curation** layer. The
upstream R pipeline produces inputs, the Atlas produces curated
review outputs, and that's the contract.
