# Atlas — interactive viewer + curation layer

This repo holds **multiple Atlases**, one per scientific workflow.
Each Atlas is a set of single-file HTML viewers backed by JS modules.

```
inversion_discovery.html      ← inversion calling (PCA scrubbers, candidate focus)
inversion_review.html         ← refinement (auto-promote, accept/reject)
inversion_catalogue.html      ← synthesis (manuscript output)
inversion_comparative.html    ← cross-species

diversity_*.html              ← future
population_*.html             ← future
assembly_*.html               ← future
```

For the full naming convention and data folder layout, see
[`data/README.md`](data/README.md).

---

## Quick start

### Run an Atlas in dev mode

```bash
cd Atlas/
python3 -m http.server 8000
# open http://localhost:8000/inversion_review.html
```

ES-module imports require a server (browsers refuse `file://` module loads
in most settings). Any static-file server works; Python's built-in is
the lowest-friction.

### Build a single-file distributable

```bash
python3 build/flatten.py inversion_review.html
# writes dist/inversion_review_flat.html
```

The flattened HTML inlines all JS modules into one `<script type="module">`
block. Drag-drop into any browser, no server needed. Use this for
sharing with collaborators or committing release builds.

### Run smoke tests

```bash
node tests/test_modular_smoke.js
```

Verifies the modular pattern (imports, exports, paths, JSON/TSV
round-trip, `flatten.py` output validity).

---

## Repo layout

```
Atlas/
├── README.md                          ← you are here
├── data/                              ← canonical data layout (see data/README.md)
│   ├── precomp/<chrom>/               ← R pipeline output, per chromosome
│   ├── cohort/                        ← cohort-level metadata
│   ├── candidates/<cid>/              ← per-candidate evidence
│   ├── comparative/                   ← cross-species
│   └── review/<workflow>/             ← atlas-side curation (writable)
│
├── shared/                            ← cross-workflow JS modules
│   └── state_io.js                    ← canonical I/O, the only file that
│                                        encodes the data folder layout
│
├── inversion_discovery.html           ← discovery phase Atlas (TBD)
├── inversion_discovery/
│   ├── main.js
│   ├── panels/
│   └── pages/
│
├── inversion_review.html              ← review phase Atlas
├── inversion_review/
│   ├── main.js
│   └── panels/                        ← (TBD: top intervals, graph, sample lines, L3 contingency, side actions)
│
├── inversion_catalogue.html           ← TBD
├── inversion_catalogue/
├── inversion_comparative.html         ← TBD
├── inversion_comparative/
│
├── build/
│   └── flatten.py                     ← inline ES modules into one HTML
│
├── tests/
│   └── test_modular_smoke.js          ← validates the modular pattern
│
└── dist/                              ← flatten.py output (gitignored)
```

---

## Conventions (locked)

1. **HTML naming**: `<workflow>_<phase>.html`. Phases: `discovery` /
   `review` / `catalogue` / `comparative` (others workflow-specific).
2. **JS module folders**: same name as the HTML they serve
   (`inversion_review/` for `inversion_review.html`).
3. **Cross-workflow modules**: `shared/`. Anything that imports from
   `shared/` must work in any Atlas.
4. **Data layout**: read-only inputs are `precomp/cohort/candidates/comparative/`;
   writes go to `review/<workflow>/`. Atlases never write outside `review/`.
5. **Imports**: relative paths only (`'../shared/state_io.js'`). No
   bundlers, no npm, no node_modules. The flatten.py inliner depends
   on this constraint.
6. **No build step at dev time**. You edit a JS file, refresh the
   browser, see the change. `flatten.py` runs only when producing
   distribution HTMLs.

---

## What lives elsewhere

- **R pipeline** (the upstream code that produces `data/precomp/`,
  `data/cohort/`, `data/candidates/`): in `inversion_codebase_v8.5/`
  and on LANTA. Not part of this repo.
- **Manuscript drafts**: separate Drive folder.
- **HPC scratch / intermediate files**: never committed.

The Atlas is the **interactive viewer + curation** layer. The R
pipeline produces inputs; the Atlas writes curation outputs to
`data/review/<workflow>/`. That's the contract.
