# Atlas merge bundle — phase_8 + phase_X + turn129 review docs

**For:** Quentin Andres
**Date:** 2026-05-04
**Source chats:** phase_8 (comparative-fragility) + phase_X (U-shape) + turn129 review bundle

This tarball drops directly into your existing `Atlas/` folder. Nothing
overwrites; everything is additive. Two new modules + three docs files are
introduced.

## What's in the bundle

```
Atlas/
├── _scripts/
│   ├── phase_8_comparative_breakpoint_fragility/   ← NEW (this bundle)
│   └── phase_X_ushape_evolution/                   ← NEW (this bundle)
│
├── handoffs_to_implement/
│   ├── HANDOFF_2026-05-04_phase_8_phase_X.md       ← NEW (this file)
│   └── turn129_bundle/                             ← NEW (preserved as-is)
│       ├── README.md
│       ├── plans/
│       │   ├── 00_master_plan.md
│       │   └── 01_priority_triage.md
│       ├── patches/   (P1.1 .. P6.1, 13 files)
│       └── specs/     (S1 .. S5, 5 files)
│
├── specs_todo/
│   ├── from_turn129/                               ← NEW (S6, S7 from turn129)
│   │   ├── README.md
│   │   ├── S6_dosage_heatmap_streaming_viewer.md
│   │   └── S7_karyotype_breakpoint_internal_evidence.md
│   └── SPEC_comparative_te_breakpoint_fragility__implementation_status.md  ← NEW
│
└── json/
    ├── comparative_breakpoint_fragility/           ← NEW (output landing dir)
    └── ushape_evolution/                           ← NEW (output landing dir)
```

**Not touched / not overwritten:**
- `Atlas/Inversion_atlas.html` (and `.bak` siblings) — left alone, per your
  instruction to skip the older atlas HTML.
- `Atlas/specs_todo/SPEC_comparative_te_breakpoint_fragility.md` — the
  existing forward-looking spec is left intact. A new sibling
  `*__implementation_status.md` records that the code now exists and points
  to it.
- `Atlas/server_turn1/popstats_server.py` — the U-shape endpoint is shipped
  as a sibling module (`server/ushape_endpoint.py` inside phase_X). To
  activate, paste the four-line snippet from `phase_X/README.md` into the
  server file. Doing it as a paste means the server file stays canonical and
  the endpoint module stays version-controlled separately.

## Module 1 — phase_8_comparative_breakpoint_fragility

**Purpose.** Comparative TE / repeat-density layer for Gar inversion
breakpoints. For every Gar candidate, asks (a) is the focal breakpoint
TE-enriched vs chromosome and local background; (b) when synteny mappings
exist to Mac / Pangasius / Silurus / Ictalurus / etc, do the homologous
boundaries also sit on TE-rich architecture. Outputs a JSON layer
(`comparative_breakpoint_fragility_v0.1`) and a flat TSV summary.

**Status.** Complete v0.1 — 8 step scripts, common utilities, R plotting,
end-to-end runner. Tolerant of messy TE folder structures (EDTA, GFF3, RM,
BED). Vocabulary discipline encoded in the classifier and the README §6.

**Where it goes:** `Atlas/_scripts/phase_8_comparative_breakpoint_fragility/`.
**How to run:**

```bash
cd Atlas/_scripts/phase_8_comparative_breakpoint_fragility/
# fill in config/species_manifest.tsv, config/candidate_breakpoints.tsv
# optionally config/synteny_mapping.tsv
bash run_all.sh --te_root /path/to/messy/te_outputs
```

**Atlas integration:** outputs at `Atlas/json/comparative_breakpoint_fragility/`.
The atlas's existing `_renderRepeatDensityPanel` (line ~15770 per the spec
notes) needs a small addition to read the new schema's `comparative_species`
block. That is **not** in this bundle — see the spec at
`Atlas/specs_todo/SPEC_comparative_te_breakpoint_fragility.md` for the
contract.

## Module 2 — phase_X_ushape_evolution

**Purpose.** Per-candidate evolutionary-shape classification. For each
inversion, computes windowed dXY, FST, π_HOMO1, π_HOMO2, AFD across the
inversion plus flanks; classifies the resulting profile shape; outputs
`ushape_evolution_v1.json` plus per-candidate plots.

**Status.** Complete v0.1 — 7 CLI step scripts, 3 sourceable R libraries,
1 FastAPI server endpoint module, 1 atlas-side JS renderer, 1 plain-language
class-rules doc, 1 YAML config, end-to-end runner.

**Architecture.** Server-primary (popstats `/api/ushape/candidate`), CLI
batch as fallback. Both paths share the same R libs in
`phase_X_ushape_evolution/R/`. The server endpoint reuses `region_popstats`
for the actual FST/dXY/π compute; only the post-processing (zoning,
scoring, classification, JSON shaping) is new.

**Where it goes:** `Atlas/_scripts/phase_X_ushape_evolution/`.
**How to activate the server endpoint:** paste the four-line snippet from
`phase_X_ushape_evolution/README.md` into `Atlas/server_turn1/popstats_server.py`
near the existing endpoints. Add `ushape_r_dir:` to
`popstats_server.config.yaml`. Restart the server.

**Atlas integration:** Drop `phase_X_ushape_evolution/atlas_js/ushape_renderer.js`
into `Atlas/js/`, add the script tag, add the five empty container `<div>`s
(IDs listed in the README) wherever the surfaces should appear. The
renderer no-ops gracefully if a container is missing, so this can be done
incrementally.

**Important biology note.** The classifier uses an explicit **age gate**:
candidates that fail an oldness threshold are labelled `young_weak_divergence`
without a shape verdict. This deviates from the original prompt (which had
`oldness_score` and `youngness_score` as parallel axes) and matches the
Berdan/Guerrero figure 2 model where shape is only informative on old
inversions. If you'd rather have shape called on young candidates too, flip
`oldness_min` low in `config_ushape.yaml`.

## What's NOT shipped here

- **Atlas HTML edits.** No `Inversion_atlas.html` modifications anywhere.
  Both modules expose container-IDs that the atlas can opt into when ready.
- **Server modifications to `popstats_server.py`.** The U-shape endpoint
  ships as a sibling module; activating it is a four-line paste documented
  in the phase_X README.
- **Real LANTA data.** No genome-wide window stats path is hardcoded. U03
  falls back to flank-only when missing.
- **The 24/28 Clair3 chromosome restriction.** Both modules iterate over
  whatever chromosomes the inputs supply — pass only the 24 you have, the
  rest cleanly drop.
- **Patches P1–P6 from turn129.** Those are still anchor-based markdown
  patch docs awaiting your decision-point sign-off (stage naming, catalogue
  location, SV merge name, turn-128 pill removal). They're preserved
  unchanged at `Atlas/handoffs_to_implement/turn129_bundle/patches/`.

## Reading order for whoever continues

1. **This file** (the one you're reading)
2. `phase_8_comparative_breakpoint_fragility/README.md` — what the comparative
   fragility module computes, vocabulary discipline, JSON schema, classification
   labels A–G
3. `phase_X_ushape_evolution/README.md` — U-shape architecture, server vs
   CLI modes, atlas integration sketch
4. `phase_X_ushape_evolution/docs/ushape_evolution_class_rules.md` — the safe
   manuscript phrasings (Reviewer 2 protection)
5. `Atlas/handoffs_to_implement/turn129_bundle/README.md` — the planning
   review bundle that motivated phase_X (and queued P1–P6 patches)
6. `Atlas/specs_todo/SPEC_comparative_te_breakpoint_fragility.md` — the
   forward-looking spec phase_8 implements

## Cross-module relationships

Both modules attach to the **candidate page** but answer different
questions:

| module                             | question                                          | placement                                     |
|------------------------------------|---------------------------------------------------|-----------------------------------------------|
| `phase_8_comparative_breakpoint_fragility` | is the breakpoint architecture fragile?    | repeat-density area + Boundaries page         |
| `phase_X_ushape_evolution`         | what does the divergence profile look like?      | Candidate page only ("Evolutionary shape")    |

They never share fields. The U-shape JSON is deliberately TE-free and the
fragility JSON is deliberately popgen-free. If you want a panel that
combines them on a single screen, do it at the atlas-rendering layer (CSS
adjacency), not at the data-schema layer.

## Quick provenance

- **phase_8** scaffold matches the suggested layout in
  `Atlas/specs_todo/SPEC_comparative_te_breakpoint_fragility.md` 1:1
  (lines 57–76 of that spec). Architecture classes A–G and age tags
  (YOUNG_POP / OLD_BP_YOUNG_INV / etc) match the prompt verbatim.
- **phase_X** rule logic deviates from the prompt by (a) gating on age
  before assigning shape (Berdan/Guerrero figure 2 dictates this), (b)
  picking Hudson FST as the default (Weir-Cockerham needs per-individual
  genotypes the server doesn't expose), (c) excluding HET samples from the
  compute (dXY/FST is between-arrangement only). All three deviations are
  flagged in the conversation that produced this bundle and in the README.

## License / attribution

Internal to `MS_Inversions_North_african_catfish`. PI: Prapansak Srisapoome.
Module authorship: Quentin Andres. Pipeline scaffolding generated with
assistance from Claude.
