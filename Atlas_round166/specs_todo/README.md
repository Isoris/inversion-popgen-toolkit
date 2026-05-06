# specs_todo/ — atlas + pipeline specs in the build queue

Specs that are approved for build but not yet fully shipped. Each spec
is self-contained, written so a future Claude session (or Quentin) can
pick it up cold.

## Lifecycle

```
specs_new_turn131/   ← drafted, awaiting review for build
       ↓
specs_todo/          ← approved, in build queue (or partially built)
       ↓
specs_done/          ← fully shipped
```

When a spec ships fully, move it to `../specs_done/` and update its
top-of-file `**Status**:` line to `SHIPPED in turn N (date)` with a
pointer to the test file.

When a spec ships *one slice* but other slices remain, leave it here
and update the `**Status**:` line to capture both — what shipped,
what's deferred.

## Adding a new spec

- One Markdown file, named `SPEC_<short_id>.md`
- Top of the file: `**Status**: forward-looking spec. Not implemented.`
- Top of the file: `**Reading order**: this spec → ...`
- Vocabulary discipline section if the spec touches manuscript-grade claims
- Tests section listing concrete test names that would gate "done"
- Open questions / deferred list at the bottom

## Inventory (2026-05-05)

### Forward-looking — no slices shipped yet

| Spec | One-line | Blocking dependency |
|---|---|---|
| `SPEC_focal_vs_background_widget.md` | Spalax-style focal-vs-bg widget for page 16 + page 11 | aggregated JSON layer optional for manuscript-locked P |
| `SPEC_ld_decay_overlay.md` | Spalax Fig. S62 — three-trace LD decay overlay | precomputed `ld_decay_genome_wide.json` optional |
| `SPEC_ncrna_density_layer.md` | rRNA/tRNA/ncRNA density panel on page 11 + flanking chips | aggregator script `aggregate_ncrna_density.R` |
| `SPEC_cross_species_dotplot.md` | Mummerplot-style dot plot with hover-enlarge + click-to-pin | mashmap pre-compute (~10s, one-shot) |
| `SPEC_inversion_age_atlas_surface.md` | Two-task atlas surface for inversion age (Task A between-inversion + Task B within-arrangement) | new HPC emitter `STEP_C01f_d_emit_age_json.py` (~80 lines) |
| `SPEC_busco_anchors.md` | BUSCO anchor layer (multispecies) | none |
| `SPEC_phylogenetic_tree_integration.md` | Tree/age polarization layer (multispecies) | tree exists per Quentin |
| `SPEC_hypothesis_registry_and_multispecies.md` | Registry data model + classification UI | none |
| `SPEC_l2_envelope_live_split.md` | Live split refinement on L2 envelopes | none |
| `SPEC_distant_band_concordance_fish_trajectory.md` | Fish-trajectory concordance pipeline | none |
| `SPEC_sliding_window_inheritance.md` | Sliding-window inheritance compute | none |
| `SPEC_band_track_extraction_and_l3_single_band_rows.md` | Per-band rows under K×K tables in L3 panel + automated band-track extraction (multi-segment / A-B-A aware) | none |
| `QR03_SPEC.md` | Q-regime detector (STEP_QR03) | STEP_QR02 outputs |

### Partially shipped — at least one slice landed

| Spec | What's shipped | What's pending |
|---|---|---|
| `SPEC_l3_het_dosage_coloring.md` | slice 1 (turn 129) | slice 2 (tracked-samples always-on, per-sample-lines extension) |
| `SPEC_lines_panel_candidate_bands.md` | slice 1 (turn 141) | slice 2 (see §9 in spec) |
| `SPEC_age_origin_panel.md` | main panel (turn 119) | §8 Q5 chip, §9 card-flip |
| `SPEC_l2_sweep_inheritance.md` | auto-promote (turn 133) + inspector (turn 134) | dispatcher refactor making inheritance compute consume L2 envelopes directly |
| `SPEC_review_surfaces_auto_and_lineages.md` | slice 0 infrastructure (turn 130) | lineage-bloc grouping tab in G-panel + review queue (now unblocked) |

### Decision-pending / meta — keep here, may never become buildable

| Spec | Why it's here |
|---|---|
| `SPEC_inheritance_unification.md` | Refactor decision document. Pulled when L2-sweep + sliding-window + fish-trajectory have all shipped *and* been used on real data. |
| `SPEC_inheritance_diagnostic_protocol.md` | Diagnostic protocol, not a code spec. Could move to `docs/`. |
| `SPEC_OVERVIEW_multispecies_architecture.md` | Integration map across the multispecies specs. Reference doc. |

## Subdirectories

- `from_turn129/` — two specs (S6 dosage streaming viewer, S7 karyotype
  breakpoint internal evidence) carried over from the turn-129 review
  bundle. See its own README.
- `later/` — `SPEC_xpehh_track.md` and similar specs that need data
  Quentin doesn't have yet (e.g. phased VCF for the 226-sample cohort).
- `_mockups/` — design mockups referenced by specs (currently the
  three SV evidence page PNGs).
