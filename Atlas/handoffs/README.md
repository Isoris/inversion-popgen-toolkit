# Atlas handoffs — directory layout

**You-are-here pointer:** the canonical current-session report is
`../HANDOFF_FULL_BUNDLE_2026-05-05.md` at repo root.

This folder organises the per-turn handoff `.md` files so sessions
don't have to wade through 47 files at root. Layout reflects code
state, not chronology — each section corresponds to a concrete
slice of `Inversion_atlas.html`'s history.

## `current/` — turns 151 → 158 (this session)

10 files. Everything in here corresponds to features verified
present in the live `Inversion_atlas.html` (73,773 LOC) and tests
that pass at **2991 / 0**.

The two parallel turn-157 deliverables are renamed for clarity:
- `turn157A_vshape_tooltip.md` — V-shape hover tooltip (track A)
- `turn157B_dosage_bridge_load.md` — dosage shim load fix + page6/7
  try/catch (track B, parallel session)

Also note: `HANDOFF_2026-05-05_turn151_session_handoff.md` is the
broader session-entry handoff, distinct from
`turn151_slab_het_coloring.md` which documents only the het-coloring
slice shipped that turn.

## `archive/2026-05-05_turn128-150/` — turns 128 → 150

25 files documenting the prior session arc. Atlas was 54k → 72k LOC
across these turns. All features shipped; kept for archaeology.

Notable bundle-handoffs in this group:
- `HANDOFF_sv_evidence_complete.md` (was `HANDOFF.md`)
- `HANDOFF_2026-05-04_FINAL.md` — SV evidence + producers
- `HANDOFF_2026-05-05_turn150_FINAL.md` — slab UV cluster modes,
  closes out the slab-track work that turns 151 then continued

## `archive/pre_turn128/` — pre-turn-128 history

11 files + `turn129_review_bundle/` subdirectory. Dates: April 30
through May 4. Atlas was 54k–62k LOC at the time. Everything here
is shipped or superseded:
- `HANDOFF_inheritance_labeling_turn1-2p.md` (was `HANDOFF.md`)
- `HANDOFF_apr04.md`, `HANDOFF_apr30.md` — pre-real-data scaffolding
- `HANDOFF_2026-05-03_*` — turn 117–120 spec implementations
- `HANDOFF_2026-05-04_phase_8_phase_X.md` — comparative-fragility +
  U-shape merge
- `HANDOFF_after_turn127.md` — pre-128 entry point
- `HANDOFF_diversity_atlas.md`, `HANDOFF_genome_atlas.md` — sister
  atlas scaffolds
- `turn129_review_bundle/` — review package (P1.1 → P6.1 patches +
  S1 → S5 specs). All P-patches verified landed in current atlas
  (canonical server URL key, `/api/dosage/chunk`, SV evidence
  skeleton, etc.). S-specs absorbed into the SV-evidence track.

## Promotion / archival rules

- **Active handoffs** (during a session) live in `current/`.
- At session close, the FULL_BUNDLE handoff at root is updated and
  the per-turn `current/` files of the closing session are moved
  to a new `archive/<date>_turn<X-Y>/` directory.
- Don't delete handoffs. Archaeology beats cleanliness when a future
  session needs to reconstruct decision history.
