# HANDOFF — next chat: build the SV Evidence page

**Date written:** 2026-05-04
**For:** the next chat that will continue this work
**From:** Quentin Andres · MS_Inversions_North_african_catfish

---

## TL;DR for the next Claude

1. **Drag-drop the tarball** `Atlas_dragdrop.tar.gz` to access the full
   working tree at `/home/claude/Atlas/` (or wherever the user uploads
   it). Everything below assumes this tree is available.
2. **Read** `Atlas/specs_todo/SPEC_sv_evidence_page.md` end to end. The
   mockups it references live at `Atlas/specs_todo/_mockups/`.
3. **Build** the SV Evidence page following §7 *Implementation order*,
   step by step. Quentin works in many parallel chats — keep
   commits small and self-contained, one §7 step per coherent unit.
4. **Don't touch** `Atlas/Inversion_atlas.html` directly unless the
   spec explicitly requires it (§4.1 tab registration, §4.2 page
   mount). All page logic lives in the new file
   `Atlas/js/atlas_sv_evidence.js`.

## Context Quentin won't repeat

Quentin is finalising a manuscript ("MS_Inversions_North_african_catfish")
analyzing 226 pure *Clarias gariepinus* hatchery samples on the LANTA
HPC (account `lt200308`). The Inversion Atlas is the visualization tool
for the inversion catalogue. He works mobile-first, terse, French-native,
pushes back on errors precisely. The atlas tabs visible in the mockups
(`Overview / Karyotype / Tiers / Boundaries / Dosage / PCA / SV Evidence
/ Coherence / Polarity / Sim Matrix`) are an **aspirational** restructure
shown in the mockup — the current `Inversion_atlas.html` has the older
numbered scheme (`1 local PCA |z|`, `4 boundaries`, `5 catalogue`, etc.).
**Don't try to align all tab labels to the mockup**. Only add the new
SV Evidence tab; leave the rest of the bar untouched. Spec §10 question
1 covers this.

Three catfish cohorts must NEVER be conflated:

1. F1 hybrid (*C. gariepinus* × *C. macrocephalus*) — genome assembly
   paper only.
2. **226-sample pure *C. gariepinus* hatchery cohort** — current
   inversion work, this is what the SV Evidence page renders.
3. Pure *C. macrocephalus* wild cohort — future paper.

The K=3 PCA clusters in the mockup (H1/H1, H1/H2, H2/H2) reflect
**hatchery broodline structure within pure C. gariepinus**, not species
admixture.

## What's already in the tarball

```
Atlas/
├── Inversion_atlas.html              ← do NOT rewrite; surgical edits only
├── Diversity_atlas.html / Genome_atlas.html / Population_atlas.html
├── atlas_server.py                   ← localhost dev server
├── js/                               ← page-specific JS modules (house pattern)
│   ├── atlas_overview.js   ← good template for atlas_sv_evidence.js
│   ├── atlas_renderers_turn4.js
│   ├── atlas_focal_vs_bg.js
│   └── ... (read these to learn the project's conventions)
├── server_turn1/popstats_server.py   ← FastAPI popgen server
├── _scripts/
│   ├── phase_8_comparative_breakpoint_fragility/   ← my recent work
│   ├── phase_X_ushape_evolution/                   ← my recent work
│   └── ... existing modules ...
├── specs_todo/
│   ├── SPEC_sv_evidence_page.md                    ← READ THIS FIRST
│   ├── _mockups/                                   ← three reference PNGs
│   ├── from_turn129/S6_dosage_heatmap_streaming_viewer.md  ← related
│   ├── from_turn129/S7_karyotype_breakpoint_internal_evidence.md  ← related
│   └── ... other specs ...
├── handoffs_to_implement/
│   ├── HANDOFF_2026-05-04_phase_8_phase_X.md       ← previous deliverable
│   ├── HANDOFF_2026-05-04_next_chat_sv_evidence.md ← THIS FILE
│   └── turn129_bundle/                             ← review docs
└── json/                             ← data layers (heavy TE folders STRIPPED for transit)
```

The four heavy TE/TSD folders (`Cgar_TE_density`, `Cmac_TE_density`,
`Cgar_TSD_density`, `Cmac_TSD_density`) are **deliberately absent** from
this drop — they live on Quentin's LANTA copy. The SV Evidence page does
not need them.

## Recent work to be aware of (so you don't redo it)

- **`phase_8_comparative_breakpoint_fragility/`** — TE/repeat density
  comparative layer for inversion breakpoints. Output schema
  `comparative_breakpoint_fragility_v0.1`. The SV Evidence page does
  **not** consume this; mentioned only so you don't think it's the same
  thing.
- **`phase_X_ushape_evolution/`** — per-candidate evolutionary-shape
  classification (U-shape vs internal peak vs flat-deep). Server endpoint
  `POST /api/ushape/candidate`, atlas JS at
  `_scripts/phase_X_ushape_evolution/atlas_js/ushape_renderer.js`. The
  SV Evidence page does **not** depend on this either; flagged so you
  don't conflate "shape" classification with "pattern label"
  classification — they answer different questions.
- **Matched-background z-scores** in phase_X are **flank-resampled, not
  genomewide**. Same applies if you ever need a similar null on the SV
  page.

## What to build (high-level checklist)

This is a compressed restatement of `SPEC_sv_evidence_page.md §7`. The
spec is authoritative; this is just to scope your first message.

- [ ] step 1 — `js/atlas_sv_evidence.js` skeleton + tab + empty page mount
- [ ] step 2 — `sv_genotype_counts_v1` loader + main SV table (sort,
      filter, pagination, FDR colours, pattern colours, TSV export)
- [ ] step 3 — locus track strip (zones, dosage, PCA, SV icons, genes,
      repeat density) — reuse existing renderers; new code = SV-icons
      track only
- [ ] step 4 — right-rail boundary summary tables + legend
- [ ] step 5 — UpSet panel (right-rail compact + main-area full)
- [ ] step 6 — sample × SV heatmap (right-rail compact + main-area full)
- [ ] step 7 — `_scripts/svgt/STEP_SV_GT_AGG_aggregate_genotype_counts.py`
      stub emitter

After each step, write a tiny test fixture under
`Atlas/tests/sv_evidence/` that the next session can use to regression-
test without real data.

## Open questions to surface to Quentin BEFORE writing code

These are §10 of the spec, restated. Send them as the FIRST message in
the build chat (one short message, three to five bullet points, easy to
respond to on mobile):

1. **Tab numbering** — accept the mockup's `4 boundaries / 5 SV evidence
   / 6 catalogue` renumber, or stay with current numbering and call it
   `5b SV evidence`?
2. **Pattern-label authority** — does the data-side emitter classify
   patterns, or should the atlas compute them from FDR + zone + group
   counts?
3. **`sv_genotype_counts` producer** — does an upstream R/Python script
   already emit per-call × group counts? If yes, in what format? If no,
   build the stub from §7 step 7 and document the schema for the data
   side.
4. **Notes column source** — free text from the producer, or wait for
   S7's cluster evidence layer?
5. **Browser-side recompute fallback** — keep (with the warning
   checkbox) or drop?

Quentin's answers to those five questions unblock the entire build.
Don't start coding until at least #1 (renumber) and #2 (label authority)
are answered — the rest can be deferred.

## House style reminders

- One JS module per page, all `window.AtlasXxx` namespaced. See
  `js/atlas_overview.js` for the template.
- No build step. The atlas is a single HTML loading `<script>` tags
  directly. Modules are vanilla ES with no bundler.
- Reuse `state.*` (the global atlas state object). Read it; only the
  scrubber writes to it.
- Empty / partial states are mandatory — never throw. The user opens the
  atlas before all layers load; render the parts you can.
- Light + dark mode: the atlas auto-themes via CSS variables. Use
  existing variable names (`var(--ink)`, `var(--bg)`, `var(--good)`,
  `var(--accent)`) — don't hardcode colours except where the spec gives
  literal hex values (FDR colour palette, pattern label colours).
- localStorage keys for persistence: `state.svEvidenceFilters` keyed by
  candidate_id. Mirror the page-11 pattern for boundary edits.
- Never invent a sample-id or candidate-id. Read them from `state.*`.

## After you finish

1. Update `Atlas/specs_todo/SPEC_sv_evidence_page.md` `Status:` line to
   `implemented v0.1`. Don't delete the spec — leave it as the design
   document.
2. Drop a sibling `SPEC_sv_evidence_page__implementation_status.md` next
   to it (matching the pattern Quentin already has for
   `SPEC_comparative_te_breakpoint_fragility__implementation_status.md`).
3. Write a HANDOFF doc dated to the day you finish, into
   `handoffs_to_implement/`, summarising what shipped, what's deferred,
   and the next page Quentin should consider building.
4. Bundle a fresh `Atlas_dragdrop.tar.gz` for Quentin to drop into LANTA.

## Last note

Quentin is on a slow mobile connection most of the time. Keep your
turn-by-turn outputs **tight**. Don't dump 800 lines of code in a single
message — incrementally build, ask for sign-off on §1 before §2.
Bundle final artifacts at the end of each major step so he can pull
them down between sessions.

Good luck. — Claude (this chat, 2026-05-04)
