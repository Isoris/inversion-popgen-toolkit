# ADR-14: Four-file atlas architecture (eventual)

**Status:** Locked, deferred implementation
**Date:** April 30, 2026
**Supersedes:** SESSION_AUDIT_apr30 §8 open question on diversity-atlas
file-vs-mode-switch
**Related:** ADR-13 (Q + Ancestry pages), SCHEMA §17, §23, §28

---

## Decision

The atlas suite will eventually consist of **four separate HTML files**,
each its own self-contained product:

| File | Scope | Audience |
|---|---|---|
| `Inversion_atlas.html` | Phase-2 inversion discovery and per-candidate refinement (current Inversion Atlas pages 1, 12, 2, 11, 3, 4, 6, 7, 8, 9, 10, 5; plus future Q page, Ancestry page, phase-4 CUSUM panels) | Manuscript reviewers; the PI running candidate review |
| `Population_atlas.html` | Samples / QC / families / heterozygosity / diversity / inbreeding (already exists as `population_atlas_v1.html`) | Breeding-program collaborators; cohort QC |
| `Diversity_atlas.html` | θπ, F_ROH, ROH atlas, ancestry-Q stratification, kinship/relatedness, family-LD — diversity-perspective lens on the same per-window data | Population genetics reviewers; comparative-diversity questions |
| `Genome_atlas.html` | Chromosome-scale assembly stats, synteny / ancestral karyotype, gene tracks, repeats / TE landscape, conserved elements, variant annotations | Genome-assembly reviewers; the F1-hybrid assembly manuscript's Phase-2 figures |

The four files share no DOM, no JS state, no localStorage namespace —
each loads its own JSONs and renders independently. They cross-link via
right-aligned tab-bar buttons (the pattern already established by the
Population Atlas link in the current Inversion Atlas).

## Promotion timeline

Given the rationale is Claude-session ergonomics, the right time to
promote is **as soon as a session is going to do non-trivial work
that's scoped to one sub-product**. Not when renderers exist.

Specifically:

**Diversity_atlas.html** → split out **before** the diversity-page
renderers get built, not after. Splitting first means every Diversity-
focused session loads ~300 KB instead of 1.9 MB. Probably the first
file to promote.

**Genome_atlas.html** → same logic. Split out at the start of any
genome-annotation work, not at the end.

**Inversion_atlas.html** → keep `atlas.html` as the working alias
through LG28 dry-run + page-3/12 stabilization. This is the file
under most active edit right now, so renaming it mid-debug invites
breakage. Rename to `Inversion_atlas.html` once page 3 and 12 are
working — that's the natural quiet point.

**Population_atlas.html** → already exists as `population_atlas_v1.html`.
Rename to `Population_atlas.html` for consistency at the same moment
`atlas.html` becomes `Inversion_atlas.html`.

The promotion step itself is mechanical: copy the relevant
`<button data-page="pageNN">` and the matching `<section id="pageNN">`
out into a new file, plus the shared CSS / scrubber / JSON-loader
infrastructure. The shared infrastructure has to be duplicated — that's
the cost — but each file then becomes self-contained and small.

## Why eventually four files (not one big atlas)

The primary reason is **Claude-session ergonomics**, not end-user
ergonomics.

The atlas HTML is currently 1.9 MB / ~39,000 lines. Every Claude
session that touches it has to load the entire file into context,
even when the work is scoped to a single sub-product (e.g. a
diversity-page tweak doesn't need the candidate-management code).
This burns tokens and increases the chance Claude conflates
unrelated sections during edits.

Splitting along natural product boundaries means future sessions
can load **just** the file relevant to the work — probably 200–400 KB
per file once split, vs 1.9 MB for the whole monolith. Inversion
Atlas changes don't have to fight Genome Atlas changes for the same
file slot in context.

Secondary benefits:

The four atlases have genuinely different scopes and audiences. A
manuscript reviewer reading the diversity argument doesn't need
inversion candidate management; a breeder reading population QC
doesn't need genome assembly stats. Per-file loading is faster on
the user side too.

The Population Atlas precedent is already established — splitting
was the right call once.

When one atlas needs a major refactor (e.g. phase-4 CUSUM adding
per-candidate sub-systems and per-carrier boundary distributions to
the Inversion Atlas), the work is contained to one file. Cross-
contamination risk between unrelated atlases drops to zero.

## Why not promote ALL of them right now

Even though the rationale favors splitting early, three reasons to
**stage** the promotion rather than do all four at once:

**LG28 dry-run is the critical path.** The Inversion Atlas is mid-debug
on the page-3/12 renderers. Renaming files in the middle of a debug
loop creates path-mismatch bugs and makes the next session waste turns
chasing "why does the runbook say `atlas.html` but the file is
`Inversion_atlas.html`." Wait for the natural quiet point.

**Diversity and Genome scaffolds are tiny.** Pages 13 and 14 today are
~50 lines of HTML each — explanatory `<div>` blocks. Splitting them out
before any real work touches them is symmetric churn for no token saving
yet. Split when the **next** session is going to add real content to
that file, not now as a preemptive reorg.

**One promotion at a time, not a batch.** Each promotion has the
mechanical risk of duplicating shared infrastructure incorrectly. Doing
four at once means four chances to ship a broken JSON-loader copy.
Sequencing them lets each go through its own smoke test.

## What changes today

- This ADR is added.
- Tab tooltips for pages 13 and 14 already say "preview / scaffold" —
  no change needed.
- Future Claude / future Quentin: do not re-litigate "should diversity
  be a tab or a separate file." It will be a separate file. Just not
  yet.

## What does NOT change today

- `atlas.html` filename. Stays as the working alias.
- Pages 13 and 14 stay as scaffold pages inside the Inversion Atlas.
- `population_atlas_v1.html` stays where it is (separate file already).
- Cross-link button at the tab bar (`→ Population Atlas`) stays.

## Promotion checklist (when the time comes)

For each atlas file getting promoted out:

1. Confirm the data-layer JSONs are produced by the cluster pipeline
   and have been smoke-tested through the scaffold page.
2. Copy the relevant section(s) out of `atlas.html` into a new file.
3. Strip out the multi-mode tab bar — each file is single-mode.
4. Rename localStorage keys to namespace by atlas (e.g.
   `pca_scrubber_v3.atlasMode` → `inversion_atlas.uiState`).
5. Add right-aligned cross-link buttons for the other three atlases
   in the tab bar of each file.
6. Update `Atlas/README.md` to reflect the four-file layout.
7. Update `JSON_CONTRACT.md` if the schema split breaks at the file
   boundary.

Do this **once per atlas**, only when that atlas's renderers are real.
Don't batch the promotions.

## Walk-back triggers

The decision should reverse only if shared infrastructure becomes a
maintenance burden — i.e. if scrubber, JSON loader, theme toggle, etc.
need synchronized changes across all four files often enough that
keeping them in sync costs more than the per-session token savings.

If that emerges, the right answer is to extract shared infrastructure
into a separate `atlas_common.js` and `atlas_common.css` loaded by all
four HTML files, **not** to re-merge into one HTML.

If a future session proposes "let's just keep everything in one big
atlas, the cross-mode toggle works fine" — point them at the token
cost. A 1.9 MB file is fine for the user but expensive for every
Claude session that has to load it. The split is for Claude's
workflow, not the user's.
