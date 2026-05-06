# Atlas integration plan — turn 129 (post-handoff scope expansion)

You sent me a wall of requirements covering UI, server, data bridges,
SV interpretation rules, indel slope, future breakpoint scoring, and
double-crossover labels. Final instruction: **"For the other stuff
write to specs."**

This is the triage. Read this first. Nothing should land in your repo
based on this turn until you confirm the plan.

---

## State-of-the-world disclaimer

My copy of the atlas is from **turn 128** (62,215 lines, 3.0 MB). Your
repo at `/mnt/c/Users/quent/Desktop/inversion-popgen-toolkit/Atlas/`
has clearly diverged — you mention `collectPopstatsTracks`,
`popgenRenderers`, `atlasServer.url`, group dock wiring, server
endpoints (`/api/health`, `/api/popstats/groupwise`, `/api/dosage/chunk`),
and a `popstats_server_cache` directory. None of that is in my
snapshot.

**I cannot reliably patch your repo from here.** Line numbers don't
match, function names may not match, and silent diffs will compound
into the "patching randomly" problem you specifically said you want
to avoid.

What I CAN do, and am doing, is:
1. Write specs and patch templates with **anchor strings** (function
   names, distinctive comments) instead of line numbers.
2. Make the patches small, atomic, and testable so you (or another
   session with the live file mounted) can apply them cleanly.
3. Be honest about what's a 1-line fix vs a multi-day re-architect.

---

## Priority triage

### Priority 1 — small, atomic fixes (do these first, ~1 hour)

1. **`const out =` reassignment in `collectPopstatsTracks`**
   ([P1.1] in patches/) — 1 line edit, real bug.
2. **Server URL canonicalization** ([P1.2]) — collapse the four
   localStorage keys into one canonical with documented fallback
   order.
3. **Server `/health` alias preservation** ([P1.3]) — ensure the
   patched alias survives re-edit.
4. **`LDSplitReq` vs `FastLDReq` naming** ([P1.4]) — pick one and
   document; remove silent compat alias.
5. **Server config startup race** ([P1.5]) — ensure
   `POPSTATS_CONFIG` env var is honored when uvicorn re-imports.

These are five small fixes. Together they unblock the popgen page
and group dock from completing wiring.

### Priority 2 — group dock from `locked_labels` (do this second)

6. **Build candidate dimensions from `locked_labels`** when
   `fish_calls` is undefined ([P2.1]).
   - Convert: 0 → H1/H1, 1 → H1/H2, 2 → H2/H2 (with `g0/g1/g2`
     fallback for K≠3).
   - Group dock should see: `per candidate · <id>`, `diploid_class@<id>`
     chips with sample counts.
   - Compute button calls `/api/popstats/groupwise` with explicit
     sample-ID lists.

This is the unblocker for the entire popgen workflow. ~50 lines of
JS in the group engine.

### Priority 3 — UI nav restructure (third)

7. **Two-row navigation** ([P3.1]) — Replace the in-row pill
   expansion (turn 128) with the cleaner two-row pattern you
   described:
   - Row 1: main section buttons only (`#atlas_main_nav`)
   - Row 2: subpages of active main (`#atlas_sub_nav`)
   - State: `activeMain`, `subButtonsByMain`
   - No appending children to the right of clicked main button

8. **Page consolidation** ([P3.2]):
   - SV evidence into Boundaries → page becomes
     "Boundaries + SV evidence"
   - Catalogue into Karyotype/Tier → page becomes
     "Karyotype / Tier / Catalogue"

   New main+sub structure I'm proposing (you confirm):

   ```
   Overview     (no subs — landing page)
   Discovery
     ├─ local PCA |z|         (page1)
     ├─ local PCA θπ          (page12)
     ├─ local PCA GHSL        (page15)
     ├─ candidate focus       (page2)
     └─ negative regions      (page19)
   Refine
     ├─ Boundaries + SV       (was page11; absorbs SV evidence)
     ├─ Karyotype / Tier / Catalogue (was page4; absorbs page3)
     └─ Windows               (page8)
   Evidence
     ├─ Popstats              (page6)
     └─ Ancestry              (page7)
   Compare
     ├─ Cross-species         (page16)
     └─ Multi-species         (page16b)
   Output
     ├─ Confirmed             (page9)
     ├─ Annotation            (page21)
     ├─ Markers               (page10)
     ├─ Stats profile         (page17)
     ├─ Marker panel          (page18)
     └─ Overview spreadsheet  (page_overview)
   Help                       (no subs)
   ```

   That's 6 main sections + Overview/Help (2 subless). Five sections
   have multiple subs — exactly the right granularity for two rows.

9. **JS file org** ([P3.3]) — production modules in `Atlas/js/`,
   tests in `Atlas/tests/`, fix script tags. Verify dependency
   order. Audit complete in `plans/03_js_org_audit.md`.

### Priority 4 — data bridges (fourth)

10. **Dosage chunk bridge** ([P4.1]) — server endpoint
    `GET /api/dosage/chunk?chrom&start&end&cap` reading
    `04_dosage_by_chr/<chrom>.{sites,dosage}.tsv.gz`. Returns the
    JSON shape the existing `renderDosageHeatmap` expects.
    Spec includes the schema reverse-engineered from the renderer's
    `state.data.dosage_chunks` reads.

11. **Theta candidate overlay** ([P4.2]) — wire
    `het_bridge_*.sample_tP_windows.tsv.gz` and
    `population_mean_tP_windows.*.tsv.gz` through a server
    endpoint. Distinct from the full θπ scrubber (which still
    needs `theta_pi_per_window`/`theta_pi_local_pca`/
    `theta_pi_envelopes`).

### Priority 5 — Boundaries + SV evidence skeleton (fifth)

12. **Page skeleton** ([P5.1]) — empty-state tabs for the 10
    sections you defined. No data wiring yet. Just the structural
    scaffold so you can iterate section by section.

13. **SV interpretation discipline baked into the renderer**
    ([P5.2]) — empty-state copy and tooltips encode the rules:
    INV ≠ BND, BND not raw read, evidence tiers, body vs boundary
    zones.

### Priority 6 — Karyotype / Tier / Catalogue combined page (sixth)

14. **Combine page4 (karyotype/tier) + page3 (catalogue)** ([P6.1])
    — single page, three views (karyotype / tier / catalogue) as
    sub-tabs or side-by-side. The catalogue table becomes the final
    candidate table with: candidate ID, karyotype structure, tier,
    confidence, boundary status, notes.

### Specs only — DO NOT IMPLEMENT THIS TURN

These are spec docs in `specs/`. No code.

- **SV evidence tables** ([S1]) — the 5-table schema you defined
  (`sv_variant_catalog.tsv`, `sv_sample_genotypes.tsv`,
  `candidate_group_membership.tsv`, `sv_group_genotype_counts.tsv`,
  `sv_group_enrichment_results.tsv`). Includes the pattern-label
  vocabulary.
- **Indel slope / indel burden layer** ([S2]) — placement on
  Boundaries page as "indel-burden support" not breakpoint proof.
  Schema for the layer JSON.
- **Karyotype-stratified breakpoint evidence scoring** ([S3]) —
  Bayesian beta-binomial framework, MAPQ0 use, evidence
  clustering. Future module spec.
- **Double-crossover / recombinant extension** ([S4]) — paired
  internal transition logic. Future module spec.
- **SV interpretation rules** ([S5]) — the discipline doc that
  governs what we can and cannot claim.

---

## What this turn ships

I'm shipping a **review bundle** — you read it, push back where I
got it wrong, then a follow-up turn lands code with anchor-based
patches against your actual file.

What's in the bundle:

```
atlas_handoff_v2_turn129.tar.gz
├── plans/
│   ├── 00_master_plan.md                       (this file)
│   ├── 01_priority_triage.md                   (P1-P6 ordered list)
│   ├── 02_nav_restructure_plan.md              (the new two-row nav)
│   ├── 03_js_org_audit.md                      (script tag audit)
│   └── 04_server_consolidation.md              (server cleanup plan)
├── patches/
│   ├── P1_1_const_out_fix.md                   (1 line)
│   ├── P1_2_server_url_canonical.md
│   ├── P1_3_health_alias.md
│   ├── P1_4_ld_naming.md
│   ├── P1_5_server_startup.md
│   ├── P2_1_group_dock_from_locked_labels.md   (~50 lines JS)
│   ├── P3_1_two_row_nav.md
│   ├── P3_2_page_consolidation.md
│   ├── P4_1_dosage_chunk_endpoint.md
│   ├── P4_2_theta_overlay_endpoint.md
│   ├── P5_1_boundaries_sv_skeleton.md
│   └── P6_1_karyo_tier_catalogue_merge.md
├── specs/
│   ├── S1_sv_evidence_tables.md
│   ├── S2_indel_slope_burden.md
│   ├── S3_breakpoint_scoring_bayes.md
│   ├── S4_double_crossover.md
│   └── S5_sv_interpretation_rules.md
└── README.md
```

Each `patches/Pn_*.md` file contains:
- Anchor strings to find the right place in your file
- The exact diff (before / after blocks)
- A test or verification step
- Line-count estimate
- Risk notes

You can pick which patches to apply. I'm not committing to landing
code without you confirming the plan first — that was the explicit
instruction.

---

## What I am NOT going to do this turn

- Write large new modules from scratch into a 70k-line file I can't
  see.
- Guess at the names of functions in your live atlas
  (`collectPopstatsTracks`, `popgenRenderers`,
  `_atlasGroupEngine_addCandidateDims`, etc.).
- Land the BUSCO integration (still deferred).
- Touch any of the four atlas HTML files (Inversion / Diversity /
  Genome / Population) without seeing them first.
- Write code for "future read-evidence algorithm" — that's spec-only
  per your instruction.

---

## Recommended next action

1. Read `plans/00_master_plan.md` (this file) and
   `plans/01_priority_triage.md`.
2. Push back on anything in the proposed page/nav layout — especially
   whether `Overview` page should exist as a landing, and whether
   `Evidence` is the right name for the popgen+ancestry stage (I
   considered `popgen` per turn 128 but you said `Evidence` in your
   message).
3. Either:
   - **(a)** confirm the plan and I prepare turn 130 to land
     P1+P2 patches against your real file (you paste the affected
     functions or run me with the file accessible), OR
   - **(b)** tell me which patches are wrong and we re-spec.

This is a deliberate planning turn. The wall of requirements deserves
a thought-through structure rather than another big drop of code.
