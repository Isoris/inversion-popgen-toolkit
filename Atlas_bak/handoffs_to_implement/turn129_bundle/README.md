# Atlas integration bundle — turn 129

This bundle is a **review package**, not a code drop. Read the
plans first, push back where I got it wrong, then a follow-up turn
lands the actual code against your real `Inversion_atlas.html`.

## Why no code this turn

You said:

> Goal: Cleanly recode the atlas integration instead of patching
> randomly.

I cannot reliably patch your repo from here. My snapshot is from
turn 128 (62,215 lines, 3.0 MB). Your repo at
`/mnt/c/Users/quent/Desktop/inversion-popgen-toolkit/Atlas/` has
clearly diverged — `collectPopstatsTracks`, `popgenRenderers`,
`atlasServer.url`, group dock wiring, server endpoints, and the
`popstats_server_cache` directory all postdate my snapshot. Line
numbers don't match. Function names may not match. Pushing code
into a copy I can't validate is exactly the "patching randomly"
you asked me to stop doing.

So this bundle:

1. Triages your wall of requirements into ordered priorities.
2. Writes each fix as an **anchor-based patch doc** — find by
   string, not by line number.
3. Specs the deferred work (SV tables, indel slope, breakpoint
   scoring, double crossover, interpretation rules).
4. Decision points are flagged for your confirmation before any
   code lands.

## Read order

1. `plans/00_master_plan.md` — full triage and rationale
2. `plans/01_priority_triage.md` — ordered list of 13 patches +
   5 specs
3. **Decision points** at end of `01_priority_triage.md` — confirm
   these before any code lands
4. Browse patches you care about
5. Browse specs for context on deferred work

## Priority order (TL;DR)

| Phase | Patches | What it unblocks |
|---|---|---|
| **P1 (5 fixes)** | const_out, server_url, health_alias, ld_naming, server_startup | Server stability, popgen panel |
| **P2** | group_dock_from_locked_labels | Popgen workflow end-to-end |
| **P3 (3)** | two_row_nav, page_consolidation, js_org | UI clean recode |
| **P4 (2)** | dosage_endpoint, theta_endpoint | Live data without R-side scripts |
| **P5** | boundaries_sv_skeleton | New combined page scaffold |
| **P6** | karyo_tier_catalogue_merge | New combined page |
| **Specs** | S1-S5 | Future SV / indel / breakpoint work |

## Decision points

These need your sign-off before I land code:

1. **Stage naming** — "Evidence" (your message) vs "popgen" (turn
   128). Pick one.
2. **Catalogue location** — merge with karyotype/tier (proposed)
   or keep on Discovery as candidate landing table?
3. **SV merge** — the page named just "Boundaries" (with SV as
   one section) or "Boundaries + SV evidence" as the page name?
4. **Turn 128 pills** — confirm fully removed in P3.1, or kept
   as alternate?

## File listing

```
plans/
  00_master_plan.md            full triage + rationale
  01_priority_triage.md        ordered patch list + decision points

patches/
  P1_1_const_out_fix.md                1-line bug fix
  P1_2_server_url_canonical.md         localStorage key consolidation
  P1_3_health_alias.md                 /health alias regression test
  P1_4_ld_naming.md                    LDSplitReq vs FastLDReq
  P1_5_server_startup.md               POPSTATS_CONFIG race fix
  P2_1_group_dock_from_locked_labels.md  group dock unblocker (~80 LOC)
  P3_1_two_row_nav.md                  two-row nav (replaces turn 128 pills)
  P3_2_page_consolidation.md           SV→Boundaries, Catalogue→Karyo
  P3_3_js_org_audit.md                 production js/ vs tests/
  P4_1_dosage_chunk_endpoint.md        /api/dosage/chunk (~120 LOC server)
  P4_2_theta_overlay_endpoint.md       /api/theta/candidate (STUB — needs find output)
  P5_1_boundaries_sv_skeleton.md       10-section page skeleton
  P6_1_karyo_tier_catalogue_merge.md   3-view internal toggle

specs/
  S1_sv_evidence_tables.md             5-table SV schema
  S2_indel_slope_burden.md             indel layer placement
  S3_breakpoint_scoring_bayes.md       Bayesian beta-binomial framework (future)
  S4_double_crossover.md               recombinant detection (future)
  S5_sv_interpretation_rules.md        discipline doc — INV ≠ raw read, etc.
```

## What to do with this

**Option A** — review the plans, confirm decision points, send
back any pushback. Next turn lands P1.1 + P1.2 + P1.3 + P1.4 +
P1.5 (the five small fixes) against your real file in one pass.

**Option B** — paste the relevant function snippets from your
real file (especially `collectPopstatsTracks` and the group engine
dim collector). I rewrite P2.1 to match your actual code and
ship it ready-to-paste in one block.

**Option C** — skip planning, just take what's already in the
patches, treat as best-effort, apply manually. Each patch doc has
anchor strings so you can grep your way to the right spot.

I'd recommend (A) for the P1 batch and (B) for P2.1 — the P2.1
patch needs to reference your actual function names to be safe.

## What's NOT in this bundle

- Any actual atlas HTML / JS / Python file modifications.
- BUSCO integration (deferred since turn 121).
- The Inversion_atlas.html.bak rotation policy (you have a .bak
  already; that's enough).
- New tests against your real file — only test templates that
  assume the patches landed.

## Notes on style

Each patch doc has the same shape:
- Risk + lines + dependencies + verification at top
- "What's wrong" — the symptom
- "The fix" — the actual change, with anchor strings
- "Verification" — how to confirm it worked
- "Test" — a `tests/test_pX_Y_*.js` template
- "Risk notes" — what could go wrong

The patches are deliberately verbose. Better to over-explain than
under-specify when I can't validate against the live file.
