# Handoff prompt for chat 9

Continuing from chat 8. Tarball is the chat-8 output. Phase 4b
understanding is complete and three fixes landed. Four findings
flagged (S, T, U, V + chat-7's 8/9 still open).

**Chat-8 fixes (all need LANTA Rscript parse-check):**

- FIX 36 (cosmetic) — 4b/STEP_C01i_b_multi_recomb.R L316: decision-rule
  comment rewritten to match the code (and the README).
- FIX 37 (crash, dormant) — 4b/STEP_C01i_b_multi_recomb.R L233: undefined
  `het_carriers` in Signal 3 hemizygous detector replaced with proper
  per-row list-column mask.
- FIX 38 v2 (silent, impactful) — 4b/lib_decompose_helpers.R L88–96:
  `%||%` precedence bug was inflating every window interval by 100 kb,
  which inflated every recombinant mosaic_length_bp by 100 kb. Against
  the real cheat24 thresholds (50 kb / 200 kb), most true gene-conversion
  events were being lost to "ambiguous". v2 stripped the dead fallback
  entirely because C01a always writes real start_bp/end_bp columns.

## Posture for chat 9

**Understand before fixing.** Chat 8 had a near-miss on this: FIX 38
v1 was committed with an incorrect impact story because I hadn't
traced what the code was actually doing. Quentin's pushback on the
arbitrary 50000 magic number surfaced that the fallback was dead
code, which changed the correct fix. Lesson internalized: when
something looks like a bug, understand its role in the real data
flow first. Trace magic numbers to their upstream source. Trace field
names to their writer. Don't assume a pattern-match on "this looks
wrong" is enough.

This means reading goes before fixing, but reading-with-fixing is
still the goal — not audit-only. When a fix IS obvious AND surgical
AND understood, apply it inline. When the understanding is partial
or the call is non-obvious, flag it for Quentin instead of guessing.

**Three-fix ceiling stays in effect** (per chat-7's handoff standard).
If more than three real, understood fixes emerge, stop at three and
note the rest for chat 10.

## Scope options for chat 9

Pick one; don't mix:

**Option A: phase 4c (C01f hypothesis tests) — RECOMMENDED.**
`STEP_C01f_hypothesis_tests.R` at 2512 lines is the biggest R script
in the tree and the central consumer of everything phase 4b produces.
Three chat-7 patches (01/02/03 in `patches/`) have already surgically
modified it: `comp_from_registry` (chat 4 FIX 20), promotion_cap
read (chat 7 FIX 31), jackknife semantics. A whole-chat audit makes
sense here because the consumer logic is where correctness of the
whole phase-4 pipeline converges.

**Option B: phase 4d + 4e combined.** Fewer lines but more
cross-script contracts (4e's `compute_candidate_status.R` reads
results from 4c AND Layer D from phase 3 AND composite flag from 4b).
Chat 7's FIX 32 modified 4e/compute_candidate_status.R and
4c/group_validation_gate.R together; chat 9 would verify the
end-to-end flow those touch.

**Option C: a design session for the deferred items.** Finding S
(STEP03 seed wiring) and Finding V (cheat24 threshold divergence)
both need Quentin's input before they're actionable. A short
decision-focused session could unblock both, then chats 10+ apply
the decisions as code changes.

My recommendation: **Option A** if you have budget for a long chat,
**Option C** if you want to close open questions quickly.

## If doing Option A (phase 4c audit)

**Reading order (understand first, then look for bugs):**

1. `AUDIT_LOG_chat8_2026-04-17.md` — state of 4b you're inheriting,
   including the four deferred findings.
2. `4c_group_validation/README.md` (if it exists) — the job 4c is
   supposed to do. If no README, look at the section headers of
   `STEP_C01f_hypothesis_tests.R` as a rough TOC.
3. `ls 4c_group_validation/` — full inventory. Note which scripts
   are main vs helpers.
4. The three patches in `4_postprocessing/patches/`:
   - `01_C01f_comp_from_registry.R` — establishes the band1/band2/band3
     vocabulary and the registry-reader contract.
   - `02_C01f_promotion_cap.R` — chat 7 FIX 31's wiring.
   - `03_C01f_jackknife_semantics.R` — the jackknife → family_linkage
     + polymorphism_class contract.
   These tell you what has already been surgically edited into C01f
   and what contracts those edits establish. **Don't re-audit the
   patches themselves** (chat 5/7 did that); audit the behaviour
   they produce when applied.
5. `STEP_C01f_hypothesis_tests.R` — navigate by section, not linearly.
   Use grep to find specific patterns. The file is ~2500 lines; a
   sequential `view` through it will exhaust tool budget.

**What to actually understand before calling anything a bug:**

- What is `comp`? Where does it come from? (4b registry, on-the-fly
  k-means, or some fallback?) What columns does each path produce?
- What are the hypothesis tests T1/T2/T7/T8/T10 testing? (Not the
  code — the biology. A one-sentence gloss per test before reading
  any of their implementations.)
- What does `compute_group_validation` actually validate? What are
  Layers A/B/C/D? What does "VALIDATED" vs "UNCERTAIN" mean in terms
  of downstream gates?
- What does the jackknife compute, and why is the output a
  classification instead of a test statistic?
- Which flat keys does C01f write vs read? Cross-reference with the
  seal writes (chat 8 Finding U identified forward-declared keys) to
  spot consumer gaps.

**Rough audit checklist once understanding is in place:**

A. Comp-from-registry + fallback chain works end-to-end for both
   registry-present and registry-empty cases. Chat 4 FIX 20 closed
   one failure mode; verify no others.

B. `compute_group_validation`'s VALIDATED promotion gate reaches its
   Layer D Fisher reads correctly (post chat-7 FIX 31). The
   `fisher_or > 5` clause (chat-7 FIX 32) agrees between 4c and 4e.

C. Jackknife overwrites `q6_family_linkage` and `q6_polymorphism_class`
   AFTER seal initializes them. Timing matters — check the execution
   order in the SLURM DAG vs the code.

D. End-to-end: does 4c's output feed 4d/4e cleanly? Column contract
   check on whatever file C01f writes (typically `hypothesis_results.tsv.gz`).

E. Forward-declared keys pattern (chat 7 Finding 8, chat 8 Finding U)
   — is C01f writing any keys that have no live reader? Flag them
   alongside the existing forward-declared set.

## Budget planning

Lessons from chats 6/7/8:

- 30% reading code + reading/verifying upstream context (C01f
  especially — don't read it linearly, navigate by section)
- 10% tracing numbers / names upstream to verify understanding
  before committing to fixes (this is the part chat 8 skimped on)
- 40% fixes + audit writeup
- 20% repack

If a fix starts to feel un-surgical — multiple touch points, unclear
downstream consumer, magic numbers you haven't traced — STOP, flag
it, move on. The three-fix ceiling is a ceiling, not a target.

## Open items carried forward

From chat 5 / chat 8 Finding S:
  STEP03 seed wiring into C01i_decompose. Design call needed.

From chat 6:
  The lib design work. Blocks the registry-block path end-to-end.

From chat 7 Findings 8 / 9:
  C01g silent skip on missing helper; 4a registry block writers.
  Both lib-pending.

From chat 8:
  - Finding T: README structure_type drift (4 listed, 7 produced)
    + decompose header comment claiming "per-window dosage track"
    when it's actually a class track. Batch for doc pass.
  - Finding U: forward-declared flat keys `q1_composite_flag`,
    `q2_decomp_quality_flags`. Consumers emerge in 4d/4e or lib chat.
  - Finding V: cheat24 threshold divergence. Real cheat24 uses
    50kb/200kb; inline fallback uses 100kb/500kb. Same candidate
    gets different event_class depending on flashlight availability.
    Needs Quentin to pick authoritative pair.

## Parse-check backlog (cumulative)

Nine R scripts + one Python script need LANTA
`Rscript -e "parse(file=...)"` + smoke tests before production:

Chat 7 (six files):
- phase_3_refine/03_statistical_tests_and_seeds.py — Python AST-parse OK, smoke-tested
- 4c/STEP_C01f_hypothesis_tests.R
- 4c/group_validation_gate.R
- 4e/compute_candidate_status.R
- 4a/STEP_C01d_candidate_scoring_wired_25_v934_registry.R
- 4a/STEP_C01g_boundary_catalog_wired_4_8_10_11_17_21_v934_registry.R
- patches/01_C01f_comp_from_registry.R (template update)

Chat 8 (two files):
- 4b/STEP_C01i_b_multi_recomb.R
- 4b/lib_decompose_helpers.R

Container has no R. Parse-check needs LANTA.

## End of chat

AUDIT_LOG_chat9_2026-04-17.md, append to SESSION_SUMMARY_2026-04-17.md
and FIXES_APPLIED_2026-04-17.md. Repack as
`inversion-popgen-toolkit_chat9_fixes_plus_phase4c_audit_2026-04-17.tar`
(or `phase4de` if you went Option B, or `design_decisions` if C).

If budget runs low on the writeup stage: abbreviated audit log is
better than no audit log. Session summary and FIXES_APPLIED appends
can be short bullet-points if needed, as long as the fix descriptions
carry enough detail for chat 10 to re-verify them.
