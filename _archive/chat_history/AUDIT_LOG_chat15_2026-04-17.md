# AUDIT_LOG_chat15_2026-04-17.md

Short audit log for chat 15. The work was narrow in scope (BK
schema canonicalisation + rename + per-key docs); full detail is in
the chat-15 section of `FIXES_APPLIED.md`.

## Session framing

Chat 15 was asked (by the chat-14 handoff) to be "the first HPC run
on LANTA" with six priority items spanning live validation,
threshold calibration, and BK schema-extraction. In practice the
chat-15 environment had no R interpreter, no LANTA access, and no
way to run any of the first five items. Chat 15 did item 6 (BK
schema-extraction) and the user added scope on scientific renaming
and per-key documentation. HPC work is now all queued for chat 16.

## Findings raised in chat 15

- **BK** (closed in-session) — three canonical phase-4b schemas
  written/rewritten; 41 `keys_extracted` directives wired;
  `compute_candidate_status.R` q2/q6 vectors updated. Validated by
  three independent checks: JSON parse + draft-07 meta-schema,
  internal `from`-path resolution, source-code cross-check. See
  `FIXES_APPLIED.md`.

- **BK-rename** (closed in-session) — 31 of 41 keys renamed from
  procedural / mechanical labels to scientific names (what the key
  measures about the genome, not what software step emitted it).
  Applied systematically via `_bk_rename.py` across the three
  schemas + `compute_candidate_status.R`. Preserved: keys whose
  existing names were already scientific (`q2_bic_gap_k3_vs_k2`,
  per-class seed counts, per-class prelim counts, `q1_composite_flag`
  which is the v10.1 canonical).

- **BK-docs** (closed in-session) — new file
  `registries/schemas/structured_block_schemas/BK_KEYS_EXPLAINED.md`
  (~500 lines). Per-key provenance / biological meaning / downstream
  consumer for every one of the 41 keys, plus intro sections on
  what Q1/Q2/Q6 mean, what a "candidate" is, what the three blocks
  characterise about an interval, and a trailing rename map with
  reasoning for each rename.

## Mistakes made and corrected

- **First-pass BK FIXES_APPLIED.md stated `q2_decomp_quality_flags`
  had been renamed to `q2_decomp_quality`.** Incorrect — they are two
  distinct keys (the first is seal-written, a comma-joined aggregate
  of validation quality_flags; the second is the silhouette-derived
  clean/noisy flag from the internal_dynamics schema). Caught while
  searching for external consumers of renamed keys. Fix: re-added
  `q2_decomp_quality_flags` to the q2 vector with a clarifying
  comment; renamed `q2_decomp_quality` to
  `q2_pca_cluster_separation_flag` in the rename pass to make the
  semantic distinction more visible at key-name level.

- **First-pass BK keys were procedurally named.** User feedback
  "the names of the keys for some irs not scientific" triggered the
  rename pass. Resolution was a systematic review of every key
  with a 3-column table (current / proposed / rationale), then
  applied via `_bk_rename.py`.

## Validation surface

Three tools were added at work root:
- `_schema_check.py` — JSON + draft-07 + `from`-path resolver.
- `_code_field_check.py` — schema-vs-source cross-check. Known
  limitation: R list() parser under-reports multi-line if/else
  values in strict mode; permissive fallback catches them. Not a
  blocker for chat-15 schemas; verify by direct grep if in doubt.
- `_bk_rename.py` — rename applier. Kept in-tree as audit trail of
  what was renamed.

All three were run to completion before tar'ing up. Final result:
41 / 41 keys structurally valid, 41 / 41 `from` paths resolve to a
source-code field.

## Open questions for chat 16

None blocking. The BK schemas have not executed against real data —
that's the main chat-16 validation. If any key is missing from
`keys.tsv` after a live run, the `BK_KEYS_EXPLAINED.md` entry for
that key + the block's writer source is where to start.

## Coverage delta

| Metric | Chat 14 end | Chat 15 end |
|---|---|---|
| Q2 | 22 / 73 = 30.1% | 59 / 91 = 64.8% |
| Registry-wide | 73 / 319 = 22.9% | ~91 / 319 = ~28.5% |

Denominator grew because the three canonical schemas exposed keys
that were being written inside JSON blocks but not tracked in
`build_key_spec()`. Meets the handoff's 50-60% BK target with a
margin.
