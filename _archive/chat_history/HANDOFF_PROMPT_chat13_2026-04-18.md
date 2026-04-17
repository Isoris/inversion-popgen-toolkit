# Handoff prompt for chat 13 — registry wiring after chat-12 audit

**Read only these four documents at the start. That's the full
context for chat 13.** Everything else is either code (read on
demand) or historical (ignore unless something here references a
specific past decision you need to reconstruct).

1. This file — `HANDOFF_PROMPT_chat13_2026-04-18.md`.
2. `SESSION_SUMMARY.md` — current-state snapshot.
3. `AUDIT_PHASE4_COHERENCE_2026-04-18.md` — architectural reference.
   Start with §6 "Recommendations for chat 13 wiring"; that is your
   checklist. §1 per-script material is reference-on-demand.
4. `AUDIT_LOG_chat12_2026-04-18.md` — findings register. Table
   format; scan for items assigned to chat 13.

**Do NOT read** unless you specifically need to:

- `_archive/chat_history/` — all pre-chat-12 handoffs and per-chat
  audit logs. They were active during their own chats; none of their
  decisions are still open questions. Earlier fixes already applied
  are summarized in `FIXES_APPLIED.md`.
- `FIXES_APPLIED.md` beyond the chat 12 section — prior chats are
  summarized in one-liners only for cross-referencing.
- Superseded code files in `inversion_modules/*/_archive_superseded/`
  and `inversion_modules/_archive/`.

---

## Posture for chat 13

Chat 12 was architectural — audit before wire. Chat 13 is mechanical
— the audit blessed the wiring paths. Don't re-open design questions
that the audit already settled. If something feels off, check the
audit's finding register first; likely it's already been logged.

All chat-12 code has parse-check passing and 83/83 unit tests passing
(50 DAG derive_R_from_regime + 33 GC detector). The rewrites are on
solid ground; treat them as stable.

---

## Wiring checklist

From audit §6. Immediate, mechanical:

- [ ] `STEP_C01j_regime_compatibility_engine.R` — redirect segment
      output via `write_block("regime_segments", ...)`. Raw
      `regime_memberships` stays as sidecar TSV.gz (too large for
      registry JSON); the block can carry the sidecar path.
- [ ] `STEP_C01l_local_structure_segments.R` — redirect via
      `write_block("local_structure_segments", ...)`.
- [ ] `STEP_C01m_distance_concordance.R` — redirect via
      `write_block("distance_concordance", ...)`.
- [ ] Confirm the `write_block("regime_sample_dag", ...)` call
      already in chat-12's `STEP_C01i_b_multi_recomb.R` round-trips
      through the registry JSON writer. If the `per_sample`
      data.table doesn't serialise cleanly, convert to
      `as.list(as.data.frame(per_sample))` before the write.
- [ ] `STEP_C01i_d_seal.R` reads the new `recomb_subgroup` field
      from multi_recomb's output; registers
      `inv_<cid>_RECOMBINANT_GC` and `inv_<cid>_RECOMBINANT_DCO`
      groups alongside the existing `inv_<cid>_RECOMBINANT`.
- [ ] `characterize_candidate.R` Q2 pulls the new keys:
      regime_segments (`q2_regime_*`), local_structure_segments
      (`q2_boundary_*`), distance_concordance, and regime_sample_dag
      (`q2_dag_*`). Audit §1.19 has the full list.

Threshold fixes (low risk, small diffs):

- [ ] **AT**: C01l `flank_bp` scales with span:
      `flank_bp = clamp(span_bp, 200000L, 500000L)`.
- [ ] **AV**: decompose `decomp_quality = silhouette >= 0.40 ?
      "clean" : "noisy"`. Non-gating annotation only.
- [ ] **BB**: `lib_ghsl_confirmation.R::resolve_karyo_bp` — add
      `stopifnot(uniqueN(annot_dt$global_window_id) == nrow(annot_dt))`
      at top.

Path verification:

- [ ] **AQ**: grep C01e for hardcoded C01j paths; update if they
      still reference pre-chat-11.5 orphan locations. Panel H is
      the specific panel that reads C01j output.
- [ ] **BD**: C01k (`STEP_C01k_annotated_simmat.R`) — confirm reads
      from the registry block locations produced by this chat's
      wiring.

Metrics:

- [ ] **BC**: run `compute_candidate_status.R::build_key_spec()`
      after wiring and compute wired-fraction. Target 85%+.

Doc:

- [ ] **AO**: C01d header line — "10 dimensions" → "12 dimensions"
      correction.
- [ ] **AY**: `lib_step03_seed_loader.R` header — note that drop-
      conflict is stricter than the chat-9 doc's "priority
      flashlight" language.

---

## Not in scope for chat 13

- D13–D15 gate-promotion in C01d. Wire the new blocks as annotation-
  only first; promote to gates after HPC calibration in chat 14.
- `boundary_v2` schema (Finding AR) — chat 14 after HPC run gives
  real boundary examples to calibrate against.
- Pattern enum extensions (AP, BE) — chat 14.
- C01j `structure_score` cutoff recalibration (AS) — chat 14.
- `min_dco_bp` recalibration (AW) — chat 14.
- Manuscript terminology (BF) — for the writeup, not the pipeline.

## Do not touch

- `registries/api/*/registry_loader.*` — mature from chat 11.
- Existing frozen schemas (regime_segments, local_structure_segments,
  distance_concordance, boundary_scan, regime_sample_dag,
  gene_conversion_tracts). ONLY ADD new schemas this chat if needed.
- `lib_recomb_combination.R`, `gene_conversion_detector.R`,
  `STEP_C01i_b_multi_recomb.R`, `STEP_C01i_decompose.R` compute
  logic. Chat 12 rewrote these; 83 tests lock behaviour. Touch only
  to fix clear wiring bugs surfaced by chat 13's integration run.

## Definition of done

1. All `write_block` calls in place for regime_segments /
   local_structure_segments / distance_concordance / regime_sample_dag.
2. Seal registers the two new recombinant subgroup names.
3. Characterize Q2 reads the new keys.
4. Threshold fixes AT / AV / BB applied.
5. Path verifications AQ / BD resolved.
6. Wired-fraction BC computed and recorded.
7. Doc fixes AO / AY applied.
8. Existing 83/83 tests still pass (no regressions).
9. `AUDIT_LOG_chat13_<date>.md` created with findings BG+ (if any).
10. `SESSION_SUMMARY.md` and `FIXES_APPLIED.md` updated in place
    (rolling docs — add chat 13 section at top, don't create dated
    copies).
11. Previous handoff + audit log archived:
    `mv HANDOFF_PROMPT_chat13_2026-04-18.md AUDIT_LOG_chat12_2026-04-18.md
    _archive/chat_history/`.
12. Create new `HANDOFF_PROMPT_chat14_<date>.md` at root for the
    HPC run.
13. Tarball as
    `inversion-popgen-toolkit_chat13_wiring_<date>.tar`.

Chat 14 then does the first HPC run on LANTA — LG12 first (strong
inversion), LG25 second (paralog-heavy control), then sweep.

---

## Housekeeping pattern going forward

Chat N:
- Reads: root-level `HANDOFF_PROMPT_chatN_*.md`, `SESSION_SUMMARY.md`,
  `AUDIT_PHASE4_COHERENCE_*.md` (or newer successor), and the
  previous chat's audit log if still at root.
- Writes: new `AUDIT_LOG_chatN_*.md` at root. Updates rolling
  `SESSION_SUMMARY.md` and `FIXES_APPLIED.md` in place. Creates
  `HANDOFF_PROMPT_chat(N+1)_*.md` at root.
- At end of chat: archives the previous chat's handoff + audit log
  to `_archive/chat_history/`.

This keeps the root at ~6 docs (handoff, session summary, fixes,
coherence audit, current chat audit log, README/LICENSE). Each new
chat only reads 3–4 documents to bootstrap. Archive is append-only.
