# `_archive/chat_history/`

Historical per-chat handoffs, audit logs, and dated rolling-doc
snapshots from chats 3 through 11.5 + chat 12 handoff-in.

**Normal working pattern: don't read these.** They describe the
state of the project during their own chats. The current state lives
at the repository root (`SESSION_SUMMARY.md`, `FIXES_APPLIED.md`,
`AUDIT_PHASE4_COHERENCE_*.md`, the current chat's audit log, and the
handoff for the next chat).

**Open these only when:**

- A root-level document references a finding or decision by its
  chat number and you need the original context (e.g. "Finding V,
  open since chat 8" → open `AUDIT_LOG_chat8_2026-04-17.md`).
- Debugging a behaviour that pre-dates the chat-12 audit and you
  suspect an earlier design choice is the cause.
- Building a manuscript methods section that needs to trace the
  pipeline's development history.

## Contents

**Per-chat audit logs (chronological):**
- `AUDIT_LOG_chat3_2026-04-17.md` — initial Phase-4 framework.
- `AUDIT_LOG_chat4_2026-04-17.md` — scoring dimensions.
- `AUDIT_LOG_chat5_2026-04-17.md` — boundary catalog.
- `AUDIT_LOG_chat6_2026-04-17.md` — cheat framework.
- `AUDIT_LOG_chat7_2026-04-17.md` — flashlight integration.
- `AUDIT_LOG_chat8_2026-04-17.md` — cheat24 drift surfaced here.
- `AUDIT_LOG_chat9_2026-04-17.md` — 4b architectural redesign.
- `AUDIT_LOG_chat10_2026-04-17.md` — phase-4e characterizer + 367-key spec.
- `AUDIT_LOG_chat11_2026-04-17.md` — registry library expansion.
- `AUDIT_LOG_chat11_5_2026-04-17.md` — C01j/C01l/C01m dispatch.

**Superseded handoff prompts:**
- `HANDOFF_PROMPT_chat9_2026-04-17.md`
- `HANDOFF_PROMPT_chat10_2026-04-17.md`
- `HANDOFF_PROMPT_chat11_2026-04-17.md`
- `HANDOFF_PROMPT_chat12_2026-04-17.md` — the instructions Chat 12
  executed.
- `HANDOFF_PROMPT_next_chat_2026-04-17.md`,
  `HANDOFF_PROMPT_next_chat_partB_2026-04-17.md` — earlier
  generic-name handoff drafts, superseded by the numbered ones.

**Dated rolling-doc snapshots:**
- `SESSION_SUMMARY_2026-04-17.md` — 59 KB snapshot covering chats
  3–11.5. A comprehensive historical narrative; most of its
  architectural content is now consolidated into the current-state
  `SESSION_SUMMARY.md` at root.
- `FIXES_APPLIED_2026-04-17.md` — 42 KB fix register for chats
  3–11.5. Full detail; the current-state `FIXES_APPLIED.md` at root
  summarises this era in one-liners and detailed entries apply only
  to the most recent chat.

## Archival convention

At the end of each chat, the previous chat's handoff and audit log
are moved here. The current chat's audit log stays at root until the
chat after it is done. The rolling docs (`SESSION_SUMMARY.md`,
`FIXES_APPLIED.md`) live exclusively at root and are updated in
place — dated snapshots are archived only at major milestones.
