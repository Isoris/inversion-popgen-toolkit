# Obsolete design docs

Three design / how-to docs that self-marked themselves as outdated
with "outdated" in the filename. Kept for history; none of them
describes the current architecture.

| Doc | Why it's here | Replacement / current equivalent |
|---|---|---|
| `HOW_TO_RUN_outdated.md` | How-to-run for the v8.5 pipeline (retired). Self-marked OUTDATED 2026-04-17. | `docs/HANDOFF.md` (current end-of-session handoff) + the live phase READMEs under `inversion_modules/phase_N_*/README.md`. |
| `PIPELINE_v8.5_ARCHITECTURE_outdated.md` | Architecture of v8.5. Self-marked outdated. | `inversion_modules/README.md` (phase_1..7 architecture). |
| `SNAKE3_v3_ARCHITECTURE_outdated_now_v5.md` | Design draft for Snake 3 (GHSL) v3, dated 2026-04-07. GHSL is now at v6. | `inversion_modules/phase_2_discovery/2e_ghsl/README.md` covers the v5→v6 implementation story. The *conceptual* "why GHSL haplotype contrast works in a hatchery broodstock" rationale (§1 of the old doc, lines 30-60) is NOT yet migrated to the live README — a useful future doc improvement. |

These were previously at `inversion_modules/docs/`, which is now removed
(it contained only these 3 "outdated" files plus `LEGACY_FOLLOWUP_BRIDGE.md`
which moved to `inversion_modules/phase_5_followup/` next to the code
it maps).
