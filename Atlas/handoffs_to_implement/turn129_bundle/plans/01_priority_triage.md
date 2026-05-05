# Priority triage — ordered list

The order matters. P1 unblocks P2. P2 unblocks the popgen workflow.
P3 reorganizes the UI without depending on P1/P2. P4 needs the
server (P1.5) to be stable. P5/P6 build on P3.

| # | Patch | Risk | Lines | Depends on |
|---|---|---|---|---|
| P1.1 | `const out` → `let out` in `collectPopstatsTracks` | trivial | 1 | none |
| P1.2 | Server URL: canonical localStorage key | low | ~30 | none |
| P1.3 | `/health` alias preservation | trivial | already done; lock in tests | none |
| P1.4 | `LDSplitReq` vs `FastLDReq` rename | low | ~10 | none |
| P1.5 | Server startup `POPSTATS_CONFIG` race | low | ~20 | none |
| P2.1 | Group dock build from `locked_labels` | medium | ~80 | P1.1 (popgen wiring needs it) |
| P3.1 | Two-row nav | medium | ~150 | none (replaces turn 128 pills) |
| P3.2 | Page consolidation (SV→Boundaries, Catalogue→Karyo) | medium | ~50 | P3.1 |
| P3.3 | JS file org audit + script tag fix | low | ~10 | none |
| P4.1 | Dosage chunk endpoint | medium | ~120 server + ~40 atlas | P1.5 stable |
| P4.2 | Theta candidate overlay endpoint | medium | ~100 server + ~40 atlas | P4.1 pattern, P1.5 stable |
| P5.1 | Boundaries + SV evidence skeleton | medium | ~200 (mostly empty-state) | P3.2 |
| P6.1 | Karyotype/Tier/Catalogue merge | medium | ~150 | P3.2 |

Total: roughly 13 patches, ~1100 lines of net code.

## Order of operations

**Session 1 (recommended):** P1.1 → P1.2 → P1.3 → P1.4 → P1.5 → P2.1
- All five P1s are atomic and don't interact.
- P2.1 depends on the group dock seeing fish_calls properly, which
  needs P1.1 (popgen wiring) so the popstats compute works after.
- Test: G key opens dock, candidate dims show up, Compute button
  hits server, server returns FST/dXY/theta, charts render.

**Session 2:** P3.1 → P3.2 → P3.3
- Nav restructure plus page consolidation. UI-only, low risk to data
  flow. P3.3 (JS file org) is independent and can land any time.

**Session 3:** P4.1 → P4.2
- Server data bridges. Once these land, the dosage heatmap and
  candidate θπ overlay come alive without needing R-side scripts.

**Session 4:** P5.1 → P6.1
- The two new combined pages. These are scaffolds — you'll iterate
  on each section over time.

## Specs (no code this turn)

- S1: SV evidence tables (5-table schema)
- S2: Indel slope / burden layer
- S3: Bayesian breakpoint scoring (future module)
- S4: Double-crossover / recombinant extension (future module)
- S5: SV interpretation rules (discipline doc)

These are deferred specifically because you said "for the other
stuff write to specs". They'll inform later sessions when you're
ready to wire SV data through.

## Decision points I need your confirmation on

Before I land any code in turn 130+:

1. **Page layout naming** — I'm proposing 6 main sections plus
   Overview/Help:
   `Overview | Discovery | Refine | Evidence | Compare | Output | Help`
   You used "Evidence" in your message but turn 128 had
   "popgen". Pick one. ("Evidence" is broader and might absorb
   future SV-specific evidence pages.)

2. **Catalogue location** — combining catalogue with
   karyotype/tier produces a heavy page. Alternative: keep
   catalogue as the **only** content of the Discovery section's
   landing tab (so Discovery's first tab is "the candidate
   table", and the four PCA tabs are subordinate). I think
   merging with karyo/tier is what you want, but flag if not.

3. **SV evidence merge** — your message has it merged with
   Boundaries. I agree. Just confirming: the page is named
   "Boundaries" (not "Boundaries + SV") and SV is a section
   within it, alongside the 9 other sections you listed.

4. **Two-row nav vs turn 128 pills** — turn 128 shipped pills that
   fold/unfold in-row. Your message rejects this:
   *"child buttons appeared to the right in the same row. This is
   bad."* — so P3.1 replaces turn 128's approach entirely. Confirm
   you want the in-row pills removed in favour of two true rows.
