# ADR-13. Q page (discovery) and Ancestry page (phase QC) — asymmetric naming

**Status:** locked April 30, 2026.
**Supersedes:** earlier sketches in this session that conflated ancestry-Q
discovery with phase-QC integration.

## Decision

Two new atlas pages, with deliberately asymmetric labeling:

| Tab label | Role | JSON consumed | Built when |
|---|---|---|---|
| **"Q"** | Per-stream discovery view of the ancestry-Q signal | `LG##_q.json` | Next, after dry-run validation |
| **"Ancestry"** | Phase-QC integration view (11-track shelf, cross-stream) | `LG##_shelf_qc.json` | After Q page lands and 4 discovery streams produce real data |

## What goes where

**"Q" page** — strictly the ancestry-Q discovery stream:
- dom_Q track (per-window dominant ancestry component)
- Δ12 track (single-scale)
- Multi-scale Δ12 track
- Local-PCA-on-Q sim_mat + |Z| profile + envelopes (when STEP_QR03 lands; deferred)

Mirrors pages 3 (GHSL) and 12 (θπ) structurally. Currently incomplete by design — STEP_QR01 + QR02 produce the tracks, but STEP_QR03 (the candidate detector) doesn't exist yet, so the page renders tracks without candidate envelopes overlaid.

**"Ancestry" page** — phase-QC / shelf QC integration view, per candidate:
- 11-track stacked composite at the candidate region (sim_mat density, Robust Z, SNPs/10kb, BEAGLE uncertainty, θπ mean+CV, dom_Q, Δ12, multi-Δ12, θπ_invgt for Hom1/Het/Hom2, Fst Hom1_Hom2)
- Per-sample karyotype call table — one row per sample, columns showing each method's classification at this candidate
- Cross-method concordance grid (the 6 pairwise contingency tables for dosage × GHSL × θπ × ancestry-Q)
- Candidate-bounded view; no chromosome-wide scrubber

The "Ancestry" label reflects the page's *purpose* (confronting ancestry confounding) not its literal contents. Most tracks shown are non-ancestry (Fst is structural, θπ is diversity, sim_mat is local-PCA), but the page exists because the QC battery's central question is "is this candidate confounded by ancestry / family structure?" — and the answer comes from looking at all 11 tracks together.

## Why asymmetric labels

Strict literal naming would call them "Ancestry-Q" and "Phase-QC integration." That's accurate but loses the reviewer-facing framing.

Reviewer reads "Ancestry" tab → immediately understands "this is where ancestry confounding gets confronted." That's the right framing for the page's job, even though the page contains many non-ancestry tracks.

If we labeled the QC page "Phase-QC integration," reviewers would have to learn what that means. If we labeled the Q-discovery page "Ancestry," they'd be confused why it shows only Q tracks and not Fst.

The current asymmetric labels prioritize reviewer comprehension over technical literal-mindedness. **Documented here so future-me doesn't "fix" the naming.**

## Why two pages, not a toggle

We considered double-clicking the tab to toggle between Q-discovery and shelf-QC views (analogous to the page 3 / 3-bis toggle). Rejected because:

1. **Different scopes.** Page 3 / 3-bis are two views of the same GHSL data. Q-discovery and shelf-QC are different data scopes (chromosome-wide per-stream vs candidate-focused cross-stream). Reusing the toggle metaphor for a structurally different relationship overloads the pattern.

2. **Different load dependencies.** Q page renders with just `LG##_q.json`. Shelf-QC needs all four discovery JSONs + the shelf QC JSON. Sharing a tab would require complex empty-state logic.

3. **Different positions in tab order.** Q page belongs in the discovery cluster (alongside pages 1, 3, 12). Shelf-QC belongs in the phase-4 / catalogue cluster (near page 4). They aren't visually adjacent in the tab bar — sharing a tab would force unnatural navigation.

4. **Different audiences in workflow.** Q page is for "is there ancestry signal at all on this chromosome?" Shelf-QC is for "for THIS candidate, do all streams agree?" Different questions, different answers, different pages.

## What this means for next session

When real data arrives:

1. **Build the Q page first.** It mirrors the existing pages 3 and 12 — known patterns, predictable scope. ~2 days of work for `export_q_to_json.R` + atlas page scaffold + 5 panel renderers (or fewer if QR03 is still pending).

2. **Build the Ancestry page second.** Cross-stream integration is harder — needs all 4 discovery JSONs available, plus a per-candidate aggregation step on the cluster (`export_shelf_qc_to_json.R`). ~3-5 days of work.

3. **Promote ancestry-Q to a "real" discovery stream when QR03 lands.** Until then, the Q page renders tracks-only with no candidate envelopes. That's an honest reflection of the upstream pipeline state — don't fake completeness with placeholder envelopes.

## Asymmetry vs other discovery streams

The four-stream discovery framing (ADR-1 + ADR-2 updated) becomes:

| Stream | Discovery page | Detector status | Hierarchy |
|---|---|---|---|
| Dosage | Page 1 | Existing | Primary |
| GHSL | Page 3 | STEP_C04b PASS-runs | Primary; |Z| + D17 secondary |
| θπ | Page 12 | STEP_TR_B \|Z\| + STEP_TR_D D17 | Both primary (no upstream calibrated detector) |
| Ancestry-Q | "Q" (NEW) | STEP_QR03 not built; QR02 produces tracks only | Tracks-only until QR03; no candidate detector |

The Q page being incomplete-by-design is acceptable because the architecture is honest about it. Better an honest tracks-only view than a faked discovery page with no real detector.

## Walk-back / reversal protection

If a future session proposes:

- "Let's combine Q and Ancestry into one page with a toggle" → reject; see "Why two pages, not a toggle" above
- "Let's rename Ancestry to Phase-QC" → reject; the asymmetric label is intentional reviewer framing
- "Let's add envelopes to the Q page using QR02 alone" → reject; without QR03 the envelopes are unjustified threshold-based detections, not validated candidates

The right walk-back trigger: **only if QR03 lands and produces calibrated candidates**, the Q page gets real envelopes. No earlier substitute.
