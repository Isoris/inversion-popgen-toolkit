# Atlas — Architecture Decision Record

A catalog of architectural decisions made during the GHSL/θπ Phase 2
build (April 2026 session). One section per decision. Each records:

- **What** was decided
- **Why** that choice over alternatives
- **What was rejected** (and why)
- **Reversals**, if any, with honest reasoning

This document is for future-you and reviewers asking "why did you do
it this way?" It does not document *implementation* — see
`RUNBOOK_produce_phase2_jsons.md` for that.

---

## ADR-1. Three-stream parallel discovery, not unified

**Decision:** Three independent inversion-discovery streams run in
parallel — dosage local PCA (page 1), GHSL local PCA (page 3), θπ
local PCA (page 12). Each produces its own candidate intervals.
Convergence across streams is evaluated post-hoc by interval overlap.

**Why:** The three streams measure different biological signals:
- Dosage local PCA sees genotype-frequency structure
- GHSL sees per-sample phased-haplotype divergence
- θπ sees per-sample nucleotide diversity

A region detected by all three is "convergent evidence" — much
stronger than any one alone. A region detected by only one is
informative in itself (sweep-only signal, recombination-suppression-
only signal, etc.). The Layer-A/B/C/D framework in §13 of the atlas
schema reflects this: same architecture, three independent layers.

**Rejected alternative:** Unified pipeline that combines signals
upstream (e.g. weighted sum of dosage + GHSL + θπ scores per
window). Discarded because it loses the orthogonality that makes
convergence interpretable. If you collapse the signals, you can't
say "GHSL finds X but dosage doesn't" — you only see the average.

---

## ADR-2. Detector hierarchy varies per stream (asymmetric)

**Decision:** The three envelope-detection methods (|Z|-threshold,
D17 cross-block scan, STEP_C04b PASS-runs) are present on each
stream where available, but their PRIMARY/SECONDARY status differs:

| Stream | \\|Z\\|-threshold | D17 cross-block | STEP_C04b PASS-runs |
|---|---|---|---|
| Dosage (page 1) | n/a | **PRIMARY** | n/a |
| GHSL (page 3) | secondary | secondary | **PRIMARY** |
| θπ (page 12) | **PRIMARY** | **PRIMARY** | n/a |

**Why:** The asymmetry isn't arbitrary — it reflects which streams
have an upstream production candidate detector calibrated on real
data. GHSL has STEP_C04b (PASS/WEAK/FAIL classifier with denominator-
confound mitigations baked in). θπ has nothing equivalent. So:

- Where a calibrated detector exists, that's the primary candidate
  set; local-PCA-derived envelopes are an orthogonal cross-check
- Where no calibrated detector exists, the local-PCA-derived
  envelopes ARE the primary set — they're the only set

**Reversal record:** Originally I removed envelope detection from
STEP_C04c entirely (turn 2 mid-session) on the grounds that
STEP_C04b already exists and adds parallel detection creates
ambiguity. Quentin pushed back: *"we can always have failed
envelopes its fine by default actually."* The fix wasn't to remove
the new envelopes; it was to ship them as SECONDARY rather than
competing with PRIMARY. Same script can produce envelopes or not;
the label (primary/secondary) is a routing decision, not a
methodology decision.

**Rejected alternative:** Always demote local-PCA envelopes to
secondary regardless of upstream availability. Discarded for θπ
because there's no primary to defer to — promoting the only
detector to "secondary" would leave the stream with no primary
candidates at all.

---

## ADR-3. Wrappers around D17, not modifications to it

**Decision:** STEP_C04d (GHSL) and STEP_TR_C (θπ) are thin Python/R
wrappers that prepare D17-compatible temp files (precomp.rds + 
sim_mat.rds), invoke `STEP_D17_multipass_L1_only_v7.R` and
`STEP_D17_multipass_L2_v8.R` as black boxes, then rename outputs to
namespaced TSVs (`<chr>_ghsl_d17L*.tsv`, `<chr>_theta_d17L*.tsv`).
**D17's 2,600 lines are not modified.**

**Why:** D17 is validated for the dosage path. Modifying it for
GHSL/θπ risks behavioral drift that affects the dosage pipeline.
Wrapping it preserves the dosage validation, applies the same
detector to the new streams without code duplication, and makes
the GHSL/θπ output schemas identical to dosage's (so downstream
consumers don't need to know which stream produced an envelope).

**Pragmatic escape hatch:** If D17's defaults need re-tuning for a
new stream, pass `--d17_l1_args "--boundary_grow_W_pct 0.005,0.02,0.05"`
to the wrapper. Wrapper forwards extra args to the D17 invocation
without touching D17.

**Rejected alternative:** Generalize D17 itself with a `--stream
dosage|ghsl|theta` flag. Discarded because (a) the stream-specific
parts are at the wrapper level, not inside D17, (b) any change to
D17 needs to re-validate the dosage path, (c) the wrapper approach
ships in 2 days; the generalization approach ships in 2 weeks.

---

## ADR-4. D17 complexity is O(N²), not O(N³)

**Decision:** D17 wrapper performance budget assumes O(N²) total
work. No optimization of D17 itself was attempted.

**Why:** I initially called D17 cubic. Empirical benchmark (turn 6)
proved this wrong:

- Step 1 (diagonal mean/sd precompute): N iterations × O(N-d) per
  diagonal = N²/2 ops total — quadratic.
- Step 2 (cross-block scoring): O(N) iterations × O(W²) per step at
  fixed W = O(N) at constant W — linear.
- Grow validator: n_peaks (~5-50 per chromosome, scales with biology
  not N) × Σ W² across grow series. When W series is percent-of-N,
  this scales as O(N²).

**Empirical wall times in R (extrapolated from numpy benchmark):**

| N | Step 1 | Grow validator (default 6 W's) | Total |
|---|---|---|---|
| 4,300 (GHSL on LG28) | ~3-6s | ~2-4s | ~5-10s |
| 16,500 (θπ fine grid) | ~40-80s | ~16-36s | ~60-120s |

**Why I'm recording this:** Future-me might re-investigate
performance and re-derive the same wrong cubic estimate. Pinning
"D17 is O(N²) at ≤2 minutes per chromosome at fine θπ grid" 
prevents that loop.

---

## ADR-5. K=3 GHSL palette deliberately mirrors dosage K=3

**Decision:** `COLORMAP_KSTRIPE[3]` reuses dosage's GROUP_COLORS
exactly: blue (HOM_REF), grey (HET), amber (HOM_INV). Higher K
values extend the palette with intermediate hues.

**Why:** Cross-method visual concordance. A sample shown blue under
dosage-K=3 AND blue under ghsl-K=3 is the same band in both
methods → convergent evidence visible at a glance, no mental
remapping. The whole point of having three orthogonal streams is
to spot agreement and disagreement; making the palettes match for
the most-likely band-count (K=3 = HOM/HET/HOM) lets the eye do that
work directly.

**Rejected alternative:** Distinct palettes per stream so streams
look different in the UI. Discarded because the goal is comparison,
not differentiation — different colors per stream would force the
user to mentally remap "blue here means stripe-1; blue there means
something else."

---

## ADR-6. Sim_mat shipped as upper-triangle-packed flat array

**Decision:** GHSL `ghsl_local_pca.sim_mat` is emitted as a flat
Float32 array of length n*(n+1)/2 with `sim_mat_format:
"upper_triangle_float32"` and `sim_mat_n` carrying the matrix size.
Atlas reader symmetrizes on load and caches the full n×n on
`state._ghSimFull`.

**Why:** The lower triangle is symmetric to the upper — storing
both wastes half the bytes. At GHSL N=4,300 on LG28, full dense
sim_mat is ~74 MB raw; upper-triangle packing brings it to ~37 MB
raw. Below the 120 MB JSON budget either way, but lighter is
better for browser load time.

**Why not banded or quantized:** GHSL fits without tricks. θπ at
fine grid (N=16,500) is ~1.1 GB raw and DOES need banding or
quantization — that's a separate decision (ADR-7) not yet made.

**Trade-off accepted:** The atlas decoder pays a one-time
symmetrization cost on JSON load. At N=4,300 this is ~18M
assignments — fast (~50ms). Cached, so subsequent renders pay zero.

---

## ADR-7. θπ sim_mat sizing — DEFERRED until real-data evidence

**Decision:** No commitment yet. Three options on the table:

- **A. Coarse grid:** sim_mat at `win50000.step10000` (N=3,300
  on LG28) → 44 MB raw, 10 MB gzip — fits without tricks.
  Loses spatial resolution.
- **B. Banded:** only ±100 windows around diagonal at fine grid
  (N=16,500) → 13 MB raw, 3 MB gzip. Captures local similarity;
  loses long-range structure (which is usually noise floor).
- **C. Int8 quantized:** full dense at fine grid → 280 MB raw,
  70 MB gzip. Borderline fit, decode-on-render overhead.

**Why deferred:** The right choice depends on what fine-grid θπ
sim_mats actually look like. If long-range cells turn out to carry
biological signal (cross-chromosome correlation patterns), banding
is wrong. If they're noise, banding is the cleanest answer. **Cannot
decide without seeing real data.**

**Status:** STEP_TR_B v5 retrofit (which adds sim_mat to θπ) is
still pending. The wrapper STEP_TR_C currently computes sim_mat
in R memory and passes it to D17 as RDS — never serialized to
JSON. So this decision is only blocking the page-12 panel
renderers, not the cluster-side detection.

---

## ADR-8. Atlas patches as one-shot Python mutators, not source-of-truth

**Decision:** The four atlas patcher scripts (`apply_page3_scaffold.py`,
`apply_ghsl_color_mode.py`, `apply_page3_panels.py`,
`apply_enrichment_merge_fix.py`) are one-shot mutators of
`pca_scrubber_v4.html`. Once applied, the patched HTML is the
canonical artifact. The patchers are kept in `_archive/` for
reference but are not part of daily development.

**Why:** The atlas HTML is hand-edited as the source of truth.
Ratcheting it through scripted patches is fragile — one syntax
change in the upstream HTML and a regex in a patcher fails. Better
to commit the patched HTML directly and move on.

**Why ship the patchers at all:** Reproducibility. If something
breaks in the patched HTML and we need to rebuild from a fresh
upstream base, the patchers document exactly what the changes were.
Also useful for code review — diff the patcher's old vs new strings
to see the shape of the change.

**Idempotency requirement:** All four patchers are idempotent
(re-running on already-patched HTML is a no-op with skip messages).
This was tested explicitly during development and is a hard
requirement — without it, a re-run could corrupt the atlas.

---

## ADR-9. Field naming — short forms in user-facing files

**Decision:** Local atlas folder uses short stream-suffixed names:
`LG28_dosage.json`, `LG28_ghsl.json`, `LG28_theta.json`. Cluster
keeps verbose phase-numbered names: `C_gar_LG28_phase2_ghsl.json`.
The `sync_jsons.sh` script bridges the two without renaming on
the cluster.

**Why:** Two audiences with different needs:
- The cluster pipeline cares about phase numbers, species, and
  pipeline stage. Long names are correct there.
- Daily atlas use cares about chromosome and stream. The user's
  brain operates on `LG28_ghsl` not `C_gar_LG28_phase2_ghsl`.

Trying to satisfy both with one name fails one. So: keep both,
bridge at the boundary.

**Why no chromosome-internal subfolders:** The user wanted simple,
flat-ish layout (the "no 999 folders" rule). At 28 chromosomes ×
3 streams = 84 files in `json/`, `ls` is still readable. Adding
`json/LG28/` with three files inside adds nesting without making
anything easier to find.

---

## ADR-10. Documentation lives next to the thing, not separately

**Decision:** Architecture docs (this file, `RUNBOOK_*.md`,
`schema_v2_addendum_*.md`) live in `_scripts/` next to the code
they describe, not in a separate `docs/` folder.

**Why:** The 1:1 mapping between code and its decision record
matters more than documentation completeness. When you open
STEP_C04c and want to know why it has |Z|-threshold envelopes
shipped as secondary, the answer should be findable in the same
folder, not in a documentation tree two levels deep.

**Why not inline in code comments:** Decision records benefit from
prose space (explaining WHY, recording reversals, cross-referencing
related decisions) that doesn't fit in code comments without making
the file unreadable. Best of both worlds: code stays readable, docs
stay close.

**`_archive/` for superseded docs:** Earlier-version decision docs
(e.g. `phase2_ghsl_arrangement_v0.md` superseded by v1) move to
`_archive/`. They're useful for archaeology if a reviewer asks "why
NOT the v0 approach?" but don't clutter daily navigation.

---

## ADR-11. Test as structural assertion, not behavioral test

**Decision:** All tests are *structural* — they verify the right
files exist, contain the right field names in the right places,
have the right CLI flags wired up. They do NOT exercise actual
runtime behavior with real data.

**Why:** I have no R installed, no LANTA access, no real data. A
"behavioral" test would require running the R pipeline against
fixtures, which I can't do. Structural tests catch the bugs I CAN
introduce: missing fields, broken anchors in patcher scripts,
wrong filename patterns, function-ordering bugs in R. They don't
catch logic bugs.

**Trade-off accepted:** The structural tests give confidence that
the pieces fit together; only a real run on LG28 confirms they
work. The runbook's dry-run section is the smoke test.

**Why this matters for documentation:** A reader sees "48/48 tests
passing" and might assume runtime correctness was verified. It
wasn't. Recording this here so the claim isn't over-interpreted.

---

## ADR-12. The "include cheap things by default" principle

**Decision:** When deciding whether to ship an optional output
field/layer/computation, default to YES if the cost is small.

**Why:** I made the same kind of mistake repeatedly during this
session — when there was tension about whether to ship something,
I defaulted to "remove it for safety" rather than "ship it and let
the user decide." Quentin's pushback on the GHSL secondary
envelopes (turn 2): *"we can always have failed envelopes its fine
by default."* This is correct as a general principle:

- Cost of overinclusion: a few KB of unused JSON, slight cognitive
  load when reading the schema
- Cost of underinclusion: have to come back later and re-derive,
  re-run the cluster pipeline, update the schema, re-emit JSONs

The asymmetry is severe. **Default to ship.** Only remove when
there's a real cost (memory, performance, schema confusion).

---

## How to read this file

If a future decision conflicts with what's recorded here, that's
fine — write a new ADR superseding the old one and link them. Don't
edit existing ADRs in place; they're the historical record. The
"Reversal record" sections within each ADR show the right pattern.

If you're a reviewer asking "why did you do X?", search for X in
this file first. If it's not here, that's a documentation gap —
file it as a question and we'll write the missing ADR.

---

*Generated end of April 2026 session. 12 ADRs covering the GHSL
page-3 build, θπ retrofit decisions, atlas patcher architecture, and
file-layout choices.*
