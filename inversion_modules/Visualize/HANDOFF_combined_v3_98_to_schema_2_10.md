# Handoff — combined (v3.98 + STEP_M06 + schema 2.10) — 2026-04-28

This handoff covers everything shipped across the last three increments
in one document, so a future chat can resume from a single read.
Detailed per-increment handoffs (`HANDOFF_v3_98.md`,
`HANDOFF_v3_99_step_m06.md`) remain available for reference.

## State at end of this chat

| Artefact | Version | Size | Status |
|---|---|---|---|
| `pca_scrubber_v3.html` | **v3.98** | 1052 KB, 23,147 lines | shipped |
| `SCHEMA_V2.md` | **v2.10** | 2,179 lines | shipped (planned-status for §10.1) |
| `STEP_M04_emit_dosage_chunks.R` | v3.95 | 24 KB | shipped (carried) |
| `STEP_M05_emit_step29_layers.R` | v3.95 | 16 KB | shipped (carried) |
| `STEP_M06_emit_boundary_evidence.R` | **post-v3.98** | 27 KB, 612 lines | shipped (this chat) |

**Test coverage end-of-chat:**
- 615/615 unit tests passing (39 new in v3.98)
- 77/77 e2e checks (boundaries lifecycle, TSV roundtrip, etc.)
- 32/32 R smoke checks for STEP_M06
- 36/36 Node cross-checks for STEP_M06 (against real `_buildBoundaryTrackScores`)
- 6/6 Monte Carlo verification checks for NULL_CONCORD

## What shipped this chat session

### Increment 1 — v3.98 NULL_CONCORD K=6 fix + v3.80 UI text rewrite

**NULL_CONCORD K=6 fix.** The contingency-table concord display
compares two clusterings' agreement rate against a random-baseline
expectation. Pre-v3.98 the table only had `K∈{2,3,4,5}` so `K=6` fell
through to a `0.30` fallback — but Monte Carlo at N=226 with 5000
trials gives `0.24`. K=6 sub-band partitions that were genuinely well
above baseline got reported as "near baseline" by ~6 percentage points.

**Fix:**

```js
// Before:
const NULL_CONCORD = { 2: 0.55, 3: 0.39, 4: 0.30, 5: 0.25 };
const nullConcord = NULL_CONCORD[K] || 0.30;

// After:
const NULL_CONCORD = { 2: 0.55, 3: 0.39, 4: 0.30, 5: 0.25, 6: 0.24 };
const nullConcord = NULL_CONCORD[K] || 0.24;
```

K=2..5 left at legacy values (within ~1.5 pp of empirical at N=226).
The comment now cites Monte Carlo / N=226 / 5000 trials / N-dependence.
`verify_null_concord_v398.js` re-derives values any time you want to audit.

**v3.80 UI text rewrite — 7 surgical edits.**

1. Sidebar candidate-mode tooltip — "draft candidate" → "draft interval; promote to add to candidate list"
2. Page 2 candidate-focus empty-state — explicit lifecycle (promote → provisional → mark confirmed) + tab number fix (page 3 → page 4 since v3.97 boundaries took tab 3)
3. Page 5 karyotype empty-state — same vocabulary + tab fix
4. Catalogue button title — "Define an inversion candidate" → "Promote the selected interval(s) to provisional candidates"
5. L3 band-continuity scope tooltips (×2) — "draft candidate-mode draft span" → "interval being assembled in candidate mode"
6. Help-page tabs table — candidate-focus + karyotype rows now describe the lifecycle
7. **NEW Help-page Vocabulary section** — 5 rows (L1/L2/L3, Interval, Candidate, Regime, Boundary zone) defining each with the v3.80 lifecycle baked in

39 new tests in `test_v398.js`.

### Increment 2 — STEP_M06_emit_boundary_evidence.R (R-side, post-v3.98)

Per-chromosome emit of the optional `boundary_evidence` layer for the
v3.97 boundary zone refinement module. Track-by-track availability —
each track independently optional based on inputs:

| Track | Required input | Algorithm |
|---|---|---|
| `fst` | `--dosage` + `--fish_regimes` | Hudson FST per window between K=3 regimes (homo1 vs homo2); ratio of mean numerator over mean denominator |
| `theta_pi_homo1`/`het`/`homo2` | `--dosage` + `--fish_regimes` | Per-window per-regime θπ with sample-size correction `n/(n-1)` |
| `discordant_pair_pileup` | `--discordant_bed` | 4-col BED (chrom/start/end/count) → midpoint binning into scan windows |
| `sv_anchors` | `--sv_vcf` (1+ files) | Parses DELLY (CT=) + Manta (INV3/INV5) for SVTYPE in {INV, DEL, BND} |

**Production usage:**

```bash
Rscript STEP_M06_emit_boundary_evidence.R \
  --candidates_registry candidates_registry.tsv \
  --chrom C_gar_LG28 \
  --out_dir /path/to/scrubber/data/ \
  --dosage C_gar_LG28.dos.tsv.gz \
  --fish_regimes fish_regime_calls.tsv \
  --samples sample_list.txt \
  --sv_vcf delly.inv.vcf manta.inv.vcf \
  --discordant_bed discordant_pairs.bed
```

Drop the resulting `<chrom>_boundary_evidence.json` into the scrubber's
enrichment loader. Open page 3 boundaries → select a candidate →
auto-propose now uses 11 contributing tracks instead of 7 (the
4 boundary_evidence-derived ones add ~15% of total weight).

**Bugs fixed from in-progress draft:** `parse_sv_vcf` was forward-
referenced before its definition (R doesn't hoist non-block-scope
function calls in script-mode evaluation). Reorganized to def-then-use.

**Validated end-to-end on synthetic data:**
- FST: 0.74 at left edge (11.95 Mb), 0.79 at right edge (14.05 Mb), 0.32 globally
- Discordant pile: 52 at left, 82 at right (background ~2)
- SV anchors: 3 in scan range (DELLY INV ×2 with CT=3to3 / 5to5, Manta_INV5 ×1); off-target Manta_INV3 at 20 Mb correctly filtered; off-chrom call correctly filtered
- Scrubber's actual `_buildBoundaryTrackScores` (loaded via `vm.createContext`) produces all 4 expected output tracks with signal landing at the correct windows

### Increment 3 — Schema 2.9 → 2.10 (within-regime subclustering)

**Status reconciliation (2.9 → 2.10):** Carried forward from previous
chats — flipped 3 layers from "planned" to "shipped" in the layer
registry table:
- `dosage_chunks`, `candidate_sample_coherence`, `candidate_marker_polarity` (R-side STEP_M04 + STEP_M05 from v3.95)
- `boundary_evidence` (R-side STEP_M06 from this chat)

No new contract surface in 2.9 — purely metadata.

**Within-regime subclustering contract (2.10):** New §10.1 sub-section
in the dosage heatmap module. The contract:

> When sub-structure exists *inside* a major regime, running density-
> based clustering naively on raw `(u, v)` across all samples conflates
> within-regime signal with the top-level HOMO_1/HET/HOMO_2 split. The
> major-regime axis dominates distances and the minor within-regime
> axis gets washed out.

**5-step algorithm (per regime separately):**
1. Extract `(u, v)` for samples in that regime
2. SVD-rotate within-regime `(u, v)` → `(u_rot, v_rot)`
3. Z-score standardize `u_rot` and `v_rot` independently
4. DBSCAN (default `eps=0.5`, `min_samples=3`) or HDBSCAN
5. Label as `<REGIME>_sub<N>` (or `<REGIME>_sub_noise` for DBSCAN noise points)

**`cluster_mode` parameter — four values:**

| Mode | Default? | Description |
|---|---|---|
| `none` | **default** | Subclustering not run. Preserves current scrubber behaviour — heatmap K=3 track shows existing solid-color cells. |
| `rotated_within_regime` | opt-in | The 5-step contract. Production mode once calibrated. |
| `raw_uv` | diagnostic | DBSCAN on raw `(u, v)` across all samples. Almost always recovers the major-regime split — diagnostic only. |
| `one_dimensional_rotated_axis` | diagnostic | Per-regime split + 1D `u_rot` only. Fallback for low-coverage where `v_rot` is noise-dominated. |

**Why default = `none`.** Subclustering is a *new* analytical contract.
Until calibrated on the 226-sample LG28 cohort (per-regime sample-size
lower bounds, eps tuning, sub-stripe color palette), the safe default
is to leave existing heatmap behaviour unchanged. Default may flip to
`rotated_within_regime` in a future schema bump after validation —
that change is explicitly out of scope for 2.10.

**Critical guardrail (in-spec, repeated three times in §10.1):**
> Subclusters are DIAGNOSTIC within-regime haplotype backgrounds, NOT
> paper-level inversion states. Manuscript text should never report
> "`HOMO_1_sub2` carries 12 samples" as a primary finding — that is a
> a description of within-regime broodline structure, not of the
> inversion itself.

**JSON contract.** Subcluster labels ride on the existing
`candidate_sample_coherence` layer (no new layer). Three new optional
columns: `coarse_group_refined`, `cluster_mode`, `u_rot`, `v_rot`.

**Heatmap rendering.** When `cluster_mode ∈ {rotated_within_regime,
one_dimensional_rotated_axis}` AND labels are present, the K=3 PCA
group track gains an inner sub-stripe (one sub-color per `_sub<N>`,
monotone within each major regime's hue). Otherwise unchanged.

**Status: planned (no scrubber implementation, no R-side emit yet).**
The contract is locked; downstream work is:
- Scrubber-side implementation of the 5-step algorithm (UI toggle, recompute button, sub-stripe rendering)
- R-side emit modification — extend `STEP_M05_emit_step29_layers.R` to add the 3 new columns, OR a new `STEP_M07_compute_subclusters.R` that does the rotation+DBSCAN cluster-side and emits to the same TSV pattern STEP29 already produces

## What's outstanding before v4.0

Three concrete blockers, ranked by gating-importance:

**1. Real-data validation on the 226-sample LG28 cohort. *Gating item.***
Everything from v3.93 onward (het-shape, dosage heatmap, hover tooltips,
boundary refinement, K=6 fix, within-regime subclustering contract)
needs real data. Synthetic inputs always saturate cleanly; ergonomic +
calibration issues only become legible on actual cohort data. **This
is a you-with-eyes-on-data task, not a code task.**

**2. Three parked R-side band-diagnostics emits:**
- `STEP_T06_emit_theta_pi_panel.R` — unlocks θπ pillar (v3.92)
- `STEP_R12_emit_roh_intervals.R` — unlocks ROH pillar (v3.92)
- `STEP_R13_emit_sample_froh.R` — unlocks FROH pillar (v3.92)

Without these, 6 of the 12 band-diagnostics metrics render as `?`,
and the `theta_pi_shift`/`ROH_confounded`/`FROH_high` flags can
never fire. Each follows the M04/M05/M06 pattern (~600 lines + smoke
+ cross-check) and takes one chat each.

**3. Marker module / phase-13.** Page 10 hides itself entirely until
`marker_panel_summary` JSONs ship (STEP_M01/M02/M03). Separate
deliverable from the scrubber.

Plus a final clean-up pass: dead-code removal, any v3.x tech debt,
and a real `CHANGELOG.md` if you want one for the manuscript Methods
section.

## Suggested entry points for next chat

In rough priority order:

1. **Real-data validation on LG28** — by far the most valuable. Drop
   STEP_M04 + STEP_M05 + STEP_M06 outputs into the scrubber, navigate
   to pages 2 (candidate-focus + dosage heatmap), 3 (boundaries), and
   the L3 cursor strip. Surface ergonomic issues, propose v3.99+ fixes.
2. **Implement schema 2.10 subclustering** — scrubber-side 5-step
   algorithm + UI toggle + sub-stripe rendering, plus a small R-side
   extension to `STEP_M05` to emit the new columns. Probably one chat.
3. **Three parked R-side emits** (T06, R12, R13) — three chats, each
   following the M04/M05/M06 pattern.

## File staging

All accumulated outputs are at `/mnt/user-data/outputs/`. A combined
tarball `scrubber_v3_98_post_m06_bundle.tar.gz` packages everything
needed for next-chat reconstruction (current scrubber + schema +
all R-side emit scripts + active handoffs + tests).

## User working style notes

(Unchanged.)

"Continue" = proceed with next agreed item. Pragmatic, multi-session
parallel chats. Strongly prefers compact non-fragmented pipelines,
original tested code verbatim, single central config file per module,
targeted fixes over rewrites, deliverables as tarballs ready to extract.

Three catfish cohorts that must NEVER be conflated:
1. F1 hybrid (*C. gariepinus* × *C. macrocephalus*) — genome assembly paper only
2. 226-sample pure *C. gariepinus* hatchery cohort — current inversion work, "MS_Inversions_North_african_catfish" — K clusters reflect broodline structure, NOT species admixture
3. Pure *C. macrocephalus* wild cohort — future paper

User's full name is Quentin Andres. Never invent a surname.
