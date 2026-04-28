# Handoff — R-side boundary_evidence unblock (post-v3.98, 2026-04-28)

## Where we are

**Scrubber head:** `pca_scrubber_v3.html` v3.98 (1052 KB, 23,147 lines) — unchanged this chat.
**Schema:** `SCHEMA_V2.md` v2.8 — unchanged.
**Tests:** 615/615 still passing + 77 e2e checks (no scrubber-side changes).

This chat shipped one R-side emit script that unblocks v3.97's boundary
zone refinement module for real-data validation on the 226-sample LG28
cohort. Once this script runs on real data and the JSON is loaded into
the scrubber, v3.97's auto-propose algorithm picks up four additional
evidence tracks (`fst_edge`, `theta_pi_step`, `discordant_pile`,
`sv_anchor`) automatically — no scrubber changes needed.

This chat was scoped explicitly to the deferred R-side work the user
parked at the start of v3.98. The in-progress draft from earlier this
session had a `parse_sv_vcf` forward-reference bug; that's now fixed.

## Work shipped this chat

### `STEP_M06_emit_boundary_evidence.R` (new, 612 lines)

Per-chromosome emit of the optional `boundary_evidence` layer for the
scrubber's v3.97 boundary zone refinement module (schema 2.8 §12). One
output file per chromosome:

```
<out_dir>/<chrom>/<chrom>_boundary_evidence.json
```

The JSON shape matches schema §12 exactly:

```jsonc
{
  "schema_version": 2,
  "_layers_present": ["boundary_evidence"],
  "boundary_evidence": [
    {
      "candidate_id": 1,
      "chrom": "C_gar_LG28",
      "scan_start_bp": 10500000,        // candidate ± scan_radius_bp, clamped
      "scan_end_bp":   15500000,
      "scan_window_bp": 5000,
      "tracks": {
        "fst": [...],                    // length = (end-start)/win
        "theta_pi_homo1": [...],
        "theta_pi_het":   [...],
        "theta_pi_homo2": [...],
        "discordant_pair_pileup": [...],
        "sv_anchors": [
          { "kind": "DELLY_INV", "ct": "3to3", "pos_bp": 11950000, "qual": 80 },
          { "kind": "Manta_INV5", "pos_bp": 14050000, "qual": 95 }
        ]
      }
    }
  ]
}
```

### Track-by-track availability — all four are independently optional

The scrubber's `_buildBoundaryTrackScores` at line 6733 reads
`data.boundary_evidence`, finds the matching candidate by `candidate_id`,
and consumes whatever tracks are present. Missing tracks are silently
skipped. STEP_M06 mirrors this: each track is emitted only when its
input source is provided.

| Track | Required input | Algorithm |
|---|---|---|
| `fst` | `--dosage` + `--fish_regimes` | Hudson FST per window between K=3 regimes (homo1 vs homo2); ratio of mean numerator over mean denominator |
| `theta_pi_homo1`/`het`/`homo2` | `--dosage` + `--fish_regimes` | Per-window per-regime θπ with sample-size correction `n/(n-1)`; NA-tolerant |
| `discordant_pair_pileup` | `--discordant_bed` | 4-col BED (chrom/start/end/count) → midpoint binning into scan windows |
| `sv_anchors` | `--sv_vcf` (1+ files) | Parses DELLY (CT=) + Manta (INV3/INV5) for SVTYPE in {INV, DEL, BND} |

### Key design decisions

- **scan_window_bp default 5 kb.** Matches schema §12. Coarser than the
  L1 PCA window step (5–10 kb) so the scrubber doesn't re-blur an
  already-smoothed input.
- **scan_radius_bp default 1.5 Mb.** Matches the scrubber's default
  scan_radius_bp. Schema's "huge candidate" expansion rule is honored:
  for candidates with `span > 3 Mb`, radius = max(1.5 Mb, span × 0.5).
- **NA encoding via -1.** Dosage TSV files often use `-1` as the NA
  sentinel (matches `selectTopMarkers` line 4949). M06 maps `-1 → NA`
  in the dosage matrix; subsequent FST/θπ computations skip NA cells
  per marker.
- **Hudson FST formula uses ratio-of-averages.** Per Hudson 1992: window
  FST = `mean(numerator)/mean(denominator)` rather than `mean(per-marker
  FST)`. This is the correct unbiased estimator for low-coverage data.
  Per-marker numerator: `(pA - pB)² - pA(1-pA)/(nA-1) - pB(1-pB)/(nB-1)`;
  denominator: `pA(1-pB) + pB(1-pA)`.
- **θπ uses sample-size correction.** `2 · p · (1-p) · n/(n-1)` where
  `p` is the regime's allele frequency at the marker and `n` is the
  sample count after NA-drop. Without the `n/(n-1)` correction, θπ
  under-estimates with small regime sizes.
- **Marker pos_bp parsed from multiple TSV formats.** The dosage TSV
  may have `pos_bp`/`pos` columns explicitly, or marker IDs in
  `<chrom>_<pos>` / `<chrom>:<pos>` form. STEP_M06 detects either.
- **SV VCF parser handles DELLY + Manta together.** Reads the standard
  8-column format, filters to `SVTYPE in {INV, DEL, BND}`, extracts
  `CT=` (DELLY connection type for inversions) and detects Manta
  `INV3`/`INV5` from the INFO field. SV calls outside scan range are
  filtered per-candidate; off-chrom calls filtered up-front.
- **Optional inputs fail soft, not hard.** Missing VCF → warn + skip
  that file. Missing BED → warn + skip the track. Missing
  fish_regimes/dosage → skip FST + θπ tracks. Only the
  `candidates_registry` is hard-required.
- **Empty `tracks{}` on a candidate → skip the candidate entry entirely.**
  Emitting an empty `boundary_evidence` row would still register the
  layer with the scrubber but provide no signal — better to omit.

### Required inputs

- `--candidates_registry <file.tsv>` — v3.97 candidate registry TSV
  (must include `candidate_id`, `chrom`, `start_bp`, `end_bp`).
- `--chrom <chrom>` — restricts to this chromosome.
- `--out_dir <dir>` — writes `<out_dir>/<chrom>/<chrom>_boundary_evidence.json`.

### Optional inputs (each enables different tracks)

- `--dosage` + `--fish_regimes` → enables `fst`, `theta_pi_homo1`,
  `theta_pi_het`, `theta_pi_homo2`
- `--sv_vcf vcf1 vcf2 ...` → enables `sv_anchors` (any number of VCFs)
- `--discordant_bed` → enables `discordant_pair_pileup`
- `--samples` → cohort sample order (defaults to dosage TSV header order)
- `--scan_radius_bp` (default 1500000)
- `--scan_window_bp` (default 5000)

### Cross-validation against the real scrubber

Both smoke + cross-check pass end-to-end:

- **`smoke_test_M06.R`** — 32/32 R-level checks. Builds a synthetic
  C_gar_LG28 candidate at 12-14 Mb with engineered FST / discordant /
  SV signal at the boundaries (allele-freq differences between regimes
  at left/right edges, hetero-elevated body markers, BED spikes, INV
  calls in two VCFs). Validates JSON shape, track lengths, signal
  landing, NA handling, off-target filtering.

- **`cross_check_M06.js`** — 36/36 Node-level checks. Loads M06's
  output JSON + the REAL `_buildBoundaryTrackScores` from
  `pca_scrubber_v3.html` via `vm.createContext`. Confirms:
  - Layer detection passes (Array.isArray + first row has candidate_id + tracks)
  - Track shape contract (lengths consistent with scan_start/end/window)
  - `fst_edge`, `theta_pi_step`, `discordant_pile`, `sv_anchor`
    all populate as Float64Array (duck-typed across VM boundary)
  - FST signal lands at the right windows (0.74 at left edge,
    0.79 at right edge — matching the synthetic edge signal)
  - `sv_anchor` fires at the exact windows containing our 3 INV calls
  - Empty `boundary_evidence` array → graceful (no crash, no tracks)
  - Missing candidate_id → graceful (no crash, no tracks)

This means: **when STEP_M06 runs on real LG28 data and the JSON is
loaded into the scrubber, v3.97's auto-propose picks up the four
additional tracks automatically.** Visiting page 3 boundaries with
a real candidate selected will use 11 contributing tracks (the 7
already-wired ones plus the 4 new) instead of just 7.

## Suggested production run

```bash
# Per-chromosome boundary_evidence (do for LG28 first, then iterate)
Rscript STEP_M06_emit_boundary_evidence.R \
  --candidates_registry  /scratch/lt200308-agbsci/.../candidates_registry.tsv \
  --chrom                C_gar_LG28 \
  --out_dir              /path/to/scrubber/data/ \
  --dosage               /scratch/lt200308-agbsci/.../C_gar_LG28.dos.tsv.gz \
  --fish_regimes         /path/to/.../C_gar_LG28_fish_regime_calls.tsv \
  --samples              /scratch/lt200308-agbsci/.../sample_list.txt \
  --sv_vcf               /scratch/lt200308-agbsci/.../module_4d_delly/C_gar_LG28.delly.inv.vcf \
                         /scratch/lt200308-agbsci/.../module_4h_manta/C_gar_LG28.manta.inv.vcf \
  --discordant_bed       /scratch/lt200308-agbsci/.../C_gar_LG28.discordant.bed \
  --scan_radius_bp       1500000 \
  --scan_window_bp       5000
```

Drop the resulting JSON into the scrubber's enrichment loader. Open
page 3 boundaries → select a promoted candidate → click auto-propose
→ the 4 additional tracks contribute to the combined score.

## Sizing on the 226-sample LG28 cohort (estimate)

LG28 ≈ 30 Mb. Per candidate (typical span 1-3 Mb + 1.5 Mb scan radius
each side) → scan range 4-6 Mb → 800-1200 windows of 5 kb. Per
candidate JSON entry:

- 4 numeric tracks × 1000 windows × 6dp + comma overhead ≈ 35 KB per cand
- SV anchors: 5-15 entries per candidate ≈ 1-2 KB

For ~50 candidates per chrom (LG28 is upper bound) → ~2 MB JSON. Light
compared to the dense GHSL panel (50-300 MB per chrom). LRU cache caps
not relevant here.

Computational cost is dominated by the FST/θπ inner loops over
markers × samples per scan window. For 50 cands × 1000 windows × 100
markers/window × 226 samples ≈ 1.1B cell reads per chrom. Should run
in ~5-15 min single-threaded in R; SLURM array per chrom (28 jobs) is
trivial parallelism if needed.

## What this unblocks (downstream)

1. **v3.97 boundaries page on real data** — auto-propose now uses 11
   contributing tracks instead of 7. The four additional tracks have
   weights summing to 0.15 (fst_edge: 0.05, theta_pi_step: 0.04,
   discordant_pile: 0.04, sv_anchor: 0.02) — modest individually but
   combining into ~15% of the total weight budget.
2. **Stronger support_class promotions** — when fst_edge / theta_pi_step
   step coherently with the existing pca_drop / ghsl_step / band_continuity_drop,
   the support_class can hit `strong` (≥3 contributing tracks at score
   ≥ 0.60) more reliably.
3. **SV-supported breakpoint_status** — with `sv_anchors` populated,
   the user can promote `breakpoint_status` from `boundary_zone_only`
   to `SV_supported` after reviewing whether DELLY/Manta calls fall
   inside the proposed zones. Currently a manual user action, but
   the data is now visible.

## Pipeline run-order notes

- **STEP_M06 is independent of STEP_M04 and STEP_M05.** Different
  candidate enrichment layers, can run in parallel.
- **Inputs must come from the same candidate set** — `candidates_registry`,
  `fish_regimes`, dosage, and SV VCFs should all reflect the same
  promoted-candidate state. If you re-promote candidates on the scrubber,
  re-run M06 (cheap) before re-loading the layer.
- **VCFs need not be filtered to in-scan ranges in advance.** STEP_M06
  filters per-candidate at read time. Pass the full chromosome VCF.
- **STEP_M06 does not require the marker catalogue.** Unlike STEP_M04
  (where partial-coverage marker_catalogue is silently dropped),
  STEP_M06 doesn't read it.

## Outstanding R-side from prior handoffs (still parked)

These don't block boundary_evidence work; listed for completeness:

- `STEP_T06_emit_theta_pi_panel.R` — unlocks θπ pillar (v3.92 band-diag).
  Note: STEP_M06 emits **per-candidate per-regime** θπ for boundary
  scoring, NOT the **per-window genome-wide** θπ that band-diagnostics
  needs. Different layer, different shape.
- `STEP_R12_emit_roh_intervals.R` — unlocks ROH pillar (v3.92).
- `STEP_R13_emit_sample_froh.R` — unlocks FROH pillar (v3.92).

## Outstanding scrubber-side (from v3.98 handoff)

- Real-data validation on the 226-sample LG28 cohort
- Per-candidate filtered TSV export buttons on page 10 (last v3.90 deferral)
- Annotation-track tooltips (heatmap left + top tracks)
- Recompute button for `purity_threshold` on candidate page (v3.81)
- L3 sub-band column captions
- K=6 coalescing UI (manual `g0a + g0b` merge)
- Per-L2 band diagnostics (currently only at ref_l2 anchor)
- Per-L2 het-shape (opt-in)
- Window-radius UI for live mode (5/10/20 toggle)

## File staging at /mnt/user-data/outputs/

- **STEP_M06_emit_boundary_evidence.R** — new, 612 lines
- **smoke_test_M06.R** — new, 351 lines, exits non-zero on any failure
- **cross_check_M06.js** — new, 278 lines, validates outputs against
  scrubber's actual `_buildBoundaryTrackScores` via vm.createContext
- **HANDOFF_v3_99_step_m06.md** — this file (R-side unblock; scrubber
  unchanged)

Scrubber and schema files unchanged this chat:

- pca_scrubber_v3.html — v3.98 (carried forward)
- SCHEMA_V2.md — v2.8 (carried forward)

## Notes for future-me / next chat

- **The 226-sample cohort is pure *C. gariepinus* hatchery broodstock.**
  K clusters reflect broodline structure, not species admixture. This
  matters for FST interpretation: high FST between K=3 regimes inside
  a candidate is an inversion signal, NOT a species/cohort split.
- **Don't conflate cohorts.** F1 hybrid (assembly paper),
  pure-gariepinus 226 (current boundary work), pure-macrocephalus wild
  (future paper).
- **The R `assembly` env at `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`
  has `data.table 1.14.10`, `jsonlite 1.8.8`** — verified during prior
  pipeline work. No new package installs needed.
- **Dosage TSV format flexibility.** STEP_M06 accepts both
  `pos_bp`-column form and marker-ID `<chrom>_<pos>`/`<chrom>:<pos>`
  form. If the actual cohort dosage uses some other format (e.g.
  bgzip-tabix .tsv.gz with a different header), this script may need
  a small extension — check the input format first.
- **Fish-regime calls TSV.** STEP_M06 expects the v3.80
  `fish_regime_calls.tsv` format (columns: `candidate_id`, `sample`,
  `regime` with values `g0`/`g1`/`g2`). The scrubber emits this TSV
  directly via the candidate-registry export button.
- **Cross-VM boundary type identity.** The cross-check uses duck-typing
  (`.BYTES_PER_ELEMENT === 8`) instead of `instanceof Float64Array`
  because functions extracted via `vm.createContext` get different
  constructor identities. Future cross-checks following this pattern
  should use the same duck-typing helper.
- **Empty boundary_evidence is graceful.** A candidate registry with
  no candidates on the chrom emits a JSON with `boundary_evidence: []`.
  The scrubber detects this as the layer being present-but-empty —
  same as no candidate matched.
