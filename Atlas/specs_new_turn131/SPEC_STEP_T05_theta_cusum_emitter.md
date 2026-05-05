# SPEC — `STEP_T05_theta_cusum.R` — per-carrier θπ CUSUM emitter

**Status**: drafted turn 132. Not yet implemented on LANTA. The atlas
hero (`#thCusumHeroPanel` on page 12, behind the `2 local PCA θπ` tab)
already paints this layer's output — see `_drawThCusumHero()` in
`Inversion_atlas.html` (turn 132 Slice 3). What's missing is the
R-side script that produces the JSON the renderer reads.

**Working directory** (LANTA):
`/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/phase_4_resolution/4b_theta_resolution/STEP_T05_theta_cusum.R`

**Builds on** (LANTA):
- `phase_4_resolution/shared_lib/lib_persample_cusum.R` — pure CUSUM
  utility (chat `487c7f04` turns 120-122). The lib was deliberately
  walked back to emit only `persample` (no aggregated modes / no
  Hartigan dip / no fitted spread classification) per Quentin's
  pushback: *"How can choose the statistical test of a distribution
  shape that we don't know what it will look like. Observe first
  empirically then correct the kernel."*

**Trigger** (Quentin, chat `487c7f04`, when shown the boundary-distribution
discovery):
> *"THATS EXACTLY RIGHT. AND SOMETIMES IN HATCHERY WITH FAMILY LD WE
> HAVE ERODED BREAKPOINTS AND NOT ERODED BREAKPOINTS. (THIS ONE I
> INVENTED BUT I THINK ITS TRUE LIKE)."*

That hatchery-LD insight is what makes the per-carrier distribution
worth visualizing rather than a single consensus boundary call:
narrow distributions = founder-block-locked breakpoint; ragged
distributions = old breakpoint with gene-conversion erosion across
lineages.

---

## 1. Why this layer exists

Page-1 dosage local-PCA has sharp breakpoints because dosage is a
binary in/out signal. **θπ is a diversity gradient, not a step** — so
the standard sim_mat → MDS → eigenvalue-Z pipeline gives diffuse
diversity boundaries that look beautiful in genome view but are useless
for declaring "the inversion runs from `start` to `end`."

The fix is **per-carrier CUSUM** on each sample's per-window θπ track.
The maximum of `|cumsum(x − mean(x))|` along the interval gives that
sample's individual best-supported changepoint. Across the cohort, the
**distribution of those 60-or-so changepoint positions** is the
informative quantity:

- **Tight cluster** (IQR < 100 kb) → one shared founder breakpoint,
  inherited as a block. Manuscript phrasing: *"the 5' breakpoint
  resolves to 14.01 Mb (IQR 70 kb across 60 carriers), consistent with
  recent founder-locked transmission."*
- **Ragged distribution** (IQR > 500 kb) → erosion. Manuscript phrasing:
  *"the 3' breakpoint shows ragged carrier-level resolution
  (IQR 1.2 Mb across 60 carriers), consistent with gene-conversion
  erosion."*
- **Multi-modal** → nested inversions or distinct ancestral haplotypes.

This is the **only signal the dosage scrubber can't produce**, which
is why the θπ page hero is CUSUM-based rather than a sim_mat copy.

---

## 2. Inputs

### 2.1 Per-sample θπ matrix

`theta_pi_per_window` JSON layer (produced by `STEP_R39`, also pending).
Shape per the atlas detector at `Inversion_atlas.html` ~line 32521:

```json
{
  "schema_version": 1,
  "layer": "theta_pi_per_window",
  "chrom": "LG28",
  "scale": "win50000.step10000",
  "n_windows": 3300,
  "n_samples": 226,
  "window_pos_bp": [start_bp_w0, start_bp_w1, ..., start_bp_wN-1],
  "samples": [
    { "sample_id": "Cga_001", "theta_pi": [<float>, ...], "n_sites": [<int>, ...] },
    ...
  ],
  "values": [<flat n_windows × n_samples row-major>]   // optional convenience
}
```

The `samples[i].theta_pi` array is the per-window per-site θπ
(scale: 50 kb window, 10 kb step — see SCHEMA §22 / atlas page-12 docs).

### 2.2 Carrier list per candidate

For each candidate interval, T05 needs the **carrier subset** (samples
whose dosage karyotype is HET or HOM_INV — non-REF). This comes from
the existing per-candidate karyotype TSV emitted by phase 2 / `STEP21`.
Path on LANTA:

```
inversion_codebase_v8.5/phase_2_discovery/.../candidate_<cand_id>_karyotype.tsv
```

Columns expected: `sample_id, karyotype_class` where
`karyotype_class ∈ {HOM_REF, HET, HOM_INV}`.

### 2.3 Candidate definitions

The candidate atlas — same TSV the dosage scrubber consumes
(`atlas_catalogue.tsv` or equivalent), with at least:
`cand_id, chrom, start_bp, end_bp`.

T05 runs once per candidate.

---

## 3. Algorithm — wraps `lib_persample_cusum.R`

For each candidate `(chrom, start_bp, end_bp)`:

### 3.1 Build the carrier × window submatrix

1. From `theta_pi_per_window`, restrict columns (windows) to those
   whose `window_pos_bp` falls in `[start_bp − flank, end_bp + flank]`
   where `flank = 0.25 × (end_bp − start_bp)` (so CUSUM has room to
   resolve the boundary against flanking baseline).
2. Restrict rows (samples) to carriers (`karyotype != HOM_REF`).
3. Result: `M[carriers × windows_in_region]`.

### 3.2 Call the shared lib

```r
source("../shared_lib/lib_persample_cusum.R")

result <- compute_persample_cusum(
  matrix       = M,
  window_pos_bp = window_pos_bp_subset,   # length = ncol(M)
  sample_ids   = rownames(M)
)
# result$persample is a data.table with one row per carrier:
#   sample_id, cp_bp, strength, asymmetry
```

`compute_persample_cusum` is the deliberately-minimal walked-back lib.
It does **NOT** compute KDE modes, Hartigan dip, spread class, or any
distributional aggregate. Per the chat-`487c7f04` walked-back design.

### 3.3 Annotate carrier rows with karyotype

Join `result$persample` against the candidate karyotype TSV to attach
each carrier's HET / HOM_INV label. The atlas renderer reads
`persample[i].karyotype` to color tick marks.

### 3.4 No further analysis

T05 explicitly does not compute:
- Mode positions
- IQR-based spread classification (atlas does this descriptively for display)
- Dip-test p-values
- Cross-stream concordance with GHSL (that's `STEP_DC06`'s job)

The per-carrier observations are the deliverable.

---

## 4. Output — `cusum_theta` JSON layer

Per-chromosome JSON (one file per chrom containing all candidates on
that chrom — matches the atlas's chrom-keyed loading model). Path:

```
phase_4_resolution/4b_theta_resolution/output/<CHROM>_cusum_theta.json
```

### 4.1 Top-level shape

```json
{
  "schema_version": 1,
  "layer": "cusum_theta",
  "chrom": "LG28",
  "source_script": "STEP_T05_theta_cusum.R",
  "input_layer": "theta_pi_per_window",
  "scale": "win50000.step10000",
  "n_candidates": 7,
  "candidates": [ <one entry per candidate, see 4.2> ]
}
```

### 4.2 Per-candidate entry

```json
{
  "cand_id": "LG28_cand_03",
  "chrom": "LG28",
  "range_bp": { "start": 14000000, "end": 16500000 },
  "n_carriers": 60,
  "n_carriers_skipped": 2,
  "skipped_reasons": ["NA_too_dense", "all_NA"],
  "persample": [
    {
      "sample_id": "Cga_017",
      "cp_bp": 14012500,
      "strength": 4.21,
      "asymmetry": -1,
      "karyotype": "HET"
    },
    ...
  ]
}
```

### 4.3 Atlas-side mapping

The atlas's existing single-`cusum_theta` model assumes one candidate's
worth of data at a time. The atlas reads
`state.data.cusum_theta.persample[]` directly — so the chrom-level JSON
above must be **flattened by the atlas loader** to the active-candidate
view: when the user selects a candidate, the loader sets

```js
state.data.cusum_theta = {
  schema_version: 1,
  layer: "cusum_theta",
  chrom: <chrom>,
  source_script: "STEP_T05_theta_cusum.R",
  range_bp: <selected candidate's range_bp>,
  persample: <selected candidate's persample[]>,
};
```

This atlas-side flatten step is **not in this spec** — it's a small
follow-up loader patch (estimated ~0.3 turns).

### 4.4 Why per-chrom file (not per-candidate file)

- Mirrors how `theta_pi_local_pca` and `theta_pi_envelopes` are
  organized (chrom-keyed)
- Avoids file-fragmentation: 226-sample × ~30 candidates per chrom × 28
  chroms = up to 840 files if per-candidate
- The atlas chrom-load orchestrator (planned in
  `specs_new_turn131/SPEC_multichrom_load_orchestrator.md`) loads one
  chrom at a time anyway

---

## 5. Edge cases the script must handle

### 5.1 Carriers with too few non-NA windows

If a carrier has fewer than 10 non-NA windows in the candidate region,
its CUSUM is unreliable. Skip the carrier and increment
`n_carriers_skipped`. Record reason in `skipped_reasons` array
(at the candidate level, not per-sample, to keep JSON small).

Threshold rationale: at 50 kb windows / 10 kb step, 10 non-NA windows
covers ~140 kb of effective callable density — minimum for the
cumsum-of-deviations signal to be meaningful.

### 5.2 Constant signal (sd(x) = 0)

`strength = max|cumsum| / sd(x)` blows up if a row's per-window θπ is
constant. Skip with reason `"constant_signal"`. Should be very rare
since θπ varies by site density even in flat regions.

### 5.3 No carriers

If `n_carriers == 0` (e.g. a candidate that's actually a sweep with
no real carriers), emit the candidate entry with empty `persample[]`.
Atlas renders as an empty hero panel (renderer already handles this —
see Slice 3 test 2b).

### 5.4 Candidate exceeds available windows

If the candidate's `[start_bp, end_bp]` extends past the chrom's
window grid (very near telomere), restrict to available windows and
emit a warning to stderr but keep the candidate.

---

## 6. Reproducibility / pinning

- Random seed: not used (CUSUM is deterministic given input matrix)
- R version: pin to 4.4.x (matches the rest of phase_4_resolution)
- Required packages: `data.table`, `jsonlite`, `optparse`. No KDE
  packages, no `diptest` — those were removed in the lib walk-back.

---

## 7. Test fixture

A small synthetic input lets the script self-test:

- 4 carriers on a 30-window region
- Carrier A: clean step at window 12 → cp_bp expected at corresponding bp
- Carrier B: same step at window 12 (tight cluster)
- Carrier C: step at window 18 (different breakpoint — second mode)
- Carrier D: no step (flat) → low strength, cp arbitrary

Expected output:
- `n_carriers = 4`
- A and B have cp_bp within one window of each other
- C cp_bp differs from A/B by ~6 windows worth of bp
- D's `strength` is the smallest

This fixture should live in
`phase_4_resolution/4b_theta_resolution/tests/synthetic_cusum_theta.R`
and be runnable via `Rscript`.

---

## 8. Companion specs (to write later)

These are the layers `cusum_theta` interlocks with. Not in scope for
this spec.

| Spec | Status | Purpose |
|------|--------|---------|
| `SPEC_STEP_R39_theta_pi_per_window_emitter.md` | TODO | The dense θπ matrix this script consumes |
| `SPEC_STEP_R02_ghsl_cusum.md` | TODO | GHSL parallel — produces `cusum_ghsl` |
| `SPEC_STEP_DC06_cusum_concordance.md` | TODO | Cross-references R02 + T05; produces `cusum_concordance` for the concord badge |
| `SPEC_atlas_loader_chrom_to_candidate_flatten.md` | TODO | The ~0.3-turn atlas patch that flattens per-chrom JSON to per-candidate `state.data.cusum_theta` |

---

## 9. Atlas integration — already in place

The atlas hero is already wired and tested. As a sanity reference for
whoever builds T05:

- **DOM**: `#thCusumHeroPanel` inside `<div id="page12">` in
  `Inversion_atlas.html`
- **Visibility**: `_refreshThetaPiPanelVisibility()` reveals the panel
  when `state.layersPresent.has('cusum_theta')` is true
- **Renderer**: `_drawThCusumHero()` reads
  `state.data.cusum_theta.persample[]` and paints
  `#thCusumStripCanvas` (per-carrier rug) +
  `#thCusumHistCanvas` (cohort histogram with median + IQR shading)
- **Spread badge**: `#thCusumSpreadBadge` — IQR-based, descriptive,
  thresholds tight<100 kb / intermediate<500 kb / ragged>500 kb
- **Concord badge**: `#thCusumConcordBadge` — reads
  `state.data.cusum_concordance` if present (from STEP_DC06); shows
  `concord —` until that layer ships

Tests covering this end of the contract:
`tests/test_turn132_page12_cusum_hero_renderer.js` (39 / 0).

---

## 10. Build estimate

~1 LANTA session. The shared lib already exists per chat `487c7f04`
turn 120-122; T05 is a thin wrapper that loops over candidates,
restricts the matrix, calls the lib, and writes JSON. The genuine
effort is in plumbing the input layer (`theta_pi_per_window`), which
requires `STEP_R39` to ship first.

**Build order**: `STEP_R39` → `STEP_T05` → atlas chrom-to-candidate
flatten patch → real data flows end-to-end into the page-12 hero.
