# Session handoff — 2026-05-01

**Status at end of session:** 16 atlas + 1 python suite all green (17 total).
All atlases shipped to `/mnt/user-data/outputs/`.

This handoff supersedes `HANDOFF_2026-04-30_to_2026-05-01.md` (which described
state up through turn 86). This session went from turn 87 → 105 (19 turns).

---

## TL;DR for tomorrow morning

1. **Detailed-mode + diamond architecture is shipped.** The 8-turn parallel
   candidate-system build (turns 87-95) is in the atlas. Sub-band visualization
   + heterozygosity histogram + 3-mode SNP density + manuscript vocabulary +
   structural-haplotype transition graph + band-reach + regime-breadth strip
   all live now (turns 96-105).

2. **Original tomorrow-morning plan still stands:** local PCA θπ on Inversion
   Atlas page 2 (`data-page="page12"`). Quentin brings R-side STEP_R39/R40/R41
   scripts + sample outputs + display preferences. Wire 4 JSON layers
   (theta_pi_per_window, theta_pi_local_pca, theta_pi_envelopes, plus
   per-sample heterozygosity for diagnostic mode).

3. **Original tomorrow-afternoon plan still stands:** Population Atlas
   relatedness page 3 (KING heatmap, family clusters from NAToRA pruning,
   1st/2nd/3rd-degree thresholds ≥0.177/≥0.0884/≥0.0442). Genome page
   (separate `Genome_atlas.html`).

4. **What's NEW for Quentin to look at first thing:** the toolbar on the
   per-sample lines panel got 3 new toggles — see "New UI" section below.
   The catalogue got a Diamond column. The karyotype page in detailed mode
   got a classification panel + heterozygosity histogram. Help page got 4
   new sections.

---

## Turns shipped this session (87 → 105)

### Architecture core (87 → 95)
- **87** — Karyotype label vocabulary resolver. Sidebar `labels: legacy / detailed`
  toggle. K=3 detailed labels: H1/H1, H1/H2, H2/H2. K=4..6 extends to H3 system.
  9 K=3 buttons across 3 sites tagged `data-vocab-band-k3` for live update.
  27 checks.
- **88** — Parallel candidate registry (data layer). `state.candidate_detailed`,
  `state.candidates_detailed`, `state.candidateList_detailed`. `state.activeMode`
  = 'default' | 'detailed'. `getActiveCandidate()` / `setActiveCandidate()`
  routing. `_pcrDeepCloneCandidate()` cross-realm safe (uses constructor.name
  for typed arrays). `initDetailedFromDefault()` 1-to-1 duplication, idempotent.
- **89** — Merge isolation. `assertSameMode(a, b)`, `buildContingencyForCandidates`,
  `mergeIsolationAudit()` for runtime safety. Existing default-mode merge code
  is auto-isolated because it only writes to default-mode slots.
- **90** — Detailed-mode classifier. `classifyDetailedCandidate(cand, het?)`
  returns per_sample + per_band + karyotype_assignment with H-system labels.
  Bimodal het detection via K-means K=2 with absolute-gap minimum (0.15) +
  z-score threshold (2.0). Source = 'precomp_het' or 'atlas_proxy' (PC residual
  fallback via `_diagComputeCandidateSuspicion` from turn 83).
- **91** — UI for activeMode toggle (separate from labelVocab). Sidebar
  `system: default / detailed` bar + cyan banner under it warning about
  operational caveat in detailed mode. First switch to detailed seeds the
  registry from default via `initDetailedFromDefault()`.
- **92** — Diamond detector. `detectDiamonds(candidate)` returns array of
  `{splitting_band, diamond_start_w, diamond_end_w, stable_bands, slanting_bands,
  strict, strict2, baseline_spread, peak_spread, peak_spread_ratio}`.
  `summarizeDiamonds()` for catalogue. Strict = ≥1 stable parallel band.
  Strict-2 = ≥2.
- **93** — Catalogue Diamond column. 3-button strictness toggle in toolbar
  (`◆ loose / ◆ strict / ◆ strict-2`). Cell rendering with diamond glyph + count.
  `_catBuildTransientCandidateFromL2()` + `_catComputeRowDiamondCount()`.
- **94** — Cyan-alpha diamond zone overlay on PC1 sub-panel. Translucent cyan
  rectangle behind diamond windows + dashed edges + "◆ split detected"
  annotation strip. Visibility heuristics: only when diamond width > 12 pixels
  AND covers >= 8% of visible mb range.
- **95** — SNP-density visualization (3-mode toggle: off / strip / shade).
  `_drawSnpDensityStrip()` paints cool→warm gradient bar above plot.
  `_snpDensityForWindow()` priority chain: precomp:n_snps → precomp:snp_count
  → precomp:lambda1 → precomp:variance_pc1 → atlas-side proxy from PC1
  variance.

### Polish (96 → 99)
- **96** — Precomp spec doc updated. Section 8 (`per_sample_het` per candidate)
  + section 9 (`n_snps` per window). Migration order updated to 9 items. Atlas
  classifier auto-reads `cand.per_sample_het` if present. Spec doc 471 → 622
  lines.
- **97** — Sub-band visualization. `_renderDetailedClassificationPanel()`
  renders per-band cards on karyotype page in detailed mode. Each card shows
  K-band swatch, label (legacy or H-system per current vocab), interpretation
  tag (hom-like / het-like / mixed / ambiguous), heterozygosity stacked bar
  (low/mid/high counts), sub-band swatches with mode centers when band is mixed.
- **98** — Heterozygosity histogram. `_renderDetailedHetHistogram()` paints an
  SVG histogram of all per-sample het values with bars colored by K-means band,
  stacked when bands overlap. 16 bins, 5 axis ticks, max-count y-label, frame.
  Sits at top of detailed panel.
- **99** — SNP density `shade` mode added (3rd toggle option after off / strip).
  `_drawSnpDensityShade()` paints translucent vertical bars across full plot
  height for low-density windows. Single-checkbox UI replaced by 3-button bar.
  Backwards-compat with legacy localStorage key.

### Manuscript + analysis primitives (100 → 105)
- **100** — Manuscript vocabulary section. Tiered terms in help page:
  structural-haplotype system → segregating structural-haplotype system →
  inversion-like structural-haplotype system → candidate polymorphic inversion
  → polymorphic inversion. Includes paste-ready methods statement and L1/L2/L3
  hierarchy mapped to new vocabulary.
- **101** — Structural-haplotype transition graph (data layer).
  `computeStructuralHaplotypeTransitionGraph()` walks adjacent L2 pairs.
  Per-boundary stats: `transition_rate`, `n_changed`, `n_samples`, `edges`
  (sorted by count desc), `position_mb`. Hotspot threshold 0.30. 12 checks.
- **102** — Transition-rate strip on lines panel. `_drawTransitionRateStrip()`
  paints red-graded bars per L2 boundary at bottom of PC1 sub-panel. Hotspots
  get a thin vertical tick across plot height. Toggle next to SNP density bar.
- **103** — Band-reach + regime-breadth data layer.
  `computeBandReachAcrossL2s(l2_indices)` projects K-means labels through
  Hungarian alignment chains, returns per-fish band-reach + per-band visit
  counts + categorical regime breadth ('narrow' / 'medium' / 'wide' /
  'no_signal'). `computeWindowedBandReachPerL2()` returns 3-L2 sliding-window
  summary per envelope. 14 checks.
- **104** — Regime-breadth strip on lines panel. `_drawRegimeBreadthStrip()`
  paints colored per-L2 segments above PC1 panel: green = narrow (clean
  segregation), amber = medium, red = wide (no clean segregation), grey =
  no_signal. Toggle next to trans-rate.
- **105** — Help-page sections added: "Band-reach + regime breadth" explaining
  the metric + "Hatchery applications" with 5 use cases + paste-ready
  manuscript paragraph for discussion section + paste-ready paragraph for
  diversity preservation rationale.

### Window-mode candidate drafts + family color (106 → 108)
- **106** — Window-mode candidate drafts. When sidebar stepMode is 1w / 5w /
  10w / Nw, ↑/↓ arrow keys now extend / shrink the L3 draft by windows
  instead of by L2 envelopes. Promoting commits a candidate with raw
  start_w / end_w that may fall **inside** an existing L2 envelope —
  effectively creating candidate-level boundaries at any window. The L2
  envelope registry is NOT modified (Option B from the design discussion).
  Schema fields `resolution: 'L2' | 'W'`, `start_w`, `end_w` were already in
  place from turn 10; turn 106 connects the arrow-key wiring to drive them.
  `_segmentDraftByCuts` and `_buildCandidateFromSegment` already handled the
  W path; they just needed the upstream draft state to reach them. Banner
  in focal pane shows `w 312…320 (9 windows)` instead of `L2 7…7 (1 L2)` in
  window mode. 14 test checks.
- **107** — Visual feedback on L2 strip for window-mode candidates. Bright
  cyan tick markers drawn at every committed candidate's start/end window
  positions on both the main Z-panel L2 row AND the minimap L2 row, but
  ONLY for candidates with `source === 'l3_draft_w'` (turn-106 tagged).
  Active draft in window mode shows GOLD tick markers (matches draft
  palette) that update live as the user extends/shrinks. The L2 segment
  blocks themselves stay unchanged because the L2 registry isn't
  authoritative anymore — Quentin's design choice: "L2 are useful for
  navigation but if we operate on window and promote to candidates then
  L2 are just for navigation." Communicates "user-defined unit boundaries
  here" without falsifying L2 envelopes. 1 additional test check (e2e
  commit verifies source='l3_draft_w' tagging).
- **108** — `color: family` mode actually colors the per-sample lines.
  Pre-existing bug: the dropdown UI + state slot + resolver hooks were all
  in place since v3.99 t14e but the lines drawing code never called the
  resolvers — every sample was painted with the same fixed grey
  `rgba(180,190,210,0.10)` regardless of mode. So switching to "family"
  did nothing visible. Fix: implemented the family branch in both
  `_resolveSampleColorByMode` and `_resolveSampleScopeColor` (uses
  existing `familyColor(si)` which reads `state.familyPalette` + sample
  `family_id`, with FAMILY_COLOR_SMALL / SINGLETON / UNMATCHED fallbacks).
  Modified the offscreen-canvas line drawing to call
  `_resolveSampleScopeColor(si, lcMode)` per sample when `lcMode === 'family'`.
  baseAlpha bumped from 0.10 to 0.25 so saturated family colors remain
  readable when 226 lines overlap. Cache key now includes `|lc=` mode for
  invalidation. The other modes (dosage, ghsl, het, theta_pi, froh,
  confounder_alert) remain stubs awaiting their respective JSON layers.

---

## New UI for Quentin to find tomorrow

### Sidebar (left of every page)
```
labels:    [legacy*]  [detailed]      ← turn 87, vocabulary toggle
system:    [default*] [detailed]      ← turn 91, candidate-system toggle
   (banner appears in detailed mode)
```
The two toggles are **independent**. 4 valid combinations: legacy+default
(default behavior), legacy+detailed (new candidates use heterozygosity but
labels stay neutral), detailed+default (default candidates with H-system
labels), detailed+detailed (full detailed mode).

### Per-sample lines panel toolbar (between K mode and K=3/K=6 buttons)
```
... lasso  SNP dens: [off*] [strip] [shade]  trans-rate ☐  regime ☐ ...
                     turn 95+99           turn 102      turn 104
```
- **strip** = thin gradient bar above plot (cool=low SNP density, warm=high)
- **shade** = translucent vertical bars across plot for low-density windows
- **trans-rate** = red bars at bottom of plot, one per L2 boundary; hotspots
  get a vertical line across plot
- **regime** = colored strip above plot, per L2: green/amber/red/grey

### Catalogue toolbar (page 5)
```
... [simple] [detailed*]  ◆ loose  ◆ strict  ◆ strict-2  ...
                          turn 93 diamond strictness
```
New `Diamond` column shows count of detected diamonds at the active strictness.

### Karyotype page (page 7) in detailed mode
A new cyan-bordered panel appears between the header and the toolbar:
- 16-bin SVG histogram of all per-sample heterozygosity values, stacked by
  K-means band
- Per-band cards: K-band swatch + label + H-class + interpretation tag
  (hom-like / het-like / mixed / ambiguous) + heterozygosity stacked bar +
  sub-band mode centers if mixed

In **default mode**, the karyotype page is unchanged.

### Help page (page 16, "about")
Four new sections in order:
1. "Manuscript vocabulary" (turn 100) — tiered terminology
2. "L1/L2/L3 hierarchy in the new vocabulary" (turn 100)
3. "Band-reach + regime breadth" (turn 103-104)
4. "Hatchery applications" (turn 105)

All four have paste-ready paragraphs in styled boxes for direct copy into the
manuscript.

---

## Test suites (17 total)

```
inversion atlas page 16 (cross-species):       35 ✓
inversion atlas page 17 (stats profile):       24 ✓
inversion atlas page 18 (marker panel):        30 ✓
stage bands + pulses:                          21 ✓
negative regions + caution:                    27 ✓
diagnostic mode:                               14 ✓ (occasionally flaky on
                                                    Math.random outlier color)
cross-page coloring:                           18 ✓
Q-ancestry:                                    26 ✓
label vocab (turn 87):                         27 ✓
multi-turn (turns 88-90, 92):                  33 ✓
ui multi-turn (turns 91, 93-95, 102):          ~30 ✓
transition graph (turn 101):                   12 ✓
band reach (turn 103):                         14 ✓
population breeding:                           28 ✓
hash routing:                                  10 ✓
hatchery health:                               28 ✓
python STEP_CS01:                               ✓
```

**Total: ~395 atlas + python checks green.**

---

## Tomorrow morning — local PCA θπ on page 2

(unchanged from previous handoff, stays the plan)

Required from Quentin:
- R-side STEP_R39/R40/R41 scripts
- Sample JSON outputs from those scripts
- Display preferences for the page (heatmap? lines? both?)

Atlas-side work to do:
- Wire 4 JSON layers per chromosome:
  - `theta_pi_per_window` — per-window θπ values (numeric grid)
  - `theta_pi_local_pca` — local PCA on θπ matrix (per-window PC1, PC2)
  - `theta_pi_envelopes` — L1/L2 envelopes computed on θπ-PCA
  - per-sample heterozygosity (becomes the input to detailed-mode classifier!)
- Add a 4-panel page on `data-page="page12"`
- Plug per-sample heterozygosity into `cand.per_sample_het` so the detailed-mode
  classifier (turn 90) goes from `atlas_proxy` to `precomp_het` source

---

## Tomorrow afternoon — Population Atlas relatedness page

KING heatmap on page 3. Family clusters from NAToRA pruning. Thresholds:
1st degree ≥ 0.177, 2nd ≥ 0.0884, 3rd ≥ 0.0442.

Then `Genome_atlas.html` as separate page.

---

## Files at /mnt/user-data/outputs/

- `Inversion_atlas.html` — 2.3 MB, 49,133 lines
- `Population_atlas.html` — 147 KB, unchanged this session
- `cross_species_breakpoints.tar.gz` — ~95 KB
  - `precomp_diagnostics_spec.md` — 622 lines (sections 1-9)
  - 18 test scripts (16 atlas + python + STEP_CS01)
  - 5 demo JSONs
- `population_atlas_breeding.tar.gz` — ~16 KB

---

## Open / deferred (not blocking tomorrow's plan)

- **🟡 RESEARCH DIRECTION (turn 108): genome-wide crossover + GC detector.**
  Quentin proposed combining contingency tables + dosage heatmaps + per-sample
  lines as a unified pattern recognition primitive. Quote:
  > "If we can succeed to have a highly reliable contingency method crossed
  > with dosage heatmaps and per sample lines we can solve all of crossovers
  > and GCs in one go for whole genome. it's pattern recognition × biology
  > (GC vs crossovers) × contingency patterns × sample lane movements ×
  > structured haplotype intervals × window SNP density as confounder ×
  > window coverage in the genome as confounder bc if no window it makes
  > linear bands. if high or low snps it distorts."
  Concrete observation that motivated this: at a single window (1w slab), the
  -1 / focal / +1 contingency tables can show "18 grey samples become blue;
  then 14 blue become grey in focal" — patterns consistent with crossovers
  or gene conversions, distinguishable by their length signature (GC = short,
  one or few windows; crossover = long, persists). Worth scoping as a major
  analytical primitive in a future session: needs a confounder-aware
  detector over (lane-change events × dosage residuals × SNP-density gating
  × coverage gating). Existing turn-101 transition_rate primitive is the
  starting infrastructure but doesn't yet resolve GC-vs-crossover or apply
  the SNP-density / coverage confounder gates.

- **🔴 BUG: slab focal panes (1w/5w/10w/Nw) use bare-bones rendering vs L2's high-quality renderer.**
  Identified end of session via 3 screenshots from Quentin (LG28).
  - L2 mode renders through `focalContentHtml()` + `drawMiniPCA()` — convex-hull
    cluster polygons, λ₁ axis labels, K-mode chips (3 / 6 / 3+6 / s / u / d),
    chip-rich invariant + K-specific metadata, conservation row, recluster-mode dropdown.
  - 1w/5w/10w/Nw modes render through `slabFocalContentHtml()` + `drawSlabMiniPCA()`
    — bare PCA dots (no polygons, no λ₁ labels), bare metric rows ("windows | K used | fam_purity"),
    no K-mode chips, no recluster dropdown.
  - Affects: 1w / 5w / 10w / N w buttons in the L3 contingency tables toolbar.
  - Layouts +Left, +Right, +L/R, ±2, Dual all inherit the bare slab renderer when
    triggered from non-L2 step modes.
  - Root cause: `slabFocalContentHtml` (line ~29167) is a parallel-but-poorer
    implementation of `focalContentHtml` (line ~29846); `drawSlabMiniPCA`
    (line ~29227) is a parallel-but-poorer implementation of `drawMiniPCA` (line ~30655).
  - **Fix scope:** at minimum 3 code paths — refactor `slabFocalContentHtml` to use
    the same `_invariantMetaInlineHtml` + `_kSpecificMetaInlineHtml` blocks; add
    convex-hull + axis-label rendering to `drawSlabMiniPCA` (or generalize
    `drawMiniPCA` to accept a window range instead of an L2 envelope); wire the
    K=3/K=6/s/u/d chip set + recluster-mode dropdown into the slab header.
  - **Why deferred:** subtle interactions with contingency-table logic,
    chip-restrict flow, and per-pane Hungarian alignment make this a "small bugs
    visible later" risk if rushed. Better as a deliberate single-purpose session.
  - **Test plan when fixing:** load a chromosome, switch to L2 mode, screenshot.
    Switch to 1w mode at the same window, screenshot. Switch to 10w mode +L/R
    layout, screenshot. All three should show: convex hulls, λ₁ axis labels,
    K-mode chips, conservation row. Cross-check against contingency tables in
    +L/R layout — they currently work in slab mode (Image 3 shows them) so the
    chip-restrict flow is fine, just the focal pane rendering is poor.

- **Per-line gradient stroke for SNP density** — `shade` mode covers the
  diagnostic need without the hot-path refactor; per-line stroke remains
  deferred until needed
- **L3 contingency validation on LG28** — Quentin HPC work, not atlas
- **Manuscript wording pass beyond help-page tooltips** — Quentin work
- **Atlas other internal `inversion-system` strings** — turn 100 added the new
  vocabulary to the help page but didn't rename other internal mentions
  (regime registry tooltips, σ verdict descriptions). Should be a single pass
  before manuscript submission, low risk
- **Real candidate JSON with diamond pattern** — atlas detector + overlay
  should be tested against the actual LG28 candidate Quentin photographed.
  Will probably reveal threshold tuning needs (e.g., `_DD_MIN_DIAMOND_WINDOWS`
  may need to drop from 3 to 2 for short candidates)

---

## Carried-forward state (unchanged from previous handoff)

- Three-cohort discipline: F₁ hybrid (assembly only) / 226 pure C. gariepinus
  hatchery (current MS) / pure C. macrocephalus wild (future). 226-sample
  cohort never F₁ hybrids; K clusters = hatchery broodline structure.
- F_ROH | H is the manuscript's novel framework (3 metrics: F_ROH | H,
  Burden | H, π_D / π_S).
- "Continue" means proceed naturally — don't ask, just go.
- Each turn ships independently with full test sweep green.
- Verdict text colored, structural data normal ink (`.fhh-good / .fhh-review
  / .fhh-caution`).
- JSDOM tests at `/tmp` with `NODE_PATH=/tmp/node_modules`. Helpers exposed
  via `window.X` near top of script (NOT end — JSDOM aborts on canvas init).
- Cross-realm safety: typed-array `instanceof` checks fail across JSDOM
  realms; use `constructor.name` instead.

