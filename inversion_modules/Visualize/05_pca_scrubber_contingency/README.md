# 05_pca_scrubber_contingency

Live PCA scrubber for L3 contingency-table debugging on the multipass
catalogue. Lets you scrub through windows on a chromosome, see L1/L2
envelope structure overlaid on sim_mat and Z, and watch K=2 / K=3
K-means contingency tables evolve **live** between the focal L2 envelope
and configurable neighbors within the same L1 parent.

Output of this folder is a **calibrated `merge_threshold` value** for
the batch L3 pipeline (`STEP_D17_multipass_L3_v1.R`, to be written next).


## Folder layout

```
05_pca_scrubber_contingency/
├─ 01_input/    L1/L2 envelopes + boundaries TSVs (LG28 sample)
├─ 02_code/     export_precomp_to_json_v3.R + pca_scrubber_v3.html
├─ 03_output/   put generated LG28.json here
└─ docs/        notes
```


## Workflow

### 1. Build the JSON

The exporter accepts multi-scale sim_mat + sample-identity layer
(bamlist + relatedness) so you can switch smoothing and color modes
in the browser without re-exporting.

```bash
Rscript 05_pca_scrubber_contingency/02_code/export_precomp_to_json_v3.R \
  --precomp        01_input/C_gar_LG28.precomp.slim.rds \
  --sim_mat_nn40   01_input/C_gar_LG28.sim_mat_nn40.rds \
  --sim_mat_nn80   01_input/C_gar_LG28.sim_mat_nn80.rds \
  --sim_mat_nn160  01_input/C_gar_LG28.sim_mat_nn160.rds \
  --sim_mat_nn320  01_input/C_gar_LG28.sim_mat_nn320.rds \
  --bamlist        ../05_ngsrelate/list_of_samples_one_per_line_same_bamfile_list.tsv \
  --pairs          ../05_ngsrelate/catfish_226_for_natora.txt \
  --theta_cutoff   0.177 \
  --samples        <optional path>/samples.tsv \
  --l1_envelopes   01_input/C_gar_LG28_d17L1_envelopes.tsv \
  --l2_envelopes   01_input/C_gar_LG28_d17L2_envelopes.tsv \
  --l2_boundaries  01_input/C_gar_LG28_d17L2_boundaries.tsv \
  --l1_boundaries  01_input/C_gar_LG28_d17L1_boundaries.tsv \
  --out            03_output/LG28.json
```

**Sample identity layer** (ported verbatim from `STEP_D16` v3.2):

| flag | role |
|---|---|
| `--bamlist` | one CGA per line, in precomp Ind* order. Maps `Ind0..IndN -> CGAxxx` for display labels. |
| `--pairs` + `--theta_cutoff` | union-find on the relatedness graph (theta ≥ cutoff). Each connected component = one family. `family_id 1` = biggest component. **RECOMMENDED** — matches the network figure exactly. Default cutoff: 0.177. |
| `--family` | fallback: 2-col TSV `cga<TAB>family_id` (natora-style). Only use if you don't have the raw pairs file. |
| `--samples` | optional 3-col TSV `ind<TAB>cga<TAB>ancestry` for ancestry coloring. |

Each scale's sim and per-distance local Z (computed exactly like
`STEP_D17c_overlay_plot_L2_v8.R::compute_local_z`) is thumbnailed at
200×200 with q05/q95 anchors. The browser renders the heatmap with
**exactly the PDF styling**.


### 2. Open the scrubber

```bash
xdg-open 05_pca_scrubber_contingency/02_code/pca_scrubber_v3.html
```

Pick `03_output/LG28.json` in the file input.


### 3. Sidebar controls

**Sim_mat smoothing** — dropdown over scales (nn40 / nn80 / nn160 / nn320).
Toggle "PDF-style" to swap between triangle-split rendering and the
legacy single-channel thumbnail.

**Where am I** — current L1 + L2 envelope info, with neighbor strip.

**L3 clustering** — K (2 or 3), aggregation method, merge τ, alpha,
min n/group, min n_windows.

**Display — color samples by**: 4-button row
- **cluster**: live K-means group, Hungarian-aligned to focal envelope
  (the colors that drive the contingency tables)
- **family**: hub families (n≥4) get distinct Tol-qualitative colors,
  small families (n=2-3) get grey, singletons (n=1) get lighter grey,
  unmatched (family_id=-1) get neutral grey. Auto-disabled when no
  family info was loaded.
- **ancestry**: by ancestry value from `--samples`. Auto-disabled when
  fewer than 2 distinct ancestry classes.
- **none**: uniform grey (read positions only)

**Tracked samples** — slider for N (0–30), then "Auto-pick: spread"
or "Auto-pick: radial". Click any dot in the main PCA to add/remove.
Chip list shows `current/N`.

**Jump** — unit dropdown (Mb / window / bp), value, Go (or Enter).
For "window", value is 1-based.

**L1 / L2 step buttons** — jump to next/prev envelope of either level.


### 4. L3 panel layout

Five buttons at the top of the L3 panel:

| button | shows |
|---|---|
| **Focal only** | Just the focal L2 — its K-means details + mini PCA |
| **+Left** | Left neighbor + Focal |
| **+Right** | Focal + Right neighbor |
| **+L/R** | Left, Focal, Right |
| **±2** | L−2, L−1, Focal, R+1, R+2 |

Each column has a **mini PCA** at the top showing samples colored by
the current color mode. In **cluster** mode, every neighbor's labels
are Hungarian-aligned against the **focal** envelope's labels (not
the immediate neighbor's), so colors are consistent across the strip.
Tracked samples carry their identity color + a halo across all minis.

The **focal** column shows:
- windows, n per group, center PC1
- power (OK / LOW_NWIN / LOW_GROUP_N)
- **fam_purity** (size-weighted mean across K-clusters of "fraction of
  samples in this cluster from its single most-common family")
- **FAM-LD?** badge when `fam_purity ≥ 0.70` — bands ARE families,
  ≠ inversion
- per-cluster top-family breakdown: `k0: F3 (12/15, 80%)` etc.

The **neighbor** columns show K×K contingency tables with concord %,
Fisher (K=2) or chi² (K=3) p-value, and MERGE / SEPARATE verdict.


### 5. Calibrate `merge_threshold`

1. Hit `n` to jump from one L2 boundary to the next.
2. Select **±2** layout to see the merging chain across 5 envelopes.
3. The dense zone is **L1_0008** (14 L2 envelopes) — start there.
4. Slide **merge τ**. The MERGE / SEPARATE verdict updates live.
5. Switch color mode to **family** to see whether bands are tracking
   family LD (FAM-LD?) instead of inversions.
6. Try K=2 if K=3 looks unstable on small envelopes.

Record the chosen value in `docs/calibration_LG28.md`.


## What the colours mean

**Sim_mat panel (PDF-style):**

| triangle | meaning | scale |
|---|---|---|
| **Lower** (j < i) | similarity | white → mint → sand → ember → maroon (q05 / mid / q95) |
| **Upper** (j > i) | local diagonal-distance Z | blue → white → red (symmetric ±z_max) |
| **Diagonal** | window-self | thin yellow line `#E8C547` |

**Envelope rectangles:**
| color | meaning |
|---|---|
| `#0042FF` blue | L1 envelope |
| `#00E676` green | L2 envelope |

**Boundary peaks on Z-track:**
| color | status |
|---|---|
| `#C8102E` red | STABLE_BLUE |
| `#E07B3F` orange | MARGINAL |
| `#999999` grey | DECAYS |
| `#7A4FBF` purple | EDGE |
| dark grey | DEDUP |

**K-means group dots (cluster color mode):**
| color | group |
|---|---|
| blue `#4fa3ff` | group 0 (lo PC1, ≈ hom1) |
| grey `#b8b8b8` | group 1 (mid PC1, ≈ het) |
| amber `#f5a524` | group 2 (hi PC1, ≈ hom2) |

**Family dots (family color mode, ported from STEP_D16):**
| tier | color |
|---|---|
| hub (n≥4) | one of 30 distinct Tol-qualitative colors, sorted by family size descending |
| small (n=2-3) | `#cbd5e1` light grey |
| singleton (n=1) | `#dde3eb` lighter grey |
| unmatched (-1) | `#94a3b8` neutral grey |


## Keyboard

| key | action |
|---|---|
| `←` `→` | step 1 window |
| Shift+`←` `→` | step 20 windows |
| `n` `p` | jump to next / previous L2 envelope |
| `space` | play/pause |


## Plumbing notes

- TSV columns are read by exact header name; the v3 exporter is
  schema-aware (parent_l1_id triggers L2 branch).
- Window indices in TSVs are 1-based; scrubber 0-based internally,
  1-based on display.
- Hungarian alignment: brute-force enumeration of K! permutations
  (K ≤ 3, ≤ 6 perms). Maximizes diagonal sum.
- Significance: Fisher's exact for K=2, chi-square (Wilson–Hilferty
  approx) for K=3.
- Verdict: MERGE if `concord ≥ merge_τ`, else SEPARATE. Significance
  is informational; concordance is primary.
- In ±2 layout, every neighbor's labels are Hungarian-aligned against
  the **focal** envelope's labels, so colors are consistent across
  the whole 5-column strip.
- PDF heatmap styling uses the exact palette from
  `STEP_D17c_overlay_plot_L2_v8.R`:
  `#F8F8F8 / #A8DBC2 / #F2DC78 / #E08838 / #7E1F1F` for similarity,
  `#2C5AA0 / #FAFAFA / #B22222` for diverging Z.
- Family logic ported verbatim from `STEP_D16_tree_node_pca_facets_v4.R`
  (v3.2): bamlist Ind→CGA remap, union-find on pairs with path-halving,
  family_id=1 is biggest component, fam_purity is size-weighted mean
  of within-cluster top-family fraction (unmatched -1 excluded from
  the ratio but kept in cluster size).


## Tested

Test files (dev environment, not shipped):
- `test_l3.js` — 8 unit tests (kmeans, alignLabels, Fisher, chi², full clusterL2)
- `test_integration.js` — integration vs real LG28 envelope structure
- `test_layout.js` — 3 tests for L3 layout helpers
- `test_v31.js` — 6 tests for PDF colors, multi-scale resolution, trackedN
- `test_v32.js` — 8 tests for family palette tiers, fam_purity, color mode dispatch

All real tests pass. PDF color anchors and family-purity arithmetic
verified against expected values from STEP_D17c and STEP_D16.


## Files in 01_input/

| file | rows | what |
|---|---|---|
| `C_gar_LG28_d17L1_envelopes.tsv` | 11 | L1 segments |
| `C_gar_LG28_d17L1_boundaries.tsv` | 24 | L1 candidate peaks |
| `C_gar_LG28_d17L1_boundary_score_curve.tsv` | 4302 | per-window L1 boundary score |
| `C_gar_LG28_d17L2_envelopes.tsv` | 48 | L2 segments inside 10 L1 parents |
| `C_gar_LG28_d17L2_boundaries.tsv` | 134 | L2 peaks |
| `C_gar_LG28_d17L2_quadrant_audit.tsv` | 134 | quadrant-test audit |
| `C_gar_LG28_d17L2_segment_stats.tsv` | 10 | per-L1 Ward stats |

The sim_mat RDS files (nn40/nn80/nn160/nn320) and the relatedness files
(bamlist, pairs) live in your project tree — pass their paths to the
exporter directly. They're too large to ship in the bundle.


## Next step

After calibrating `merge_threshold`: write `STEP_D17_multipass_L3_v1.R`
— translation of the JS verdict logic into batch R. Single-chromosome,
no SLURM array, K=2 first.
