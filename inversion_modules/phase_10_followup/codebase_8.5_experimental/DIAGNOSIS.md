# Diagnosis and Fix Plan — Inversion Follow-Up Pipeline

## Response to Prompts 1/3, 2/3, 3/3

---

## 1. Root Cause Diagnosis

### A. FIG_C10 Topology Track: appears blank/grey

**Root cause confirmed:** `geom_tile(height=0.8)` with 365 windows in a 10-inch
plot creates tiles ~0.027 inches wide. The fill color IS being applied (the code
is correct about factor levels after the v4 fix), but the tiles are so narrow that
the fill is visually indistinguishable from the grey panel background.

STEP32 confirms: "W=100 clean ~100%". So the topology classifier IS returning
`clean_3band` for most windows. The data is correct; the rendering fails.

**Fix delivered:** `STEP26_topology_fix.R` — uses `geom_rect()` with explicit
`xmin/xmax` from window coordinates, `color="grey30"` borders, and full alpha.

### B. FIG_C08 Heatmap: hex-code Position legend

**Root cause confirmed:** `columnAnnotation(Position = colorRamp2(...)(pos))` maps
200 unique genomic positions to 200 unique colors. ComplexHeatmap's legend engine
tries to show ALL of them as discrete entries, producing the hex-code wall.

The Het annotation has the SAME bug when using raw `col = list(Het = colorRamp2(...))`
because ComplexHeatmap converts the continuous mapping to discrete legend entries
when there are many unique values.

**Fix delivered:** `STEP28_heatmap_fix.R` — removes Position annotation entirely;
uses `anno_simple()` for Het which creates a proper continuous legend bar.

### C. FIG_C30/C35: HET coherence is always "intermediate" or worse

**Root cause:** The coherence metric computes agreement as "fraction of markers
where sample is closer to own group mean than to other-group mean." For true
heterokaryotypes, their dosage is ~1.0 (intermediate between 0 and 2), and the
HET group mean is also ~1.0. But the "other mean" is the average of HOMO_1 mean
(~0) and HOMO_2 mean (~2), which is also ~1.0. So HET samples are equidistant
from own and other means, giving agreement ≈ 0.5 — classified as "intermediate."

This is **expected biology, not a bug in the data**. But the threshold treats it
as ambiguous, which then makes all HET samples "peripheral" or "junk."

**Fix delivered:** `STEP29 v2` — uses group-specific thresholds:
- Homos: coherent ≥0.70, intermediate ≥0.45, discordant <0.45
- HET: coherent ≥0.55, intermediate ≥0.40, discordant <0.40

### D. FIG_C34 Core-only PCA: missing HET samples

**Root cause:** When all HET samples are "intermediate"/"peripheral" (see C above),
none pass the "core" tier. So the core-only PCA shows only HOMO_1 and HOMO_2.

**Fix:** With group-specific thresholds, HET samples with agreement ≥0.55 become
"core" and appear in the core-only PCA.

### E. FIG_C32/C36: 15,967 blocks from 36,500 markers

**Root cause:** Block definition = contiguous same-sign Δhom. With noisy ~5x WGS
dosage, adjacent markers often have small random sign flips in Δhom, creating
blocks of 2-3 markers each.

**Fix delivered:** `STEP29 v2` — L2 polarity uses a sliding window of
MIN_BLOCK_SIZE=5 markers, taking the weighted median sign. This merges tiny
noise-driven blocks into coherent longer blocks.

### F. STEP30 trajectories: all "moderate", none "stable"

**Root cause:** Window-local k-means produces labels that are not aligned across
windows. Window 1's "group 1" may correspond to HOMO_2, while window 2's "group 1"
may correspond to HOMO_1, depending on random k-means initialization and centroid
ordering.

**Current status:** This is a fundamental limitation. Fixing it properly requires
aligning per-window labels to the candidate-wide coarse groups (e.g., by comparing
each window's centroids to the candidate-wide HOMO_1/HET/HOMO_2 centroids).

**Recommendation:** Mark trajectory outputs as provisional until alignment is fixed.

---

## 2. Polarity Validity Assessment (Prompt 2)

### Is marker flipping biologically valid?

**Yes, under specific conditions.** Here is the rigorous reasoning:

The dosage at a biallelic site encodes the count of the ALT allele (0, 1, or 2).
Which allele is called "ALT" depends on the reference genome assembly. In an
inversion region, if the reference carries arrangement A, then arrangement B
carriers will have ALT-enriched sites at some positions and REF-enriched at others,
depending on the local allele structure.

The key insight: **for a true inversion, all markers within the inverted block
should show the same homokaryotype contrast direction** (HOMO_2 consistently
higher OR consistently lower than HOMO_1). If some markers show the opposite
direction, it is because the reference allele encoding is inconsistent — the
ALT allele sometimes tags arrangement A and sometimes arrangement B.

Flipping = re-orienting so that the same arrangement is consistently "high dosage"
across all markers. This is an **encoding convention**, not a biological alteration.

### When is flipping valid?

1. When the contrast is learned from stable samples (core homokaryotypes)
2. When the flipping rule is marker-level or block-level, not candidate-global
3. When the raw and harmonized versions are both preserved and displayed
4. When weak/ambiguous markers are flagged, not silently flipped

### When is flipping dangerous?

1. When polarity is learned from the same samples used to display the result
   (circularity risk — addressed by infer-on-core-homos, apply-to-all)
2. When nearly all markers get flipped (suggests the reference may carry the
   minor arrangement, which is fine but should be noted)
3. When blocks are so small that flipping is just noise-smoothing

### Anti-circularity safeguards

The current code infers polarity from ALL samples' group means, including the
samples that were assigned to groups based on the same PCA that polarity is
trying to validate. This is mildly circular.

**Better approach:** Infer polarity from core homokaryotypes only (the most
clearly assigned samples), then apply to all samples. This breaks the circularity
because core samples are defined by PCA geometry + coherence, not by polarity.

The L2 block consensus further stabilizes by requiring neighborhood agreement.

---

## 3. Proposed Figure Architecture (Prompt 3)

### Family A — Candidate Overview
- FIG_C02: PCA by het ✓ (working)
- FIG_C03: PCA by ancestry ✓ (working)
- FIG_C05: Summary card ✓ (working, pattern label fixed)
- FIG_C06B: Group counts ✓ (working)
- FIG_C07: Het boxplot ✓ (count label fixed)
- FIG_C07B: Het ridges ✓ (working)
- FIG_C16: Marker density ✓ (working)

### Family B — Stripe Assignment
- FIG_C04: Rotated PCA by subcluster ✓ (working, but DBSCAN over-fragments)
- FIG_C14: DBSCAN diagnostics ✓ (working, needs density contour alternate)
- NEW: PCA with density contours / ellipses per group
- NEW: PCA with centroid + boundary overlay

### Family C — Window Behavior
- FIG_C10: Topology track ✗ → FIX: use geom_rect with borders
- FIG_C09: Faceted local PCA ✓ (working)
- FIG_C11: Window agreement heatmap ✓ (working)
- NEW: Multiscale comparison panel

### Family D — Coherence + Quality
- FIG_C30: Coherence scatter ✓ (working, but HET threshold needs fix)
- FIG_C31: Quality tiers ✓ (working, same fix)
- FIG_C34: Core vs all PCA ✗ → FIX: include HET core
- FIG_C35: Coherence by stripe ✓ (working, needs HET threshold fix)

### Family E — Polarity + Harmonization
- FIG_C32: Polarity track ✓ (working, needs L1/L2 separation)
- FIG_C36: Block summary ✓ (working, needs L2 blocks)
- FIG_C08: Marker heatmap ✗ → FIX: remove Position legend
- Harmonized heatmap ✓ (working well with STEP35)
- NEW: Raw vs harmonized 2×2 quartet
- NEW: Δhom support track with confidence shading

### Family F — Synthesis
- Compact composite ✓ (working)
- Intermediate composite ✓ (working)
- NEW: Full diagnostic composite with all families

---

## 4. Implementation Priority Order

### Phase 1 — Fix broken fundamentals (this delivery)
1. ✅ STEP29 v2: group-specific coherence thresholds, L1/L2 polarity
2. ✅ FIG_C10 topology fix: geom_rect with borders
3. ✅ FIG_C08 heatmap fix: remove Position legend, use anno_simple for Het

### Phase 2 — Downstream consequences of STEP29 fix
4. Re-run STEP34 anchor resolution (will now include HET core samples)
5. Re-run STEP35 harmonized heatmaps (will use L2 polarity)
6. Re-run STEP25B interpretation table (will reflect corrected coherence)

### Phase 3 — New debug figures
7. Raw vs harmonized 2×2 quartet panel
8. PCA with density contours
9. Multiscale comparison panel
10. Δhom support track with L1/L2 comparison

### Phase 4 — Trajectory alignment fix
11. Align per-window k-means labels to candidate-wide coarse groups
12. Re-score trajectories with aligned labels
