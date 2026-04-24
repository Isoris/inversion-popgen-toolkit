# MODULE_5B — Inversion Follow-up

Sample-structure-first candidate interpretation framework. Takes discovery outputs (MODULE_5A) and characterizes each inversion candidate through multi-layer analysis.

## Core engine (current_followup/, STEP20–STEP36)

### New sample-structure-first steps

**STEP20 — Window Profile + Vector/Profile Engine** (v4, 587 lines)

The central method. For each candidate:
- Returns to raw marker-level dosage from STEP08 (bypasses STEP09 compression)
- Maintains three parallel raw marker-vector encodings as first-class objects:
  - `candidate_raw_minor_vectors.tsv.gz` — continuous minor dosage
  - `candidate_raw_major_vectors.tsv.gz` — continuous major dosage
  - `candidate_012_marker_matrix.tsv.gz` — discretized 0/1/2
- Computes per-sample profile strings (012, AB encoding)
- Computes pairwise distances independently per encoding (Manhattan minor, Manhattan major, Hamming 012)
- Compares each sample to group consensus profiles
- Assigns quality tiers (core / peripheral / ambiguous) based on reference match margin
- Groups windows by induced sample structure (ARI-based clustering)

Two modes: MODE 1 (dosage-only, works now) / MODE 2 (dosage + Clair3, compatible outputs).

**STEP20b — Anchor Geometry Module** (411 lines)

Anchor-centered (u,v) transformation:
- Dense core points define left/middle/right anchors
- All samples projected into anchor space
- Per-sample: u, v, center_distance, broad_axis_position, left_affinity, right_affinity, branch_balance, ambiguity_score
- Diagnostic splits: VIEW A (all), VIEW B (homo-only), VIEW C (het-only)

**STEP20c — Gradient-Marker Experiments** (324 lines)

Per-marker Spearman correlation with anchor-derived gradients:
- MODE A: all samples × {u, center_distance, branch_balance}
- MODE B: HET-only × same gradients
- Marker classification: core_inversion_axis / het_center_distance / branch_specific / uninformative

### Existing follow-up steps (enhanced by richer upstream)

```
STEP21   Stripe-aware state assignment (GMM/k-means + rotation)
STEP22   Within-stripe DBSCAN subclustering
STEP23   Ancestry/Q association
STEP24   Het + marker support
STEP25B  Master interpretation table v2
STEP29   Coherence + L1/L2 polarity
STEP30   Window-state trajectories
STEP31   Diagnostic figures
STEP32   Multiscale windows
STEP33   Within-stripe analysis + AB asymmetry
STEP34   Anchor-based HET decomposition
STEP35   Multi-mode harmonized heatmaps
STEP36   Clair3 local-signature support (when available)
```

### Figure steps

```
STEP26   15-figure catalogue | STEP27  Faceted PCA | STEP28  Marker heatmap
STEP31   Diagnostic figures  | STEP35  Harmonized heatmaps
```

## Legacy follow-up (legacy_followup/, STEP14–STEP18)

Classic follow-up backbone. STEP17c exports the shared contrast-group interface for LD/FST/HOBS.

## Running

```bash
# Single candidate, full pipeline
bash run_candidate_followup_v6.sh 42 full

# All candidates via SLURM
sbatch run_candidate_followup.slurm

# Core steps only (STEP10c → STEP20 → STEP20b → STEP20c → STEP21–25)
bash run_candidate_followup_v6.sh all core
```

## Key outputs per candidate

| Step | Output | Description |
|------|--------|-------------|
| STEP20 | `candidate_raw_minor_vectors.tsv.gz` | First-class raw marker vectors (minor) |
| STEP20 | `candidate_raw_major_vectors.tsv.gz` | First-class raw marker vectors (major) |
| STEP20 | `candidate_012_marker_matrix.tsv.gz` | Discretized 012 marker matrix |
| STEP20 | `candidate_pairwise_profile_distance_*.tsv.gz` | 3 distance matrices (per encoding) |
| STEP20 | `candidate_best_profile_assignment.tsv` | Quality tiers from reference matching |
| STEP20b | `candidate_anchor_projection.tsv` | Per-sample geometric features |
| STEP20b | `candidate_split_reanalysis_summary.tsv` | Homo-only / het-only diagnostics |
| STEP20c | `candidate_marker_gradient_all.tsv` | Marker × gradient correlations |
| STEP20c | `candidate_top_gradient_markers.tsv` | Classified top markers |
