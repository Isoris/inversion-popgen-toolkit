# MODULE 2B — Population Structure & Ancestry Visualization

## Methods

### Extended ancestry inference (A01–A03)

To comprehensively characterize population substructure in the hatchery cohort, ancestry proportions were estimated across a wider parameter space than the initial MODULE_2A discovery pass. NGSadmix was run for K = 2 through K = 20 with 5 independent random seeds per K value, using BEAGLE genotype likelihood files at two thinning distances (500 bp and 1,000 bp) for both the full 226-sample set and the pruned unrelated set (81 samples after first-degree kinship removal). Per-chromosome analyses were additionally run at thin-500 for the full sample set to assess chromosome-level structure consistency.

Task lists were generated programmatically (A01) and submitted as SLURM array jobs for NGSadmix (A02) and evalAdmix (A03) with dependency chaining to ensure evalAdmix ran only after its corresponding NGSadmix job completed. EvalAdmix was run with 5 EM iterations and a tolerance of 0.05 to assess model fit by computing pairwise residual correlation matrices between observed and expected genotype likelihoods.

### Best-seed selection and ancestry export (A04)

For each K value within each analysis scope (whole-genome or per-chromosome) and thinning panel, the optimal seed was selected by a hierarchical decision rule: (1) highest log-likelihood, (2) lowest mean absolute off-diagonal evalAdmix residual, (3) lowest maximum absolute residual, and (4) lowest seed number as deterministic tie-break. If evalAdmix outputs were missing for some seeds, selection proceeded on log-likelihood alone with a warning flag.

A stable 20-color palette was assigned to ancestry clusters, and per-sample dominant cluster assignments with full Q-vectors were exported as the canonical ancestry reference table. Selected best-seed outputs were auto-registered to the project sample registry for cross-module consumption.

### K-hierarchy merge tree (A05)

To characterize the nested founder lineage structure of the hatchery population, a K-hierarchy merge tree was constructed by matching ancestry components between consecutive K values. For each transition from K to K−1, components were matched by maximum Q-vector overlap to identify which clusters merge. The resulting tree reveals the hierarchical relationships among inferred founder lineages, from coarse (K = 2) to fine (K = 20) resolution.

### Admixture composite figure (B01)

A single tall composite figure was produced per thinning panel, displaying from top to bottom: a relatedness classification strip, followed by evalAdmix residual strip and admixture barplot pairs for each K = 2 through K = 20. Samples were ordered by dominant ancestry at the selected best K, with the best-K row highlighted. This figure provides a complete overview of how population structure emerges across increasing model complexity.

### Kinship network visualization (B02)

Pairwise kinship relationships were visualized as a network graph using ggraph, with nodes colored by dominant ancestry at best K, sized by degree (number of related pairs), and darkness scaled by maximum Q-value. Edges were drawn as solid lines for duplicate/MZ twin pairs (θ > 0.354) and dashed lines for first-degree pairs (θ > 0.177). Connected component hulls provided subtle grouping.

### Kinship heatmap (B03)

A 226 × 226 pairwise theta matrix was displayed as an annotated heatmap, with samples ordered by PCAngsd hierarchical clustering. Annotation bars along both axes indicated dominant ancestry at best K and retained/pruned status. Horizontal reference lines marked KING-based theta thresholds for relatedness classification.

### PCA with kinship overlay (B04)

Principal components (PC1 vs PC2 from PCAngsd) were plotted with points colored by dominant ancestry and sized by Q-value. First-degree kinship edges were overlaid as line segments connecting related sample pairs. Convex hulls delineated ancestry groups, providing a joint view of genetic structure and familial relationships.

### K-hierarchy Sankey diagram (B05)

The merge tree from A05 was visualized as an alluvial/Sankey diagram with K = 2 at the top and K = 20 at the bottom. Each node represented an ancestry component at a given K, with width proportional to sample count. Flows between adjacent K values showed component merges, colored by component identity at best K.

### Relatedness summary (B06)

A four-panel composite summarized the kinship pruning process: (A) theta distribution histogram with threshold lines for each relatedness class, (B) category barplot with pair counts per class, (C) admixture barplot at best K with relatedness strip, and (D) before/after sample counts showing the pruning outcome.

### evalAdmix diagnostic (B07)

Per-K diagnostic panels displayed the upper-triangle residual correlation heatmap and a per-sample burden barplot (mean absolute off-diagonal residual). These diagnostics guided best-K selection and identification of poorly modeled samples or ancestry components.

## Key Parameters Summary

| Parameter | Value | Used in |
|-----------|-------|---------|
| K range | 2–20 | A01–A05 |
| Seeds per K | 5 | A01–A03 |
| NGSadmix -minMaf | 0.05 | A02 |
| evalAdmix -nIts | 5 | A03 |
| evalAdmix -misTol | 0.05 | A03 |
| PCAngsd eigenvectors | 6 | A04 |
| PCAngsd EM iterations | 250 | A04 |
| KING theta thresholds | 0.354 / 0.177 / 0.0884 / 0.0442 | B02, B03, B06 |
| Best-seed rule | loglik → mean\|resid\| → max\|resid\| → seed | A04 |
| Thinning panels | 500, 1000 bp | A01–A04 |
| Sample sets | all (226), pruned (81) | A01–A04 |
