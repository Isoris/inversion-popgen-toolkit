#!/usr/bin/env python3
"""
Logic-parity self-test of the STEP40 classifier, ported to NumPy.

Rationale: we don't have R in this sandbox, but we can verify that the
classifier's *logic* (feature extraction + thresholded decision tree)
behaves correctly on synthetic composite/simple/technical/demographic
scenarios. If this passes, the R version is very likely correct too —
because the algorithms are mathematically identical (only the packaging
differs).

If any synthetic scenario fails, that's a real signal that either the
threshold set or the feature definition needs tuning.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple

rng = np.random.default_rng(42)

# ──────────────────────────────────────────────────────────────────────────
# PARAMETERS — must mirror STEP40_architecture_classifier.R PARAMS block
# ──────────────────────────────────────────────────────────────────────────
PARAMS = dict(
    min_markers_for_classification = 30,
    min_samples_per_group           = 5,
    top_marker_variance_fraction    = 0.5,
    marker_family_max_k             = 4,
    subinterval_n_bins              = 3,
    min_block_markers               = 5,
    thr_simple_axis_strong          = 0.35,
    thr_simple_axis_weak            = 0.15,
    thr_boundary_consistency_strong = 0.70,
    thr_boundary_consistency_weak   = 0.45,
    thr_monotonic_strong            = 0.75,
    thr_monotonic_weak              = 0.50,
    thr_within_substructure_high    = 0.30,
    thr_subinterval_discord_high    = 0.25,
    thr_technical_high              = 0.50,
    thr_hwe_p_ok                    = 0.01,
    thr_family_ld_corr_high         = 0.60,
)

# ──────────────────────────────────────────────────────────────────────────
# FEATURE HELPERS
# ──────────────────────────────────────────────────────────────────────────

def safe_var(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.var(x, ddof=1)) if x.size >= 2 else 0.0


def safe_corr(x, y):
    x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3: return np.nan
    xv, yv = x[ok], y[ok]
    if np.std(xv) == 0 or np.std(yv) == 0: return np.nan
    return float(np.corrcoef(xv, yv)[0, 1])


def mean_impute(M):
    """Column-wise mean imputation (ignores NaN)."""
    M = M.copy()
    col_means = np.nanmean(M, axis=0)
    col_means[~np.isfinite(col_means)] = 0.0
    for j in range(M.shape[1]):
        bad = ~np.isfinite(M[:, j])
        if bad.any():
            M[bad, j] = col_means[j]
    return M


def pc1_var_explained(M):
    """M: n_rows × n_cols. Return fraction of total variance on PC1."""
    if M is None or M.ndim != 2 or M.shape[0] < 3 or M.shape[1] < 2:
        return np.nan
    M = mean_impute(M)
    col_var = np.var(M, axis=0, ddof=1)
    keep = col_var > 1e-10
    if keep.sum() < 2: return np.nan
    M = M[:, keep]
    try:
        # Center, SVD
        Mc = M - M.mean(axis=0, keepdims=True)
        _, s, _ = np.linalg.svd(Mc, full_matrices=False)
        sd2 = (s ** 2) / (M.shape[0] - 1)
        total = np.sum(col_var[keep])
        return float(sd2[0] / total) if total > 0 else np.nan
    except np.linalg.LinAlgError:
        return np.nan


def monotonic_fraction(grp_means):
    """grp_means: n_markers × 3 matrix of (HOMO_1, HET, HOMO_2) means."""
    h1, he, h2 = grp_means[:, 0], grp_means[:, 1], grp_means[:, 2]
    ok = np.isfinite(h1) & np.isfinite(he) & np.isfinite(h2)
    if ok.sum() < 5: return np.nan
    up = (h1 <= he) & (he <= h2) & ((h2 - h1) > 0.05)
    down = (h1 >= he) & (he >= h2) & ((h1 - h2) > 0.05)
    return float(np.mean((up | down)[ok]))


def within_class_substructure(X_samp_mark, groups, min_n=5):
    """X_samp_mark: n_samples × n_markers. Mean PC1 VE inside each karyotype."""
    scores = []
    for g in ("HOMO_1", "HET", "HOMO_2"):
        idx = np.where(groups == g)[0]
        if len(idx) >= min_n:
            ve = pc1_var_explained(X_samp_mark[idx, :])
            if np.isfinite(ve):
                scores.append(ve)
    return float(np.mean(scores)) if scores else np.nan


def cluster_marker_families(X_mark_samp, max_k=4, top_n=200):
    """Return (k, cluster_labels, silhouette). Polarity-flip to common sign."""
    n_mark = X_mark_samp.shape[0]
    if n_mark < 10:
        return dict(k=1, cluster=np.ones(n_mark, int), silhouette=np.nan,
                    keep_idx=np.arange(n_mark))
    mv = np.nan_to_num(np.var(X_mark_samp, axis=1, ddof=1))
    keep = np.argsort(-mv)[:min(top_n, n_mark)]
    M = X_mark_samp[keep].copy()
    # mean-impute rows
    for j in range(M.shape[1]):
        bad = ~np.isfinite(M[:, j])
        if bad.any():
            M[bad, j] = np.nanmean(M[:, j]) if np.any(np.isfinite(M[:, j])) else 0.0
    M = np.nan_to_num(M)

    # Polarity-flip rows negative-correlated with row 0 (to collapse polarity-
    # inverted copies of the same partition; keep genuinely different
    # partitions distinct — they'll have low correlation in either direction)
    ref = M[0]
    for i in range(M.shape[0]):
        cr = safe_corr(M[i], ref)
        if np.isfinite(cr) and cr < 0:
            M[i] = 2 - M[i]

    # Correlation matrix between rows
    if np.var(M) == 0:
        return dict(k=1, cluster=np.ones(n_mark, int), silhouette=np.nan,
                    keep_idx=keep)
    R = np.corrcoef(M)
    R = np.nan_to_num(R, nan=0.0)
    D = 1 - R  # correlation distance

    # Hierarchical clustering (average linkage)
    try:
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        # scipy wants condensed dist
        np.fill_diagonal(D, 0)
        D = (D + D.T) / 2
        D_sym = np.clip(D, 0, None)
        cond = squareform(D_sym, checks=False)
        Z = linkage(cond, method="average")
    except Exception:
        return dict(k=1, cluster=np.ones(n_mark, int), silhouette=np.nan,
                    keep_idx=keep)

    best_k, best_sil, best_cl = 1, -np.inf, np.ones(len(keep), int)
    for k in range(2, max_k + 1):
        cl = fcluster(Z, t=k, criterion="maxclust")
        _, counts = np.unique(cl, return_counts=True)
        if counts.min() < 3: continue
        # Silhouette on distance matrix
        sil = silhouette_on_D(D_sym, cl)
        if np.isfinite(sil) and sil > best_sil:
            best_sil, best_k, best_cl = sil, k, cl

    if best_sil < 0.15:
        best_k, best_cl = 1, np.ones(len(keep), int)

    full_cl = np.ones(n_mark, int)
    full_cl[keep] = best_cl
    return dict(k=best_k, cluster=full_cl, silhouette=best_sil, keep_idx=keep)


def silhouette_on_D(D, labels):
    labels = np.asarray(labels)
    n = len(labels)
    uniq = np.unique(labels)
    if len(uniq) < 2: return np.nan
    a, b = np.zeros(n), np.full(n, np.inf)
    for i in range(n):
        own = (labels == labels[i])
        own[i] = False
        a[i] = D[i, own].mean() if own.any() else 0
        for c in uniq:
            if c == labels[i]: continue
            oth = (labels == c)
            if oth.any():
                b[i] = min(b[i], D[i, oth].mean())
    denom = np.maximum(a, b); denom[denom == 0] = 1e-10
    return float(np.mean((b - a) / denom))


def boundary_consistency(X_mark_samp, n_markers_use=100):
    if X_mark_samp.shape[0] < 5 or X_mark_samp.shape[1] < 5: return np.nan
    mv = np.nan_to_num(np.var(X_mark_samp, axis=1, ddof=1))
    top = np.argsort(-mv)[:min(n_markers_use, X_mark_samp.shape[0])]
    M = X_mark_samp[top].copy().astype(float)
    ref = M[0].copy()
    for i in range(1, M.shape[0]):
        cr = safe_corr(M[i], ref)
        if np.isfinite(cr) and cr < 0:
            M[i] = 2 - M[i]
    M = np.nan_to_num(M)
    # Spearman = Pearson on ranks
    from scipy.stats import rankdata
    R = np.apply_along_axis(rankdata, 1, M)
    C = np.corrcoef(R)
    np.fill_diagonal(C, np.nan)
    return float(np.nanmean(np.abs(C)))


def subinterval_discordance(X_mark_samp, positions, n_bins=3, min_m=10):
    n_mark = X_mark_samp.shape[0]
    if n_mark < n_bins * min_m: return np.nan
    lo, hi = np.min(positions), np.max(positions)
    if hi <= lo: return np.nan
    edges = np.linspace(lo, hi, n_bins + 1)
    bin_id = np.digitize(positions, edges[1:-1])  # 0..n_bins-1
    scores = np.full((X_mark_samp.shape[1], n_bins), np.nan)
    for b in range(n_bins):
        mi = np.where(bin_id == b)[0]
        if len(mi) < min_m: continue
        sub = X_mark_samp[mi].copy()
        for j in range(sub.shape[1]):
            bad = ~np.isfinite(sub[:, j])
            if bad.any():
                sub[bad, j] = np.nanmean(sub[:, j]) if np.any(np.isfinite(sub[:, j])) else 0.0
        sub = np.nan_to_num(sub)
        v = np.var(sub, axis=1, ddof=1)
        kk = v > 1e-10
        if kk.sum() < 2: continue
        try:
            Mc = sub[kk].T - sub[kk].T.mean(axis=0, keepdims=True)
            U, S, Vt = np.linalg.svd(Mc, full_matrices=False)
            scores[:, b] = U[:, 0] * S[0]
        except np.linalg.LinAlgError:
            continue
    ok_bins = [b for b in range(n_bins) if np.isfinite(scores[:, b]).sum() >= 5]
    if len(ok_bins) < 2: return np.nan
    # Sign-align to ref bin
    ref = ok_bins[0]
    for b in ok_bins[1:]:
        cr = safe_corr(scores[:, b], scores[:, ref])
        if np.isfinite(cr) and cr < 0:
            scores[:, b] = -scores[:, b]
    cors = []
    for i, bi in enumerate(ok_bins):
        for bj in ok_bins[i+1:]:
            cr = safe_corr(scores[:, bi], scores[:, bj])
            if np.isfinite(cr): cors.append(cr)
    if not cors: return np.nan
    return 1 - float(np.mean(cors))


def hwe_gof(n_h1, n_het, n_h2):
    n = n_h1 + n_het + n_h2
    if n < 10: return dict(p_value=np.nan, chisq=np.nan, p_allele=np.nan)
    p = (2 * n_h1 + n_het) / (2 * n); q = 1 - p
    e = np.array([n * p * p, 2 * n * p * q, n * q * q])
    o = np.array([n_h1, n_het, n_h2])
    keep = e > 0
    if keep.sum() < 2:
        return dict(p_value=np.nan, chisq=np.nan, p_allele=p)
    chi = np.sum((o[keep] - e[keep]) ** 2 / e[keep])
    from scipy.stats import chi2
    pv = 1 - chi2.cdf(chi, df=1)
    return dict(p_value=float(pv), chisq=float(chi), p_allele=float(p))


def technical_score(sites_pos, interval_start, interval_end,
                     gap_intervals=None, density_bin_bp=50000):
    span = interval_end - interval_start
    if span <= 0: return 0.0
    n_bins = max(1, int(np.ceil(span / density_bin_bp)))
    edges = np.linspace(interval_start, interval_end, n_bins + 1)
    bin_ids = np.digitize(sites_pos, edges[1:-1])
    counts = np.bincount(bin_ids, minlength=n_bins)[:n_bins]
    dropout = float(np.mean(counts == 0))
    nz = counts[counts > 0]
    cv = float(np.std(nz) / np.mean(nz)) if nz.size > 2 else 0.0
    cv_n = min(cv / 2, 1.0)
    gap_frac = 0.0
    if gap_intervals:
        ov = 0
        for (gs, ge) in gap_intervals:
            lo2, hi2 = max(interval_start, gs), min(interval_end, ge)
            if hi2 > lo2: ov += hi2 - lo2
        gap_frac = ov / span
    return min(0.5 * dropout + 0.3 * cv_n + 0.2 * gap_frac, 1.0)


def family_ld_correlation(sample_scores, kinship_pcs):
    if kinship_pcs is None or not kinship_pcs: return np.nan
    # kinship_pcs expected format: {"PC1": {sample_id: value}, "PC2": {...}}
    # Use the first PC's samples as the reference set of samples
    first_pc = next(iter(kinship_pcs.values()))
    common = [s for s in sample_scores if s in first_pc]
    if len(common) < 10: return np.nan
    s = np.array([sample_scores[c] for c in common])
    max_abs = 0
    for pc_name, pc_vec in kinship_pcs.items():
        if pc_name == "sample": continue
        pc_common = np.array([pc_vec[c] for c in common if c in pc_vec])
        if len(pc_common) < 10: continue
        cr = safe_corr(s[:len(pc_common)], pc_common)
        if np.isfinite(cr) and abs(cr) > max_abs:
            max_abs = abs(cr)
    return max_abs


# ──────────────────────────────────────────────────────────────────────────
# MAIN FEATURE EXTRACTION + CLASSIFICATION
# ──────────────────────────────────────────────────────────────────────────

def compute_features(X, positions, groups, sample_names,
                      interval_start, interval_end,
                      gap_intervals=None, kinship_pcs=None):
    """X: markers × samples. Returns dict of features."""
    # Informative markers
    mv = np.nan_to_num(np.var(X, axis=1, ddof=1))
    thr = np.quantile(mv, 1 - PARAMS["top_marker_variance_fraction"])
    info = np.where(mv >= thr)[0]
    if len(info) < 10:
        info = np.argsort(-mv)[:min(30, X.shape[0])]
    Xi = X[info]
    pos_i = positions[info]

    # F1: simple axis
    f1 = pc1_var_explained(Xi.T)

    # F2: boundary consistency
    f2 = boundary_consistency(Xi, n_markers_use=min(100, Xi.shape[0]))

    # F3: marker family k
    fam = cluster_marker_families(Xi, max_k=PARAMS["marker_family_max_k"])
    f3_k, f3_sil = fam["k"], fam["silhouette"]

    # F4: subinterval discordance
    f4 = subinterval_discordance(Xi, pos_i, n_bins=PARAMS["subinterval_n_bins"])

    # F5: within-class substructure (using full X)
    f5 = within_class_substructure(X.T, np.asarray(groups),
                                    min_n=PARAMS["min_samples_per_group"])

    # F6: HWE
    n_h1 = int(np.sum(np.asarray(groups) == "HOMO_1"))
    n_het = int(np.sum(np.asarray(groups) == "HET"))
    n_h2 = int(np.sum(np.asarray(groups) == "HOMO_2"))
    f6 = hwe_gof(n_h1, n_het, n_h2)["p_value"]

    # F7: monotonic fraction
    grp_means = np.full((X.shape[0], 3), np.nan)
    for j, g in enumerate(("HOMO_1", "HET", "HOMO_2")):
        idx = np.where(np.asarray(groups) == g)[0]
        if len(idx) >= PARAMS["min_samples_per_group"]:
            grp_means[:, j] = np.nanmean(X[:, idx], axis=1)
    f7 = monotonic_fraction(grp_means)

    # F8: polarity block count
    delta = np.nan_to_num(grp_means[:, 2] - grp_means[:, 0])
    hw = PARAMS["min_block_markers"] // 2
    sm_sign = np.ones(len(delta), int)
    for i in range(len(delta)):
        lo2, hi2 = max(0, i - hw), min(len(delta) - 1, i + hw)
        local = delta[lo2:hi2 + 1]
        w = np.abs(local)
        pos_w = np.sum(w[local > 0]); neg_w = np.sum(w[local < 0])
        sm_sign[i] = 1 if pos_w >= neg_w else -1
    f8 = int(np.sum(np.abs(np.diff(sm_sign)) > 0) + 1)

    # F9: technical
    f9 = technical_score(positions, interval_start, interval_end,
                          gap_intervals=gap_intervals)

    # F10: family-LD correlation
    # Use global PC1 on informative markers
    Mi = mean_impute(Xi)
    v = np.var(Mi, axis=1, ddof=1)
    kk = np.where(v > 1e-10)[0]
    if len(kk) >= 2:
        Mc = Mi[kk].T - Mi[kk].T.mean(axis=0, keepdims=True)
        try:
            U, S, _ = np.linalg.svd(Mc, full_matrices=False)
            pc_global = U[:, 0] * S[0]
        except np.linalg.LinAlgError:
            pc_global = np.full(X.shape[1], np.nan)
    else:
        pc_global = np.full(X.shape[1], np.nan)
    pc_map = {sample_names[i]: pc_global[i] for i in range(len(sample_names))}
    f10 = family_ld_correlation(pc_map, kinship_pcs) if kinship_pcs else np.nan

    return dict(
        simple_axis_score=f1,
        boundary_consistency=f2,
        marker_family_k=f3_k,
        marker_family_silhouette=f3_sil,
        subinterval_discordance=f4,
        within_class_substructure=f5,
        hwe_p_value=f6,
        monotonic_fraction=f7,
        polarity_block_count=f8,
        technical_score=f9,
        family_ld_correlation=f10,
        n_markers=X.shape[0], n_samples=X.shape[1],
        n_hom1=n_h1, n_het=n_het, n_hom2=n_h2,
    )


def classify_architecture(f):
    def isT(v):
        try: return bool(np.isfinite(v)) and bool(v)
        except Exception: return False

    # Technical gate
    if f["technical_score"] is not None and np.isfinite(f["technical_score"]) \
            and f["technical_score"] >= PARAMS["thr_technical_high"]:
        return dict(cls="TECHNICAL", tier="D",
                    reason=f"technical_score={f['technical_score']:.2f} >= {PARAMS['thr_technical_high']:.2f}")

    # Demographic gate
    if (f["family_ld_correlation"] is not None
            and np.isfinite(f["family_ld_correlation"])
            and f["family_ld_correlation"] >= PARAMS["thr_family_ld_corr_high"]
            and np.isfinite(f["monotonic_fraction"])
            and f["monotonic_fraction"] < PARAMS["thr_monotonic_weak"]):
        return dict(cls="DEMOGRAPHIC_BLOCK", tier="C",
                    reason=f"family_ld={f['family_ld_correlation']:.2f}; monotonic={f['monotonic_fraction']:.2f}")

    sas = np.isfinite(f["simple_axis_score"]) and f["simple_axis_score"] >= PARAMS["thr_simple_axis_strong"]
    saw = np.isfinite(f["simple_axis_score"]) and f["simple_axis_score"] >= PARAMS["thr_simple_axis_weak"]
    bcs = np.isfinite(f["boundary_consistency"]) and f["boundary_consistency"] >= PARAMS["thr_boundary_consistency_strong"]
    bcw = np.isfinite(f["boundary_consistency"]) and f["boundary_consistency"] >= PARAMS["thr_boundary_consistency_weak"]
    ms  = np.isfinite(f["monotonic_fraction"]) and f["monotonic_fraction"] >= PARAMS["thr_monotonic_strong"]
    mw  = np.isfinite(f["monotonic_fraction"]) and f["monotonic_fraction"] >= PARAMS["thr_monotonic_weak"]
    hwe_ok = np.isfinite(f["hwe_p_value"]) and f["hwe_p_value"] >= PARAMS["thr_hwe_p_ok"]
    wh  = np.isfinite(f["within_class_substructure"]) and f["within_class_substructure"] >= PARAMS["thr_within_substructure_high"]
    sih = np.isfinite(f["subinterval_discordance"]) and f["subinterval_discordance"] >= PARAMS["thr_subinterval_discord_high"]
    mult = f["marker_family_k"] >= 2

    # SIMPLE_STRONG — strict: no substructure, no multi-family, no subinterval
    if sas and bcs and ms and not mult and not wh and not sih:
        tier = "A" if hwe_ok else "B"
        return dict(cls="SIMPLE_STRONG", tier=tier,
                    reason="axis + boundary + monotonic all strong; single family")

    # COMPOSITE_OVERLAP — strongest composite signal: distinct sample partitions
    # across the interval (left vs right disagree) + multiple marker families
    if sih and mult:
        return dict(cls="COMPOSITE_OVERLAP", tier="B",
                    reason=f"subinterval={f['subinterval_discordance']:.2f} + families={f['marker_family_k']}")

    # COMPOSITE_INTERNAL — main axis exists but within-class substructure
    # (one karyotype is internally heterogeneous) OR multiple marker families
    # without left/right disagreement. Check BEFORE SIMPLE_WEAK.
    if saw and (wh or mult):
        return dict(cls="COMPOSITE_INTERNAL", tier="B",
                    reason=f"axis={f['simple_axis_score']:.2f} within={f['within_class_substructure']:.2f} families={f['marker_family_k']}")

    # SIMPLE_WEAK — 3-state visible but low contrast / noisy (fallthrough)
    if saw and mw and not mult and not sih:
        return dict(cls="SIMPLE_WEAK", tier="B",
                    reason=f"axis={f['simple_axis_score']:.2f} monotonic={f['monotonic_fraction']:.2f}")

    return dict(cls="UNRESOLVED", tier="C", reason="mixed evidence")


# ──────────────────────────────────────────────────────────────────────────
# SCENARIO GENERATORS
# ──────────────────────────────────────────────────────────────────────────

N_SAMP = 226
N_MARK = 300


def make_clean_matrix(gts, n_markers=N_MARK, noise=0.1, signal=1.0):
    X = np.zeros((n_markers, len(gts)))
    for s, g in enumerate(gts):
        X[:, s] = g + rng.normal(0, noise, n_markers)
    X = 1 + (X - 1) * signal
    return np.clip(X, 0, 2)


def gts_60_106_60():
    return np.array([0]*60 + [1]*106 + [2]*60)


def scenario_simple_strong():
    gts = gts_60_106_60()
    return dict(X=make_clean_matrix(gts, noise=0.08, signal=0.9),
                gts=gts, expected="SIMPLE_STRONG")


def scenario_composite_internal():
    gts = gts_60_106_60()
    X = make_clean_matrix(gts, noise=0.1, signal=0.8)
    # Split HOMO_1 into two subpops at 80 random markers
    h1 = np.where(gts == 0)[0]
    half = h1[::2]
    marks = rng.choice(N_MARK, 80, replace=False)
    X[np.ix_(marks, half)] += 0.9
    X = np.clip(X, 0, 2)
    return dict(X=X, gts=gts, expected="COMPOSITE_INTERNAL")


def scenario_composite_overlap():
    gts1 = gts_60_106_60()
    gts2 = rng.permutation(gts1)
    n_left = N_MARK // 2
    X_left = make_clean_matrix(gts1, n_markers=n_left, noise=0.1, signal=0.85)
    X_right = make_clean_matrix(gts2, n_markers=N_MARK - n_left, noise=0.1, signal=0.85)
    X = np.vstack([X_left, X_right])
    return dict(X=X, gts=gts1, expected="COMPOSITE_OVERLAP")


def scenario_demographic():
    n_fam = 5
    fam_id = np.array([i % n_fam for i in range(N_SAMP)])
    rng.shuffle(fam_id)
    gts = rng.permutation(gts_60_106_60())
    X = rng.uniform(0.4, 1.6, (N_MARK, N_SAMP))
    for m in range(N_MARK):
        fam_eff = rng.normal(0, 0.5, n_fam)
        X[m] += fam_eff[fam_id]
    X = np.clip(X, 0, 2)
    kpc1 = fam_id.astype(float) + rng.normal(0, 0.2, N_SAMP)
    return dict(X=X, gts=gts, expected="DEMOGRAPHIC_BLOCK", kinship_pc1=kpc1)


def scenario_technical():
    gts = gts_60_106_60()
    X = make_clean_matrix(gts, noise=0.2, signal=0.4)
    return dict(X=X, gts=gts, expected="TECHNICAL", patchy=True)


# ──────────────────────────────────────────────────────────────────────────
# RUNNER
# ──────────────────────────────────────────────────────────────────────────

def run_scenario(name, scen, interval=(15_000_000, 18_000_000), seed=None):
    # Reseed per-scenario so each test is deterministic regardless of what
    # other scenarios have consumed from the shared RNG beforehand.
    global rng
    if seed is not None:
        rng = np.random.default_rng(seed)
    print()
    print("=" * 70)
    print(f"Scenario: {name}   (expected = {scen['expected']})")
    print("=" * 70)

    X = scen["X"]; gts = scen["gts"]
    groups = np.array(["HOMO_1", "HET", "HOMO_2"])[gts]
    sample_names = [f"S{i:03d}" for i in range(X.shape[1])]

    interval_start, interval_end = interval

    if scen.get("patchy"):
        # Technical scenario: all SNPs packed into 2 tight clusters,
        # leaving the large majority of the interval empty.
        # The interval spans 3 Mb; SNPs occupy only ~300 kb of it.
        half = N_MARK // 2
        positions = np.sort(np.concatenate([
            rng.uniform(interval_start,              interval_start + 150_000, half),
            rng.uniform(interval_end - 150_000,      interval_end,             N_MARK - half),
        ]))
    else:
        positions = np.sort(rng.uniform(interval_start, interval_end, N_MARK))

    kpcs = None
    if "kinship_pc1" in scen:
        kpcs = {"PC1": {sample_names[i]: float(scen["kinship_pc1"][i]) for i in range(N_SAMP)},
                "PC2": {sample_names[i]: float(rng.normal()) for i in range(N_SAMP)}}

    feat = compute_features(X, positions, groups, sample_names,
                             interval_start, interval_end, kinship_pcs=kpcs)
    for k, v in feat.items():
        if isinstance(v, float):
            print(f"  {k:30s} {v:.4f}" if np.isfinite(v) else f"  {k:30s} NA")
        else:
            print(f"  {k:30s} {v}")
    cls = classify_architecture(feat)
    got = cls["cls"]
    status = "MATCH" if got == scen["expected"] else "MISMATCH"
    print(f"\nResult:   {cls['cls']}-{cls['tier']}   [{cls['reason']}]")
    print(f"Expected: {scen['expected']}   -> {status}")
    return got == scen["expected"], feat, cls


if __name__ == "__main__":
    print("#" * 72)
    print("#  STEP40 classifier self-test (NumPy port)")
    print("#" * 72)

    results = []
    results.append(("simple_strong",       *run_scenario("SIMPLE_STRONG",       scenario_simple_strong())))
    results.append(("composite_internal",  *run_scenario("COMPOSITE_INTERNAL",  scenario_composite_internal())))
    results.append(("composite_overlap",   *run_scenario("COMPOSITE_OVERLAP",   scenario_composite_overlap())))
    results.append(("demographic",         *run_scenario("DEMOGRAPHIC_BLOCK",   scenario_demographic())))
    results.append(("technical",           *run_scenario("TECHNICAL",           scenario_technical())))

    print()
    print("#" * 72)
    print("#  SUMMARY")
    print("#" * 72)
    n_match = sum(1 for r in results if r[1])
    print(f"Accuracy: {n_match}/{len(results)} scenarios matched")
    for name, ok, _, cls in results:
        print(f"  {name:25s} {'PASS' if ok else 'FAIL'}  (got {cls['cls']}-{cls['tier']})")
