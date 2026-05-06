// =============================================================================
// atlas_turn7.js — Turn 7 of chat A
// =============================================================================
//
// Two deliverables:
//
//   A. Page 7 (ancestry) per-group stratified Δ12 — companion to turn 6's
//      page 6 wiring. Uses the popstats live response when present, falls
//      back to static state.data.ancestry_window cohort-mean otherwise.
//
//   B. Q09b shelf-LD test — port of phase_5_qc_triage/R/q09b_shelf_ld_check.R
//      to atlas-side JS. Inputs: dosage chunk + per-sample karyotype labels +
//      shelf bp interval + N_BINS. Outputs: per-bin per-sample arrangement
//      score, between-bin Pearson correlation matrix, summary verdict.
//
// API (window.popgenTurn7):
//   .runShelfLDTest({ chunk, sample_invgt, shelf_start_bp, shelf_end_bp,
//                     n_bins?, diagnostic_loose?, diagnostic_strong? })
//   .computePerSnpStats(chunk, sample_invgt)
//   .scorePerBin(...)
//   .pearsonCorrelationMatrix(matrix)
//   .verdictFromCorrMatrix(corMat, n_bins)
//   .makeShelfLDPanel(opts)        — DOM element with controls + heatmap
//                                     canvas + summary line. Mounts onto
//                                     page 3 (candidate focus).
//   .makePage7AncestryStratifier() — DOM element to mount on page 7
//   .wirePage7DeltaStratified()    — convenience installer
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenTurn7 = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Q09b — shelf LD test
  // ===========================================================================
  //
  // Input shape (all required):
  //   chunk = {
  //     samples: string[],
  //     markers: [{marker_id, pos_bp, missingness?, diagnostic_score?}, ...],
  //     dosage:  number[][]   // dosage[m][s], -1 or null = NA
  //   }
  //   sample_invgt: Map<sample_id, 'Hom1'|'Het'|'Hom2'>   (string ID → label)
  //   shelf_start_bp, shelf_end_bp:  inclusive bounds
  //   n_bins (default 20)
  //   diagnostic_loose (default 0.3)  — threshold for being included in bin score
  //   diagnostic_strong (default 0.6) — reported in summary, not used in score
  //
  // Output:
  //   {
  //     ok: true,
  //     n_shelf_snps, n_diag_loose, n_diag_strong,
  //     per_snp: [{ marker_id, bp, mb, p_hom1, p_hom2, p_het, diagnostic, bin }],
  //     score_per_bin: number[][],       // [n_samples][n_bins], NaN where empty
  //     samples_in_score: string[],      // ordering used in score_per_bin
  //     cor_mat: number[][],             // [n_bins][n_bins], NaN when undefined
  //     bin_centers_mb: number[],
  //     summary: {
  //       within_first_half_corr, within_second_half_corr,
  //       between_halves_corr,    overall_mean_corr,
  //       verdict: 'SINGLE_INVERSION' | 'LIKELY_MULTIPLE_ARRANGEMENTS' |
  //                'AMBIGUOUS_OR_PARTIAL_LD'
  //     }
  //   }

  function runShelfLDTest(opts) {
    opts = opts || {};
    const {
      chunk, sample_invgt,
      shelf_start_bp, shelf_end_bp,
    } = opts;
    const n_bins = opts.n_bins || 20;
    const diag_loose  = (opts.diagnostic_loose  != null) ? opts.diagnostic_loose  : 0.3;
    const diag_strong = (opts.diagnostic_strong != null) ? opts.diagnostic_strong : 0.6;

    if (!chunk || !Array.isArray(chunk.markers) || !Array.isArray(chunk.dosage) ||
        !Array.isArray(chunk.samples)) {
      return { ok: false, error: 'invalid chunk shape' };
    }
    if (!(sample_invgt instanceof Map) && typeof sample_invgt !== 'object') {
      return { ok: false, error: 'sample_invgt must be a Map or {sample_id: label}' };
    }
    const groupOf = _normalizeInvgt(sample_invgt);
    if (!isFinite(shelf_start_bp) || !isFinite(shelf_end_bp) ||
        shelf_end_bp <= shelf_start_bp) {
      return { ok: false, error: 'invalid shelf interval' };
    }

    // 1. Index shelf SNPs
    const shelfSnpIdx = [];
    for (let mi = 0; mi < chunk.markers.length; mi++) {
      const m = chunk.markers[mi];
      if (!m || m.pos_bp == null) continue;
      if (m.pos_bp < shelf_start_bp || m.pos_bp > shelf_end_bp) continue;
      shelfSnpIdx.push(mi);
    }
    if (shelfSnpIdx.length === 0) {
      return { ok: false, error: 'no shelf SNPs' };
    }

    // 2. Map sample columns to group indices
    const hom1_cols = [], hom2_cols = [], het_cols = [];
    for (let si = 0; si < chunk.samples.length; si++) {
      const id = chunk.samples[si];
      const grp = groupOf.get(id);
      if (grp === 'Hom1' || grp === 'HOM1' || grp === 'H1/H1') hom1_cols.push(si);
      else if (grp === 'Hom2' || grp === 'HOM2' || grp === 'H2/H2') hom2_cols.push(si);
      else if (grp === 'Het'  || grp === 'HET'  || grp === 'H1/H2') het_cols.push(si);
    }
    if (hom1_cols.length < 2 || hom2_cols.length < 2) {
      return { ok: false, error: 'need ≥2 Hom1 and ≥2 Hom2 samples to score diagnostic SNPs' };
    }

    // 3. Per-SNP allele frequencies
    const per_snp = [];
    for (const mi of shelfSnpIdx) {
      const row = chunk.dosage[mi];
      const p_hom1 = _meanDosageOverCols(row, hom1_cols) / 2;
      const p_hom2 = _meanDosageOverCols(row, hom2_cols) / 2;
      const p_het  = (het_cols.length > 0)
                     ? _meanDosageOverCols(row, het_cols) / 2 : null;
      const diagnostic = (isFinite(p_hom1) && isFinite(p_hom2))
                          ? Math.abs(p_hom1 - p_hom2) : null;
      const bp = chunk.markers[mi].pos_bp;
      per_snp.push({
        marker_id: chunk.markers[mi].marker_id || ('m' + mi),
        chunk_idx: mi,
        bp, mb: bp / 1e6,
        p_hom1, p_hom2, p_het, diagnostic,
      });
    }

    // 4. Bin assignment
    const bin_edges_bp = [];
    for (let i = 0; i <= n_bins; i++) {
      bin_edges_bp.push(shelf_start_bp + (i / n_bins) * (shelf_end_bp - shelf_start_bp));
    }
    const bin_centers_mb = [];
    for (let i = 0; i < n_bins; i++) {
      bin_centers_mb.push((bin_edges_bp[i] + bin_edges_bp[i + 1]) / 2 / 1e6);
    }
    for (const s of per_snp) {
      let bi = Math.floor(((s.bp - shelf_start_bp) / (shelf_end_bp - shelf_start_bp)) * n_bins);
      if (bi >= n_bins) bi = n_bins - 1;
      if (bi < 0)       bi = 0;
      s.bin = bi;
    }

    // 5. Score per bin per sample
    // For each bin, collect diagnostic SNPs (diag > diag_loose), polarize each
    // by sign(p_hom2 - p_hom1) so dosage on the "Hom2-like" side is positive,
    // average across SNPs per sample.
    const n_samples = chunk.samples.length;
    const score_per_bin = [];
    for (let si = 0; si < n_samples; si++) {
      score_per_bin.push(new Array(n_bins).fill(NaN));
    }
    for (let b = 0; b < n_bins; b++) {
      const snps = per_snp.filter(s => s.bin === b &&
                    s.diagnostic != null && s.diagnostic > diag_loose);
      if (snps.length < 5) continue;
      // For each sample, mean of polarized dosage across diagnostic SNPs
      for (let si = 0; si < n_samples; si++) {
        let sum = 0, n = 0;
        for (const s of snps) {
          const v = chunk.dosage[s.chunk_idx][si];
          if (v == null || !isFinite(v) || v < 0) continue;   // NA
          const pol = Math.sign(s.p_hom2 - s.p_hom1);
          if (pol === 0) continue;   // tied, no polarity
          sum += pol * v;
          n++;
        }
        score_per_bin[si][b] = (n > 0) ? sum / n : NaN;
      }
    }

    // 6. Pairwise Pearson over bins
    const cor_mat = pearsonCorrelationMatrix(score_per_bin);

    // 7. Summary
    const summary = verdictFromCorrMatrix(cor_mat, n_bins);

    const n_diag_loose = per_snp.filter(s =>
      s.diagnostic != null && s.diagnostic > diag_loose).length;
    const n_diag_strong = per_snp.filter(s =>
      s.diagnostic != null && s.diagnostic > diag_strong).length;

    return {
      ok: true,
      n_shelf_snps: per_snp.length,
      n_diag_loose, n_diag_strong,
      per_snp,
      score_per_bin,
      samples_in_score: chunk.samples.slice(),
      cor_mat,
      bin_centers_mb,
      summary,
    };
  }

  function _normalizeInvgt(input) {
    if (input instanceof Map) return input;
    const m = new Map();
    if (input && typeof input === 'object') {
      for (const k of Object.keys(input)) m.set(String(k), input[k]);
    }
    return m;
  }

  function _meanDosageOverCols(row, cols) {
    if (!row || cols.length === 0) return NaN;
    let sum = 0, n = 0;
    for (const c of cols) {
      const v = row[c];
      if (v == null || !isFinite(v) || v < 0) continue;
      sum += v; n++;
    }
    return n > 0 ? sum / n : NaN;
  }

  // ---------------------------------------------------------------------------
  // Pearson correlation matrix with pairwise complete observations
  // ---------------------------------------------------------------------------
  // Input: matrix[n_rows][n_cols] (in Q09b: rows = samples, cols = bins)
  // Output: cor[n_cols][n_cols] (correlation between every pair of columns)
  //         NaN for pairs where < 2 paired-non-NA rows exist.

  function pearsonCorrelationMatrix(M) {
    if (!Array.isArray(M) || M.length === 0) return [];
    const n_cols = M[0].length;
    const n_rows = M.length;
    const cor = [];
    for (let i = 0; i < n_cols; i++) {
      const row = new Array(n_cols).fill(NaN);
      cor.push(row);
    }
    for (let a = 0; a < n_cols; a++) {
      cor[a][a] = 1.0;
      for (let b = a + 1; b < n_cols; b++) {
        let sumA = 0, sumB = 0, sumAB = 0, sumA2 = 0, sumB2 = 0, n = 0;
        for (let r = 0; r < n_rows; r++) {
          const va = M[r][a], vb = M[r][b];
          if (va == null || vb == null || isNaN(va) || isNaN(vb)) continue;
          sumA += va; sumB += vb;
          sumAB += va * vb;
          sumA2 += va * va; sumB2 += vb * vb;
          n++;
        }
        if (n < 2) continue;   // NaN
        const meanA = sumA / n, meanB = sumB / n;
        const cov = sumAB / n - meanA * meanB;
        const varA = sumA2 / n - meanA * meanA;
        const varB = sumB2 / n - meanB * meanB;
        if (varA <= 0 || varB <= 0) continue;
        cor[a][b] = cor[b][a] = cov / Math.sqrt(varA * varB);
      }
    }
    return cor;
  }

  function verdictFromCorrMatrix(cor, n_bins) {
    const mid = Math.floor(n_bins / 2);
    const a_bins = [], b_bins = [];
    for (let i = 0; i < mid; i++) a_bins.push(i);
    for (let i = mid; i < n_bins; i++) b_bins.push(i);
    function meanOfPairs(rows, cols, requireUpper) {
      let sum = 0, n = 0;
      for (const r of rows) {
        for (const c of cols) {
          if (requireUpper && c <= r) continue;
          const v = cor[r][c];
          if (v == null || isNaN(v)) continue;
          sum += v; n++;
        }
      }
      return n > 0 ? sum / n : NaN;
    }
    const within_a  = meanOfPairs(a_bins, a_bins, true);
    const within_b  = meanOfPairs(b_bins, b_bins, true);
    const between   = meanOfPairs(a_bins, b_bins, false);

    let sumAll = 0, nAll = 0;
    for (let i = 0; i < n_bins; i++) {
      for (let j = i + 1; j < n_bins; j++) {
        const v = cor[i][j];
        if (v != null && !isNaN(v)) { sumAll += v; nAll++; }
      }
    }
    const overall = nAll > 0 ? sumAll / nAll : NaN;

    let verdict;
    if (!isFinite(between)) {
      verdict = 'AMBIGUOUS_OR_PARTIAL_LD';
    } else if (between > 0.7) {
      verdict = 'SINGLE_INVERSION';
    } else if (between < 0.3) {
      verdict = 'LIKELY_MULTIPLE_ARRANGEMENTS';
    } else {
      verdict = 'AMBIGUOUS_OR_PARTIAL_LD';
    }

    return {
      within_first_half_corr:  within_a,
      within_second_half_corr: within_b,
      between_halves_corr:     between,
      overall_mean_corr:       overall,
      verdict,
    };
  }

  // ===========================================================================
  // makeShelfLDPanel — DOM panel with heatmap + summary line, mountable on
  // page 3 (candidate focus) or anywhere else
  // ===========================================================================

  function makeShelfLDPanel(opts) {
    opts = opts || {};
    if (typeof document === 'undefined') return null;
    _injectCSS();
    const root = document.createElement('div');
    root.setAttribute('data-popgen-shelf-ld', 'root');

    // Header
    const header = document.createElement('div');
    header.setAttribute('data-popgen-shelf-ld', 'header');
    const title = document.createElement('span');
    title.setAttribute('data-popgen-shelf-ld', 'title');
    title.textContent = 'Q09b shelf-LD test';
    header.appendChild(title);

    // Run button
    const runBtn = document.createElement('button');
    runBtn.setAttribute('data-popgen-shelf-ld', 'run');
    runBtn.textContent = '▶ run';
    runBtn.title = 'Run Q09b on the focal candidate';
    header.appendChild(runBtn);

    root.appendChild(header);

    // Status / summary line
    const status = document.createElement('div');
    status.setAttribute('data-popgen-shelf-ld', 'status');
    status.textContent = '— click ▶ run to compute —';
    root.appendChild(status);

    // Heatmap canvas
    const canvas = document.createElement('canvas');
    canvas.setAttribute('data-popgen-shelf-ld', 'heatmap');
    canvas.width = 360; canvas.height = 360;
    canvas.style.display = 'block';
    canvas.style.margin = '6px auto';
    root.appendChild(canvas);

    // Verdict line
    const verdict = document.createElement('div');
    verdict.setAttribute('data-popgen-shelf-ld', 'verdict');
    root.appendChild(verdict);

    // Wire run button — pulls the chunk + invgt from atlas state at click time
    runBtn.addEventListener('click', async () => {
      status.textContent = 'collecting inputs...';
      const inputs = (opts.inputProvider || _defaultInputProvider)();
      if (!inputs || !inputs.chunk) {
        status.textContent = 'no dosage chunk available';
        return;
      }
      status.textContent = 'computing...';
      const result = runShelfLDTest({
        chunk: inputs.chunk,
        sample_invgt: inputs.sample_invgt,
        shelf_start_bp: inputs.shelf_start_bp,
        shelf_end_bp: inputs.shelf_end_bp,
        n_bins: opts.n_bins || 20,
      });
      if (!result.ok) {
        status.textContent = 'error: ' + result.error;
        verdict.textContent = '';
        verdict.removeAttribute('data-verdict');
        return;
      }
      status.textContent = result.n_shelf_snps + ' shelf SNPs · ' +
                           result.n_diag_loose + ' loose / ' +
                           result.n_diag_strong + ' strong diagnostic';
      _drawHeatmap(canvas, result.cor_mat, result.bin_centers_mb);
      const v = result.summary.verdict;
      verdict.textContent =
        v + '  ·  within-A=' + _f2(result.summary.within_first_half_corr) +
        '  within-B=' + _f2(result.summary.within_second_half_corr) +
        '  between=' + _f2(result.summary.between_halves_corr);
      verdict.setAttribute('data-verdict', v);
    });

    return root;
  }

  function _f2(x) { return (x == null || isNaN(x)) ? 'NA' : x.toFixed(2); }

  function _defaultInputProvider() {
    if (typeof window === 'undefined' || !window.state) return null;
    const s = window.state;
    if (!s.candidate || !s.candidate.id) return null;
    const cand = s.candidate;
    // Need shelf interval — use candidate's bp range
    const start_bp = (typeof cand.start_bp === 'number') ? cand.start_bp
                    : (typeof cand.start_mb === 'number') ? Math.round(cand.start_mb * 1e6)
                    : null;
    const end_bp   = (typeof cand.end_bp === 'number') ? cand.end_bp
                    : (typeof cand.end_mb === 'number') ? Math.round(cand.end_mb * 1e6)
                    : null;
    if (start_bp == null || end_bp == null) return null;

    // Get the dosage chunk covering this region. The atlas cache exposes
    // chunks via window.popgenDosage if a helper is registered; else we
    // look at state.data.dosage_chunks.chunks (which is just an index of URLs
    // — so this works only when the chunk is already pre-fetched and cached
    // somewhere we can reach).
    let chunk = null;
    if (typeof window.popgenDosage !== 'undefined' &&
        typeof window.popgenDosage.getCachedChunk === 'function') {
      chunk = window.popgenDosage.getCachedChunk(start_bp, end_bp);
    }
    if (!chunk && s.data && s.data._dosage_chunk_inline) {
      chunk = s.data._dosage_chunk_inline;
    }
    if (!chunk) return null;

    // Build sample_invgt from cand.fish_calls
    const sample_invgt = new Map();
    if (Array.isArray(cand.fish_calls)) {
      const labels = (cand.k === 6 || cand._system === 'detailed')
        ? ['H1/H1', 'H1/H2', 'H2/H2', 'H1/H3', 'H2/H3', 'H3/H3']
        : ['Hom1', 'Het', 'Hom2'];
      for (const fc of cand.fish_calls) {
        if (!fc || fc.regime == null || fc.regime < 0) continue;
        const id = fc.sample_id || (Array.isArray(s.data.samples)
          ? s.data.samples[fc.sample_idx] : String(fc.sample_idx));
        const lbl = labels[fc.regime] || ('regime_' + fc.regime);
        sample_invgt.set(String(id), lbl);
      }
    }

    return { chunk, sample_invgt, shelf_start_bp: start_bp, shelf_end_bp: end_bp };
  }

  // Render a square correlation heatmap onto the canvas with diverging
  // red→white→blue scale (matches the R ggplot version).
  function _drawHeatmap(canvas, corMat, binCentersMb) {
    if (!canvas || !canvas.getContext) return;
    if (!Array.isArray(corMat) || corMat.length === 0) return;
    const ctx = canvas.getContext('2d');
    const N = corMat.length;
    const W = canvas.width, H = canvas.height;
    const padL = 60, padB = 30, padT = 10, padR = 10;
    const plotW = W - padL - padR, plotH = H - padT - padB;
    const cellW = plotW / N, cellH = plotH / N;

    ctx.clearRect(0, 0, W, H);

    // Heatmap cells
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const v = corMat[i][j];
        const x = padL + j * cellW;
        const y = padT + i * cellH;
        if (v == null || isNaN(v)) {
          ctx.fillStyle = '#e6e8ec';
        } else {
          ctx.fillStyle = _redBlueRamp(v);
        }
        ctx.fillRect(x, y, cellW + 0.5, cellH + 0.5);
      }
    }

    // Frame
    ctx.strokeStyle = '#555';
    ctx.lineWidth = 1;
    ctx.strokeRect(padL + 0.5, padT + 0.5, plotW, plotH);

    // Axis ticks at start/middle/end
    ctx.fillStyle = '#333';
    ctx.font = '10px ui-monospace, monospace';
    ctx.textAlign = 'center';
    if (binCentersMb && binCentersMb.length) {
      const ticks = [0, Math.floor(N / 2), N - 1];
      for (const t of ticks) {
        const cx = padL + (t + 0.5) * cellW;
        const lbl = binCentersMb[t].toFixed(1) + ' Mb';
        ctx.fillText(lbl, cx, padT + plotH + 14);
        ctx.save();
        ctx.translate(padL - 6, padT + (t + 0.5) * cellH);
        ctx.rotate(-Math.PI / 2);
        ctx.fillText(lbl, 0, 0);
        ctx.restore();
      }
    }
    // Color-scale label
    ctx.textAlign = 'left';
    ctx.fillText('r∈[-1,1]: red→white→blue', padL, padT - 2);
  }

  // Diverging red→white→blue ramp for r ∈ [-1, 1]
  function _redBlueRamp(r) {
    if (r == null || isNaN(r)) return '#e6e8ec';
    const t = Math.max(-1, Math.min(1, r));
    if (t >= 0) {
      // white (255,255,255) → blue (31,78,121)
      const f = t;
      const R = Math.round(255 + f * (31  - 255));
      const G = Math.round(255 + f * (78  - 255));
      const B = Math.round(255 + f * (121 - 255));
      return 'rgb(' + R + ',' + G + ',' + B + ')';
    } else {
      // white (255,255,255) → red (192,80,77)
      const f = -t;
      const R = Math.round(255 + f * (192 - 255));
      const G = Math.round(255 + f * (80  - 255));
      const B = Math.round(255 + f * (77  - 255));
      return 'rgb(' + R + ',' + G + ',' + B + ')';
    }
  }

  // ===========================================================================
  // Page 7 stratified Δ12 wiring — small, mostly mirrors turn 6
  // ===========================================================================

  // Returns a getData()-shaped resolver that produces multi-line stratified
  // Δ12 from the live ancestry response. Single-line cohort-mean Δ12 is the
  // fallback (state.data.ancestry_window.delta12).

  function makePage7AncestryStratifier() {
    return function (s) {
      if (!s) return null;
      // Try live first
      if (s.popstatsLive && s.popstatsLive.lastAncestryResponse) {
        const renderers = (typeof window !== 'undefined' && window.popgenRenderers)
          ? window.popgenRenderers : null;
        if (renderers && renderers.adaptAncestryResponse) {
          const adapted = renderers.adaptAncestryResponse(s.popstatsLive.lastAncestryResponse);
          // Flatten K → a multi-line track of "stratified Δ12" if the
          // response carries delta12 per group. Otherwise emit Q1 per group.
          if (adapted && adapted.ancestry_q1_per_group) {
            return adapted.ancestry_q1_per_group;
          }
        }
      }
      // Fallback: cohort-mean Δ12 from state.data.ancestry_window
      if (s.data && s.data.ancestry_window && Array.isArray(s.data.ancestry_window.delta12)) {
        const win = s.data.ancestry_window;
        const mb = Array.isArray(win.center_mb) ? win.center_mb
                  : (s.data.windows ? s.data.windows.map(w => w.center_mb)
                                    : []);
        return {
          mb,
          series: [{ name: 'cohort Δ12', color: '#2c7a39', values: win.delta12 }],
          values: win.delta12,
        };
      }
      return null;
    };
  }

  // Auto-mount on page 7. The atlas's existing renderAncestryPage looks for
  // the Δ12 track in its scalar stack; we patch the resolver to use ours
  // when popgenLive is available.
  function wirePage7DeltaStratified() {
    if (typeof window === 'undefined') return;
    if (window.popgenTurn7Page7Wired) return;
    const stratifier = makePage7AncestryStratifier();
    // Provide a hook for the atlas to look up before falling back to its
    // native getData()
    window.popgenPage7DeltaResolver = stratifier;
    // Trigger rerender on mode/Compute changes
    window.addEventListener('popgen:mode-changed', () => {
      if (typeof window.renderAncestryPage === 'function') {
        try { window.renderAncestryPage(); } catch (_) {}
      }
    });
    window.popgenTurn7Page7Wired = true;
  }

  // ===========================================================================
  // CSS injection
  // ===========================================================================

  const SHELF_CSS = `
    [data-popgen-shelf-ld="root"] {
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      border-radius: 4px;
      padding: 8px 10px;
      margin: 8px 0;
      font-family: var(--mono, ui-monospace, monospace);
      font-size: 11px; color: var(--ink, #0e1116);
    }
    [data-popgen-shelf-ld="header"] {
      display: flex; align-items: center; gap: 8px;
      padding-bottom: 4px;
      border-bottom: 1px solid var(--rule, #c8cdd2);
      margin-bottom: 6px;
    }
    [data-popgen-shelf-ld="title"] {
      flex: 1; font-weight: 600; font-size: 11px;
      letter-spacing: 0.04em; text-transform: uppercase;
      color: var(--ink-dim, #555e69);
    }
    [data-popgen-shelf-ld="run"] {
      padding: 3px 10px;
      border: 1px solid var(--accent, #f5a524);
      background: var(--accent, #f5a524);
      color: #0e1116; cursor: pointer;
      border-radius: 4px;
      font-family: inherit; font-size: 10px; font-weight: 600;
    }
    [data-popgen-shelf-ld="run"]:hover { background: #d99119; }
    [data-popgen-shelf-ld="status"] {
      font-size: 10px; color: var(--ink-dim, #555e69);
      padding: 2px 0;
    }
    [data-popgen-shelf-ld="verdict"] {
      padding: 4px 6px; margin-top: 4px;
      border-radius: 3px;
      font-size: 10px;
      background: var(--panel-2, #f0f1f3);
    }
    [data-popgen-shelf-ld="verdict"][data-verdict="SINGLE_INVERSION"] {
      background: rgba(60, 192, 138, 0.15); color: #1b5b3e;
      border: 1px solid #3cc08a;
    }
    [data-popgen-shelf-ld="verdict"][data-verdict="LIKELY_MULTIPLE_ARRANGEMENTS"] {
      background: rgba(224, 85, 92, 0.15); color: #6e1b1f;
      border: 1px solid #e0555c;
    }
    [data-popgen-shelf-ld="verdict"][data-verdict="AMBIGUOUS_OR_PARTIAL_LD"] {
      background: rgba(245, 165, 36, 0.15); color: #6d4707;
      border: 1px solid #f5a524;
    }
  `;

  function _injectCSS() {
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-shelf-ld-css')) return;
    const el = document.createElement('style');
    el.id = 'popgen-shelf-ld-css';
    el.textContent = SHELF_CSS;
    document.head.appendChild(el);
  }

  // ===========================================================================
  // Public exports
  // ===========================================================================

  return {
    // Q09b
    runShelfLDTest,
    pearsonCorrelationMatrix,
    verdictFromCorrMatrix,
    makeShelfLDPanel,

    // Page 7
    makePage7AncestryStratifier,
    wirePage7DeltaStratified,

    // Internals
    _internals: {
      _meanDosageOverCols, _redBlueRamp, _normalizeInvgt,
      _drawHeatmap, _defaultInputProvider,
    },
  };
}));
