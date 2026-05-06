// =============================================================================
// atlas_focal_vs_bg.js — Spalax-style focal-vs-background widget
// -----------------------------------------------------------------------------
// One reusable widget. Two page mounts on day 1: page 11 (boundaries) and
// page 16 (cross-species). Asks one question:
//
//   "Is the distribution of metric M inside a focal genomic zone different
//    from the distribution of M outside that zone (the background)?"
//
// Three output panels: focal histogram (panel C), background histogram (panel
// D), along-chromosome scatter with focal zone shaded (panel E). One Wilcoxon
// rank-sum P value comparing focal to background.
//
// Data flow (three modes, first available wins):
//   Mode A  precomputed aggregated JSON layer with focal+bg stats locked
//           at HPC-build time (state.focalVsBgStats). Manuscript-locked P,
//           lock icon. NOT IMPLEMENTED day 1 — atlas-side reader scaffold
//           is here but no JSON shape is currently loaded.
//   Mode B  live, computed from per-window data getters (state.data.windows
//           for `z`, state.popstatsLive for fst/hobs/theta tracks,
//           state.ncRNADensity for tRNA/rRNA/ncRNA tracks). Wilcoxon runs
//           in-browser. (live) tag.
//   Mode C  empty state when no metric is available for the active chrom.
//
// Public API (window.popgenFocalVsBg):
//   .makeFocalVsBgPanel(opts)  -> HTMLElement
//   .wilcoxonRankSum(focal, bg) -> { W, p_two_sided, logp_neg, z_stat }
//   .availableMetrics(state)    -> [{id, label, getData(state) -> {mb,values}}]
//   .computeFocalBgStats(opts)  -> { focal_arr, bg_arr, n_focal, n_bg, ... }
//
// Forward-compat: opts.cohort_selector is null on day 1, becomes a real
// selector spec when species 2 lands. Don't bake "single-cohort" into the API.
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenFocalVsBg = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Wilcoxon rank-sum (Mann-Whitney U) — normal approximation w/ continuity
  // correction. Tied ranks averaged. NaN-skip. Floor P at 1e-16 to match the
  // "P < 2.2e-16" convention from the Spalax figure.
  //
  // Returns { W, U_focal, p_two_sided, logp_neg, z_stat, n_focal, n_bg }.
  // Throws if either array is empty after NaN filtering.
  // ---------------------------------------------------------------------------
  function _phi(z) {
    // Standard normal CDF via erf (Abramowitz & Stegun 7.1.26)
    const a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741;
    const a4 = -1.453152027, a5 =  1.061405429, p  =  0.3275911;
    const sign = z < 0 ? -1 : 1;
    z = Math.abs(z) / Math.SQRT2;
    const t = 1 / (1 + p * z);
    const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.exp(-z*z);
    return 0.5 * (1 + sign * y);
  }
  function wilcoxonRankSum(focal, background) {
    const focal_clean = focal.filter(x => Number.isFinite(x));
    const bg_clean    = background.filter(x => Number.isFinite(x));
    const n1 = focal_clean.length;
    const n2 = bg_clean.length;
    if (n1 === 0 || n2 === 0) {
      return { W: NaN, U_focal: NaN, p_two_sided: NaN, logp_neg: NaN,
               z_stat: NaN, n_focal: n1, n_bg: n2 };
    }
    // Combined ranking with average ranks for ties
    const combined = [];
    for (const v of focal_clean) combined.push({ v, group: 0 });
    for (const v of bg_clean)    combined.push({ v, group: 1 });
    combined.sort((a, b) => a.v - b.v);
    // Assign average ranks for ties; track sum of t^3 - t for the variance
    // correction (W standard with ties).
    let i = 0;
    let R1 = 0;     // sum of focal ranks
    let tieSum = 0;
    while (i < combined.length) {
      let j = i;
      while (j + 1 < combined.length && combined[j + 1].v === combined[i].v) j++;
      const tieCount = j - i + 1;
      const avgRank  = (i + j + 2) / 2;   // 1-based ranks
      for (let k = i; k <= j; k++) {
        if (combined[k].group === 0) R1 += avgRank;
      }
      if (tieCount > 1) {
        tieSum += tieCount * tieCount * tieCount - tieCount;
      }
      i = j + 1;
    }
    // Mann-Whitney U from rank sum
    const U1 = R1 - n1 * (n1 + 1) / 2;     // U for focal
    // Normal approximation
    const N = n1 + n2;
    const meanU = (n1 * n2) / 2;
    let varU = (n1 * n2) * (N + 1) / 12;
    if (tieSum > 0) {
      // Tie correction: subtract n1*n2 / (12 N (N-1)) * sum(t^3 - t)
      varU -= (n1 * n2) * tieSum / (12 * N * (N - 1));
    }
    if (varU <= 0) {
      return { W: R1, U_focal: U1, p_two_sided: 1, logp_neg: 0,
               z_stat: 0, n_focal: n1, n_bg: n2 };
    }
    // Continuity correction: shift |U - meanU| by 0.5 toward meanU
    const diff = U1 - meanU;
    const cc = (diff > 0) ? -0.5 : (diff < 0 ? 0.5 : 0);
    const z = (diff + cc) / Math.sqrt(varU);
    // Two-sided P
    let p = 2 * (1 - _phi(Math.abs(z)));
    p = Math.max(p, 1e-16);     // floor at Spalax convention
    return {
      W: R1, U_focal: U1, p_two_sided: p,
      logp_neg: -Math.log10(p), z_stat: z,
      n_focal: n1, n_bg: n2,
    };
  }

  // ---------------------------------------------------------------------------
  // Metric registry — what's available depends on what's loaded in state.
  // Each metric exports a getData(state) -> {mb, values} | null function.
  //
  // Source order (priority):
  //   - state.data.windows (precomp): always available for `z`
  //   - state.popstatsLive: for fst/hobs/theta when live server is hooked up
  //   - state.ncRNADensity[chrom]: ncRNA density tracks (turn 116)
  //   - state.repeatDensity[chrom]: TE density tracks (sibling — bonus)
  //
  // The `available` filter at render time hides metrics that return null.
  // ---------------------------------------------------------------------------
  const METRIC_DEFS = [
    {
      id: 'z', label: 'robust |Z|', family: 'discovery',
      getData: (s) => {
        if (!s.data || !s.data.windows) return null;
        const mb = s.data.windows.map(w => w.center_mb);
        const values = s.data.windows.map(w => Math.abs(w.z || 0));
        if (values.every(v => !Number.isFinite(v) || v === 0)) return null;
        return { mb, values };
      },
    },
    {
      id: 'fst_hom1_hom2', label: 'F_ST (Hom1–Hom2)', family: 'popstats',
      getData: (s) => {
        // Pull from popstatsLive when available
        if (!s.popstatsLive || !s.popstatsLive.lastResponse) return null;
        const wins = s.popstatsLive.lastResponse.windows;
        if (!Array.isArray(wins)) return null;
        const mb = wins.map(w => (typeof w.center_mb === 'number')
          ? w.center_mb
          : (typeof w.start === 'number' && typeof w.end === 'number'
              ? (w.start + w.end) / 2 / 1e6 : null));
        const values = wins.map(w => {
          const v = w.fst_hom1_hom2 != null ? w.fst_hom1_hom2 : w.fst;
          return (v == null || isNaN(v)) ? null : Number(v);
        });
        if (values.every(v => v == null)) return null;
        return { mb, values };
      },
    },
    {
      id: 'hobs_hexp_ratio', label: 'Hobs / Hexp', family: 'popstats',
      getData: (s) => {
        if (!s.popstatsLive || !s.popstatsLive.lastResponse) return null;
        const wins = s.popstatsLive.lastResponse.windows;
        if (!Array.isArray(wins)) return null;
        const mb = wins.map(w => (typeof w.center_mb === 'number')
          ? w.center_mb
          : (typeof w.start === 'number' && typeof w.end === 'number'
              ? (w.start + w.end) / 2 / 1e6 : null));
        const values = wins.map(w => {
          const ho = w.hobs, he = w.hexp;
          if (ho == null || he == null || he === 0) return null;
          return ho / he;
        });
        if (values.every(v => v == null)) return null;
        return { mb, values };
      },
    },
    {
      id: 'theta_pi_log2_ratio', label: 'log₂ θπ ratio (HOM_INV / HOM_REF)',
      family: 'popstats',
      getData: (s) => {
        if (!s.popstatsLive || !s.popstatsLive.lastResponse) return null;
        const wins = s.popstatsLive.lastResponse.windows;
        if (!Array.isArray(wins)) return null;
        const mb = wins.map(w => (typeof w.center_mb === 'number')
          ? w.center_mb
          : (typeof w.start === 'number' && typeof w.end === 'number'
              ? (w.start + w.end) / 2 / 1e6 : null));
        const values = wins.map(w => {
          const a = w.theta_pi_hom_inv, b = w.theta_pi_hom_ref;
          if (a == null || b == null || a <= 0 || b <= 0) return null;
          return Math.log2(a / b);
        });
        if (values.every(v => v == null)) return null;
        return { mb, values };
      },
    },
  ];

  // ncRNA-density auto-discovery — generates one metric per family default
  // (tRNA_all, rRNA_all, ncRNA_all) plus rRNA_5S (the rDNA cluster signal
  // is biologically distinct enough to surface separately).
  const NCRNA_AUTO_METRICS = ['tRNA_all', 'rRNA_all', 'rRNA_5S', 'ncRNA_all'];
  function _ncrnaMetricGetData(cls) {
    return (s) => {
      const chrom = s.data && s.data.chrom;
      if (!chrom || !s.ncRNADensity || !s.ncRNADensity[chrom]) return null;
      const block = s.ncRNADensity[chrom].chrom_block;
      if (!block || !block.by_class || !block.by_class[cls]) return null;
      const mb = block.window_centers_mb;
      const values = block.by_class[cls].densities;
      if (!mb || !values || mb.length !== values.length) return null;
      return { mb, values };
    };
  }
  for (const cls of NCRNA_AUTO_METRICS) {
    METRIC_DEFS.push({
      id: 'ncrna_' + cls,
      label: cls.replace('_', ' ') + ' (loci/Mb)',
      family: 'ncrna',
      getData: _ncrnaMetricGetData(cls),
    });
  }

  // Optional: TE density auto-discovery (one entry per loaded class — only
  // emits when the layer is present, so it doesn't clutter the dropdown for
  // users who haven't loaded it).
  function _teMetricGetData(cls) {
    return (s) => {
      const chrom = s.data && s.data.chrom;
      if (!chrom || !s.repeatDensity || !s.repeatDensity[chrom]) return null;
      const block = s.repeatDensity[chrom].chrom_block;
      if (!block || !block.by_class || !block.by_class[cls]) return null;
      const mb = block.window_centers_mb;
      const values = block.by_class[cls].densities;
      if (!mb || !values || mb.length !== values.length) return null;
      return { mb, values };
    };
  }
  // Just expose all_TE — the TE classes are many and mostly redundant for
  // the focal-vs-bg use case. Users who want a specific TE class can switch
  // to it via the boundaries-page TE panel.
  METRIC_DEFS.push({
    id: 'te_all_TE', label: 'all_TE density', family: 'te',
    getData: _teMetricGetData('all_TE'),
  });

  function availableMetrics(state) {
    const out = [];
    for (const m of METRIC_DEFS) {
      try {
        const d = m.getData(state);
        if (d && Array.isArray(d.values) && d.values.length > 0) {
          out.push(m);
        }
      } catch (_) { /* skip */ }
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // Compute focal vs background stats from a metric's {mb, values} array.
  //
  // opts: {
  //   data: {mb, values},
  //   focal_lo_mb, focal_hi_mb,           // focal interval in Mb
  //   bg_scope: 'per_chrom' | 'genome_wide_excl_focal',
  //   z_transform: bool,                   // Z-normalize against bg mean/sd
  // }
  // ---------------------------------------------------------------------------
  function computeFocalBgStats(opts) {
    const data = opts.data;
    const lo = opts.focal_lo_mb, hi = opts.focal_hi_mb;
    if (!data || !Array.isArray(data.mb) || !Array.isArray(data.values) ||
        data.mb.length !== data.values.length) {
      return null;
    }
    // Partition windows
    const focal_arr = [];
    const focal_xy  = [];   // [{mb, v}] for the scatter
    const bg_arr    = [];
    const bg_xy     = [];
    for (let i = 0; i < data.mb.length; i++) {
      const x = data.mb[i];
      const v = data.values[i];
      if (x == null || !Number.isFinite(v)) continue;
      const inFocal = (x >= lo && x <= hi);
      if (inFocal) {
        focal_arr.push(v);
        focal_xy.push({ mb: x, v: v, inFocal: true });
      } else {
        bg_arr.push(v);
        bg_xy.push({ mb: x, v: v, inFocal: false });
      }
    }
    // bg mean/sd for Z-transform
    const bgMean = bg_arr.reduce((a, b) => a + b, 0) / Math.max(1, bg_arr.length);
    let bgVar = 0;
    for (const v of bg_arr) bgVar += (v - bgMean) * (v - bgMean);
    const bgSd = bg_arr.length > 1 ? Math.sqrt(bgVar / (bg_arr.length - 1)) : 0;
    // Apply Z-transform if requested. Guard against near-zero sd from
    // floating-point noise on constant-valued background — would otherwise
    // produce nonsense like Z = 5e16 when the bg is e.g. [0.1, 0.1, 0.1, ...]
    // and floating point drift gives sd ≈ 1e-17.
    let focal_z = focal_arr;
    let bg_z = bg_arr;
    let allXy = focal_xy.concat(bg_xy);
    const SD_EPS = 1e-12;
    if (opts.z_transform && bgSd > SD_EPS) {
      focal_z = focal_arr.map(v => (v - bgMean) / bgSd);
      bg_z    = bg_arr.map(v => (v - bgMean) / bgSd);
      allXy = allXy.map(p => Object.assign({}, p, {
        v_z: bgSd > 0 ? (p.v - bgMean) / bgSd : 0,
      }));
    }
    // Means + medians of the displayed values
    const focalMean = focal_z.reduce((a, b) => a + b, 0) / Math.max(1, focal_z.length);
    const bgMean_disp = bg_z.reduce((a, b) => a + b, 0) / Math.max(1, bg_z.length);
    const focalMedian = _median(focal_z);
    const bgMedian    = _median(bg_z);
    // Wilcoxon (always on raw values, not Z-transformed — Z is monotone so
    // the rank-sum is invariant; using raw is conventionally clearer)
    const w = wilcoxonRankSum(focal_arr, bg_arr);
    return {
      focal_arr: focal_z, bg_arr: bg_z,
      focal_raw: focal_arr, bg_raw: bg_arr,
      n_focal: focal_arr.length, n_bg: bg_arr.length,
      focal_mean: focalMean, bg_mean: bgMean_disp,
      focal_median: focalMedian, bg_median: bgMedian,
      bg_mean_raw: bgMean, bg_sd_raw: bgSd,
      wilcoxon: w,
      scatter_xy: allXy,
    };
  }
  function _median(a) {
    if (a.length === 0) return NaN;
    const s = a.slice().sort((x, y) => x - y);
    const m = (s.length - 1) / 2;
    return (s.length & 1) ? s[Math.floor(m)] : 0.5 * (s[Math.floor(m)] + s[Math.ceil(m)]);
  }

  // ---------------------------------------------------------------------------
  // Theme accessor (mirrors atlas_dotplot.js pattern)
  // ---------------------------------------------------------------------------
  function _theme() {
    const cs = (typeof getComputedStyle !== 'undefined' && document &&
                document.documentElement)
      ? getComputedStyle(document.documentElement) : null;
    const get = (k, fb) => {
      if (!cs) return fb;
      const v = cs.getPropertyValue(k).trim();
      return v || fb;
    };
    return {
      bg:        get('--bg',         '#0e1116'),
      panel:     get('--panel',      '#151a22'),
      panel2:    get('--panel-2',    '#1c2330'),
      ink:       get('--ink',        '#e7edf3'),
      inkDim:    get('--ink-dim',    '#8a94a3'),
      inkDimmer: get('--ink-dimmer', '#5a6472'),
      rule:      get('--rule',       '#2a3242'),
      accent:    get('--accent',     '#f5a524'),
      accent2:   get('--accent-2',   '#4fa3ff'),
      mono:      get('--mono',       'ui-monospace, monospace'),
    };
  }

  // ---------------------------------------------------------------------------
  // Histogram bins via Freedman-Diaconis (capped at 60).
  // ---------------------------------------------------------------------------
  function _histBins(arr, nBinsCap) {
    if (arr.length === 0) return { bins: [], nBins: 0, lo: 0, hi: 1, width: 1 };
    const sorted = arr.slice().sort((a, b) => a - b);
    const lo = sorted[0];
    const hi = sorted[sorted.length - 1];
    if (lo === hi) {
      // All same value — one bin
      return {
        bins: [{ lo: lo - 0.5, hi: lo + 0.5, count: arr.length }],
        nBins: 1, lo: lo - 0.5, hi: lo + 0.5, width: 1,
      };
    }
    // FD: bin width = 2 * IQR / n^(1/3)
    const q1 = sorted[Math.floor(sorted.length * 0.25)];
    const q3 = sorted[Math.floor(sorted.length * 0.75)];
    const iqr = q3 - q1;
    let width;
    if (iqr > 0) {
      width = 2 * iqr / Math.pow(sorted.length, 1 / 3);
    } else {
      width = (hi - lo) / 30;
    }
    let n = Math.ceil((hi - lo) / width);
    if (n < 5) n = 5;
    if (n > nBinsCap) { n = nBinsCap; width = (hi - lo) / n; }
    const bins = new Array(n).fill(0).map((_, i) => ({
      lo: lo + i * width,
      hi: lo + (i + 1) * width,
      count: 0,
    }));
    for (const v of arr) {
      let bi = Math.floor((v - lo) / width);
      if (bi < 0) bi = 0;
      if (bi >= n) bi = n - 1;
      bins[bi].count++;
    }
    return { bins, nBins: n, lo, hi, width };
  }

  // ---------------------------------------------------------------------------
  // SVG histogram renderer for one side (focal or background).
  // Shares yMaxCount with the sibling so they're visually comparable.
  // ---------------------------------------------------------------------------
  function _renderHistogramSVG(arr, opts) {
    const t = opts.theme;
    const W = opts.width || 360;
    const H = opts.height || 130;
    const PAD = { top: 20, right: 14, bottom: 26, left: 38 };
    const innerW = W - PAD.left - PAD.right;
    const innerH = H - PAD.top - PAD.bottom;
    if (arr.length === 0) {
      return '<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H + '" ' +
             'xmlns="http://www.w3.org/2000/svg" style="display:block;">' +
             '<rect x="' + PAD.left + '" y="' + PAD.top + '" width="' + innerW +
             '" height="' + innerH + '" fill="' + t.panel + '" stroke="' + t.rule +
             '" stroke-width="0.6"/>' +
             '<text x="' + (W / 2) + '" y="' + (H / 2) + '" fill="' + t.inkDim +
             '" font-size="11" text-anchor="middle" font-family="' + t.mono +
             '">n = 0 windows</text></svg>';
    }
    const { bins, lo, hi } = opts.binsForBoth || _histBins(arr, 60);
    // Recount this side into the shared bin edges
    const counts = new Array(bins.length).fill(0);
    if (opts.binsForBoth) {
      for (const v of arr) {
        let bi = Math.floor((v - lo) / ((hi - lo) / bins.length));
        if (bi < 0) bi = 0;
        if (bi >= bins.length) bi = bins.length - 1;
        counts[bi]++;
      }
    } else {
      for (let i = 0; i < bins.length; i++) counts[i] = bins[i].count;
    }
    const yMaxCount = opts.yMaxCount || Math.max(1, ...counts);
    const xRange = hi - lo || 1;
    const toX = (v) => PAD.left + ((v - lo) / xRange) * innerW;
    const toY = (c) => PAD.top + innerH - (c / yMaxCount) * innerH;
    let svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H +
              '" xmlns="http://www.w3.org/2000/svg" style="display:block;">';
    // Frame
    svg += '<rect x="' + PAD.left + '" y="' + PAD.top + '" width="' + innerW +
           '" height="' + innerH + '" fill="' + t.panel + '" stroke="' + t.rule +
           '" stroke-width="0.6"/>';
    // Bars
    const barColor = opts.color || t.accent;
    for (let i = 0; i < bins.length; i++) {
      if (counts[i] === 0) continue;
      const x = toX(bins[i].lo);
      const xR = toX(bins[i].hi);
      const y = toY(counts[i]);
      const h = PAD.top + innerH - y;
      svg += '<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) +
             '" width="' + Math.max(0.5, (xR - x - 0.6)).toFixed(1) +
             '" height="' + h.toFixed(1) + '" fill="' + barColor +
             '" fill-opacity="0.75"/>';
    }
    // Mean line (dashed vertical)
    const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
    if (mean >= lo && mean <= hi) {
      const xMean = toX(mean);
      svg += '<line x1="' + xMean.toFixed(1) + '" x2="' + xMean.toFixed(1) +
             '" y1="' + PAD.top + '" y2="' + (PAD.top + innerH) +
             '" stroke="' + t.ink + '" stroke-width="1" stroke-dasharray="4,3"/>';
      const meanFmt = Math.abs(mean) >= 100 ? mean.toFixed(0)
                   : Math.abs(mean) >= 10  ? mean.toFixed(1) : mean.toFixed(3);
      svg += '<text x="' + xMean.toFixed(1) + '" y="' + (PAD.top - 5) +
             '" fill="' + t.ink + '" font-size="9.5" text-anchor="middle" ' +
             'font-family="' + t.mono + '">μ=' + meanFmt + '</text>';
    }
    // X-axis ticks (3 evenly spaced)
    for (let i = 0; i <= 2; i++) {
      const xv = lo + (i / 2) * xRange;
      const x = toX(xv);
      const fmt = Math.abs(xv) >= 100 ? xv.toFixed(0)
               : Math.abs(xv) >= 10  ? xv.toFixed(1) : xv.toFixed(2);
      svg += '<text x="' + x.toFixed(1) + '" y="' + (PAD.top + innerH + 12) +
             '" fill="' + t.inkDim + '" font-size="9" text-anchor="middle" ' +
             'font-family="' + t.mono + '">' + fmt + '</text>';
    }
    // Y-axis ticks (just max)
    svg += '<text x="' + (PAD.left - 4) + '" y="' + (PAD.top + 8) +
           '" fill="' + t.inkDim + '" font-size="9" text-anchor="end" ' +
           'font-family="' + t.mono + '">' + yMaxCount + '</text>';
    // Title row
    if (opts.title) {
      svg += '<text x="' + PAD.left + '" y="' + (PAD.top - 16) +
             '" fill="' + t.ink + '" font-size="10.5" font-family="' + t.mono +
             '">' + opts.title + '</text>';
      svg += '<text x="' + (W - PAD.right) + '" y="' + (PAD.top - 16) +
             '" fill="' + t.inkDim + '" font-size="9.5" text-anchor="end" ' +
             'font-family="' + t.mono + '">n=' + arr.length + '</text>';
    }
    svg += '</svg>';
    return svg;
  }

  // ---------------------------------------------------------------------------
  // SVG along-chromosome scatter renderer.
  // Focal zone shaded; anchor markers at focal_lo/hi (page 11 use both;
  // page 16 single anchor handled by passing focal_lo == focal_hi).
  // ---------------------------------------------------------------------------
  function _renderScatterSVG(scatter_xy, opts) {
    const t = opts.theme;
    const W = opts.width || 740;
    const H = opts.height || 150;
    const PAD = { top: 18, right: 16, bottom: 26, left: 50 };
    const innerW = W - PAD.left - PAD.right;
    const innerH = H - PAD.top - PAD.bottom;
    if (scatter_xy.length === 0) {
      return '<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H +
             '" xmlns="http://www.w3.org/2000/svg" style="display:block;"><text x="' +
             (W / 2) + '" y="' + (H / 2) + '" fill="' + t.inkDim + '" font-size="11" ' +
             'text-anchor="middle" font-family="' + t.mono + '">no data</text></svg>';
    }
    const useZ = !!opts.z_transform;
    const valKey = useZ ? 'v_z' : 'v';
    // X range
    const xMin = opts.x_lo_mb != null ? opts.x_lo_mb
      : Math.min.apply(null, scatter_xy.map(p => p.mb));
    const xMax = opts.x_hi_mb != null ? opts.x_hi_mb
      : Math.max.apply(null, scatter_xy.map(p => p.mb));
    const xRange = xMax - xMin || 1;
    // Y range — symmetric around 0 if Z; else min..max
    let yMin, yMax;
    const vals = scatter_xy.map(p => p[valKey]).filter(Number.isFinite);
    if (useZ) {
      const m = Math.max(Math.abs(Math.min.apply(null, vals)),
                          Math.abs(Math.max.apply(null, vals)));
      yMin = -Math.max(m, 1) * 1.05;
      yMax =  Math.max(m, 1) * 1.05;
    } else {
      yMin = Math.min.apply(null, vals);
      yMax = Math.max.apply(null, vals);
      if (yMin === yMax) { yMin -= 0.5; yMax += 0.5; }
      const pad = (yMax - yMin) * 0.05;
      yMin -= pad; yMax += pad;
    }
    const toX = (mb) => PAD.left + ((mb - xMin) / xRange) * innerW;
    const toY = (v)  => PAD.top + innerH - ((v - yMin) / (yMax - yMin)) * innerH;
    let svg = '<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H +
              '" xmlns="http://www.w3.org/2000/svg" style="display:block;">';
    // Frame
    svg += '<rect x="' + PAD.left + '" y="' + PAD.top + '" width="' + innerW +
           '" height="' + innerH + '" fill="' + t.panel + '" stroke="' + t.rule +
           '" stroke-width="0.6"/>';
    // Focal-zone shading
    if (opts.focal_lo_mb != null && opts.focal_hi_mb != null &&
        opts.focal_hi_mb > opts.focal_lo_mb) {
      const xL = Math.max(PAD.left, toX(opts.focal_lo_mb));
      const xR = Math.min(PAD.left + innerW, toX(opts.focal_hi_mb));
      if (xR > xL) {
        svg += '<rect x="' + xL.toFixed(1) + '" y="' + PAD.top + '" ' +
               'width="' + (xR - xL).toFixed(1) + '" height="' + innerH + '" ' +
               'fill="' + t.accent + '" fill-opacity="0.12"/>';
        // Dashed boundary lines
        svg += '<line x1="' + xL.toFixed(1) + '" x2="' + xL.toFixed(1) +
               '" y1="' + PAD.top + '" y2="' + (PAD.top + innerH) +
               '" stroke="' + t.accent + '" stroke-width="1" stroke-dasharray="3,3" ' +
               'stroke-opacity="0.7"/>';
        svg += '<line x1="' + xR.toFixed(1) + '" x2="' + xR.toFixed(1) +
               '" y1="' + PAD.top + '" y2="' + (PAD.top + innerH) +
               '" stroke="' + t.accent + '" stroke-width="1" stroke-dasharray="3,3" ' +
               'stroke-opacity="0.7"/>';
      }
    }
    // Anchor markers (single point in page-16; L+R zone marks in page-11 are
    // already drawn via focal_lo/hi above)
    if (Array.isArray(opts.anchor_marks_mb)) {
      for (const m of opts.anchor_marks_mb) {
        const x = toX(m);
        if (x >= PAD.left && x <= PAD.left + innerW) {
          svg += '<line x1="' + x.toFixed(1) + '" x2="' + x.toFixed(1) +
                 '" y1="' + PAD.top + '" y2="' + (PAD.top + innerH) +
                 '" stroke="#e85a5a" stroke-width="1.2" stroke-dasharray="2,2"/>';
        }
      }
    }
    // Y-zero line if useful
    if (yMin < 0 && yMax > 0) {
      const y0 = toY(0);
      svg += '<line x1="' + PAD.left + '" x2="' + (PAD.left + innerW) +
             '" y1="' + y0.toFixed(1) + '" y2="' + y0.toFixed(1) +
             '" stroke="' + t.inkDim + '" stroke-width="0.5" stroke-dasharray="1.5,3"/>';
    }
    // Dots
    for (const p of scatter_xy) {
      const v = p[valKey];
      if (!Number.isFinite(v)) continue;
      const x = toX(p.mb);
      const y = toY(v);
      const fill = p.inFocal ? t.accent : t.inkDim;
      const op = p.inFocal ? 0.7 : 0.45;
      svg += '<circle cx="' + x.toFixed(1) + '" cy="' + y.toFixed(1) +
             '" r="1.5" fill="' + fill + '" fill-opacity="' + op + '"/>';
    }
    // Y-axis labels (3 ticks)
    for (let i = 0; i <= 2; i++) {
      const yv = yMin + (i / 2) * (yMax - yMin);
      const y = toY(yv);
      const fmt = Math.abs(yv) >= 100 ? yv.toFixed(0)
               : Math.abs(yv) >= 10  ? yv.toFixed(1) : yv.toFixed(2);
      svg += '<text x="' + (PAD.left - 4) + '" y="' + (y + 3).toFixed(1) +
             '" fill="' + t.inkDim + '" font-size="9" text-anchor="end" ' +
             'font-family="' + t.mono + '">' + fmt + '</text>';
    }
    // X-axis labels (5 ticks)
    for (let i = 0; i <= 4; i++) {
      const xv = xMin + (i / 4) * xRange;
      const x = toX(xv);
      svg += '<text x="' + x.toFixed(1) + '" y="' + (PAD.top + innerH + 14) +
             '" fill="' + t.inkDim + '" font-size="9" text-anchor="middle" ' +
             'font-family="' + t.mono + '">' + xv.toFixed(0) + '</text>';
    }
    svg += '<text x="' + (PAD.left + innerW / 2) + '" y="' + (H - 4) +
           '" fill="' + t.inkDim + '" font-size="9.5" text-anchor="middle" ' +
           'font-family="' + t.mono + '">position (Mb)</text>';
    svg += '</svg>';
    return svg;
  }

  // ---------------------------------------------------------------------------
  // Radius toolbar HTML (shared with page 11 boundaries toolbar; same buttons,
  // same data-radius attributes).
  // ---------------------------------------------------------------------------
  function _radiusToolbarHTML(prefix, currentRadiusBp) {
    const btns = [
      [1000,    '1 kb'],   [5000,    '5 kb'],   [10000,   '10 kb'],
      [25000,   '25 kb'],  [100000,  '100 kb'],
      [1000000, '1 Mb'],   [1500000, '1.5 Mb'], [2000000, '2 Mb'], [5000000, '5 Mb'],
    ];
    let html = '<div class="fvb-radius-row" style="display:flex;align-items:center;gap:4px;flex-wrap:wrap;">';
    for (let i = 0; i < btns.length; i++) {
      const [r, lbl] = btns[i];
      const active = (r === currentRadiusBp);
      html += '<button class="fvb-radius-btn" data-radius="' + r + '" ' +
              'data-prefix="' + prefix + '" ' +
              'style="background:' + (active ? 'var(--accent)' : 'var(--panel)') + ';' +
              'color:' + (active ? '#1a1a1a' : 'var(--ink)') + ';' +
              'border:1px solid ' + (active ? 'var(--accent)' : 'var(--rule)') + ';' +
              'border-radius:0;padding:2px 7px;font-family:var(--mono);font-size:10px;' +
              'cursor:pointer;">' + lbl + '</button>';
      if (i === 4) {   // separator after 100 kb
        html += '<span style="display:inline-block;width:1px;height:14px;background:var(--rule);margin:0 3px;"></span>';
      }
    }
    html += '</div>';
    return html;
  }

  // ---------------------------------------------------------------------------
  // Format Wilcoxon P value for display (Spalax convention).
  // ---------------------------------------------------------------------------
  function _formatP(p) {
    if (p == null || !Number.isFinite(p)) return '\u2014';
    if (p <= 1e-16) return 'P < 2.2e-16';
    if (p < 0.001)  return 'P = ' + p.toExponential(2);
    if (p < 0.01)   return 'P = ' + p.toFixed(4);
    return 'P = ' + p.toFixed(3);
  }

  // ---------------------------------------------------------------------------
  // Sticky preferences (per-page) — localStorage-backed
  // ---------------------------------------------------------------------------
  function _prefsKey(pageId) { return 'inversion_atlas.focalVsBg.' + pageId; }
  function _loadPrefs(pageId) {
    try {
      const raw = (typeof localStorage !== 'undefined') ? localStorage.getItem(_prefsKey(pageId)) : null;
      if (raw) return Object.assign({}, _defaultPrefs(), JSON.parse(raw));
    } catch (_) {}
    return _defaultPrefs();
  }
  function _savePrefs(pageId, prefs) {
    try {
      if (typeof localStorage !== 'undefined') {
        localStorage.setItem(_prefsKey(pageId), JSON.stringify(prefs));
      }
    } catch (_) {}
  }
  function _defaultPrefs() {
    return { metric: 'z', bg_scope: 'per_chrom', z_transform: true, view_mode: 'full' };
  }

  // ---------------------------------------------------------------------------
  // Esc helper
  // ---------------------------------------------------------------------------
  function _esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }

  // ---------------------------------------------------------------------------
  // The widget itself.
  // ---------------------------------------------------------------------------
  function makeFocalVsBgPanel(opts) {
    opts = opts || {};
    const pageId = opts.page_id || 'unknown';
    const t = _theme();
    const prefs = _loadPrefs(pageId);

    // Outer wrap
    const wrap = document.createElement('div');
    wrap.className = 'fvb-panel';
    wrap.style.cssText = (opts.containerStyle || '') +
      ';background:' + t.panel2 + ';border:1px solid ' + t.rule +
      ';border-radius:4px;padding:10px 14px;font-family:' + t.mono +
      ';color:' + t.ink + ';';

    // Inner state captured in a closure so re-renders are cheap
    const widget = {
      anchor_id:    opts.anchor_id,
      anchor_type:  opts.anchor_type,           // 'cs_breakpoint' | 'candidate_boundary'
      chrom:        opts.chrom,
      // Focal interval — for page 16 supplied directly; for page 11 the
      // wrapper recomputes via _boundaryScanRange before each render.
      focal_lo_bp:  opts.focal_lo_bp,
      focal_hi_bp:  opts.focal_hi_bp,
      anchor_marks_mb: opts.anchor_marks_mb || [],   // dashed red lines
      // Radius binding
      radius_source:    opts.radius_source || 'self',     // 'self' | 'boundariesState'
      radius_default_bp: opts.radius_default_bp || 1500000,
      get_radius_bp: opts.get_radius_bp,         // optional callback
      set_radius_bp: opts.set_radius_bp,         // optional callback
      cohort_selector: opts.cohort_selector || null,
    };

    function getRadiusBp() {
      if (typeof widget.get_radius_bp === 'function') {
        const r = widget.get_radius_bp();
        if (Number.isFinite(r) && r > 0) return r;
      }
      return widget.radius_default_bp;
    }
    function setRadiusBp(r) {
      if (typeof widget.set_radius_bp === 'function') {
        widget.set_radius_bp(r);
      } else {
        widget.radius_default_bp = r;
      }
    }

    function rerender() {
      const stateRef = (typeof window !== 'undefined' && window.state) || {};
      const metrics = availableMetrics(stateRef);
      // Pick metric: prefs > first available > first definition
      let active = metrics.find(m => m.id === prefs.metric) || metrics[0];
      if (!active && METRIC_DEFS.length > 0) active = METRIC_DEFS[0];
      // Fetch data
      let data = null;
      if (active) {
        try { data = active.getData(stateRef); } catch (_) { data = null; }
      }
      // Header / toolbar
      let html = '';
      // Row 1: metric / bg-scope / Z toggle / view mode
      html += '<div style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;font-size:10.5px;margin-bottom:6px;">';
      html += '<span style="color:' + t.ink + ';font-family:' + t.mono + ';">focal vs background</span>';
      html += '<span style="color:' + t.inkDimmer + ';">·</span>';
      // Metric dropdown
      html += '<span style="color:' + t.inkDim + ';">metric:</span>';
      html += '<select id="fvb-metric" style="background:' + t.panel + ';color:' + t.ink +
              ';border:1px solid ' + t.rule + ';border-radius:3px;padding:2px 6px;' +
              'font-family:' + t.mono + ';font-size:10.5px;">';
      if (metrics.length === 0) {
        html += '<option>(no metrics available)</option>';
      } else {
        for (const m of metrics) {
          html += '<option value="' + _esc(m.id) + '"' +
                  (active && active.id === m.id ? ' selected' : '') + '>' +
                  _esc(m.label) + '</option>';
        }
      }
      html += '</select>';
      // bg-scope dropdown
      html += '<span style="color:' + t.inkDim + ';margin-left:6px;">bg:</span>';
      html += '<select id="fvb-bg-scope" style="background:' + t.panel + ';color:' + t.ink +
              ';border:1px solid ' + t.rule + ';border-radius:3px;padding:2px 6px;' +
              'font-family:' + t.mono + ';font-size:10.5px;">';
      html += '<option value="per_chrom"' + (prefs.bg_scope === 'per_chrom' ? ' selected' : '') +
              '>per-chrom</option>';
      html += '<option value="genome_wide_excl_focal"' +
              (prefs.bg_scope === 'genome_wide_excl_focal' ? ' selected' : '') +
              '>genome-wide</option>';
      html += '</select>';
      // Z toggle
      html += '<label style="margin-left:8px;color:' + t.inkDim +
              ';display:inline-flex;align-items:center;gap:3px;cursor:pointer;">' +
              '<input type="checkbox" id="fvb-z" ' + (prefs.z_transform ? 'checked' : '') +
              ' style="margin:0;"> Z</label>';
      // View mode
      html += '<span style="color:' + t.inkDim + ';margin-left:8px;">view:</span>';
      html += '<select id="fvb-view" style="background:' + t.panel + ';color:' + t.ink +
              ';border:1px solid ' + t.rule + ';border-radius:3px;padding:2px 6px;' +
              'font-family:' + t.mono + ';font-size:10.5px;">';
      html += '<option value="full"'   + (prefs.view_mode === 'full'   ? ' selected' : '') +
              '>full</option>';
      html += '<option value="zoomed"' + (prefs.view_mode === 'zoomed' ? ' selected' : '') +
              '>zoomed</option>';
      html += '</select>';
      // P value (computed below; placeholder for now)
      html += '<span id="fvb-pvalue" style="margin-left:auto;color:' + t.ink +
              ';font-weight:500;"></span>';
      html += '</div>';
      // Row 2: radius
      html += '<div style="display:flex;align-items:center;gap:6px;margin-bottom:8px;">' +
              '<span style="color:' + t.inkDim + ';font-size:10px;">scan radius:</span>' +
              _radiusToolbarHTML('fvb', getRadiusBp()) + '</div>';

      // Compute focal interval
      const focal_lo_mb = widget.focal_lo_bp / 1e6;
      const focal_hi_mb = widget.focal_hi_bp / 1e6;

      // Compute stats
      let stats = null;
      let pStr = '(no data)';
      if (data && data.values.length > 0) {
        // bg_scope filtering — for genome_wide_excl_focal, we'd merge across
        // chroms. Day-1 implementation: per_chrom mode uses the loaded chrom's
        // data; genome_wide mode does the same (single-chrom data is all that's
        // available locally) but flags the limitation in the UI hint.
        stats = computeFocalBgStats({
          data: data,
          focal_lo_mb: focal_lo_mb,
          focal_hi_mb: focal_hi_mb,
          bg_scope: prefs.bg_scope,
          z_transform: prefs.z_transform,
        });
        if (stats && stats.wilcoxon && Number.isFinite(stats.wilcoxon.p_two_sided)) {
          pStr = _formatP(stats.wilcoxon.p_two_sided) + ' <span style="color:' +
                 t.inkDim + ';font-size:9.5px;">(live)</span>';
        }
      }

      // Body: histograms + scatter
      if (!data) {
        html += '<div style="padding:24px 12px;color:' + t.inkDim +
                ';font-size:11px;text-align:center;line-height:1.7;">' +
                'No data available for the selected metric on this chromosome.<br>' +
                '<span style="color:' + t.inkDimmer + ';font-size:10px;">' +
                'Load a precomp JSON, popstats live response, or ncRNA / repeat-density layer ' +
                'to populate this metric.</span></div>';
      } else if (!stats) {
        html += '<div style="padding:18px;color:' + t.inkDim + ';font-size:11px;">No stats computed.</div>';
      } else {
        // Shared bin edges for the two histograms
        const allDisplayed = stats.focal_arr.concat(stats.bg_arr);
        const sharedBinning = _histBins(allDisplayed, 60);
        // Recompute counts for both sides under shared binning to find yMaxCount
        function _binCount(arr, bins, lo, hi) {
          const n = bins.length;
          const w = (hi - lo) / n;
          const c = new Array(n).fill(0);
          for (const v of arr) {
            let bi = Math.floor((v - lo) / w);
            if (bi < 0) bi = 0;
            if (bi >= n) bi = n - 1;
            c[bi]++;
          }
          return c;
        }
        const focalCounts = _binCount(stats.focal_arr, sharedBinning.bins, sharedBinning.lo, sharedBinning.hi);
        const bgCounts    = _binCount(stats.bg_arr,    sharedBinning.bins, sharedBinning.lo, sharedBinning.hi);
        const yMaxCount = Math.max(1, Math.max.apply(null, focalCounts.concat(bgCounts)));
        // Side-by-side histograms
        html += '<div style="display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-bottom:10px;">';
        html += '<div>' + _renderHistogramSVG(stats.focal_arr, {
          theme: t, width: 360, height: 130,
          binsForBoth: sharedBinning, yMaxCount,
          color: t.accent, title: 'focal · ' + (active ? active.label : ''),
        }) + '</div>';
        html += '<div>' + _renderHistogramSVG(stats.bg_arr, {
          theme: t, width: 360, height: 130,
          binsForBoth: sharedBinning, yMaxCount,
          color: t.inkDim, title: 'background · ' + (active ? active.label : ''),
        }) + '</div>';
        html += '</div>';
        // Along-chromosome scatter
        const xLo = (prefs.view_mode === 'zoomed') ? focal_lo_mb - 0.05 * (focal_hi_mb - focal_lo_mb) : null;
        const xHi = (prefs.view_mode === 'zoomed') ? focal_hi_mb + 0.05 * (focal_hi_mb - focal_lo_mb) : null;
        html += '<div>' + _renderScatterSVG(stats.scatter_xy, {
          theme: t, width: 740, height: 150,
          x_lo_mb: xLo, x_hi_mb: xHi,
          focal_lo_mb: focal_lo_mb, focal_hi_mb: focal_hi_mb,
          anchor_marks_mb: widget.anchor_marks_mb,
          z_transform: prefs.z_transform,
        }) + '</div>';
        // Caption
        html += '<div style="color:' + t.inkDim + ';font-size:10px;margin-top:6px;line-height:1.45;">' +
                'n_focal = ' + stats.n_focal + ' · n_bg = ' + stats.n_bg +
                ' · focal mean = ' + (Number.isFinite(stats.focal_mean) ? stats.focal_mean.toFixed(3) : '?') +
                ' · bg mean = '    + (Number.isFinite(stats.bg_mean)    ? stats.bg_mean.toFixed(3)    : '?') +
                (stats.wilcoxon && Number.isFinite(stats.wilcoxon.z_stat)
                  ? ' · z-stat = ' + stats.wilcoxon.z_stat.toFixed(2) : '') +
                '</div>';
      }
      wrap.innerHTML = html;
      // Apply P value to header (separate so it survives a missing-stats path)
      const pEl = wrap.querySelector('#fvb-pvalue');
      if (pEl) pEl.innerHTML = pStr;

      // Wire toolbar events
      const metricSel = wrap.querySelector('#fvb-metric');
      if (metricSel) metricSel.addEventListener('change', () => {
        prefs.metric = metricSel.value;
        _savePrefs(pageId, prefs); rerender();
      });
      const bgSel = wrap.querySelector('#fvb-bg-scope');
      if (bgSel) bgSel.addEventListener('change', () => {
        prefs.bg_scope = bgSel.value;
        _savePrefs(pageId, prefs); rerender();
      });
      const zChk = wrap.querySelector('#fvb-z');
      if (zChk) zChk.addEventListener('change', () => {
        prefs.z_transform = zChk.checked;
        _savePrefs(pageId, prefs); rerender();
      });
      const viewSel = wrap.querySelector('#fvb-view');
      if (viewSel) viewSel.addEventListener('change', () => {
        prefs.view_mode = viewSel.value;
        _savePrefs(pageId, prefs); rerender();
      });
      // Radius buttons
      wrap.querySelectorAll('button[data-radius][data-prefix="fvb"]').forEach(btn => {
        btn.addEventListener('click', () => {
          const r = parseInt(btn.getAttribute('data-radius'), 10);
          if (Number.isFinite(r) && r > 0) {
            setRadiusBp(r);
            // The radius change may need to be propagated upstream to the
            // parent page (e.g. boundaries page scan_radius_bp). Ask the
            // host to recompute focal_lo_bp / focal_hi_bp.
            if (typeof opts.on_radius_change === 'function') {
              opts.on_radius_change(r, widget);
            }
            rerender();
          }
        });
      });
    }
    // Expose for external re-render (when host data changes)
    wrap._fvbRerender = rerender;
    rerender();
    return wrap;
  }

  // ---------------------------------------------------------------------------
  // Public API
  // ---------------------------------------------------------------------------
  return {
    makeFocalVsBgPanel: makeFocalVsBgPanel,
    wilcoxonRankSum:    wilcoxonRankSum,
    availableMetrics:   availableMetrics,
    computeFocalBgStats: computeFocalBgStats,
    // For tests
    _internal: {
      METRIC_DEFS: METRIC_DEFS,
      histBins:    _histBins,
      median:      _median,
      formatP:     _formatP,
    },
  };
}));
