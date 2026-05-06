// =============================================================================
// atlas_compare.js — Turn 10 of chat A
// =============================================================================
//
// Set-arithmetic helper for comparing two (or more) sample-id groups.
// Fills a gap that's been called out repeatedly: cross-LG comparisons of
// homozygous arrangements are NOT divergence questions — they're set-overlap
// questions, because there's no shared coordinate system across chromosomes
// to define a meaningful between-chromosome Fst.
//
//   "Is HOM2 at LG28 the same set of individuals as HOM2 at LG14?"
//   "Does the family-pruned HOM2 still cover the same shelf carriers?"
//   "Of the 60 LG28 HOM2 carriers, how many overlap with the 47 LG14 HOM2?"
//
// All of these collapse to set arithmetic + a Fisher's exact independence
// test against the cohort universe.
//
// API (window.popgenCompare):
//   .compareGroups(A, B, opts)
//     → {
//         n_A, n_B, n_intersect, n_only_A, n_only_B, n_union,
//         jaccard, overlap_coef, dice,
//         expected_intersect_under_indep,
//         odds_ratio, fisher_p_two_tailed,
//         relation: 'DISJOINT' | 'NESTED' | 'IDENTICAL' | 'OVERLAPPING',
//         within_universe: number,   // |A∪B| ≤ |universe|
//         universe_size: number      // either provided in opts or |A∪B|
//       }
//
//   .compareNGroups(groupsByName, opts)
//     → {
//         names: string[],
//         sizes: number[],
//         intersect_matrix: number[N][N],
//         jaccard_matrix:   number[N][N],
//         overlap_matrix:   number[N][N],
//       }
//
//   .resolveSlotsToIdSets(slotNames)  — bridge to turn 2's slot system
//   .makeComparePanel(opts)           — DOM for inline comparison view
//   .compareSlotsByName(nameA, nameB) — convenience: lookup slot by name + compare
//
// Pure-function semantics: compareGroups doesn't touch atlas state at all.
// Universe size for Fisher's exact defaults to |A∪B|; pass opts.universe_size
// to use the full cohort instead (the right choice for cross-LG questions on
// the SAME cohort).
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenCompare = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ===========================================================================
  // Engine accessors
  // ===========================================================================

  function _engine() {
    if (typeof window !== 'undefined' && window.popgen) return window.popgen;
    if (typeof globalThis !== 'undefined' && globalThis.popgen) return globalThis.popgen;
    return null;
  }
  function _atlasState() {
    if (typeof window !== 'undefined' && window.state) return window.state;
    if (typeof globalThis !== 'undefined' && globalThis.state) return globalThis.state;
    return null;
  }

  function _cohortSize() {
    const s = _atlasState();
    if (s && s.data && Array.isArray(s.data.samples)) return s.data.samples.length;
    if (s && s.data && typeof s.data.n_samples === 'number') return s.data.n_samples;
    return null;
  }

  // ===========================================================================
  // Set arithmetic
  // ===========================================================================

  function _toIdSet(input) {
    if (input == null) return new Set();
    if (input instanceof Set) return new Set(input);
    if (Array.isArray(input)) return new Set(input.map(String));
    return new Set();
  }

  function _classifyRelation(n_A, n_B, n_intersect, n_union) {
    if (n_intersect === 0) return 'DISJOINT';
    if (n_A === n_B && n_A === n_intersect) return 'IDENTICAL';
    if (n_intersect === n_A || n_intersect === n_B) return 'NESTED';
    return 'OVERLAPPING';
  }

  // ---------------------------------------------------------------------------
  // Log factorials cached up to ceiling needed
  // ---------------------------------------------------------------------------

  let _logFactCache = [0];   // logFact[0] = 0 = log(1)
  function _logFact(n) {
    if (n < 0) return NaN;
    if (n < _logFactCache.length) return _logFactCache[n];
    let last = _logFactCache[_logFactCache.length - 1];
    for (let i = _logFactCache.length; i <= n; i++) {
      last += Math.log(i);
      _logFactCache.push(last);
    }
    return last;
  }

  // log of binomial(n, k)
  function _logBinom(n, k) {
    if (k < 0 || k > n) return -Infinity;
    return _logFact(n) - _logFact(k) - _logFact(n - k);
  }

  // Hypergeometric log-pmf: P(K = k | N, K_white, n_drawn)
  function _logHypergeom(N, K, n, k) {
    if (k < Math.max(0, n - (N - K))) return -Infinity;
    if (k > Math.min(n, K))            return -Infinity;
    return _logBinom(K, k) + _logBinom(N - K, n - k) - _logBinom(N, n);
  }

  // Two-tailed Fisher's exact for a 2x2 table:
  //
  //                in A      not in A
  //     in B        k         n_B - k
  //     not in B  n_A - k    N - n_A - n_B + k
  //
  // P(K = j) sum for all j with P(K = j) ≤ P(K = k_obs).
  function _fisherTwoTailed(N, n_A, n_B, k_obs) {
    if (N <= 0 || n_A < 0 || n_B < 0 || n_A > N || n_B > N) return null;
    const lo = Math.max(0, n_A + n_B - N);
    const hi = Math.min(n_A, n_B);
    if (k_obs < lo || k_obs > hi) return null;
    const log_obs = _logHypergeom(N, n_A, n_B, k_obs);
    if (!isFinite(log_obs)) return null;
    // Sum probabilities ≤ observed using stable log-sum-exp
    // Add a small epsilon for floating-point comparison robustness
    const EPS = 1e-12;
    const tail_logs = [];
    for (let j = lo; j <= hi; j++) {
      const lj = _logHypergeom(N, n_A, n_B, j);
      if (lj <= log_obs + EPS) tail_logs.push(lj);
    }
    if (tail_logs.length === 0) return 1.0;
    // logsumexp
    const max_l = Math.max(...tail_logs);
    let sum = 0;
    for (const l of tail_logs) sum += Math.exp(l - max_l);
    const log_p = max_l + Math.log(sum);
    return Math.min(1, Math.exp(log_p));
  }

  // Odds ratio with Haldane-Anscombe continuity correction (adds 0.5 to each
  // cell when any cell is 0). Returns null if the universe is too small to
  // even define the table.
  function _oddsRatio(N, n_A, n_B, k_obs) {
    if (N - n_A - n_B + k_obs < 0) return null;   // table is invalid
    let a = k_obs;
    let b = n_B - k_obs;
    let c = n_A - k_obs;
    let d = N - n_A - n_B + k_obs;
    if (a === 0 || b === 0 || c === 0 || d === 0) {
      a += 0.5; b += 0.5; c += 0.5; d += 0.5;
    }
    if (b * c === 0) return Infinity;
    return (a * d) / (b * c);
  }

  // ===========================================================================
  // compareGroups — the headline function
  // ===========================================================================

  function compareGroups(A, B, opts) {
    opts = opts || {};
    const setA = _toIdSet(A);
    const setB = _toIdSet(B);
    const n_A = setA.size;
    const n_B = setB.size;

    // Intersection
    const inter = new Set();
    for (const id of setA) if (setB.has(id)) inter.add(id);
    const n_intersect = inter.size;
    const n_only_A    = n_A - n_intersect;
    const n_only_B    = n_B - n_intersect;
    const n_union     = n_A + n_B - n_intersect;

    // Universe size for Fisher's exact independence test:
    // - If user passed a number → use it
    // - If user passed 'cohort' or omitted with cohort context → atlas cohort size
    // - Otherwise → |A∪B| (degenerate case; means we can't reject the null in any
    //   meaningful way, but we still report a nominal P)
    let universe_size;
    if (typeof opts.universe_size === 'number' && opts.universe_size > 0) {
      universe_size = opts.universe_size;
    } else if (opts.universe_size === 'cohort') {
      universe_size = _cohortSize() || n_union;
    } else if (Array.isArray(opts.universe) || opts.universe instanceof Set) {
      universe_size = _toIdSet(opts.universe).size;
    } else {
      universe_size = _cohortSize() || n_union;
    }
    if (universe_size < n_union) universe_size = n_union;   // can't be smaller

    const jaccard      = (n_union > 0)              ? n_intersect / n_union              : null;
    const overlap_coef = (Math.min(n_A, n_B) > 0)   ? n_intersect / Math.min(n_A, n_B)   : null;
    const dice         = ((n_A + n_B) > 0)          ? (2 * n_intersect) / (n_A + n_B)    : null;

    const expected = (universe_size > 0) ? (n_A * n_B) / universe_size : null;
    const odds = _oddsRatio(universe_size, n_A, n_B, n_intersect);
    const fisher_p = _fisherTwoTailed(universe_size, n_A, n_B, n_intersect);
    const relation = _classifyRelation(n_A, n_B, n_intersect, n_union);

    return {
      n_A, n_B, n_intersect, n_only_A, n_only_B, n_union,
      jaccard, overlap_coef, dice,
      expected_intersect_under_indep: expected,
      odds_ratio: odds,
      fisher_p_two_tailed: fisher_p,
      relation,
      universe_size,
      // Sample IDs — useful for downstream chaining
      ids_intersect: Array.from(inter).sort(),
      ids_only_A:    Array.from(setA).filter(id => !setB.has(id)).sort(),
      ids_only_B:    Array.from(setB).filter(id => !setA.has(id)).sort(),
    };
  }

  // ===========================================================================
  // compareNGroups — generalize to N groups, pairwise matrices
  // ===========================================================================

  function compareNGroups(groupsByName, opts) {
    opts = opts || {};
    const names = Object.keys(groupsByName);
    const N = names.length;
    if (N < 2) return null;
    const sets = names.map(n => _toIdSet(groupsByName[n]));
    const sizes = sets.map(s => s.size);
    const intersect_matrix = [];
    const jaccard_matrix   = [];
    const overlap_matrix   = [];
    for (let i = 0; i < N; i++) {
      const interRow = new Array(N).fill(0);
      const jacRow   = new Array(N).fill(0);
      const ovRow    = new Array(N).fill(0);
      for (let j = 0; j < N; j++) {
        if (i === j) {
          interRow[j] = sizes[i];
          jacRow[j]   = 1;
          ovRow[j]    = 1;
        } else {
          let inter = 0;
          for (const id of sets[i]) if (sets[j].has(id)) inter++;
          const uni = sizes[i] + sizes[j] - inter;
          interRow[j] = inter;
          jacRow[j]   = uni > 0 ? inter / uni : 0;
          const minSize = Math.min(sizes[i], sizes[j]);
          ovRow[j]    = minSize > 0 ? inter / minSize : 0;
        }
      }
      intersect_matrix.push(interRow);
      jaccard_matrix.push(jacRow);
      overlap_matrix.push(ovRow);
    }
    return { names, sizes, intersect_matrix, jaccard_matrix, overlap_matrix };
  }

  // ===========================================================================
  // Bridges to turn-2 slot system
  // ===========================================================================

  function resolveSlotsToIdSets(slotNames) {
    const eng = _engine();
    if (!eng) return {};
    const slots = (typeof eng.getSlots === 'function') ? eng.getSlots() : [];
    const out = {};
    for (const name of slotNames) {
      const idx = slots.findIndex(s => s && s.name === name);
      if (idx < 0) continue;
      try {
        out[name] = eng.resolveSlotByIdx(idx);
      } catch (_) {}
    }
    return out;
  }

  function compareSlotsByName(nameA, nameB, opts) {
    const r = resolveSlotsToIdSets([nameA, nameB]);
    if (!r[nameA] || !r[nameB]) return null;
    return compareGroups(r[nameA], r[nameB], opts);
  }

  // ===========================================================================
  // DOM panel — used in the overview's comparison view + as a slot tooltip
  // ===========================================================================

  let _injectedCSS = false;
  function _injectCSS() {
    if (_injectedCSS) return;
    if (typeof document === 'undefined') return;
    if (document.getElementById('popgen-compare-css')) {
      _injectedCSS = true; return;
    }
    const style = document.createElement('style');
    style.id = 'popgen-compare-css';
    style.textContent = COMPARE_CSS;
    document.head.appendChild(style);
    _injectedCSS = true;
  }

  // makeComparePanel({A, B, nameA, nameB, opts}) returns DOM with:
  //   - Header line: A ⟷ B
  //   - 2x2 venn-like glyph (CSS-only): two overlapping circles with counts
  //   - Stats grid: n_A, n_B, intersect, jaccard, overlap, dice, fisher p
  //   - Verdict line classifying the relation
  function makeComparePanel(opts) {
    opts = opts || {};
    if (typeof document === 'undefined') return null;
    _injectCSS();

    const A = opts.A, B = opts.B;
    const nameA = opts.nameA || 'A';
    const nameB = opts.nameB || 'B';
    const result = compareGroups(A, B, opts);

    const root = document.createElement('div');
    root.setAttribute('data-popgen-compare', 'panel');

    // Header
    const header = document.createElement('div');
    header.setAttribute('data-popgen-compare', 'header');
    header.textContent = nameA + ' ⟷ ' + nameB;
    root.appendChild(header);

    // Venn glyph
    const venn = _makeVennGlyph(result, nameA, nameB);
    root.appendChild(venn);

    // Stats grid
    const grid = document.createElement('div');
    grid.setAttribute('data-popgen-compare', 'stats');
    function row(label, value, kind) {
      const r = document.createElement('div');
      r.setAttribute('data-popgen-compare', 'stat-row');
      if (kind) r.setAttribute('data-kind', kind);
      const l = document.createElement('span');
      l.setAttribute('data-popgen-compare', 'stat-label');
      l.textContent = label;
      const v = document.createElement('span');
      v.setAttribute('data-popgen-compare', 'stat-value');
      v.textContent = value;
      r.appendChild(l); r.appendChild(v);
      grid.appendChild(r);
    }
    row('|A|',                 String(result.n_A));
    row('|B|',                 String(result.n_B));
    row('|A ∩ B|',             String(result.n_intersect), 'highlight');
    row('|A \\ B|',            String(result.n_only_A));
    row('|B \\ A|',            String(result.n_only_B));
    row('|A ∪ B|',             String(result.n_union));
    row('Jaccard',             _f3(result.jaccard));
    row('overlap coef.',       _f3(result.overlap_coef));
    row('Dice',                _f3(result.dice));
    row('expected (indep.)',   _f1(result.expected_intersect_under_indep));
    row('odds ratio',          (result.odds_ratio === Infinity) ? '∞' : _f2(result.odds_ratio));
    row('Fisher P (2-tail)',   _formatP(result.fisher_p_two_tailed));
    row('universe N',          String(result.universe_size));
    root.appendChild(grid);

    // Verdict
    const verdict = document.createElement('div');
    verdict.setAttribute('data-popgen-compare', 'verdict');
    verdict.setAttribute('data-relation', result.relation);
    verdict.textContent = _verdictText(result);
    root.appendChild(verdict);

    return root;
  }

  function _makeVennGlyph(r, nameA, nameB) {
    const wrap = document.createElement('div');
    wrap.setAttribute('data-popgen-compare', 'venn');
    // Just two labelled circles in CSS, with the three counts overlaid as
    // text. No SVG — keeps the bundle small.
    const left = document.createElement('div');
    left.setAttribute('data-popgen-compare', 'venn-left');
    const center = document.createElement('div');
    center.setAttribute('data-popgen-compare', 'venn-center');
    const right = document.createElement('div');
    right.setAttribute('data-popgen-compare', 'venn-right');

    left.innerHTML   = '<b>' + r.n_only_A + '</b><br><span>' + nameA + '</span>';
    center.innerHTML = '<b>' + r.n_intersect + '</b><br><span>both</span>';
    right.innerHTML  = '<b>' + r.n_only_B + '</b><br><span>' + nameB + '</span>';

    wrap.appendChild(left);
    wrap.appendChild(center);
    wrap.appendChild(right);
    return wrap;
  }

  function _verdictText(r) {
    if (r.relation === 'IDENTICAL') {
      return 'IDENTICAL — same set of individuals';
    }
    if (r.relation === 'DISJOINT') {
      return 'DISJOINT — no shared individuals (P = ' + _formatP(r.fisher_p_two_tailed) + ')';
    }
    if (r.relation === 'NESTED') {
      const subset = (r.n_intersect === r.n_A) ? 'A ⊆ B' : 'B ⊆ A';
      return 'NESTED — ' + subset + ' (Jaccard ' + _f3(r.jaccard) + ')';
    }
    // OVERLAPPING — describe relative to expected
    if (r.expected_intersect_under_indep != null) {
      const ratio = r.n_intersect / r.expected_intersect_under_indep;
      if (ratio > 1.5) return 'OVERLAPPING — observed ' + r.n_intersect +
        ' is ' + ratio.toFixed(1) + '× expected ' + _f1(r.expected_intersect_under_indep) +
        ' (P = ' + _formatP(r.fisher_p_two_tailed) + ')';
      if (ratio < 0.5) return 'OVERLAPPING — observed ' + r.n_intersect +
        ' is ' + ratio.toFixed(2) + '× expected ' + _f1(r.expected_intersect_under_indep) +
        ' — under-overlap (P = ' + _formatP(r.fisher_p_two_tailed) + ')';
    }
    return 'OVERLAPPING — Jaccard ' + _f3(r.jaccard) +
      ', P = ' + _formatP(r.fisher_p_two_tailed);
  }

  function _f1(x) { return (x == null || isNaN(x)) ? 'NA' : x.toFixed(1); }
  function _f2(x) { return (x == null || isNaN(x)) ? 'NA' : x.toFixed(2); }
  function _f3(x) { return (x == null || isNaN(x)) ? 'NA' : x.toFixed(3); }
  function _formatP(p) {
    if (p == null || isNaN(p)) return 'NA';
    if (p < 1e-10) return '< 1e-10';
    if (p < 1e-3)  return p.toExponential(2);
    return p.toFixed(3);
  }

  // ===========================================================================
  // CSS
  // ===========================================================================

  const COMPARE_CSS = `
    [data-popgen-compare="panel"] {
      display: flex; flex-direction: column;
      background: var(--panel, #fbfcfd);
      border: 1px solid var(--rule, #c8cdd2);
      border-radius: 4px;
      padding: 10px 12px;
      font-family: var(--mono, ui-monospace, monospace);
      font-size: 11px; color: var(--ink, #0e1116);
      max-width: 360px;
    }
    [data-popgen-compare="header"] {
      font-weight: 700; letter-spacing: 0.05em; text-transform: uppercase;
      font-size: 11px; color: var(--ink-dim, #555e69);
      padding-bottom: 6px;
      border-bottom: 1px solid var(--rule, #c8cdd2);
      margin-bottom: 8px;
      text-align: center;
    }
    [data-popgen-compare="venn"] {
      display: grid;
      grid-template-columns: 1fr 1fr 1fr;
      gap: 0;
      padding: 8px 0;
      align-items: center;
      justify-items: center;
      text-align: center;
      font-size: 10px;
    }
    [data-popgen-compare="venn-left"],
    [data-popgen-compare="venn-center"],
    [data-popgen-compare="venn-right"] {
      width: 70px; height: 70px;
      border-radius: 50%;
      display: flex; flex-direction: column; align-items: center; justify-content: center;
      line-height: 1.2;
      font-size: 9px; color: var(--ink, #0e1116);
    }
    [data-popgen-compare="venn-left"] {
      background: rgba(31, 78, 121, 0.18);
      border: 1.5px solid #1f4e79;
      margin-right: -22px; z-index: 1;
    }
    [data-popgen-compare="venn-center"] {
      background: rgba(245, 165, 36, 0.30);
      width: 50px; height: 50px;
      border: 1.5px solid var(--accent, #f5a524);
      z-index: 2;
      font-weight: 700;
    }
    [data-popgen-compare="venn-right"] {
      background: rgba(60, 192, 138, 0.18);
      border: 1.5px solid #3cc08a;
      margin-left: -22px; z-index: 1;
    }
    [data-popgen-compare="venn"] b {
      font-size: 16px; font-weight: 700; color: var(--ink, #0e1116);
    }
    [data-popgen-compare="venn"] span {
      font-size: 8px; color: var(--ink-dim, #555e69); margin-top: 2px;
    }
    [data-popgen-compare="stats"] {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 2px 16px;
      padding: 6px 0;
      border-top: 1px solid var(--rule, #c8cdd2);
      border-bottom: 1px solid var(--rule, #c8cdd2);
      margin-bottom: 8px;
    }
    [data-popgen-compare="stat-row"] {
      display: flex; justify-content: space-between;
      padding: 1px 0;
      font-variant-numeric: tabular-nums;
      font-size: 10px;
    }
    [data-popgen-compare="stat-row"][data-kind="highlight"] {
      background: rgba(245, 165, 36, 0.10);
      font-weight: 600;
    }
    [data-popgen-compare="stat-label"] {
      color: var(--ink-dim, #555e69);
    }
    [data-popgen-compare="stat-value"] {
      color: var(--ink, #0e1116);
    }
    [data-popgen-compare="verdict"] {
      padding: 6px 8px;
      border-radius: 3px;
      font-size: 10px;
      text-align: center;
    }
    [data-popgen-compare="verdict"][data-relation="IDENTICAL"] {
      background: rgba(60, 192, 138, 0.18);
      color: #1b5b3e;
      border: 1px solid #3cc08a;
      font-weight: 600;
    }
    [data-popgen-compare="verdict"][data-relation="NESTED"] {
      background: rgba(60, 192, 138, 0.10);
      color: #1b5b3e;
      border: 1px solid #3cc08a;
    }
    [data-popgen-compare="verdict"][data-relation="DISJOINT"] {
      background: rgba(224, 85, 92, 0.13);
      color: #6e1b1f;
      border: 1px solid #e0555c;
    }
    [data-popgen-compare="verdict"][data-relation="OVERLAPPING"] {
      background: rgba(245, 165, 36, 0.13);
      color: #6d4707;
      border: 1px solid var(--accent, #f5a524);
    }
  `;

  // ===========================================================================
  // Public exports
  // ===========================================================================

  return {
    compareGroups,
    compareNGroups,
    resolveSlotsToIdSets,
    compareSlotsByName,
    makeComparePanel,
    _internals: {
      _logFact, _logBinom, _logHypergeom,
      _fisherTwoTailed, _oddsRatio,
      _classifyRelation,
      _toIdSet,
    },
  };
}));
