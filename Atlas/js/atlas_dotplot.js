// =============================================================================
// atlas_dotplot.js — cross-species dot plot widget (Quentin's manuscript)
// -----------------------------------------------------------------------------
// Self-contained module. Zero external dependencies. Renders a synteny dot
// plot from either:
//   (a) wfmash synteny_blocks (already inside cs_breakpoints_v1.json v2), or
//   (b) mashmap multi-resolution grid (dotplot_mashmap_v1.json).
//
// Visual idiom: mummerplot-flavour, but modernised for the atlas.
//   - Each chromosome pair is a sub-rectangle inside the global axes.
//     Boxes drawn with a 1px rule line; chrom labels in the gutter.
//   - Each alignment is a line segment (start_q, start_t) -> (end_q, end_t)
//     within its chrom-pair box.
//   - Color encodes percent identity on a perceptual ramp that lives in the
//     atlas's accent family (low PI = blue, high PI = amber-red).
//   - Strand encoded by *texture*, not color: forward = solid α=0.85,
//     reverse = solid α=0.55 with a 1.5px dotted overlay.
//   - Layout: chroms reordered on both axes so the heaviest synteny falls on
//     the main diagonal (LIS-flavour layout, ported from the Perl ref).
//
// Display modes (state machine):
//   MINI ──hover──> ENLARGED ──leave──> MINI
//
// (Earlier versions had a PINNED state on click; removed per Quentin's
// preference for a pure-hover model — open on hover, close on un-hover,
// no separate pin step.)
//
// Public API (window.popgenDotplot):
//   .makeDotplotPanel(opts) -> HTMLElement
//       opts: {
//         containerStyle?,                // optional inline CSS string for the wrapper
//         data: {
//           species_query:  { name, haplotype },
//           species_target: { name, haplotype },
//           chrom_lengths_query:  { LG28: 38500000, ... },
//           chrom_lengths_target: { Mac_27: 41200000, ... },
//           // Either one source OR both:
//           wfmash_blocks?: [ { gar_chr, gar_start, gar_end, mac_chr, mac_start,
//                                mac_end, strand, block_size_bp,
//                                mapping_quality, pi? } ],
//           mashmap_resolutions?: [ { label, segment_length_bp, pi_min,
//                                      n_alignments, alignments: [
//                                        { q_chr, q_start, q_end, t_chr,
//                                          t_start, t_end, strand, pi } ] } ],
//         },
//         defaultResolution?,             // resolution label, default 'wfmash'
//         miniSize?,                      // px, default 280
//         enlargedSize?,                  // px, default 720
//         onAlignmentClick?,              // (alignment, ev) => void
//       }
//   .colormap(pi)        -> "#rrggbb"
//   .layoutChroms(...)   -> { qOrder, tOrder, qFlipped, tFlipped }
//
// =============================================================================

(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) {
    module.exports = factory();
  } else {
    root.popgenDotplot = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Color ramp — perceptual, on-theme (atlas accent family)
  // 5 stops from low PI (cool slate-blue) to high PI (warm amber-red).
  // PI ∈ [0.85, 1.0]; values outside clamp to ramp ends.
  // ---------------------------------------------------------------------------
  const RAMP_STOPS = [
    { pi: 0.85, hex: '#3a5a78' },  // muted slate
    { pi: 0.90, hex: '#3f7a4d' },  // sage
    { pi: 0.94, hex: '#c08a00' },  // amber (atlas accent base)
    { pi: 0.97, hex: '#c2541d' },  // burnt orange
    { pi: 1.00, hex: '#a8201a' },  // brick red
  ];

  function _hex2rgb(h) {
    h = h.replace('#', '');
    return [parseInt(h.slice(0, 2), 16),
            parseInt(h.slice(2, 4), 16),
            parseInt(h.slice(4, 6), 16)];
  }
  function _rgb2hex(r, g, b) {
    const t = (x) => Math.round(Math.max(0, Math.min(255, x)))
      .toString(16).padStart(2, '0');
    return '#' + t(r) + t(g) + t(b);
  }

  function colormap(pi) {
    if (!Number.isFinite(pi)) return RAMP_STOPS[0].hex;
    const lo = RAMP_STOPS[0].pi;
    const hi = RAMP_STOPS[RAMP_STOPS.length - 1].pi;
    if (pi <= lo) return RAMP_STOPS[0].hex;
    if (pi >= hi) return RAMP_STOPS[RAMP_STOPS.length - 1].hex;
    // Find bracketing stops
    for (let i = 0; i < RAMP_STOPS.length - 1; i++) {
      const a = RAMP_STOPS[i];
      const b = RAMP_STOPS[i + 1];
      if (pi >= a.pi && pi <= b.pi) {
        const t = (pi - a.pi) / (b.pi - a.pi);
        const ar = _hex2rgb(a.hex);
        const br = _hex2rgb(b.hex);
        return _rgb2hex(
          ar[0] + (br[0] - ar[0]) * t,
          ar[1] + (br[1] - ar[1]) * t,
          ar[2] + (br[2] - ar[2]) * t,
        );
      }
    }
    return RAMP_STOPS[RAMP_STOPS.length - 1].hex;
  }

  // ---------------------------------------------------------------------------
  // Chrom layout — LIS-flavour fattest-alignment diagonal.
  // Algorithm (paraphrased from Perl LayoutIDs/SpanXwY):
  //   1. For each query chrom Q, find the dominant target chrom T(Q):
  //      the target carrying the most aligned bp from Q.
  //   2. Symmetric for targets.
  //   3. Order queries by total aligned bp (descending). For each Q,
  //      assign T(Q) the next available target slot. If T(Q)'s dominant
  //      query is NOT this Q, defer.
  //   4. Append unplaced chroms at the end.
  //   5. Flip a target's orientation if its dominant alignment is
  //      reverse-strand.
  // Returns { qOrder: [chrom...], tOrder: [chrom...],
  //          qFlipped: Set, tFlipped: Set }
  // qFlipped is always empty (queries don't flip — convention); tFlipped
  // marks target chroms drawn reversed (high coord at bottom of the box).
  // ---------------------------------------------------------------------------
  function layoutChroms(alignments, qLens, tLens) {
    if (!alignments || alignments.length === 0) {
      return {
        qOrder: Object.keys(qLens || {}).sort(),
        tOrder: Object.keys(tLens || {}).sort(),
        qFlipped: new Set(),
        tFlipped: new Set(),
      };
    }
    // Aggregate aligned bp per (q, t) pair, separately for forward and reverse
    const pairCov = new Map();      // "q__t" -> { fwd: bp, rev: bp, total: bp }
    const qTot = new Map();
    const tTot = new Map();
    for (const a of alignments) {
      const span = Math.max(1, a.q_end - a.q_start);
      const k = a.q_chr + '__' + a.t_chr;
      const cur = pairCov.get(k) || { fwd: 0, rev: 0, total: 0 };
      if (a.strand === '+') cur.fwd += span;
      else cur.rev += span;
      cur.total += span;
      pairCov.set(k, cur);
      qTot.set(a.q_chr, (qTot.get(a.q_chr) || 0) + span);
      tTot.set(a.t_chr, (tTot.get(a.t_chr) || 0) + span);
    }
    // For each chrom, dominant counterpart (max total)
    const qDom = new Map();   // q -> { t, span, fwd, rev }
    const tDom = new Map();   // t -> { q, span, fwd, rev }
    for (const [k, v] of pairCov) {
      const [q, t] = k.split('__');
      const cq = qDom.get(q);
      if (!cq || v.total > cq.span) qDom.set(q, { t, span: v.total, fwd: v.fwd, rev: v.rev });
      const ct = tDom.get(t);
      if (!ct || v.total > ct.span) tDom.set(t, { q, span: v.total, fwd: v.fwd, rev: v.rev });
    }
    // Order queries by total aligned bp, descending
    const qSorted = Array.from(qTot.keys()).sort((a, b) => qTot.get(b) - qTot.get(a));
    // Assign target order by walking qSorted; T(Q) goes next if its own
    // dominant query is Q (mutual). Otherwise defer the target.
    const qOrder = [];
    const tOrder = [];
    const tPlaced = new Set();
    const tFlipped = new Set();
    for (const q of qSorted) {
      qOrder.push(q);
      const dom = qDom.get(q);
      if (!dom) continue;
      const t = dom.t;
      if (tPlaced.has(t)) continue;
      // Mutual? T's dominant query == this q?
      const tdq = tDom.get(t);
      if (tdq && tdq.q === q) {
        tOrder.push(t);
        tPlaced.add(t);
        // Flip target if reverse dominates
        if (dom.rev > dom.fwd) tFlipped.add(t);
      }
    }
    // Place remaining targets — by descending tTot, in order
    const remaining = Array.from(tTot.keys())
      .filter(t => !tPlaced.has(t))
      .sort((a, b) => tTot.get(b) - tTot.get(a));
    for (const t of remaining) {
      tOrder.push(t);
      tPlaced.add(t);
      const tdq = tDom.get(t);
      if (tdq && tdq.rev > tdq.fwd) tFlipped.add(t);
    }
    // Append query chroms with no alignments at all (so axes show every chrom)
    for (const q of Object.keys(qLens || {})) {
      if (qOrder.indexOf(q) === -1) qOrder.push(q);
    }
    for (const t of Object.keys(tLens || {})) {
      if (tOrder.indexOf(t) === -1) tOrder.push(t);
    }
    return { qOrder, tOrder, qFlipped: new Set(), tFlipped };
  }

  // ---------------------------------------------------------------------------
  // Cumulative-offset helper — given chrom order + lengths, produce a Map
  // chrom -> { offset_bp, length_bp } so a (chrom, pos_bp) coord maps to a
  // global axis position via offset + (flipped ? length - pos : pos).
  // ---------------------------------------------------------------------------
  function _cumulOffsets(order, lens, flipped) {
    const m = new Map();
    let off = 0;
    for (const c of order) {
      const L = lens[c] || 0;
      m.set(c, { offset_bp: off, length_bp: L, flipped: flipped.has(c) });
      off += L;
    }
    return { offsets: m, total_bp: off };
  }

  // ---------------------------------------------------------------------------
  // Theme accessor — read CSS vars from the document root so the dot plot
  // tracks the active atlas theme (dark/light/ivory) automatically.
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
      panel3:    get('--panel-3',    '#2a3242'),
      ink:       get('--ink',        '#e7edf3'),
      inkDim:    get('--ink-dim',    '#8a94a3'),
      inkDimmer: get('--ink-dimmer', '#5a6472'),
      rule:      get('--rule',       '#2a3242'),
      accent:    get('--accent',     '#f5a524'),
      mono:      get('--mono',       'ui-monospace, monospace'),
    };
  }

  // ---------------------------------------------------------------------------
  // Pick the active alignments array based on resolution selection.
  // ---------------------------------------------------------------------------
  function _pickResolution(data, resolutionLabel) {
    if (resolutionLabel === 'wfmash' || resolutionLabel === 'wfmash_synteny') {
      // Normalize wfmash blocks to the unified shape
      if (!data.wfmash_blocks) return null;
      return data.wfmash_blocks.map(b => ({
        q_chr: b.gar_chr, q_start: b.gar_start, q_end: b.gar_end,
        t_chr: b.mac_chr, t_start: b.mac_start, t_end: b.mac_end,
        strand: b.strand,
        // wfmash block has mapping_quality but not always pi; synthesize a
        // reasonable PI from match content if available, else default 0.97
        // (high — wfmash blocks are typically very high identity).
        pi: (typeof b.pi === 'number') ? b.pi : 0.97,
        block_size_bp: b.block_size_bp,
        mapping_quality: b.mapping_quality,
      }));
    }
    if (!data.mashmap_resolutions) return null;
    const r = data.mashmap_resolutions.find(x => x.label === resolutionLabel);
    if (!r) return null;
    return r.alignments;
  }

  // ---------------------------------------------------------------------------
  // Resolution dropdown options builder.
  // ---------------------------------------------------------------------------
  function _resolutionOptions(data) {
    const opts = [];
    if (data.wfmash_blocks && data.wfmash_blocks.length > 0) {
      opts.push({ label: 'wfmash synteny (precise)', value: 'wfmash' });
    }
    if (data.mashmap_resolutions) {
      for (const r of data.mashmap_resolutions) {
        opts.push({
          label: 'mashmap · ' + r.label.replace('_', ' · ') +
                 ' · n=' + (r.n_alignments || (r.alignments || []).length),
          value: r.label,
        });
      }
    }
    return opts;
  }

  // ---------------------------------------------------------------------------
  // Build SVG markup for the dot plot at a given size.
  // size: integer px (square plot). Returns an outer string of SVG.
  // ---------------------------------------------------------------------------
  function _renderDotplotSVG(data, alignments, layout, size, opts) {
    const t = _theme();
    const PAD = { top: 22, right: 18, bottom: 56, left: 56 };
    const W = size;
    const H = size;
    const innerW = W - PAD.left - PAD.right;
    const innerH = H - PAD.top - PAD.bottom;

    const qInfo = _cumulOffsets(layout.qOrder, data.chrom_lengths_query,  layout.qFlipped);
    const tInfo = _cumulOffsets(layout.tOrder, data.chrom_lengths_target, layout.tFlipped);
    const qTotal = qInfo.total_bp || 1;
    const tTotal = tInfo.total_bp || 1;

    // Coord mappers (bp → svg px). Y increases downward in SVG, but
    // dot-plot convention is target-chrom-1 at top, so we map directly.
    const xOf = (chr, pos) => {
      const e = qInfo.offsets.get(chr);
      if (!e) return null;
      const p = e.flipped ? (e.length_bp - pos) : pos;
      return PAD.left + (e.offset_bp + p) / qTotal * innerW;
    };
    const yOf = (chr, pos) => {
      const e = tInfo.offsets.get(chr);
      if (!e) return null;
      const p = e.flipped ? (e.length_bp - pos) : pos;
      return PAD.top + (e.offset_bp + p) / tTotal * innerH;
    };

    let svg = '';
    svg += '<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W + '" height="' + H + '" ' +
           'xmlns="http://www.w3.org/2000/svg" ' +
           'style="display:block;background:' + t.panel2 + ';font-family:' + t.mono + ';">';

    // Plot area background (slightly contrasted from the panel)
    svg += '<rect x="' + PAD.left + '" y="' + PAD.top + '" ' +
           'width="' + innerW + '" height="' + innerH + '" ' +
           'fill="' + t.panel + '" stroke="' + t.rule + '" stroke-width="0.8"/>';

    // Chrom-pair grid (vertical lines at query chrom boundaries; horizontal
    // at target chrom boundaries). Subtle — we want the alignments to pop.
    let cumQ = 0;
    for (const c of layout.qOrder) {
      cumQ += (data.chrom_lengths_query[c] || 0);
      if (cumQ >= qTotal) break;
      const x = PAD.left + cumQ / qTotal * innerW;
      svg += '<line x1="' + x.toFixed(1) + '" y1="' + PAD.top + '" ' +
             'x2="' + x.toFixed(1) + '" y2="' + (PAD.top + innerH) + '" ' +
             'stroke="' + t.rule + '" stroke-width="0.4" stroke-opacity="0.6"/>';
    }
    let cumT = 0;
    for (const c of layout.tOrder) {
      cumT += (data.chrom_lengths_target[c] || 0);
      if (cumT >= tTotal) break;
      const y = PAD.top + cumT / tTotal * innerH;
      svg += '<line x1="' + PAD.left + '" y1="' + y.toFixed(1) + '" ' +
             'x2="' + (PAD.left + innerW) + '" y2="' + y.toFixed(1) + '" ' +
             'stroke="' + t.rule + '" stroke-width="0.4" stroke-opacity="0.6"/>';
    }

    // Alignment line segments — sort by PI ascending so high-PI hits draw on
    // top (the eye picks up the strong matches over the weak ones).
    const sorted = alignments.slice().sort((a, b) =>
      (a.pi || 0) - (b.pi || 0));

    // For dense plots, batch into one path-per-color-bin to keep DOM small.
    // 16 PI bins is plenty.
    const nBins = 16;
    const binsFwd = Array.from({ length: nBins }, () => []);
    const binsRev = Array.from({ length: nBins }, () => []);
    for (const a of sorted) {
      const x1 = xOf(a.q_chr, a.q_start);
      const x2 = xOf(a.q_chr, a.q_end);
      const y1 = yOf(a.t_chr, a.t_start);
      const y2 = yOf(a.t_chr, a.t_end);
      if (x1 == null || x2 == null || y1 == null || y2 == null) continue;
      const pi = Number.isFinite(a.pi) ? a.pi : 0.95;
      const rampLo = 0.85, rampHi = 1.0;
      const f = Math.max(0, Math.min(1, (pi - rampLo) / (rampHi - rampLo)));
      const bin = Math.min(nBins - 1, Math.floor(f * nBins));
      const seg = 'M' + x1.toFixed(1) + ' ' + y1.toFixed(1) +
                  ' L' + x2.toFixed(1) + ' ' + y2.toFixed(1);
      ((a.strand === '+') ? binsFwd : binsRev)[bin].push(seg);
    }
    // Render reverse first (so forward is on top — main diagonal stays
    // legible when both are present in the same chrom-pair box).
    for (let b = 0; b < nBins; b++) {
      if (binsRev[b].length === 0) continue;
      const piMid = 0.85 + (b + 0.5) / nBins * 0.15;
      const stroke = colormap(piMid);
      // Solid pass (lower opacity), then dotted overlay
      svg += '<path d="' + binsRev[b].join(' ') + '" ' +
             'fill="none" stroke="' + stroke + '" stroke-width="' + (size > 400 ? 1.2 : 0.9) + '" ' +
             'stroke-opacity="0.55" stroke-linecap="round"/>';
      svg += '<path d="' + binsRev[b].join(' ') + '" ' +
             'fill="none" stroke="' + stroke + '" stroke-width="' + (size > 400 ? 1.2 : 0.9) + '" ' +
             'stroke-opacity="0.85" stroke-dasharray="1.5,2.2" stroke-linecap="round"/>';
    }
    for (let b = 0; b < nBins; b++) {
      if (binsFwd[b].length === 0) continue;
      const piMid = 0.85 + (b + 0.5) / nBins * 0.15;
      const stroke = colormap(piMid);
      svg += '<path d="' + binsFwd[b].join(' ') + '" ' +
             'fill="none" stroke="' + stroke + '" stroke-width="' + (size > 400 ? 1.2 : 0.9) + '" ' +
             'stroke-opacity="0.85" stroke-linecap="round"/>';
    }

    // Chrom labels (only when size is large enough — mini panel skips them)
    const showLabels = size >= 480;
    if (showLabels) {
      // X-axis (query) labels
      let cum = 0;
      for (const c of layout.qOrder) {
        const L = data.chrom_lengths_query[c] || 0;
        const cx = PAD.left + (cum + L / 2) / qTotal * innerW;
        cum += L;
        svg += '<text x="' + cx.toFixed(1) + '" y="' + (PAD.top + innerH + 14) + '" ' +
               'fill="' + t.inkDim + '" font-size="9" text-anchor="middle" ' +
               'transform="rotate(-45 ' + cx.toFixed(1) + ' ' + (PAD.top + innerH + 14) + ')">' +
               _esc(c) + '</text>';
      }
      // Y-axis (target) labels
      cum = 0;
      for (const c of layout.tOrder) {
        const L = data.chrom_lengths_target[c] || 0;
        const cy = PAD.top + (cum + L / 2) / tTotal * innerH;
        cum += L;
        const flip = layout.tFlipped.has(c) ? '*' : '';
        svg += '<text x="' + (PAD.left - 6) + '" y="' + (cy + 3).toFixed(1) + '" ' +
               'fill="' + t.inkDim + '" font-size="9" text-anchor="end">' +
               _esc(flip + c) + '</text>';
      }
    } else {
      // Mini: just show "query" and "target" axis labels in the gutter
      svg += '<text x="' + (PAD.left + innerW / 2) + '" y="' + (H - 8) + '" ' +
             'fill="' + t.inkDim + '" font-size="10" text-anchor="middle" ' +
             'font-style="italic">' + _esc(data.species_query.name || 'query') + '</text>';
      svg += '<text transform="translate(14, ' + (PAD.top + innerH / 2) +
             ') rotate(-90)" fill="' + t.inkDim + '" font-size="10" ' +
             'text-anchor="middle" font-style="italic">' +
             _esc(data.species_target.name || 'target') + '</text>';
    }

    // Color legend (bottom-right, small) — only on enlarged
    if (showLabels) {
      const legW = 120, legH = 8;
      const legX = W - PAD.right - legW;
      const legY = PAD.top - 14;
      const gradId = 'dpgrad' + Math.floor(Math.random() * 1e9);
      svg += '<defs><linearGradient id="' + gradId + '" x1="0" x2="1">';
      for (const s of RAMP_STOPS) {
        const t01 = (s.pi - 0.85) / (1.0 - 0.85);
        svg += '<stop offset="' + (t01 * 100).toFixed(0) + '%" stop-color="' + s.hex + '"/>';
      }
      svg += '</linearGradient></defs>';
      svg += '<rect x="' + legX + '" y="' + legY + '" width="' + legW +
             '" height="' + legH + '" fill="url(#' + gradId + ')" stroke="' +
             t.rule + '" stroke-width="0.5"/>';
      svg += '<text x="' + legX + '" y="' + (legY - 3) + '" fill="' + t.inkDim +
             '" font-size="8.5">% identity</text>';
      svg += '<text x="' + legX + '" y="' + (legY + legH + 9) + '" fill="' +
             t.inkDimmer + '" font-size="8" text-anchor="start">85</text>';
      svg += '<text x="' + (legX + legW) + '" y="' + (legY + legH + 9) + '" fill="' +
             t.inkDimmer + '" font-size="8" text-anchor="end">100</text>';
    }

    svg += '</svg>';
    return svg;
  }

  // ---------------------------------------------------------------------------
  // Tooltip rendering — find the alignment closest to (x, y) in the SVG
  // and return a small HTML tooltip string. Used in enlarged mode.
  // ---------------------------------------------------------------------------
  function _findNearestAlignment(alignments, layout, data, sizePx, mouseX, mouseY) {
    const PAD = { top: 22, right: 18, bottom: 56, left: 56 };
    const innerW = sizePx - PAD.left - PAD.right;
    const innerH = sizePx - PAD.top - PAD.bottom;
    const qInfo = _cumulOffsets(layout.qOrder, data.chrom_lengths_query,  layout.qFlipped);
    const tInfo = _cumulOffsets(layout.tOrder, data.chrom_lengths_target, layout.tFlipped);
    const qTotal = qInfo.total_bp || 1;
    const tTotal = tInfo.total_bp || 1;
    const xOf = (chr, pos) => {
      const e = qInfo.offsets.get(chr);
      if (!e) return null;
      const p = e.flipped ? (e.length_bp - pos) : pos;
      return PAD.left + (e.offset_bp + p) / qTotal * innerW;
    };
    const yOf = (chr, pos) => {
      const e = tInfo.offsets.get(chr);
      if (!e) return null;
      const p = e.flipped ? (e.length_bp - pos) : pos;
      return PAD.top + (e.offset_bp + p) / tTotal * innerH;
    };
    let best = null;
    let bestD = Infinity;
    for (const a of alignments) {
      const x1 = xOf(a.q_chr, a.q_start), x2 = xOf(a.q_chr, a.q_end);
      const y1 = yOf(a.t_chr, a.t_start), y2 = yOf(a.t_chr, a.t_end);
      if (x1 == null || y1 == null) continue;
      // Distance from point to the line segment (x1,y1)-(x2,y2)
      const dx = x2 - x1, dy = y2 - y1;
      const L2 = dx * dx + dy * dy;
      let t01 = L2 > 0 ? ((mouseX - x1) * dx + (mouseY - y1) * dy) / L2 : 0;
      t01 = Math.max(0, Math.min(1, t01));
      const px = x1 + t01 * dx, py = y1 + t01 * dy;
      const d = Math.hypot(mouseX - px, mouseY - py);
      if (d < bestD) { bestD = d; best = a; }
    }
    return (bestD < 6) ? best : null;
  }

  // ---------------------------------------------------------------------------
  // HTML escaping + bp formatter
  // ---------------------------------------------------------------------------
  function _esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  }
  function _fmtBp(bp) {
    if (bp == null || !Number.isFinite(bp)) return '?';
    if (bp >= 1e6) return (bp / 1e6).toFixed(2) + ' Mb';
    if (bp >= 1e3) return (bp / 1e3).toFixed(0) + ' kb';
    return bp + ' bp';
  }

  // ---------------------------------------------------------------------------
  // The widget itself.
  // ---------------------------------------------------------------------------
  function makeDotplotPanel(opts) {
    opts = opts || {};
    const data = opts.data;
    if (!data) {
      const empty = document.createElement('div');
      empty.style.cssText = 'padding:20px;color:var(--ink-dim);font-family:var(--mono);font-size:11px;';
      empty.textContent = 'No dot plot data loaded.';
      return empty;
    }
    const t = _theme();
    const miniSize = opts.miniSize || 280;
    const enlargedSize = opts.enlargedSize || 720;
    const resOptions = _resolutionOptions(data);
    let resolution = opts.defaultResolution ||
      (resOptions[0] ? resOptions[0].value : 'wfmash');

    // Compute layout once per (data, resolution) pair — cached
    const layoutCache = new Map();
    const computeLayoutFor = (resLabel) => {
      const cached = layoutCache.get(resLabel);
      if (cached) return cached;
      const aln = _pickResolution(data, resLabel) || [];
      const lay = layoutChroms(aln, data.chrom_lengths_query, data.chrom_lengths_target);
      const v = { alignments: aln, layout: lay };
      layoutCache.set(resLabel, v);
      return v;
    };

    // ----- Outer wrapper -----
    const wrap = document.createElement('div');
    wrap.className = 'dotplot-panel';
    wrap.style.cssText = (opts.containerStyle || '') +
      ';background:' + t.panel2 + ';border:1px solid ' + t.rule +
      ';border-radius:4px;padding:10px 12px;font-family:' + t.mono +
      ';color:' + t.ink + ';';

    // ----- Header -----
    const header = document.createElement('div');
    header.style.cssText = 'display:flex;align-items:center;gap:10px;flex-wrap:wrap;margin-bottom:8px;font-size:11px;';
    const title = document.createElement('span');
    title.style.cssText = 'color:' + t.ink + ';';
    title.innerHTML = 'Dot plot &middot; <i>' + _esc(data.species_query.name) +
      '</i> &times; <i>' + _esc(data.species_target.name) + '</i>';
    header.appendChild(title);

    const sep = document.createElement('span');
    sep.style.cssText = 'color:' + t.inkDimmer + ';';
    sep.textContent = '·';
    header.appendChild(sep);

    const resLabel = document.createElement('span');
    resLabel.style.cssText = 'color:' + t.inkDim + ';';
    resLabel.textContent = 'resolution:';
    header.appendChild(resLabel);

    const resSelect = document.createElement('select');
    resSelect.style.cssText = 'background:' + t.panel + ';color:' + t.ink +
      ';border:1px solid ' + t.rule + ';border-radius:3px;padding:2px 6px;font-family:' +
      t.mono + ';font-size:11px;';
    for (const o of resOptions) {
      const optEl = document.createElement('option');
      optEl.value = o.value;
      optEl.textContent = o.label;
      if (o.value === resolution) optEl.selected = true;
      resSelect.appendChild(optEl);
    }
    header.appendChild(resSelect);

    // Hint text on the right
    const hint = document.createElement('span');
    hint.style.cssText = 'color:' + t.inkDimmer + ';font-size:10px;margin-left:auto;';
    hint.textContent = 'hover to enlarge';
    header.appendChild(hint);

    wrap.appendChild(header);

    // ----- Mini SVG holder -----
    const miniHolder = document.createElement('div');
    miniHolder.style.cssText = 'display:flex;justify-content:center;cursor:pointer;';
    wrap.appendChild(miniHolder);

    const renderMini = () => {
      const v = computeLayoutFor(resolution);
      miniHolder.innerHTML = _renderDotplotSVG(data, v.alignments, v.layout, miniSize, opts);
    };
    renderMini();
    resSelect.addEventListener('change', () => {
      resolution = resSelect.value;
      renderMini();
      // If popup is open, re-render that too
      if (popupState.open) renderEnlarged();
    });

    // ----- Popup (created lazily) -----
    // Pure hover model: open on mouseenter, close on mouseleave (with a
    // small grace period for moving cursor between mini and popup).
    const popupState = { open: false, el: null, tipEl: null };

    const closePopup = () => {
      if (!popupState.open) return;
      if (popupState.el && popupState.el.parentNode) {
        popupState.el.parentNode.removeChild(popupState.el);
      }
      popupState.el = null;
      popupState.tipEl = null;
      popupState.open = false;
    };

    const renderEnlarged = () => {
      const v = computeLayoutFor(resolution);
      // Build (or reuse) popup element
      if (!popupState.el) {
        const backdrop = document.createElement('div');
        backdrop.className = 'dotplot-popup-backdrop';
        backdrop.style.cssText =
          'position:fixed;inset:0;background:rgba(10,12,14,0.45);' +
          'z-index:9998;display:flex;align-items:center;justify-content:center;' +
          'pointer-events:none;';   // pointer-events flipped on pin
        const popup = document.createElement('div');
        popup.className = 'dotplot-popup';
        popup.style.cssText =
          'position:relative;background:' + t.panel2 + ';' +
          'border:1px solid ' + t.rule + ';border-radius:6px;' +
          'box-shadow:0 12px 40px rgba(0,0,0,0.5),0 2px 6px rgba(0,0,0,0.3);' +
          'padding:14px 16px 16px;font-family:' + t.mono +
          ';color:' + t.ink + ';pointer-events:auto;';
        // Header inside popup
        const ph = document.createElement('div');
        ph.style.cssText = 'display:flex;align-items:center;gap:10px;margin-bottom:8px;font-size:11px;';
        ph.innerHTML =
          '<span><i>' + _esc(data.species_query.name) + '</i> &times; ' +
          '<i>' + _esc(data.species_target.name) + '</i> &middot; ' +
          _esc(resolution) + '</span>';
        const closeBtn = document.createElement('button');
        closeBtn.type = 'button';
        closeBtn.title = 'Close (Esc)';
        closeBtn.textContent = '×';
        closeBtn.style.cssText =
          'margin-left:auto;background:transparent;border:1px solid ' + t.rule +
          ';color:' + t.ink + ';font-size:18px;line-height:1;width:26px;height:26px;' +
          'border-radius:3px;cursor:pointer;font-family:' + t.mono + ';';
        closeBtn.addEventListener('click', closePopup);
        ph.appendChild(closeBtn);
        popup.appendChild(ph);
        // Plot area
        const plotHolder = document.createElement('div');
        plotHolder.className = 'dotplot-popup-plot';
        plotHolder.style.cssText = 'position:relative;';
        popup.appendChild(plotHolder);
        // Tooltip
        const tip = document.createElement('div');
        tip.className = 'dotplot-tip';
        tip.style.cssText =
          'position:absolute;pointer-events:none;background:' + t.panel +
          ';color:' + t.ink + ';border:1px solid ' + t.rule + ';border-radius:3px;' +
          'padding:5px 8px;font-size:10.5px;line-height:1.45;font-family:' +
          t.mono + ';display:none;z-index:1;white-space:nowrap;' +
          'box-shadow:0 2px 8px rgba(0,0,0,0.35);';
        plotHolder.appendChild(tip);
        // Caption (status line)
        const caption = document.createElement('div');
        caption.className = 'dotplot-caption';
        caption.style.cssText = 'color:' + t.inkDim +
          ';font-size:10px;margin-top:8px;line-height:1.4;';
        caption.textContent = (v.alignments.length ? v.alignments.length : 0) +
          ' alignments · ' + v.layout.qOrder.length + ' query × ' +
          v.layout.tOrder.length + ' target chroms · * = target reverse-flipped';
        popup.appendChild(caption);
        backdrop.appendChild(popup);
        document.body.appendChild(backdrop);
        popupState.el = backdrop;
        popupState.tipEl = tip;
        popupState.plotEl = plotHolder;
        popupState.popupInner = popup;
        popupState.captionEl = caption;
        popupState.headerEl = ph;
        // Mouse move for tooltip
        plotHolder.addEventListener('mousemove', (ev) => {
          const rect = plotHolder.getBoundingClientRect();
          const mx = ev.clientX - rect.left;
          const my = ev.clientY - rect.top;
          const v2 = computeLayoutFor(resolution);
          const hit = _findNearestAlignment(v2.alignments, v2.layout, data,
            enlargedSize, mx, my);
          if (hit) {
            const piStr = (hit.pi != null) ? (hit.pi * 100).toFixed(1) + '%' : '?';
            tip.innerHTML =
              '<div><b>' + _esc(hit.q_chr) + '</b> ' +
              _fmtBp(hit.q_start) + '–' + _fmtBp(hit.q_end) +
              ' <span style="color:' + t.inkDim + ';">×</span> ' +
              '<b>' + _esc(hit.t_chr) + '</b> ' +
              _fmtBp(hit.t_start) + '–' + _fmtBp(hit.t_end) + '</div>' +
              '<div style="color:' + t.inkDim + ';">' +
              'strand ' + _esc(hit.strand) + ' · PI ' + piStr +
              ' · block ' + _fmtBp(hit.q_end - hit.q_start) + '</div>';
            tip.style.display = 'block';
            const tipW = tip.offsetWidth || 200;
            const tipH = tip.offsetHeight || 40;
            let tx = mx + 12, ty = my + 12;
            if (tx + tipW > enlargedSize) tx = mx - tipW - 12;
            if (ty + tipH > enlargedSize) ty = my - tipH - 12;
            tip.style.left = tx + 'px';
            tip.style.top  = ty + 'px';
          } else {
            tip.style.display = 'none';
          }
        });
        plotHolder.addEventListener('mouseleave', () => {
          tip.style.display = 'none';
        });
        // Optional: alignment click handler
        if (typeof opts.onAlignmentClick === 'function') {
          plotHolder.addEventListener('click', (ev) => {
            const rect = plotHolder.getBoundingClientRect();
            const mx = ev.clientX - rect.left;
            const my = ev.clientY - rect.top;
            const v2 = computeLayoutFor(resolution);
            const hit = _findNearestAlignment(v2.alignments, v2.layout, data,
              enlargedSize, mx, my);
            if (hit) opts.onAlignmentClick(hit, ev);
          });
        }
        // Esc to close (defensive safety hatch — also closes if user
        // tabs focus away or otherwise loses the hover trail)
        const onEsc = (ev) => { if (ev.key === 'Escape') closePopup(); };
        document.addEventListener('keydown', onEsc);
        // Stash so close removes
        popupState._onEsc = onEsc;
      }
      // Render the SVG
      popupState.plotEl.innerHTML = '';                 // clear (preserves tip — re-add)
      popupState.plotEl.appendChild(popupState.tipEl);
      const svgString = _renderDotplotSVG(data, v.alignments, v.layout, enlargedSize, opts);
      popupState.plotEl.insertAdjacentHTML('afterbegin', svgString);
      // Update caption
      popupState.captionEl.textContent =
        v.alignments.length + ' alignments · ' +
        v.layout.qOrder.length + ' query × ' + v.layout.tOrder.length +
        ' target chroms · * = target reverse-flipped';
      popupState.open = true;
    };

    // Hover wiring on the mini panel
    let hoverTimer = null;
    miniHolder.addEventListener('mouseenter', () => {
      // Open immediately on hover (no pin step — pure hover model).
      if (!popupState.open) renderEnlarged();
    });
    miniHolder.addEventListener('mouseleave', () => {
      // Small grace period so the user can move into the popup itself
      // without the popup closing under them.
      if (hoverTimer) clearTimeout(hoverTimer);
      hoverTimer = setTimeout(() => closePopup(), 80);
    });
    // If the user moves into the popup, cancel the close timer; if they
    // move out of both mini and popup, schedule close.
    document.addEventListener('mouseover', (ev) => {
      if (!popupState.open) return;
      if (popupState.el && popupState.el.contains(ev.target)) {
        if (hoverTimer) { clearTimeout(hoverTimer); hoverTimer = null; }
      } else if (!miniHolder.contains(ev.target)) {
        if (hoverTimer) clearTimeout(hoverTimer);
        hoverTimer = setTimeout(() => closePopup(), 80);
      }
    });

    return wrap;
  }

  // ---------------------------------------------------------------------------
  // Public API
  // ---------------------------------------------------------------------------
  return {
    makeDotplotPanel: makeDotplotPanel,
    colormap:         colormap,
    layoutChroms:     layoutChroms,
    // For tests
    _internal: {
      cumulOffsets:           _cumulOffsets,
      pickResolution:         _pickResolution,
      resolutionOptions:      _resolutionOptions,
      findNearestAlignment:   _findNearestAlignment,
      RAMP_STOPS:             RAMP_STOPS,
    },
  };
}));
