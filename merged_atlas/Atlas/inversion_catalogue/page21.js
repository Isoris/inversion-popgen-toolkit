// Atlas/inversion_catalogue/page21.js
// =============================================================================
// page21 — Annotation cockpit (per-sample-lines canvas with cursor-driven
// candidate selection)
// (`<div id="page21">` contains canvas `annoCockpitCanvas` + footer panels)
//
// NOTE on handoff doc: HANDOFF_BATCH_3.md describes page21 as "Manual
// karyotype groups list" — that label is incorrect. The HTML region
// (legacy lines 7822–7859) is the annotation cockpit. Line ranges are
// authoritative; this extraction follows the line ranges.
//
// Source: legacy/Inversion_atlas.html lines 46938–47616 (JS body)
//                                    lines  7822–7859  (HTML shell)
// Public entry: refreshAnnotationCockpit()  (legacy line 47590)
//
// External dependencies:
//   TODO_MISSING(_gatherActiveCandidatesForInheritance)
//                                — legacy line 41196; gathers the active
//                                  candidate list for the current chrom.
//                                  Likely promoted to shared/ at merge.
//   TODO_MISSING(_wireCandidateHaplotypeAnnotations)
//                                — wires the band-annotation UI under
//                                  the canvas. Probably owned by Batch 1
//                                  (page2 candidate focus).
//   TODO_MISSING(candidateHaplotypeAnnotationsHtml)
//                                — HTML builder for the haplotype-
//                                  annotation panel (sibling of above).
//   TODO_MISSING(computeTrackedLinkageProjection)
//                                — legacy line ~46912; computes per-
//                                  candidate band purity for tracked
//                                  fish. Likely shared/ candidate.
//
//   global `state`               — reads:
//                                    state.data
//                                    state.tracked
//                                    state.cockpitCursor (initialized here)
// =============================================================================

const state = (typeof window !== 'undefined' && window.state) ? window.state : {};

// ---------------------------------------------------------------------------
// VERBATIM extraction from legacy lines 46970–47616.
// (The legacy `window.foo = foo` re-export block at the very end is replaced
// with ES-module exports below; otherwise the body is unchanged.)
// ---------------------------------------------------------------------------

const _ACK_PAD = { l: 60, r: 30, t: 30, b: 50 };
const _ACK_LINE_ALPHA = 0.35;
const _ACK_HIGHLIGHT_ALPHA = 0.85;
const _ACK_CURSOR_COLOR = '#f5a524';
const _ACK_BAND_PALETTE = ['#4fa3ff','#b8b8b8','#f5a524','#3cc08a','#e0555c','#b07cf7'];

function _ackBandColor(k) {
  if (k == null || k < 0) return '#666';
  return _ACK_BAND_PALETTE[k % _ACK_BAND_PALETTE.length];
}

function _ackEnsureState() {
  const _state = (typeof window !== 'undefined' && window.state) ? window.state : state;
  if (!_state.cockpitCursor) {
    _state.cockpitCursor = { mb: null, candidate_id: null };
  }
  return _state;
}

// Find the candidate under a given mb position. Returns the candidate
// from the items list (sorted by start_bp) or null.
function _annoCockpitCandidateAtCursor(mb, items) {
  if (mb == null || !Array.isArray(items)) return null;
  for (const it of items) {
    const lo = it.start_bp / 1e6;
    const hi = it.end_bp / 1e6;
    if (mb >= lo && mb <= hi) return it;
  }
  return null;
}

// Get the chromosome's mb extent — try state.data.chrom_len_bp if present,
// else fall back to max end_bp of items + a small margin. If neither is
// available, return null.
function _annoCockpitChromExtent(items) {
  const _state = (typeof window !== 'undefined' && window.state) ? window.state : state;
  const data = _state && _state.data;
  if (data && typeof data.chrom_len_bp === 'number' && data.chrom_len_bp > 0) {
    return { mbMin: 0, mbMax: data.chrom_len_bp / 1e6 };
  }
  // Use windows if present
  if (data && Array.isArray(data.windows) && data.windows.length > 0) {
    let lastEnd = null;
    for (let i = data.windows.length - 1; i >= 0; i--) {
      const w = data.windows[i];
      if (w && typeof w.end_bp === 'number') { lastEnd = w.end_bp; break; }
    }
    if (lastEnd != null) return { mbMin: 0, mbMax: lastEnd / 1e6 };
  }
  // Fallback to candidate range
  if (Array.isArray(items) && items.length > 0) {
    let maxBp = 0;
    for (const it of items) if (it.end_bp > maxBp) maxBp = it.end_bp;
    return { mbMin: 0, mbMax: maxBp / 1e6 * 1.05 };
  }
  return null;
}

function _annoCockpitDraw() {
  const canvas = document.getElementById('annoCockpitCanvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;
  const _state = _ackEnsureState();
  const items = _gatherActiveCandidatesForInheritance();
  const extent = _annoCockpitChromExtent(items);

  // Adjust canvas resolution to match its CSS-displayed size
  // (so trajectories render at the visible width)
  const dpr = (typeof window !== 'undefined' && window.devicePixelRatio) || 1;
  const cssW = canvas.clientWidth || 1200;
  const cssH = canvas.clientHeight || 500;
  if (canvas.width !== Math.floor(cssW * dpr)) canvas.width = Math.floor(cssW * dpr);
  if (canvas.height !== Math.floor(cssH * dpr)) canvas.height = Math.floor(cssH * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  // Clear
  ctx.fillStyle = '#1c2231';
  ctx.fillRect(0, 0, cssW, cssH);

  if (!extent || items.length === 0) {
    ctx.fillStyle = '#7a8398';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('no candidates on this chromosome yet', cssW / 2, cssH / 2);
    return;
  }

  const pad = _ACK_PAD;
  const plotW = cssW - pad.l - pad.r;
  const plotH = cssH - pad.t - pad.b;
  const mbMin = extent.mbMin;
  const mbMax = extent.mbMax;
  const mbToX = mb => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;

  // Find candidate under cursor
  const cursor = _state.cockpitCursor;
  const activeItem = _annoCockpitCandidateAtCursor(cursor.mb, items);
  if (activeItem) cursor.candidate_id = activeItem.id;
  else cursor.candidate_id = null;

  // 1. Background candidate-region shading. Each candidate gets a faint
  //    rectangle. The candidate under the cursor gets a brighter one.
  for (const it of items) {
    const xLo = mbToX(it.start_bp / 1e6);
    const xHi = mbToX(it.end_bp / 1e6);
    const w = xHi - xLo;
    const isActive = activeItem && it.id === activeItem.id;
    ctx.fillStyle = isActive
      ? 'rgba(245, 165, 36, 0.18)'
      : 'rgba(168, 177, 196, 0.06)';
    ctx.fillRect(xLo, pad.t, w, plotH);
    // Outline
    ctx.strokeStyle = isActive ? 'rgba(245, 165, 36, 0.55)' : 'rgba(168, 177, 196, 0.20)';
    ctx.lineWidth = isActive ? 1.5 : 0.5;
    ctx.strokeRect(xLo + 0.5, pad.t + 0.5, w - 1, plotH - 1);
    // Label "I#" centered above the rect
    ctx.fillStyle = isActive ? '#f5a524' : '#7a8398';
    ctx.font = isActive ? 'bold 10px sans-serif' : '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.fillText(`I${it.seq_num}`, (xLo + xHi) / 2, pad.t - 2);
  }

  // turn 2l: linkage shading from state.tracked. When the user has
  // selected a band (digit key) or lassoed on another page, draw alpha-
  // shaded intervals per candidate keyed to its inheritance group of
  // the dominant band among the lassoed fish. Same logic as turn 2i's
  // _drawTrackedLinkageStrip, but tuned to cockpit geometry.
  const trackedFish = _state && _state.tracked;
  if (Array.isArray(trackedFish) && trackedFish.length >= 3
      && typeof computeTrackedLinkageProjection === 'function') {
    const proj = computeTrackedLinkageProjection(trackedFish);
    if (proj && Array.isArray(proj.per_candidate)) {
      for (const rec of proj.per_candidate) {
        if (rec.purity < 0.30) continue;
        const xLo = mbToX(rec.start_bp / 1e6);
        const xHi = mbToX(rec.end_bp / 1e6);
        const ww = xHi - xLo;
        if (ww < 1) continue;
        const alpha = rec.purity * rec.purity * 0.50;
        const hex = rec.inh_group_color || '#7a8398';
        const rr = parseInt(hex.slice(1, 3), 16);
        const gg = parseInt(hex.slice(3, 5), 16);
        const bb = parseInt(hex.slice(5, 7), 16);
        ctx.fillStyle = `rgba(${rr}, ${gg}, ${bb}, ${alpha.toFixed(3)})`;
        ctx.fillRect(xLo, pad.t, ww, plotH);
        // Purity label inside the rect at the bottom
        if (ww >= 36 && rec.purity >= 0.50) {
          ctx.font = 'bold 9px sans-serif';
          ctx.textAlign = 'center';
          ctx.textBaseline = 'bottom';
          ctx.fillStyle = `rgba(${rr}, ${gg}, ${bb}, 0.95)`;
          ctx.fillText(
            `b${rec.dominant_band} · ${(rec.purity * 100).toFixed(0)}%`,
            (xLo + xHi) / 2, pad.t + plotH - 3
          );
        }
      }
    }
  }

  // 2. Per-sample lines, colored by band-of-active-candidate.
  // We render trajectories as thin polylines across the whole chromosome.
  // Each fish takes its band assignment from activeItem; if no active
  // candidate, fish are all gray.
  const data = _state.data;
  const samples = (data && Array.isArray(data.samples)) ? data.samples : [];
  const N = samples.length;
  const windows = (data && Array.isArray(data.windows)) ? data.windows : [];
  if (windows.length > 0 && N > 0) {
    // Find PC1 range for y-scale
    let pc1Min = +Infinity, pc1Max = -Infinity;
    for (const w of windows) {
      if (w && w.pca && w.pca.pc1) {
        for (let s = 0; s < N; s++) {
          const v = w.pca.pc1[s];
          if (isFinite(v)) {
            if (v < pc1Min) pc1Min = v;
            if (v > pc1Max) pc1Max = v;
          }
        }
      }
    }
    if (!isFinite(pc1Min)) { pc1Min = -1; pc1Max = 1; }
    const yPadFrac = 0.08;
    const yRange = pc1Max - pc1Min;
    const yLo = pc1Min - yRange * yPadFrac;
    const yHi = pc1Max + yRange * yPadFrac;
    const pcToY = v => pad.t + plotH - ((v - yLo) / (yHi - yLo)) * plotH;

    // Build per-fish color (band of active candidate)
    const fishColor = new Array(N).fill('rgba(140, 150, 165, 0.25)');
    if (activeItem && activeItem.labels) {
      for (let s = 0; s < N; s++) {
        const k = activeItem.labels[s];
        if (k >= 0) fishColor[s] = _ackBandColor(k);
      }
    }

    // For efficiency, draw all fish per band as a single path
    if (activeItem) {
      const K = activeItem.K;
      const byBand = Array.from({ length: K }, () => []);
      for (let s = 0; s < N; s++) {
        const k = activeItem.labels[s];
        if (k >= 0) byBand[k].push(s);
      }
      for (let k = 0; k < K; k++) {
        ctx.strokeStyle = _ackBandColor(k);
        ctx.globalAlpha = _ACK_LINE_ALPHA;
        ctx.lineWidth = 0.8;
        for (const s of byBand[k]) {
          ctx.beginPath();
          let started = false;
          for (let wi = 0; wi < windows.length; wi++) {
            const wn = windows[wi];
            if (!wn || !wn.pca || !wn.pca.pc1) continue;
            const v = wn.pca.pc1[s];
            if (!isFinite(v)) { started = false; continue; }
            // mid-bp of the window — fall back to wi proportional if no bp data
            let mb;
            if (typeof wn.start_bp === 'number' && typeof wn.end_bp === 'number') {
              mb = (wn.start_bp + wn.end_bp) / 2 / 1e6;
            } else {
              mb = mbMin + (wi / Math.max(1, windows.length - 1)) * (mbMax - mbMin);
            }
            const x = mbToX(mb);
            const y = pcToY(v);
            if (!started) { ctx.moveTo(x, y); started = true; } else { ctx.lineTo(x, y); }
          }
          ctx.stroke();
        }
      }
      ctx.globalAlpha = 1;
    } else {
      // No active candidate — render all fish in flat gray
      ctx.strokeStyle = 'rgba(140, 150, 165, 0.25)';
      ctx.lineWidth = 0.6;
      for (let s = 0; s < N; s++) {
        ctx.beginPath();
        let started = false;
        for (let wi = 0; wi < windows.length; wi++) {
          const wn = windows[wi];
          if (!wn || !wn.pca || !wn.pca.pc1) continue;
          const v = wn.pca.pc1[s];
          if (!isFinite(v)) { started = false; continue; }
          let mb;
          if (typeof wn.start_bp === 'number' && typeof wn.end_bp === 'number') {
            mb = (wn.start_bp + wn.end_bp) / 2 / 1e6;
          } else {
            mb = mbMin + (wi / Math.max(1, windows.length - 1)) * (mbMax - mbMin);
          }
          const x = mbToX(mb);
          const y = pcToY(v);
          if (!started) { ctx.moveTo(x, y); started = true; } else { ctx.lineTo(x, y); }
        }
        ctx.stroke();
      }
    }
  } else {
    // No window data — show placeholder
    ctx.fillStyle = '#7a8398';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('per-window PC1 data not loaded — only candidate regions shown', cssW / 2, cssH / 2);
  }

  // 3. Cursor line
  if (cursor.mb != null) {
    const cx = mbToX(cursor.mb);
    if (cx >= pad.l && cx <= pad.l + plotW) {
      ctx.strokeStyle = _ACK_CURSOR_COLOR;
      ctx.globalAlpha = _ACK_HIGHLIGHT_ALPHA;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(cx, pad.t);
      ctx.lineTo(cx, pad.t + plotH);
      ctx.stroke();
      ctx.globalAlpha = 1;
      // mb label above
      ctx.fillStyle = _ACK_CURSOR_COLOR;
      ctx.font = 'bold 10px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText(`${cursor.mb.toFixed(2)} Mb`, cx, pad.t - 14);
    }
  }

  // 4. X-axis: bp/Mb tick marks
  ctx.strokeStyle = 'rgba(168, 177, 196, 0.3)';
  ctx.lineWidth = 0.5;
  ctx.fillStyle = '#a8b1c4';
  ctx.font = '9px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const span = mbMax - mbMin;
  let step = 1;
  if (span > 200) step = 20;
  else if (span > 100) step = 10;
  else if (span > 50) step = 5;
  else if (span > 20) step = 2;
  for (let mb = Math.ceil(mbMin / step) * step; mb <= mbMax; mb += step) {
    const x = mbToX(mb);
    ctx.beginPath();
    ctx.moveTo(x, pad.t + plotH);
    ctx.lineTo(x, pad.t + plotH + 4);
    ctx.stroke();
    ctx.fillText(`${mb}`, x, pad.t + plotH + 6);
  }
  ctx.fillText('Mb', pad.l + plotW / 2, pad.t + plotH + 26);
}

function _annoCockpitUpdateReadouts() {
  const _state = _ackEnsureState();
  const items = _gatherActiveCandidatesForInheritance();
  const cursor = _state.cockpitCursor;
  const cohortReadout = document.getElementById('annoCohortReadout');
  if (cohortReadout) cohortReadout.textContent = `${items.length} candidate${items.length === 1 ? '' : 's'}`;
  const cursorReadout = document.getElementById('annoCursorReadout');
  if (cursorReadout) {
    if (cursor.mb == null) {
      cursorReadout.textContent = 'cursor: — Mb · no candidate';
    } else {
      const at = _annoCockpitCandidateAtCursor(cursor.mb, items);
      const candText = at
        ? `<b style="color:#f5a524;">I${at.seq_num}</b> · ${at.id} · ${(at.start_bp/1e6).toFixed(2)}–${(at.end_bp/1e6).toFixed(2)} Mb · K=${at.K}`
        : '<i>(between candidates)</i>';
      cursorReadout.innerHTML = `cursor: <b>${cursor.mb.toFixed(2)} Mb</b> · ${candText}`;
    }
  }
  // turn 2l: the right footer pane is now the haplotype-annotation panel
  // (full vocab dropdowns + auto-fill + linkage-aware band buttons). Only
  // re-render when the cursor moves to a DIFFERENT candidate, so the user
  // doesn't lose focus while typing custom labels. Cursor moves within
  // the same candidate (or in a gap) are no-ops here.
  const at = (cursor.mb != null) ? _annoCockpitCandidateAtCursor(cursor.mb, items) : null;
  const newCandId = at ? at.id : null;
  if (newCandId !== _state._cockpitLastHapCandId) {
    _state._cockpitLastHapCandId = newCandId;
    _annoCockpitRenderHapPanel();
  }
}

function _annoCockpitOnKey(ev) {
  const _state = _ackEnsureState();
  const items = _gatherActiveCandidatesForInheritance();
  if (items.length === 0) return;
  const extent = _annoCockpitChromExtent(items);
  if (!extent) return;
  const cursor = _state.cockpitCursor;

  // Step size: small step = 1% of span; large step = jump to next candidate
  const span = extent.mbMax - extent.mbMin;
  const smallStep = span * 0.005;     // 0.5% of chromosome span

  if (ev.key === 'ArrowRight' || ev.key === 'ArrowLeft') {
    ev.preventDefault();
    const dir = ev.key === 'ArrowRight' ? 1 : -1;
    if (cursor.mb == null) {
      cursor.mb = (extent.mbMin + extent.mbMax) / 2;
    } else if (ev.shiftKey) {
      // Jump to next/prev candidate boundary
      const sortedBoundaries = [];
      for (const it of items) {
        sortedBoundaries.push(it.start_bp / 1e6);
        sortedBoundaries.push(it.end_bp / 1e6);
      }
      sortedBoundaries.sort((a, b) => a - b);
      let target = null;
      if (dir > 0) {
        for (const b of sortedBoundaries) if (b > cursor.mb + 1e-9) { target = b; break; }
      } else {
        for (let i = sortedBoundaries.length - 1; i >= 0; i--) {
          if (sortedBoundaries[i] < cursor.mb - 1e-9) { target = sortedBoundaries[i]; break; }
        }
      }
      if (target != null) cursor.mb = target;
    } else {
      cursor.mb = Math.max(extent.mbMin, Math.min(extent.mbMax, cursor.mb + dir * smallStep));
    }
    _annoCockpitDraw();
    _annoCockpitUpdateReadouts();
    return;
  }
  if (ev.key === 'Escape') {
    ev.preventDefault();
    cursor.mb = null;
    cursor.candidate_id = null;
    // turn 2l: ESC also clears tracked selection
    _state.tracked = [];
    _annoCockpitDraw();
    _annoCockpitUpdateReadouts();
    return;
  }

  // turn 2l: digit keys 0-9 select that band of the candidate under
  // the cursor and set state.tracked to its members. Triggers the
  // linkage shading on the canvas + populates the linkage table panel.
  if (/^[0-9]$/.test(ev.key)) {
    const at = _annoCockpitCandidateAtCursor(cursor.mb, items);
    if (!at) return;
    const k = parseInt(ev.key, 10);
    if (k >= at.K) return;
    ev.preventDefault();
    // Build members for band k
    const members = [];
    for (let s = 0; s < at.labels.length; s++) {
      if (at.labels[s] === k) members.push(s);
    }
    if (members.length === 0) return;
    _state.tracked = members;
    _state.cockpitSelectedBand = { candidate_id: at.id, band: k, n: members.length };
    _annoCockpitDraw();
    _annoCockpitUpdateReadouts();
    _annoCockpitRenderLinkagePanel();
    _annoCockpitRenderHapPanel();
    return;
  }
}

function _annoCockpitOnClick(ev) {
  const canvas = ev.currentTarget || document.getElementById('annoCockpitCanvas');
  if (!canvas) return;
  const _state = _ackEnsureState();
  const items = _gatherActiveCandidatesForInheritance();
  const extent = _annoCockpitChromExtent(items);
  if (!extent) return;
  const rect = canvas.getBoundingClientRect();
  const x = ev.clientX - rect.left;
  const cssW = rect.width;
  const pad = _ACK_PAD;
  const plotW = cssW - pad.l - pad.r;
  if (x < pad.l || x > pad.l + plotW) return;
  const frac = (x - pad.l) / plotW;
  const mb = extent.mbMin + frac * (extent.mbMax - extent.mbMin);
  _state.cockpitCursor.mb = mb;
  _annoCockpitDraw();
  _annoCockpitUpdateReadouts();
  canvas.focus();
}

// turn 2l: render the linkage panel (left footer pane). Reuses turn 2j's
// table format adapted for the cockpit context (no jump-to-page-1 button
// since we ARE the cockpit).
function _annoCockpitRenderLinkagePanel() {
  const panel = document.getElementById('annoLinkagePanel');
  if (!panel) return;
  const _state = _ackEnsureState();
  const tracked = _state.tracked;
  const sel = _state.cockpitSelectedBand;
  if (!Array.isArray(tracked) || tracked.length < 3) {
    panel.innerHTML = '<div style="color:var(--ink-dimmer);">No band selected. Move cursor onto a candidate, press <b>0–9</b> to select a band.</div>';
    return;
  }
  const proj = (typeof computeTrackedLinkageProjection === 'function')
    ? computeTrackedLinkageProjection(tracked) : null;
  if (!proj || proj.per_candidate.length === 0) {
    panel.innerHTML = '<div style="color:var(--ink-dimmer);">No projection data.</div>';
    return;
  }
  // Header
  let header = '';
  if (sel) {
    const swatch = _ackBandColor(sel.band);
    header = `
      <div style="font-size:11px;color:#e8edf6;margin-bottom:6px;">
        <span style="display:inline-block;width:10px;height:10px;background:${swatch};border-radius:1px;vertical-align:middle;margin-right:4px;"></span>
        <b>${sel.candidate_id}</b> · band ${sel.band} · ${sel.n} fish projected
      </div>
    `;
  } else {
    header = `
      <div style="font-size:11px;color:#e8edf6;margin-bottom:6px;">
        ${tracked.length} fish projected (lassoed)
      </div>
    `;
  }
  // Table
  let tbl = '<table style="border-collapse:collapse;font-family:var(--mono);font-size:10px;width:100%;">';
  tbl += '<thead><tr style="color:var(--ink-dimmer);border-bottom:1px solid var(--rule);">';
  tbl += '<th style="text-align:left;padding:2px 4px;">I#</th>';
  tbl += '<th style="text-align:left;padding:2px 4px;">span (Mb)</th>';
  tbl += '<th style="text-align:center;padding:2px 4px;">dom</th>';
  tbl += '<th style="text-align:right;padding:2px 4px;">purity</th>';
  tbl += '<th style="text-align:center;padding:2px 4px;">grp</th>';
  tbl += '</tr></thead><tbody>';
  const isSelf = (recId) => sel && String(recId) === String(sel.candidate_id);
  for (const r of proj.per_candidate) {
    const purityPct = (r.purity * 100).toFixed(0);
    const purityColor = r.purity >= 0.90 ? '#3cc08a'
                      : r.purity >= 0.70 ? '#88c45e'
                      : r.purity >= 0.50 ? '#f5a524'
                      : '#7a8398';
    const inhTag = r.inh_group_id != null
      ? `<span style="display:inline-block;padding:1px 4px;background:${r.inh_group_color};color:#1c2231;border-radius:2px;font-weight:bold;font-size:9px;">g${r.inh_group_id}</span>`
      : '<span style="color:var(--ink-dimmer);">—</span>';
    const rowBg = isSelf(r.candidate_id) ? 'rgba(245, 165, 36, 0.08)' : 'transparent';
    const swatch = _ackBandColor(r.dominant_band);
    tbl += `<tr style="border-bottom:1px solid var(--rule);background:${rowBg};">`;
    tbl += `<td style="padding:2px 4px;color:#a8b1c4;">I${r.seq_num}</td>`;
    tbl += `<td style="padding:2px 4px;color:var(--ink-dimmer);">${(r.start_bp/1e6).toFixed(1)}-${(r.end_bp/1e6).toFixed(1)}</td>`;
    tbl += `<td style="padding:2px 4px;text-align:center;"><span style="display:inline-block;width:8px;height:8px;background:${swatch};border-radius:1px;vertical-align:middle;"></span> b${r.dominant_band}</td>`;
    tbl += `<td style="padding:2px 4px;text-align:right;color:${purityColor};font-weight:bold;">${purityPct}%</td>`;
    tbl += `<td style="padding:2px 4px;text-align:center;">${inhTag}</td>`;
    tbl += '</tr>';
  }
  tbl += '</tbody></table>';
  // Clear button
  tbl += `
    <div style="margin-top:8px;">
      <button id="annoCockpitClearBtn"
              style="font-family:var(--mono);font-size:10px;padding:2px 6px;background:var(--panel-2);border:1px solid var(--rule);color:var(--ink);border-radius:2px;cursor:pointer;">
        clear (Esc)
      </button>
    </div>
  `;
  panel.innerHTML = header + tbl;
  const clearBtn = document.getElementById('annoCockpitClearBtn');
  if (clearBtn) clearBtn.addEventListener('click', () => {
    _state.tracked = [];
    _state.cockpitSelectedBand = null;
    _annoCockpitDraw();
    _annoCockpitUpdateReadouts();
    _annoCockpitRenderLinkagePanel();
    _annoCockpitRenderHapPanel();
  });
}

// turn 2l: render the haplotype-annotation panel (right footer pane).
// Uses turn 2f's candidateHaplotypeAnnotationsHtml + wires it into the
// cockpit's container. Updates whenever the cursor moves to a different
// candidate.
function _annoCockpitRenderHapPanel() {
  const panel = document.getElementById('annoCandidateInfo');
  if (!panel) return;
  const _state = _ackEnsureState();
  const items = _gatherActiveCandidatesForInheritance();
  const cursor = _state.cockpitCursor;
  const at = _annoCockpitCandidateAtCursor(cursor.mb, items);
  if (!at) {
    panel.innerHTML = '<div style="color:var(--ink-dimmer);">Move the cursor over a candidate to see its haplotype labels here.</div>';
    return;
  }
  // Need the actual candidate object from state (not the items projection,
  // which is a slim view). Look up by id.
  const isDetailed = _state.activeMode === 'detailed';
  const src = isDetailed ? _state.candidates_detailed : _state.candidates;
  const cand = src && src[at.id];
  if (!cand) {
    panel.innerHTML = '<div style="color:#e0555c;">candidate object missing</div>';
    return;
  }
  // Build a compact header with candidate ID + band sizes
  const counts = new Array(at.K).fill(0);
  for (let s = 0; s < at.labels.length; s++) {
    const k = at.labels[s];
    if (k >= 0 && k < at.K) counts[k]++;
  }
  let header = `<div style="font-size:11px;color:#e8edf6;margin-bottom:6px;">`;
  header += `<b>I${at.seq_num}</b> · ${at.id} · ${(at.start_bp/1e6).toFixed(2)}–${(at.end_bp/1e6).toFixed(2)} Mb · K=${at.K}</div>`;
  header += `<div style="display:flex;gap:6px;flex-wrap:wrap;margin-bottom:8px;">`;
  for (let k = 0; k < at.K; k++) {
    const swatch = _ackBandColor(k);
    const isSelected = _state.cockpitSelectedBand
      && _state.cockpitSelectedBand.candidate_id === at.id
      && _state.cockpitSelectedBand.band === k;
    const border = isSelected ? '1.5px solid #f5a524' : '1px solid var(--rule)';
    header += `<button class="anno-band-pick" data-band-idx="${k}"
       title="Press ${k} to select" style="display:inline-flex;align-items:center;gap:4px;font-family:var(--mono);font-size:10px;padding:2px 6px;background:var(--panel-2);border:${border};color:var(--ink);border-radius:2px;cursor:pointer;">`;
    header += `<span style="display:inline-block;width:8px;height:8px;background:${swatch};border-radius:1px;"></span>`;
    header += `${k}: ${counts[k]}</button>`;
  }
  header += `</div>`;
  // Reuse turn 2f's annotation HTML (vocab + dropdowns + auto-fill button)
  const hapHtml = (typeof candidateHaplotypeAnnotationsHtml === 'function')
    ? candidateHaplotypeAnnotationsHtml(cand) : '';
  panel.innerHTML = header + hapHtml;
  // Wire the hap panel using the new sectionEl parameter (turn 2l refactor)
  const section = panel.querySelector('#hapLabelsSection');
  if (section && typeof _wireCandidateHaplotypeAnnotations === 'function') {
    try { _wireCandidateHaplotypeAnnotations(cand, section); }
    catch (e) { console.warn('[annoCockpit hap]', e && e.message); }
  }
  // Wire clickable band buttons (synth digit-key behavior)
  panel.querySelectorAll('.anno-band-pick').forEach(btn => {
    btn.addEventListener('click', () => {
      const k = parseInt(btn.getAttribute('data-band-idx'), 10);
      if (Number.isNaN(k) || k >= at.K) return;
      const members = [];
      for (let s = 0; s < at.labels.length; s++) {
        if (at.labels[s] === k) members.push(s);
      }
      if (members.length === 0) return;
      _state.tracked = members;
      _state.cockpitSelectedBand = { candidate_id: at.id, band: k, n: members.length };
      _annoCockpitDraw();
      _annoCockpitUpdateReadouts();
      _annoCockpitRenderLinkagePanel();
      _annoCockpitRenderHapPanel();
    });
  });
}

let _annoCockpitInited = false;
function _annoCockpitInit() {
  if (_annoCockpitInited) return;
  const canvas = document.getElementById('annoCockpitCanvas');
  if (!canvas) return;
  canvas.addEventListener('click', _annoCockpitOnClick);
  canvas.addEventListener('keydown', _annoCockpitOnKey);
  _annoCockpitInited = true;
}

// Public API: refresh the cockpit. Called when:
//   - the user navigates to page21 (tab click)
//   - candidates added/removed
//   - chromosome changes
function refreshAnnotationCockpit() {
  _ackEnsureState();
  _annoCockpitInit();
  const items = _gatherActiveCandidatesForInheritance();
  const empty = document.getElementById('annotationCockpitEmpty');
  const body = document.getElementById('annotationCockpitBody');
  if (items.length === 0) {
    if (empty) empty.style.display = 'block';
    if (body) body.style.display = 'none';
    return;
  }
  if (empty) empty.style.display = 'none';
  if (body) body.style.display = 'block';
  _annoCockpitDraw();
  _annoCockpitUpdateReadouts();
}


// ---------------------------------------------------------------------------
// End verbatim extraction.
// (Original legacy code had `window.foo = foo` re-exports here; ES-module
//  exports below replace those.)
// ---------------------------------------------------------------------------

// Public entry — refresh the cockpit. Called when:
//   - the user navigates to page21 (tab click)
//   - candidates added/removed
//   - chromosome changes
export { refreshAnnotationCockpit };

// Internal helpers re-exported for the merge chat / tests.
export {
  _annoCockpitInit,
  _annoCockpitDraw,
  _annoCockpitOnKey,
  _annoCockpitOnClick,
  _annoCockpitUpdateReadouts,
  _annoCockpitCandidateAtCursor,
  _annoCockpitChromExtent,
  _annoCockpitRenderLinkagePanel,
  _annoCockpitRenderHapPanel,
  _ackBandColor,
  _ackEnsureState,
};

// Constants other modules may want to read.
export {
  _ACK_PAD,
  _ACK_LINE_ALPHA,
  _ACK_HIGHLIGHT_ALPHA,
  _ACK_CURSOR_COLOR,
  _ACK_BAND_PALETTE,
};

export const __MODULE_ID__ = 'inversion_catalogue/page21';
