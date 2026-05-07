// Atlas/inversion_discovery/page1.js
//
// Local PCA on dosage. Sim_mat heatmap, robust |Z|, per-sample lines,
// K-means PCA, L3 contingency. Per HANDOFF_BATCH_1 this is THE big page.
//
// Extracted LITERALLY from legacy/Inversion_atlas.html.
// Per HANDOFF_BATCH_1: state is now passed as first arg (was a global in legacy).
// All other unresolved references are flagged TODO_MISSING — the merge chat
// will resolve them once all 5 batches land.
//
// Entry points (in extraction order):
  // drawSim() — legacy lines 31331-31638
  // drawSimMini() — legacy lines 31650-31767
  // drawZ() — legacy lines 31839-32460
  // drawLinesPanel() — legacy lines 34894-35744
  // drawPCA() — legacy lines 35749-35949
  // drawAnchorStrip() — legacy lines 36147-36256
  // renderL3Panel() — legacy lines 48685-49193
  // renderL3PanelSlab() — legacy lines 49209-49478
  // renderL3PanelScaleStability() — legacy lines 12833-12912
  // updateWinLabel() — legacy lines 51738-51742
  // setCur() — legacy lines 51747-51792
  // autoPickRadial() — legacy lines 51862-51892
  // applyData() — legacy lines 54476-54690
  // onSimClick() — legacy lines 52067-52111
  // onZClick() — legacy lines 52112-52202
  // onPCAClick() — legacy lines 52203-52244
  // togglePlay() — legacy lines 52399-52417
  // cycleKAside() — legacy lines 56537-56565
  // renderTrackedList() — legacy lines 52004-52062
  // renderManualGroupsList() — legacy lines 48152-48191
  // buildLinesPanelCheckboxes() — legacy lines 33016-33117
  // buildLinesPanel() — legacy lines 33334-33610
  // buildTrackPanels() — legacy lines 32805-32848
  // drawTracks() — legacy lines 32860-32872
  // refreshLinesColorMode() — legacy lines 33165-33193
  // setLinesPanelCandidateBands() — legacy lines 34002-34007
//

import { contextFromState, clusterL2, ClusterCache } from '../shared/per_l2_cluster.js';
import { hetRateColor } from '../shared/het_rate.js';
import { alignLabels, hungarianChainProjection, concordanceMatrix } from '../shared/hungarian.js';
import { buildContingency, computeARI, computeNMI, cramersV } from '../shared/contingency.js';
import { kmeans1D, kmeans2D, silhouette1D, adaptiveK1D } from '../shared/kmeans.js';

// ---------------------------------------------------------------------------
// TODO_MISSING — functions referenced but not defined in this module or shared/*.
// (These are called by the extracted bodies below. The merge chat will either
//  pull them from another batch, or hoist them into shared/ if they are reused.)
// ---------------------------------------------------------------------------
// TODO_MISSING(_assignCandidateLanes)
// TODO_MISSING(_buildJumpMask)
// TODO_MISSING(_drawBandTraceStrip)
// TODO_MISSING(_drawDiamondOverlay)
// TODO_MISSING(_drawInheritanceLabelsStrip)
// TODO_MISSING(_drawLineageStrip)
// TODO_MISSING(_drawRegimeBreadthStrip)
// TODO_MISSING(_drawSnpDensityShade)
// TODO_MISSING(_drawSnpDensityStrip)
// TODO_MISSING(_drawTrackedLinkageStrip)
// TODO_MISSING(_drawTransitionRateStrip)
// TODO_MISSING(_drawWRow)
// TODO_MISSING(_drawWinNavLane)
// TODO_MISSING(_ensureCsOverlayIndex)
// TODO_MISSING(_paintCandidateBands)
// TODO_MISSING(_refreshScreeInset)
// TODO_MISSING(_resolveSampleScopeColor)
// TODO_MISSING(_vColor)
// TODO_MISSING(_wRowBand)
// TODO_MISSING(_winNavBand)
// TODO_MISSING(allSampleIdx)
// TODO_MISSING(currentMbRange)
// TODO_MISSING(dataset)
// TODO_MISSING(drawCandidateBar)
// TODO_MISSING(drawRect)
// TODO_MISSING(escapeHtml)
// TODO_MISSING(families)
// TODO_MISSING(fitCanvas)
// TODO_MISSING(flushRun)
// TODO_MISSING(formatTrackVal)
// TODO_MISSING(getActiveSimScale)
// TODO_MISSING(getL2Cluster)
// TODO_MISSING(getLinesGrid)
// TODO_MISSING(getLinesSignAt)
// TODO_MISSING(getLinesValuesAt)
// TODO_MISSING(getPCRender)
// TODO_MISSING(getSampleColor)
// TODO_MISSING(hubs)
// TODO_MISSING(jittered)
// TODO_MISSING(layer)
// TODO_MISSING(mbAt)
// TODO_MISSING(niceTicks)
// TODO_MISSING(recomputeAnchorConcord)
// TODO_MISSING(samples)
// TODO_MISSING(simColor)
// TODO_MISSING(simColorPDF)
// TODO_MISSING(strokeSamplePath)
// TODO_MISSING(strokeSamplePathStyled)
// TODO_MISSING(themeColor)
// TODO_MISSING(toPx)
// TODO_MISSING(toPy)
// TODO_MISSING(toX)
// TODO_MISSING(toY)
// TODO_MISSING(trackedColor)
// TODO_MISSING(withAlpha)
// TODO_MISSING(xOfWin)
// TODO_MISSING(zColorPDF)

// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT — state.<slot> references that are NOT in
// shared/state.js SLOT_REGISTRY. Some may be ad-hoc geometry caches the
// legacy added imperatively (state._simGeom, state._zGeom, etc.); some may
// be slot names we need to register. Merge chat decides per-slot.
// ---------------------------------------------------------------------------
// TODO_MISSING_SLOT(state._l3CacheFp)
// TODO_MISSING_SLOT(state._l3CacheRendered)
// TODO_MISSING_SLOT(state._lineageComputeScheduled)
// TODO_MISSING_SLOT(state._simGeom)
// TODO_MISSING_SLOT(state._simMinimapGeom)
// TODO_MISSING_SLOT(state.ancestryPalette)
// TODO_MISSING_SLOT(state.cacheKey)
// TODO_MISSING_SLOT(state.candidateMode)
// TODO_MISSING_SLOT(state.colorMode)
// TODO_MISSING_SLOT(state.compareUnit)
// TODO_MISSING_SLOT(state.crossSpecies)
// TODO_MISSING_SLOT(state.hubFamilies)
// TODO_MISSING_SLOT(state.l2GroupCache)
// TODO_MISSING_SLOT(state.l3Draft)
// TODO_MISSING_SLOT(state.l3Mode)
// TODO_MISSING_SLOT(state.l3ReclusterMode)
// TODO_MISSING_SLOT(state.l3SecondaryMetric)
// TODO_MISSING_SLOT(state.scaleStabilityPanes)
// TODO_MISSING_SLOT(state.schemaVersion)
// TODO_MISSING_SLOT(state.secondaryL2)
// TODO_MISSING_SLOT(state.simInMinimap)
// TODO_MISSING_SLOT(state.singletonFamilyIds)
// TODO_MISSING_SLOT(state.smallFamilyIds)
// TODO_MISSING_SLOT(state.trackedN)
// TODO_MISSING_SLOT(state.windowToL1)
// TODO_MISSING_SLOT(state.windowToL2)
// TODO_MISSING_SLOT(state.zColorMode)
// TODO_MISSING_SLOT(state.zHighlightThr)
// TODO_MISSING_SLOT(state.zValueMode)

// ---------------------------------------------------------------------------
// Extracted bodies
// ---------------------------------------------------------------------------

// --- drawSim() — legacy lines 31331-31638 ---
export function drawSim(state) {
  const canvas = document.getElementById('simCanvas');
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  if (!state.data) return;
  const d = state.data;

  // ---- Compute centered square geometry ----
  // The heatmap must be a perfect square so the diagonal reads at 45°
  // (otherwise visual structure inside L2 envelopes is unreadable).
  // Reserve a margin for axis text + a side legend strip.
  const padTop    = 24;
  const padBottom = 22;
  const padLeft   = 50;
  const padRight  = 130;   // leaves room for the legend strip
  const availW = Math.max(50, w - padLeft - padRight);
  const availH = Math.max(50, h - padTop - padBottom);
  const side   = Math.max(50, Math.min(availW, availH));
  const simX0  = padLeft + Math.floor((availW - side) / 2);
  const simY0  = padTop  + Math.floor((availH - side) / 2);
  const simX1  = simX0 + side;
  const simY1  = simY0 + side;
  // Persist geometry for the click handler
  state._simGeom = { x0: simX0, y0: simY0, x1: simX1, y1: simY1, side };

  const scale = getActiveSimScale();
  if (scale && scale.sim && scale.n > 0) {
    const N = scale.n;
    const sim = scale.sim;
    const zArr = scale.z;
    const usePdf = state.pdfStyle && zArr;
    const img = ctx.createImageData(N, N);

    if (usePdf) {
      const q_lo = scale.q_lo, q_hi = scale.q_hi, z_max = scale.z_max;
      for (let k = 0; k < sim.length; k++) {
        const i = Math.floor(k / N), j = k - i * N;
        let r, g, b;
        if (i === j) {
          r = 232; g = 197; b = 71;  // diagonal yellow
        } else if (j < i) {
          [r, g, b] = simColorPDF(sim[k], q_lo, q_hi);
        } else {
          [r, g, b] = zColorPDF(zArr[k], z_max);
        }
        img.data[k*4]=r; img.data[k*4+1]=g; img.data[k*4+2]=b; img.data[k*4+3]=255;
      }
    } else {
      let mn = Infinity, mx = -Infinity;
      for (let i = 0; i < sim.length; i++) {
        if (sim[i] < mn) mn = sim[i]; if (sim[i] > mx) mx = sim[i];
      }
      const rng = Math.max(1e-9, mx - mn);
      for (let i = 0; i < sim.length; i++) {
        const v = (sim[i] - mn) / rng;
        const [r, g, b] = simColor(v);
        img.data[i*4]=r; img.data[i*4+1]=g; img.data[i*4+2]=b; img.data[i*4+3]=255;
      }
    }
    const off = document.createElement('canvas');
    off.width = N; off.height = N;
    off.getContext('2d').putImageData(img, 0, 0);
    // v4 turn 72: adaptive smoothing. When the source matrix N is comparable
    // to or larger than the display side, keep crisp pixels (false). When the
    // source is significantly smaller (e.g. 200×200 thumbnail blit to 800×800
    // panel), enable high-quality smoothing so the upscale uses bilinear/cubic
    // interpolation instead of nearest-neighbor. Threshold: smooth when the
    // upscale ratio exceeds 1.5×.
    if (side / N > 1.5) {
      ctx.imageSmoothingEnabled = true;
      ctx.imageSmoothingQuality = 'high';
    } else {
      ctx.imageSmoothingEnabled = false;
    }
    ctx.drawImage(off, simX0, simY0, side, side);

    // Frame around the heatmap
    ctx.strokeStyle = 'rgba(120,128,140,0.45)';
    ctx.lineWidth = 1;
    ctx.strokeRect(simX0 + 0.5, simY0 + 0.5, side - 1, side - 1);
  } else {
    ctx.fillStyle = themeColor('panel-2');
    ctx.fillRect(simX0, simY0, side, side);
  }

  // Window-coord -> pixel coord (mapped INSIDE the square)
  const Nw = d.n_windows;
  const toPx = (wIdx) => simX0 + (wIdx + 0.5) * side / Nw;
  const toPy = (wIdx) => simY0 + (wIdx + 0.5) * side / Nw;

  // L1 envelope rectangles
  if (Array.isArray(d.l1_envelopes)) {
    const curL1 = state.windowToL1 ? state.windowToL1[state.cur] : -1;
    ctx.lineWidth = 1;
    d.l1_envelopes.forEach((e, i) => {
      const x0 = toPx(e._s0), x1 = toPx(e._e0);
      const y0 = toPy(e._s0), y1 = toPy(e._e0);
      const isCur = (i === curL1);
      ctx.strokeStyle = isCur ? 'rgba(0,66,255,1.0)' : 'rgba(0,66,255,0.45)';
      ctx.lineWidth = isCur ? 2 : 1;
      if (isCur) {
        ctx.fillStyle = 'rgba(0,66,255,0.10)';
        ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
      }
      ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
    });
  }

  // L2 envelope rectangles
  if (Array.isArray(d.l2_envelopes)) {
    const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
    d.l2_envelopes.forEach((e, i) => {
      const x0 = toPx(e._s0), x1 = toPx(e._e0);
      const y0 = toPy(e._s0), y1 = toPy(e._e0);
      const isCur = (i === curL2);
      const isPinned = (i === state.secondaryL2 && state.secondaryL2 != null && i !== curL2);
      if (isPinned) {
        // Magenta highlight for the 2nd pinned focal
        ctx.strokeStyle = 'rgba(232,121,249,1.0)';  // #e879f9
        ctx.lineWidth = 2.5;
        ctx.fillStyle = 'rgba(232,121,249,0.10)';
        ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
      } else {
        ctx.strokeStyle = isCur ? 'rgba(0,230,118,1.0)' : 'rgba(0,230,118,0.55)';
        ctx.lineWidth = isCur ? 2 : 1;
        if (isCur) {
          ctx.fillStyle = 'rgba(0,230,118,0.12)';
          ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
        }
      }
      ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
    });
  }

  // v3.99 turn 13 ask 2: confirmed candidates render as green-outlined zones
  // in the sim_mat, with their short ID label on top. The same dark-green
  // hue used for the catalogue confirmed-row indicator (rgba(60,192,138)).
  // Each confirmed candidate is drawn as an outlined box spanning its
  // bp interval, with the candidate ID rendered in dark green just above
  // the box. Provisional candidates do NOT get this treatment — only
  // candidates that the user has explicitly marked confirmed.
  if (Array.isArray(state.candidateList) && state.candidateList.length > 0) {
    ctx.save();
    for (const cand of state.candidateList) {
      if (!cand || !cand.confirmed) continue;
      // Map bp range to window indices (nearest-neighbor); guard against missing
      // start_bp / end_bp on malformed candidates.
      if (!isFinite(cand.start_bp) || !isFinite(cand.end_bp)) continue;
      let s0 = -1, e0 = -1;
      for (let wi = 0; wi < d.n_windows; wi++) {
        const w0 = d.windows[wi];
        if (!w0) continue;
        if (s0 < 0 && w0.end_bp >= cand.start_bp) s0 = wi;
        if (w0.start_bp <= cand.end_bp) e0 = wi;
      }
      if (s0 < 0 || e0 < s0) continue;
      const x0 = toPx(s0), x1 = toPx(e0);
      const y0 = toPy(s0), y1 = toPy(e0);
      // Outline + soft-green fill (tint matches catalogue confirmed-row alpha)
      ctx.strokeStyle = 'rgba(60,192,138,0.85)';
      ctx.fillStyle   = 'rgba(60,192,138,0.10)';
      ctx.lineWidth = 1.5;
      ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
      ctx.strokeRect(x0, y0, x1 - x0, y1 - y0);
      // ID label centered above the top edge. Short form strips chrom prefix
      // (e.g. "C_gar_LG28_d17L2_0001_03" → "d17L2_0001_03"); falls back to
      // full ID for non-conventional names.
      const fullId = String(cand.id || '?');
      const shortId = fullId.replace(/^[^_]+_[^_]+_[^_]+_/, '');
      // Position label just above the upper-left corner of the box. If the
      // box is at the very top of the matrix, render the label INSIDE the
      // top-left corner instead so it doesn't get clipped above the matrix.
      const labelY = (y0 - 4 < simY0 + 8) ? y0 + 11 : y0 - 4;
      ctx.font = 'bold 9.5px ui-monospace, monospace';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'alphabetic';
      // Background pill so the label is readable over heatmap data.
      const labelW = ctx.measureText(shortId).width;
      const padX = 3, padY = 2;
      ctx.fillStyle = 'rgba(60,192,138,0.18)';
      ctx.fillRect(x0 - padX, labelY - 9, labelW + padX * 2, 9 + padY);
      ctx.fillStyle = 'rgba(40,140,90,1.0)';
      ctx.fillText(shortId, x0, labelY);
    }
    ctx.restore();
  }

  // v4 turn 114b: cross-species breakpoint red-cross overlay.
  // For each cs-breakpoint that maps to a window in the loaded chromosome,
  // draw a red cross (axis-aligned, ~14 px arms) at (toPx(wi), toPy(wi)) —
  // i.e., on the diagonal at the breakpoint's window position. This
  // marks where on the sim_mat the breakpoint sits, complementing the
  // page-16 catalogue. Drawn before the orange cursor crosshair so the
  // cursor stays visually on top.
  try {
    const csIdx = (typeof _ensureCsOverlayIndex === 'function')
      ? _ensureCsOverlayIndex() : null;
    if (csIdx && csIdx.bps.length > 0) {
      ctx.save();
      ctx.strokeStyle = 'rgba(232, 90, 90, 0.85)';
      ctx.lineWidth = 2;
      ctx.lineCap = 'round';
      const armPx = 7;  // half-arm length → 14 px total
      for (const e of csIdx.bps) {
        if (e.win == null || e.win < 0) continue;
        const cx = toPx(e.win);
        const cy = toPy(e.win);
        if (!Number.isFinite(cx) || !Number.isFinite(cy)) continue;
        // Clamp so cross stays visible even at edges
        if (cx < simX0 - armPx || cx > simX1 + armPx) continue;
        if (cy < simY0 - armPx || cy > simY1 + armPx) continue;
        ctx.beginPath();
        ctx.moveTo(cx - armPx, cy - armPx);
        ctx.lineTo(cx + armPx, cy + armPx);
        ctx.moveTo(cx + armPx, cy - armPx);
        ctx.lineTo(cx - armPx, cy + armPx);
        ctx.stroke();
      }
      ctx.restore();
    }
  } catch (_) { /* fail-soft */ }

  // Crosshair (clipped to inside the square)
  const xc = toPx(state.cur), yc = toPy(state.cur);
  ctx.strokeStyle = '#f5a524';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(Math.round(xc) + 0.5, simY0); ctx.lineTo(Math.round(xc) + 0.5, simY1);
  ctx.moveTo(simX0, Math.round(yc) + 0.5); ctx.lineTo(simX1, Math.round(yc) + 0.5);
  ctx.stroke();

  // Mb tick labels along the bottom (every ~5 Mb)
  const wins = d.windows;
  if (wins && wins.length > 0) {
    const mbMin = wins[0].center_mb, mbMax = wins[wins.length - 1].center_mb;
    const ticks = niceTicks(mbMin, mbMax, 6);
    ctx.fillStyle = themeColor('ink-dim');
    ctx.font = '10px ui-monospace, monospace';
    ctx.textAlign = 'center';
    for (const mb of ticks) {
      // Find the nearest window
      const t = (mb - mbMin) / Math.max(1e-9, mbMax - mbMin);
      const x = simX0 + t * side;
      ctx.fillText(mb.toFixed(0), x, simY1 + 14);
    }
    ctx.textAlign = 'right';
    ctx.fillText('Mb', simX1, simY1 + 14);
    // Y-axis label
    ctx.save();
    ctx.translate(simX0 - 8, simY0 + side / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText('window index', 0, 0);
    ctx.restore();
  }

  // ---- Side legend strip ----
  // Show: scale label, q_lo / q_hi, z_max, and small color ramps
  if (scale) {
    const lx = simX1 + 14;
    const lw = Math.max(30, Math.min(padRight - 20, w - lx - 8));
    let ly = simY0 + 4;
    ctx.fillStyle = themeColor('ink');
    ctx.font = 'bold 11px ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.fillText(scale.label || 'scale', lx, ly);
    ly += 16;

    if (state.pdfStyle && scale.z) {
      // Sim ramp (lower triangle palette)
      ctx.font = '10px ui-monospace, monospace';
      ctx.fillStyle = themeColor('ink-dim');
      ctx.fillText('similarity', lx, ly);
      ly += 4;
      const rampW = Math.min(lw, 90), rampH = 8;
      for (let i = 0; i < rampW; i++) {
        const v = i / (rampW - 1);
        const [r, g, b] = simColorPDF(v, scale.q_lo, scale.q_hi);
        ctx.fillStyle = `rgb(${r},${g},${b})`;
        ctx.fillRect(lx + i, ly, 1, rampH);
      }
      ly += rampH + 11;
      ctx.fillStyle = themeColor('ink-dim');
      ctx.fillText(`q05 ${scale.q_lo.toFixed(2)}`, lx, ly);
      ly += 11;
      ctx.fillText(`q95 ${scale.q_hi.toFixed(2)}`, lx, ly);
      ly += 16;

      // Z ramp (upper triangle palette)
      ctx.fillText('local |Z|', lx, ly);
      ly += 4;
      for (let i = 0; i < rampW; i++) {
        const t = i / (rampW - 1);
        const z = (t * 2 - 1) * scale.z_max;
        const [r, g, b] = zColorPDF(z, scale.z_max);
        ctx.fillStyle = `rgb(${r},${g},${b})`;
        ctx.fillRect(lx + i, ly, 1, rampH);
      }
      ly += rampH + 11;
      ctx.fillText(`±${scale.z_max.toFixed(1)}`, lx, ly);
      ly += 16;
    }
    ctx.fillStyle = themeColor('ink-dimmer');
    ctx.fillText(`${scale.n}×${scale.n}`, lx, ly);
    ly += 11;
    ctx.fillText(`thumb`, lx, ly);
  }
}

// --- drawSimMini() — legacy lines 31650-31767 ---
export function drawSimMini(state) {
  const canvas = document.getElementById('simMinimapCanvas');
  if (!canvas) return;
  const fit = fitCanvas(canvas);
  if (!fit) return;
  const { ctx, w, h } = fit;
  ctx.clearRect(0, 0, w, h);
  if (!state.data) return;
  const d = state.data;

  // Square heatmap centered in the available space. Smaller padding than
  // drawSim because the minimap is much smaller and we want to maximize
  // the heatmap area. No legend strip.
  const pad = 4;
  const availW = Math.max(20, w - 2 * pad);
  const availH = Math.max(20, h - 2 * pad);
  const side = Math.min(availW, availH);
  const x0 = pad + Math.floor((availW - side) / 2);
  const y0 = pad + Math.floor((availH - side) / 2);
  state._simMinimapGeom = { x0, y0, side };

  const scale = getActiveSimScale();
  if (scale && scale.sim && scale.n > 0) {
    const N = scale.n;
    const sim = scale.sim;
    const zArr = scale.z;
    const usePdf = state.pdfStyle && zArr;
    const img = ctx.createImageData(N, N);
    if (usePdf) {
      const q_lo = scale.q_lo, q_hi = scale.q_hi, z_max = scale.z_max;
      for (let k = 0; k < sim.length; k++) {
        const i = Math.floor(k / N), j = k - i * N;
        let r, g, b;
        if (i === j) { r = 232; g = 197; b = 71; }
        else if (j < i) { [r, g, b] = simColorPDF(sim[k], q_lo, q_hi); }
        else { [r, g, b] = zColorPDF(zArr[k], z_max); }
        img.data[k*4] = r; img.data[k*4+1] = g; img.data[k*4+2] = b; img.data[k*4+3] = 255;
      }
    } else {
      let mn = Infinity, mx = -Infinity;
      for (let i = 0; i < sim.length; i++) {
        if (sim[i] < mn) mn = sim[i]; if (sim[i] > mx) mx = sim[i];
      }
      const rng = Math.max(1e-9, mx - mn);
      for (let i = 0; i < sim.length; i++) {
        const v = (sim[i] - mn) / rng;
        const [r, g, b] = simColor(v);
        img.data[i*4] = r; img.data[i*4+1] = g; img.data[i*4+2] = b; img.data[i*4+3] = 255;
      }
    }
    const off = document.createElement('canvas');
    off.width = N; off.height = N;
    off.getContext('2d').putImageData(img, 0, 0);
    // v4 turn 72: adaptive smoothing — see main sim_mat draw above.
    if (side / N > 1.5) {
      ctx.imageSmoothingEnabled = true;
      ctx.imageSmoothingQuality = 'high';
    } else {
      ctx.imageSmoothingEnabled = false;
    }
    ctx.drawImage(off, x0, y0, side, side);
    ctx.strokeStyle = 'rgba(120,128,140,0.35)';
    ctx.lineWidth = 1;
    ctx.strokeRect(x0 + 0.5, y0 + 0.5, side - 1, side - 1);
  } else {
    ctx.fillStyle = themeColor('panel-2');
    ctx.fillRect(x0, y0, side, side);
  }

  const Nw = d.n_windows;
  const toPx = (wi) => x0 + (wi + 0.5) * side / Nw;
  const toPy = (wi) => y0 + (wi + 0.5) * side / Nw;

  // L1 envelopes — outlined only (no fill, since space is tight)
  if (Array.isArray(d.l1_envelopes)) {
    const curL1 = state.windowToL1 ? state.windowToL1[state.cur] : -1;
    d.l1_envelopes.forEach((e, i) => {
      const xa = toPx(e._s0), xb = toPx(e._e0);
      const ya = toPy(e._s0), yb = toPy(e._e0);
      const isCur = (i === curL1);
      ctx.strokeStyle = isCur ? 'rgba(0,80,255,0.95)' : 'rgba(0,80,255,0.35)';
      ctx.lineWidth = isCur ? 1.5 : 0.7;
      ctx.strokeRect(xa, ya, xb - xa, yb - ya);
    });
  }
  // L2 envelopes — green outlines
  if (Array.isArray(d.l2_envelopes)) {
    const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
    d.l2_envelopes.forEach((e, i) => {
      // _s0/_e0 are computed by buildEnvelopeIndices; defensively recompute
      // from start_w/end_w if missing (mini may render before that runs)
      const s0 = (e._s0 != null) ? e._s0 : Math.max(0, (e.start_w | 0) - 1);
      const e0 = (e._e0 != null) ? e._e0 : Math.max(s0, (e.end_w | 0) - 1);
      const xa = toPx(s0), xb = toPx(e0);
      const ya = toPy(s0), yb = toPy(e0);
      const isCur = (i === curL2);
      ctx.strokeStyle = isCur ? 'rgba(0,230,118,0.95)' : 'rgba(0,230,118,0.35)';
      ctx.lineWidth = isCur ? 1.5 : 0.6;
      ctx.strokeRect(xa, ya, xb - xa, yb - ya);
    });
  }

  // Current scrubber position — orange crosshair on diagonal
  const cur = state.cur;
  if (cur != null && cur >= 0 && cur < Nw) {
    const cx = toPx(cur), cy = toPy(cur);
    ctx.strokeStyle = 'rgba(245,165,36,0.95)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    // Vertical line through current column
    ctx.moveTo(cx, y0);
    ctx.lineTo(cx, y0 + side);
    // Horizontal line through current row
    ctx.moveTo(x0, cy);
    ctx.lineTo(x0 + side, cy);
    ctx.stroke();
  }
}

// --- drawZ() — legacy lines 31839-32460 ---
export function drawZ(state) {
  const canvas = document.getElementById('zCanvas');
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  if (!state.data) {
    // v3.99 turn 14e ask 2: when no data is loaded, draw a centered hint
    // text instead of a silent blank canvas. In compact mode the right
    // column was completely empty on first open with no clue what to do.
    // The hint mirrors the sidebar emptyState language so the user knows
    // both surfaces point to the same action.
    ctx.fillStyle = themeColor('ink-dim');
    ctx.font = '12px ui-monospace, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.fillText('Load a precomp JSON in the sidebar (Browse…) to begin.', w / 2, h / 2 - 8);
    ctx.font = '10.5px ui-monospace, monospace';
    ctx.fillStyle = themeColor('ink-dimmer') || themeColor('ink-dim');
    ctx.fillText('Robust |Z| profile renders here once data is loaded.', w / 2, h / 2 + 10);
    return;
  }
  const d = state.data;
  // v3.61: when collapsed, render a COMPACT view that keeps the L1/L2 zone
  // bars + boundary arrows + a yellow candidate strip visible. Previously
  // collapsing the Z panel hid everything (including the L1/L2 zone bars,
  // which the user relies on for orientation while scrubbing). Compact
  // mode is the new default for the collapse state.
  if (state.zCollapsed) {
    const pad = { l: 44, r: 16, t: 2, b: 1 };
    const zoneH = 14;
    const peaksH = 5;
    // v3.99 t14e+ continue: candidate bar above L1/L2 (was below boundaries
    // in old strip placement). Slightly thinner here than in non-collapsed
    // mode since vertical space is constrained; ★ glyph still fits at h=5.
    // v4 turn 13 (Deliverable B): when overlapping candidates require N>1
    // lanes, the bar's TOTAL height is N × candBarH so each lane gets the
    // same per-lane height as the single-lane case.
    const candBarH = 5;
    const candGap  = 1;
    const _candLanes = (typeof _assignCandidateLanes === 'function')
      ? _assignCandidateLanes(state.candidateList || []).n_lanes : 1;
    const candBarTotal = candBarH * _candLanes;
    const zoneTop  = pad.t + candBarTotal + candGap;
    const plotW = w - pad.l - pad.r;
    const _mbR = currentMbRange();
    const mbMin = _mbR.mbMin, mbMax = _mbR.mbMax;
    const toX = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;
    const xOfWin = (wi) => toX(d.windows[wi].center_mb);
    // v3.99 t14e+ continue: Candidate bar — pending grey + confirmed gold ★
    // sitting ABOVE the L1/L2 zone bars
    drawCandidateBar(ctx, d, toX, pad.t, candBarTotal);
    // v4 turn 1 ask 5: "cand" label to the left of the candidate strip,
    // mirroring the L1/L2 label pattern below. Quentin: "the candidate
    // track is missing L3 or Candidate label to the left." Using "cand"
    // (4 chars) so it fits the same pad.l offset as "L1"/"L2" without
    // overlapping the y-axis tick labels. Color is amber (--accent /
    // rgba(245,165,36)) — matches the candidate bar's confirmed gold
    // and is the recognizable "candidate" hue throughout the scrubber.
    ctx.fillStyle = 'rgba(245, 165, 36, 0.95)';
    ctx.font = '8px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText('cand', pad.l - 4, pad.t + candBarTotal / 2);
    // L1 zone bars (top half of zoneH, deeper blue — v3.99 turn 13 ask 3)
    if (Array.isArray(d.l1_envelopes)) {
      const y0 = zoneTop, hzone = zoneH / 2 - 1;
      const curL1 = state.windowToL1 ? state.windowToL1[state.cur] : -1;
      d.l1_envelopes.forEach((e, i) => {
        const x0 = xOfWin(e._s0), x1 = xOfWin(e._e0);
        ctx.fillStyle = (i === curL1) ? 'rgba(48,116,200,0.95)' : 'rgba(48,116,200,0.50)';
        ctx.fillRect(x0, y0, Math.max(1, x1 - x0), hzone);
      });
      // v4 turn 4 ask 2: dark-red 1px separator at each fragment's right edge
      // so adjacent L1 envelopes are visible. Without these, neighbors at
      // alpha=0.50 visually merge into one bar.
      ctx.fillStyle = 'rgba(140, 29, 36, 0.95)';
      d.l1_envelopes.forEach((e, i) => {
        if (i === d.l1_envelopes.length - 1) return;
        const x1 = xOfWin(e._e0);
        ctx.fillRect(x1 - 0.5, y0, 1, hzone);
      });
      ctx.fillStyle = 'rgba(48,116,200,0.95)';
      ctx.font = '9px ui-monospace, monospace';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText('L1', pad.l - 4, y0 + hzone / 2);
    }
    // L2 zone bars (bottom half, deeper green — v3.99 turn 13 ask 3)
    if (Array.isArray(d.l2_envelopes)) {
      const y0 = zoneTop + zoneH / 2 + 1, hzone = zoneH / 2 - 1;
      const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
      d.l2_envelopes.forEach((e, i) => {
        const x0 = xOfWin(e._s0), x1 = xOfWin(e._e0);
        ctx.fillStyle = (i === curL2) ? 'rgba(34,160,80,0.95)' : 'rgba(34,160,80,0.50)';
        ctx.fillRect(x0, y0, Math.max(1, x1 - x0), hzone);
      });
      // v4 turn 4 ask 2: dark-red 1px separator (same as L1)
      ctx.fillStyle = 'rgba(140, 29, 36, 0.95)';
      d.l2_envelopes.forEach((e, i) => {
        if (i === d.l2_envelopes.length - 1) return;
        const x1 = xOfWin(e._e0);
        ctx.fillRect(x1 - 0.5, y0, 1, hzone);
      });
      // v4 turn 107: window-mode candidate boundary ticks on minimap L2 row
      if (Array.isArray(state.candidateList)) {
        ctx.fillStyle = 'rgba(60, 223, 255, 0.95)';
        for (const c of state.candidateList) {
          if (!c || c.source !== 'l3_draft_w') continue;
          if (c.start_w != null) {
            const xs = xOfWin(c.start_w);
            ctx.fillRect(xs - 0.5, y0 - 1, 1, hzone + 2);
          }
          if (c.end_w != null) {
            const xe = xOfWin(c.end_w);
            ctx.fillRect(xe + 0.5, y0 - 1, 1, hzone + 2);
          }
        }
      }
      if (state.l3Draft && state.l3Draft.resolution === 'W' &&
          state.l3Draft.start_w != null && state.l3Draft.end_w != null) {
        ctx.fillStyle = 'rgba(245, 196, 58, 0.95)';   // gold (draft palette)
        const xs = xOfWin(state.l3Draft.start_w);
        const xe = xOfWin(state.l3Draft.end_w);
        ctx.fillRect(xs - 0.5, y0 - 1, 1, hzone + 2);
        ctx.fillRect(xe + 0.5, y0 - 1, 1, hzone + 2);
      }
      ctx.fillStyle = 'rgba(34,160,80,0.95)';
      ctx.font = '9px ui-monospace, monospace';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText('L2', pad.l - 4, y0 + hzone / 2);
    }
    // v4 turn 33 (Ask A): always-on windows-navigation lane, sits between
    // the L2 zone bar and the candidate-mode W-row. Visual rhythm is
    // cand → L1 → L2 → win-nav → W-row(if cand mode) → boundary arrows.
    const _navBandC = _winNavBand({ collapsed: true, zoneTop, zoneH });
    if (_navBandC) {
      _drawWinNavLane(ctx, d, toX, xOfWin, _navBandC, pad.l);
    }
    const _navExtraC = _navBandC ? (_navBandC.h + _navBandC.gap) : 0;
    // v4 turn 10: W-row (per-window candidate-draft strip), collapsed mode.
    // Sits directly below the L2 zone bar; only visible when in candidate mode.
    // v4 turn 33: when nav-lane is visible we pass an inflated zoneH so the
    // W-row band lands below it instead of stacking on top.
    const _wBandC = _wRowBand({ collapsed: true, zoneTop,
                                  zoneH: zoneH + _navExtraC });
    if (_wBandC) {
      _drawWRow(ctx, d, toX, xOfWin, _wBandC);
    }
    // Boundary arrows (just below zone bars)
    if (Array.isArray(d.l2_boundaries)) {
      // v4 turn 10: when W-row is visible, push the arrow band down by its
      // height + gap so it doesn't overlap the W-row.
      // v4 turn 33: also push down by the nav-lane's height + gap.
      const _wExtra = _wBandC ? (_wBandC.h + _wBandC.gap) : 0;
      const _vertExtra = _navExtraC + _wExtra;
      const yArrowTop  = zoneTop + zoneH + _vertExtra + 1;
      const yArrowBase = zoneTop + zoneH + _vertExtra + peaksH - 1;
      const halfW = 3;
      for (const b of d.l2_boundaries) {
        const wi = (b.boundary_w | 0) - 1;
        if (wi < 0 || wi >= d.n_windows) continue;
        const x = toX(d.windows[wi].center_mb);
        const col = STATUS_COLOR[b.validation_status] || '#666';
        ctx.fillStyle = col;
        ctx.beginPath();
        ctx.moveTo(x, yArrowTop);
        ctx.lineTo(x - halfW, yArrowBase);
        ctx.lineTo(x + halfW, yArrowBase);
        ctx.closePath();
        ctx.fill();
      }
    }
    // Cursor — extends from pad.t (top of candidate bar) through L1/L2 down
    // to the bottom of the boundary-arrow band.
    // v4 turn 10: include W-row in cursor extent when W-row is visible.
    // v4 turn 33: include nav-lane too.
    const _wExtraCur = _wBandC ? (_wBandC.h + _wBandC.gap) : 0;
    const _vertExtraCur = _navExtraC + _wExtraCur;
    const xCur = toX(d.windows[state.cur].center_mb);
    ctx.strokeStyle = '#f5a524';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(xCur + 0.5, pad.t);
    ctx.lineTo(xCur + 0.5, zoneTop + zoneH + _vertExtraCur + peaksH);
    ctx.stroke();
    return;
  }
  // v3.57: pad.t reduced from 18 → 6 because the panel-label header is now
  // in normal flow ABOVE the canvas (flex column on #zPanel), not absolute
  // overlay. Zone bars sit at canvas-top with just 6px breathing room.
  // v3.99 turn 14e: bumped 6 → 14 because turn 14d hid the panel-label
  // header strip entirely (moved explanation to a panel-level tooltip and
  // ⚙/▼ buttons to canvas bottom-left). With the strip gone, the canvas
  // is flush with the top of the panel and 6 px isn't enough breathing
  // room — the L1/L2 zone bars at y=pad.t looked like they were overflowing
  // into adjacent panels above. 14 px restores proper visual separation.
  const pad = { l: 44, r: 16, t: 14, b: 22 };
  // Reserve top strip for L1/L2 zone bars
  const zoneH = 14;
  // v3.61: reserve a thin band BELOW the zone bars for L2 boundary arrows
  // (previously rendered as filled circles at the BOTTOM of the plot area on
  // the dashed yellow line). Putting them just under the zone bars makes them
  // read as "boundary markers tied to the L1/L2 envelope coloring above"
  // rather than getting confused with the z-score scatter.
  const peaksH = 8;
  // v3.99 t14e+ continue: reserve a band ABOVE the L1/L2 zone bars for the
  // candidate bar (drawCandidateBar). Pending candidates render grey here,
  // confirmed candidates render shiny gold + ★. 2 px gap separates this
  // strip from the L1 bar below it. Total top-band cost: candBarH + 2 px.
  const candBarH = 7;
  const candGap  = 2;
  // v4 turn 13 (Deliverable B): when overlapping candidates require N>1
  // lanes, the bar's TOTAL height is N × candBarH so each lane gets the
  // same per-lane height as the single-lane case. zoneTop and plotH both
  // adjust to keep L1/L2/W-row/peaks/scatter visually in the same place
  // relative to the cand bar's BOTTOM.
  const _candLanes = (typeof _assignCandidateLanes === 'function')
    ? _assignCandidateLanes(state.candidateList || []).n_lanes : 1;
  const candBarTotal = candBarH * _candLanes;
  // Effective top of L1 bar shifts down by (candBarTotal + candGap)
  const zoneTop = pad.t + candBarTotal + candGap;
  // v4 turn 33 (Ask A): always-on nav-windows lane between L2 and W-row.
  // _navExtra captures its height + gap, and the W-row computation uses
  // an inflated zoneH so it lands BELOW the nav-lane.
  const _navBandEx = _winNavBand({ collapsed: false, zoneTop, zoneH });
  const _navExtra  = _navBandEx ? (_navBandEx.h + _navBandEx.gap) : 0;
  // v4 turn 10: reserve a per-window draft strip between the L2 zone bar and
  // the boundary-arrow band when in candidate mode. _wExtra captures both the
  // row height and its top gap; everything below it (boundary arrows, plot,
  // toY ladder) shifts down by this amount so the layout doesn't collide.
  // _wRowBand returns null when the row is not visible — in that case
  // _wExtra evaluates to 0 and the layout is identical to pre-v4-turn-10.
  // v4 turn 33: pass zoneH + _navExtra so W-row lands below the nav-lane.
  const _wRowEx = _wRowBand({ collapsed: false, zoneTop,
                                zoneH: zoneH + _navExtra });
  const _wExtra = _wRowEx ? (_wRowEx.h + _wRowEx.gap) : 0;
  // Total vertical reservation between L2 bottom and the boundary-arrow top.
  const _vertExtra = _navExtra + _wExtra;
  const plotW = w - pad.l - pad.r;
  const plotH = h - pad.t - pad.b - candBarTotal - candGap - zoneH - _vertExtra - peaksH;

  const mbs = d.windows.map(w0 => w0.center_mb);
  // v3.41: signed-or-absolute Z values per state.zValueMode
  const zMode_val = state.zValueMode || 'abs';
  const zs = (zMode_val === 'signed')
    ? d.windows.map(w0 => +(w0.z || 0))
    : d.windows.map(w0 => Math.abs(w0.z || 0));
  const _mbR = currentMbRange();
  const mbMin = _mbR.mbMin, mbMax = _mbR.mbMax;
  // Compute zMax / zMin from the visible window subset only.
  let zMax = 0, zMin = 0;
  for (let i = 0; i < d.n_windows; i++) {
    const cm = mbs[i];
    if (cm < mbMin || cm > mbMax) continue;
    if (zs[i] > zMax) zMax = zs[i];
    if (zMode_val === 'signed' && zs[i] < zMin) zMin = zs[i];
  }
  zMax = Math.max(3, Math.ceil(zMax));
  if (zMode_val === 'signed') {
    // Symmetric range so 0 sits in the middle. Use max(|zMax|, |zMin|).
    const r = Math.max(zMax, Math.abs(Math.floor(zMin)));
    zMax = r;
    zMin = -r;
  } else {
    zMin = 0;
  }

  const toX = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;
  // v3.61: z-scatter starts at zoneTop + zoneH + peaksH (below boundary-arrows)
  // v3.99 t14e+ continue: zoneTop = pad.t + candBarH + candGap (was pad.t)
  // v4 turn 10: + _wExtra so the W-row sits between the L2 bar and the arrows.
  // v4 turn 33: switched to _vertExtra (which includes the always-on nav-lane).
  const toY = (z)  => zoneTop + zoneH + _vertExtra + peaksH + plotH - ((z - zMin) / (zMax - zMin)) * plotH;
  const xOfWin = (wi) => toX(d.windows[wi].center_mb);

  // L1 zone bars (top half of zoneH, deeper blue per Quentin's request).
  // v3.99 turn 13 ask 3: colors darkened to be more legible against the
  // panel background; old 79,163,255 was washed-out and the inactive band
  // at 0.35 alpha was nearly invisible. New rgb(48,116,200) is the same hue
  // family but with more saturation. Active state stays brighter via alpha.
  if (Array.isArray(d.l1_envelopes)) {
    const y0 = zoneTop, hzone = zoneH / 2 - 1;
    const curL1 = state.windowToL1 ? state.windowToL1[state.cur] : -1;
    d.l1_envelopes.forEach((e, i) => {
      const x0 = xOfWin(e._s0), x1 = xOfWin(e._e0);
      ctx.fillStyle = (i === curL1) ? 'rgba(48,116,200,0.95)' : 'rgba(48,116,200,0.50)';
      ctx.fillRect(x0, y0, Math.max(1, x1 - x0), hzone);
    });
    // v4 turn 4 ask 2: dark-red 1px separator at each fragment's right edge
    // so adjacent L1 envelopes are visually distinct.
    ctx.fillStyle = 'rgba(140, 29, 36, 0.95)';
    d.l1_envelopes.forEach((e, i) => {
      if (i === d.l1_envelopes.length - 1) return;
      const x1 = xOfWin(e._e0);
      ctx.fillRect(x1 - 0.5, y0, 1, hzone);
    });
    // v3.99 turn 13 ask 3: "L1" label in the y-axis margin so users can
    // tell what the blue bar means. Right-aligned to pad.l - 4 so it sits
    // just inside the bar's left edge.
    ctx.fillStyle = 'rgba(48,116,200,0.95)';
    ctx.font = '9px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText('L1', pad.l - 4, y0 + hzone / 2);
  }
  // L2 zone bars (bottom half, deeper green)
  if (Array.isArray(d.l2_envelopes)) {
    const y0 = zoneTop + zoneH / 2 + 1, hzone = zoneH / 2 - 1;
    const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
    d.l2_envelopes.forEach((e, i) => {
      const x0 = xOfWin(e._s0), x1 = xOfWin(e._e0);
      ctx.fillStyle = (i === curL2) ? 'rgba(34,160,80,0.95)' : 'rgba(34,160,80,0.50)';
      ctx.fillRect(x0, y0, Math.max(1, x1 - x0), hzone);
    });
    // v4 turn 4 ask 2: dark-red 1px separator (same as L1)
    ctx.fillStyle = 'rgba(140, 29, 36, 0.95)';
    d.l2_envelopes.forEach((e, i) => {
      if (i === d.l2_envelopes.length - 1) return;
      const x1 = xOfWin(e._e0);
      ctx.fillRect(x1 - 0.5, y0, 1, hzone);
    });
    // v4 turn 107: bright cyan tick markers on the L2 strip at the start/end
    // window of every committed candidate whose boundaries don't align with
    // the underlying L2 envelopes. These are window-mode candidates from
    // turn 106 (source='l3_draft_w'). Visually communicates "user-defined
    // unit boundaries here that L2 doesn't reflect" — the L2 segments
    // themselves stay unchanged because the L2 registry isn't authoritative
    // anymore (Quentin's design choice for turn 106-107: L2 = navigation,
    // candidates = real units).
    if (Array.isArray(state.candidateList)) {
      ctx.fillStyle = 'rgba(60, 223, 255, 0.95)';
      for (const c of state.candidateList) {
        if (!c) continue;
        // Only draw ticks for window-mode candidates (turn-106 tagged) OR
        // candidates whose start_w / end_w don't fall on an L2 boundary.
        const isWindowMode = c.source === 'l3_draft_w';
        if (!isWindowMode) continue;
        if (c.start_w != null) {
          const xs = xOfWin(c.start_w);
          ctx.fillRect(xs - 1, y0 - 1, 2, hzone + 2);
        }
        if (c.end_w != null) {
          const xe = xOfWin(c.end_w);
          ctx.fillRect(xe, y0 - 1, 2, hzone + 2);
        }
      }
    }
    // Active draft in window mode — show its current ticks too (as the user
    // is extending/shrinking with arrows, the cyan ticks update live).
    if (state.l3Draft && state.l3Draft.resolution === 'W' &&
        state.l3Draft.start_w != null && state.l3Draft.end_w != null) {
      ctx.fillStyle = 'rgba(245, 196, 58, 0.95)';   // gold (matches draft palette)
      const xs = xOfWin(state.l3Draft.start_w);
      const xe = xOfWin(state.l3Draft.end_w);
      ctx.fillRect(xs - 1, y0 - 1, 2, hzone + 2);
      ctx.fillRect(xe,     y0 - 1, 2, hzone + 2);
    }
    // v3.99 turn 13 ask 3: "L2" label, same pattern as L1
    ctx.fillStyle = 'rgba(34,160,80,0.95)';
    ctx.font = '9px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText('L2', pad.l - 4, y0 + hzone / 2);
  }

  // Threshold lines + |Z| points
  ctx.strokeStyle = 'rgba(42,50,66,0.6)';
  for (let i = 0; i <= 4; i++) {
    // v3.61: gridlines start at top of z-scatter (after zone bars + peaks band)
    // v3.99 t14e+ continue: zoneTop replaces pad.t (candidate bar offset)
    // v4 turn 10: + _wExtra so gridlines stay aligned with toY when W-row visible.
    // v4 turn 33: switched to _vertExtra.
    const y = zoneTop + zoneH + _vertExtra + peaksH + (i / 4) * plotH;
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + plotW, y); ctx.stroke();
  }
  const ths = [
    { z: 2.5, color: 'rgba(224,85,92,0.45)' },
    { z: 1.8, color: 'rgba(60,192,138,0.40)' },
    { z: 1.2, color: 'rgba(245,165,36,0.40)' },
  ];
  ctx.setLineDash([4, 4]);
  for (const th of ths) {
    if (th.z > zMax) continue;
    ctx.strokeStyle = th.color;
    let y = toY(th.z);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + plotW, y); ctx.stroke();
    // In signed mode, also draw the negative mirror so deviations in either
    // direction are flagged
    if (zMode_val === 'signed' && -th.z >= zMin) {
      y = toY(-th.z);
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + plotW, y); ctx.stroke();
    }
  }
  ctx.setLineDash([]);
  // Solid baseline at zero (signed mode) — important visual reference
  if (zMode_val === 'signed') {
    ctx.strokeStyle = themeColor('ink-dim');
    ctx.lineWidth = 1;
    const y0 = toY(0);
    ctx.beginPath();
    ctx.moveTo(pad.l, y0); ctx.lineTo(pad.l + plotW, y0);
    ctx.stroke();
  }
  // v3.50: y-axis tick labels are now generated by niceTicks() below in the
  // "Axes labels" block. The previous implementation hard-coded zMin/mid/zMax,
  // and was then overwritten by an integer ladder loop. Both paths consolidated.
  // Y-axis label (mode indicator) at top-left
  ctx.fillStyle = themeColor('ink-dim');
  ctx.font = '9px ui-monospace, monospace';
  ctx.textAlign = 'left';
  ctx.fillText(zMode_val === 'signed' ? 'Z' : '|Z|', 4, zoneTop + zoneH + _vertExtra + peaksH + 9);

  // ---------------------------------------------------------------------------
  // Z dot coloring — switchable mode (state.zColorMode):
  //   'bands'     band thresholds (current default; preserves existing colors)
  //   'gradient'  continuous viridis-ish from cool low |Z| to hot high |Z|
  //   'zone'      bright if window is in current focal L2; dim otherwise
  //   'highlight' bright if |Z| >= state.zHighlightThr; dim below
  // ---------------------------------------------------------------------------
  const zMode = state.zColorMode || 'bands';
  const curL1 = state.windowToL1 ? state.windowToL1[state.cur] : -1;
  const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
  const zThr  = (state.zHighlightThr != null && isFinite(state.zHighlightThr))
              ? state.zHighlightThr : 1.8;

  for (let i = 0; i < d.n_windows; i++) {
    const z = zs[i];                              // signed or abs depending on mode
    const zAbs = Math.abs(z);                     // for thresholding (significance is symmetric)
    let color;
    if (zMode === 'gradient') {
      // Continuous gradient: 0..zMax mapped through HSL hue 220→0 (blue→red)
      // with brightness scaling so very low-|Z| is muted.
      const t = Math.max(0, Math.min(1, zAbs / zMax));
      const hue = (1 - t) * 220;          // 220° (blue) at low, 0° (red) at high
      const sat = 70;                      // %
      const lit = 35 + 25 * t;             // 35→60% lightness
      const alpha = 0.30 + 0.70 * t;       // dim at bottom, bright at top
      color = `hsla(${hue.toFixed(0)},${sat}%,${lit}%,${alpha.toFixed(2)})`;
    } else if (zMode === 'zone') {
      // In current focal L2 → bright green (matches L2 zone bar palette).
      // In current L1 (but not L2) → muted blue. Outside both → very dim grey.
      // v3.99 turn 14e+: colors aligned to the darker turn-13 palette so the
      // dots match the L1/L2 zone-bar colors above. Old rgba(0,230,118)
      // and rgba(79,163,255) were the pre-turn-13 washed-out values.
      const inL2 = (curL2 >= 0) && (state.windowToL2 && state.windowToL2[i] === curL2);
      const inL1 = (curL1 >= 0) && (state.windowToL1 && state.windowToL1[i] === curL1);
      if (inL2)      color = 'rgba(34,160,80,0.85)';
      else if (inL1) color = 'rgba(48,116,200,0.55)';
      else           color = 'rgba(180,180,180,0.18)';
    } else if (zMode === 'highlight') {
      // Bright above threshold, dim below. Threshold itself shown as a
      // dashed horizontal line for clarity (drawn after the dot loop).
      if (zAbs >= zThr) {
        color = 'rgba(245,165,36,0.92)';   // bright amber
      } else {
        color = 'rgba(170,170,170,0.22)';  // muted grey
      }
    } else {
      // Default: 'bands' — palette ported from the R z-score profile plot
      // (v3.60). Five bins that spread the high-|Z| range across three colors
      // so peaks visibly differentiate. Bins use 1.5/2/3/4 breakpoints, which
      // are visualization-only and intentionally distinct from the biological
      // S1L/S1M/S1S thresholds (1.2/1.8/2.5) drawn as dashed reference lines
      // elsewhere in the panel. Compared on |Z| so signed mode still
      // highlights extreme negative deviations.
      if      (zAbs >= 4)   color = 'rgba(26, 47, 112, 0.95)';   // deep navy
      else if (zAbs >= 3)   color = 'rgba(45, 87, 160, 0.90)';   // rich blue
      else if (zAbs >= 2)   color = 'rgba(79, 127, 184, 0.80)';  // mid blue
      else if (zAbs >= 1.5) color = 'rgba(166, 197, 229, 0.65)'; // pale ice blue
      else                  color = 'rgba(200, 200, 200, 0.32)'; // muted grey (z<1.5)
    }
    ctx.fillStyle = color;
    const x = toX(mbs[i]), y = toY(z);
    ctx.beginPath(); ctx.arc(x, y, 1.4, 0, Math.PI * 2); ctx.fill();
  }
  // Draw the highlight threshold as a dashed amber line in 'highlight' mode
  if (zMode === 'highlight' && zThr > 0 && zThr <= zMax) {
    ctx.setLineDash([4, 4]);
    ctx.strokeStyle = 'rgba(245,165,36,0.65)';
    ctx.lineWidth = 1;
    const yT = toY(zThr);
    ctx.beginPath(); ctx.moveTo(pad.l, yT); ctx.lineTo(pad.l + plotW, yT); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(245,165,36,0.92)';
    ctx.font = '9px ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.fillText(`|Z|≥${zThr.toFixed(1)}`, pad.l + 4, yT - 3);
  }

  // L2 boundary peaks colored by validation_status
  // v3.61: rendered as small UPWARD triangles in the dedicated peaksH band
  // between the L2 zone bar and the z-scatter (was filled circles at the
  // bottom of the plot area, on the dashed yellow line). Triangles point up
  // so they read as "boundary marker tied to the L1/L2 envelope above".
  // v4 turn 10: W-row (per-window candidate-draft strip), expanded mode.
  // Sits directly below the L2 zone bar; only visible when in candidate mode.
  // _wRowEx was computed up front so plotH already accounts for this band.
  // v4 turn 33: nav-lane (always-on) is drawn FIRST so it sits between L2
  // and the W-row. _navBandEx was reserved up front so plotH/toY/arrows
  // already account for it via _vertExtra.
  if (_navBandEx) {
    _drawWinNavLane(ctx, d, toX, xOfWin, _navBandEx, pad.l);
  }
  if (_wRowEx) {
    _drawWRow(ctx, d, toX, xOfWin, _wRowEx);
  }
  if (Array.isArray(d.l2_boundaries)) {
    // v4 turn 10: + _wExtra so arrows sit below the W-row, not on top of it.
    // v4 turn 33: switched to _vertExtra so arrows sit below the always-on nav-lane too.
    const yArrowTop  = zoneTop + zoneH + _vertExtra + 2;          // pointed tip touches just below zone bars (or W-row)
    const yArrowBase = zoneTop + zoneH + _vertExtra + peaksH - 1; // triangle base
    const halfW      = 3;
    for (const b of d.l2_boundaries) {
      const wi = (b.boundary_w | 0) - 1;
      if (wi < 0 || wi >= d.n_windows) continue;
      const x = toX(d.windows[wi].center_mb);
      const col = STATUS_COLOR[b.validation_status] || '#666';
      ctx.fillStyle = col;
      ctx.beginPath();
      ctx.moveTo(x, yArrowTop);                   // tip
      ctx.lineTo(x - halfW, yArrowBase);          // base-left
      ctx.lineTo(x + halfW, yArrowBase);          // base-right
      ctx.closePath();
      ctx.fill();
    }
  }
  // v3.99 t14e+ continue: Candidate bar moved ABOVE the L1/L2 zone bars per
  // Quentin's request — pending candidates render grey, confirmed render
  // shiny gold + ★. Drawn at pad.t (top of canvas plot area, before the
  // L1/L2 zone bars at zoneTop). The old strip-below-boundary-arrows
  // location has been removed.
  // v4 turn 13: candBarTotal = candBarH × n_lanes (lane-stacking when
  // overlapping candidates exist; identical to candBarH when n_lanes=1).
  drawCandidateBar(ctx, d, toX, pad.t, candBarTotal);
  // v4 turn 1 ask 5: "cand" label to the left of the candidate strip,
  // mirroring the L1/L2 label pattern. Same amber color, same right-
  // aligned at pad.l - 4 anchor. See the compact-mode call site above
  // for full rationale.
  ctx.fillStyle = 'rgba(245, 165, 36, 0.95)';
  ctx.font = '8px ui-monospace, monospace';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  ctx.fillText('cand', pad.l - 4, pad.t + candBarTotal / 2);

  // Axes labels
  // v3.50: previous code rendered an integer label at EVERY z from 0..zMax,
  // which became unreadable for zMax > ~10 (e.g. zMax=51 → 52 stacked
  // labels). Use niceTicks for ~5 well-spaced labels at round numbers.
  ctx.fillStyle = themeColor('ink-dim');
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'right';
  const zTicks = niceTicks(zMin, zMax, 5);
  for (const z of zTicks) {
    const y = toY(z);
    // Stay below the zone-bar + peaks strip (v3.61)
    // v3.99 t14e+ continue: use zoneTop instead of pad.t
    // v4 turn 10: + _wExtra so guard tracks toY in W-row mode.
    if (y < zoneTop + zoneH + _wExtra + peaksH - 2) continue;
    // Use 1 decimal for fractional ticks (e.g. 2.5), no decimals for integers
    const lbl = (Math.abs(z - Math.round(z)) < 1e-9) ? z.toFixed(0) : z.toFixed(1);
    ctx.fillText(lbl, pad.l - 5, y + 3);
  }
  ctx.textAlign = 'center';
  const mbTicks = niceTicks(mbMin, mbMax, 6);
  for (const mb of mbTicks) {
    const x = toX(mb);
    // v4 turn 10: + _wExtra so mb-axis labels stay aligned with plot bottom.
    ctx.fillText(mb.toFixed(0), x, zoneTop + zoneH + _wExtra + peaksH + plotH + 14);
  }
  // v3.99 turn 13 ask 5: chromosome name + "(Mb)" label moved from the
  // right edge to the left edge of the X-axis.
  // v3.99 turn 14e: anchored further left at x=2 (was pad.l=44) so the
  // label sits in the y-axis margin area starting near the panel's left
  // edge instead of starting at where the plot area begins. Frees the
  // pad.l..pad.l+plotW zone for the leftmost mb tick label without
  // collision.
  ctx.textAlign = 'left';
  ctx.fillText(d.chrom + ' (Mb)', 2, h - 4);

  // v4 turn 114a: cross-species breakpoint overlay lines.
  // Renders a thin red dashed vertical line at each cs-breakpoint Gar
  // genomic position, spanning the full plot region. Distinct from the
  // existing cursor (orange solid) and L1/L2 boundaries (blue/green solid).
  // Source data: state.crossSpecies.breakpoints filtered by chrom (built
  // by _ensureCsOverlayIndex).
  try {
    const csIdx = (typeof _ensureCsOverlayIndex === 'function')
      ? _ensureCsOverlayIndex() : null;
    if (csIdx && csIdx.bps.length > 0) {
      ctx.save();
      ctx.strokeStyle = '#e85a5a';
      ctx.lineWidth = 1.2;
      ctx.setLineDash([4, 3]);
      ctx.globalAlpha = 0.85;
      const yTop = pad.t;
      const yBot = zoneTop + zoneH + _wExtra + peaksH + plotH;
      for (const e of csIdx.bps) {
        if (e.mb < mbMin || e.mb > mbMax) continue;
        const xx = toX(e.mb);
        if (!Number.isFinite(xx)) continue;
        ctx.beginPath();
        ctx.moveTo(xx + 0.5, yTop);
        ctx.lineTo(xx + 0.5, yBot);
        ctx.stroke();
      }
      ctx.restore();
    }
  } catch (_) { /* fail-soft: overlay is decorative */ }

  // Current cursor
  // v3.99 t14e+ continue: cursor extends from pad.t (top — through candidate
  // bar) all the way down to plot bottom, so it visually links the candidate
  // strip with the L1/L2 bars and the z-scatter as one continuous indicator.
  // v4 turn 10: + _wExtra so the cursor reaches plot bottom in W-row mode too.
  const xCur = toX(d.windows[state.cur].center_mb);
  ctx.strokeStyle = '#f5a524';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(xCur + 0.5, pad.t);
  ctx.lineTo(xCur + 0.5, zoneTop + zoneH + _wExtra + peaksH + plotH);
  ctx.stroke();
}

// --- drawLinesPanel() — legacy lines 34894-35744 ---
export function drawLinesPanel(state) {
  if (!state.data) return;
  const container = document.getElementById('linesCanvasContainer');
  if (!container || typeof container.querySelectorAll !== 'function') return;
  const subs = container.querySelectorAll('.lines-subpanel');
  if (!subs || subs.length === 0) return;

  const d = state.data;
  const nWin = d.n_windows;
  const nS = d.n_samples;
  if (nWin < 2 || nS === 0) return;

  const trackedSet = new Set(state.tracked);
  const mbs = d.windows.map(w0 => w0.center_mb);
  const _mbR = currentMbRange();
  const mbMin = _mbR.mbMin, mbMax = _mbR.mbMax;

  for (const sub of subs) {
    const source = sub.dataset.linesSource;
    const cv = sub.querySelector('canvas');
    if (!cv) continue;
    const { ctx, w, h } = fitCanvas(cv);
    ctx.clearRect(0, 0, w, h);

    // v3.99 turn 6: pad.t reduced 18 → 6 for tighter compact layout, per
     // Quentin's request. MUST match _computeLinesLassoSamples pad above.
    const pad = { l: 44, r: 16, t: 6, b: 8 };
    const plotW = w - pad.l - pad.r;
    const plotH = h - pad.t - pad.b;
    if (plotW <= 0 || plotH <= 0) continue;

    // Resolve the source's native grid. PC sources use the precomp grid;
    // GHSL sources (het, ghsl_div_<scale>) use the GHSL panel grid which has
    // its own n_windows and bp positions. The X-axis is always the precomp's
    // [mbMin, mbMax] so all sub-panels stay vertically aligned along the
    // chromosome regardless of which grid each source uses.
    const sourceGrid = getLinesGrid(source);   // null = precomp grid
    let nGrid;
    let mbAt;   // function: gridIdx -> mb position
    if (sourceGrid) {
      nGrid = sourceGrid.length;
      mbAt = (gi) => sourceGrid[gi] / 1e6;
    } else {
      nGrid = nWin;
      mbAt = (gi) => mbs[gi];
    }
    if (nGrid < 2) {
      ctx.fillStyle = themeColor('dim');
      ctx.font = '12px ui-monospace, monospace';
      ctx.textAlign = 'center';
      ctx.fillText(`(${source.toUpperCase()} not available in this dataset)`, w / 2, h / 2);
      continue;
    }

    // v3.99 turn 7 perf: per-source cache for the expensive y-matrix +
    // y-bounds + x-by-grid. None of these depend on state.cur — they're
    // pure functions of (source data, source sign, plotW, plotH, mb range,
    // pad). On window-by-window stepping the inputs don't change, so we
    // can reuse the entire computation. Cache key includes everything that
    // affects the values; on cache hit we just re-stroke the polylines on
    // top of the freshly cleared canvas.
    if (!state.__linesCache) state.__linesCache = {};
    const cacheKey = [
      'v1',
      source,
      nGrid,
      mbMin.toFixed(6), mbMax.toFixed(6),
      plotW, plotH,
      pad.l, pad.r, pad.t, pad.b,
      // sign array is per-grid-step; the data it derives from rarely
      // changes during stepping, but include a coarse marker so flips
      // (sign-align toggle) invalidate. state.flipPC1 controls whether
      // the per-window sign-flip is applied to the PC1 source.
      state.flipPC1 ? 'sa' : '',
      // data identity — if the user loads a different chromosome the cache
      // must miss
      state.data && state.data.chrom ? state.data.chrom : '',
    ].join('|');
    let cached = state.__linesCache[source];
    if (!cached || cached._key !== cacheKey) {
      // Compute Y range across all samples and all windows for this source
      let yMin = Infinity, yMax = -Infinity;
      let validData = false;
      for (let gi = 0; gi < nGrid; gi++) {
        const vals = getLinesValuesAt(gi, source);
        if (!vals) continue;
        const sign = getLinesSignAt(gi, source);
        for (let si = 0; si < nS; si++) {
          const v = vals[si] * sign;
          if (!isFinite(v)) continue;
          if (v < yMin) yMin = v;
          if (v > yMax) yMax = v;
          validData = true;
        }
      }
      if (!validData) {
        // Source not present — render an info message (not cached)
        ctx.fillStyle = themeColor('dim');
        ctx.font = '12px ui-monospace, monospace';
        ctx.textAlign = 'center';
        ctx.fillText(`(${source.toUpperCase()} not available in this dataset)`, w / 2, h / 2);
        continue;
      }
      const yPad = (yMax - yMin) * 0.05 || 0.01;
      yMin -= yPad; yMax += yPad;
      const toY = (v) => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;
      const xByGrid = new Float32Array(nGrid);
      for (let gi = 0; gi < nGrid; gi++) {
        xByGrid[gi] = pad.l + ((mbAt(gi) - mbMin) / (mbMax - mbMin)) * plotW;
      }
      const yMatrix = new Array(nS);
      for (let si = 0; si < nS; si++) {
        const arr = new Float32Array(nGrid);
        for (let i = 0; i < nGrid; i++) arr[i] = NaN;
        yMatrix[si] = arr;
      }
      for (let gi = 0; gi < nGrid; gi++) {
        const vals = getLinesValuesAt(gi, source);
        if (!vals) continue;
        const sign = getLinesSignAt(gi, source);
        for (let si = 0; si < nS; si++) {
          const v = vals[si] * sign;
          if (isFinite(v)) yMatrix[si][gi] = toY(v);
        }
      }
      cached = { _key: cacheKey, yMin, yMax, xByGrid, yMatrix };
      state.__linesCache[source] = cached;
    }
    const yMin = cached.yMin;
    const yMax = cached.yMax;
    const xByGrid = cached.xByGrid;
    const yMatrix = cached.yMatrix;
    // toY still defined for downstream cursor-overlay code that needs to
    // forward-map a y value not in the cache (rare — kept for compatibility
    // with the geometry stash + click-inspect popover).
    const toY = (v) => pad.t + plotH - ((v - yMin) / (yMax - yMin)) * plotH;
    const toX = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;

    // v3.87: stash this sub-canvas's geometry on state so the click-inspect
    // popover can forward-map per-sample values into canvas y to find the
    // closest tracked fish at a click position. Keyed by source name.
    if (!state.__linesGeom) state.__linesGeom = {};
    state.__linesGeom[source] = {
      pad: { l: pad.l, r: pad.r, t: pad.t, b: pad.b },
      plotW, plotH, w, h,
      mbMin, mbMax, yMin, yMax,
      nGrid, sourceGrid: !!sourceGrid,
    };

    // turn 141 Slice 1: per-candidate vertical band highlights.
    // SPEC §2.4 — bands paint FIRST so they're pure background; lines,
    // tracking highlights, the cursor crosshair, and the inheritance
    // strip all stack on top. Gated on state.linesPanelCandidateBands
    // (default ON, persisted). Helper handles all filtering / clipping
    // / palette assignment. Wrapped in try/catch so a misbehaving
    // candidate (e.g. a malformed bp value) can't take out the whole
    // panel — band drawing is purely additive and safe to skip.
    if (state.linesPanelCandidateBands !== false &&
        typeof _paintCandidateBands === 'function') {
      try {
        _paintCandidateBands(ctx, {
          pad, plotW, plotH, toX, mbMin, mbMax,
          candidates: Array.isArray(state.candidateList) ? state.candidateList : [],
          chrom: (state.data && state.data.chrom) ? state.data.chrom : null,
        });
      } catch (e) {
        // Defensive: never let band-drawing failure block the rest of
        // the lines panel from rendering.
        if (typeof console !== 'undefined' && console.warn) {
          console.warn('[lines/candBands] paint failed:', e && e.message);
        }
      }
    }

    // Frame
    ctx.strokeStyle = themeColor('rule');
    ctx.lineWidth = 1;
    ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, plotW, plotH);

    // Y-axis ticks (3 ticks: min, mid, max)
    // v3.79: switched from themeColor('dim') (--ink-dim, light-grey) to
    // themeColor('ink') (--ink, white in dark mode / dark-grey in light) for
    // full contrast — the dim tone was reading as black/illegible against
    // dark backgrounds in some configurations.
    ctx.fillStyle = themeColor('ink');
    ctx.font = '9px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const yMid = (yMin + yMax) / 2;
    for (const v of [yMin + (yMax - yMin) * 0.05, yMid, yMax - (yMax - yMin) * 0.05]) {
      const yp = toY(v);
      ctx.fillText(formatTrackVal(v), pad.l - 4, yp);
      ctx.strokeStyle = 'rgba(42,50,66,0.4)';
      ctx.beginPath();
      ctx.moveTo(pad.l, yp + 0.5); ctx.lineTo(pad.l + plotW, yp + 0.5);
      ctx.stroke();
    }

    // v3.99 turn 7 perf: xByGrid + yMatrix were previously built here on
    // every draw. They're now provided by state.__linesCache[source] above
    // (rebuilt only when inputs change). Helpers below (strokeSamplePath,
    // strokeSamplePathStyled) close over the cached references.

    // Helper to stroke a polyline with NaN gap handling
    function strokeSamplePath(si) {
      const ys = yMatrix[si];
      let started = false;
      ctx.beginPath();
      for (let gi = 0; gi < nGrid; gi++) {
        const y = ys[gi];
        if (!isFinite(y)) { started = false; continue; }
        const x = xByGrid[gi];
        if (!started) { ctx.moveTo(x, y); started = true; }
        else { ctx.lineTo(x, y); }
      }
      ctx.stroke();
    }

    // v3.85: per-segment-styled stroke for tracked fish whose K=3 vote inside
    // the active candidate differs from their consensus regime. Inputs:
    //   si:        sample index
    //   jumpMask:  Uint8Array(nGrid), 1 = this window falls inside an L2
    //              where this fish jumped parents
    //   normalStyle: { strokeStyle, lineWidth, dash }
    //   jumpStyle:   { strokeStyle, lineWidth, dash }
    // The polyline is broken into runs of constant mask value, each stroked
    // with its own style. NaN gaps still break the line. A small overlap of
    // 1 vertex between adjacent runs keeps the visual transition smooth.
    function strokeSamplePathStyled(si, jumpMask, normalStyle, jumpStyle) {
      const ys = yMatrix[si];
      // Walk gi, accumulate runs
      let runStart = -1;
      let runMask = 0;
      const flushRun = (giEnd) => {
        if (runStart < 0) return;
        const style = runMask ? jumpStyle : normalStyle;
        ctx.save();
        ctx.strokeStyle = style.strokeStyle;
        ctx.lineWidth = style.lineWidth;
        if (style.dash) ctx.setLineDash(style.dash); else ctx.setLineDash([]);
        ctx.beginPath();
        let started = false;
        // Include 1 extra vertex on either end (where finite) so segments
        // visually meet. Use clamped indices.
        const giA = Math.max(0, runStart - (runMask ? 0 : 0));
        const giB = Math.min(nGrid - 1, giEnd + (runMask ? 0 : 0));
        for (let gi = giA; gi <= giB; gi++) {
          const y = ys[gi];
          if (!isFinite(y)) { started = false; continue; }
          const x = xByGrid[gi];
          if (!started) { ctx.moveTo(x, y); started = true; }
          else { ctx.lineTo(x, y); }
        }
        ctx.stroke();
        ctx.restore();
        runStart = -1;
      };
      for (let gi = 0; gi < nGrid; gi++) {
        const m = jumpMask[gi];
        if (runStart < 0) { runStart = gi; runMask = m; continue; }
        if (m !== runMask) {
          flushRun(gi);   // close previous run at gi-1, but include gi as overlap
          runStart = gi - 1;   // start next run one vertex earlier so segments meet
          if (runStart < 0) runStart = 0;
          runMask = m;
        }
      }
      flushRun(nGrid - 1);
    }

    // Draw untracked first as a single batched stroke per sample (cheap)
    // v3.99 turn 8 perf: the untracked-gray cloud (~200-220 polylines × nGrid
    // segments) is a pure function of (yMatrix, xByGrid, trackedSet, w, h).
    // None of those depend on state.cur. Cache an offscreen canvas with the
    // pre-stroked cloud and blit it on each frame instead of re-stroking.
    // Cache key includes the trackedSet hash since picking a sample to track
    // moves it from the untracked-gray bucket to the colored-overlay bucket.
    // v4 turn 108: cache key also includes state.linesColorMode so switching
    // from kmeans to family (or any other mode) invalidates the cache and
    // triggers a re-stroke with per-sample colors.
    let trackedHashLocal = 0;
    for (const si of trackedSet) trackedHashLocal += si | 0;
    trackedHashLocal = (trackedHashLocal * 31) ^ (trackedSet.size << 16);
    const lcMode = state.linesColorMode || 'kmeans';
    const overlayCacheKey = cacheKey + '|tr=' + trackedHashLocal +
                            '|wh=' + w + 'x' + h + '|lc=' + lcMode;
    if (!cached.bgCanvas || cached.bgCanvasKey !== overlayCacheKey) {
      // (Re)build offscreen canvas. Use a plain canvas (matches main 2D
      // context) sized to the visible CSS pixels of the live canvas — so
      // drawImage maps 1:1 without scaling artifacts.
      const off = (typeof OffscreenCanvas !== 'undefined')
        ? new OffscreenCanvas(Math.max(1, w | 0), Math.max(1, h | 0))
        : (function() {
            const c = document.createElement('canvas');
            c.width = Math.max(1, w | 0);
            c.height = Math.max(1, h | 0);
            return c;
          })();
      const offCtx = off.getContext('2d');
      offCtx.lineWidth = 0.6;
      // v4 turn 108: per-sample coloring. For modes that produce one color
      // per sample (family, F_ROH, kmeans-as-stable-band), call the resolver
      // per sample. For 'kmeans' we keep the legacy grey-cloud behavior
      // because the per-window lane assignments are conveyed by the WALKING
      // of each line through the band y-positions, not by line color (every
      // sample is the same grey). For 'family', each fish gets its family
      // color across all windows; alpha bumped to 0.25 so saturated colors
      // remain readable when 226 lines overlap.
      const usePerSampleColor = (lcMode === 'family');
      const baseAlpha = usePerSampleColor ? 0.25 : 0.10;
      const defaultStroke = `rgba(180,190,210,${baseAlpha.toFixed(3)})`;
      offCtx.strokeStyle = defaultStroke;
      // v4 turn 126: track how many samples got a real per-sample color.
      // If we're in family mode but every sample falls back to the default
      // stroke (because family_id isn't loaded on samples), the visual
      // result is identical to kmeans mode — Quentin reported "color by
      // family doesn't color." A notice gets drawn on the live canvas
      // below to make this state visible instead of silently no-op'ing.
      let _perSampleColorHits = 0;
      for (let si = 0; si < nS; si++) {
        if (trackedSet.has(si)) continue;
        // Pick per-sample stroke color when in a per-sample-coloring mode.
        if (usePerSampleColor && typeof _resolveSampleScopeColor === 'function') {
          const c = _resolveSampleScopeColor(si, lcMode);
          if (c) {
            _perSampleColorHits++;
            // c may be 'rgb(r,g,b)' or '#rrggbb' — wrap with alpha via
            // withAlpha if available, else pass through.
            offCtx.strokeStyle = (typeof withAlpha === 'function')
              ? withAlpha(c, baseAlpha)
              : c;
          } else {
            offCtx.strokeStyle = defaultStroke;
          }
        }
        const ys = yMatrix[si];
        let started = false;
        offCtx.beginPath();
        for (let gi = 0; gi < nGrid; gi++) {
          const y = ys[gi];
          if (!isFinite(y)) { started = false; continue; }
          const x = xByGrid[gi];
          if (!started) { offCtx.moveTo(x, y); started = true; }
          else { offCtx.lineTo(x, y); }
        }
        offCtx.stroke();
      }
      cached.bgCanvas = off;
      cached.bgCanvasKey = overlayCacheKey;
      // Stash the diagnostic so the live-canvas pass below can draw a notice.
      // We stash it on `cached` so subsequent frames (when the cache is
      // re-blitted instead of re-built) can still surface the notice.
      cached.bgFamilyMissing = (usePerSampleColor && _perSampleColorHits === 0);
      cached.bgFamilyDistinct = 0;
      if (usePerSampleColor) {
        // Count distinct non-fallback colors so the notice can hint at it.
        const distinct = new Set();
        for (let si = 0; si < nS; si++) {
          if (typeof _resolveSampleScopeColor !== 'function') break;
          const c = _resolveSampleScopeColor(si, lcMode);
          if (c) distinct.add(c);
        }
        cached.bgFamilyDistinct = distinct.size;
      }
    }
    // Blit the pre-rendered untracked cloud (single drawImage call instead
    // of ~220 strokes × nGrid segments). Massive speedup on stepping.
    ctx.drawImage(cached.bgCanvas, 0, 0);

    // v4 turn 126: family-mode "no data" notice. When the user picks
    // "color: family" but the loaded JSON doesn't have family_id on samples,
    // every line falls back to the default grey stroke — looking identical
    // to kmeans mode and producing Quentin's "color by family doesn't color"
    // report. Drawn ONLY on PC1 sub-panel (one notice, not per-subpanel) and
    // ONLY when the diagnostic flagged us as no-hits. Notice is small,
    // amber-tinted, top-right of the plot so it doesn't obscure data.
    if (source === 'pc1' && cached.bgFamilyMissing) {
      ctx.save();
      const msg = 'family mode: no family data loaded — drag-drop ngsRelate JSON';
      ctx.font = '10px ui-monospace, monospace';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'top';
      const textW = ctx.measureText(msg).width;
      const noticeX = w - 8;
      const noticeY = pad.t + 4;
      // Faint amber backdrop so the text reads against the line cloud
      ctx.fillStyle = 'rgba(245, 165, 36, 0.15)';
      ctx.fillRect(noticeX - textW - 6, noticeY - 2, textW + 8, 14);
      ctx.strokeStyle = 'rgba(245, 165, 36, 0.55)';
      ctx.lineWidth = 0.8;
      ctx.strokeRect(noticeX - textW - 6 + 0.5, noticeY - 2 + 0.5, textW + 8 - 1, 14 - 1);
      ctx.fillStyle = 'rgba(245, 165, 36, 0.95)';
      ctx.fillText(msg, noticeX - 2, noticeY);
      ctx.restore();
    }

    // v3.76: Lasso mode for the PC1 sub-canvas — suppress tracked-color so
    // the user can pick fresh by trajectory. Tracked samples still render but
    // in dim grey (same treatment as untracked, slightly brighter so they're
    // identifiable). Lassoed samples render in accent color so the user sees
    // the live selection. Other sub-canvases ignore lasso mode entirely.
    const isPC1 = (source === 'pc1');
    const lassoOnHere = isPC1 && state.linesLassoActive;
    const lassoSet = lassoOnHere ? new Set(state.linesLassoSelected || []) : null;
    if (lassoOnHere) {
      // Tracked samples → dim grey (slightly brighter than untracked so still readable)
      ctx.lineWidth = 0.8;
      ctx.strokeStyle = 'rgba(200,210,225,0.18)';
      for (const si of state.tracked) {
        if (lassoSet && lassoSet.has(si)) continue;   // will render in accent below
        strokeSamplePath(si);
      }
      // Lassoed samples → accent color (gold)
      if (lassoSet && lassoSet.size > 0) {
        ctx.lineWidth = 1.4;
        ctx.strokeStyle = '#f5a524';
        for (const si of lassoSet) strokeSamplePath(si);
      }
    } else {
      // Tracked lines (full color, full opacity) — original behavior
      // v3.85: when an active candidate has fish_calls and this is the PC1
      // sub-canvas, build a per-fish jumpMask: 1 = this grid index falls
      // inside an L2 where the fish's K=3 vote differs from its consensus
      // regime. Tracked fish with no jumps stroke the normal way; jumpers
      // get strokeSamplePathStyled which breaks the polyline at L2 boundaries
      // and renders the "jumped" segments in dashed red.
      const cand = (isPC1 && state.candidate) ? state.candidate : null;
      const candFishCalls = (cand && Array.isArray(cand.fish_calls))
        ? cand.fish_calls : null;
      const candL2s = (cand && Array.isArray(cand.l2_indices))
        ? cand.l2_indices : null;
      // Precompute per-L2 mb-range for the candidate (so we can map gi → L2idx
      // by mb position without rebuilding the lookup per fish).
      // v3.86: when the lines panel is rendering on the precomp grid (i.e.
      // sourceGrid is null and nGrid === nWin), gi IS the window index, so
      // we compare integer window ranges directly via env._s0/_e0 — both
      // faster (no float ops) and more accurate at L2 boundaries (no rounding
      // edge cases when boundary mb falls exactly on a window center).
      // For non-precomp sources (GHSL etc.) the grid is mb-positioned and we
      // fall back to the mb-range path.
      let l2BpRanges = null;     // mb-based, used when sourceGrid != null
      let l2WinRanges = null;    // window-index-based, used when sourceGrid == null
      let gridSlotByGi = null;   // Int16Array(nGrid): slot index (0..candL2s.length-1)
                                 // covering this gi, or -1 if outside all L2s.
                                 // Pre-bucketed once so the per-fish loop is O(nGrid),
                                 // not O(nGrid * nL2). Int16 (not Int8) because some
                                 // candidates merge >127 L2s on long chromosomes.
      if (candL2s && candL2s.length > 0) {
        const envs = d.l2_envelopes || [];
        if (sourceGrid) {
          l2BpRanges = candL2s.map(li => {
            const env = envs[li];
            return env ? [env.start_bp / 1e6, env.end_bp / 1e6] : null;
          });
          // Pre-bucket gi → slot by mb position
          gridSlotByGi = new Int16Array(nGrid).fill(-1);
          for (let gi = 0; gi < nGrid; gi++) {
            const mb = mbAt(gi);
            for (let p = 0; p < l2BpRanges.length; p++) {
              const r = l2BpRanges[p];
              if (!r) continue;
              if (mb >= r[0] && mb <= r[1]) { gridSlotByGi[gi] = p; break; }
            }
          }
        } else {
          // Precomp grid: compare integer window indices.
          l2WinRanges = candL2s.map(li => {
            const env = envs[li];
            return env ? [env._s0, env._e0] : null;
          });
          gridSlotByGi = new Int16Array(nGrid).fill(-1);
          for (let gi = 0; gi < nGrid; gi++) {
            for (let p = 0; p < l2WinRanges.length; p++) {
              const r = l2WinRanges[p];
              if (!r) continue;
              if (gi >= r[0] && gi <= r[1]) { gridSlotByGi[gi] = p; break; }
            }
          }
        }
      }
      // Helper: build the mask for a given fish_call entry. Returns null when
      // no jumping (so caller can use the cheap path). Uses the pre-bucketed
      // gridSlotByGi so the lookup is O(1) per gi.
      function _buildJumpMask(fc) {
        if (!fc || !candFishCalls || !candL2s || !gridSlotByGi) return null;
        if (fc.subband_stability == null) return null;
        if (fc.subband_stability >= 1) return null;
        if (!Array.isArray(fc.votes) || fc.regime < 0) return null;
        const mask = new Uint8Array(nGrid);
        let anyJumped = false;
        for (let gi = 0; gi < nGrid; gi++) {
          const slot = gridSlotByGi[gi];
          if (slot < 0) continue;
          const vote = fc.votes[slot];
          if (vote >= 0 && vote !== fc.regime) {
            mask[gi] = 1;
            anyJumped = true;
          }
        }
        return anyJumped ? mask : null;
      }
      const normalStyle_jump = { strokeStyle: '#e0555c', lineWidth: 1.6, dash: [4, 3] };
      ctx.lineWidth = 1.4;
      for (const si of state.tracked) {
        const fc = candFishCalls ? candFishCalls[si] : null;
        const mask = _buildJumpMask(fc);
        if (mask) {
          // Mixed: normal segments in tracked color, jumped segments in dashed red
          strokeSamplePathStyled(si, mask,
            { strokeStyle: trackedColor(si), lineWidth: 1.4, dash: null },
            normalStyle_jump
          );
        } else {
          ctx.strokeStyle = trackedColor(si);
          strokeSamplePath(si);
        }
      }
    }

    // Crosshair at current window. state.cur indexes into the precomp grid;
    // we draw it at the corresponding mb position so all panels' crosshairs
    // line up vertically with the precomp-driven page layout.
    const curMb = (state.cur >= 0 && state.cur < nWin) ? mbs[state.cur] : NaN;
    if (isFinite(curMb)) {
      const xc = toX(curMb);
      ctx.strokeStyle = '#f5a524';
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(Math.round(xc) + 0.5, pad.t);
      ctx.lineTo(Math.round(xc) + 0.5, pad.t + plotH);
      ctx.stroke();
    }

    // v3.84: active-candidate overlay on the PC1 sub-panel. When a saved
    // candidate is active (state.candidate set) AND it has fish_calls (i.e.
    // committed via the v3.80+ flow), draw:
    //   1. A translucent amber band shading the candidate's bp-span horizontally
    //   2. Vertical dashed lines at each supporting L2's right boundary
    //   3. A small ⚠ glyph at the right edge of the plot for each tracked fish
    //      whose subband_stability < 1 (i.e. crossed K=3 parents within this
    //      candidate). The glyph sits at the fish's y-position at the last
    //      window of the candidate span — tells the user at a glance which
    //      tracked samples are regime-jumpers within this candidate.
    if (isPC1 && state.candidate && state.candidate.l2_indices &&
        state.candidate.l2_indices.length > 0) {
      const cand = state.candidate;
      const envs = d.l2_envelopes || [];
      // v3.91: reset marker hit list at the top of the overlay block. The
      // populated path inside (if (lastGi >= 0)) overwrites this; clearing
      // here keeps the list in sync when the candidate scrolls out of view
      // or all tracked fish have stability >= 1.
      cv.__candJumperMarkers = [];
      state.__candJumperMarkers = [];
      // Map candidate bp-span to screen x. Use start_bp/end_bp directly so
      // the band aligns with the chromosome coords on the existing axis.
      const candStartMb = cand.start_bp / 1e6;
      const candEndMb = cand.end_bp / 1e6;
      // Clip to visible mb range
      const visStartMb = Math.max(candStartMb, mbMin);
      const visEndMb = Math.min(candEndMb, mbMax);
      if (visEndMb > visStartMb) {
        const xLo = toX(visStartMb);
        const xHi = toX(visEndMb);
        // Translucent band
        ctx.save();
        ctx.fillStyle = 'rgba(245,196,58,0.06)';
        ctx.fillRect(xLo, pad.t, xHi - xLo, plotH);
        // Top + bottom edges
        ctx.strokeStyle = 'rgba(245,196,58,0.40)';
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        ctx.moveTo(xLo, pad.t + 0.5); ctx.lineTo(xHi, pad.t + 0.5);
        ctx.moveTo(xLo, pad.t + plotH - 0.5); ctx.lineTo(xHi, pad.t + plotH - 0.5);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.restore();
        // L2 boundary verticals (between supporting intervals; skip first
        // and last since those are the band edges)
        ctx.save();
        ctx.strokeStyle = 'rgba(245,196,58,0.25)';
        ctx.lineWidth = 1;
        ctx.setLineDash([2, 4]);
        for (let i = 0; i < cand.l2_indices.length - 1; i++) {
          const env = envs[cand.l2_indices[i]];
          if (!env) continue;
          const boundaryMb = env.end_bp / 1e6;
          if (boundaryMb < mbMin || boundaryMb > mbMax) continue;
          const xb = Math.round(toX(boundaryMb)) + 0.5;
          ctx.beginPath();
          ctx.moveTo(xb, pad.t);
          ctx.lineTo(xb, pad.t + plotH);
          ctx.stroke();
        }
        ctx.restore();
        // Jumper markers — tracked fish with subband_stability < 1
        if (Array.isArray(cand.fish_calls) && state.tracked.length > 0) {
          // Get y-source for last visible window in the candidate span. The
          // mbAt(gi) function gives mb position; find the gridIdx whose
          // mb is closest to (and ≤) visEndMb.
          let lastGi = -1;
          for (let gi = 0; gi < nGrid; gi++) {
            if (mbAt(gi) <= visEndMb) lastGi = gi;
            else break;
          }
          if (lastGi >= 0) {
            const valsLast = (typeof getLinesValuesAt === 'function')
              ? getLinesValuesAt(lastGi, source) : null;
            const xMarker = toX(visEndMb) + 6;   // just outside the band
            // v3.91: stash marker positions for click hit-testing. One entry
            // per drawn ⚠ glyph; carries the fish index, marker center on
            // canvas, and the first-jump L2 within this candidate (used to
            // jump the scrubber when the marker is clicked).
            const markerHits = [];
            ctx.save();
            for (const si of state.tracked) {
              const fc = cand.fish_calls[si];
              if (!fc || fc.subband_stability == null) continue;
              if (fc.subband_stability >= 1) continue;
              if (!valsLast) continue;
              const y0 = valsLast[si];
              if (!isFinite(y0)) continue;
              const yp = toY(y0);
              // Small triangle glyph (warning ⚠ in geometry form so it
              // renders consistently across fonts).
              ctx.fillStyle = '#e0555c';
              ctx.strokeStyle = themeColor('bg');
              ctx.lineWidth = 1.5;
              ctx.beginPath();
              ctx.moveTo(xMarker, yp - 4);
              ctx.lineTo(xMarker - 4, yp + 3);
              ctx.lineTo(xMarker + 4, yp + 3);
              ctx.closePath();
              ctx.stroke();   // dark outline first for contrast
              ctx.fill();
              // Exclamation tick inside
              ctx.fillStyle = themeColor('bg');
              ctx.fillRect(xMarker - 0.5, yp - 1, 1, 3);
              // v3.91: find the first L2 in cand.l2_indices where this fish's
              // subband_path's K=3 parent disagrees with fc.regime — that's
              // a "jump". Records the L2 envelope index + center window so a
              // click on this marker can scrub there. Falls back to the last
              // L2 of the candidate if subband_path is missing or no jump
              // detected (shouldn't happen since stability < 1 by definition).
              let jumpL2Idx = -1, jumpWinIdx = -1;
              if (Array.isArray(fc.subband_path) && fc.regime != null) {
                const path = fc.subband_path;
                for (let p = 0; p < path.length; p++) {
                  const tok = path[p];
                  if (typeof tok !== 'string' || tok.length < 2 || tok[0] !== 'g') continue;
                  const parentDigit = tok.charCodeAt(1) - 48;   // '0'..'9' → 0..9
                  if (parentDigit < 0 || parentDigit > 9) continue;
                  if (parentDigit !== fc.regime) {
                    jumpL2Idx = cand.l2_indices[p];
                    break;
                  }
                }
              }
              if (jumpL2Idx < 0 && cand.l2_indices.length > 0) {
                jumpL2Idx = cand.l2_indices[cand.l2_indices.length - 1];
              }
              if (jumpL2Idx >= 0 && envs[jumpL2Idx]) {
                const env2 = envs[jumpL2Idx];
                const s0 = env2._s0, e0 = env2._e0;
                if (Number.isFinite(s0) && Number.isFinite(e0) && e0 >= s0) {
                  jumpWinIdx = (s0 + e0) >> 1;
                }
              }
              markerHits.push({
                si: si,
                x: xMarker,
                y: yp,
                jump_l2_idx: jumpL2Idx,
                jump_win_idx: jumpWinIdx,
              });
            }
            ctx.restore();
            // Stash on canvas + state. Canvas stash is what the click handler
            // consults (it has the cv element directly). State stash makes
            // the data inspectable from console / tests.
            cv.__candJumperMarkers = markerHits;
            state.__candJumperMarkers = markerHits;
          }
        }
      }
    } else if (isPC1) {
      // v3.91: no active candidate (or empty l2_indices) on the PC1 panel —
      // clear any stale marker hit-test entries so unmodified clicks fall
      // through to the existing click-to-jump behavior.
      cv.__candJumperMarkers = [];
      state.__candJumperMarkers = [];
    }


    // turn 123: tracked-linkage shading. When the user has lassoed fish
    // (state.tracked), shade each candidate's bp range on the lines panel
    // by the dominant band's inheritance group, with alpha keyed to purity.
    // Drawn FIRST so subsequent overlays (diamond, transitions, lines) sit on top.
    if (isPC1 && typeof _drawTrackedLinkageStrip === 'function') {
      try {
        _drawTrackedLinkageStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // v4 turn 94: Diamond zone overlay. When the user is zoomed (visible
    // mb range smaller than the candidate's full span), and the active
    // candidate has detectable diamond patterns, draw a translucent cyan
    // background rectangle behind the diamond windows on the PC1 sub-panel,
    // with a small "split detected" annotation strip on top.
    //
    // Only on PC1 sub-panel and only when zoomed enough that the diamond
    // is visually distinguishable (otherwise it looks like a thin smear).
    if (isPC1 && typeof _drawDiamondOverlay === 'function') {
      try {
        _drawDiamondOverlay(ctx, pad, plotW, plotH, mbMin, mbMax, w, h);
      } catch (_) {}
    }

    // v4 turn 95: SNP-density strip on PC1 sub-panel (only when activated
    // via toolbar button + the user is zoomed). Thin gradient bar above the
    // plot showing where SNP density is low (PC1 loses resolution) vs high.
    if (isPC1 && typeof _drawSnpDensityStrip === 'function') {
      try {
        _drawSnpDensityStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // v4 turn 99: SNP-density shade on PC1 sub-panel — translucent vertical
    // bars across plot height for low-density windows. Drawn on top of lines
    // because the alpha is low (≤ 0.18) and the visual "darkening" of low-density
    // columns is the intended effect. Activated via toolbar button (mode='shade').
    if (isPC1 && typeof _drawSnpDensityShade === 'function') {
      try {
        _drawSnpDensityShade(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // v4 turn 102: structural-haplotype transition-rate strip on PC1 sub-panel.
    // Per-L2-boundary transition_rate (fraction of fish that change band)
    // shown as red-graded bars at the bottom of the plot. Hotspots (rate >=
    // 0.30) get a thin vertical tick across the plot height.
    if (isPC1 && typeof _drawTransitionRateStrip === 'function') {
      try {
        _drawTransitionRateStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // v4 turn 104: regime-breadth strip on PC1 sub-panel — categorical
    // per-L2 classification (narrow/medium/wide/no_signal) drawn as a thin
    // colored bar at the top of the plot. Green = narrow (clean segregation),
    // red = wide (no clean segregation), amber = mixed, grey = no signal.
    if (isPC1 && typeof _drawRegimeBreadthStrip === 'function') {
      try {
        _drawRegimeBreadthStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // turn 117: inheritance group labels strip on PC1 sub-panel — small
    // "I1·3g" labels above each candidate region. I1 = sequential candidate
    // number on this chromosome by start_bp; 3g = number of inheritance
    // groups (computed via inheritanceGroupClustering on confirmed
    // candidates' band labels). Placed above the regime-breadth strip.
    if (isPC1 && typeof _drawInheritanceLabelsStrip === 'function') {
      try {
        _drawInheritanceLabelsStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // turn 130 Slice 2: lineage strip — per-L2 dominant lineage
    // (computed from the fish-trajectory clustering, runLineageCompute).
    // Sibling to the regime-breadth strip. Toggle: state.linesLineageStripOn.
    if (isPC1 && typeof _drawLineageStrip === 'function') {
      try {
        _drawLineageStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // turn 161 — band-trace strip (Slice 4 UI). Per-L2 stacked bar of
    // band_fractions for state.bandTraceFishSet, with a regime-color
    // top stripe. Sits above the lineage strip so it doesn't displace
    // the existing diagnostic. Off by default; user toggles via the
    // lines header. Observation-only — no inversion-call markers in
    // this slice (manuscript framing: report co-segregation, do not
    // interpret).
    if (isPC1 && typeof _drawBandTraceStrip === 'function') {
      try {
        _drawBandTraceStrip(ctx, pad, plotW, plotH, mbMin, mbMax);
      } catch (_) {}
    }

    // v3.76: Lasso rectangle overlay on the PC1 panel. Live drag uses
    // linesLassoRect (filled translucent gold + dashed border); after pointer-up
    // the rect persists in linesLassoCommitted (lighter dashed border only) so
    // the user sees what was selected until Confirm or Clear.
    if (isPC1 && state.linesLassoActive) {
      const live = state.linesLassoRect;
      const committed = state.linesLassoCommitted;
      const drawRect = (rect, fillStyle, strokeStyle, lineDash) => {
        if (!rect) return;
        const x0 = Math.min(rect.x0, rect.x1);
        const y0 = Math.min(rect.y0, rect.y1);
        const x1 = Math.max(rect.x0, rect.x1);
        const y1 = Math.max(rect.y0, rect.y1);
        ctx.save();
        ctx.fillStyle = fillStyle;
        ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
        ctx.strokeStyle = strokeStyle;
        ctx.lineWidth = 1.2;
        ctx.setLineDash(lineDash);
        ctx.strokeRect(x0 + 0.5, y0 + 0.5, x1 - x0, y1 - y0);
        ctx.restore();
      };
      if (live) {
        drawRect(live, 'rgba(245,165,36,0.10)', '#f5a524', [4, 3]);
      } else if (committed) {
        drawRect(committed, 'rgba(245,165,36,0.05)', 'rgba(245,165,36,0.6)', [2, 4]);
      }
    }

    // v4 turn 114a: cross-species breakpoint overlay vertical lines.
    // Renders a thin red dashed vertical line at each cs-breakpoint Gar
    // genomic position, spanning the subpanel's plot region. Distinct from
    // L1/L2 boundaries (blue/green) and lasso (gold). Per-subpanel loop so
    // every sub gets the overlay regardless of source (PC1, PC2, GHSL,
    // dosage, etc.). Source data: state.crossSpecies.breakpoints filtered
    // by chrom (built lazily by _ensureCsOverlayIndex).
    try {
      const csIdx = (typeof _ensureCsOverlayIndex === 'function')
        ? _ensureCsOverlayIndex() : null;
      if (csIdx && csIdx.bps.length > 0) {
        ctx.save();
        ctx.strokeStyle = '#e85a5a';
        ctx.lineWidth = 1.2;
        ctx.setLineDash([4, 3]);
        ctx.globalAlpha = 0.85;
        const yTop = pad.t;
        const yBot = pad.t + plotH;
        for (const e of csIdx.bps) {
          if (e.mb < mbMin || e.mb > mbMax) continue;
          const xx = pad.l + ((e.mb - mbMin) / (mbMax - mbMin)) * plotW;
          if (!Number.isFinite(xx)) continue;
          ctx.beginPath();
          ctx.moveTo(xx + 0.5, yTop);
          ctx.lineTo(xx + 0.5, yBot);
          ctx.stroke();
        }
        ctx.restore();
      }
    } catch (_) { /* fail-soft */ }
  }
}

// --- drawPCA() — legacy lines 35749-35949 ---
export function drawPCA(state) {
  const canvas = document.getElementById('pcaCanvas');
  const { ctx, w, h } = fitCanvas(canvas);
  ctx.clearRect(0, 0, w, h);
  if (!state.data) return;
  document.getElementById('emptyState').style.display = 'none';
  const d = state.data;
  const cur = state.cur;
  const trailStart = Math.max(0, cur - state.trailN);

  // v3.25: which two PCs to plot (default PC1×PC2). PC1 keeps its sign-flip
  // rule (signX); other PCs render in raw orientation. The analytics path
  // (k-means, sample identity, lockedLabels) still uses canonical PC1×PC2
  // via getPC() — only the visual axis changes here.
  const [axisX, axisY] = state.viewControls.pcaXY;

  // Range across trail span
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (let wi = trailStart; wi <= cur; wi++) {
    const { x: xs, y: ys, signX, signY } = getPCRender(wi, axisX, axisY);
    if (!xs || !ys) continue;
    const samples = (wi === cur) ? allSampleIdx() : state.tracked;
    for (const si of samples) {
      const x = xs[si] * signX, y = ys[si] * signY;
      if (x < xMin) xMin = x; if (x > xMax) xMax = x;
      if (y < yMin) yMin = y; if (y > yMax) yMax = y;
    }
  }
  const xPad = (xMax - xMin) * 0.08 || 0.01;
  const yPad = (yMax - yMin) * 0.08 || 0.01;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;

  // v3.73: top padding tightened from 32 → 10 — there's no chart title or
  // axis label above the plot frame, just the YMax tick label which sits at
  // pad.t+4 inside the frame. The previous 32px reserved blank space that
  // pushed the data-point cluster down and squeezed the visible area. Bottom
  // padding 38 stays (X-axis tick labels + "PC1" label sit there).
  // v3.99 turn 14d ask 3: right padding reduced 200 → 16. Previously the
  // 200-px gutter reserved space for an external legend, but the per-band
  // legend is rendered INSIDE the plot frame at the top-right (line ~14391
  // computes lx = pad.l + plotW - legendW - 3), so the right gutter was
  // pure empty space. Reducing it lets the scatter use the full canvas
  // width. Hit-detection in onPCAClick() must use the same value.
  const pad = { l: 50, r: 16, t: 10, b: 38 };
  const plotW = w - pad.l - pad.r, plotH = h - pad.t - pad.b;
  const toX = v => pad.l + ((v - xMin) / (xMax - xMin)) * plotW;
  const toY = v => pad.t + (1 - (v - yMin) / (yMax - yMin)) * plotH;

  // Frame + grid
  ctx.strokeStyle = themeColor('rule');
  ctx.strokeRect(pad.l + 0.5, pad.t + 0.5, plotW, plotH);
  ctx.strokeStyle = 'rgba(42,50,66,0.5)';
  ctx.setLineDash([2, 3]);
  for (let i = 1; i < 5; i++) {
    const gx = pad.l + (i / 5) * plotW;
    ctx.beginPath(); ctx.moveTo(gx, pad.t); ctx.lineTo(gx, pad.t + plotH); ctx.stroke();
    const gy = pad.t + (i / 5) * plotH;
    ctx.beginPath(); ctx.moveTo(pad.l, gy); ctx.lineTo(pad.l + plotW, gy); ctx.stroke();
  }
  ctx.setLineDash([]);

  // Axis labels
  // v4 turn 2 ask 5: include %variance contribution in the axis labels.
  // For local PCA, each window has lam1 and lam2 (the two eigenvalues of
  // the local covariance matrix). The %variance shown is the fraction
  // each eigenvalue contributes to the top-2-PC basis: lam_k / (lam1+lam2).
  // This is locally interpretable — "of the variance captured by PC1+PC2,
  // PC1 explains XX%" — which is the only honest percentage we can give
  // from per-window data alone (full-dimensional variance is not in JSON).
  // Falls back to bare "PC1"/"PC2" when eigenvalues are missing or invalid.
  // Quentin: "you add (%) variance in axis text in y and x bc now it just
  // says PC1 which is not enough."
  let _pcaPC1Label = 'PC1';
  let _pcaPC2Label = 'PC2';
  if (state.data && state.data.windows && state.cur >= 0
      && state.cur < state.data.windows.length) {
    const _wObj = state.data.windows[state.cur];
    const _l1 = _wObj && _wObj.lam1;
    const _l2 = _wObj && _wObj.lam2;
    if (_l1 != null && isFinite(_l1) && _l2 != null && isFinite(_l2)
        && (_l1 + _l2) > 1e-12) {
      const _sum = _l1 + _l2;
      const _p1 = (100 * _l1 / _sum).toFixed(1);
      const _p2 = (100 * _l2 / _sum).toFixed(1);
      _pcaPC1Label = `PC1 (λ₁ ${_p1}%)`;
      _pcaPC2Label = `PC2 (λ₂ ${_p2}%)`;
    }
  }
  ctx.fillStyle = themeColor('ink-dim');
  ctx.font = '10px ui-monospace, monospace';
  ctx.textAlign = 'center';
  ctx.fillText(_pcaPC1Label, pad.l + plotW / 2, h - 12);
  ctx.save(); ctx.translate(14, pad.t + plotH / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText(_pcaPC2Label, 0, 0); ctx.restore();
  ctx.textAlign = 'right';
  ctx.fillText(yMax.toFixed(3), pad.l - 5, pad.t + 4);
  ctx.fillText(yMin.toFixed(3), pad.l - 5, pad.t + plotH);
  ctx.textAlign = 'center';
  ctx.fillText(xMin.toFixed(3), pad.l, pad.t + plotH + 14);
  ctx.fillText(xMax.toFixed(3), pad.l + plotW, pad.t + plotH + 14);

  // Determine L2 group labels for cluster coloring.
  // Priority: state.lockedLabels (manually frozen) > current L2's K-means labels.
  const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
  let groupLabels = null;
  if (state.colorMode === 'cluster') {
    if (state.lockedLabels) {
      groupLabels = state.lockedLabels;
    } else if (curL2 >= 0) {
      const cl = getL2Cluster(curL2);
      if (cl && cl.labels) groupLabels = cl.labels;
    }
  }

  // Non-tracked samples (current window) — use user-selected PC axes
  const cR = getPCRender(cur, axisX, axisY);
  const pc1c = cR.x, pc2c = cR.y, signXc = cR.signX, signYc = cR.signY;
  const trackedSet = new Set(state.tracked);
  // Cache current-window screen-space positions for lasso readback.
  // Stored as Float32Array pairs (x, y) indexed by sample idx.
  const _pcaScreenXY = new Float32Array(d.n_samples * 2);
  if (pc1c && pc2c) {
    for (let si = 0; si < d.n_samples; si++) {
      const x = toX(pc1c[si] * signXc);
      const y = toY(pc2c[si] * signYc);
      _pcaScreenXY[si * 2]     = x;
      _pcaScreenXY[si * 2 + 1] = y;
      if (trackedSet.has(si)) continue;
      const baseCol = getSampleColor(si, state.colorMode, groupLabels);
      const col = withAlpha(baseCol, state.colorMode === 'cluster' ? 0.45 : 0.7);
      ctx.fillStyle = col;
      ctx.beginPath(); ctx.arc(x, y, 2.8, 0, Math.PI * 2); ctx.fill();
    }
  }
  // Expose for lasso handler. Plot bounds also stored so lasso can clip
  // its drag rectangle to the data area.
  state.__pcaScreenXY = _pcaScreenXY;
  state.__pcaPlotRect = { x: pad.l, y: pad.t, w: plotW, h: plotH };

  // Trails (tracked samples)
  if (state.trailOn && state.tracked.length > 0 && state.trailN > 0) {
    // turn 128: when neither PCA axis is PC1, the trail lines wander in 2D
    // (they're no longer ~horizontal as in the default pc1×pc2 view because
    // PC2/3/4 vary a lot per window). Drop the line alpha so the cluttered
    // wander doesn't visually dominate the scatter. The dot fade below
    // already implies temporal direction; the line is just connective tissue.
    const _trailNeitherIsPc1 = (axisX !== 'pc1' && axisY !== 'pc1');
    const _trailLineAlpha = _trailNeitherIsPc1 ? 0.32 : 0.85;
    for (const si of state.tracked) {
      const col = trackedColor(si);
      ctx.strokeStyle = col; ctx.lineWidth = 1.3;
      ctx.globalAlpha = _trailLineAlpha;
      ctx.beginPath();
      let first = true;
      for (let wi = trailStart; wi <= cur; wi++) {
        const r = getPCRender(wi, axisX, axisY);
        if (!r.x || !r.y) continue;
        const x = toX(r.x[si] * r.signX), y = toY(r.y[si] * r.signY);
        if (first) { ctx.moveTo(x, y); first = false; } else ctx.lineTo(x, y);
      }
      ctx.stroke();
      for (let wi = trailStart; wi < cur; wi++) {
        const r = getPCRender(wi, axisX, axisY);
        if (!r.x || !r.y) continue;
        const t = (wi - trailStart) / Math.max(1, cur - trailStart);
        ctx.fillStyle = col; ctx.globalAlpha = 0.15 + 0.5 * t;
        const x = toX(r.x[si] * r.signX), y = toY(r.y[si] * r.signY);
        ctx.beginPath(); ctx.arc(x, y, 2.0, 0, Math.PI * 2); ctx.fill();
      }
    }
    ctx.globalAlpha = 1;
  }

  // Tracked dots — outline color = group, fill = sample identity, label = CGA
  if (pc1c && pc2c) for (const si of state.tracked) {
    const col = trackedColor(si);
    const x = toX(pc1c[si] * signXc), y = toY(pc2c[si] * signYc);
    // Group halo
    if (groupLabels) {
      const gcol = ['#4fa3ff', '#b8b8b8', '#f5a524'][groupLabels[si]] || '#666';
      ctx.strokeStyle = gcol; ctx.lineWidth = 3;
      ctx.beginPath(); ctx.arc(x, y, 9, 0, Math.PI * 2); ctx.stroke();
    }
    ctx.fillStyle = col;
    ctx.strokeStyle = themeColor('bg');
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.arc(x, y, 6, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
    // Label with text shadow for readability
    const name = d.samples[si].cga || d.samples[si].ind;
    ctx.font = 'bold 11px ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.strokeStyle = themeColor('bg'); ctx.lineWidth = 3;
    ctx.strokeText(name, x + 9, y + 4);
    ctx.fillStyle = col;
    ctx.fillText(name, x + 9, y + 4);
  }
  // turn 120: refresh the scree inset on every PCA draw. Cheap (pure SVG
  // string write to an absolutely-positioned div, no canvas, no layout).
  // The renderer handles the off/on toggle and the empty-state internally.
  try { _refreshScreeInset(); } catch (e) { /* fail-soft */ }
}

// --- drawAnchorStrip() — legacy lines 36147-36256 ---
export function drawAnchorStrip(state) {
  const canvas = document.getElementById('anchorStripCanvas');
  if (!canvas) return;
  const fit = fitCanvas(canvas);
  if (!fit) return;
  const { ctx, w, h } = fit;
  ctx.clearRect(0, 0, w, h);
  if (!state.data) return;
  const N = state.data.n_windows;
  // v3.59: defensive recompute. If samples are tracked but concord is missing
  // or stale (e.g. data just loaded, or some upstream caller forgot to call
  // recomputeAnchorConcord), do it inline so the strip never silently shows
  // the "track samples to enable" hint while tracking IS actually active.
  if (state.tracked.length > 0 &&
      (!state.anchorConcord ||
       state.anchorConcord.length !== N ||
       !state.trackingAnchor)) {
    if (typeof recomputeAnchorConcord === 'function') {
      try { recomputeAnchorConcord(); } catch (_) {}
    }
  }
  const concord = state.anchorConcord;
  // v3.58: padding matches drawLines (l:44, r:16) so the per-window x-coords
  // line up exactly with the PC1/PC2 plot area above. Previously was l:30,
  // r:8 which shifted the strip a few pixels left of the lines plot.
  const padL = 44, padR = 16;
  const plotW = Math.max(10, w - padL - padR);

  // Background
  ctx.fillStyle = themeColor('panel');
  ctx.fillRect(0, 0, w, h);

  if (state.tracked.length === 0) {
    // No tracked samples yet — render generic hint
    ctx.fillStyle = themeColor('ink-dimmer');
    ctx.font = '10px ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.fillText('anchor concord — track samples to enable', 8, h / 2 + 3);
    return;
  }
  if (!state.trackingAnchor) {
    // v3.59: tracked but anchor failed to capture (no L2 envelopes on the
    // chrom, or some other edge case). Better diagnostic than the generic hint.
    ctx.fillStyle = themeColor('ink-dimmer');
    ctx.font = '10px ui-monospace, monospace';
    ctx.textAlign = 'left';
    const noL2 = !Array.isArray(state.data.l2_envelopes) ||
                 state.data.l2_envelopes.length === 0;
    ctx.fillText(noL2
      ? 'anchor concord — no L2 envelopes on this chromosome'
      : 'anchor concord — anchor capture failed (try re-tracking inside an L2 zone)',
      8, h / 2 + 3);
    return;
  }
  if (!concord) {
    // Recompute returned with concord=null somehow. Render gray.
    ctx.fillStyle = themeColor('ink-dimmer');
    ctx.font = '10px ui-monospace, monospace';
    ctx.textAlign = 'left';
    ctx.fillText('anchor concord — computing…', 8, h / 2 + 3);
    return;
  }

  // Draw one vertical strip per window
  const stripH = h - 4;
  const stripY = 2;
  for (let wi = 0; wi < N; wi++) {
    const x0 = padL + Math.floor((wi / N) * plotW);
    const x1 = padL + Math.floor(((wi + 1) / N) * plotW);
    ctx.fillStyle = _vColor(concord[wi]);
    ctx.fillRect(x0, stripY, Math.max(1, x1 - x0), stripH);
  }

  // Y-axis label on the left
  ctx.fillStyle = themeColor('ink-dim');
  ctx.font = '9px ui-monospace, monospace';
  ctx.textAlign = 'left';
  ctx.fillText('V', 4, h / 2 + 3);

  // Vertical orange line at current scrubber position
  const curX = padL + Math.floor(((state.cur + 0.5) / N) * plotW);
  ctx.strokeStyle = '#f5a524';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(curX + 0.5, 0);
  ctx.lineTo(curX + 0.5, h);
  ctx.stroke();

  // Vertical green line at anchor position
  if (state.trackingAnchor.winIdx != null && state.trackingAnchor.winIdx >= 0) {
    const ax = padL + Math.floor(((state.trackingAnchor.winIdx + 0.5) / N) * plotW);
    ctx.strokeStyle = 'rgba(0,230,118,0.85)';
    ctx.lineWidth = 1.2;
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    ctx.moveTo(ax + 0.5, 0);
    ctx.lineTo(ax + 0.5, h);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  // Current V value in tiny text on the right
  const vCur = concord[state.cur];
  if (isFinite(vCur)) {
    ctx.fillStyle = themeColor('ink-dim');
    ctx.font = '9.5px ui-monospace, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(`V=${vCur.toFixed(2)}`, w - 4, h / 2 + 3);
  }
}

// --- renderL3Panel() — legacy lines 48685-49193 ---
export function renderL3Panel(state) {
  // v4.1: install delegated click handler for the L3 metric-cycle chip,
  // once per document lifetime. We listen on document.body and check the
  // event target's class so we don't re-attach on every L3 re-render.
  if (!window._v41HandlerInstalled && typeof document !== 'undefined') {
    window._v41HandlerInstalled = true;
    document.body.addEventListener('click', function (ev) {
      // Find the closest .l3-metric-chip (event might fire on inner span)
      let el = ev.target;
      while (el && el !== document.body) {
        if (el.classList && el.classList.contains('l3-metric-chip')) break;
        el = el.parentNode;
      }
      if (!el || el === document.body) return;
      // Cycle through cramer → nmi → ami → ari → cramer
      const order = ['cramer', 'nmi', 'ami', 'ari'];
      const labels = { cramer: "Cramér's V", nmi: 'NMI', ami: 'AMI', ari: 'ARI' };
      const tips = {
        cramer: 'Effect size of association in contingency table; N-independent. Does NOT account for label permutations — answers "are these variables associated?" not "are these the same partition?"',
        nmi: 'Normalized Mutual Information (Strehl-Ghosh, geometric-mean variant). [0,1]. Standard for partition comparison; permutation-invariant. Slightly sensitive to imbalanced cluster sizes.',
        ami: 'Adjusted Mutual Information (Vinh-Epps-Bailey 2010). Corrects NMI for chance agreement under hypergeometric null. Best metric for imbalanced K-means partitions like the hatchery 180/30/16 split.',
        ari: 'Adjusted Rand Index (Hubert-Arabie 1985). [-1,1] with 0 = chance. Counts agreeing sample-pairs across both partitions, adjusted for expectation under random labelling. Robust to cluster-size imbalance.',
      };
      // Different bands per metric
      function bandFor(key, v) {
        const t = { cramer: [0.10, 0.30, 0.50], nmi: [0.30, 0.50, 0.75],
                    ami: [0.20, 0.40, 0.65], ari: [0.20, 0.40, 0.65] }[key];
        if (v >= t[2]) return 'large';
        if (v >= t[1]) return 'medium';
        if (v >= t[0]) return 'small';
        return 'none';
      }
      const cur = el.dataset.metricActive || 'cramer';
      const nextIdx = (order.indexOf(cur) + 1) % order.length;
      const next = order[nextIdx];
      const v = parseFloat(el.dataset[next]);
      // Update state + persist
      state.l3SecondaryMetric = next;
      try { localStorage.setItem('pca_scrubber_v3.l3SecondaryMetric', next); } catch (_) {}
      // Update chip in-place — read the inner spans by querySelector
      el.dataset.metricActive = next;
      el.title = `${tips[next]}\n\nClick to cycle: Cramér's V → NMI → AMI → ARI → Cramér's V`;
      const lblSpan = el.querySelector('.lbl');
      const valSpan = el.querySelector('.val');
      if (lblSpan) lblSpan.innerHTML = `${labels[next]} <span class="dim" style="font-size:9px; opacity:0.6;">▸</span>`;
      if (valSpan) valSpan.innerHTML = `${v.toFixed(2)} <span class="dim" style="font-size:10px;">${bandFor(next, v)}</span>`;
      ev.stopPropagation();
    }, true);
    // Also restore from localStorage on first render
    try {
      const saved = localStorage.getItem('pca_scrubber_v3.l3SecondaryMetric');
      if (saved && ['cramer', 'nmi', 'ami', 'ari'].includes(saved)) {
        state.l3SecondaryMetric = saved;
      }
    } catch (_) {}
    // v4.1: wire the L3 re-clustering mode picker. Restore from localStorage,
    // set the dropdown's value, attach change listener that persists. The
    // handler does NOT yet trigger a re-render of the contingency — that's
    // post-v4.1 work. Today the dropdown is a state slot only.
    try {
      const savedMode = localStorage.getItem('pca_scrubber_v3.l3ReclusterMode');
      const validModes = ['kmeans-K3', 'kmeans-K6', 'distance-uv', 'uv-rotated',
                          'uv-denoise', 'uv-dbscan', 'uv-dist-rank', 'uv-dist-fuzzy'];
      if (savedMode && validModes.includes(savedMode)) {
        state.l3ReclusterMode = savedMode;
      }
    } catch (_) {}
    const reSel = document.getElementById('l3ReclusterSel');
    if (reSel) {
      reSel.value = state.l3ReclusterMode;
      reSel.addEventListener('change', function () {
        const v = reSel.value;
        // Only accept enabled modes; disabled options shouldn't reach here
        // (the disabled attr blocks selection in browsers) but guard anyway.
        const validModes = ['kmeans-K3', 'kmeans-K6', 'distance-uv', 'uv-rotated',
                            'uv-denoise', 'uv-dbscan', 'uv-dist-rank', 'uv-dist-fuzzy'];
        if (validModes.includes(v)) {
          state.l3ReclusterMode = v;
          try { localStorage.setItem('pca_scrubber_v3.l3ReclusterMode', v); } catch (_) {}
          // v4.1 continue: trigger an L3 re-render so the contingency rebuilds
          // against the new mode. Three call sites pick up the change:
          // (1) the focal-vs-neighbor contingency (line ~14328 in
          // _renderL3Fingerprint via the offset != 0 branch), (2) the
          // conservation calculation (in the focal pane's stats), (3) any
          // tracked-samples readouts that read from cmp.table. Cache hits
          // for the new mode are O(1) after first compute, so re-renders
          // are fast even on chromosomes with many L2s.
          if (typeof renderL3Panel === 'function') {
            renderL3Panel();
          }
        }
      });
    }
  }

  // v4 turn 32: install per-pane toolbar delegation once. Fires for clicks
  // on K-mode / color-mode buttons and changes on the recluster select
  // inside any of the L3 panes' header toolbars.
  if (typeof _wireL3PaneToolsDelegation === 'function') {
    try { _wireL3PaneToolsDelegation(); } catch (_) {}
  }

  const d = state.data;
  if (!d) return;

  // v4 turn 12a (Deliverable D): scale-stability mode short-circuits both
  // the L2 path AND the slab path. When state.l3Mode === 'scale_stability'
  // the carousel shows three different SCALES (fine / medium / coarse) on
  // the same focal window, with cross-scale fuse/split detection.
  if (state.l3Mode === 'scale_stability') {
    return renderL3PanelScaleStability();
  }

  // v3.45: slab mode short-circuits the L2-based path. compareUnit ≠ 'L2'
  // means each pane represents a slab of W windows, not an L2 envelope.
  // Sub-step A: render only the focal slab pane; neighbor panes show
  // a placeholder (sub-step B will wire them).
  if (state.compareUnit && state.compareUnit !== 'L2') {
    return renderL3PanelSlab();
  }

  // v3.99 turn 7 perf: short-circuit when nothing that affects pane content
  // has changed since the last render. The L3 panes are a pure function of
  // (curL2, layout, color mode, K mode, K, agg, merge τ, α, min n/grp,
  //  lockedLabels, l3Detailed, secondaryL2, candidate-draft state, tracked
  //  samples for the detailed readout, bcScope). state.cur within the same
  //  L2 doesn't affect anything — the panes show envelope-level clusters,
  //  not cursor-relative state. Skipping the rebuild on same-L2 steps avoids
  //  ~5 mini-PCA canvas redraws + the DOM teardown of body.innerHTML='' +
  //  full pane reconstruction. Big payoff for window-by-window scrubbing.
  const cur = state.cur;
  const curL2 = state.windowToL2 ? state.windowToL2[cur] : -1;
  const fp = _renderL3Fingerprint(curL2);
  if (state._l3CacheFp === fp && state._l3CacheRendered === true) {
    return;   // same render — DOM unchanged
  }
  state._l3CacheFp = fp;
  state._l3CacheRendered = true;

  const layoutKey = state.l3Layout || 'leftright';
  const layout = L3_LAYOUTS[layoutKey] || L3_LAYOUTS.leftright;
  const body = document.getElementById('l3Body');
  body.style.gridTemplateColumns = layout.cols;
  body.innerHTML = '';

  // Meta line
  const metaEl = document.getElementById('l3Meta');
  if (curL2 < 0) {
    metaEl.textContent = '— scrub into an L2 envelope —';
    body.innerHTML = '<div class="l3-col"><div class="dim" style="padding:12px;">scroll into an L2 envelope</div></div>';
    return;
  }
  metaEl.innerHTML = `K=${state.k} · agg=${state.aggMethod} · merge τ=${state.mergeThr.toFixed(2)} · α=${state.alpha.toFixed(3)} · min n/grp=${state.minNGroup}`;

  // ---- DUAL LAYOUT: two pinned L2 envelopes + middle comparison column ----
  if (layoutKey === 'dual' && state.secondaryL2 != null) {
    const aIdx = curL2;
    const bIdx = state.secondaryL2;
    if (aIdx < 0 || bIdx < 0 || aIdx === bIdx) {
      body.innerHTML = '<div class="l3-col"><div class="dim" style="padding:12px;">Dual layout needs two distinct L2 envelopes pinned. Scroll to one, click 📌, then scroll to the other.</div></div>';
      return;
    }
    const aEnv = d.l2_envelopes[aIdx];
    const bEnv = d.l2_envelopes[bIdx];
    const aCl  = getL2Cluster(aIdx);
    const bCl  = getL2Cluster(bIdx);
    const cmpAB = compareL2Pair(aIdx, bIdx);

    // ---- Column A: focal A ----
    const colA = document.createElement('div');
    colA.className = 'l3-col';
    const h3A = document.createElement('h3');
    h3A.classList.add('focal');
    h3A.innerHTML = `🅰️ <b>${shortId(aEnv.candidate_id)}</b> <span class="dim" style="font-weight:400;">${aEnv.n_windows}W · sim ${fmt(aEnv.mean_sim)}</span>`;
    colA.appendChild(h3A);
    const miniA = document.createElement('canvas');
    miniA.className = 'mini-pca';
    colA.appendChild(miniA);
    const contentA = document.createElement('div');
    contentA.className = 'ct-content';
    contentA.innerHTML = focalContentHtml(aCl, aEnv, aIdx);
    colA.appendChild(contentA);
    body.appendChild(colA);

    // ---- Column M (middle): A vs B contingency ----
    const colM = document.createElement('div');
    colM.className = 'l3-col';
    const h3M = document.createElement('h3');
    h3M.innerHTML = `🅰️ vs 🅱️ <span class="dim" style="font-weight:400;">contingency</span>`;
    colM.appendChild(h3M);
    // No mini for middle
    const contentM = document.createElement('div');
    contentM.className = 'ct-content';
    contentM.innerHTML = ctHtml(cmpAB, +1);  // offset +1 just for arrow display
    colM.appendChild(contentM);
    body.appendChild(colM);

    // ---- Column B: focal B (pinned) ----
    const colB = document.createElement('div');
    colB.className = 'l3-col';
    const h3B = document.createElement('h3');
    h3B.classList.add('focal');
    h3B.innerHTML = `🅱️ <b>${shortId(bEnv.candidate_id)}</b> <span class="dim" style="font-weight:400;">${bEnv.n_windows}W · sim ${fmt(bEnv.mean_sim)}</span>`;
    colB.appendChild(h3B);
    const miniB = document.createElement('canvas');
    miniB.className = 'mini-pca';
    colB.appendChild(miniB);
    const contentB = document.createElement('div');
    contentB.className = 'ct-content';
    contentB.innerHTML = focalContentHtml(bCl, bEnv, bIdx);
    colB.appendChild(contentB);
    body.appendChild(colB);

    // Draw both mini PCAs after layout, both Hungarian-aligned to A so colors
    // are consistent across the two-column strip (matches the rest of the UI).
    // In the A/B dual-pinned layout we treat A as offset=0 (focal) and B as
    // offset=+1 for palette purposes, so the user can still distinguish panes
    // when l3ColorMode is 'unique' or 'dual'.
    requestAnimationFrame(() => {
      drawMiniPCA(miniA, aIdx, aCl ? aCl.labels : null,
        { paneOffset: 0, focalLabels: null });
      const aCl_full = aCl;
      const aFocalLabels = aCl_full ? (aCl_full.fixedKLabels || aCl_full.labels) : null;
      drawMiniPCA(miniB, bIdx, alignedLabelsTo(aIdx, bIdx),
        { paneOffset: +1, focalLabels: aFocalLabels });
    });
    return;
  }

  // Build one column per offset
  // For 'dual' mode we need focal-aligned labels per neighbor (i.e. how the
  // focal would label this neighbor's samples). That's exactly what
  // alignLabels(focal_lab, neighbor_lab) gives us — the inverse perm — so we
  // pass each neighbor its own labels for fill, and the focal labels (cf)
  // re-keyed to neighbor sample index for the ring. Since labels are arrays
  // indexed by sample index in BOTH L2s (sample identity is preserved across
  // L2s), focal's labels work directly as the "what would focal say about
  // this sample" reference for any pane.
  //
  // l3KMode controls per-column section rendering:
  //   'k3' (default): one section per column at state.k (typically 3)
  //   'k6': one section per column forced to K=6
  //   'both': two sections per column — first at K=3, then K=6 — stacked
  const l3kMode = state.l3KMode || 'k3';
  // Decide which K(s) to render. Always include state.k as the primary K when
  // mode is 'k3' (so adaptive-K mode flows through correctly), use exactly 6
  // for 'k6', and both 3 and 6 for 'both'.
  const ksToRender = (l3kMode === 'k6') ? [6]
                   : (l3kMode === 'both') ? [state.k, 6]
                   : [state.k];

  layout.offsets.forEach((offset, colIdx) => {
    const isFocal = (offset === 0);
    const l2idx   = l2NeighborAt(curL2, offset);

    const col = document.createElement('div');
    col.className = 'l3-col';
    // v3.66: tag the column with a gold backdrop when this L2 is part of the
    // active candidate-mode draft. Solid gold = "all joins are MERGE";
    // dashed gold border = "draft includes a non-MERGE join" (warning).
    if (typeof l2DraftStatus === 'function') {
      const dStatus = l2DraftStatus(l2idx);
      if (dStatus === 'in')      col.classList.add('in-candidate');
      else if (dStatus === 'in-warn') col.classList.add('in-candidate-warn');
    }

    const h3 = document.createElement('h3');
    if (isFocal) h3.classList.add('focal');
    let titlePrefix = '';
    if (offset < 0) titlePrefix = (offset === -1 ? '←' : '←←');
    else if (offset > 0) titlePrefix = (offset === +1 ? '→' : '→→');
    else titlePrefix = '◆';
    // v4 turn 32: per-pane toolbar mirrors of the global L3 controls. Sits
    // on the right of the h3 (flex layout below). Tools render even on
    // empty panes (l2idx == null) so layout stays stable across panes.
    const paneToolsHtml = (typeof _l3PaneHeaderToolsHtml === 'function')
      ? _l3PaneHeaderToolsHtml() : '';
    if (l2idx == null) {
      h3.innerHTML = `<span class="l3-pane-title">${titlePrefix} <b>—</b></span>${paneToolsHtml}`;
    } else {
      const env = d.l2_envelopes[l2idx];
      h3.innerHTML = `<span class="l3-pane-title">${titlePrefix} <b>${shortId(env.candidate_id)}</b> <span class="dim" style="font-weight:400;">${env.n_windows}W · sim ${fmt(env.mean_sim)}</span></span>${paneToolsHtml}`;
    }
    col.appendChild(h3);

    if (l2idx == null) {
      // Empty placeholder column — single mini + message, no per-K split.
      const miniCanvas = document.createElement('canvas');
      miniCanvas.className = 'mini-pca';
      col.appendChild(miniCanvas);
      const content = document.createElement('div');
      content.className = 'ct-content';
      content.innerHTML = `<div class="dim" style="padding: 4px 0;">no neighbor at offset ${offset}<br>(boundary of L1 parent)</div>`;
      col.appendChild(content);
      drawMiniPCAEmpty(miniCanvas);
      body.appendChild(col);
      return;
    }

    // For each K in ksToRender, append a (mini-PCA + content) section.
    // Each section is independently rendered at its own K. In 'both' mode,
    // a small label badge identifies the K above each mini-PCA so the user
    // can tell K=3 from K=6 at a glance.
    //
    // v3.52: when only one K is rendered, append " · K=N" to the h3 caption
    // (saves a vertical line). Standalone badge only shown in multi-K mode
    // and only between sections (not before the first one).
    if (ksToRender.length === 1) {
      const K = ksToRender[0];
      // Append ' · K=N' to the existing h3 inline
      const dimSpan = h3.querySelector('span.dim');
      if (dimSpan) {
        dimSpan.innerHTML += ` · K=${K}`;
      }
    }
    // v3.72: in K=both mode the invariant meta block (windows/span/SNPs/density)
    // is the same for K=3 and K=6 — render it once at the top of the focal
    // column, then pass skipInvariantMeta:true to focalContentHtml in each
    // per-K section so it doesn't repeat. Only do this for the focal column
    // (offset === 0) since neighbor columns don't render the meta block.
    const hoistInvariantMeta = (ksToRender.length > 1) && isFocal;
    if (hoistInvariantMeta) {
      const env_focal = d.l2_envelopes[l2idx];
      // Use any K's cluster to read cl.nW (it's K-independent — number of
      // windows in the L2). Prefer state.k to keep semantics with old code.
      const clForStats = getL2Cluster(l2idx) || getL2ClusterAt(l2idx, ksToRender[0]);
      if (clForStats) {
        const stats = _l2InvariantStats(clForStats, env_focal, l2idx);
        const sharedHeader = document.createElement('div');
        sharedHeader.style.cssText = 'padding: 4px 0; border-bottom: 1px dashed var(--rule); margin-bottom: 2px;';
        sharedHeader.innerHTML = _invariantMetaInlineHtml(stats);
        col.appendChild(sharedHeader);
      }
    }
    ksToRender.forEach((K, sectionIdx) => {
      // K-label badge: only when multi-K rendering AND not the first section
      // (the first section's K is implied by the h3 caption append above —
      // wait, in multi-K we don't append to h3 — so first section needs the
      // badge too in multi-K mode). Show the badge in multi-K mode for all
      // sections. Use minimal vertical space.
      if (ksToRender.length > 1) {
        const kBadge = document.createElement('div');
        kBadge.style.cssText = 'font-family: var(--mono); font-size: 9.5px; ' +
          'color: var(--ink-dim); margin: 2px 0 0 0; ' +
          'border-top: ' + (sectionIdx > 0 ? '1px dashed var(--rule)' : 'none') + '; ' +
          'padding-top: ' + (sectionIdx > 0 ? '4px' : '0') + ';';
        kBadge.textContent = `K=${K}`;
        col.appendChild(kBadge);
      }

      // Pick the labels and helpers based on K.
      // For K === state.k, use the existing primary cluster path (preserves
      // adaptive-K, family-purity, coherence, and other diagnostics seen in
      // focalContentHtml). For K !== state.k, fall through to the per-K path.
      const useExisting = (K === state.k);
      const focalCl = useExisting
        ? getL2Cluster(curL2)
        : getL2ClusterAt(curL2, K);
      const focalLabels_forRing = focalCl
        ? (focalCl.fixedKLabels || focalCl.labels)
        : null;

      // v3.99 turn 6: for FOCAL pane, render the karyotype chips (per group /
      // center PC1) ABOVE the mini-PCA, in a tight small-font row that
      // matches the K=N badge size (~9.5px). Previously these chips lived
      // BELOW the plot inside `content`, which forced the user's eye to
      // travel down past the scatter to read them. Per Quentin's request:
      // "the karyotype info should be above the pca plot not below". For
      // neighbor panes the chips don't apply — they show contingency
      // tables instead via ctHtml.
      let focalChipsAbove = null;
      if (isFocal) {
        const cl = useExisting ? getL2Cluster(l2idx) : getL2ClusterAt(l2idx, K);
        const chipsHtml = _kSpecificMetaInlineHtml(cl, l2idx);
        if (chipsHtml) {
          focalChipsAbove = document.createElement('div');
          // Tight small font, matches K=N badge (9.5px). Override
          // the .meta-inline default size via inline style on the wrapper.
          focalChipsAbove.style.cssText = 'font-size: 9.5px; line-height: 1.2; ' +
            'padding: 2px 0; margin: 0;';
          focalChipsAbove.innerHTML = chipsHtml;
          col.appendChild(focalChipsAbove);
        }
      }

      const miniCanvas = document.createElement('canvas');
      miniCanvas.className = 'mini-pca';
      col.appendChild(miniCanvas);
      // v4 turn 6: spotlight click — fires after drawMiniPCA so __l3_render
      // is populated before the user can interact.
      if (typeof _setupL3MiniClick === 'function') _setupL3MiniClick(miniCanvas);

      // v4 turn 37 (Ask C step 2): band-selector strip below the focal
      // mini-PCA. Visible only on the FOCAL pane and only when
      // state.candidateMode === true. Returns '' otherwise so non-focal
      // panes and out-of-cand-mode views are unchanged.
      if (isFocal && typeof _l3BandSelectorHtml === 'function') {
        try {
          // Pull per-band sample counts from the focal cluster for hover
          // tooltips. Fail-soft if the cluster object is missing.
          let bandCounts = null;
          try {
            const cl = useExisting ? getL2Cluster(l2idx) : getL2ClusterAt(l2idx, K);
            if (cl) {
              const labels = cl.fixedKLabels || cl.labels;
              if (labels && labels.length) {
                bandCounts = new Array(K).fill(0);
                for (let s = 0; s < labels.length; s++) {
                  const lbl = labels[s] | 0;
                  if (lbl >= 0 && lbl < K) bandCounts[lbl]++;
                }
              }
            }
          } catch (_) {}
          const selHtml = _l3BandSelectorHtml(K, bandCounts);
          if (selHtml) {
            const sel = document.createElement('div');
            sel.innerHTML = selHtml;
            // unwrap one level so the .l3-band-selector div is the direct
            // child of col (no extra wrapping div)
            const inner = sel.firstChild;
            if (inner) col.appendChild(inner);
          }
          if (typeof _wireBandSelectorClicks === 'function') {
            _wireBandSelectorClicks();
          }
        } catch (_) { /* fail-soft — never break focal pane */ }
      }

      const content = document.createElement('div');
      content.className = 'ct-content';
      col.appendChild(content);

      if (isFocal) {
        // Focal column: K-means details + focal mini-PCA
        const cl = useExisting ? getL2Cluster(l2idx) : getL2ClusterAt(l2idx, K);
        const env = d.l2_envelopes[l2idx];
        // v3.72: when in K=both mode skip the invariant meta in each K's
        // section (it was hoisted to the top of the column above)
        // v3.99 turn 6: also skip kSpecific meta — already rendered above
        // the mini-PCA via focalChipsAbove.
        content.innerHTML = focalContentHtml(cl, env, l2idx,
          { skipInvariantMeta: hoistInvariantMeta, skipKSpecificMeta: true });
        const alignedLabels = useExisting
          ? alignedLabelsTo(curL2, l2idx)
          : alignedLabelsTo_atK(curL2, l2idx, K);
        requestAnimationFrame(() => drawMiniPCA(miniCanvas, l2idx,
          alignedLabels,
          { paneOffset: 0, focalLabels: null }));
      } else {
        // Neighbor column: contingency vs FOCAL + neighbor mini-PCA.
        // v4.1 continue: when state.l3ReclusterMode is set to a non-default
        // mode (kmeans-K6 or distance-uv), route through compareL2Pair_byMode
        // which uses the alternative clustering. Default 'kmeans-K3' falls
        // through to the existing useExisting / _atK branch — no behavioral
        // change for users who don't touch the picker.
        const reMode = state.l3ReclusterMode || 'kmeans-K3';
        let cmp;
        if (reMode !== 'kmeans-K3') {
          cmp = (offset < 0)
            ? compareL2Pair_byMode(l2idx, curL2, reMode)
            : compareL2Pair_byMode(curL2, l2idx, reMode);
        } else {
          cmp = useExisting
            ? ((offset < 0) ? compareL2Pair(l2idx, curL2) : compareL2Pair(curL2, l2idx))
            : ((offset < 0) ? compareL2Pair_atK(l2idx, curL2, K) : compareL2Pair_atK(curL2, l2idx, K));
        }
        // v3.65: compute alignedLabels BEFORE the ctHtml call so we can pass
        // them in (used by tracked-samples readout — neighbor's labels need
        // to be in focal's K-frame for cross-pane comparability).
        const alignedLabels = useExisting
          ? alignedLabelsTo(curL2, l2idx)
          : alignedLabelsTo_atK(curL2, l2idx, K);
        content.innerHTML = ctHtml(cmp, offset, alignedLabels);
        // v4 turn 6: stash spotlight-resolution metadata on the column so
        // _applySpotlightHighlights can find the (r, c) cell for each
        // sample without re-running the alignment. cmp.perm is the
        // permutation `bestPerm` from alignLabels — bestPerm[r] = right
        // raw label that aligns to left raw label r; invPerm maps a raw
        // right label to its aligned position. We compute invPerm here
        // to keep the post-processor cheap.
        const _invPerm = new Array(K);
        if (cmp && cmp.perm) {
          for (let r = 0; r < K; r++) _invPerm[cmp.perm[r]] = r;
        }
        col.__l3_meta = {
          K, offset,
          leftIdx:  (offset < 0) ? l2idx : curL2,
          rightIdx: (offset < 0) ? curL2 : l2idx,
          invPerm: _invPerm,
        };
        requestAnimationFrame(() => drawMiniPCA(miniCanvas, l2idx,
          alignedLabels,
          { paneOffset: offset, focalLabels: focalLabels_forRing }));
      }
    });

    body.appendChild(col);
  });
  // v4 turn 6: apply cross-pane sample-spotlight highlights to the just-
  // rendered contingency tables. Walks each .l3-col, finds the cell where
  // state.spotlight (and/or state.tracked when spotlightTrackedAll is on)
  // sits, and tags it with the spotlight CSS class. PCA dot rendering
  // is handled inside drawMiniPCA — drawMiniPCA reads state.spotlight
  // directly, no post-processing needed there.
  if (typeof _applySpotlightHighlights === 'function') {
    try { _applySpotlightHighlights(curL2); } catch (_) {}
  }
}

// --- renderL3PanelSlab() — legacy lines 49209-49478 ---
export function renderL3PanelSlab(state) {
  const d = state.data;
  if (!d) return;
  const body = document.getElementById('l3Body');
  if (!body) return;
  const layoutKey = state.l3Layout || 'leftright';
  const layout = L3_LAYOUTS[layoutKey] || L3_LAYOUTS.leftright;
  body.style.gridTemplateColumns = layout.cols;
  body.innerHTML = '';

  const halfW = compareUnitHalfW();
  const cur = state.cur;
  const range = (halfW != null) ? slabRange(cur, halfW) : null;
  const W = (halfW != null) ? (2 * halfW + 1) : 0;

  // v3.55: honor state.l3KMode in slab mode (parity with L2 mode).
  // 'k3' → render at state.k; 'k6' → render at K=6; 'both' → render K=3 then K=6 stacked.
  const l3kMode = state.l3KMode || 'k3';
  const ksToRender = (l3kMode === 'k6') ? [6]
                   : (l3kMode === 'both') ? [state.k, 6]
                   : [state.k];

  // Toolbar meta — show all Ks rendered when multiple
  const metaEl = document.getElementById('l3Meta');
  if (metaEl) {
    if (range) {
      const startMb = d.windows[range[0]].center_mb;
      const endMb   = d.windows[range[1]].center_mb;
      const kStr = ksToRender.length === 1 ? `K=${ksToRender[0]}` : `K=${ksToRender.join('+')}`;
      metaEl.innerHTML =
        `slab=${W}w · windows ${range[0]+1}–${range[1]+1} · ` +
        `${startMb.toFixed(2)}–${endMb.toFixed(2)} Mb · ` +
        `${kStr} · merge τ=${state.mergeThr.toFixed(2)}`;
    } else {
      metaEl.textContent = '— invalid slab —';
    }
  }
  if (!range) {
    body.innerHTML = '<div class="l3-col"><div class="dim" style="padding:12px;">invalid slab range</div></div>';
    return;
  }

  // Render each pane defined by the layout. v3.55: each pane iterates over
  // ksToRender, producing one (mini-PCA + content) section per K. In single-K
  // mode this matches the v3.46 behavior; in multi-K mode the sections stack
  // vertically with a thin separator.
  layout.offsets.forEach((offset, paneIdx) => {
    const isFocal = (offset === 0);
    const col = document.createElement('div');
    col.className = 'l3-col';

    // Pane header (one per pane, K-agnostic)
    const head = document.createElement('div');
    head.style.cssText = 'padding: 6px 10px; font-family: var(--mono); font-size: 11px;' +
                        ' display: flex; gap: 8px; align-items: center;' +
                        ' border-bottom: 1px solid var(--rule);';
    let offsetSlab = null;
    // turn 148: per-pane toolbar mirrors of the global L3 controls (parity
    // with L2 mode line 44592). Sits to the right of the title text via
    // margin-left:auto inside the existing .l3-pane-tools span. Renders even
    // on out-of-range panes so layout stays stable across panes.
    const paneToolsHtml = (typeof _l3PaneHeaderToolsHtml === 'function')
      ? _l3PaneHeaderToolsHtml() : '';
    if (isFocal) {
      // Inline K=N suffix when single-K (v3.53 convention)
      const kSuffix = (ksToRender.length === 1) ? ` · K=${ksToRender[0]}` : '';
      head.innerHTML = `<span class="l3-pane-title" style="font-weight: 600;">◆ slab focal` +
                       ` <span class="dim" style="font-weight:400;">w${range[0]+1}–w${range[1]+1} (${W}w)${kSuffix}</span></span>` +
                       paneToolsHtml;
    } else {
      const arrow = offset < 0 ? '←' : '→';
      offsetSlab = slabRangeOffset(cur, halfW, offset);
      const kSuffix = (ksToRender.length === 1) ? ` · K=${ksToRender[0]}` : '';
      const lblText = offsetSlab
        ? `slab ${offset > 0 ? '+' : ''}${offset} · w${offsetSlab[0]+1}–w${offsetSlab[1]+1}${kSuffix}`
        : `slab ${offset} · out of range`;
      head.innerHTML = `<span class="l3-pane-title"><span class="dim">${arrow}</span> ` +
                       `<span style="font-weight: 500;">${lblText}</span></span>` +
                       paneToolsHtml;
    }
    col.appendChild(head);

    // Edge cases (neighbor only) — handle once before the K loop
    if (!isFocal) {
      if (!offsetSlab) {
        const oot = document.createElement('div');
        oot.className = 'dim';
        oot.style.cssText = 'padding: 16px 12px; font-size: 11px;';
        oot.textContent = 'out of range';
        col.appendChild(oot);
        body.appendChild(col);
        return;
      }
      if (offsetSlab[0] === offsetSlab[1] && offsetSlab[0] === range[0] && range[0] === range[1]) {
        const oot = document.createElement('div');
        oot.className = 'dim';
        oot.style.cssText = 'padding: 16px 12px; font-size: 11px;';
        oot.textContent = 'identical to focal at edge';
        col.appendChild(oot);
        body.appendChild(col);
        return;
      }
    }

    // For each K, append (badge + mini-PCA + content) to the column
    ksToRender.forEach((K, sectionIdx) => {
      // K-badge: only in multi-K mode (single-K's K is in the head already)
      if (ksToRender.length > 1) {
        const kBadge = document.createElement('div');
        kBadge.style.cssText = 'font-family: var(--mono); font-size: 9.5px; ' +
          'color: var(--ink-dim); padding: 4px 10px 0 10px; ' +
          'border-top: ' + (sectionIdx > 0 ? '1px dashed var(--rule)' : 'none') + '; ' +
          'margin-top: ' + (sectionIdx > 0 ? '4px' : '0') + ';';
        kBadge.textContent = `K=${K}`;
        col.appendChild(kBadge);
      }

      // turn 148: karyotype chips ABOVE the mini-PCA on the focal pane
      // (parity with L2 mode line 44687-44699). _kSpecificMetaInlineHtml works
      // off cl.n_per_group + cl.centers — both available on slab clusters.
      // Pass l2idx=null to skip the sub-band annotation chip (which only
      // applies to L2 envelopes in active draft context).
      if (isFocal) {
        const clFocal = getSlabClusterAt(range[0], range[1], K);
        if (typeof _kSpecificMetaInlineHtml === 'function') {
          const chipsHtml = _kSpecificMetaInlineHtml(clFocal, null);
          if (chipsHtml) {
            const focalChipsAbove = document.createElement('div');
            focalChipsAbove.style.cssText = 'font-size: 9.5px; line-height: 1.2; ' +
              'padding: 2px 10px; margin: 0;';
            focalChipsAbove.innerHTML = chipsHtml;
            col.appendChild(focalChipsAbove);
          }
        }
      }

      // Mini-PCA canvas
      const miniWrap = document.createElement('div');
      miniWrap.style.cssText = 'padding: 6px;';
      const mini = document.createElement('canvas');
      mini.style.cssText = 'display: block; width: 100%; height: 110px; cursor: crosshair; background: var(--panel-2);';
      miniWrap.appendChild(mini);
      col.appendChild(miniWrap);
      // turn 148: spotlight click — fires after drawSlabMiniPCA has populated
      // canvas.__l3_render. Same handler as L2 mode (parity with line 44707).
      if (typeof _setupL3MiniClick === 'function') _setupL3MiniClick(mini);

      // turn 148: band-selector strip on focal pane in candidate mode (parity
      // with L2 mode line 44713-44744). _l3BandSelectorHtml is K-and-counts
      // only — no L2 dependency — so it works for slabs too.
      if (isFocal && typeof _l3BandSelectorHtml === 'function') {
        try {
          let bandCounts = null;
          try {
            const cl = getSlabClusterAt(range[0], range[1], K);
            if (cl) {
              const labels = cl.fixedKLabels || cl.labels;
              if (labels && labels.length) {
                bandCounts = new Array(K).fill(0);
                for (let s = 0; s < labels.length; s++) {
                  const lbl = labels[s] | 0;
                  if (lbl >= 0 && lbl < K) bandCounts[lbl]++;
                }
              }
            }
          } catch (_) {}
          const selHtml = _l3BandSelectorHtml(K, bandCounts);
          if (selHtml) {
            const sel = document.createElement('div');
            sel.innerHTML = selHtml;
            const inner = sel.firstChild;
            if (inner) col.appendChild(inner);
          }
          if (typeof _wireBandSelectorClicks === 'function') {
            _wireBandSelectorClicks();
          }
        } catch (_) { /* fail-soft — never break focal pane */ }
      }

      // Content block: focal shows cluster info; neighbors show contingency
      const content = document.createElement('div');
      content.className = 'ct-content';
      content.style.cssText = 'padding: 6px 10px;';
      col.appendChild(content);

      if (isFocal) {
        const cl = getSlabClusterAt(range[0], range[1], K);
        content.innerHTML = slabFocalContentHtml(cl, range, K);
        requestAnimationFrame(() => drawSlabMiniPCA(mini, range, cl ? cl.labels : null));
      } else {
        // ctHtml expects focal-on-rows for offset>=0, focal-on-cols for offset<0.
        // turn 149: when state.l3ReclusterMode is set to a non-default mode,
        // route through compareSlabPair_byMode (parity with L2 line 44832).
        // The byMode variant handles kmeans-K6 directly and falls back to
        // kmeans-K3 with a `fellBack: true` flag for U/V modes (which don't
        // yet have slab implementations). Default kmeans-K3 keeps the simple
        // path so unaffected users see no behavioral change.
        const reMode = state.l3ReclusterMode || 'kmeans-K3';
        let cmp;
        if (reMode !== 'kmeans-K3') {
          cmp = (offset < 0)
            ? compareSlabPair_byMode(offsetSlab, range, reMode)
            : compareSlabPair_byMode(range, offsetSlab, reMode);
        } else {
          cmp = (offset < 0)
            ? compareSlabPair(offsetSlab, range, K)
            : compareSlabPair(range, offsetSlab, K);
        }
        // v3.65: compute alignedLabels first so we can pass to ctHtml for
        // tracked-samples readout (slab variant). When reMode forces a
        // different K (e.g. kmeans-K6 against l3KMode='k3' sections), the
        // aligned labels are still at the section K — visual K mismatch
        // between the mini-PCA palette (K=3) and the contingency (K=6) is
        // the same trade-off L2 mode accepts.
        const alignedLabels = alignedSlabLabelsTo(range, offsetSlab, K);
        const neighborCl = getSlabClusterAt(offsetSlab[0], offsetSlab[1], K);
        if (!cmp) {
          content.innerHTML = '<div class="dim">low power · slab K-means failed</div>';
        } else {
          // turn 149: surface fallback notice when a U/V mode was requested
          // but slab compare fell back to kmeans-K3 (no slab-aware U/V
          // implementation yet). Prepended to the contingency content so
          // the user knows the visible table doesn't reflect the dropdown
          // setting.
          let prefix = '';
          if (cmp.fellBack && cmp.requestedMode) {
            prefix = `<div class="dim" style="font-size:9.5px; line-height:1.3; ` +
                     `padding:2px 4px; margin-bottom:4px; ` +
                     `border-left: 2px solid var(--accent); padding-left:6px;" ` +
                     `title="The U/V recluster modes (uv-rotated, uv-denoise, uv-dbscan, uv-dist-rank, uv-dist-fuzzy) require a slab-aware rotation pipeline that isn't shipped yet. ` +
                     `Slab compare falls back to kmeans-K3 for these modes. Use kmeans-K3 or kmeans-K6 for accurate slab-mode reclustering, or switch the compare unit back to L2 for U/V modes.">` +
                     `↩ slab fallback: <b>${cmp.requestedMode}</b> not available for slabs · using kmeans-K3` +
                     `</div>`;
          }
          content.innerHTML = prefix + ctHtml(cmp, offset, alignedLabels);
        }
        // turn 148: stash spotlight-resolution metadata on the column (parity
        // with L2 mode line 44798-44807). Slab variant uses leftRange /
        // rightRange instead of L2 indices and sets isSlab:true so the
        // post-pass can route to getSlabClusterAt.
        // turn 149: meta.K must reflect the K that was actually used in the
        // contingency, which may differ from the section K when reMode forced
        // a different K (e.g. kmeans-K6). Read from cmp.K.
        const ctK = (cmp && cmp.K) ? cmp.K : K;
        const _invPerm = new Array(ctK);
        if (cmp && cmp.perm) {
          for (let r = 0; r < ctK; r++) _invPerm[cmp.perm[r]] = r;
        }
        col.__l3_meta = {
          K: ctK, offset,
          isSlab: true,
          leftRange:  (offset < 0) ? offsetSlab : range,
          rightRange: (offset < 0) ? range      : offsetSlab,
          invPerm: _invPerm,
        };
        requestAnimationFrame(() => drawSlabMiniPCA(mini, offsetSlab,
          alignedLabels || (neighborCl ? neighborCl.labels : null)));
      }
    });

    body.appendChild(col);
  });

  // turn 148: apply cross-pane sample-spotlight highlights to the just-rendered
  // contingency tables (parity with L2 mode line 44822). _applySpotlightHighlights
  // detects slab columns via meta.isSlab and routes to getSlabClusterAt.
  if (typeof _applySpotlightHighlights === 'function') {
    try { _applySpotlightHighlights(null); } catch (_) {}
  }
}

// --- renderL3PanelScaleStability() — legacy lines 12833-12912 ---
export function renderL3PanelScaleStability(state) {
  const d = state.data;
  const body = document.getElementById('l3Body');
  const metaEl = document.getElementById('l3Meta');
  if (!body) return;

  if (!d) {
    body.innerHTML = '<div class="l3-col"><div class="dim" style="padding:12px;">no data</div></div>';
    return;
  }

  const result = _scaleStabilityCompute();
  if (metaEl) {
    metaEl.innerHTML = 'scale-stability mode · focal window ' + (state.cur + 1) +
      ' · verdict <b>' + result.verdict + '</b>';
  }

  // Layout: 3 panes + 2 inter-pane "gap" columns (numeric contingency tables)
  // Total 5 columns: pane | gap | pane | gap | pane
  body.style.gridTemplateColumns =
    'minmax(180px,1fr) minmax(120px,160px) minmax(180px,1fr) minmax(120px,160px) minmax(180px,1fr)';
  body.innerHTML = '';

  // Append a column helper
  function appendCol(html, classes) {
    const col = document.createElement('div');
    col.className = (classes || 'l3-col');
    col.innerHTML = html;
    body.appendChild(col);
    return col;
  }

  for (let i = 0; i < 3; i++) {
    const cfg  = state.scaleStabilityPanes[i];
    const pane = result.panes[i] || { ok: false, why: 'missing' };
    const label = _scaleStabilityPaneLabel(i, cfg, pane);
    let paneHtml = '<div class="ss-pane-header">' +
                   '<span class="ss-pane-icon">' + ['Σ₁','Σ₂','Σ₃'][i] + '</span>' +
                   '<span>' + label.replace(/^Σ[₁₂₃]\s/, '') + '</span>';
    if (pane.lowConf) paneHtml += '<span class="ss-low-conf">low conf</span>';
    paneHtml += '</div>';
    if (!pane.ok) {
      paneHtml += '<div class="dim" style="padding:12px;font-family:var(--mono);font-size:11px;">— ' +
                  (pane.why || 'unavailable') + ' —</div>';
    } else {
      const range = pane.range ? ('w' + (pane.range[0] + 1) + '–w' + (pane.range[1] + 1)) : '';
      const counts = pane.n_per_group ? pane.n_per_group.join(' · ') : '';
      paneHtml += '<div style="padding:10px 12px;font-family:var(--mono);font-size:11px;line-height:1.6;">' +
                  '<div class="dim">range: ' + range + '</div>' +
                  '<div class="dim">cluster sizes: ' + counts + '</div>' +
                  '</div>';
    }
    appendCol(paneHtml);

    // Gap (between this pane and the next, but not after the last)
    if (i < 2) {
      const pw = result.pairwise[i];
      let gapHtml = '<div class="ss-gap-label">' + ['Σ₁→Σ₂','Σ₂→Σ₃'][i] + '</div>';
      // v4 turn 12b: SVG Sankey ribbons sit ABOVE the numeric contingency
      // table. Visual headline (fan-in / parallel ribbons / reshuffle) reads
      // first, the numeric backup follows.
      if (pw && pw.ok) {
        gapHtml += '<div class="ss-sankey">' + _ssSankeyHtml(pw.contingency) + '</div>';
      }
      gapHtml += pw && pw.ok ? _ssContingencyTableHtml(pw.contingency) :
                               '<div class="ss-event-row" style="color:var(--ink-dimmer);">—</div>';
      gapHtml += _ssEventRowHtml(pw);
      appendCol(gapHtml, 'ss-gap');
    }
  }

  // Verdict row spans all columns
  const verdictRow = document.createElement('div');
  verdictRow.className = 'ss-verdict-row';
  const cssCls = _scaleStabilityVerdictCss(result.verdict);
  verdictRow.innerHTML =
    '<span class="ss-verdict-pill ' + cssCls + '">' + result.verdict + '</span>' +
    '<span class="ss-verdict-text">' + _scaleStabilityVerdictText(result.verdict) + '</span>';
  body.appendChild(verdictRow);
}

// --- updateWinLabel() — legacy lines 51738-51742 ---
export function updateWinLabel(state) {
  const w = state.data.windows[state.cur];
  document.getElementById('winIdx').textContent = state.cur;
  document.getElementById('winBp').innerHTML = `· <b>${w.center_mb.toFixed(3)} Mb</b>`;
}

// --- setCur() — legacy lines 51747-51792 ---
export function setCur(state, i) {
  if (!state.data) return;
  state.cur = Math.max(0, Math.min(state.data.n_windows - 1, i | 0));
  document.getElementById('scrubber').value = state.cur;
  updateWinLabel();
  drawSim(); drawZ(); drawTracks(); drawLinesPanel(); drawPCA();
  // v3.51: keep the minimap's orange crosshair in sync with the scrubber
  if (state.simInMinimap && typeof drawSimMini === 'function') {
    try { drawSimMini(); } catch (_) {}
  }
  // v3.52: anchor concord strip — orange cursor line follows scrubber
  if (typeof drawAnchorStrip === 'function') {
    try { drawAnchorStrip(); } catch (_) {}
  }
  // v3.71: concord V badge (above per-sample-lines header) follows the scrubber
  if (typeof _updateConcordBadge === 'function') {
    try { _updateConcordBadge(); } catch (_) {}
  }
  updateSidebarInfo();
  renderZoneBlock();
  renderL3Panel();
  // v3.94: live dosage heatmap follows cursor (debounced; no-op when closed)
  if (typeof redrawCursorHeatmap === 'function') {
    try { redrawCursorHeatmap(); } catch (_) {}
  }
  // v4 turn 1 ask 1: keep the Windows page table .cur highlight in sync
  // with the scrubber even when the user is on page 1 (or any other page).
  // Auto-scroll the highlighted row into view + brief orange flash.
  // No-ops cheaply if the windows page DOM hasn't been built yet, or if
  // the current window's row is filtered out / not in the visible 2000.
  if (typeof _refreshWinSumCurRow === 'function') {
    try { _refreshWinSumCurRow(); } catch (_) {}
  }
  // v4 turn 31: redraw the windows-page strip too so its orange cursor
  // line follows the scrubber when ←/→ arrow keys move state.cur. Cheap
  // no-op when the strip canvas isn't visible / hasn't been built.
  if (typeof drawWinSumStrip === 'function') {
    try { drawWinSumStrip(); } catch (_) {}
  }
  // v4 turn 22: keep the boundaries-page repeat density panel cursor in
  // sync with the global scrubber cursor. No-op when not on the boundaries
  // page or when no repeat density is loaded for the current chrom.
  if (typeof _renderRepeatDensityPanel === 'function') {
    try { _renderRepeatDensityPanel(); } catch (_) {}
  }
}

// --- autoPickRadial() — legacy lines 51862-51892 ---
export function autoPickRadial(state, n) {
  if (!state.data) return;
  if (n == null || !isFinite(n)) n = state.trackedN;
  n = Math.max(0, Math.min(state.data.n_samples, n | 0));
  if (n === 0) { state.tracked = []; renderTrackedList(); drawLinesPanel(); drawPCA(); renderL3Panel(); return; }
  const { pc1, pc2, sign } = getPC(state.cur);
  const N = state.data.n_samples;
  let mx = 0, my = 0;
  for (let i = 0; i < N; i++) { mx += pc1[i] * sign; my += pc2[i]; }
  mx /= N; my /= N;
  const picks = []; const used = new Set();
  for (let k = 0; k < n; k++) {
    const target = (k / n) * Math.PI * 2;
    let bestI = -1, bestScore = -Infinity;
    for (let i = 0; i < N; i++) {
      if (used.has(i)) continue;
      const dx = pc1[i] * sign - mx, dy = pc2[i] - my;
      const r = Math.sqrt(dx * dx + dy * dy);
      const ang = Math.atan2(dy, dx);
      let dA = Math.abs(ang - target);
      if (dA > Math.PI) dA = 2 * Math.PI - dA;
      const score = r - 0.5 * r * (dA / Math.PI);
      if (score > bestScore) { bestScore = score; bestI = i; }
    }
    if (bestI >= 0) { picks.push(bestI); used.add(bestI); }
  }
  state.tracked = picks;
  renderTrackedList();
  drawPCA();
  renderL3Panel();
}

// --- applyData() — legacy lines 54476-54690 ---
export function applyData(state, data) {
  // Cross-chrom safety: candidate from a different chromosome is meaningless.
  // Clear it before swapping data so refreshes downstream see no stale state.
  if (state.candidate && state.candidate.chrom !== data.chrom) {
    state.candidate = null;
  }
  // Schema detection (v3.20). state.layersPresent is the source of truth
  // for conditional UI rendering across the rest of the scrubber.
  const det = detectSchemaAndLayers(data);
  state.schemaVersion = det.schemaVersion;
  state.layersPresent = det.layers;
  console.log(`[scrubber] loaded ${data.chrom} · schema v${det.schemaVersion} · layers: [${Array.from(det.layers).sort().join(', ')}]`);
  state.data = data;
  // turn 130 Slice 2: lineage compute is per-chromosome (the concordance
  // matrix and lineage labels only mean something within one chrom). When
  // data swaps, drop the cached result so the next paint re-triggers
  // compute on the new chromosome's L2 inventory.
  if (typeof invalidateLineageCache === 'function') {
    try { invalidateLineageCache(); } catch (_) {}
  }
  // turn 161: same for the band-trace cache (per-chromosome, per-fish-set).
  // Also drop the fish-set itself — sample indices are per-chrom and
  // generally don't transfer to a new chromosome's data shape.
  if (typeof _bandTraceClearCache === 'function') {
    try { _bandTraceClearCache(); } catch (_) {}
  }
  state.bandTraceFishSet = null;
  state._lineageComputeScheduled = false;
  // v3.99 turn 14e: if simInMinimap was restored from localStorage at
  // startup but the visual .active class was held back (because data
  // wasn't loaded yet), reconcile now that data has landed.
  if (typeof _reapplyMinimapActiveOnDataLoad === 'function') {
    _reapplyMinimapActiveOnDataLoad();
  }
  // v3.99 t14e+ continue: refresh page 12 (θπ scrubber) layer-status
  // indicators in the empty-state. Each row's status flips from ⚪ to 🟢
  // when that layer is detected.
  if (typeof _refreshThetaPiLayerStatus === 'function') {
    try { _refreshThetaPiLayerStatus(); } catch (_) {}
  }
  // v4 turn 132 Slice 2: also flip panel visibility — empty-state hides
  // and per-layer panels reveal as their required layers arrive.
  if (typeof _refreshThetaPiPanelVisibility === 'function') {
    try { _refreshThetaPiPanelVisibility(); } catch (_) {}
  }
  // v4 turn 132 Slice 3: paint the CUSUM hero panel from cusum_theta if
  // present. Visibility wiring above already revealed/hid the panel; this
  // draws into its canvases. Other renderers (sim_mat, |Z|, lines, PCA,
  // L3) ship in later slices.
  if (typeof _drawThCusumHero === 'function') {
    try { _drawThCusumHero(); } catch (_) {}
  }
  // v4 turn 132 Slice 5: paint the per-sample θπ lines panel from
  // theta_pi_per_window if present. Single-source (no PC1/PC2/GHSL/het
  // stacking like page 1), no lasso, no caching — minimum viable mirror.
  if (typeof _drawThLinesPanel === 'function') {
    try { _drawThLinesPanel(); } catch (_) {}
  }
  // v4 turn 132 Slice 6a/6b: paint sim_mat heatmap + |Z| waveform from
  // theta_pi_local_pca if present. Mirrors page 1's drawSim/drawZ
  // minimum-viable subset — no L1/L2 overlays (need theta_pi_envelopes),
  // no click-to-jump, no PDF-style triangle split.
  if (typeof _drawThSimMatPanel === 'function') {
    try { _drawThSimMatPanel(); } catch (_) {}
  }
  if (typeof _drawThZPanel === 'function') {
    try { _drawThZPanel(); } catch (_) {}
  }
  // v4 turn 132 Slice 7a/7b: paint envelope anchor strip + PC1×PC2 scatter
  // from theta_pi_envelopes / theta_pi_local_pca.
  if (typeof _drawThAnchorStripPanel === 'function') {
    try { _drawThAnchorStripPanel(); } catch (_) {}
  }
  if (typeof _drawThPcaPanel === 'function') {
    try { _drawThPcaPanel(); } catch (_) {}
  }
  if (typeof _refreshGhslLayerStatus === 'function') {
    try { _refreshGhslLayerStatus(); } catch (_) {}
  }
  // Load saved candidate list for this chromosome from localStorage.
  // Each chromosome has its own list (cross-chrom labels are meaningless).
  if (typeof loadCandidateList === 'function') loadCandidateList();
  // turn 133 Slice 1 follow-up: chrom-load hook for L2-sweep.
  // Cache key is chrom-prefixed so the previous chrom's result wouldn't
  // re-serve, but explicit invalidation is cleaner. Then if the toggle
  // is on, run sweep + auto-promote against THIS chrom's L2 inventory
  // (which is now fresh in state.data, with candidates already loaded).
  // Order matters: candidates first → sweep can dedupe against them →
  // auto-promotes land in candidateList. Without this hook, switching
  // chrom with the toggle on would leave the sweep stale until the user
  // toggled it off-and-on.
  if (typeof invalidateL2SweepCache === 'function') {
    try { invalidateL2SweepCache(); } catch (_) {}
  }
  if (state.l2SweepEnabled
      && typeof runL2SweepInheritance === 'function'
      && typeof _autoPromoteFromSweep === 'function') {
    try {
      const sweepRes = runL2SweepInheritance({ force: true });
      if (sweepRes) _autoPromoteFromSweep(sweepRes);
    } catch (e) {
      console.warn('[l2sweep] applyData hook failed:', e && e.message);
    }
  }
  // turn 134 Slice 2: refresh inspector trigger-button availability after
  // the chrom-load sweep (or its absence). The button is disabled until
  // state.l2SweepResult is populated.
  if (typeof _refreshL2SweepInspectBtnAvailability === 'function') {
    try { _refreshL2SweepInspectBtnAvailability(); } catch (_) {}
  }
  // v3.99 turn 10: load saved catalogue favorites for this chromosome.
  // catState may not exist yet at very first call (it's defined later in the
  // file), so guard with typeof.
  if (typeof catState !== 'undefined' && typeof _loadCatFavorites === 'function') {
    catState.favorites = _loadCatFavorites();
  }
  // Load saved view controls (PCA axis selection etc.) and reconcile against
  // the data we just loaded — drops PC3/PC4 if not available, keeps PC1×PC2.
  loadViewControls();
  reconcileViewControlsForData();
  state.cur = 0;
  state.tracked = [];
  state.ancestryPalette = {};
  state.l2GroupCache = null;
  state.cacheKey = null;
  // v3.99 turn 7 perf: clear render caches whenever a new dataset loads
  if (typeof _l3CacheInvalidate === 'function') _l3CacheInvalidate();
  if (typeof _linesCacheInvalidate === 'function') _linesCacheInvalidate();
  buildIndexes();
  computePC1Signs();
  populateSimScales();
  buildFamilyPalette();
  refreshColorModeBar();
  refreshBandPickBar();
  refreshPcaAxisBar();
  // Manual groups: reload from localStorage now that we know the chrom
  if (typeof _mgRefreshOnDataLoad === 'function') _mgRefreshOnDataLoad();
  // v4 turn 128 (AS1): active samples — restore the saved CGA list for
  // this cohort and refresh the badge text. AS1 is purely scaffolding;
  // no other atlas function reads state.activeSampleSet yet.
  if (typeof loadActiveSamples === 'function') loadActiveSamples();
  if (typeof refreshActiveSamplesBadge === 'function') refreshActiveSamplesBadge();
  if (typeof buildLinesPanelCheckboxes === 'function') buildLinesPanelCheckboxes();
  if (typeof buildLinesPanel === 'function') buildLinesPanel();
  // v3.99 turn 14e+ continue: revalidate the lines coloring mode against
  // the layers we just discovered. If the user previously selected, e.g.,
  // 'theta_pi' on a different JSON that had the layer, but this JSON
  // doesn't, the picker falls back to 'kmeans' silently.
  if (typeof refreshLinesColorMode === 'function') refreshLinesColorMode();
  state.secondaryL2 = null;
  refreshPinUI();
  state.lockedLabels = null;
  state.lockedRefL2 = null;
  if (typeof refreshLockBtn === 'function') refreshLockBtn();
  buildTrackPanels();
  document.getElementById('scrubber').max = data.n_windows - 1;
  document.getElementById('scrubber').disabled = false;
  const l1c = (data.l1_envelopes || []).length;
  const l2c = (data.l2_envelopes || []).length;
  const bc  = (data.l2_boundaries || []).length;
  const scaleSummary = data.sim_scales && Object.keys(data.sim_scales).length > 0
    ? `<br><span class="dim">sim scales: ${Object.keys(data.sim_scales).join(', ')} (default: ${state.simScale})</span>`
    : '';
  let famSummary = '';
  if (data.family_source && data.family_source !== 'none') {
    const nHub = state.hubFamilies.length;
    const nSmall = state.smallFamilyIds.size;
    const nSing = state.singletonFamilyIds.size;
    const src = data.family_source === 'pairs' ? `pairs θ≥${data.theta_cutoff}` : 'natora --family';
    famSummary = `<br><span class="dim">families (${src}): ${nHub} hubs (n≥4), ${nSmall} small, ${nSing} singletons</span>`;
  }
  const sampleNames = data.samples && data.samples[0] && data.samples[0].cga !== data.samples[0].ind ? 'CGA' : 'Ind';
  const pc2Note = data.has_pc2 === false
    ? `<br><span class="dim" style="color: var(--accent);">PC2: jittered (slim precomp — no per-sample PC2)</span>`
    : '';
  const trackList = (data.tracks && Object.keys(data.tracks).length > 0)
    ? `<br><span class="dim">tracks: ${Object.keys(data.tracks).join(', ')}</span>`
    : '';
  document.getElementById('dataStatus').innerHTML =
    `<b>${data.chrom}</b><br>${data.n_windows} W · ${data.n_samples} samples (${sampleNames})<br>` +
    `<span class="dim">L1 envelopes: ${l1c}, L2 envelopes: ${l2c}, L2 peaks: ${bc}</span>` +
    scaleSummary + famSummary + pc2Note + trackList;
  // v3.99 t14e+ continue: surface optional species name in header. Used by
  // the manuscript prompt template too. Italicized as a scientific name
  // when present. Falls back to omission when absent.
  const speciesPart = (typeof data.species === 'string' && data.species.trim())
    ? ` · <i>${data.species.trim()}</i>`
    : '';
  document.getElementById('headerMeta').innerHTML =
    `<b>${data.chrom}</b>${speciesPart} · <b>${data.n_windows}</b> W · <b>${data.n_samples}</b> samples`;
  // Schema badge
  const schemaBadge = document.getElementById('schemaBadge');
  if (schemaBadge) {
    const layerNames = listLayers();
    schemaBadge.textContent = `schema v${state.schemaVersion} · ${layerNames.length} layer${layerNames.length === 1 ? '' : 's'}`;
    schemaBadge.title = `Schema version: v${state.schemaVersion}\nLayers: ${layerNames.join(', ')}\n\nUse + load enrichment to add layers from cluster phases 6+`;
    schemaBadge.className = 'v' + state.schemaVersion;
    schemaBadge.style.display = 'inline-block';
  }
  renderTrackedList();
  refreshCandidateUI();         // re-render page 2 metadata + tab dot
  if (typeof refreshCandidateListUI === 'function') refreshCandidateListUI();
  // v3.90: activate (or hide) the marker page based on whether a phase-13
  // marker layer was loaded. Re-render its content after activation so it
  // reflects whatever subset of {summary, catalogue, primers} arrived.
  if (typeof _refreshMarkerPageActivation === 'function') _refreshMarkerPageActivation();
  if (typeof renderMarkerPage === 'function') renderMarkerPage();
  setCur(0);
  // v4 turn 73f: persist this chromosome to IndexedDB so it survives page
  // reloads / cross-atlas navigation. Async, fire-and-forget; failures log
  // to console but don't block the UI.
  if (typeof _idbPersistChrom === 'function') {
    try { _idbPersistChrom(data); } catch (_) { /* fail-soft */ }
  }
}

// --- onSimClick() — legacy lines 52067-52111 ---
export function onSimClick(state, evt) {
  if (!state.data) return;
  const rect = document.getElementById('simCanvas').getBoundingClientRect();
  const px = evt.clientX - rect.left;
  const py = evt.clientY - rect.top;
  const g = state._simGeom;
  if (!g) return;
  // Click outside the heatmap square: ignore
  if (px < g.x0 || px > g.x1 || py < g.y0 || py > g.y1) return;
  // v4 turn 114c (remaining): cs-breakpoint click-to-jump on the diagonal
  // red-cross overlay (drawSim renders these via 114b). Check 2-D distance
  // to (toPx(bp.win), toPy(bp.win)) before falling back to the existing
  // x-axis scrub. Tolerance is _CS_BP_HIT_TOL_PX from the centralized
  // helper. Done before the x-axis fallback so a click ON a cross goes
  // to that cross, not to whatever window is at that x-fraction.
  if (typeof _ensureCsOverlayIndex === 'function') {
    const csIdx = _ensureCsOverlayIndex();
    if (csIdx && csIdx.bps.length > 0) {
      const Nw = state.data.n_windows;
      // Replicate drawSim's mapping exactly (it lives inside drawSim's
      // closure so we can't reuse it; re-derive from state._simGeom).
      const _toPx = (wIdx) => g.x0 + (wIdx + 0.5) * g.side / Nw;
      const _toPy = (wIdx) => g.y0 + (wIdx + 0.5) * g.side / Nw;
      // Allow a slightly larger tolerance here because the cross arms
      // extend ±7 px from center (see drawSim 114b armPx=7), so a click
      // near the tip of an arm is a legitimate target.
      const SIM_TOL = 7;
      const hit = _csBpHitTest2D(
        csIdx.bps, px, py,
        (bp) => (bp.win >= 0 && bp.win < Nw)
          ? { x: _toPx(bp.win), y: _toPy(bp.win) }
          : null,
        SIM_TOL
      );
      if (hit) {
        _csBpJumpToWindow(hit);
        return;
      }
    }
  }
  // Prefer x-axis (horizontal) → window index. The heatmap is symmetric, so
  // either axis works, but x is the natural "scrub through chromosome" gesture.
  const frac = (px - g.x0) / g.side;
  setCur(Math.round(frac * (state.data.n_windows - 1)));
}

// --- onZClick() — legacy lines 52112-52202 ---
export function onZClick(state, evt) {
  if (!state.data) return;
  const canvas = document.getElementById('zCanvas');
  const rect = canvas.getBoundingClientRect();
  const pad = { l: 44, r: 16 };
  const x = evt.clientX - rect.left;
  const y = evt.clientY - rect.top;     // v4 turn 10: track y for W-row hit-test
  const plotW = rect.width - pad.l - pad.r;
  const frac = Math.max(0, Math.min(1, (x - pad.l) / plotW));
  const d = state.data;
  const mbMin = d.windows[0].center_mb;
  const mbMax = d.windows[d.n_windows - 1].center_mb;
  const targetMb = mbMin + frac * (mbMax - mbMin);
  let bestI = 0, bestD = Infinity;
  for (let i = 0; i < d.n_windows; i++) {
    const dd = Math.abs(d.windows[i].center_mb - targetMb);
    if (dd < bestD) { bestD = dd; bestI = i; }
  }
  // v4 turn 10: if the click landed inside the W-row, treat it as a window-
  // resolution draft edit instead of a setCur jump. The W-row's y-band is
  // computed from the same layout constants drawZ uses.
  // We need to mirror drawZ's collapsed/expanded layout choice here:
  //   collapsed mode uses padT=2, candBarH=5, candGap=1, zoneH=14
  //   expanded mode uses padT=14, candBarH=7, candGap=2, zoneH=14
  // _wRowBand returns null when the W-row isn't visible, so we fall through
  // to setCur in those cases.
  // v4 turn 13 (Deliverable B): account for multi-lane candBar height
  // (candBarTotal = candBarH × n_lanes). The candidate-click hit-test runs
  // BEFORE the W-row test so clicks on stacked candidate bars route to
  // candidate selection rather than setCur.
  const collapsed = !!state.zCollapsed;
  const padT     = collapsed ? 2 : 14;
  const candBarH = collapsed ? 5 : 7;
  const candGap  = collapsed ? 1 : 2;
  const _candLanes = (typeof _assignCandidateLanes === 'function')
    ? _assignCandidateLanes(state.candidateList || []).n_lanes : 1;
  const candBarTotal = candBarH * _candLanes;
  const zoneTop  = padT + candBarTotal + candGap;
  const zoneH    = 14;
  // Mb→px function matching drawZ's
  const toX_click = (mb) => pad.l + ((mb - mbMin) / (mbMax - mbMin)) * plotW;
  // v4 turn 13: candidate hit-test (lane-aware) — sets active candidate on hit
  if (typeof _candidateAtClick === 'function') {
    const hit = _candidateAtClick(x, y, padT, candBarTotal, toX_click, d);
    if (hit) {
      state.candidate = hit;
      // v4 turn 56: route through helper so persistence stays consistent
      // with the prev/next nav path (_navigateToCandidate).
      _persistActiveCandidate(hit.id || '');
      // Trigger a re-render so the new active candidate's downstream views
      // (candidate-focus tab, scale-stability "candidate" scale) update.
      if (typeof drawZ === 'function') drawZ();
      if (typeof renderL3Panel === 'function') renderL3Panel();
      if (typeof renderCandidateMetadata === 'function') renderCandidateMetadata();
      return;
    }
  }
  // v4 turn 33 (Ask A): nav-lane click test runs BEFORE the W-row test.
  // Its band sits between the L2 zone bar and the W-row. The W-row's zoneH
  // is inflated by the nav-lane's height + gap so its hit-test lands below
  // the nav-lane (matching the drawZ layout).
  const _navBandClick = (typeof _winNavBand === 'function')
    ? _winNavBand({ collapsed, zoneTop, zoneH }) : null;
  const _navExtraClick = _navBandClick ? (_navBandClick.h + _navBandClick.gap) : 0;
  if (typeof _winNavHandleClick === 'function' &&
      _winNavHandleClick(y, bestI, { collapsed, zoneTop, zoneH })) return;
  if (_wRowHandleClick(y, bestI, { collapsed, zoneTop,
                                    zoneH: zoneH + _navExtraClick })) return;
  // v4 turn 114c (remaining): cs-breakpoint click-to-jump. The red dashed
  // vertical lines drawn by drawZ at toX_click(bp.mb) take priority over
  // the generic setCur(bestI) fallback so users can land on a breakpoint
  // window even when the click isn't pixel-perfect on the dashed line.
  // Hit-test runs AFTER candidate / nav-lane / W-row tests because those
  // are higher-precedence interactive zones; users dragging a candidate
  // boundary or W-row scrubbing shouldn't be hijacked by a cs-bp click.
  if (typeof _ensureCsOverlayIndex === 'function') {
    const csIdx = _ensureCsOverlayIndex();
    if (csIdx && csIdx.bps.length > 0) {
      const _toXBp = (bp) => {
        if (bp.mb < mbMin || bp.mb > mbMax) return NaN;
        return toX_click(bp.mb);
      };
      const hit = _csBpHitTestXFromList(csIdx.bps, x, _CS_BP_HIT_TOL_PX, _toXBp);
      if (hit) {
        _csBpJumpToWindow(hit);
        return;
      }
    }
  }
  setCur(bestI);
}

// --- onPCAClick() — legacy lines 52203-52244 ---
export function onPCAClick(state, evt) {
  if (!state.data) return;
  const canvas = document.getElementById('pcaCanvas');
  const rect = canvas.getBoundingClientRect();
  const px = evt.clientX - rect.left, py = evt.clientY - rect.top;
  const d = state.data, cur = state.cur;
  const trailStart = Math.max(0, cur - state.trailN);
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (let wi = trailStart; wi <= cur; wi++) {
    const { pc1, pc2, sign } = getPC(wi);
    const samples = (wi === cur) ? allSampleIdx() : state.tracked;
    for (const si of samples) {
      const x = pc1[si] * sign, y = pc2[si];
      if (x < xMin) xMin = x; if (x > xMax) xMax = x;
      if (y < yMin) yMin = y; if (y > yMax) yMax = y;
    }
  }
  const xPad = (xMax - xMin) * 0.08 || 0.01, yPad = (yMax - yMin) * 0.08 || 0.01;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;
  const w = rect.width, h = rect.height;
  // v3.73: must match drawPCA's pad constants. If they diverge, click-to-pick
  // hit detection will be off by the difference.
  // v3.99 turn 14d ask 3: r reduced 200 → 16 (matches drawPCA above).
  const pad = { l: 50, r: 16, t: 10, b: 38 };
  const plotW = w - pad.l - pad.r, plotH = h - pad.t - pad.b;
  if (px < pad.l || px > w - pad.r || py < pad.t || py > h - pad.b) return;
  const toX = v => pad.l + ((v - xMin) / (xMax - xMin)) * plotW;
  const toY = v => pad.t + (1 - (v - yMin) / (yMax - yMin)) * plotH;
  const { pc1, pc2, sign } = getPC(cur);
  let bestI = -1, bestD = Infinity;
  for (let si = 0; si < d.n_samples; si++) {
    const x = toX(pc1[si] * sign), y = toY(pc2[si]);
    const dd = (x - px) * (x - px) + (y - py) * (y - py);
    if (dd < bestD) { bestD = dd; bestI = si; }
  }
  if (bestI >= 0 && bestD < 400) {
    const i = state.tracked.indexOf(bestI);
    if (i >= 0) state.tracked.splice(i, 1);
    else if (state.tracked.length < state.trackedN) state.tracked.push(bestI);
    renderTrackedList(); drawLinesPanel(); drawPCA(); renderL3Panel();
  }
}

// --- togglePlay() — legacy lines 52399-52417 ---
export function togglePlay(state) {
  const btn = document.getElementById('playBtn');
  if (state.playing) {
    clearInterval(state.playTimer);
    state.playing = false;
    btn.textContent = '▶ Play';
    btn.classList.remove('playing');
  } else {
    state.playing = true;
    btn.textContent = '❚❚ Pause';
    btn.classList.add('playing');
    state.playTimer = setInterval(() => {
      if (!state.data) return;
      let next = state.cur + 2;
      if (next >= state.data.n_windows) next = 0;
      setCur(next);
    }, 80);
  }
}

// --- cycleKAside() — legacy lines 56537-56565 ---
export function cycleKAside(state) {
  const cur = state.k;
  const i = _K_CYCLE_ORDER.indexOf(cur);
  const next = (i < 0)
    ? _K_CYCLE_ORDER[0]
    : _K_CYCLE_ORDER[(i + 1) % _K_CYCLE_ORDER.length];
  state.kMode = 'fixed';
  state.k = next;
  state.l2GroupCache = null; state.cacheKey = null;
  // Mirror to sidebar kSelect so both UIs stay aligned
  const _kSel = document.getElementById('kSelect');
  if (_kSel) _kSel.value = String(next);
  if (typeof refreshBandPickBar === 'function') refreshBandPickBar();
  if (state.trackingAnchor) state.trackingAnchor = null;
  if (typeof recomputeAnchorConcord === 'function') {
    try { recomputeAnchorConcord(); } catch (_) {}
  }
  if (typeof drawPCA === 'function') drawPCA();
  if (typeof drawLinesPanel === 'function') drawLinesPanel();
  if (typeof renderZoneBlock === 'function') renderZoneBlock();
  if (typeof renderL3Panel === 'function') renderL3Panel();
  if (typeof drawAnchorStrip === 'function') {
    try { drawAnchorStrip(); } catch (_) {}
  }
  if (typeof _updateConcordBadge === 'function') {
    try { _updateConcordBadge(); } catch (_) {}
  }
  if (typeof _syncTrackedCompactUI === 'function') _syncTrackedCompactUI();
}

// --- renderTrackedList() — legacy lines 52004-52062 ---
export function renderTrackedList(state) {
  // v4 turn 73c: refactored to dual-write into both #trackedList (sidebar)
  // AND #trackedListCompact (compact panel). Same chip-building logic; each
  // chip is built fresh per container because a DOM node can only belong to
  // one parent.
  const elSidebar = document.getElementById('trackedList');
  const elCompact = document.getElementById('trackedListCompact');
  if (elSidebar) elSidebar.innerHTML = '';
  if (elCompact) elCompact.innerHTML = '';
  const ci = document.getElementById('trackedCountInfo');
  const ni = document.getElementById('trackedNInfo');
  if (ci) ci.textContent = state.tracked.length;
  if (ni) ni.textContent = state.trackedN;
  // v3.70: keep the compact-mode tracked-samples panel in sync
  if (typeof _syncTrackedCompactUI === 'function') {
    try { _syncTrackedCompactUI(); } catch (e) {}
  }
  // v3.71: refresh the concord badge after recomputeAnchorConcord runs below
  // (we call it after the recompute so the badge reflects the freshest data)
  // v3.52: tracked set changed — refresh anchor + recompute concord
  if (typeof recomputeAnchorConcord === 'function') {
    try { recomputeAnchorConcord(); } catch (e) {}
  }
  if (typeof drawAnchorStrip === 'function') {
    try { drawAnchorStrip(); } catch (e) {}
  }
  // v3.71: badge update after the recompute
  if (typeof _updateConcordBadge === 'function') {
    try { _updateConcordBadge(); } catch (e) {}
  }
  if (!state.data) return;
  // Pre-compute per-sample spread for the current L2 once.
  const curL2 = state.windowToL2 ? state.windowToL2[state.cur] : -1;
  const sd = curL2 >= 0 ? sampleSpreadL2(curL2) : null;

  // v4 turn 73c: helper to build a chip for one sample. Called twice (once
  // for sidebar, once for compact) so each container gets its own DOM node.
  const buildChip = (si) => {
    const s = state.data.samples[si];
    const chip = document.createElement('span');
    chip.className = 'tag on';
    let label = s.cga || s.ind;
    if (sd && isFinite(sd[si])) {
      label += ` <span style="opacity:0.6; font-size:10px;">σ${sd[si].toFixed(3)}</span>`;
      if (sd[si] > 0.05) chip.style.borderColor = 'var(--bad)';
    }
    chip.innerHTML = label;
    chip.onclick = () => {
      state.tracked = state.tracked.filter(x => x !== si);
      renderTrackedList(); drawLinesPanel(); drawPCA(); renderL3Panel();
    };
    return chip;
  };

  state.tracked.forEach(si => {
    if (elSidebar) elSidebar.appendChild(buildChip(si));
    if (elCompact) elCompact.appendChild(buildChip(si));
  });
}

// --- renderManualGroupsList() — legacy lines 48152-48191 ---
export function renderManualGroupsList(state) {
  const containers = [];
  if (typeof document !== 'undefined') {
    const sidebar = document.getElementById('manualGroupsList');
    const compact = document.getElementById('manualGroupsListCompact');
    // turn 135 Slice 1 (SPEC_g_panel_unified_groups.md): popup re-host.
    // The G-panel manual tab body holds a #manualGroupsListPopup div.
    // When the popup is open, this renderer also fills it so the
    // single source of truth (state.manualGroups) drives all three
    // surfaces. When the popup isn't open the lookup returns null
    // and the loop just skips it (existing pattern).
    const popup   = document.getElementById('manualGroupsListPopup');
    if (sidebar) containers.push(sidebar);
    if (compact) containers.push(compact);
    if (popup)   containers.push(popup);
  }
  if (containers.length === 0) return;
  const groups = state.manualGroups || [];
  if (groups.length === 0) {
    const emptyHtml = '<div class="mg-empty">No groups yet — pick samples or grab a K-band.</div>';
    for (const box of containers) box.innerHTML = emptyHtml;
    return;
  }
  let html = '';
  for (const g of groups) {
    const pinClass = g.scope === 'cohort' ? 'mg-pin pinned' : 'mg-pin';
    const pinTitle = g.scope === 'cohort'
      ? 'Pinned to cohort — follows you across chromosomes (click to unpin)'
      : 'Per-chromosome — click to pin to cohort';
    html += `<div class="mg-row" data-mgid="${g.id}">` +
      `<span class="mg-swatch" style="background:${g.color}"></span>` +
      `<span class="mg-name" contenteditable="true" spellcheck="false" ` +
      `data-mgid="${g.id}">${escapeHtml(g.name)}</span>` +
      `<span class="mg-count">n=${g.members.length}</span>` +
      `<button class="${pinClass}" data-mgid="${g.id}" title="${pinTitle}">📌</button>` +
      `<button class="mg-del" data-mgid="${g.id}" title="Remove this group">×</button>` +
      `</div>`;
  }
  for (const box of containers) box.innerHTML = html;
}

// --- buildLinesPanelCheckboxes() — legacy lines 33016-33117 ---
export function buildLinesPanelCheckboxes(state) {
  const wrap = document.getElementById('linesYsourceCheckboxes');
  if (!wrap || typeof wrap.appendChild !== 'function') return;
  if ('innerHTML' in wrap) wrap.innerHTML = '';
  const avail = availablePCs();   // ['pc1', 'pc2', ...]
  const sources = avail.slice();
  // Stage 3: append GHSL panel sources when the panel layer is loaded.
  // 'het' — alias for divergence at the primary scale (per-sample het rate
  //          per phased-snp denominator; from STEP_C04_snake3_ghsl_v6.R).
  // 'ghsl_div_s<scale>' — one per available rolling scale on the panel.
  if (state.data && state.data.ghsl_panel && state.data.ghsl_panel.div_roll) {
    const panel = state.data.ghsl_panel;
    sources.push('het');
    const scales = panel.scales || Object.keys(panel.div_roll);
    for (const scale of scales) {
      // scale comes through as 's10', 's20', etc. — keep the prefix as is.
      sources.push(`ghsl_div_${scale}`);
    }
  }
  const selected = new Set(state.viewControls.linesYsources);

  for (const src of sources) {
    const id = `linesYsrc_${src}`;
    const lab = document.createElement('label');
    lab.style.cssText = 'display: flex; align-items: center; gap: 4px; cursor: pointer;';
    const cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.id = id;
    cb.value = src;
    cb.checked = selected.has(src);
    cb.addEventListener('change', () => {
      const checkedNow = Array.from(wrap.querySelectorAll('input[type="checkbox"]:checked'))
        .map(el => el.value);
      // Always keep at least one source selected; if user unchecked the last
      // one, re-check the one they just toggled to keep the panel populated.
      if (checkedNow.length === 0) {
        cb.checked = true;
        return;
      }
      // Preserve order from the sources array (PCs first, then GHSL chips)
      const ordered = sources.filter(s => checkedNow.includes(s));
      setLinesYsources(ordered);
      // If linked, the PCA selector also updates — refresh its UI
      if (state.viewControls.linked) {
        if (typeof refreshPcaAxisBar === 'function') refreshPcaAxisBar();
        if (typeof drawPCA === 'function') drawPCA();
      }
      buildLinesPanel();
      drawLinesPanel();
    });
    lab.appendChild(cb);
    const span = document.createElement('span');
    // Keep PC labels short ('PC1'); GHSL labels readable ('het', 'div_s50')
    if (/^pc[1-4]$/.test(src)) {
      span.textContent = src.toUpperCase();
      // v3.99 turn 13 ask 4: per-PC tooltip explaining --npc availability so
      // the verbose visible note can be hidden. The text is built once we
      // know `avail.length` and `sources.length`, so we set the title
      // attribute on each PC label after the loop. For now, leave a
      // placeholder; the tooltip is patched below.
      lab.dataset.pcTooltip = '1';
    } else if (src === 'het') {
      span.textContent = 'het';
    } else if (/^ghsl_div_s\d+$/.test(src)) {
      span.textContent = `div_${src.replace(/^ghsl_div_/, '')}`;
    } else {
      span.textContent = src;
    }
    lab.appendChild(span);
    wrap.appendChild(lab);
  }

  // v3.99 turn 13 ask 4: build the npc-availability text once, then attach it
  // as a tooltip to every PC label (instead of a separate visible inline note
  // that pushed the toolbar to wrap onto a second line). The visible
  // #linesYsourceNote span is now hidden via display:none — kept in DOM so
  // older code paths that reference it don't crash.
  let npcTooltip = '';
  {
    const ghslExtras = sources.length - avail.length;
    if (ghslExtras > 0) {
      npcTooltip = `${avail.length} PCs · ${ghslExtras} GHSL sources in this dataset`;
    } else {
      npcTooltip = avail.length <= 2
        ? `${avail.length} PCs in this dataset (run with --npc 4 for more)`
        : `${avail.length} PCs available`;
    }
  }
  wrap.querySelectorAll('label[data-pc-tooltip]').forEach(lab => {
    lab.title = npcTooltip;
  });

  // Note about availability — v3.99 turn 13 ask 4: hidden by default. The
  // text is now in the per-PC label tooltips above. Element kept so the
  // initial-render code path at line ~19842 (which sets a different message
  // pre-buildLinesPanel) still finds the node without erroring.
  const note = document.getElementById('linesYsourceNote');
  if (note) {
    note.textContent = '';
    note.style.display = 'none';
  }
}

// --- buildLinesPanel() — legacy lines 33334-33610 ---
export function buildLinesPanel(state) {
  // Construct one sub-canvas per source in state.viewControls.linesYsources.
  // Recreates all sub-canvases on every call (cheap; sources rarely change).
  const container = document.getElementById('linesCanvasContainer');
  const panel = document.getElementById('linesPanel');
  if (!container || !panel) return;
  if (typeof container.appendChild !== 'function') return;   // test-shim safety
  if (!state.data) {
    if (panel.style) panel.style.display = 'none';
    return;
  }
  if (panel.style) panel.style.display = '';
  if ('innerHTML' in container) container.innerHTML = '';
  const sources = state.viewControls.linesYsources.slice();
  // v3.58: container flex/sizing is now in CSS (#linesPanel #linesCanvasContainer).
  // We only need to ensure display:flex flex-direction:column for sub-canvas
  // stacking; the parent #linesPanel is a flex column that gives this container
  // its share of the available height via min-height:0 + flex:1 1 auto.
  if (container.style) {
    container.style.display = 'flex';
    container.style.flexDirection = 'column';
  }
  for (const src of sources) {
    const sub = document.createElement('div');
    sub.className = 'lines-subpanel';
    sub.dataset.linesSource = src;
    // Each subpanel gets equal share via flex:1; min-height:0 lets it shrink
    sub.style.cssText = `position: relative; flex: 1 1 0; min-height: 0; border-bottom: 1px solid var(--rule, #2a3242);`;
    const cv = document.createElement('canvas');
    cv.style.cssText = 'display: block; width: 100%; height: 100%; cursor: crosshair;';
    cv.dataset.linesSource = src;
    sub.appendChild(cv);
    // turn 2p: install inheritance-pill tooltip on PC1 canvas. Pills are
    // only ever drawn on the PC1 sub-panel (per turn 2c), so we skip the
    // other sources to keep the handler load minimal.
    if (src === 'pc1' && typeof _wireInheritancePillTooltip === 'function') {
      try { _wireInheritancePillTooltip(cv); } catch (_) {}
    }
    // turn 162 — band-trace strip tooltip. Same gating as the inheritance
    // pill (PC1 only), since the strip itself is also PC1-only. Idempotent
    // via the canvas dataset marker, so re-running drawLinesPanel after
    // each chrom switch (which destroys + recreates these canvases) just
    // re-attaches handlers to fresh nodes.
    if (src === 'pc1' && typeof _wireBandTraceTooltip === 'function') {
      try { _wireBandTraceTooltip(cv); } catch (_) {}
    }
    const lbl = document.createElement('div');
    lbl.className = 'lines-subpanel-label';
    lbl.style.cssText = 'position: absolute; top: 4px; left: 8px; font-size: 11px; color: var(--dim, #888); pointer-events: none; font-family: ui-monospace, monospace; font-weight: 600;';
    lbl.textContent = src.toUpperCase();
    sub.appendChild(lbl);
    // Click-to-jump: same gesture as Z panel
    cv.addEventListener('click', e => {
      // v3.76: in lasso mode on the PC1 sub-canvas, the pointer handlers
      // below own the gesture; the click event still fires after a no-motion
      // pointerup, but we suppress it to avoid double-acting (the pointer-up
      // logic already clears any tiny rect). dataset.linesSource === 'pc1'
      // is the gating condition because lasso only lives on that source.
      if (state.linesLassoActive && cv.dataset.linesSource === 'pc1') {
        // If the down was treated as a click (no drag), let it fall through;
        // we set _lassoSwallowClick on pointerup when there was a drag.
        if (cv.__lassoSwallowClick) {
          cv.__lassoSwallowClick = false;
          e.preventDefault();
          e.stopPropagation();
          return;
        }
      }
      // v3.91: unmodified click on a ⚠ jumper marker (PC1 panel only) →
      // jump scrubber to the L2 where the fish first crossed K=3 parents
      // within the active candidate, then open the same popover the v3.87
      // shift+click path uses (anchored at the marker, fish forced to the
      // marker's si). The hit list `cv.__candJumperMarkers` is populated
      // during draw inside the candidate-overlay block. Hit radius 7 px.
      // No-op (falls through) when no marker hit, no candidate, or shift is
      // held (shift+click still routes to the existing line-inspect logic).
      if (!e.shiftKey && cv.dataset.linesSource === 'pc1' &&
          Array.isArray(cv.__candJumperMarkers) && cv.__candJumperMarkers.length > 0) {
        const rect0 = cv.getBoundingClientRect();
        const cx = e.clientX - rect0.left;
        const cy = e.clientY - rect0.top;
        const HIT_R = 7;
        let bestM = null, bestD2 = HIT_R * HIT_R;
        for (const m of cv.__candJumperMarkers) {
          const dx = m.x - cx, dy = m.y - cy;
          const d2 = dx * dx + dy * dy;
          if (d2 <= bestD2) { bestD2 = d2; bestM = m; }
        }
        if (bestM) {
          if (Number.isFinite(bestM.jump_win_idx) && bestM.jump_win_idx >= 0 &&
              typeof setCur === 'function') {
            setCur(bestM.jump_win_idx);
          }
          if (typeof _showFishInspectPopover === 'function') {
            _showFishInspectPopover(e, cv, bestM.si);
          }
          e.preventDefault();
          e.stopPropagation();
          return;
        }
      }
      // v3.87: shift+click on PC1 → inspect nearest tracked fish for the
      // active candidate. Shows a popover with the fish's votes, subband_path,
      // regime, confidence, and subband_stability. No-op when no candidate
      // is active or the click misses all tracked-fish lines. Without shift,
      // falls through to the existing click-to-jump behavior.
      if (e.shiftKey && cv.dataset.linesSource === 'pc1' &&
          typeof _showFishInspectPopover === 'function') {
        const handled = _showFishInspectPopover(e, cv);
        if (handled) {
          e.preventDefault();
          e.stopPropagation();
          return;
        }
      }
      const rect = cv.getBoundingClientRect();
      const pad = { l: 44, r: 16 };
      const x = e.clientX - rect.left;
      const plotW = rect.width - pad.l - pad.r;
      const frac = Math.max(0, Math.min(1, (x - pad.l) / plotW));
      // v3.99 turn 3: match the Z-panel's Mb-based mapping. Previously this
      // used window-index-fraction (`frac * (wins.length - 1)`), which is
      // wrong when windows are non-uniform in Mb across the chromosome —
      // the click would jump to a window at the wrong genomic position
      // relative to where the user clicked, and clicks at the same x
      // position across the Z-panel and lines-panel would land on
      // different windows. Now both panels use the same Mb axis: convert
      // the click fraction to an Mb position via currentMbRange(), then
      // find the window with the closest center_mb.
      const d = state.data;
      const _mbR = (typeof currentMbRange === 'function')
        ? currentMbRange()
        : { mbMin: d.windows[0].center_mb, mbMax: d.windows[d.n_windows - 1].center_mb };
      // v4 turn 114c (remaining): cs-breakpoint click-to-jump. Same red
      // dashed lines that drawLinesPanel paints across each subpanel get
      // a click-target priority over the generic Mb-frac → setCur fallback,
      // so a click near a cs-bp line lands on that breakpoint's window
      // exactly. Re-derive the per-subpanel toX from rect + pad here
      // (drawLinesPanel uses the same pad constants).
      if (typeof _ensureCsOverlayIndex === 'function') {
        const csIdx = _ensureCsOverlayIndex();
        if (csIdx && csIdx.bps.length > 0) {
          const _mbMinL = _mbR.mbMin, _mbMaxL = _mbR.mbMax;
          const _toX = (bp) => {
            if (bp.mb < _mbMinL || bp.mb > _mbMaxL) return NaN;
            return pad.l + ((bp.mb - _mbMinL) / (_mbMaxL - _mbMinL)) * plotW;
          };
          const hit = _csBpHitTestXFromList(csIdx.bps, x, _CS_BP_HIT_TOL_PX, _toX);
          if (hit) {
            _csBpJumpToWindow(hit);
            return;
          }
        }
      }
      const targetMb = _mbR.mbMin + frac * (_mbR.mbMax - _mbR.mbMin);
      let bestWin = 0, bestD = Infinity;
      for (let i = 0; i < d.n_windows; i++) {
        const dd = Math.abs(d.windows[i].center_mb - targetMb);
        if (dd < bestD) { bestD = dd; bestWin = i; }
      }
      if (typeof setCur === 'function') setCur(bestWin);
    });
    // v4 turn 114d: cs-breakpoint hover-glow on this subpanel. Same toX
    // mapping as the click hit-test above; tolerance is _CS_BP_HOVER_TOL_PX
    // (one px more generous than _CS_BP_HIT_TOL_PX so the glow engages
    // slightly before a click would commit). Idempotent via the wiring
    // helper's internal flag — buildLinesPanel rebuilds subpanels on every
    // call, but each new canvas is a fresh DOM node so the flag is fresh
    // too.
    if (typeof _wireCsBpHoverOnCanvas === 'function') {
      _wireCsBpHoverOnCanvas(cv, (canvas, evt, csIdx) => {
        if (!state.data) return null;
        const rect = canvas.getBoundingClientRect();
        const x = evt.clientX - rect.left;
        const pad = { l: 44, r: 16 };
        const plotW = rect.width - pad.l - pad.r;
        if (plotW <= 0) return null;
        const _mbR = (typeof currentMbRange === 'function')
          ? currentMbRange()
          : { mbMin: state.data.windows[0].center_mb,
              mbMax: state.data.windows[state.data.n_windows - 1].center_mb };
        const _toX = (bp) => {
          if (bp.mb < _mbR.mbMin || bp.mb > _mbR.mbMax) return NaN;
          return pad.l + ((bp.mb - _mbR.mbMin) / (_mbR.mbMax - _mbR.mbMin)) * plotW;
        };
        return _csBpHitTestXFromList(csIdx.bps, x, _CS_BP_HOVER_TOL_PX, _toX);
      });
    }
    // v3.76: lasso pointer handlers — only on the PC1 sub-canvas. Drag a
    // rectangle to select samples whose trajectory passes through it. On
    // pointer-up with non-trivial motion, compute the lassoed sample set
    // and store it; the Confirm button then commits to state.tracked.
    if (src === 'pc1') {
      let dragging = false;
      let downX = 0, downY = 0;
      cv.addEventListener('pointerdown', e => {
        if (!state.linesLassoActive) return;
        // We own the gesture. Capture pointer; track screen→canvas mapping.
        dragging = true;
        const rect = cv.getBoundingClientRect();
        downX = e.clientX - rect.left;
        downY = e.clientY - rect.top;
        state.linesLassoRect = { x0: downX, y0: downY, x1: downX, y1: downY };
        state.linesLassoCommitted = null;
        cv.setPointerCapture(e.pointerId);
        e.preventDefault();
      });
      cv.addEventListener('pointermove', e => {
        if (!dragging) return;
        const rect = cv.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        if (state.linesLassoRect) {
          state.linesLassoRect.x1 = x;
          state.linesLassoRect.y1 = y;
        }
        // Live redraw (cheap — single subcanvas)
        if (typeof drawLinesPanel === 'function') drawLinesPanel();
      });
      cv.addEventListener('pointerup', e => {
        if (!dragging) return;
        dragging = false;
        try { cv.releasePointerCapture(e.pointerId); } catch (_) {}
        if (!state.linesLassoRect) return;
        const r = state.linesLassoRect;
        const dx = Math.abs(r.x1 - r.x0), dy = Math.abs(r.y1 - r.y0);
        // Negligible motion → treat as click (let click handler run normally,
        // jumps to that window). Clear rect so no overlay sticks around.
        if (dx < 4 && dy < 4) {
          state.linesLassoRect = null;
          if (typeof drawLinesPanel === 'function') drawLinesPanel();
          return;
        }
        // Real drag → commit rect, compute lasso. Suppress the trailing click.
        cv.__lassoSwallowClick = true;
        state.linesLassoCommitted = {
          x0: Math.min(r.x0, r.x1), y0: Math.min(r.y0, r.y1),
          x1: Math.max(r.x0, r.x1), y1: Math.max(r.y0, r.y1),
        };
        state.linesLassoRect = null;
        const lassoed = _computeLinesLassoSamples(cv, state.linesLassoCommitted);
        state.linesLassoSelected = lassoed;
        if (typeof _updateLinesLassoUI === 'function') _updateLinesLassoUI();
        if (typeof drawLinesPanel === 'function') drawLinesPanel();
      });
      cv.addEventListener('pointercancel', () => {
        dragging = false;
        state.linesLassoRect = null;
        if (typeof drawLinesPanel === 'function') drawLinesPanel();
      });
      // v3.91: cursor affordance — show a pointer cursor when hovering over
      // a ⚠ jumper marker so the user knows it's clickable. Cheap: just walks
      // the marker hit list (≤ |state.tracked|) on each move event and flips
      // cv.style.cursor. Doesn't fire when dragging (lasso owns the gesture).
      cv.addEventListener('mousemove', e => {
        if (dragging) return;
        const hits = cv.__candJumperMarkers;
        if (!Array.isArray(hits) || hits.length === 0) {
          if (cv.style.cursor === 'pointer') cv.style.cursor = '';
          return;
        }
        const r2 = cv.getBoundingClientRect();
        const cx = e.clientX - r2.left;
        const cy = e.clientY - r2.top;
        const HIT_R = 7;
        let hovering = false;
        for (const m of hits) {
          const dx = m.x - cx, dy = m.y - cy;
          if (dx * dx + dy * dy <= HIT_R * HIT_R) { hovering = true; break; }
        }
        const want = hovering ? 'pointer' : '';
        if (cv.style.cursor !== want) cv.style.cursor = want;
      });
    }
    container.appendChild(sub);
  }
}

// --- buildTrackPanels() — legacy lines 32805-32848 ---
export function buildTrackPanels(state) {
  const container = document.getElementById('tracksContainer');
  if (!container) return;
  container.innerHTML = '';
  if (!state.data || !state.data.tracks) return;
  const labels = Object.keys(state.data.tracks);
  if (labels.length === 0) return;
  for (const label of labels) {
    const trk = state.data.tracks[label];
    const panel = document.createElement('div');
    panel.className = 'track-panel';
    panel.dataset.trackLabel = label;
    const cv = document.createElement('canvas');
    panel.appendChild(cv);
    const lbl = document.createElement('div');
    lbl.className = 'track-label';
    lbl.innerHTML = `<b>${label}</b>`;
    panel.appendChild(lbl);
    const rng = document.createElement('div');
    rng.className = 'track-range';
    if (isFinite(trk.min) && isFinite(trk.max)) {
      rng.textContent = `${formatTrackVal(trk.min)} – ${formatTrackVal(trk.max)}`;
    }
    panel.appendChild(rng);
    // Click-to-jump (same gesture as Z panel)
    panel.addEventListener('click', e => {
      const rect = cv.getBoundingClientRect();
      const pad = { l: 44, r: 16 };
      const x = e.clientX - rect.left;
      const plotW = rect.width - pad.l - pad.r;
      const frac = Math.max(0, Math.min(1, (x - pad.l) / plotW));
      const wins = state.data.windows;
      const mbMin = wins[0].center_mb, mbMax = wins[wins.length - 1].center_mb;
      const targetMb = mbMin + frac * (mbMax - mbMin);
      let bestI = 0, bestD = Infinity;
      for (let i = 0; i < wins.length; i++) {
        const dd = Math.abs(wins[i].center_mb - targetMb);
        if (dd < bestD) { bestD = dd; bestI = i; }
      }
      setCur(bestI);
    });
    container.appendChild(panel);
  }
}

// --- drawTracks() — legacy lines 32860-32872 ---
export function drawTracks(state) {
  if (!state.data || !state.data.tracks) return;
  const container = document.getElementById('tracksContainer');
  if (!container) return;
  const panels = container.querySelectorAll('.track-panel');
  for (const panel of panels) {
    const label = panel.dataset.trackLabel;
    const cv = panel.querySelector('canvas');
    const trk = state.data.tracks[label];
    if (!trk || !trk.values) continue;
    drawOneTrack(cv, trk, label);
  }
}

// --- refreshLinesColorMode() — legacy lines 33165-33193 ---
export function refreshLinesColorMode(state) {
  const sel = (typeof document !== 'undefined' && document.getElementById)
    ? document.getElementById('linesColorModeSelect')
    : null;
  // Validate state slot — fall back to kmeans if invalid or its layer gone
  const valid = _LINES_COLOR_MODES.some(m => m.id === state.linesColorMode);
  if (!valid || !_isLinesColorModeAvailable(state.linesColorMode)) {
    state.linesColorMode = 'kmeans';
  }
  // Sync the <select>: enable/disable each option based on layer availability
  if (sel && sel.options) {
    for (let i = 0; i < sel.options.length; i++) {
      const opt = sel.options[i];
      const def = _LINES_COLOR_MODES.find(m => m.id === opt.value);
      if (!def) continue;
      const avail = _isLinesColorModeAvailable(def.id);
      opt.disabled = !avail;
      // Refresh title (tooltip) to reflect current availability
      if (avail) {
        opt.title = `Color by ${def.label} — source layer present.`;
      } else if (def.layer) {
        opt.title = `Needs ${def.layer} JSON layer (currently not loaded).`;
      } else {
        opt.title = `Color by ${def.label}.`;
      }
    }
    sel.value = state.linesColorMode;
  }
}

// --- setLinesPanelCandidateBands() — legacy lines 34002-34007 ---
export function setLinesPanelCandidateBands(state, b) {
  const _state = (typeof window !== 'undefined' && window.state) ? window.state : state;
  _state.linesPanelCandidateBands = !!b;
  try { localStorage.setItem(_LINES_PANEL_CAND_BANDS_KEY, b ? '1' : '0'); } catch (_) {}
  if (typeof drawLinesPanel === 'function') drawLinesPanel();
}
