// =============================================================================
// atlas_js/ushape_renderer.js
// =============================================================================
// U-shape evolution module renderer for Inversion_atlas.html.
//
// Exports a single global `UShapeRenderer` namespace; everything else is
// closure-scoped. Idempotent: calling init() twice is safe.
//
// Render surfaces (each is opt-in; the renderer no-ops if the host element
// is absent — so this file is safe to load before page-specific HTML is
// ready):
//
//   1. Candidate page  — card "Evolutionary shape"  → #ushape_card
//   2. Page 6 popstats — track group "Arrangement divergence profile"
//                       → #ushape_popstats_panel
//   3. Page 4 axis chip  → #ushape_axis_chip
//   4. Page 17 stats     → #ushape_stats_table
//   5. Gallery filters   → #ushape_filter_select
//   6. Tooltip           → automatic via .ushape-pill[data-tip]
//
// Data flow:
//
//   - Server-primary mode (default):
//       UShapeRenderer.fetchForCandidate(candidateId, groupsObj, intervalObj)
//       → POSTs to popstats server /api/ushape/candidate
//       → caches the response in a candidateId → block map
//       → calls renderAll() with the block
//
//   - Static-JSON-fallback mode:
//       UShapeRenderer.loadStaticLayer(url)  // ushape_evolution_v1.json
//       → splits into candidate blocks; renders the currently active one
//
// If a candidate has no data, every surface shows a graceful "U-shape
// evolution layer unavailable for this candidate" message — the atlas
// continues working.
// =============================================================================

(function () {
  "use strict";

  const STATE = {
    server_url: null,           // set via init()
    cache: new Map(),           // candidate_id -> candidate block
    static_layer: null,         // parsed ushape_evolution_v1.json (fallback)
    last_candidate_id: null,
    fetch_in_flight: new Map(), // candidate_id -> Promise (de-dup)
  };

  const PILL_LABELS = {
    neutral_like_U_shape:                "Neutral-like U-shape",
    locally_adapted_internal_peak_like:  "Internal peak / locally-adapted-like",
    locally_adapted_breakpoint_like:     "Breakpoint-adapted-like",
    flat_high_deep_structural_haplotype: "Flat deep haplotype",
    young_weak_divergence:               "Young / weak divergence",
    asymmetric_edge:                     "Asymmetric edge",
    complex_mixed:                       "Complex mixed",
    insufficient_data:                   "Insufficient data",
  };

  const PILL_COLORS = {
    neutral_like_U_shape:                "#5e81ac",
    locally_adapted_internal_peak_like:  "#bf616a",
    locally_adapted_breakpoint_like:     "#d08770",
    flat_high_deep_structural_haplotype: "#4c566a",
    young_weak_divergence:               "#a3be8c",
    asymmetric_edge:                     "#ebcb8b",
    complex_mixed:                       "#b48ead",
    insufficient_data:                   "#888888",
  };

  // ---- helpers --------------------------------------------------------------
  const fmt = (x, n = 2) =>
    (x === null || x === undefined || !isFinite(x)) ? "—" : Number(x).toFixed(n);

  const $ = (sel, root) => (root || document).querySelector(sel);
  const setHTML = (el, html) => { if (el) el.innerHTML = html; };
  const setUnavailable = (el, why) => setHTML(el,
    `<div class="ushape-unavailable" style="color:#888;font-style:italic">
       U-shape evolution layer unavailable${why ? `: ${why}` : ""}.
     </div>`);

  // ---- API: init / load -----------------------------------------------------
  function init(opts) {
    STATE.server_url = (opts && opts.server_url) || null;
    if (opts && opts.static_layer_url) {
      loadStaticLayer(opts.static_layer_url).catch((e) =>
        console.warn("[ushape] static layer load failed:", e));
    }
    document.addEventListener("mouseover", onTooltipHover, true);
  }

  async function loadStaticLayer(url) {
    const r = await fetch(url, { cache: "no-cache" });
    if (!r.ok) throw new Error("load failed: " + r.status);
    const j = await r.json();
    if (j.format_version !== "ushape_evolution_v1") {
      console.warn("[ushape] unexpected format_version:", j.format_version);
    }
    STATE.static_layer = j;
    if (Array.isArray(j.candidates)) {
      for (const c of j.candidates) STATE.cache.set(c.candidate_id, c);
    }
  }

  // ---- API: fetch ONE candidate from the server -----------------------------
  async function fetchForCandidate(candidate_id, groups, interval) {
    if (STATE.cache.has(candidate_id)) return STATE.cache.get(candidate_id);
    if (STATE.fetch_in_flight.has(candidate_id))
      return STATE.fetch_in_flight.get(candidate_id);
    if (!STATE.server_url) {
      console.warn("[ushape] no server_url; cannot fetch", candidate_id);
      return null;
    }
    const body = {
      candidate_id,
      chrom: interval.chrom,
      start_bp: interval.start_bp,
      end_bp: interval.end_bp,
      left_breakpoint:  interval.left_breakpoint  || interval.start_bp,
      right_breakpoint: interval.right_breakpoint || interval.end_bp,
      groups,
    };
    const p = (async () => {
      try {
        const r = await fetch(STATE.server_url.replace(/\/+$/, "")
                              + "/api/ushape/candidate", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(body),
        });
        if (!r.ok) throw new Error("HTTP " + r.status);
        const j = await r.json();
        const block = j.candidate || j;
        STATE.cache.set(candidate_id, block);
        return block;
      } finally {
        STATE.fetch_in_flight.delete(candidate_id);
      }
    })();
    STATE.fetch_in_flight.set(candidate_id, p);
    return p;
  }

  // ---- API: render all surfaces for one candidate ---------------------------
  function renderAll(candidate_id) {
    STATE.last_candidate_id = candidate_id;
    const block = STATE.cache.get(candidate_id) || null;
    renderCandidateCard(block);
    renderPopstatsPanel(block);
    renderAxisChip(block);
    renderStatsTable(block);
    renderGalleryOptions();   // rebuild filter from full cache
  }

  // ---- 1. Candidate card ----------------------------------------------------
  function renderCandidateCard(block) {
    const root = $("#ushape_card");
    if (!root) return;
    if (!block) return setUnavailable(root);
    const cls = block.classification || {};
    const sc = block.shape_scores || {};
    const pillKey = cls.primary_class || "insufficient_data";
    const tip = encodeURIComponent(JSON.stringify({
      reason: cls.reason || "",
      n_homo1: block.groups_used && block.groups_used.n_homo1,
      n_homo2: block.groups_used && block.groups_used.n_homo2,
      n_het:   block.groups_used && block.groups_used.n_het,
      U: sc.dxy_u_score,
      inflank: sc.dxy_inside_flank_ratio,
      peak: sc.dxy_internal_peak_score,
      oldness: sc.oldness_score,
      flags: cls.flags || [],
    }));
    setHTML(root, `
      <div class="ushape-card" style="font-family:inherit">
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:6px">
          <span style="font-weight:600">Evolutionary shape</span>
          <span class="ushape-pill" data-tip="${tip}"
                style="padding:2px 8px;border-radius:10px;color:white;
                       background:${PILL_COLORS[pillKey] || "#888"};font-size:0.85em">
            ${PILL_LABELS[pillKey] || pillKey}
          </span>
          <span style="color:#888;font-size:0.85em">conf: ${cls.confidence || "—"}</span>
        </div>
        <table style="font-size:0.85em">
          <tr><td>U-score (dXY)</td><td>${fmt(sc.dxy_u_score)}</td>
              <td style="padding-left:16px">inside / flank dXY</td><td>${fmt(sc.dxy_inside_flank_ratio)}</td></tr>
          <tr><td>internal peak</td><td>${fmt(sc.dxy_internal_peak_score)}</td>
              <td style="padding-left:16px">oldness</td><td>${fmt(sc.oldness_score)}</td></tr>
          <tr><td>local adaptation (internal)</td><td>${fmt(sc.local_adaptation_internal_score)}</td>
              <td style="padding-left:16px">breakpoint adaptation</td><td>${fmt(sc.breakpoint_adaptation_score)}</td></tr>
        </table>
      </div>`);
  }

  // ---- 2. Page 6 popstats panel  --------------------------------------------
  // We draw a tiny inline SVG with five tracks aligned on a shared x-axis.
  function renderPopstatsPanel(block) {
    const root = $("#ushape_popstats_panel");
    if (!root) return;
    if (!block || !Array.isArray(block.windows) || block.windows.length === 0)
      return setUnavailable(root);

    const W = block.windows;
    const x0 = W[0].start_bp, x1 = W[W.length - 1].end_bp;
    const W_PX = 800, H_PX = 70, GAP = 8;
    const xpx = (x) => ((x - x0) / Math.max(1, (x1 - x0))) * W_PX;

    const tracks = [
      { key: "dxy_homo1_homo2", title: "dXY (HOMO_1 vs HOMO_2)", color: "#bf616a" },
      { key: "fst_homo1_homo2", title: "FST",                    color: "#5e81ac" },
      { key: "pi_homo1",        title: "π HOMO_1",                color: "#a3be8c" },
      { key: "pi_homo2",        title: "π HOMO_2",                color: "#88c0d0" },
      { key: "allele_freq_delta", title: "|Δp|",                  color: "#d08770" },
    ];

    const zoneFill = {
      left_flank: "#f2f2f2", left_edge: "#fde0dd",
      center: "#ffffff", right_edge: "#fde0dd", right_flank: "#f2f2f2",
    };

    // group consecutive windows by zone for shading
    const zoneRects = [];
    for (let i = 0; i < W.length; i++) {
      const z = W[i].zone, s = W[i].start_bp, e = W[i].end_bp;
      const last = zoneRects[zoneRects.length - 1];
      if (last && last.zone === z && last.end === s)
        last.end = e;
      else
        zoneRects.push({ zone: z, start: s, end: e });
    }

    const trackSvg = (t, idx) => {
      const ys = W.map(w => Number(w[t.key])).filter(v => isFinite(v));
      const vmax = ys.length ? Math.max(...ys, 1e-12) : 1;
      const path = W.map((w, i) => {
        const x = xpx((w.start_bp + w.end_bp) / 2);
        const v = Number(w[t.key]);
        const y = isFinite(v) ? (H_PX - (v / vmax) * (H_PX - 4) - 2) : null;
        return (y === null) ? null : `${i === 0 ? "M" : "L"}${x.toFixed(1)},${y.toFixed(1)}`;
      }).filter(Boolean).join(" ");
      const rects = zoneRects.map(r =>
        `<rect x="${xpx(r.start).toFixed(1)}" y="0"
               width="${(xpx(r.end) - xpx(r.start)).toFixed(1)}" height="${H_PX}"
               fill="${zoneFill[r.zone] || "#fff"}" />`).join("");
      return `
        <div style="margin-bottom:${GAP}px">
          <div style="font-size:0.8em;color:#444">${t.title}
            <span style="float:right;color:#888">max=${fmt(vmax, 4)}</span></div>
          <svg width="${W_PX}" height="${H_PX}" style="display:block">
            ${rects}
            <path d="${path}" fill="none" stroke="${t.color}" stroke-width="1.4"/>
          </svg>
        </div>`;
    };

    setHTML(root, `
      <div style="font-weight:600;margin-bottom:6px">Arrangement divergence profile</div>
      ${tracks.map(trackSvg).join("")}
      <div style="font-size:0.75em;color:#888;margin-top:4px">
        Zones: <span style="background:${zoneFill.left_flank};padding:0 4px">flank</span>
        <span style="background:${zoneFill.left_edge};padding:0 4px">edge</span>
        <span style="background:${zoneFill.center};padding:0 4px;border:1px solid #ccc">center</span>
        <span style="background:${zoneFill.right_edge};padding:0 4px">edge</span>
        <span style="background:${zoneFill.right_flank};padding:0 4px">flank</span>
      </div>`);
  }

  // ---- 3. Page 4 axis chip --------------------------------------------------
  function renderAxisChip(block) {
    const root = $("#ushape_axis_chip");
    if (!root) return;
    if (!block) return setUnavailable(root);
    const cls = block.classification || {};
    const sc = block.shape_scores || {};
    let label = PILL_LABELS[cls.primary_class] || "—";
    let detail = "";
    switch (cls.primary_class) {
      case "neutral_like_U_shape":
        detail = `U-shape: ${fmt(sc.dxy_u_score)}× edge/center`; break;
      case "locally_adapted_internal_peak_like":
        detail = `Internal peak: ${fmt(sc.dxy_internal_peak_score)}×`; break;
      case "locally_adapted_breakpoint_like":
        detail = `Edge dominance: ${fmt(sc.breakpoint_adaptation_score)}`; break;
      case "flat_high_deep_structural_haplotype":
        detail = `Flat-deep: dXY ${fmt(sc.dxy_inside_flank_ratio)}× flank`; break;
      case "young_weak_divergence":
        detail = "Young / weak"; break;
      case "asymmetric_edge":
        detail = `Asym log2=${fmt(sc.abs_log2_asymmetry)}`; break;
      default: detail = "";
    }
    setHTML(root, `
      <span class="ushape-pill"
            style="padding:2px 8px;border-radius:10px;color:white;
                   background:${PILL_COLORS[cls.primary_class] || "#888"}">
        ${label}
      </span>
      <span style="margin-left:6px;color:#444;font-size:0.85em">${detail}</span>`);
  }

  // ---- 4. Page 17 stats table ----------------------------------------------
  function renderStatsTable(block) {
    const root = $("#ushape_stats_table");
    if (!root) return;
    if (!block) return setUnavailable(root);
    const cls = block.classification || {};
    const sc = block.shape_scores || {};
    const mb = block.matched_background || {};
    const rows = [
      ["primary class",                 PILL_LABELS[cls.primary_class] || cls.primary_class || "—"],
      ["secondary class",               cls.secondary_class || "—"],
      ["confidence",                    cls.confidence || "—"],
      ["dXY inside / flank",            fmt(sc.dxy_inside_flank_ratio)],
      ["dXY U-score",                   fmt(sc.dxy_u_score)],
      ["dXY internal peak",             fmt(sc.dxy_internal_peak_score)],
      ["FST inside / flank",            fmt(sc.fst_inside_flank_ratio)],
      ["FST internal peak",             fmt(sc.fst_internal_peak_score)],
      ["π HOMO_1 / HOMO_2",             fmt(sc.pi_ratio_homo1_homo2)],
      ["oldness score",                 fmt(sc.oldness_score)],
      ["youngness score",               fmt(sc.youngness_score)],
      ["local adaptation (internal)",   fmt(sc.local_adaptation_internal_score)],
      ["neutral U-shape score",         fmt(sc.neutral_u_shape_score)],
      ["breakpoint adaptation score",   fmt(sc.breakpoint_adaptation_score)],
      ["unsupervised cluster",          cls.cluster_label || cls.unsupervised_cluster || "—"],
      ["bg z(dXY inside)",              fmt(mb.z_dxy_inside)],
      ["bg z(FST inside)",              fmt(mb.z_fst_inside)],
      ["bg z(U-score)",                 fmt(mb.z_dxy_u_score)],
      ["bg z(internal peak)",           fmt(mb.z_dxy_internal_peak)],
      ["bg mode",                       mb.bg_mode || "—"],
    ];
    setHTML(root, `
      <div style="font-weight:600;margin-bottom:6px">U-shape evolution statistics</div>
      <table style="font-size:0.85em;border-collapse:collapse">
        ${rows.map(([k, v]) => `
          <tr><td style="padding:1px 12px 1px 0;color:#555">${k}</td>
              <td><strong>${v}</strong></td></tr>`).join("")}
      </table>`);
  }

  // ---- 5. Gallery filter dropdown ------------------------------------------
  function renderGalleryOptions() {
    const sel = $("#ushape_filter_select");
    if (!sel) return;
    const counts = {};
    for (const block of STATE.cache.values()) {
      const k = (block.classification && block.classification.primary_class) || "insufficient_data";
      counts[k] = (counts[k] || 0) + 1;
    }
    const opts = Object.keys(PILL_LABELS).map(k => {
      const n = counts[k] || 0;
      return `<option value="${k}">${PILL_LABELS[k]} (${n})</option>`;
    }).join("");
    sel.innerHTML = `<option value="">— any class —</option>${opts}`;
  }

  function filterCandidatesByClass(klass) {
    const out = [];
    for (const [cid, block] of STATE.cache) {
      const k = (block.classification && block.classification.primary_class) || null;
      if (!klass || k === klass) out.push(cid);
    }
    return out;
  }

  // ---- 6. Tooltip on hover --------------------------------------------------
  let tipEl = null;
  function ensureTipEl() {
    if (tipEl) return tipEl;
    tipEl = document.createElement("div");
    tipEl.className = "ushape-tooltip";
    Object.assign(tipEl.style, {
      position: "fixed", zIndex: 9999, display: "none",
      background: "rgba(40,40,40,0.95)", color: "#fff",
      padding: "8px 10px", borderRadius: "4px",
      fontSize: "0.8em", maxWidth: "320px", pointerEvents: "none",
      lineHeight: "1.35",
    });
    document.body.appendChild(tipEl);
    document.addEventListener("mouseout", (e) => {
      if (e.target.classList && e.target.classList.contains("ushape-pill")) {
        tipEl.style.display = "none";
      }
    }, true);
    return tipEl;
  }
  function onTooltipHover(e) {
    const t = e.target;
    if (!t || !t.classList || !t.classList.contains("ushape-pill")) return;
    const raw = t.getAttribute("data-tip");
    if (!raw) return;
    let info; try { info = JSON.parse(decodeURIComponent(raw)); } catch { return; }
    const html = `
      <div><strong>Reason:</strong> ${info.reason || "—"}</div>
      <div style="margin-top:4px">
        <span>U: ${fmt(info.U)}</span> ·
        <span>inside/flank: ${fmt(info.inflank)}</span> ·
        <span>peak: ${fmt(info.peak)}</span> ·
        <span>oldness: ${fmt(info.oldness)}</span>
      </div>
      <div style="margin-top:4px">
        n: HOMO_1=${info.n_homo1 ?? "?"} · HET=${info.n_het ?? "?"} · HOMO_2=${info.n_homo2 ?? "?"}
      </div>
      ${(info.flags || []).length ?
        `<div style="margin-top:4px;color:#fbb">flags: ${info.flags.join(", ")}</div>` : ""}
      <div style="margin-top:4px;color:#aaa">
        Class: shape verdict, not a causal claim.
      </div>`;
    const tip = ensureTipEl();
    tip.innerHTML = html;
    tip.style.display = "block";
    tip.style.left = (e.clientX + 12) + "px";
    tip.style.top  = (e.clientY + 12) + "px";
  }

  // ---- export --------------------------------------------------------------
  window.UShapeRenderer = {
    init,
    loadStaticLayer,
    fetchForCandidate,
    renderAll,
    filterCandidatesByClass,
    PILL_LABELS, PILL_COLORS,
    _state: STATE, // for debugging
  };
})();
