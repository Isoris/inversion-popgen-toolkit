// shared/het_rate.js
// =====================================================================
// Heterozygosity color ramp — pure, used by every panel that paints
// per-sample het signal.
//
// Extracted from legacy Inversion_atlas.html (turn 165 close binary).
// Source line refs:
//   _HET_RAMP        line 16569
//   _hetRateColor    line 16575
//
// NOT extracted (state-coupled — extracted in inversion_discovery step):
//   _computeHetRateForL2     line 16735  (reads state.data.l2_envelopes)
//   _computeHetRateForRange  (reads state.__hetRateCache, dosage_chunks index)
//   _computeHetRateForSlab   line 16752
//   _getHetRateCache         line 16612  (touches state.__hetRateCache)
//   _invalidateHetRateCache  line 16619
//
// Those compute helpers belong to inversion_discovery because they
// depend on the dosage chunk LRU and the discovery-side state.data
// layout. They'll move during the discovery extraction step.
// =====================================================================

/**
 * Diverging warm-cold-warm ramp anchored at 0.5 (the diallelic
 * expected-HET rate). cold/warm endpoints are colorbrewer RdBu's
 * extreme tones; neutral is near-white so 0.5 reads as "expected".
 */
export const HET_RAMP = Object.freeze({
  cold:    [0x21, 0x66, 0xAC],   // 0.0 — all-homo (REF or INV)
  neutral: [0xF7, 0xF7, 0xF7],   // 0.5 — expected HET (midpoint)
  warm:    [0xB2, 0x18, 0x2B],   // 1.0 — pure-het (anomaly)
});

/**
 * Map a het rate ∈ [0, 1] to a CSS color string.
 *   - rate <= 0.5  → cold..neutral interpolation
 *   - rate >= 0.5  → neutral..warm interpolation
 *   - non-finite   → 'var(--ink-dimmer)' sentinel for "no calls"
 *
 * Out-of-range inputs are clamped to [0, 1] so they never produce
 * undefined RGB.
 *
 * @param {number} rate
 * @returns {string}  e.g. 'rgb(247,247,247)' or 'var(--ink-dimmer)'
 */
export function hetRateColor(rate) {
  if (rate == null || !Number.isFinite(rate)) {
    return 'var(--ink-dimmer)';
  }
  const r = Math.max(0, Math.min(1, rate));
  let rgb;
  if (r <= 0.5) {
    const t = r / 0.5;
    rgb = [
      Math.round(HET_RAMP.cold[0] + (HET_RAMP.neutral[0] - HET_RAMP.cold[0]) * t),
      Math.round(HET_RAMP.cold[1] + (HET_RAMP.neutral[1] - HET_RAMP.cold[1]) * t),
      Math.round(HET_RAMP.cold[2] + (HET_RAMP.neutral[2] - HET_RAMP.cold[2]) * t),
    ];
  } else {
    const t = (r - 0.5) / 0.5;
    rgb = [
      Math.round(HET_RAMP.neutral[0] + (HET_RAMP.warm[0] - HET_RAMP.neutral[0]) * t),
      Math.round(HET_RAMP.neutral[1] + (HET_RAMP.warm[1] - HET_RAMP.neutral[1]) * t),
      Math.round(HET_RAMP.neutral[2] + (HET_RAMP.warm[2] - HET_RAMP.neutral[2]) * t),
    ];
  }
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

// ---------------------------------------------------------------------
// Console-debug exposures (preserves legacy `window._hetRateColor`)
// ---------------------------------------------------------------------
if (typeof window !== 'undefined') {
  window._HET_RAMP      = HET_RAMP;
  window._hetRateColor  = hetRateColor;
}
