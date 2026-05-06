# P5.1 — Boundaries + SV evidence page skeleton

## Risk: medium
## Lines changed: ~200 (mostly empty-state HTML)
## Depends on: P3.2 (the page must already be merged into page11)
## Verification: page11 has 10 sections, each with empty-state copy if data missing

---

## What you said

The Boundaries + SV evidence page should include:

1. candidate coordinate track
2. left/right statistical boundary zones
3. internal transition zones if present
4. local-PCA / dosage / GHSL / theta edge tracks
5. dosage heatmap, samples ordered by H_class and family
6. SV evidence table and markers
7. MAPQ0 / ambiguity track placeholder
8. future split/clipped/discordant-mate evidence clusters
9. group association table
10. final classification

Plus the discipline:
- Don't claim exact molecular breakpoint from dosage/theta/PCA.
- BND ≠ raw read; INV is VCF-level.
- Cluster SVs by zone and event/mate relationships.
- Tier 1: BND/INV/split-read/assembly. Tier 2: DEL/DUP/INS near
  boundary with group association. Tier 3: SV inside body
  associated with groups. Tier 4: family/repeat artifact.

## Skeleton

This is a scaffold. Each section starts as empty-state with a
description and the data dependency. As bridges land (P4.1, P4.2,
SV ingest), each section becomes alive.

```html
<!-- turn 129 P5.1: Boundaries + SV evidence page skeleton.
     Ten sections corresponding to the boundary-evidence stack
     defined in the design discussion. Each section starts in
     empty-state and fills in as bridges/layers ship. -->
<div class="page" id="page11">

  <header class="page-header">
    <h2>Boundaries + SV evidence</h2>
    <p class="page-subtitle">
      Statistical transition zones from cohort signal, plus structural
      variant evidence. <b>Not exact molecular breakpoints</b> —
      those need split reads, long reads, or assembly junctions.
    </p>
  </header>

  <!-- Section 1: candidate coordinate track -->
  <section class="boundary-section" id="bs-coord">
    <h3>1 · Candidate coordinate track</h3>
    <div class="section-body" id="bs-coord-body">
      <!-- shows: chrom, [start_bp, end_bp], width, candidate id, K, family -->
      <div class="empty-hint">No candidate active. Promote one from page 1
        (local PCA |z|) or page 3 (candidate focus).</div>
    </div>
  </section>

  <!-- Section 2: left/right statistical boundary zones -->
  <section class="boundary-section" id="bs-zones">
    <h3>2 · Statistical boundary zones</h3>
    <div class="section-body" id="bs-zones-body">
      <!-- existing E/F (auto-propose) zones from page11 v0 -->
      <div class="empty-hint">Activate a candidate to compute boundary zones.</div>
    </div>
  </section>

  <!-- Section 3: internal transition zones -->
  <section class="boundary-section" id="bs-internal">
    <h3>3 · Internal transition zones</h3>
    <div class="section-body" id="bs-internal-body">
      <!-- if any per-window dosage/theta jumps inside the candidate body, list -->
      <div class="empty-hint">No internal transitions detected for this
        candidate (or candidate not yet active).</div>
    </div>
  </section>

  <!-- Section 4: edge-evidence tracks (PCA / dosage / GHSL / theta) -->
  <section class="boundary-section" id="bs-edge-tracks">
    <h3>4 · Edge-evidence tracks</h3>
    <div class="section-body" id="bs-edge-tracks-body">
      <!-- four mini-tracks aligned to candidate region:
           - local PCA |z| edge slope
           - dosage edge slope (needs P4.1 bridge)
           - GHSL edge slope
           - θπ edge slope (needs P4.2 bridge) -->
      <div class="empty-hint">Activate a candidate to render edge tracks.</div>
    </div>
  </section>

  <!-- Section 5: dosage heatmap (rows = samples, cols = SNPs) -->
  <section class="boundary-section" id="bs-dosage">
    <h3>5 · Dosage heatmap</h3>
    <p class="section-hint">
      Rows: samples ordered by H_class then family. Columns: SNPs across
      candidate region. Colour: dosage 0/1/2/missing.
      <b>Not</b> a one-band collapse — full sample × marker matrix.
    </p>
    <div class="section-body" id="bs-dosage-body">
      <!-- existing renderDosageHeatmap mounts here. Once P4.1 lands,
           it auto-fills from /api/dosage/chunk. -->
      <div class="empty-hint">Dosage panel not wired yet. See patch P4.1.</div>
    </div>
  </section>

  <!-- Section 6: SV evidence table -->
  <section class="boundary-section" id="bs-sv-table">
    <h3>6 · SV evidence</h3>
    <p class="section-hint">
      <b>Interpretation rules</b> (see specs/S5):
      INV is a VCF-level call, not raw read data.
      BND is a breakend variant call.
      Raw read evidence (split reads, discordant mates) lives in BAMs.
      DEL/DUP/INS near boundary = secondary context, not breakpoint proof.
      SV inside body = linked marker, not breakpoint proof.
    </p>
    <div class="section-body" id="bs-sv-table-body">
      <!-- table from sv_variant_catalog.tsv + sv_sample_genotypes.tsv
           grouped by zone (left_boundary / right_boundary / body /
           left_flank / right_flank). Columns: sv_id, type, caller,
           zone, distance_to_edge, tier (1-4). -->
      <div class="empty-hint">SV catalog not loaded. Drop
        <code>sv_variant_catalog.tsv</code> +
        <code>sv_sample_genotypes.tsv</code> via the load-enrichment
        flow (spec S1).</div>
    </div>
  </section>

  <!-- Section 7: MAPQ0 / ambiguity track placeholder -->
  <section class="boundary-section" id="bs-mapq0">
    <h3>7 · MAPQ0 / mappability</h3>
    <p class="section-hint">
      MAPQ0 density and mappability tracks. Useful for distinguishing
      "real" boundary-zone SV signal from repeat/ambiguity artifacts.
      <b>Not direct breakpoint proof</b> by itself.
    </p>
    <div class="section-body" id="bs-mapq0-body">
      <div class="empty-hint">MAPQ0 layer not yet loaded. Future
        bridge will read from per-chrom mappability + MAPQ0 windows.</div>
    </div>
  </section>

  <!-- Section 8: future split/clipped/discordant-mate evidence -->
  <section class="boundary-section" id="bs-read-evidence">
    <h3>8 · Read-evidence clusters (planned)</h3>
    <p class="section-hint">
      Clustered split-read / soft-clipped / discordant-mate-pair
      counts per candidate boundary zone, stratified by H_class
      with coverage normalization.
      <b>Future module</b> — see specs/S3 (Bayesian breakpoint scoring).
    </p>
    <div class="section-body" id="bs-read-evidence-body">
      <div class="empty-hint">Pending bridge to BAM-derived per-zone
        evidence clusters. Spec S3.</div>
    </div>
  </section>

  <!-- Section 9: group association table -->
  <section class="boundary-section" id="bs-group-assoc">
    <h3>9 · Group association</h3>
    <p class="section-hint">
      For each SV in the boundary zone: do its genotypes correlate
      with H_class membership? Pattern label assigned (canonical
      breakpoint marker / dominant marker / het-specific / etc.) per
      spec S1's pattern vocabulary.
    </p>
    <div class="section-body" id="bs-group-assoc-body">
      <div class="empty-hint">Pending sv_group_enrichment_results.tsv
        (spec S1, table 5).</div>
    </div>
  </section>

  <!-- Section 10: final classification -->
  <section class="boundary-section" id="bs-classification">
    <h3>10 · Final classification</h3>
    <p class="section-hint">
      Rolled-up label for this candidate's boundary status:
      <code>breakpoint_resolved_candidate</code>,
      <code>breakpoint_supported_interval</code>,
      <code>boundary_support_signal</code>,
      <code>linked_internal_marker</code>,
      <code>mapq0_repeat_ambiguity</code>,
      <code>family_artifact</code>, or
      <code>unresolved_noise</code>.
    </p>
    <div class="section-body" id="bs-classification-body">
      <div class="empty-hint">Pending all upstream sections. Default
        when nothing is loaded: <code>unresolved_noise</code>.</div>
    </div>
  </section>

</div>
```

## CSS (minimal)

```css
.boundary-section {
  border-bottom: 1px solid var(--rule);
  padding: 10px 22px 14px;
}
.boundary-section:last-child { border-bottom: none; }
.boundary-section h3 {
  margin: 0 0 4px;
  font-family: var(--mono); font-size: 12px;
  font-weight: 500;
  color: var(--ink);
  letter-spacing: 0.05em;
}
.boundary-section .section-hint {
  margin: 0 0 8px;
  font-size: 11px; color: var(--ink-dim);
  line-height: 1.5;
}
.boundary-section .empty-hint {
  padding: 16px;
  background: var(--panel-2);
  border: 1px dashed var(--rule);
  border-radius: 3px;
  font-size: 11px; color: var(--ink-dim);
  font-family: var(--mono);
  text-align: center;
}
```

## Migration of existing page11

If the current `page11` already has a renderer for the
"E/F auto-propose boundary zones" UI, that becomes Section 2's body
content. Anchor: search for `boundary` or `zoneL` / `zoneR` in the
current page11 renderer.

For each currently-existing piece of page11 functionality, place it
under the right section heading. Don't rewrite the renderers; just
move the mount points.

## Verification

1. Apply P3.2 (consolidation) first.
2. Apply P5.1.
3. Reload. Click Refine → Boundaries + SV.
4. Page shows 10 numbered sections with headings.
5. Without any data loaded, every section shows its empty-hint.
6. Once a candidate is active, sections 1-2 (and possibly 4-5) come
   alive; 6-10 stay empty until SV/MAPQ0 layers load.

## Test

```js
// tests/test_p5_1_boundary_skeleton.js
const fs = require('fs');
const html = fs.readFileSync('Inversion_atlas.html', 'utf8');

let pass = 0, fail = 0;
const ok = (n, c) => { if (c) { pass++; console.log('PASS', n); } else { fail++; console.log('FAIL', n); } };

// Page exists
ok('page11 declared', /<div class="page" id="page11">/.test(html));

// All 10 sections present
const sectionIds = ['bs-coord', 'bs-zones', 'bs-internal',
  'bs-edge-tracks', 'bs-dosage', 'bs-sv-table', 'bs-mapq0',
  'bs-read-evidence', 'bs-group-assoc', 'bs-classification'];
for (const id of sectionIds) {
  ok(`section ${id} present`,
     new RegExp(`id="${id}"`).test(html));
  ok(`section body ${id}-body present`,
     new RegExp(`id="${id}-body"`).test(html));
}

// Discipline copy present
ok('mentions BND not raw read',
   /BND[^.]*VCF-level call/i.test(html));
ok('mentions INV is VCF-level',
   /INV[^.]*VCF-level/i.test(html));
ok('mentions tier system or evidence rules',
   /Tier|interpretation rules/i.test(html));

// Empty-hint pattern present
ok('empty-state copy uses .empty-hint class',
   /class="empty-hint"/.test(html));

console.log(`\n${pass}/${pass+fail}`);
process.exit(fail ? 1 : 0);
```

## What this patch does NOT do

- Does NOT implement section renderers — those land per-section
  later as bridges arrive.
- Does NOT load any SV or MAPQ0 data — the schemas live in spec S1.
- Does NOT introduce the future read-evidence module — that's spec S3.

This is purely the structural scaffold so you can iterate on
sections one at a time without re-architecting the page each time.
