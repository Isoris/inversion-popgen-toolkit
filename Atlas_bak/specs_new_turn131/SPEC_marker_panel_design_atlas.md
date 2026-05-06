# SPEC — Marker panel design atlas page

**Status**: drafted turn 130 final session. Already substantially
spec'd in chat `819d8454` ("Live PCA visualization with contingency
tables") and `1b4d8e12` ("Catfish chromosome cross-species analysis").
Atlas-side scrubber already has a marker page scaffold (turn 30-ish).
This spec consolidates everything for next-chat continuation.

**Trigger** (Quentin, chat `819d8454`):
> *"Add one page/module: Marker page / PCR assay page. For each
> inversion candidate, the atlas should show: Candidate ID, Chr / start
> / end / length, Regime counts, Best diagnostic marker set, Primer
> pairs, Expected amplicon sizes, Which regime each marker identifies,
> Confidence / specificity, Warnings."*
>
> Plus from chat `1b4d8e12` (Cgar × Cmac F1 hybrid breeding context):
> *"Cross-species negative control: Bighead catfish DNA/genome —
> Confirms the marker is not just amplifying the maternal species/
> background."*

The atlas-side marker page is a **review surface** for marker panels
designed by an upstream LANTA pipeline. Atlas computes some scores
live from a `variant_afs.json` layer; LANTA emits the canonical panel
TSV.

---

## 1. Tier system (rebuilt around private indels)

Following Quentin's redesign in chat `1b4d8e12`:

| Tier | Criteria |
|---|---|
| **Tier 1** | Private indel/SNP with clean dosage (AF_REF ≤ 0.02, AF_HET ∈ [0.25, 0.75], AF_INV ≥ 0.80) AND positive/negative/het controls identified AND bighead specificity confirmed |
| **Tier 2** | Multi-marker haplotype panel (≥3 linked markers) OR strong-tag (AF_REF ≤ 0.05, AF_INV ≥ 0.70) with controls but imperfect het signal |
| **Tier 3** | Breakpoint PCR candidate (demoted from Tier 1 because breakpoint precision is uncertain at our 9× coverage) |
| **Tier 4** | Exploratory (association exists, controls incomplete) |

## 2. Atlas-side scoring

When `variant_afs.json` is loaded:

```
private_score   = AF_INV - AF_REF
dosage_score    = 1 - |AF_HET - 0.5| * 2
indel_size_pen  = (size in [10, 30] bp) ? 0 : penalty for too small/large
final_score     = private_score + dosage_score - indel_size_pen
tier            = thresholds(final_score, controls_present, bighead_status)
```

Atlas computes these live in JS so the user can adjust thresholds
interactively.

## 3. Page layout (per candidate)

```
Marker Panel — Inversion {seq_num} on {chrom}
─────────────────────────────────────────────────────────────────
Span: {start_mb}–{end_mb} Mb · Karyotype: REF/{n} HET/{n} INV/{n}

[Variant table — sortable, filterable]
  marker_id | pos | type | size | AF_REF | AF_HET | AF_INV | tier | warnings

  MRK_LG28_15123456_INDEL_5    SNP   1bp  0.00  0.49  1.00   1   ✓
  MRK_LG28_15234567_INDEL_8    DEL  18bp  0.02  0.50  0.95   1   ⚠ flank repeat
  MRK_LG28_15345678_INDEL_2    SNP   1bp  0.05  0.42  0.85   2   ✓
  ...

[Panel suggestions]
  Panel A (3 markers, all Tier 1)        → diagnostic: REF / HET / INV
  Panel B (5 markers, Tier 1+2 mix)      → robust to local mutation
  Panel C (3 markers, breakpoint-anchored) → Tier 3 (validation pending)

[Controls table]
  positive INV/INV samples:    {atlas-derived from karyotype calls}
  negative REF/REF samples:    {atlas-derived}
  het samples:                 {atlas-derived}
  cross-species negative:      {user-uploaded; bighead reference DNA}

[Primer design feasibility]
  marker_id | preferred_assay | local_var_density | warnings
  MRK_...   | KASP_or_TaqMan  | low               | -
  MRK_...   | fragment_analysis | medium          | indel size 18bp
  ...
```

## 4. Output files (the panel TSVs)

Per chat `819d8454`, three TSVs:

```
candidate_marker_catalogue.tsv  → one row per marker per candidate
candidate_marker_primers.tsv    → one row per primer pair
candidate_marker_panel_summary.tsv → one row per panel (panel = ≥3 markers)
```

Schemas detailed in chat `819d8454`. Reproduced in
`marker_panel_schema_v1.md` next to this spec.

## 5. Sample suggestion logic (controls)

Atlas auto-derives positive/negative samples from `state.candidateList`
karyotype calls:

- **Positive INV/INV**: top 3 samples by INV-arrangement support score
  (highest dosage in inverted markers, lowest in reference markers).
- **Negative REF/REF**: top 3 samples by REF-arrangement support.
- **Het**: top 3 samples with closest-to-50% het rate.

User can override via `marker_panel.json` `controls_override` block.
Cross-species controls (bighead) are user-uploaded only — atlas can't
derive them.

## 6. Implementation slices

### Slice 1 — atlas-side scoring + page (~1 turn)
- `variant_afs_v1` JSON layer schema + detection
- Live JS scoring function
- Per-candidate page section (likely on existing page-2 candidate focus,
  new sub-tab "markers")

### Slice 2 — control-sample suggestion (~0.5 turn)
- Auto-derive positive/negative/het from karyotype calls
- Override mechanism for user-uploaded controls

### Slice 3 — panel suggestion algorithm (~0.5 turn)
- Greedy panel builder: pick 3-5 highest-scoring markers spread across
  the candidate interval
- Multiple panels suggested with different tier mixes

### Slice 4 — export (~0.3 turn)
- `candidate_marker_catalogue.tsv` and friends as part of manuscript
  bundle (`SPEC_manuscript_bundle_export.md`)

## 7. Open questions

1. **Primer design itself — atlas or LANTA?** Atlas can show feasibility
   warnings (local variant density, indel size) but can't actually
   design primers. Recommend: LANTA-side `primer3` pipeline emits
   primer pairs, atlas displays them.
2. **Specificity vs sensitivity tradeoff for tier thresholds**: AF_INV
   ≥ 0.80 is strict. Some real markers will fall just below. Allow
   user-tunable thresholds.
3. **Bighead specificity check**: requires bighead reference genome
   alignment. This is upstream of the atlas; atlas just consumes the
   `bighead_status` field per marker.

## 8. Dependencies

- LANTA-side variant AF computation per candidate per K-band.
- LANTA-side primer3 wrapper (existing).
- Bighead reference genome alignment of markers (existing per the
  `catfish-synteny-toolkit`).
- Per-candidate karyotype calls (existing).

## 9. What this is NOT

- Not a primer design tool itself. The atlas displays designs computed
  upstream.
- Not a wet-lab validation framework. The atlas reports predicted
  specificity; the lab confirms.
- Not a marker quality assurance suite. Tier thresholds are signals,
  not guarantees.

## 10. Cross-references

- `SPEC_per_candidate_breeding_readiness_card.md` consumes the marker
  panel TSV for the per-candidate breeding card.
- `SPEC_manuscript_bundle_export.md` includes the panel TSVs in
  Supplementary Data.
- Chat `819d8454` for the original feature design.
- Chat `1b4d8e12` for the cross-species (Cgar × Cmac F1) tier rebuild.
- Chat `58bd3863` (Clair3 marker handoff) for the upstream marker
  classification framework.
