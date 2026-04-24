# Session audit — 2026-04-23 (Phase 4 framework v5 → v6 → v7 + KDE fix)

Append to the cumulative SESSION_AUDIT. This chat ran across multiple
back-and-forth turns that produced three framework iterations
(v5, v6, v7) plus one targeted statistical fix (HSM patch for
fragment-mode estimation). This doc records what was produced, what
was superseded, and what is shipping.

---

## Timeline

### Part 1 — v5 new-evidence layer (continuation from prior session)

Framework v5 added three existing-but-unwired capabilities:
- `cheat29b_assembled_junction.py` — parses CONSENSUS / HOMLEN / HOMSEQ
  from PRECISE DELLY and Manta SV records. Previously these fields were
  written by the SV callers but not consumed downstream.
- `bnd_sided_support.py` — single-sided BND rescue. DELLY BND records
  at one side of an inversion often lack a mate on the other side due
  to repetitive boundary regions. v5 uses the single-sided evidence
  rather than discarding.
- `cross_species_bridge.py` — v5 version. Reads between-species SV
  events (STEP_02) and flank coherence (STEP_09c), joins to the per-
  candidate schema. Writes `q5_conservation_class`.

All three have schemas in `schemas/`. All three shipped.

### Part 2 — v6 Kuang-inspired cross-species scheme (honest about Dollo)

User's direct feedback: "I would not dare to use dollo but if you
think its ok we do both. I thought dollo would find many incongruous
results." This was correct. Dollo parsimony is brittle for inversions
(forbids re-gain after loss, which is exactly the assumption real
recurrence violates).

v6 restructured Q5:
- PRIMARY: direct between-species SV-event overlap (Kuang 2026 style)
- SECONDARY: flank gene-family coherence (existing STEP_09c)
- TERTIARY: tree-parsimony polarization with confidence
- CROSS-CHECK ONLY: Dollo, reported as agree/disagree with tree

Ship: `cross_species_bridge_v6.py` + `synteny_v6.schema.json`.
v5's `cross_species_bridge.py` + `synteny_dollo.schema.json` are now
deprecated (still functional, just not preferred).

Six conservation classes now: `species_specific_polymorphism`,
`species_specific_rearranged`, `shared_ancestral`, `ASEAN_diagnostic`,
`unclassified_complex`, `untested`.

### Part 3 — v6 interior-structure diagnostic (WRONG APPROACH, superseded)

I proposed `interior_structure_diagnostic.py` — a Python script that
computed per-group row-SDs, PELT changepoints, and a Hartigan dip proxy
on the sample × position dosage matrix, to classify interior complexity.

User pushed back: "its been months I have [this] in my head" — and
uploaded the `breakpoint_pipeline` tarball showing that the per-carrier
ancestral fragment scan + KDE mode + bootstrap CI pipeline was already
designed and mostly built, with explicit rationale in METHODOLOGY.md.

My v6 interior diagnostic was a coarser re-derivation of the same
signal.

Also proposed DUP-TRP/INV-DUP framing for the LG28 interior banding.
User pushed back: "that looks like cancer research to me but not fit
for catfish inversions." Correct. Human cancer aCGH template doesn't
apply to 2.89 Mb polymorphic inversion on BEAGLE 0-2 dosage.

**Superseded by v7.** `interior_structure_diagnostic.py` and
`interior_structure.schema.json` are NOT in the final bundle.

### Part 4 — v7 correction (DUP-TRP gate + bridge to existing pipeline)

v7 fixed the v6 mistakes:

1. **Gate DUP-TRP/INV-DUP as out-of-scope, don't reconstruct.**
   New terminal label `complex_rearrangement_out_of_scope`. Fires when
   interior_recomb_fraction > 70%. Candidates receiving it are retained
   in the catalog but excluded from mechanism/age/conservation. One
   honest manuscript sentence.

2. **Bridge to existing breakpoint_pipeline.**
   `bp_pipeline_bridge.py` reads the three existing output TSVs
   (candidate_breakpoints_consensus.tsv, candidate_ancestral_fragments.tsv,
   candidate_breakpoints_per_method.tsv) and writes derived Phase 4
   keys. No new biology, no new statistics. Pure wiring.

3. **v7 assigner with graduated interior suffixes.**
   `clean_simple` / `edge_recombinants_only` → no suffix
   `bimodal_boundary_signal` → `_with_partial_interior_recombination`
   `deep_interior_recombinants` → `_with_deep_recombinants`
   `complex_rearrangement_out_of_scope` → terminal label

Label count: 17 total (16 regular + 1 terminal gate).

Ship: `bp_pipeline_bridge.py`, `fragment_distribution.schema.json`,
`assign_structural_class_v7.py`.

Smoke-tested:
- Clean LG28-like data → `supported_balanced_inversion_NAHR_like_hypothesis`
  with appropriate suffix
- >70% recombinant → `complex_rearrangement_out_of_scope`

### Part 5 — KDE fix (HSM patch for fragment-mode estimation)

User: "I read somewhere that also the method of kde is not fit so much
it over smooths. Can you double check."

Verified via web search. User's memory is correct. Silverman's rule
(`bw.nrd0` in R) assumes Gaussian density; it's known to oversmooth
multimodal and skewed distributions. The per-carrier fragment-boundary
distribution is specifically a sharp peak + recombinant tail — exactly
the shape Silverman handles worst.

References located:
- Nielsen, CSwR §2: explicit statement that Silverman oversmooths
  multimodal densities
- statsmodels PR #1271: worked examples of Silverman oversmoothing
  on bimodal/skewed data
- Sheather-Jones 1991, Bickel 2002a, Bickel & Frühwirth 2006:
  alternative mode estimators for asymmetric data
- Wikipedia KDE article: bandwidth selection for heavy-tailed
  distributions is documented difficulty

Recommended fix: half-sample mode (HSM; Bickel & Frühwirth 2006,
*Comp Stat & Data Anal* 50:3500-3530).
- Bandwidth-free
- Breakdown point 0.5
- Low-bias on asymmetric distributions
- On CRAN as `modeest::hsm()`
- Fast O(n log n), comparable to KDE

Ship: `02_ancestral_fragments_hsm_patch.R` — drop-in replacement for
`summarize_fragment_side()` in user's existing
`breakpoint_pipeline/02_ancestral_fragments.R`.

Minor issues also flagged in review doc (not fatal, just noted):
1. Current code double-truncates bandwidth (trims to 5-95% before
   calling bw.nrd0 which internally robustifies via IQR — redundant).
2. Current bootstrap uses fixed bandwidth across resamples → CI
   narrower than a full bootstrap would give. HSM has no bandwidth so
   this goes away.

---

## Final bundle inventory

Bundle: `phase4_chat_final/` (6 scripts + 5 schemas + 3 docs + 1 launcher)

### Scripts (6 files)
| File | Source | Role |
|---|---|---|
| `cheat29b_assembled_junction.py` | v5 | Parse PRECISE SV junction forensics |
| `bnd_sided_support.py` | v5 | Single-sided BND rescue |
| `cross_species_bridge_v6.py` | v6 | Kuang-style cross-species annotation |
| `bp_pipeline_bridge.py` | v7 | Bridge to user's breakpoint_pipeline |
| `assign_structural_class_v7.py` | v7 | Final label with interior gating |
| `02_ancestral_fragments_hsm_patch.R` | kde_fix | HSM mode estimator patch |

### Schemas (5 files)
| File | Block type | Primary producer |
|---|---|---|
| `mechanism_assembled.schema.json` | `mechanism_assembled` | cheat29b |
| `bnd_sided_support.schema.json` | `bnd_sided_support` | bnd_sided_support |
| `flank_coherence.schema.json` | `flank_coherence` | cross_species_bridge |
| `synteny_v6.schema.json` | `synteny_v6` | cross_species_bridge_v6 |
| `fragment_distribution.schema.json` | `fragment_distribution` | bp_pipeline_bridge |

### Launchers (1 file)
| File | Role |
|---|---|
| `run_phase4_v5_new_evidence.sh` | Orchestrator for v5 new-evidence layer |

### Docs (3 files)
| File | Purpose |
|---|---|
| `phase4_framework_v7_correction.md` | v6→v7 rationale; DUP-TRP gating decision |
| `toolkit_audit.md` | Merge / combine / simplify opportunities |
| `kde_mode_estimation_review.md` | HSM recommendation + literature basis |

### Explicitly NOT shipped (superseded or mistake)
- `interior_structure_diagnostic.py` (v6) — superseded by user's
  `02_ancestral_fragments.R` + `bp_pipeline_bridge.py`
- `interior_structure.schema.json` (v6) — replaced by
  `fragment_distribution.schema.json`
- `assign_structural_class_v5.py` (v5) — replaced by v7
- `assign_structural_class_v6.py` (v6) — replaced by v7
- `cross_species_bridge.py` (v5) — replaced by v6
- Any DUP-TRP / INV-DUP reconstruction code — gated out as
  out-of-scope per user decision

---

## Key decisions recorded this session

### D21 — Dollo parsimony is a cross-check, never primary

User pushed back on Dollo. Agreed. Primary cross-species annotation
uses direct between-species SV-event overlap (Kuang 2026 Molecular
Ecology scheme) + flank gene-family coherence + tree polarization.
Dollo reported only as agreement/disagreement with tree.

### D22 — DUP-TRP/INV-DUP architecture is out-of-scope

User pushed back on my DUP-TRP framing. Agreed. Cancer-genomics
template does not fit catfish hatchery polymorphic inversions.
Honest handling is gating, not reconstruction. Terminal label
`complex_rearrangement_out_of_scope` added. One-sentence
manuscript framing: excluded from mechanism/age/conservation,
retained in catalog for transparency.

### D23 — User's breakpoint_pipeline IS the interior signal

Do not write new interior-complexity diagnostics. The existing
per-carrier ancestral fragment scan + KDE mode + bootstrap CI
pipeline already produces the signal. Phase 4 consumes it via
the bridge script `bp_pipeline_bridge.py`.

### D24 — HSM replaces Silverman-bandwidth KDE for mode estimation

Silverman's rule oversmooths the fragment-boundary distribution
shape. Half-sample mode (Bickel & Frühwirth 2006) is bandwidth-free,
robust, published, CRAN-available. Switch primary estimator; keep
KDE as secondary for sensitivity.

### D25 — Wiring comes before new modules (D19 reaffirmed)

Session audit §D19 from 2026-04-20: next session must be wiring, not
module-building. I walked into the trap once in Part 3 (v6 interior
diagnostic). v7 corrected this. Bundle is now a clean wiring package.
No more new modules in next session.

### D26 — LG28 is test case, not manuscript flagship

User explicit: "LG28 its not main result its just the first one I
focused for testing. But that inversion on chromosome LG28 is small...
Its just that I choose this one for testing but in the manuscript I
probably put a better and bigger one for the example. And that LG28
its not that simple inside is not clean it looks like complicated."

Honest manuscript strategy: LG28 as methods validation case
(pipeline correctly labels it as `_with_interior_complexity` rather
than forcing it into `_simple` = credibility win). Pick a cleaner
candidate for flagship.

---

## Unresolved / deferred

1. **Kuang-style boundary annotations (F4 from v6 framework)** —
   LOESS repeat density ±100 kb, dyad symmetry at breakpoints,
   centromeric/satellite hits, StainedGlass-style internal alignment.
   ~300 lines R, not written. Deferred past the wiring session.

2. **Internal alignment BLAT wrapper (F5)** — ~150 lines. Deferred.

3. **Toolkit merges (D.1 of toolkit_audit.md)** — cheat14+cheat27,
   cheat29+cheat29b dispatcher, STEP_C01l+m, existence_layer_b merge.
   Not urgent. Do after public release planning.

4. **Toolkit wiring gaps (D.2)** — 4 small wiring gaps flagged:
   breakpoint_evidence_audit Q7B, cheat30 GDS in Q5, v7
   `q_overall_structural_class` as Axis 5, gene_conversion_detector
   output in Q2. ~50 lines total. Ship with main wiring session.

5. **v9→v10 alias cleanup (D.3 simplify 1)** — 15 alias keys still in
   build_key_spec for backwards compat. v9 writers gone. Not urgent.

6. **phase_12_cargo consumer wiring** — BREEDING_A/C/D scripts need
   `final_label.json` as input. Mechanical connection. Ship with
   main wiring session.

7. **Real diptest vs proxy** — `bp_pipeline_bridge.py` uses a
   Python stdlib proxy for Hartigan dip. For manuscript, switch to
   R `diptest::dip.test()` on the per-side fragment distributions.
   One-liner to add in R; trivial.

---

## Runtime validation status

Everything is smoke-tested on synthetic data. Runtime validation on
real LANTA outputs waits for:
- HPC back up (user said Monday 2026-04-27)
- `breakpoint_pipeline/run_pipeline.sh` run on LG28 with HSM patch
  installed
- Full Phase 4 pipeline run on the 226-sample cohort with v7 assigner

No R runtime tests were run this session (I don't have R in sandbox,
and real LG28 data is not in sandbox). Python scripts are ast.parse
clean and integration-tested on synthetic data.

---

## File-naming conventions used

- New scripts keep `_v5`, `_v6`, `_v7` suffix for audit trail
- Patch files named `<original>_hsm_patch.R` — indicates what they
  replace
- Schemas use singular snake_case block-type names matching the
  `block_type` field inside the JSON

---

## End of audit for this chat
