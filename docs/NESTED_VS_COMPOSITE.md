# NESTED_VS_COMPOSITE — what "nested composition" does and doesn't claim

Short note to settle repeated confusion about two scripts that share the
name "nested composition".

## The two scripts

**Engine — `unified_ancestry/engines/nested_composition/internal_ancestry_composition.py`**

Generic interval-internal ancestry classifier. Given any parent region
(whole chromosome, MDS candidate, C01d inversion candidate, arbitrary
interval) + per-window ancestry labels from Engine B's local-Q cache,
characterizes the interior composition:

- per-sample dominant label blocks and switches
- fragmentation score
- internal entropy
- structural type: `homogeneous | dominant_plus_secondary |
  two_block_composite | continuous_gradient | multi_block_fragmented |
  diffuse_mixed`

The engine knows nothing about inversions. It works on any interval. Renamed
from `nested_composition.py` 2026-04-24 to match the schema/block_type it
feeds (`internal_ancestry_composition`).

**Wrapper — `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py`**

Inversion-pipeline caller. Feeds C01d inversion candidates through the
engine + layers an inversion-specific interpretation on top
(`composite_flag` ∈ `{clean, maybe_composite, likely_composite,
unknown_no_engine_b}`). Writes a Tier-2 `internal_ancestry_composition`
block per candidate via the v10 registry. Keeps the "nested_composition"
filename for back-compat with `run_phase4b.sh` and orchestration.

Before the 2026-04-24 dedup pass, the wrapper imported from a local vendored
`nested_composition_core.py` that was byte-identical to the engine. Four
copies of the same algorithm existed. The wrapper now imports from the
engine via `ANCESTRY_ENGINE_DIR`-locating `sys.path` insert, similar to how
it already finds the registry API.

## What "nested" actually means

In data-structure terms, "nested" is correct:

```
parent interval (chrom / candidate / arbitrary region)
└── Engine B local-Q windows nested inside it
    ├── window 1: assigned_pop = K2
    ├── window 2: assigned_pop = K2
    ├── window 3: assigned_pop = K5
    ├── window 4: assigned_pop = K5
    └── window 5: assigned_pop = K7
```

The engine asks: "is the ancestry composition of the windows nested inside
this parent interval simple or internally structured?"

Examples on the test cases from the design discussion:

| Label sequence         | Engine verdict             |
|------------------------|----------------------------|
| `K2 K2 K2 K2 K2`       | `homogeneous`              |
| `K2 K2 K2 K5 K5 K5`    | `two_block_composite` or `continuous_gradient` |
| `K2 K5 K2 K7 K5 K2 K7` | `multi_block_fragmented` / `diffuse_mixed` |

## What "nested" does NOT mean

In inversion biology, "nested" connotes **nested inversion**
(inversion-inside-inversion). That is a STRONGER biological claim than this
analysis supports. A `likely_composite` verdict from the engine can also be
explained by any of:

- Double crossover in a single inversion
- Gene-conversion tract
- Badly-merged candidate (two real candidates mis-joined by C01d)
- Family-LD contamination bleeding into the candidate's interior
- Old haplotype background
- Sub-inversion
- Local ancestry mixture unrelated to the candidate
- Label boundary ambiguity in Engine B

Interpreting a positive composite signal as a nested inversion requires
cross-checks with SV calls, paired breakpoints, phased-HET dose, and
karyotype evidence. That composite interpretation lives at or after
`compute_candidate_status.R`, not in this engine.

## How both pieces get used

The engine is a library. It gets called twice in the full pipeline, on
different parent intervals:

**Pass 1 — chromosomes as parents.**

`unified_ancestry/run_full_pipeline.sh` Step 4 runs the engine CLI
(`internal_ancestry_composition.py`) with `--parents` set to the chromosome
or MDS-candidate TSV. Output: `nested_composition.tsv` +
`nested_composition_summary.tsv`. Tells the pipeline: which chromosomes
or large regions have fine-grained local ancestry structure that isn't
population-wide? This is the backdrop denominator.

**Pass 2 — C01d inversion candidates as parents.**

`STEP_C01i_c_nested_composition.py` (the Phase 7 wrapper) imports the
engine as a library, runs `analyze_parent_composition` per candidate, adds
a `composite_flag`, writes a registry block per candidate. Tells the
pipeline: which candidates are internally composite?

A candidate flagged composite in Pass 2 is weaker evidence if Pass 1 showed
its chromosome has broad background composite structure — the candidate
might just be overlapping that structure rather than having intrinsic
internal composition. Reading Pass 1 and Pass 2 together is the correct
interpretation; reading Pass 2 alone is where the "nested inversion"
confusion comes from.

## Naming (what is and isn't renamed)

| File                                                                                | Renamed? | Why |
|-------------------------------------------------------------------------------------|----------|-----|
| `unified_ancestry/engines/nested_composition/internal_ancestry_composition.py`      | YES, was `nested_composition.py` | Matches the schema and block_type it feeds; makes the engine's role (generic, not inversion-specific) clear in the filename. |
| `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py` | NO | Kept for back-compat with `run_phase4b.sh` and orchestration scripts. The filename is still "nested" but the docstring (and this note) say exactly what the wrapper claims. |
| `inversion_modules/phase_7_karyotype_groups/proposal/nested_composition_core.py`    | DELETED | Was a byte-identical duplicate of the engine. Wrapper now imports from the engine. |
| `inversion_modules/phase_2_discovery/_archive/analysis/nested_composition.py`       | DELETED | Archive duplicate, unused. |
| `unified_ancestry/engines/nested_composition/` (directory)                          | NO | Directory name kept for now; would require path updates in config and launchers. Renaming is low value, high breakage — deferred. |

## If you want to change the verdict after reading this

The `composite_flag` thresholds (see wrapper docstring) aren't sacred:

- `clean`: ≥80% of samples homogeneous OR dominant_plus_secondary
- `likely_composite`: >20% two_block_composite OR >20% multi_block_fragmented
- `maybe_composite`: everything else
- `unknown_no_engine_b`: Engine B Q cache was unavailable

At 9× coverage with the 226-sample hatchery cohort, the 20% threshold for
`likely_composite` may be too sensitive (family structure can easily push
20% of samples into two-block territory without a real nested event).
If the first full run produces too many `likely_composite` flags, raise
to 30% or 40% and rerun. The thresholds live at lines ~150-165 of the
wrapper in `derive_composite_flag()`.
