# `fuzzy_merge_abandoned/` — retired region-merge engine

`STEP_C01b_2_region_merge.R` (formerly `STEP_C01b_2_merge.R`) was the
merge step that consolidated multi-tier seeded regions from C01b_1 into
merged candidate regions before scoring, using fuzzy max-min composition
of two parallel relations (membership compatibility + geometric
continuity), optionally adjusted by PHASE_01C landscape outputs.

## Why it was dropped

On test chromosomes the 1D tail-to-head fuzzy merge overmerged across
real boundaries that were visible in the 2D sim_mat:

- **LG19** — regions at 8–13 Mb and 25–30 Mb merged across the
  low-similarity gap between them
- **LG28** — regions at 7–8, 12–14, and 15–17 Mb merged into a single
  ~10 Mb block that did not reflect the sim_mat structure
- Carriers visible in the upper-triangle heatmap (raw sim_mat) did not
  match the merged-region PA — the merge connected groups that did not
  share carriers

Root cause: the membership relation was computed from soft PC1-band
k-means (k=3) which picked up family overlap between regions that
were structurally distinct. Family LD in a 226-sample hybrid hatchery
population makes most pairs of regions "compatible" on membership,
even when they are not structurally related. The geometric relation
alone was not enough to reject those pairs.

The landscape-integration patch (v8.5) reduced but did not eliminate
the problem — clear_hard boundaries from PHASE_01C penalised merges
but did not veto them, and soft / inner / diffuse boundaries remained
bridgeable.

## What replaced it

Two parallel boundary detectors running on the same `precomp/*.rds`:

- **Staircase detector** at `phase_2/2d_candidate_detection/`
  (vote-based per-window traversal → block registry +
  boundary-vote table, no seeds required)
- **PHASE_01C landscape detector** at `phase_2/2c_precomp/`
  (row-profile clustering → block registry + classified boundaries +
  blue-cross diagnosis)

Plus the seeded regions from `STEP_C01b_1_seeded_regions.R`, which are
now consumed directly at `phase_4/4a/STEP_C01d_candidate_scoring`
(`--cores_dir`) as internal evidence on candidates whose edges came
from the boundary detectors above.

This removes the 1D-merge step entirely. See
`phase_2/2c_precomp/README.md` "Architectural note — why no
`2d_seeded_merge`" for the full rationale.

## Why the file is kept

Design history — the fuzzy composition approach is documented in
several conversation threads and its retirement is not obvious from the
final pipeline alone. Keeping the script here (instead of deleting it)
lets future readers see what was tried and why it did not work.
