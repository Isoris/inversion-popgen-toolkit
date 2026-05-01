# Handoff — Diversity Atlas (next session)

**Status**: ⚪ scaffold only · ready to start building real renderers
**File**: `Atlas/Diversity_atlas.html` (33 KB, 1 page, accent green)
**Trigger to start**: HPC under maintenance for jobs but login node is up,
so all scripts can be run interactively against existing outputs that
landed before maintenance. Most of the data this atlas needs is **already
on disk** from MODULE_3.

---

## Why this atlas can ship today

The Diversity Atlas is the only sibling whose primary inputs are already
on disk. MODULE_3 finished its heterozygosity/ROH/θπ runs before the
maintenance window. The θπ scrubber data (per-window, in-progress) lives
on a different track — that goes into the Inversion Atlas's page 12, not
here. The Diversity Atlas is for **chromosome- and sample-level** views
of diversity that already exist as TSVs.

Concretely, the data is at:

```
/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/het_roh/
  02_heterozygosity/
    03_theta/
      *.win500000.step500000.pestPG       ← genome-scale θπ tracks (28)
      multiscale/*.pestPG                  ← 3 scales × 28 chroms
    04_summary/
      genomewide_heterozygosity.tsv        ← per-sample mean H
  03_ngsF_HMM/                             ← per-sample ROH calls
  04_roh_summary/
    per_sample_froh.tsv                    ← F_ROH per sample
    per_sample_roh_bins.tsv                ← short/medium/long counts
    per_chr_roh_summary.tsv                ← chrom × sample ROH map
  09_final_tables/
    master_summary.tsv                     ← joined per-sample table
```

All of these are ready to be transformed into atlas-side JSONs.

---

## What the atlas should show — recommended first cut

Five tabs, in priority order. Tab 1 alone makes the atlas useful even
if the rest never ship.

### Tab 1 — Sample diversity table (highest priority, easiest to build)

A single table, 226 rows, columns:
- sample_id, family, broodline (from `samples_qcpass.txt` + ancestry-Q)
- genome-wide H (from `genomewide_heterozygosity.tsv`)
- F_ROH (from `per_sample_froh.tsv`)
- ROH count by bin (short/medium/long)
- Largest single ROH (Mb)

Sortable by any column. Highlight the per-column outliers (top/bottom 5
percentile in colored stripes). Click a row to open the per-sample
view (tab 2). This is essentially the table that already exists as
`master_summary.tsv`; the work is just JSON-ifying it and rendering.

**Atlas file size**: ~30 KB JSON. Trivial.
**Time estimate**: 1–2 hours including R script + atlas renderer.

### Tab 2 — Per-sample diversity profile (medium priority)

Click a sample, see:
- 28-chromosome θπ ribbon (one row per chrom, length-scaled)
- Cohort percentile band (5–95) shaded behind the sample's line
- ROH segments overlaid as vertical strips (bin-colored)
- Click a chromosome → drill in (tab 3)

Reads `*.win500000.step500000.pestPG` for the cohort's mean ± percentile,
plus the per-sample θπ which we'd need to compute. **Actually** the
per-sample θπ isn't in MODULE_3's output by default — that needs an
extra realSFS+thetaStat pass per individual. This tab is therefore
gated on either:

- (a) reusing the SAF files already in `02_heterozygosity/01_saf/` (these
  are per-sample, so per-sample thetaStat IS feasible — just hasn't been
  run as a batch step), OR
- (b) deferring tab 2 to phase 2 of the atlas and shipping tab 1 + 3 + 4 + 5.

Recommended: option (b). Tab 2 needs new compute. Tab 1, 3, 4, 5 don't.

### Tab 3 — Chromosome ribbon (high priority, medium effort)

Pick a chromosome, see:
- Full-chrom θπ track at the 500 kb / 500 kb scale (cohort mean)
- Switchable to multiscale (5 kb, 10 kb, 50 kb)
- Region-of-interest callouts: top/bottom 5 percentile windows highlighted
- A second track below: ROH heat (sample × position, binned) — black
  pixel = many samples in ROH at that position; reveals genome-wide
  recombination cold spots, founder sweeps, etc.

This is the "diversity perspective" version of the Inversion Atlas's
page 1 — same chromosome-scale view but answering "is this region
diverse or depleted across the cohort" rather than "is there a
candidate inversion here". Reads `*.pestPG` + `per_chr_roh_summary.tsv`.

**Time estimate**: 4–6 hours. Mostly the ROH heatmap rendering — θπ
track is straightforward.

### Tab 4 — Family / broodline diversity comparison (medium priority)

A box-and-strip plot, x-axis = family or broodline (toggle), y-axis =
genome-wide H (or F_ROH, toggle), one dot per sample. Highlights
broodlines with elevated inbreeding or unusual diversity. This is
the "is broodline X depleted?" question for the breeding-program
collaborators.

Quentin's K=8 NGSadmix structure means broodlines are well-defined —
the assignment lives in the existing pop_struct outputs.

**Time estimate**: 2–3 hours.

### Tab 5 — About / glossary (low priority, last)

Standard "about this atlas" page. What's a θπ. What's an ROH bin.
Where the data came from. How this differs from the Inversion Atlas.
Mirror Population_atlas's tab 7 conventions.

**Time estimate**: 1 hour.

---

## Suggested direction — pick one of these for the next session

### Path A: ship tab 1 only (90 minutes, very low risk)

- Login node: write `STEP_DA01_emit_sample_diversity_json.R` that joins
  `master_summary.tsv` + `per_sample_froh.tsv` + ancestry-Q assignments
  → emits `Atlas/json/diversity_sample_table.json` (~30 KB).
- Atlas side: replace `Diversity_atlas.html`'s page1 scaffold with a
  real table renderer reading that JSON. Sortable, filterable, click-
  to-detail (the click-to-detail is a no-op until tab 2 ships).

This gives Quentin an immediately useful tool: a single sortable
table of 226 samples with diversity metrics. Pair with the R script
+ a runbook entry and that's a complete "MVP" of the Diversity Atlas.

### Path B: ship tab 1 + tab 3 (one full afternoon, medium risk)

Adds the chromosome-ribbon view. More valuable but more complex
because the ROH heatmap rendering is non-trivial (sample × position).
Reads existing `per_chr_roh_summary.tsv`.

### Path C: ship tab 1 + tab 4 (afternoon, low risk)

Tab 1 + family/broodline boxplot. Lower compute requirements than
path B. Better for "show breeding-program collaborators what we have"
demos. Strong story-telling angle.

**My recommendation**: Path A in the upcoming session. Ship tab 1, get
the file-format / JSON schema / atlas wiring right. Path B or C in a
follow-up once tab 1 is solid.

---

## What the JSON schema for tab 1 should look like

Suggested minimum (for tab 1):

```json
{
  "tool": "diversity_atlas_v1",
  "schema_version": 1,
  "generated_at": "2026-04-30T10:30:00Z",
  "cohort": "C_gariepinus_226_hatchery",
  "n_samples": 226,
  "samples": [
    {
      "sample_id": "...",
      "family": "...",
      "broodline": "K1",
      "ancestry_top": "K1",
      "ancestry_top_q": 0.92,
      "is_pruned81": true,
      "h_genomewide": 0.000497,
      "froh": 0.082,
      "roh_n_short": 142,
      "roh_n_medium": 18,
      "roh_n_long": 3,
      "roh_longest_mb": 14.2
    },
    ...
  ],
  "metadata": {
    "h_units": "per_site",
    "roh_short_max_kb": 1000,
    "roh_medium_max_kb": 5000,
    "froh_callable_genome_mb": 950
  }
}
```

Same conventions as the Inversion Atlas precomp JSONs (`tool`,
`schema_version`, `cohort`, `generated_at` at top level).

---

## What this handoff is NOT

- Not a θπ scrubber spec — that's the Inversion Atlas page 12 work.
- Not for the per-window dosage data — that's also Inversion Atlas.
- Not for the C. macrocephalus wild cohort — that's a future paper.
- Not for the F1 hybrid genome assembly — that's the Genome Atlas.

The Diversity Atlas is specifically about **per-sample and per-chromosome
diversity summaries** of the 226 hatchery cohort, derived from MODULE_3
outputs.

---

## Open questions for the next session

1. Is the per-sample θπ (tab 2) compute worth running, or punt to phase 2?
2. Should tab 1 link the F1 hybrid genome assembly numbers anywhere,
   or stay strictly cohort-scoped?
3. Do we want to embed the ancestry-Q stack-bar (per-sample) on tab 1,
   or keep tab 1 numerical-only and let tab 2 own the visual breakdown?
4. The Inversion Atlas's "promote to candidate" workflow has no analog
   here — should clicking a row do anything beyond opening the per-
   sample profile? (E.g. flag-as-outlier, write to a notes file?)

These questions are for the start of the next session, not for now.
