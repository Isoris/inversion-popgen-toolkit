# CODE_REVIEW_TO_FIGURE_OUT.md — running list of concerns

Not bugs. Not TODOs. Things where the code does what it says on the tin
but the *premise* feels shaky and deserves a second look before we trust
the downstream results. Append to the bottom; keep oldest at top.

Each entry:

- **Location** — file + line anchors so the concern is re-findable.
- **Concern** — the thing that bugs me.
- **Why it matters** — what downstream conclusion this could poison.
- **Current state** — what the code actually does today.
- **Resolution options** — what we could do about it, without picking.

---

## §1 — SV-prior seeding in `STEP_C01i_decompose.R`

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R`
- `seeded_kmeans()` L164-219 — the seeding + fallback logic
- sv_prior loader: L143-151 (sources `flashlight_v2/utils/...`)
- Tier-1 seed combiner: L115-133 (`lib_step03_seed_loader.R`,
  `combine_tier1_seeds` — sv_prior ∪ STEP03 with conflict drop)
- Main per-candidate call: ~L390 with `min_seeds_per_class = opt$min_seeds_per_class`

**Concern**

Seeding k-means from SV calls to recover the HOM_REF / HET / HOM_INV
trichotomy is a clever idea *conditional on the SV priors being well-
calibrated*. My sense from the previous runs is that they aren't —
DELLY/Manta misplaced the LG28 breakpoint by ~200 kb in the genome
assembly paper, and the per-chromosome sv_prior sets have had a lot
of noise/inconsistency across methods. I don't have a clean quality
score for sv_prior that I trust today. So the whole "seeded" path is
built on an input whose quality we haven't validated at the right
granularity for this use.

**Why it matters**

If sv_prior is systematically noisy, the seeded k-means will inherit
that noise as **structural bias** on the centroids, not random noise
that evens out across samples:

- If sv_prior over-calls HOM_INV (which is what a liberal SV caller
  does on a ~3 Mb inversion candidate), `center_inv` is pulled toward
  the HET distribution, the three centroids compress, and more real
  HET samples get assigned to HOM_INV.
- If sv_prior under-calls (strict caller), HOM_INV centroid anchors on
  only the *cleanest* inversions and splits the true HOM_INV cluster
  against the remainder.
- The silhouette scores come out looking *better* in the seeded case
  than the unsupervised case because the centroids are anchored, not
  because the clustering is more correct. So `silhouette > 0.25`
  stops being a trustworthy gate when seeding is active.
- The recombinant detection downstream (C01i_b) compares per-window
  class calls against the candidate-level decompose class. If the
  candidate-level call is biased by seeding, the per-window
  "deviation" signals shift correspondingly, and the recombinant
  gate inherits that.

The current `q6_group_validation = UNCERTAIN` default (capped to
UNCERTAIN unless quality flags clear it) mitigates this *if* we never
promote. But as soon as we start promoting, sv_prior noise propagates
into the promoted karyotypes.

**Current state**

- `seeded_kmeans()` only has a **quantity** guard, not a quality one:
  it falls back to unsupervised if
  `n_seeds < 10 || length(seed_ref) < 3 || length(seed_inv) < 2`
  (L170). There's no check on seed consistency, no check on caller
  agreement, no check on SV call support strength.
- sv_prior is treated as **Tier-1 truth** by the seed combiner, on
  equal footing with STEP03. Conflicts between sv_prior and STEP03
  are dropped (per the chat-12 "conflict drop" rule), but that's
  "agreement = trust" logic, which is right direction but doesn't
  catch systematic per-caller bias.
- The fallback path — `decomp_status = "no_seeding"` — is the
  conservative option that *does* exist (chat-12 removed the
  unsupervised fallback when seeding fails). But it only kicks in
  when there are too few seeds, not when the seeds are dubious.
- Quality flags on the output include `low_silhouette` and
  `weak_k3_fit`, but I don't emit a `seeded_centroid_drift` flag
  or similar — nothing warns downstream that seeding pulled the
  centroids hard.

**Resolution options** (not picking)

1. **Calibration study on LG28** — we know the karyotype (60/106/60)
   from independent phase-based analysis. Run decompose on LG28 with
   seeding on vs seeding off. Compare class assignments. If seeding
   flips > ~5% of samples relative to the trusted phase-based call,
   that's evidence seeding introduces bias at the candidate we know
   best. Cheap, low-stakes sanity check.
2. **Add a seed quality score** — score each sv_prior seed set by
   internal agreement across callers (DELLY/Manta/Clair3) and by
   supporting evidence counts. Require per-class mean seed quality
   to pass a threshold before seeding is allowed; fall back to
   `no_seeding` otherwise. This is the "don't trust what hasn't
   earned trust" path.
3. **Make seeded results advisory, not promotable** — leave seeding
   on for the cluster-shape heuristics but explicitly cap
   `q6_group_validation = UNCERTAIN` on any candidate where
   `seeded = TRUE` until a validation pipeline confirms the seed was
   honest. This is the "keep the clever idea, fence it" path.
4. **Orthogonal-signal validation before using seeded output** — for
   each candidate, cross-check the seeded class calls against a
   signal that doesn't depend on SV priors: e.g. GHSL haplotype
   divergence within-sample, or Fst contrast between called classes.
   Only emit `seeded_ok` if the orthogonal signal concurs.

**What would change my mind**

If a quick LG28-on-LANTA comparison (seeded vs unseeded vs trusted
phase-based) showed class-call agreement > ~95% in both seeded and
unseeded paths, the concern largely evaporates — seeding is then
just giving us faster convergence on the same answer, and sv_prior
quality matters less for this step than I'm imagining. I think this
is testable in one SLURM job once LANTA is back.

---

## §2 — Stale `sv_prior_dir` default in `lib_decompose_helpers.R`

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R` L39-40
- `resolve_decomp_paths()` function

**Concern**

The default resolution is:

```r
sv_prior_dir = args$sv_prior_dir %||% Sys.getenv("FLASHLIGHT_DIR",
  file.path(base, "flashlight_v2/cache"))
```

The `flashlight_v2/cache` path is outdated — it reflects an earlier
layout. Current running runs override via the `FLASHLIGHT_DIR` env var
or the `--sv_prior_dir` arg, so this doesn't actually bite in practice,
but the baked-in fallback is a footgun: anyone running without either
override silently reads from a path that isn't the current priors dir.

**Why it matters**

- Silent divergence between what people *think* is being loaded and
  what actually is. If the old path still exists with old data on
  LANTA, loads succeed but return stale priors — indistinguishable
  from correct behavior at the surface.
- Anyone copying `resolve_decomp_paths()` as a template for other
  scripts inherits the stale default.

**Current state**

- Env override works.
- CLI arg override works.
- Fallback is stale.

**Resolution options** (not picking)

1. Update the fallback to the current canonical priors directory.
2. Remove the fallback entirely — require either `--sv_prior_dir`
   or `FLASHLIGHT_DIR`. `stop()` with a clear error if neither.
3. Keep the fallback but have `resolve_decomp_paths()` emit a
   `message()` whenever the fallback fires, so stale-path usage is
   at least audible in job logs.

**What would change my mind**

Confirming the current canonical path on LANTA. Once we know the
right answer, option 1 is trivial.

---

## §3 — Stale `q_cache_dir` default in `lib_decompose_helpers.R`

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R` L41-42
- `resolve_decomp_paths()` function

**Concern**

Same shape as §2:

```r
q_cache_dir = args$q_cache_dir %||% Sys.getenv("Q_CACHE_DIR",
  file.path(base, "unified_ancestry/local_Q"))
```

The `unified_ancestry/local_Q` path is also out of date. Current
Engine B runs write the local-Q caches to a different layout (per
README §34: `<chr>.<sample_set>.local_Q_samples.tsv.gz` under K-split
subdirs).

**Why it matters**

- This is the same "silent stale fallback" risk as §2 — jobs that
  don't set `Q_CACHE_DIR` or `--q_cache_dir` silently read whatever
  is at the stale path, or nothing.
- **Direct blast radius on the newly-wired band_composition field**:
  the `compute_ancestry_div_hom_ref_vs_hom_inv` call lands at
  `STEP_C01i_decompose.R` L630-635 and is fed by `--local_q_dir`. If
  someone falls back to `resolve_decomp_paths()$q_cache_dir` instead
  of passing `--local_q_dir` explicitly, they'll hit the stale path.
  The helper's graceful-NA behavior will hide it (key comes back
  NA, looks normal), so no error surfaces — but every candidate gets
  NA for ancestry divergence, which LOOKS like "data not available"
  when it's really "data available at a different path I didn't
  find".

**Current state**

- `STEP_C01i_decompose.R` doesn't use `resolve_decomp_paths()` for
  local_q_dir today — it takes `--local_q_dir` directly. So the
  staleness is dormant for decompose specifically.
- But other scripts DO call `resolve_decomp_paths()`, and any that
  inherit `q_cache_dir` from it inherit the stale default.

**Resolution options** (not picking)

Same as §2 — update the path, remove the fallback, or make the
fallback audible.

**What would change my mind**

Same as §2: confirming the canonical path on LANTA.

---

## §4 — `try_load_sv_prior` reads the precomp priors, NOT phase-3 priors

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R` L217-225
- `try_load_sv_prior()` function

**Concern**

```r
try_load_sv_prior <- function(chr, sv_prior_dir) {
  fl_path <- file.path(sv_prior_dir, paste0(chr, ".rds"))
  if (!file.exists(fl_path)) return(NULL)
  tryCatch(readRDS(fl_path), error = function(e) {
    message("[lib] sv_prior load failed for ", chr, ": ", conditionMessage(e))
    NULL
  })
}
```

The `sv_prior` this loads is the **precomp-stage** SV prior, not the
higher-quality phase-3 prior. These are two different signals with
two different quality levels:

- **Precomp sv_prior** — computed during the 2c precomp phase from
  raw SV caller output (DELLY/Manta merged). Used early in the
  pipeline as a cheap coarse signal. Known to be noisier and to miss
  breakpoints by hundreds of kb.
- **Phase-3 sv_prior** — produced further down the pipeline with
  refinement steps on top. Better-calibrated but only available
  after phase 3 has run.

The name `sv_prior` and the function name `try_load_sv_prior` don't
distinguish which one this is; the implementation quietly uses the
precomp version.

**Why it matters**

This compounds §1 (the seeding concern). When C01i_decompose calls
`try_load_sv_prior` and feeds the result into `seeded_kmeans()`, it
is seeding k-means with the **noisier** of the two available priors,
even though the **better** prior exists elsewhere in the pipeline by
the time decompose runs. The seeded clustering's centroids are
anchored on the worse data. The quality-of-signal objection from §1
is sharper here because we *could be* using the better signal.

Separately, this is a naming hazard. Anyone reading decompose and
seeing `sv_prior` assumes the prior used is the phase-3 one unless
they trace the loader.

**Current state**

- The precomp prior is what gets loaded.
- `register_C01i_frequency_block`, `sv_prior_mode` tagging, and
  everything downstream that says "sv_prior" is referring to the
  precomp prior.
- There's no path in the codebase today for decompose to consume
  the phase-3 prior.

**Resolution options** (not picking)

1. **Rename for honesty** — `try_load_sv_prior_precomp()` (and the
   field `sv_prior_mode` → `sv_prior_precomp_mode`, the dir
   `sv_prior_dir` → `sv_prior_precomp_dir`). Surface the provenance
   in names. No behavior change; just no more assumption that
   "sv_prior" means "the best SV prior we have".
2. **Switch the source** — add a phase-3 prior loader and have
   decompose prefer phase-3 when it's available, fall back to
   precomp otherwise. Much bigger change; depends on phase-3
   being guaranteed available by the time decompose runs (which
   it should be, but confirm pipeline order).
3. **Hybrid** — load both, log which one was used per candidate,
   and include the source in the `sv_prior_mode` field so
   downstream can see "this candidate was seeded on precomp
   priors" vs "...on phase-3 priors".

**What would change my mind**

Measuring precomp-vs-phase-3 agreement on LG28. If the two priors
give ~95%+ the same seed set for LG28 samples, the distinction is
academic and option 1 (rename only) is enough. If they diverge a
lot, option 2 or 3 is warranted.

This connects directly to §1's proposed LG28 calibration study —
extending that study to compare precomp-seeded vs phase-3-seeded
vs unseeded vs trusted phase-based would answer both concerns in
one SLURM job.

---

## §5 — GC detector driver duplicates phase 6's diagnostic-SNP work (but not its tract-detection)

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/RUN_gene_conversion_detector.R` (driver wired 2026-04-24 chat D)
- `inversion_modules/phase_7_karyotype_groups/proposal/gene_conversion_detector.R`
  - `diagnostic_snp_mask()` L199-219
- `inversion_modules/phase_6_breakpoint_refinement/01_dosage_signal.R`
  - `compute_informative_markers()` L133-170
  - `BP01_PARAMS` L111-120

**Concern**

The GC detector and phase 6's `01_dosage_signal.R` BOTH compute a
diagnostic / informative SNP mask from the same inputs (per-candidate
dosage + karyotype groups). The driver I wrote for the detector
re-loads the per-chromosome dosage files and re-derives the mask
from scratch. Phase 6 has already done this work — and written the
result to `dosage_informative_markers.tsv.gz`.

The architecturally right home for GC detection is as a post-pass on
phase 6's output, not a standalone driver after decompose.

**Why it matters**

- I/O duplication: per-chromosome dosage files are large. Doing the
  load twice per candidate is wasteful on LANTA.
- Subtle divergence risk: the two sites compute effectively the same
  mask but on different *scales*, so changes only hit one.
  - Phase 6: `delta = h2_mean - h1_mean` on **raw dosage** (0–2
    scale); `informative` iff `|delta| >= 1.0`.
  - GC detector: `af_ref = colMeans(...)/2` on **allele frequency**
    (0–1 scale); `diagnostic` iff `|af_ref - af_inv| >= 0.5`.
  - Algebraically equivalent (0.5 AF = 1.0 dosage), but if someone
    tunes one threshold they have to remember the other exists.
- The detector's `snp_info` parameter is accepted but the driver
  leaves it NULL (no depth_z or per-SNP call-rate available from
  the BEAGLE-majmin dosage files). Phase 6's input path could
  plausibly enrich the mask with its own QC signals, if they
  matter.

**Important scale nuance — what IS and ISN'T shareable**

The two scripts are NOT doing the same thing at the tract-length
layer, and the review should not pretend otherwise:

| Step | Phase 6 (`01_dosage_signal`) | GC detector |
|---|---|---|
| Diagnostic-SNP mask | `compute_informative_markers` | `diagnostic_snp_mask` |
| Length gate | `min_core_size = 20 SNPs`, `max_extend = 5000 SNPs` | `min_run = 2`, `max_flagged_snps = 10`, `max_span_bp = 20 kb` |
| Purpose | Locate the **inversion** (Mb-scale block) | Locate **short arrangement-discordant tracts** inside it |
| Algorithm | Correlation core-block + correlation-extend | Per-SNP run-length with tolerance |
| Scale | Thousands of SNPs | 2–10 SNPs |

Phase 6 explicitly looks for the LARGE structure; the GC detector
explicitly looks for SHORT discordant runs inside.

**The GC detector's length gates are set by TWO DIFFERENT logics
that should not be confused:**

1. **Lower bound — literature-justified (biology).** `min_run = 2`
   and the Harris 1993 references (18–774 bp typical, ~5 kb
   atypical upper) justify the detection *floor*. The chat-12
   rewrite of this detector was motivated specifically because
   the old windowed detector had a ~60 kb floor — ABOVE real GC
   biology, so genuine short tracts were invisible. We care about
   the lower bound being low enough to catch real biology. This
   is where the Harris citations live.

2. **Upper bound — operational handoff (NOT biology).**
   `max_span_bp = 20000` and `max_flagged_snps = 10` are NOT
   claiming real GC tracts stop at 20 kb. They are an arbitrary
   routing boundary: things shorter than this go to the GC
   detector; things longer go to C01j's recombinant machinery
   (which handles double crossovers and single crossovers at the
   Mb scale). The detector's own header makes this explicit:
   "Longer events are recombinants, handled by C01j — discard
   here." Nothing about that number is a biological ceiling; it's
   a "don't double-count the same event in both pipelines" line.

   **Historical context — the upper bound came from the dosage
   heatmap, not the literature.** The chat-12 rewrite was
   motivated by wanting to match what's visually separable in
   the sample × SNP dosage heatmap:
   - SHORT HORIZONTAL STRIPE (one sample, 1–3 SNP columns) →
     GC-like tract, goes to this detector
   - VERTICAL STRIPES (many samples, one column) → non-diagnostic
     SNP or paralog mismapping, filtered by diagnostic_snp_mask
   - WIDE HORIZONTAL FLIPS (one sample, many SNP columns) →
     recombinant, handled by C01j
   The 20 kb number is "at our SNP density (~1.6 kb/SNP at 9×),
   anything wider than ~12 SNPs is clearly a recombinant-sized
   stripe flip rather than a fleck." It's a heatmap-visual
   threshold operationalized in bp, not a biological tract-length
   claim.

The honest framing for the upper bound: we don't know, from
sequence alone, whether a given short arrangement-discordant tract
is gene conversion, a short double crossover, paralog mismapping,
or coincidental IBS. The manuscript wording ("arrangement-discordant
IBS tracts consistent with gene conversion") acknowledges this. The
20 kb cutoff is a pipeline routing decision, not a biological
claim, and shouldn't be justified in the manuscript as if it were.

So the consolidation story is targeted, not wholesale:

- **Share** Steps 1–2 (QC mask + diagnostic-SNP mask). Same math,
  same data, currently duplicated.
- **Do NOT share** Steps 3–7 (per-SNP flagging, run-length
  aggregation, length gate, confidence, direction, summary). Those
  operate at a different scale and are the actual IP of the GC
  detector.

**Aside — why not swap the two detection algorithms entirely?**

Reasonable question: phase 6 uses correlation-core-block + extend;
the GC detector uses per-SNP run-length with tolerance. Could
phase 6 use run-length and the GC detector use correlation?

No. The two algorithms are paired with their tasks by algorithmic
structure, not historical accident:

- Phase 6 is looking for a **cross-sample** signal (every HOM_INV
  sample deviates from every HOM_REF sample at the same markers,
  across the whole inversion). Correlation-core-block extracts
  exactly this — marker-to-marker correlation across all samples.
  Per-SNP run-length would scan each sample independently and then
  have to aggregate, throwing away the cross-sample correlation
  structure that makes the signal separable from noise.
- GC is a **single-sample** phenomenon — one carrier's short tract
  of opposite-arrangement-looking genotypes. Per-SNP run-length
  scans each sample independently, which is the right shape.
  Correlation-core-block averages across samples, which is
  precisely what erases single-sample deviations (one GC tract in
  one sample out of 226 is a 1/226 perturbation at that marker —
  invisible to cross-sample correlation).

Swapping would put each algorithm on the task it's worst at. The
*sharing* in Steps 1–2 is about shared inputs, not shared
algorithms. The algorithms correctly stay where they are.

**Current state**

- Driver works, JSON output feeds C01i_b correctly, registry block
  emits the 3 schema keys. Scanner state includes the new WIRED
  contract. No correctness issue in the current setup.
- Duplication is latent — it'll bite when someone changes one
  threshold without the other, or when the BEAGLE-majmin QC is
  revisited.

**Resolution options** (not picking)

1. **Refactor the detector's Step 2 out** as a standalone function
   `build_diagnostic_snp_mask(X, groups, delta_min_af = 0.5)` that
   phase 6 ALSO calls (replacing its inline `compute_informative_markers`
   delta computation with a thin wrapper). Both sites then share the
   same mask implementation; thresholds can still be different but
   they're in one place.
2. **Move the GC detector's invocation into phase 6's orchestrator**.
   After `01_dosage_signal.R` has computed `X` and the informative
   mask, call `detect_cohort_gene_conversion_v2()` on the same
   `X` + an explicit diagnostic-mask argument. Writes the
   `gene_conversion_tracts.json` file in the same per-candidate
   output dir phase 6 already uses. The standalone driver from chat
   D stays as a re-run tool for when you want to redo GC without
   redoing phase 6.
3. **Nothing**. Accept the duplication as the cost of keeping the
   two phases decoupled.

**What would change my mind**

Measuring on LG28: does phase 6's `informative` mask agree with
the detector's `diagnostic_snp_mask` on the exact same sample-class
assignment? They SHOULD — one just divides by 2 — but a test
would rule out surprise divergences from NA handling or sample
intersection differences.

---

## §6 — Why breakpoint refinement before final decompose? Ordering is unchallenged, possibly inherited

**Location**
- `inversion_modules/phase_6_breakpoint_refinement/01_dosage_signal.R`
  (runs before Phase 7)
- `inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R`
  (runs after)
- C01d `register_C01d_keys` in `utils/registry_key_helpers.R` — upstream
  rough k-means on sim_mat PC1 that phase 6 currently consumes

**Concern**

The pipeline is **breakpoint-first, karyotype-second**. Phase 6
refines breakpoints using a set of karyotype groups that come from
phase 4's rough candidate scoring, not from phase 7's proper
decompose. The proper decompose runs AFTER phase 6 and never feeds
back.

There's circularity under the ordering:
- Good breakpoints want good groups (so `delta_i` is a real
  HOM_REF-vs-HOM_INV contrast, not a noisy one).
- Good groups want good breakpoints (so the sample-by-SNP matrix
  being clustered is tight to inverted sequence, not diluted by
  flanks).

The current order resolves the cycle by doing a **coarse pass
early**: phase 4's rough k-means → phase 6 breakpoint refinement
on those rough groups → phase 7 proper decompose on the refined
interval. The pipeline executes the cycle exactly once, starting
from the coarse end.

**Why it matters**

If phase 7's decompose would give meaningfully different groups
from phase 4's rough k-means, then phase 6 is computing breakpoints
from the wrong groups. Breakpoint localization quality depends
directly on group purity (`delta_i` is a group-mean contrast).

The LG28 breakpoint is currently reported with CI widths in the
low-kb range — tight enough that group-purity perturbations at the
few-sample level could shift the breakpoint noticeably. If three
samples flip from HET to HOM_INV between phase 4's call and phase
7's decompose (plausible given §1's concerns about seeded vs
unseeded k-means), the dosage contrast at marginally-informative
SNPs near the breakpoint edges is going to shift, and so will
the correlation-extend stopping point.

**Current state**

- Phase 6 depends on groups via `reg$samples$get_groups_for_candidate`
  (01_dosage_signal L343). Those groups are whatever's in
  sample_registry at phase-6 runtime. In the live pipeline that's
  phase 4's call (possibly augmented by whatever C01d writes).
- Phase 7 decompose does its own k-means on PC1 loadings over the
  phase-6-refined interval. It does NOT rewrite the breakpoints.
- There is NO second pass — phase 6 runs once, phase 7 runs once,
  and if they disagree on groups, the disagreement is invisible.

**Why the current order might be there (honest guesses, none confirmed)**

1. **Empirical:** someone measured and found group calls don't
   change enough between phase 4 and phase 7 to move breakpoints.
   Breakpoint-first is then fine. I don't see this documented.
2. **Historical:** phase 6 was built first; phase 7 was added
   later; nobody re-derived ordering from first principles.
   Plausible given the chat-12 / chat-14 rewrites of phase 7 with
   phase 6 sitting fixed.
3. **Dependency convenience:** breakpoint-first lets phase 6
   consume phase 4's rough groups directly, without waiting on a
   phase-7 output that doesn't exist when phase 6 starts. An
   ordering that would require phase 7 → phase 6 would need
   decompose to accept a ROUGH interval (which it can; it's group
   assignment, not breakpoint localization) and then feed its
   groups back up to phase 6.

My gut says (2) + (3) more than (1), but it's genuinely unclear
from the code alone.

**Resolution options** (not picking)

1. **Measure, then decide.** For LG28, run the pipeline as-is and
   record the group call from C01d (phase 4) AND the group call
   from C01i_decompose (phase 7). Count disagreements. If
   disagreements are < ~3 samples, (1) is retroactively confirmed
   and the ordering is fine. If more, there's a real issue.
2. **Two-pass breakpoint refinement.** Run phase 6 → phase 7 as
   usual, then re-run phase 6 with phase 7's groups. Compare
   breakpoints between the two phase 6 runs. If they shift > 2 kb,
   the second pass is buying something. The infrastructure cost is
   low (phase 6 is not slow) and this gives a breakpoint-stability
   error bar almost for free.
3. **Decompose-first architecture.** Decompose on the rough
   interval (phase 4 output). Run phase 6 with decompose's groups.
   One breakpoint pass, from the best-available groups. Cleanest,
   but requires proving decompose is tolerant enough of rough
   intervals that skipping phase 6's contribution to group quality
   doesn't hurt decompose. Likely fine based on how PC1 dominates
   inside the inversion, but needs checking.

**What would change my mind**

Option 1 alone. On LG28 where we have a trusted phase-based
karyotype, compare C01d's group call to C01i_decompose's group
call sample-by-sample. If they agree within ~1-2%, the current
ordering is empirically fine even if it's not first-principles
justified. If they diverge more, phase 6's breakpoints are
anchored on a suboptimal signal and option 2 or 3 is worth the
effort.

This is the same LG28 SLURM job that would settle §1 and §4. One
run, four answers: (§1) seeded vs unseeded decompose class calls,
(§4) precomp vs phase-3 SV priors, (§5) phase-6 mask vs
GC-detector mask on the same samples, (§6) phase-4 vs phase-7
group calls. Bundle them.

---

## §7 — GC detector flagging assumes fixed dosage polarity; diagnostic-mask does not  [RESOLVED 2026-04-24 chat D]

**Status: FIXED.** The polarity-aware flagging is now the shipped
behavior in `scan_one_sample` (L274+ after the chat-D rewrite).
`diagnostic_snp_mask` carries af_ref/af_inv as attributes; `v2` pulls
them off and threads them into `scan_one_sample`, which compares each
sample's dosage to the per-SNP group means rather than to fixed 0.5 /
1.5 thresholds. The fixed-threshold code path is retained as a fallback
for callers that don't pass af_ref/af_inv (with a warning), to keep
back-compat for any edge caller. Test coverage added in tests 11
(all-reverse-polarity cohort, asserts 0 false tracts) and 12 (mixed
polarity with one real tract) of `test_gc_detector.R`. Test 13 asserts
the new diagnostic_snps block field is populated.

Still open for empirical follow-up on LG28: measure the fraction of
diagnostic SNPs that are reverse-polarity on the real cohort, and
compare tract counts under the old vs new rule — if the old rule was
producing many false positives in practice, the new tract counts
should drop meaningfully. Record the number for the Methods write-up.

Original concern text below for context:

---

**Location**
- `inversion_modules/phase_7_karyotype_groups/proposal/gene_conversion_detector.R`
  - `scan_one_sample()` L241-288 (flagging)
  - `diagnostic_snp_mask()` L199-219 (mask)

**Concern**

The detector has two stages that disagree on whether polarity matters.

**Stage 1 — diagnostic_snp_mask (polarity-agnostic, correct):**

```r
af_ref <- colMeans(dosage_mat[ref_idx, ]) / 2
af_inv <- colMeans(dosage_mat[inv_idx, ]) / 2
diag_delta <- abs(af_ref - af_inv)
# diagnostic iff |delta| >= 0.5
```

This correctly treats SNPs where HOM_REF is at dosage ≈ 2 and SNPs
where HOM_REF is at dosage ≈ 0 as EQUALLY diagnostic, provided the
two groups disagree strongly. Good.

**Stage 2 — scan_one_sample (fixed polarity, incorrect in general):**

```r
if (baseline_class == "HOM_REF") {
  state[d <= 0.5]  <- "consistent"   # assumes HOM_REF → low dosage
  state[d >= 1.5]  <- "flag"         # assumes HOM_INV-looking = high dosage
  ...
}
```

This ASSUMES every diagnostic SNP has HOM_REF at low dosage and
HOM_INV at high dosage. That is true when ANGSD's major/minor call
polarises the data so the common-cohort-allele is the major — which
is typical when the REF arrangement is much more frequent than INV.

It is NOT true when the two arrangements are comparably common (e.g.,
LG28 at 60/106/60 on a pure C. gariepinus hatchery cohort). At a
diagnostic SNP where HOM_REF samples carry what ANGSD called the
"minor" allele (because the overall cohort allele balance was close
to 50/50 and rolled the other way), HOM_REF samples have dosage ≈ 2
legitimately. The detector's `scan_one_sample` would flag all of them
as "INV-looking" = FALSE POSITIVE.

**Why it matters**

- In balanced cohorts, an unknown fraction of diagnostic SNPs may
  have reverse-polarity. Every HOM_REF sample would produce a
  spurious "tract" at every reverse-polarity diagnostic SNP.
- This is exactly the kind of bug that inflates tract counts
  quietly. Manuscript tables reporting "N samples with GC-like
  tracts" would be over-counting by an amount that depends on the
  cohort balance, without any signal that anything's wrong.
- The detector's per-candidate diagnostic mask means the polarity
  is locked in at mask-build time. But the flagging stage doesn't
  read the polarity — it reads the raw dosage with fixed thresholds.

**Current state**

- Code at L274-284 is as described above.
- No polarity-alignment check exists between mask and scan.
- No test in test_gc_detector.R exercises reverse-polarity SNPs
  (worth checking — the synthetic fixtures may all be same-polarity
  by accident).

**Resolution options** (not picking)

1. **Make flagging polarity-aware.** Compute per-SNP group means
   once (already done in diagnostic_snp_mask — just pass them
   through). For each sample at each diagnostic SNP, flag if its
   dosage is closer to the OTHER group's mean than to its own. This
   is the principled fix. Small code change; touches only the
   flagging thresholds, not the tract-detection algorithm.
2. **Pre-polarise dosages at mask build time.** For each diagnostic
   SNP, sign-flip dosages so HOM_REF is always at low end and
   HOM_INV at high end. Then fixed thresholds ≤ 0.5 / ≥ 1.5 are
   correct again. This is less invasive but couples the mask and
   scan stages more tightly.
3. **Measure first.** Run the detector on LG28 with both the
   current fixed-threshold rule and a polarity-aware rule. Report
   tract counts. If they match within a few percent, the current
   code is fine in practice. If they diverge, fix per option 1.

**What would change my mind**

Empirical measurement on a balanced cohort (LG28 at 60/106/60 is
the obvious test case). If the fraction of reverse-polarity
diagnostic SNPs is small or if their spurious flags don't survive
the run-length gate (e.g., reverse-polarity SNPs are scattered such
that they rarely form consecutive runs of ≥ 2), the real-world
impact may be minimal. But this needs measuring, not assuming.

Note: this concern surfaced from re-reading the detector during the
DESIGN_RATIONALE writeup. It was NOT caught by the chat-12 rewrite
or chat-17's tests. It's the kind of bug that hides because the
diagnostic-mask code reads plausibly, the scan code reads plausibly,
and their incompatibility is only visible when you trace an
individual SNP with unusual polarity through both stages.

