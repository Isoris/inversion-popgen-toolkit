# HANDOFF — 2026-05-06 — band-tracking + age-method specs (next-chat queue)

**Date**: 2026-05-06
**Atlas main file in production**: `Inversion_atlas.html` (76,034 lines,
turn 164 baseline — `Atlas_turn164_2026-05-05_tar.gz`)
**Working tree from this chat**: `Atlas_band_consensus_2026-05-05.tar.gz`
(modular split — independent of the production atlas; only adds files
under `Atlas/shared/band_tracking/` and `Atlas/tests/`)
**Project**: `MS_Inversions_North_african_catfish` — 226-sample pure
*C. gariepinus* hatchery cohort, LANTA HPC.

---

## 0. What this chat actually shipped

### 0.1 Code (modular Atlas tree, 991/991 tests passing)

Five-tarball progression. Each is **additive** on top of the previous —
unzip them in order over a fresh `Atlas/` tree and you get the full
state. Each tarball ships into `Atlas/shared/band_tracking/` and
`Atlas/tests/`; **none touch `Inversion_atlas.html`** (the production
atlas).

| # | Tarball | What it adds | Tests added |
|---|---|---|---|
| 1 | `Atlas_FINAL_2026-05-05.tar.gz` | merge handoff: 4 sub-atlas HTMLs, slot registry, escape parity | 754 |
| 2 | `Atlas_band_tracking_2026-05-05.tar.gz` | single-band tracking + het skeleton + iv calls | +90 → 844 |
| 3 | `Atlas_band_layers_2026-05-05.tar.gz` | trajectory + projection + karyotype-model combiner | +73 → 917 |
| 4 | `Atlas_band_consensus_2026-05-05.tar.gz` | vote evidence + per-band view + bruteforce + 6 classes | +74 → 991 |
| 5 | (this chat — specs only) | spec amendment + new BUSCO 4D spec, no code | (no test change) |

**Test count: 991/991 pass on synthetic data. Zero tests on real catfish.**

### 0.2 Specs (this chat, no code patches into atlas)

Two new spec docs in `Atlas/specs_todo/`:

1. **`SPEC_inversion_age_atlas_surface_AMENDMENT.md`** — drops absolute-My
   displays from the existing dXY age spec (turn 117), replaces with
   ranking voice for Methods 1 + 2, adds Row C (Method 3, BUSCO 4D
   three-μ brackets) and Row D (supporting region-popstats traces).
2. **`SPEC_busco_4d_age_brackets.md`** — defines Method 3 end-to-end:
   new LANTA producer `STEP_C01f_e_emit_busco_4d_age.py`, JSON shape,
   atlas display, manuscript bundle integration, voice rules.

Both are in this tarball. **No code shipped this turn.**

---

## 1. The age-method framework (final, locked)

After multi-turn discussion concluding that all dating methods are
unreliable at hatchery Ne ≈ 20 with sweepstakes reproduction (millions
of eggs / spawn), we settled on three methods reported **alongside each
other** so the reader can triangulate:

- **Method 1** (relative between-inversions): rank inversions on each
  chrom by between-arrangement dXY. Page 5 catalogue column. **No My.**
- **Method 2** (relative within-arrangement): rank arrangements within
  one inversion by within-class π. Page 3 panel Row B. **No My.**
- **Method 3** (absolute via BUSCO 4D + three-μ brackets): dXY at 4D
  sites in BUSCO genes inside the inversion, divided by 2μ for three μ
  values (1×10⁻⁹, 3×10⁻⁹ Liu 2023, 9×10⁻⁹). Page 3 panel Row C. **My
  reported as a triple, never alone.**

Plus the **already-shipped** GDS panel (cheat30, turn 119): OLD/INTER/
YOUNG bins + single/recurrent origin, on page 3.

Plus **Row D supporting stats** (NEW): dXY trace, θπ per arrangement,
Tajima's D shape (Guerrero 2012 age proxy classes), all per-window over
the inversion interval. From `region_popstats.c`, displayed via a new
`region_popstats_v1.json` joiner pipeline (`STEP_C01f_f`, ~80 lines).

**The atlas displays four readings of the same locus**: GDS bin (already
shipped), Method 1 ranking, Method 2 ranking, Method 3 bracket. Plus
diagnostic stats (Row D). When all agree → high confidence. When they
disagree → diagnostic (selection, recent introgression, method-specific
breakdown).

---

## 2. What the next chat should do — three options

### Option A — Implement the age layer (Methods 1 + 2 + 3 + Row D)

Pure atlas-side work, ~880 LOC patch into `Inversion_atlas.html`.
No HPC dependency for Methods 1 + 2 + Row D (LANTA already produces
the inputs). Method 3 ships with empty state until LANTA runs the new
`STEP_C01f_e` pipeline.

**Risk**: 880 LOC into a 76k-line monolith without browser-side
verification. The two specs are detailed enough to make the work
mechanical, but a function rename in turns 165+ could trip me. **Mitigate
by asking Quentin whether turn 164 is still current before patching.**

**Order** (per AMENDMENT §8):
1. `region_popstats_v1` loader (~50 LOC)
2. Page-3 Row D supporting stats panel (~150 LOC)
3. `inversion_age_v1` loader (~80 LOC)
4. Page-3 Row A Method 1 (~40 LOC)
5. Page-3 Row B Method 2 (~60 LOC)
6. Page-3 Row C Method 3 (~80 LOC, gated empty state)
7. Page-5 "rel age" column (~80 LOC)
8. Manuscript bundle (~80 LOC)
9. Help-page row (~10 LOC)
10. Tests (~250 LOC, ~35 assertions)

### Option B — Run the band-tracking pipeline on real LG28

The 991-test synthetic suite has zero contact with real fish. The
LG28 prototype is at 15.115–18.005 Mb, karyotype 60/106/60. Wire
`consensus_partition` to the existing precomp JSON for that locus,
look at the output. Adjust thresholds if needed. Decision-relevant
data for everything that comes after.

### Option C — Write `het_location.js`

Parked module from turns 9–10. Inputs all available now
(`kt_infer_macro_band_groups`, `consensus_partition`, `het_track_skeleton`).
Seven-class het-position classifier per the contract negotiated in
turn 14. ~400 LOC + ~400 LOC tests, fresh module, no integration risk.
Clearly specified, mechanical to write.

### My recommendation, in priority order

1. **B first.** The synthetic suite is detailed but disconnected from
   real data. Running it on LG28 tells us whether the
   thresholds (τ=0.7, projectionMinAgree=0.4, ambiguousBandThr=0.5,
   etc.) work on hatchery data. Five-minute exercise on Quentin's
   side; saves 10× more debugging downstream.
2. **A second.** Manuscript-relevant; the spec is the most thoroughly
   drafted item in the queue. Before patching, confirm turn 164 is
   the current production atlas (or get the latest tarball).
3. **C third.** Useful but needs A and B done first to be calibrated
   against real data.

If the next chat is short on time and Quentin says "pick", lean B
(empirical) over A (more code), because Option B might reveal threshold
issues that should be fixed before A's tests are written.

---

## 3. Open decisions still pending Quentin's input

1. **Atlas baseline.** Is `Atlas_turn164_2026-05-05_tar.gz` still the
   current production atlas, or has Quentin shipped turns 165+ since?
   The age-layer patch (Option A) needs this confirmed before code is
   written. **If turns 165+ exist, ask for the latest tarball.**

2. **JSON fixture source.** For Option A's tests, do we have a real
   `inversion_age_v1.json` from a `STEP_C01f_d_emit_age_json.py` run
   on LANTA, or do we synthesize one matching the spec? (Synthesizing
   is fine for renderer tests but may diverge from real-data shape.)

3. **HPC pipeline order.** Method 3 (BUSCO 4D) requires
   `STEP_C01f_e_emit_busco_4d_age.py` to actually run on LANTA. The
   atlas can ship with empty-state Row C, but the manuscript story
   needs Method 3 numbers. Quentin's call when to run that pipeline.

4. **Tajima-shape classifier source.** Row D displays a Tajima D shape
   class from C01f_c §4e. Confirmed C01f_c emits this? If not, also
   need that step's output to be exposed in the per-candidate JSON
   (covered in `SPEC_busco_4d_age_brackets.md` §10's joiner).

5. **`het_location.js` per-layer iteration.** Original turn-14 contract
   said "one verdict per layer for MULTI_LAYER_STRUCTURE cases." If
   Option C is picked, confirm whether per-layer iteration is built
   in from the start or deferred to a follow-up.

---

## 4. Verbatim block from this chat

The user requested this section be preserved literally in the handoff,
as a record of where we paused before the spec turn:

> "Continue" is ambiguous here, and I'd rather ask than guess wrong.
> Five honest possibilities for what you might mean, in roughly
> descending likelihood:
>
> **1. "Now write Pass B."** Earlier I split the work into Pass A
> (foundation) and Pass B (multi-layer detection + extra QC). I ended
> up doing both in Pass A — the consensus tarball already has all six
> classes, all five QC fields, multi-layer detection, the lot. There's
> no Pass B left to write under the original plan.
>
> **2. "Now write `het_location.js`."** This is the parked module from
> turns 9–10. Its inputs are now all available: macro-state groupings
> (`kt_infer_macro_band_groups`), partition consensus + QC
> (`consensus_partition`), het skeleton tracking (`het_track_skeleton`).
> Writing it would tie everything together into the seven-class het
> location classifier. Earlier I argued for deferring this until you'd
> seen real outputs from the lower layers; I'd still argue that, but
> if you want to move forward I can do it.
>
> **3. "Now write `dos_build_macro_plus_het_model`."** The dosage model
> that consumes macro-state labels + het location and emits per-sample
> karyotype calls under the new framework. Replaces or supplements
> `iv_call_samples_from_skeleton`. Same reasoning as #2 for why I'd
> argue to wait.
>
> **4. "Now write the visualization."** Earlier you said the
> per-sample-lines panel "directly relates to" the band-centric voter
> view, "for visualizing." A page module or HTML artifact that
> *renders* `consensus_partition` output — probably one panel per
> band showing the visiting/excluding voters across the chromosome,
> plus a locus-level partition strip showing top-N partitions with
> the QC bar, plus the per-sample-lines panel re-colored by
> macro-state. This would let you actually look at what the algorithm
> finds.
>
> **5. "Just keep going on whatever's next."** You're tired and want
> me to pick. In that case I'd pick #4 — the visualization. The math
> is in much better shape than my ability to know which threshold
> defaults match real data. Until you can *see* the algorithm's
> output on a real chromosome, every additional algorithm I write
> risks being tuned for synthetic fixtures rather than fish.
>
> A sixth possibility worth flagging: **none of the above**. Look at
> where we are right now — over multiple turns we've gone from "the
> merge handoff" to "macro-state grouping" to "voter consensus with
> multi-layer detection." That's a lot of unrun code. The single most
> valuable thing might be **stopping** and running what exists on a
> real LG28 candidate. If I keep adding layers without ground truth,
> I'm building an increasingly elaborate structure on assumptions I
> can't verify. The deliverables in `/mnt/user-data/outputs/` are 991
> tests on synthetic data; zero of them are on real catfish.
>
> My recommendation, in priority order:
>
> 1. **Pause coding. Run the current pipeline on one real candidate**
>    (e.g., the LG28 prototype: 15.115–18.005 Mb, 60/106/60 karyo
>    from the registry). That tells us whether the thresholds need
>    tuning, whether the multi-layer signal you described shows up
>    where you expect, and whether the voter view is actually
>    informative on data that wasn't built to make it informative.
>
> 2. **If you don't want to pause, write the visualization (#4).**
>    It's the lowest-risk way to make progress — you can't see
>    anything go wrong on real data until something is rendered.
>
> 3. **If you want to keep adding algorithm, write `het_location.js`
>    (#2).** It's the most clearly-specified parked item and it
>    consumes outputs that are now all available.
>
> What do you want me to do? If you say "pick" I'll do #4.

After that exchange, Quentin chose to push for the age specs instead.
This handoff records the choice point so the next chat can re-enter
the parked work (visualization / het_location / dos_model) when ready.

---

## 5. Constants worth carrying forward

(In addition to whatever the next chat's system context already has.)

### Three catfish cohorts — never conflate

1. F₁ hybrid (*C. gariepinus* × *C. macrocephalus*) — genome assembly
   paper only.
2. **226-sample pure *C. gariepinus* hatchery cohort on LANTA — current
   inversion work.** K clusters reflect hatchery broodline structure,
   NOT species admixture.
3. Pure *C. macrocephalus* wild cohort — future paper.

### KaryoConstants

- NAToRA-pruned 81 unrelated samples
- NGSadmix K=8 on pruned set
- `mergeThr 0.85`, `minNGroup 5`, `minNWin 5`, `alpha 0.001`
- `state.k 3`, `purity 0.80`, `ambiguous 0.5`
- KING/Manichaikul thresholds: 1st≥0.177, 2nd≥0.0884, 3rd≥0.0442
- Reference: `fClaHyb_Gar_LG.fa` (28 pseudochromosomes, ~964 Mb)

### LG28 prototype

- 15.115 – 18.005 Mb, 28-band assembly
- karyotype 60/106/60 (HOM_REF / HET / HOM_INV)
- shelf Fst ~0.308 vs flanks ~0.04
- The locus to test the band-tracking pipeline on first

### Mutation rate (final, locked)

- **μ_year = 3×10⁻⁹/site/year**
- Liu et al. 2023, *Mar Life Sci Technol*, electric catfish
  (*Malapterurus electricus*), r8s-estimated from Siluriformes
  phylogeny calibrated against teleost fossils
- doi:10.1007/s42995-023-00197-8
- Used **per year directly** — generation time NOT assumed
- Method 3 brackets this with μ_low = 1×10⁻⁹ and μ_high = 9×10⁻⁹

### Band-consensus algorithm thresholds (defaults)

- trajectory τ = 0.7 (initial groups)
- projection split threshold = 0.4 (visited_jaccard below → split)
- projection merge threshold = 0.85 visited_jaccard AND 0.85 purity_cosine
- primary partner threshold = 0.7
- secondary partner = [0.4, 0.7]
- ambiguous-band threshold = 0.5 hidden_regime_residual
- δ for top-N adaptive = 0.10
- K_cap for bruteforce = 10
- hard-fail K = 13

These are synthetic-data defaults. **Real-data tuning is exactly the
question Option B (run on LG28) would answer.**

---

## 6. Files in this tarball

```
Atlas/specs_todo/
├── SPEC_inversion_age_atlas_surface_AMENDMENT.md   (NEW, this chat)
└── SPEC_busco_4d_age_brackets.md                   (NEW, this chat)

Atlas/HANDOFF_2026-05-06_age_specs_and_band_tracking.md  (this file)

(plus the existing band_consensus tree from the prior tarball:
 Atlas/shared/band_tracking/ and Atlas/tests/, unchanged)
```

The two specs are placed in `specs_todo/` (per the existing convention:
spec'd but not implemented). When the age layer ships into the
production atlas, move both specs to `specs_done/` after archiving.

---

## 7. Honest framing

**What's solid:**
- Two age specs are the most thoroughly drafted items in the queue.
  Method 3 has a complete LANTA-side pipeline definition; both atlas
  surfaces are pinned to specific JSON schemas; voice discipline rules
  are explicit.
- The 991-test band-tracking suite gives a strong synthetic foundation.
- All decisions about absolute-My (drop) vs three-μ (keep with brackets)
  are recorded with reasoning, so future readers won't relitigate.

**What's risky:**
- Zero contact with real catfish data. Every threshold is synthetic.
- The age layer's atlas-side patch (Option A) is large (~880 LOC into
  76k). High value; non-zero patching risk.
- Method 3 depends on a LANTA pipeline that doesn't exist yet
  (`STEP_C01f_e_emit_busco_4d_age.py`). The atlas can ship Row C with
  empty state, but the manuscript needs the numbers.
- `region_popstats_v1.json` (Row D supporting stats) needs a small new
  joiner script (`STEP_C01f_f`, ~80 lines) on LANTA. Specified but
  not yet written.
- The chat is exhausted. Threshold tuning, real-data review, multi-
  layer detection on real LG28, and the parked het_location /
  dos_model / visualization items all wait for the next chat.

**What's queued:**
- **Critical**: Quentin runs band-tracking pipeline on real LG28
  candidate. Output: confirmation or threshold adjustments.
- Atlas-side age layer patch (Option A above, ~880 LOC, ~250 LOC tests).
- LANTA-side `STEP_C01f_e` (BUSCO 4D pipeline) and `STEP_C01f_f`
  (region_popstats joiner). Both small Python scripts.
- Then: `het_location.js`, `dos_build_macro_plus_het_model`,
  band-tracking visualization page.

End of handoff.
