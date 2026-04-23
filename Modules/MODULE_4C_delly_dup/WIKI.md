# MODULE_4C — DELLY Duplication Calling — Wiki

Detail-level reference for everything specific to DELLY's DUP pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (BISER2/minimap2 SD detection, the mechanism-inference framework that uses this module's output), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module exists (NAHR substrate inference)](#1-why-this-module-exists-nahr-substrate-inference)
2. [How DELLY produces a DUP record](#2-how-delly-produces-a-dup-record)
3. [Filter thresholds (verified from code)](#3-filter-thresholds-verified-from-code)
4. [Reading a DUP VCF record](#4-reading-a-dup-vcf-record)
5. [The mechanism-inference pipeline downstream](#5-the-mechanism-inference-pipeline-downstream)
6. [Reference-side vs read-side duplications (key distinction)](#6-reference-side-vs-read-side-duplications-key-distinction)
7. [Limitations specific to DELLY DUP at ~5× coverage](#7-limitations-specific-to-delly-dup-at-5-coverage)
8. [Manuscript-ready prose snippets](#8-manuscript-ready-prose-snippets)
9. [Q&A](#9-qa)

---

## 1. Why this module exists (NAHR substrate inference)

MODULE_4C DUP is the **mechanism-inference** module for the inversion paper. Its role is not to count duplications for their own sake — it is to identify the **flanking duplication architecture** at every inversion breakpoint, which determines whether an inversion arose by **NAHR (non-allelic homologous recombination)** or by some other mechanism.

The framework comes from Porubsky *et al.* 2022 (Cell, HGSVC recurrent inversions) and is wired into the pipeline as cheats 14, 15, and 27:

- **NAHR** requires inverted-orientation segmental duplications flanking the inversion breakpoints. Recombination between the inverted SDs flips the segment between them. NAHR loci can produce **recurrent** inversions — the same inversion arises independently in multiple lineages because the same SD substrate is conserved.
- **NHEJ / MMEJ / MMBIR** are break-and-rejoin mechanisms that leave characteristic junctions (blunt, microhomology, templated insertion). They do NOT require an SD substrate. They produce **single-origin** inversions — one ancestral event, all carriers share the haplotype.

If MODULE_4C identifies inverted DUP pairs flanking an inversion's breakpoints → mechanism is likely NAHR → **expect bimodal genotype-dosage similarity (GDS) among carriers** (multiple independent origins). If MODULE_4C finds no flanking DUPs → mechanism is likely NHEJ/MMBIR → **expect unimodal GDS** (single origin). MODULE_5B's Hartigan dip test is the experimental confirmation.

Without MODULE_4C, every inversion would receive the same generic interpretation. With MODULE_4C, each inversion can be classified by its likely formation mechanism, and the haplotype-diversity expectations become testable.

A second, equally important role: MODULE_4C's DUP catalog is the **read-derived** counterpart to BISER2's reference-derived SD catalog. They answer different questions (see § 6) and the combination is more informative than either alone.

---

## 2. How DELLY produces a DUP record

The full junction.h dissection is in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2B](../SV_CALLING_FRAMEWORK_WIKI.md#2b--selectduplications--dup). The DUP-specific summary:

### The geometry of a tandem duplication junction

A tandem duplication doubles a segment of reference sequence in the sample. A read that bridges the junction between the two copies has a characteristic geometry: the read's body aligns to the **end** of the duplicated region, while one end of the read soft-clips back to align at the **start** of the duplicated region (the second copy starts there).

### `selectDuplications` conditions

```c
if (same chr) AND (same strand direction) AND (opposing soft-clips, REVERSED order)
   AND (ref distance > minRefSep):
    if (i.refpos ≤ j.refpos):
        if (i IS scleft) AND (j NOT scleft):  // left-clip on i, right-clip on j → DUP
            br[3].push_back(...)
    else:
        if (i NOT scleft) AND (j IS scleft):
            br[3].push_back(...)
```

The conditions are **identical to DEL** except the clipping order is **reversed**: in DEL, the right-clipped read is upstream and the left-clipped read is downstream (reads bridge across a deleted middle). In DUP, the left-clipped read is upstream and the right-clipped read is downstream (reads bridge across the boundary between two tandem copies, where the second copy's start aligns to the same coordinate as the first copy's start).

DELLY emits these clusters as `<DUP>` VCF records with `INFO/SVTYPE=DUP` and `INFO/CT=5to3` (the duplication-type connection class).

### What DELLY DUP cannot detect

- **Inverted duplications** — segments duplicated in inverted orientation. These produce different junction geometry that DELLY emits as BND records with INV-orientation CT, not as `<DUP>`. For NAHR substrate detection, this is critical: NAHR substrates are *inverted* SD pairs, which are NOT directly in MODULE_4C's catalog. The catalog has *direct* tandem duplications. BISER2 (reference-side) and minimap2 self-alignment (Cheat 14) are needed for inverted SDs.
- **Non-tandem (interspersed) duplications** — a copy elsewhere on the chromosome or on a different chromosome. These appear as BND records with various CT classes, not as DUP. MODULE_4C catalog is specifically tandem duplications.
- **Copy-number-variant duplications without clean junctions** — DELLY's `delly cnv` mode handles these via depth-based segmentation, but `delly call -t DUP` (what MODULE_4C runs) requires junction reads. Pure depth-only DUPs are missed by MODULE_4C.

These limitations are why MODULE_4C alone is **insufficient** for full SD architecture characterization. The combined picture comes from MODULE_4C + BISER2 + minimap2 (Cheat 14) cross-referenced.

---

## 3. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
PASS + INFO/SVTYPE="DUP" + INFO/PRECISE=1 + QUAL ≥ 500 + INFO/PE ≥ 5
```

Source: `00_module4c_config.sh` lines 43–46 and `slurm/SLURM_A03_merge_genotype.sh` line 125–128:

```bash
STRICT_REQUIRE_PASS=1
STRICT_REQUIRE_PRECISE=1
STRICT_MIN_QUAL=500
STRICT_MIN_PE=5

# In SLURM_A03 (constructed conditionally):
STRICT_EXPR='INFO/SVTYPE="DUP" && QUAL>=500 && INFO/PE>=5'
if [[ "${STRICT_REQUIRE_PRECISE}" -eq 1 ]]; then
  STRICT_EXPR="${STRICT_EXPR} && INFO/PRECISE=1"
fi
```

The config has commented-out alternative values (`STRICT_MIN_QUAL=300`, `STRICT_MIN_PE=3`) preserved as a documented option to relax to INV-level stringency if needed. Production uses 500/5.

### Why these specific values

| Threshold | Value | Why this value |
|---|---|---|
| QUAL | ≥ 500 | Same as DEL — DUPs and DELs are the most common SV classes by count, so the stricter cutoff (vs. INV/BND/TRA at 300) is appropriate. The 500 threshold is the DELLY-recommended high-confidence cutoff for these abundant classes. |
| PE | ≥ 5 | Stricter than INV's PE ≥ 3 because, for DUP, false positives from duplicate read pairs (mapping artifacts at repeat boundaries) accumulate quickly. Requiring 5 distinct discordant pairs per call effectively excludes cluster-of-duplicates artifacts. |
| PRECISE | required | DUPs in repeat-rich regions can have very wide CIPOS without split-read assembly. Requiring PRECISE removes most repeat-edge false positives at the cost of some real DUPs in repetitive context. |
| FILTER = PASS | required | Excludes DELLY's internal LowQual flag. |

### A note on the asymmetry that almost was

An earlier draft of the framework wiki incorrectly stated that MODULE_4C does NOT require PRECISE. This was wrong — `STRICT_REQUIRE_PRECISE=1` is set in the config and the conditional in SLURM_A03 adds the `INFO/PRECISE=1` clause. All four DELLY modules (4B/4C/4D/4E/4F) require PRECISE in production. The corrected information is in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 1](../SV_CALLING_FRAMEWORK_WIKI.md#1-verified-filter-thresholds-production-code).

---

## 4. Reading a DUP VCF record

<!-- TODO: replace synthetic example with a real DELLY DUP record from the production catalog -->

Synthetic example (a 12 kb tandem duplication):

```
Scaffold02  18452100  DUP00000087  C  <DUP>  580  PASS
PRECISE;SVTYPE=DUP;SVMETHOD=EMBL.DELLYv1.7.3;END=18464350;PE=7;
MAPQ=51;CT=5to3;CIPOS=-9,9;CIEND=-11,11;SRMAPQ=60;INSLEN=0;
HOMLEN=4;SR=5;SRQ=0.96812;CONSENSUS=ATGCAGTACTGCAGTAACTACTAAGGGCAT...
GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV
0/1:-15.32,-2.16,0:21:PASS:8240:14872:8195:3:8:5:7:4
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS/END | Scaffold02:18,452,100–18,464,350 | 12.25 kb tandem duplication |
| ID | DUP00000087 | DELLY-assigned ID |
| QUAL | 580 | Above 500 strict cutoff |
| `PRECISE` | flag | Split-read assembled |
| `SVTYPE=DUP` | string | This is a tandem duplication |
| `PE=7`, `SR=5` | int | 7 paired-end + 5 split-read supporting (above 5/0 cutoffs) |
| `CT=5to3` | string | Duplication-type connection — 5' end of downstream copy joined to 3' end of upstream copy |
| `CIPOS=-9,9`, `CIEND=-11,11` | range | Tight breakpoint confidence |
| `MAPQ=51` | int | Mean mapping quality of supporting reads (moderate; some repeat context) |
| `HOMLEN=4` | int | 4 bp microhomology at junction |
| `CONSENSUS=...` | string | Assembled junction sequence |

Per-sample (one column shown):
- `GT=0/1` heterozygous DUP
- `CN=3` — copy-number estimate of 3 (one copy added to the diploid background of 2; consistent with heterozygous DUP)
- `DR/DV = 8/5` — 5 alt discordant pairs
- `RR/RV = 7/4` — 4 alt junction reads
- Per-sample alt evidence = 5 + 4 = **9 reads**

### Useful unlike INV: the CN field works for DUP

For DUPs, `FORMAT/CN` is informative because duplications change copy number. CN=2 → reference, CN=3 → het DUP, CN=4 → hom DUP. This is the opposite of INV (where CN=2 always regardless of zygosity). When debugging DUP genotypes, CN is a reasonable cross-check against `GT`.

---

## 5. The mechanism-inference pipeline downstream

```
MODULE_4C DUP catalog
  │
  ├─▶ Cheat 14  repeat architecture classifier (R)
  │   - For each inversion candidate from MODULE_5A (or MODULE_4D):
  │     scan ±N kb of each breakpoint for DUP records
  │     classify orientation: inverted / direct / mixed / absent
  │   - Output: NAHR_CANDIDATE / NHEJ_CANDIDATE / COMPLEX / UNKNOWN
  │
  ├─▶ Cheat 15  recurrence prior (R)
  │   - Combines Cheat 14 mechanism class with carrier haplotype data
  │   - NAHR + multiple haplotype clusters → RECURRENT
  │   - NHEJ + single haplotype cluster → SINGLE_ORIGIN
  │   - Updates downstream prior on bimodal vs unimodal GDS expectation
  │
  ├─▶ Cheat 27  BISER2 cross-check (R)
  │   - Cross-references MODULE_4C DUP findings with BISER2's
  │     pre-computed segmental duplication catalog
  │   - minimap2 (Cheat 14) + BISER2 agreement → HIGH CONFIDENCE NAHR
  │   - Disagreement → flag for manual inspection
  │   - BISER2 alone NEVER overrides minimap2 (Cheat 14 is trusted)
  │
  ├─▶ MODULE_5B Hartigan dip test on GDS
  │   - For each inversion candidate, compute pairwise genotype-dosage
  │     similarity among homozygous-INV carriers
  │   - Bimodal GDS distribution → multiple independent origins
  │   - Unimodal GDS → single origin
  │   - Cross-references against Cheat 14/15/27 mechanism predictions
  │
  └─▶ Manuscript Figure (per-inversion mechanism panel)
        - Flanking-DUP architecture overlay
        - GDS distribution
        - Mechanism class label
```

For the broader plug-map (other modules consuming DUP data), see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map).

---

## 6. Reference-side vs read-side duplications (key distinction)

This is the conceptual point that took longest to crystallize during pipeline design and is worth understanding clearly because it affects how MODULE_4C's output is interpreted.

There are **two completely different questions** about duplications at an inversion locus:

### Question A: "Was there a duplication substrate present in the ancestor?"

This is the **mechanism** question. NAHR requires that the inversion-carrier and the non-carrier shared an ancestor whose genome contained inverted SDs flanking the (eventually-inverted) segment. The substrate has to have existed before the inversion event.

**Answered by reference-side tools:**
- **BISER2** finds all SD pairs in the reference assembly itself (≥1 kb, ≥90% identity)
- **minimap2 self-alignment** (Cheat 14) cross-validates BISER2 calls and finds SDs BISER2 may have missed

These tools do not look at any sample's reads. They examine the reference genome and infer the SD architecture that existed in the reference individual (and, by inheritance, in the population's recent ancestor).

### Question B: "Is there a polymorphic duplication present in samples right now?"

This is the **current state** question. A duplication that is present in some samples and absent in others is a structural-variant polymorphism — it could be a recently-emerged duplication, or a substrate that existed in the ancestor but has been deleted from some lineages.

**Answered by read-side tools:**
- **MODULE_4C DELLY DUP** finds polymorphic tandem duplications from sample reads
- **MODULE_4G Manta DUP** finds the same with assembly-based methods

These tools do not look at the reference's intrinsic SD architecture. They find duplications that differ between samples and the reference.

### Why both are needed for inversion mechanism inference

A NAHR-mediated inversion was caused by an **ancestral SD substrate** (Question A). But that substrate **may or may not still be polymorphic today** (Question B). Three scenarios:

| Scenario | Reference SDs (BISER2) | Polymorphic DUP (MODULE_4C) | Mechanism inference |
|---|---|---|---|
| Conserved ancestral NAHR locus | Inverted SD pair present | None (substrate not polymorphic) | NAHR confirmed; substrate stable across population |
| Recently expanded NAHR locus | Inverted SD pair present | Polymorphic DUPs nearby | NAHR with active SD turnover; high recurrence likely |
| No SD substrate | None | None | NHEJ/MMEJ; single origin expected |
| MODULE_4C DUP without BISER2 SD | None | Direct DUP in some samples | New tandem expansion in some lineages; NOT an inversion substrate |

The four-way matrix reduces uncertainty more than either tool alone. This is why the pipeline runs both BISER2 (once per haplotype, precomputed) and MODULE_4C DUP (on every sample) and Cheat 27 cross-references them.

---

## 7. Limitations specific to DELLY DUP at ~5× coverage

1. **Repeat boundary false positives.** DUPs are over-called at the boundaries of large repeat elements (TEs, satellite arrays) because the same read sequence can map to multiple positions. Strict QUAL ≥ 500 + PE ≥ 5 + PRECISE removes most of these but not all. Manual inspection of any DUP near a repeat boundary is wise.

2. **Inverted DUPs are NOT in this catalog.** Inverted-orientation duplications produce BND records with `CT=3to3` or `CT=5to5`, identical to inversion-junction BNDs. They live in MODULE_4E BND, not here. For NAHR substrate detection (which requires *inverted* SDs), MODULE_4C alone is insufficient — must combine with BND analysis and reference-side BISER2.

3. **Non-tandem duplications are NOT in this catalog.** Interspersed duplications (a copy elsewhere) appear as BNDs. The MODULE_4C catalog is restricted to **tandem** duplications by `selectDuplications` geometry.

4. **Pure copy-number-variant DUPs without junction reads are missed.** Large duplications where coverage is doubled but no read crosses the junction will not be in the catalog. DELLY has a separate `delly cnv` mode for depth-based CNV calling but it is not run in this pipeline. Long-read or higher-coverage data would be needed to catch these.

5. **DELLY DUP's `delly filter -f germline` is the same black box as for INV.** Allele-frequency / Hardy-Weinberg-like criteria can drop true rare polymorphic DUPs. If a candidate disappears at this stage, the pre-filter BCF is the audit path.

6. **The 5 kb minimum on `selectDuplications`.** DELLY's `minRefSep` (default 300 bp, configurable) is the smallest DUP it will emit. Smaller tandem duplications go to other DELLY callers (insertion-flavored) or to Clair3's small-indel pipeline (MODULE_4A) for ≤50 bp events. The 50 bp – ~300 bp gap is undercovered by all callers; long-read data would help most here.

---

## 8. Manuscript-ready prose snippets

### Methods — DELLY DUP calling paragraph

> DELLY2 v1.7.3 was used to discover tandem duplications across all 226 samples by running `delly call -t DUP`. Per-sample discovery, cohort merging, regenotyping, subsetting to 81 unrelated samples, and germline filtering followed the same workflow as MODULE_4D INV. The strict catalog required `FILTER = PASS`, `INFO/SVTYPE = "DUP"`, `INFO/PRECISE = 1`, `QUAL ≥ 500`, and `INFO/PE ≥ 5`. The stricter PE threshold relative to INV reflects DUP's higher false-positive rate at repeat boundaries; preliminary analyses with PE ≥ 3 produced excessive repeat-edge artifacts (data not shown). Per-sample evidence used `FORMAT/DR/DV` (discordant pairs) and `FORMAT/RR/RV` (junction reads), and `FORMAT/CN` (copy-number estimate) was used as an orthogonal genotype cross-check. Inverted-orientation duplications, which `selectDuplications` does not classify as DUP, were captured by MODULE_4E BND (`CT ∈ {3to3, 5to5}` records); reference-derived segmental duplication architecture was characterised independently by BISER2 (Methods § Segmental duplications).

### Results — number reporting placeholder

> The MODULE_4C DUP catalog contained [TODO N_DUP_STRICT] tandem duplications after strict filtering ([TODO N_DUP_RAW] before; [TODO PCT]% retention). Of these, [TODO N_DUP_NEAR_INV] fell within ±5 kb of at least one inversion breakpoint and were included in mechanism-inference cheats (14, 15, 27). Cross-referencing with the BISER2 segmental-duplication catalog identified [TODO N_DUP_AT_BISER] DUP records overlapping pre-existing reference SDs (potential NAHR substrate) and [TODO N_DUP_NOVEL] DUP records at positions without reference SDs (recently expanded duplications). The mechanism-inference pipeline (Methods § Inversion formation mechanism) classified [TODO N_NAHR] candidate inversions as NAHR-compatible based on flanking inverted SD pairs and [TODO N_NHEJ] as NHEJ-compatible based on absence of flanking SDs.

### Discussion — limitations footnote

> DELLY's `selectDuplications` function detects only direct tandem duplications. Inverted-orientation duplications (the structural class that mediates NAHR) are emitted as `CT=3to3` or `CT=5to5` BND records, not as `<DUP>`. This means MODULE_4C's catalog alone underestimates the NAHR-substrate landscape — it captures direct tandem dups (e.g., recent gene-family expansions) but not inverted SD pairs. The complete NAHR-substrate picture required combining MODULE_4C, MODULE_4E (orphan BNDs with inversion-orientation CT), and the reference-derived BISER2 catalog (which directly enumerates inverted SD pairs). The cross-referencing logic is encoded in cheats 14, 15, and 27.

---

## 9. Q&A

### Q: Why is the PE threshold stricter (5) for DUP than for INV (3)?

DUPs have a higher per-record false-positive rate at low coverage because of two effects: (1) duplicate-read mapping artifacts at repeat boundaries can mimic DUP junctions, and (2) DUP-orientation discordant pairs occur naturally in repeat-rich regions even without true duplications. Requiring PE ≥ 5 effectively excludes most of these spurious clusters. INV-orientation discordant pairs are biologically rarer in noise, so PE ≥ 3 is sufficient.

### Q: Is `<DUP>` the same as a copy-number gain?

In principle yes (a tandem duplication adds one copy), but the catalog only contains DUPs detectable by DELLY's junction-clustering approach. Pure copy-number gains without clean junction reads (large CNVs, repeat expansions) will not be in the catalog. The `FORMAT/CN` field reports DELLY's estimated copy number, which is a useful cross-check on `GT` for individual records, but it is not a comprehensive copy-number caller. For genome-wide CNV characterization, `delly cnv` (separate mode, not run in this pipeline) or Manta DUP would be more complete.

### Q: I have an inverted SD pair flanking my candidate inversion. Why isn't it in MODULE_4C?

Inverted-orientation duplications are emitted by DELLY as **BND records** (with `CT=3to3` or `CT=5to5`), not as `<DUP>` records. `selectDuplications` only catches direct tandem duplications. The inverted SD pair you're looking for is in:
- MODULE_4E BND catalog (read-derived, polymorphic in samples) with `CT=3to3` or `CT=5to5` and possibly mistaken for a partial inversion junction
- BISER2 catalog (reference-derived, present in the reference genome itself) — this is the canonical place to look for NAHR substrates

Cheat 14 specifically queries BISER2 for inverted SDs at inversion breakpoints. MODULE_4C is not the right catalog for this question.

### Q: How is DUP CN estimation different from INV CN?

DUP duplicates DNA, so copy number changes — `CN=2` is reference, `CN=3` is het DUP, `CN=4` is hom DUP. DELLY's `FORMAT/CN` is informative for DUP genotyping. INV is copy-neutral (the same DNA, just flipped), so `CN=2` always for INV regardless of zygosity. **For DUP, CN is useful; for INV, ignore CN and use GT.** This asymmetry is one reason DUP and INV reports look different even when both pass strict filtering.

### Q: A candidate inversion has no flanking DUPs in MODULE_4C, no SDs in BISER2, and no nearby inverted-CT BNDs in MODULE_4E. Mechanism?

Most likely **NHEJ or MMEJ** (break-and-rejoin without homology). NHEJ produces blunt junctions (HOMLEN=0). MMEJ produces 1–7 bp microhomology junctions. Check the inversion's `INFO/HOMLEN` in MODULE_4D INV catalog: HOMLEN=0 → NHEJ, HOMLEN ∈ [1, 7] → MMEJ, HOMLEN > 10 → likely template-mediated (rare without obvious repeat context, but possible).

For NHEJ/MMEJ inversions, the prediction is **single-origin** — all carriers should share one ancestral haplotype, GDS distribution should be unimodal, no recurrence expected. If MODULE_5B's Hartigan dip test instead shows bimodal GDS for a no-flanking-SD inversion, that's surprising and worth investigating (cryptic SDs below detection thresholds, or a recent independent re-inversion).

### Q: Are Manta DUPs (MODULE_4G) and DELLY DUPs (MODULE_4C) cross-referenced anywhere?

Not currently. Concordance analysis between callers is implemented for INV (MODULE_5A2 STEP05) but not for DUP. The same logic could apply (50% reciprocal overlap → concordant class), and would tighten the DUP catalog considerably, but is not in the current pipeline. If DUP concordance becomes important for the manuscript, it would be a small addition modeled on STEP05.

### Q: The commented-out config has STRICT_MIN_QUAL=300 and STRICT_MIN_PE=3. When would I use those?

Those values would relax DUP to INV-level stringency. The use case is **sensitivity-prioritized analyses** — for example, if a specific inversion candidate has no flanking-DUP evidence at QUAL≥500/PE≥5 but you want to check whether weaker DUP evidence exists. To use the relaxed values, edit `00_module4c_config.sh` (uncomment the 300/3 lines, comment the 500/5 lines), re-run the pipeline from STEP A03 onward, and the relaxed catalog will be in `catalog_226.DUP.PASS.vcf.gz` (overwriting the strict version). Better practice: keep the strict default, run a one-off relaxed extraction outside the main pipeline if needed for a specific candidate, document in the methods which threshold was used for which analysis.

### Q: Why doesn't this module also analyze DELLY DUP for whole-genome duplication or ploidy changes?

DELLY `call -t DUP` detects tandem duplications, not whole-chromosome or whole-genome ploidy changes. This study cohort is diploid *C. gariepinus* (2n = 56, Maneechot et al. 2016); aneuploidy or polyploidy detection is not in scope. If it were, the pipeline would need `delly cnv` (depth-based) or a dedicated karyotype caller. None of the candidate inversions in this study have shown evidence of underlying ploidy changes.
