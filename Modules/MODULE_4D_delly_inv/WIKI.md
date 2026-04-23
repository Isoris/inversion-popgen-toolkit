# MODULE_4D — DELLY Inversion Calling — Wiki

Detail-level reference for everything specific to DELLY's INV pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (how DELLY classifies SVs in general, BISER2/minimap2, population binomial, DELLY×Manta concordance), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module is central to the inversion paper](#1-why-this-module-is-central-to-the-inversion-paper)
2. [How DELLY produces an INV record](#2-how-delly-produces-an-inv-record)
3. [Filter thresholds (verified from code)](#3-filter-thresholds-verified-from-code)
4. [Reading an INV VCF record](#4-reading-an-inv-vcf-record)
5. [What the output feeds downstream](#5-what-the-output-feeds-downstream)
6. [Limitations specific to DELLY INV at ~5× coverage](#6-limitations-specific-to-delly-inv-at-5-coverage)
7. [Manuscript-ready prose snippets](#7-manuscript-ready-prose-snippets)
8. [Q&A](#8-qa)

---

## 1. Why this module is central to the inversion paper

DELLY INV is the **primary direct SV-caller evidence** for every inversion candidate the paper claims. The framing chain is:

- MODULE_5A discovers inversion candidates from population-level signal — the 4-layer independence framework combining Layer A (local PCA / dosage), Layer C (GHSL within-sample haplotype divergence), and the integrative scoring that combines them. These are *statistical* candidates — regions of the genome where the population data behaves as if there is an inversion.
- MODULE_4D INV provides *molecular* evidence — actual split reads and discordant pairs at the inferred breakpoints, classifiable by orientation as 3to3 (left) or 5to5 (right) inversion junctions.
- MODULE_5A2 then tests whether the molecular and statistical evidence agree. Concordance is the highest-confidence class.

Without MODULE_4D, every inversion claim in the paper would rest on population-genetic inference alone. With MODULE_4D, claims have direct read-level support at the breakpoints. This is the difference between *"the population data is consistent with an inversion here"* and *"there is an inversion here, with split-read evidence at both breakpoints, in N samples"*.

A second, equally important role: DELLY INV calls with **no MODULE_5A signal** are not failures — they are candidates for **low-frequency inversions below local PCA's detection limit**. At ~5× coverage, local PCA needs roughly a population-frequency floor to discriminate karyotypes; DELLY can call rare or even singleton inversions that local PCA cannot see. These DELLY-only INV calls feed into the final candidate set as a separate class.

---

## 2. How DELLY produces an INV record

The full general explanation of DELLY's split-read classifier lives in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2C](../SV_CALLING_FRAMEWORK_WIKI.md#2c--selectinversions--inv-3to3--5to5). The INV-specific summary is below.

### The two junctions of an inversion

An inverted segment between coordinates `bp1` and `bp2` creates exactly two breakpoint junctions in the reference frame:

- **Left breakpoint (3to3 / tail-to-tail)** — the 3' end of the upstream non-inverted segment is joined to the 3' end of the (now reversed) inverted segment. Reads spanning this junction map in opposing strand orientations with both soft-clips on the **right** (3') side.
- **Right breakpoint (5to5 / head-to-head)** — the 5' end of the (now reversed) inverted segment is joined to the 5' end of the downstream non-inverted segment. Reads spanning this junction map in opposing strand orientations with both soft-clips on the **left** (5') side.

### How DELLY's `selectInversions` populates them

In `src/junction.h`:

```c
// For each pair of split-read alignments at the same locus
if (same chr) AND (different strand direction) AND (agreeing soft-clip side)
   AND (ref distance > minRefSep):
    if (i.scleft)  →  br[1].push_back(...)   // 5to5 = right breakpoint
    else           →  br[0].push_back(...)   // 3to3 = left breakpoint
```

The single boolean check `i.scleft` is what splits the two breakpoints into separate pre-classification buckets. After all reads are bucketed, `delly call -t INV` performs the **pairing step**: each `3to3` junction is matched to the nearest `5to5` junction within size limits to produce a `<INV>` VCF record.

### What happens when pairing fails

If only one of the two junctions has enough supporting reads to pass DELLY's internal thresholds, the unpaired junction is **demoted to a BND record** (with `CT=3to3` or `CT=5to5`) instead of being promoted to an INV. These demoted junctions live in MODULE_4E's BND catalog and are the substrate for MODULE_5A2 STEP06's **orphan rescue** pipeline — it pairs leftover BNDs across both DELLY and Manta to recover inversions both INV-typers missed.

### Comparison to Manta

Manta represents inversions natively as **two BND records sharing an `EVENT` tag**, with `INFO/INV3` and `INFO/INV5` flags marking the orientation. The script `convertInversion.py` (distributed with Manta) merges each `INV3+INV5` pair into a single `<INV>` VCF record before cohort merging. Pairs that fail to merge stay as BNDs in the raw pre-conversion VCF and feed into the same orphan-rescue pipeline.

---

## 3. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
PASS + INFO/PRECISE=1 + QUAL ≥ 300 + INFO/PE ≥ 3
```

Source: `00_module4d_config.sh` lines 70–71 and `slurm/SLURM_A03_merge_genotype.sh` line 123:

```bash
STRICT_EXPR='INFO/SVTYPE="INV" && INFO/PRECISE=1 && QUAL>=300 && INFO/PE>=3'
```

### Why these specific values

| Threshold | Value | Why this value |
|---|---|---|
| QUAL | ≥ 300 | DELLY's internal QUAL combines paired-end + split-read + mapping-quality scores. 300 is the documented threshold above which DELLY considers a call high-confidence. The DELLY paper uses 300 as the recommended INV cutoff. |
| PE | ≥ 3 | Three or more paired-end reads supporting the inversion. At ~5× coverage, requiring 4 or 5 would drop a substantial fraction of true heterozygous inversions where only one allele provides discordant pairs. |
| PRECISE | required | INV calls without split-read confirmation (`PRECISE=0`, i.e. `IMPRECISE`) often have very wide CIPOS (±500 bp or more), making downstream breakpoint validation unreliable. |
| FILTER = PASS | required | Excludes DELLY's own LowQual flag. |

Why **not** higher PE? Looping the threshold (PE 3 → 4 → 5) at the discovery filter would change the candidate set itself, not just evidence strength. The principled approach is: call once at the inclusive threshold (PE ≥ 3), then let MODULE_5A2 STEP02 re-extract per-sample BAM evidence and STEP03's Fisher / Cochran–Armitage tests decide which candidates are statistically supported. The seed qualification at STEP03 already requires `≥60% INV-group support`, `concordance ≥ 0.80`, `≥5 INV samples` — that is where stringency lives, not at the discovery threshold.

If a candidate's robustness to evidence stringency is in doubt, filter the STEP02 evidence table to `alt ≥ 4` or `alt ≥ 5` and re-run STEP03. Costs nothing, gives you the sensitivity curve per-candidate without polluting the discovery set.

---

## 4. Reading an INV VCF record

A real example from this pipeline:

```
Scaffold03  7600386  INV00001054  T  <INV>  209  PASS
PRECISE;SVTYPE=INV;SVMETHOD=EMBL.DELLYv0.8.3;END=9802651;PE=2;MAPQ=10;
CT=5to5;CIPOS=-10,10;CIEND=-10,10;SRMAPQ=60;INSLEN=0;HOMLEN=9;
SR=4;SRQ=0.98913;CONSENSUS=GGGTGCT...
GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV
1/1:-18.132,-2.4413,0:24:PASS:25114:131670:78351:3:0:2:0:9
```

### Field-by-field

| Field | Value | Meaning |
|---|---|---|
| CHROM | Scaffold03 | Reference contig |
| POS | 7600386 | Left breakpoint (1-based) |
| ID | INV00001054 | DELLY-assigned record ID |
| REF / ALT | T / `<INV>` | Symbolic ALT — actual content not specified |
| QUAL | 209 | DELLY's combined quality score (this one is **below the 300 strict cutoff** — would be filtered out) |
| FILTER | PASS | Passed DELLY's internal filter |
| `PRECISE` | flag | Split-read confirmed |
| `SVTYPE=INV` | symbolic | This is an inversion |
| `END=9802651` | int | Right breakpoint |
| `PE=2` | int | Paired-end support: 2 reads (**below PE≥3 strict cutoff**) |
| `MAPQ=10` | int | Mean mapping quality of supporting reads (low — repetitive region likely) |
| `CT=5to5` | string | This is the **right** breakpoint of the inversion (head-to-head junction) |
| `CIPOS=-10,10` | range | 95% confidence interval on POS (±10 bp — tight because PRECISE) |
| `CIEND=-10,10` | range | 95% confidence interval on END (±10 bp) |
| `SRMAPQ=60` | int | Mapping quality of split reads specifically (high — split reads themselves map well) |
| `INSLEN=0` | int | No inserted sequence at the junction |
| `HOMLEN=9` | int | 9 bp of microhomology at the junction (bridges blunt end + small overlap) |
| `SR=4` | int | Split-read support: 4 reads spanning the junction |
| `SRQ=0.98913` | float | Split-read consensus alignment quality (0–1; ≥0.8 needed for PRECISE) |
| `CONSENSUS=...` | string | Assembled breakpoint sequence — can be BLAT'd to confirm orientation |

### Per-sample FORMAT fields

For each sample (one column per BAM), DELLY reports:

| Format key | Value in example | Meaning |
|---|---|---|
| `GT` | `1/1` | Genotype: homozygous inverted |
| `GL` | `-18.132,-2.4413,0` | Genotype likelihoods (log10) for `0/0`, `0/1`, `1/1` |
| `GQ` | `24` | Genotype quality |
| `FT` | `PASS` | Per-sample filter |
| `RCL` | `25114` | Read count, left flanking (reference signal) |
| `RC` | `131670` | Read count, in the inversion (informative for copy-number, not orientation) |
| `RCR` | `78351` | Read count, right flanking |
| `CN` | `3` | Copy-number estimate (DELLY's call — not always correct for INV which is copy-neutral) |
| `DR` | `0` | Discordant **R**eference paired-end reads |
| `DV` | `2` | Discordant **V**ariant paired-end reads (alt-supporting) |
| `RR` | `0` | Junction **R**eference reads |
| `RV` | `9` | Junction **V**ariant reads (alt-supporting) |

**Per-sample alt evidence** = `DV + RV` = `2 + 9` = **11 reads**. This is what MODULE_5A2 STEP02 re-extracts per sample, and what enters the binomial population-confidence annotation.

### What this record actually tells you

- An inversion is supported between Scaffold03:7,600,386 and 9,802,651 (~2.2 Mb)
- Both breakpoints have ±10 bp confidence (tight)
- The right breakpoint shows the 5to5 head-to-head junction
- The CT=5to5 record is one half — the matching `CT=3to3` record at the left breakpoint should appear as a separate VCF row, paired into this INV by `delly call -t INV`
- The single sample shown is homozygous for the inversion with strong evidence (11 alt reads)
- The 9 bp microhomology at the junction suggests an MMEJ-mediated formation mechanism (not NAHR which would show flanking duplications instead)
- **This particular record fails strict filter** (QUAL=209 < 300, PE=2 < 3) — it would not enter the strict catalog. It might still enter the relaxed catalog and be revisited via MODULE_5A2 BAM re-evidence.

---

## 5. What the output feeds downstream

The INV catalog is the most heavily consumed output in the entire MODULE_4 family.

```
MODULE_4D INV catalog
  │
  ├─▶ MODULE_5A2 STEP01  candidate extraction
  │   - QUAL ≥ 300, PE ≥ 3 (matches strict filter)
  │   - Size window 5 kb ≤ SVLEN ≤ 20 Mb
  │   - Output: unified candidate table with caller='delly'
  │
  ├─▶ MODULE_5A2 STEP02  per-sample BAM re-evidence
  │   - For each candidate, pysam scan ±300 bp at each breakpoint
  │   - Five evidence classes: discordant, split, mate-link, soft-clip, shared
  │   - Composite score = n_disc + 2×n_split + n_mate + 3×n_shared
  │
  ├─▶ MODULE_5A2 STEP03  Fisher + Cochran–Armitage seed qualification
  │   - Three tests per candidate (Fisher INV-vs-non, Fisher INV-vs-REF, CA trend)
  │   - Bonferroni P < 0.05, ≥60% INV support, ≥0.80 concordance, ≥5 INV samples
  │   - Qualified candidates → SEEDS for MODULE_5A C01i decomposition
  │
  ├─▶ MODULE_5A2 STEP05  DELLY×Manta concordance
  │   - 50% reciprocal overlap test against MODULE_4G Manta INV catalog
  │   - Six classes: concordant / delly_only / manta_only / bnd_rescue_* / orphan
  │
  ├─▶ MODULE_5A   regime / candidate matching
  │   - 50% reciprocal overlap against the multi-layer consensus candidate regions (Layers A+C)
  │   - Tags whether each population candidate has SV-caller anchor
  │
  └─▶ Manuscript Figure (per-candidate breakpoint validation panel)
      - DELLY split-read pile-up at both breakpoints
      - INV3/INV5 orientation track
```

For the specific cross-module wiring (cheats, framework methods, plug-map for the whole family), see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map).

---

## 6. Limitations specific to DELLY INV at ~5× coverage

1. **PE-supported-only INV calls are common.** At low coverage, getting split reads across both breakpoints is hard. DELLY frequently produces `IMPRECISE` INV calls where only PE evidence supports both junctions. The strict PRECISE requirement removes ~30–50% of such calls. Some are real inversions whose breakpoints we will miss until long-read data is available.

2. **Repeat-rich breakpoints get downgraded.** When one or both breakpoints fall in a repeat element (TE family, segmental duplication), DELLY's split-read aligner often fails to assemble a clean junction. These calls become BND records and require MODULE_4E's orphan-rescue pathway to recover. Inversions formed by NAHR — which by definition have flanking repeats — are particularly affected.

3. **Large inversions (>5 Mb) may not pair correctly.** DELLY's INV pairing has internal size limits. Very large inversions can produce two separate BNDs that DELLY does not recognize as a pair. Cross-validation with Manta and BND rescue handles this in MODULE_5A2 STEP06, but the candidate may be missed by DELLY's own INV catalog entirely.

4. **Heterozygous inversions in single samples.** A heterozygous inversion produces only the alt-supporting reads from one of the two haplotypes — at ~5× total coverage, that means ~2.5× per haplotype, meaning often only 1–2 reads supporting the alt. PE ≥ 3 already excludes some real heterozygous singletons. The population-binomial annotation partly compensates by looking for accumulation across samples ([`../SV_CALLING_FRAMEWORK_WIKI.md` § population binomial](../SV_CALLING_FRAMEWORK_WIKI.md#5-filter-attrition-through-the-pipeline-upsetr)).

5. **`delly filter -f germline` is a black box.** Documented as applying allele-frequency, Hardy-Weinberg, and genotype-quality filters, but the exact internal criteria are not publicly specified per filter step. Treated as a trusted upstream tool; for any single candidate that disappears at this stage, manual re-inspection of the pre-filter BCF is the only audit path.

6. **CN=2 always reported regardless of inversion orientation.** Inversions are copy-neutral, so DELLY's `FORMAT/CN` field is uninformative for INV. Don't use it for genotyping. Use `GT` (which DELLY infers from `DV/DR` and `RV/RR`) instead, or recompute genotype from per-sample alt evidence in MODULE_5A2 STEP02.

---

## 7. Manuscript-ready prose snippets

### Methods — DELLY INV calling paragraph

> DELLY2 v1.7.3 (Rausch *et al.*, 2012) was used for inversion discovery across all 226 *C. gariepinus* samples. Per-sample `delly call -t INV` was run on duplicate-marked BAMs against the *C. gariepinus* (Gar) subgenome reference (`fClaHyb_Gar_LG.fa`, 28 pseudochromosomes, ~964 Mb — the *gariepinus* haplotype extracted from the haplotype-resolved F₁ hybrid assembly) using the empirical exclusion BED constructed from regions with <500 callable base pairs per 50-kb bin and unconditional 50-kb chromosome-end masking. Per-sample calls were merged (`delly merge`), regenotyped against all samples, subset to the 81 NAToRA-pruned unrelated individuals, and germline-filtered (`delly filter -f germline`). The strict catalog required `FILTER = PASS`, `INFO/PRECISE = 1`, `QUAL ≥ 300`, and `INFO/PE ≥ 3`. Per-sample alt evidence was computed as the sum of discordant variant pairs and variant junction reads (`FORMAT/DV + FORMAT/RV`).

### Results — number reporting placeholder

> DELLY identified [TODO N_RAW] candidate inversions across the 226 samples in the cohort merge. After germline filtering and the strict per-record cutoff (PASS + PRECISE + QUAL ≥ 300 + PE ≥ 3), [TODO N_STRICT] inversions remained ([TODO PCT]% of the raw merged set; per-type attrition profile in Supplementary Figure SX). [TODO N_CONCORDANT] of these were concordant with Manta INV calls within a 50% reciprocal overlap window (MODULE_5A2 STEP05); the remaining [TODO N_DELLYONLY] were DELLY-specific calls. The orphan-BND rescue pipeline (MODULE_5A2 STEP06) recovered an additional [TODO N_RESCUE] inversion candidates from paired 3to3+5to5 breakpoints that did not survive either caller's INV-typing step.

### Discussion — limitations footnote

> Sensitivity for heterozygous inversions at ~5× coverage is bounded by per-sample read depth. Each heterozygote contributes alt-supporting reads from only one chromatid, so even at the inclusive PE ≥ 3 threshold, singleton heterozygous inversions are systematically under-detected. The population-binomial annotation (see Methods § Confidence annotation) partly mitigates this by accumulating evidence across carriers, but inversions present in a single individual remain at the limit of detection.

---

## 8. Q&A

Open-ended questions that came up during pipeline development. Add new ones as they arise; the repo-level Q&A file (`Modules/QA.md`, planned) will aggregate across modules.

### Q: Why PE ≥ 3 and not 4 or 5?

A short answer: at ~5× total coverage, raising the threshold drops a meaningful fraction of true heterozygous singletons. The detailed reasoning is in [§ 3](#3-filter-thresholds-verified-from-code). The principled separation between "discovery threshold" and "evidence stringency" lives in MODULE_5A2 STEP02/STEP03 — re-filter the per-sample evidence table there rather than at the discovery filter.

### Q: Should I loop the discovery threshold over PE 3, 4, 5 to test robustness?

No. Looping the discovery threshold gives you three different candidate sets to reconcile, which is a downstream nightmare. Call once at the inclusive threshold and re-stratify the per-sample evidence in STEP02. See [§ 3](#3-filter-thresholds-verified-from-code).

### Q: What is the difference between `CT=3to3` and `CT=5to5`?

They are the two breakpoints of the same inversion. `3to3` is the left breakpoint (tail-to-tail junction), `5to5` is the right breakpoint (head-to-head junction). A complete INV record has both; `delly call -t INV` pairs them. See [§ 2](#2-how-delly-produces-an-inv-record) for the geometry and [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2C](../SV_CALLING_FRAMEWORK_WIKI.md#2c--selectinversions--inv-3to3--5to5) for the full junction.h dissection.

### Q: An inversion is in my MODULE_5A candidate list but not in MODULE_4D. Should I trust it?

Possibly yes — see [§ 1](#1-why-this-module-is-central-to-the-inversion-paper) for why DELLY-INV-negative population candidates can still be real (low-frequency inversions below local PCA's effective limit go in the opposite direction). The right diagnostic is MODULE_4E (BND orphan rescue). If MODULE_4E shows paired CT=3to3 + CT=5to5 BNDs at the candidate coordinates, the inversion is real but DELLY's INV-pairing step missed it. If neither MODULE_4D nor MODULE_4E shows anything, the population signal is more likely a recombination-suppression artifact (selection sweep, balanced rearrangement) than a true inversion.

### Q: An inversion is in my MODULE_4D catalog but not in MODULE_5A. Real or artifact?

Often real but at low population frequency. MODULE_5A's population-signal-based detection requires enough karyotype carriers to discriminate a band in the local PCA. Inversions at frequency <5–10% in 226 samples may not produce a discriminable population pattern. The DELLY-only INV calls feed into the final candidate set as a separate "rare/private inversion" class. To distinguish from artifacts, check (i) whether the call is concordant with Manta (MODULE_5A2 STEP05), (ii) whether MODULE_5A2 STEP02 finds per-sample BAM evidence in the apparent carriers, and (iii) whether the breakpoints fall in repeat-rich or assembly-error-prone regions.

### Q: What does `IMPRECISE` mean — is the inversion fake?

No. `IMPRECISE` means DELLY found discordant paired-end evidence but could not assemble a split-read consensus across the junction. The inversion likely exists; only the **breakpoint coordinates** are uncertain (look at the wide CIPOS/CIEND values). At ~5× coverage many real inversions are IMPRECISE because there isn't enough depth to assemble both junctions. The strict catalog requires PRECISE because downstream BAM-evidence extraction needs reliable coordinates; IMPRECISE INV records are kept in the relaxed catalog and revisited only when MODULE_5A flags the population-signal region.

### Q: DELLY says CN=2 for my homozygous inversion but it should be CN=2 anyway. Is the genotype correct?

Yes — inversions are copy-neutral by definition (the same DNA is present, just flipped), so CN=2 is correct for a diploid sample regardless of inversion zygosity. The genotype itself comes from `GT`, which DELLY infers from the `DV/DR` and `RV/RR` ratios, not from copy-number depth. For INV calls specifically, ignore `FORMAT/CN` — it tells you nothing about the inversion.

### Q: How is HOMLEN at the breakpoint useful?

`HOMLEN=N` reports N bp of microhomology between the two sides of the junction. Length distribution tells you about the formation mechanism: 0 bp (blunt) → NHEJ; 1–7 bp → MMEJ (the dominant mechanism for de novo inversions, peak at 2–3 bp per Sultana 2017); >10 bp → likely template-mediated (FoSTeS/MMBIR or NAHR). For MODULE_5A's mechanism inference (Cheat 14, Cheat 27), HOMLEN is one of the key features distinguishing NAHR-mediated recurrent inversions from single-origin events.
