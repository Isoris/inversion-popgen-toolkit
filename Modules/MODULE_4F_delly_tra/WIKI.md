# MODULE_4F — DELLY Translocation (TRA) Calling — Wiki

Detail-level reference for everything specific to DELLY's TRA pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (junction.h `selectTranslocations` logic, the four CT classes for inter-chromosomal events), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module exists (false-positive control + complex events)](#1-why-this-module-exists-false-positive-control--complex-events)
2. [How DELLY identifies TRA records (BND filtered to CHR2 ≠ CHROM)](#2-how-delly-identifies-tra-records-bnd-filtered-to-chr2--chrom)
3. [The four inter-chromosomal CT classes](#3-the-four-inter-chromosomal-ct-classes)
4. [Filter thresholds (verified from code)](#4-filter-thresholds-verified-from-code)
5. [Reading a TRA VCF record](#5-reading-a-tra-vcf-record)
6. [What the output feeds downstream](#6-what-the-output-feeds-downstream)
7. [Limitations specific to DELLY TRA at ~5× coverage](#7-limitations-specific-to-delly-tra-at-5-coverage)
8. [Manuscript-ready prose snippets](#8-manuscript-ready-prose-snippets)
9. [Q&A](#9-qa)

---

## 1. Why this module exists (false-positive control + complex events)

MODULE_4F TRA has the **smallest direct contribution** to the inversion paper of any MODULE_4 module, but it has two specific roles:

**Role 1: False-positive control for inversion candidates.** A region of the genome that shows population-level inversion signal (MODULE_5A) and DELLY INV evidence (MODULE_4D) might in fact be an artifact of a reciprocal translocation. Reciprocal translocations create patterns of long-range linkage disequilibrium that can mimic inversion signals in local PCA, especially when the translocation involves segments of similar length on different chromosomes. If MODULE_4F shows a translocation with breakpoints near a candidate inversion's coordinates, the inversion call should be re-evaluated.

`[CONFIRM:` whether MODULE_5A2 or MODULE_5B currently cross-checks inversion candidates against TRA breakpoints. The data is available; the join may or may not be wired. `]`

**Role 2: Complex rearrangement detection.** Some real inversions are embedded inside reciprocal translocations or other complex rearrangements (e.g., inv-trans events common in cancer; rare but possible in germline). When this happens, the inversion's breakpoints will be flanked by translocation breakpoints in MODULE_4F. Detecting this allows the manuscript to flag the inversion as part of a complex event rather than a simple inversion.

In the catfish hybrid context, the expected number of true reciprocal translocations is small. Catfish hybrids have a stable karyotype across the cohort (chromosomal-level rearrangements do not segregate); MODULE_4F's catalog mostly contains low-frequency private translocations and call artifacts. The high-confidence TRA records are useful as cross-checks even though they are not numerous.

---

## 2. How DELLY identifies TRA records (BND filtered to CHR2 ≠ CHROM)

DELLY does not have a dedicated `delly call -t TRA` mode. **TRA is a derived class** — it is BND records (output of `delly call -t BND`) filtered to keep only those where the two breakends are on **different chromosomes**.

The MODULE_4F pipeline:

```
1. SLURM_A02 runs `delly call -t BND` per sample
   (DELLY_CALL_TYPE="BND" in 00_module4f_config.sh)

2. SLURM_A03 merges, regenotypes, applies strict filter

3. SLURM_A03 then applies the TRA-specific filter:
   bcftools view -i 'INFO/CHR2!=CHROM' input.bcf > tra_only.bcf
   (line 197 in SLURM_A03_merge_genotype.sh)

4. The output catalog contains only inter-chromosomal events
```

This means MODULE_4F and MODULE_4E share the **same upstream `delly call -t BND` invocation conceptually**, but they keep complementary subsets:
- MODULE_4E BND keeps `CHR2 == CHROM` (intra-chromosomal — inversion-rescue substrate)
- MODULE_4F TRA keeps `CHR2 != CHROM` (inter-chromosomal — translocations)

The two modules currently run independently (each does its own `delly call -t BND`), which is wasteful but operationally simple. A future optimization could share the discovery step.

### `selectTranslocations` in junction.h

The full dissection is in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2E](../SV_CALLING_FRAMEWORK_WIKI.md#2e--selecttranslocations--bnd--tra). Brief recap: `selectTranslocations` clusters split-read pairs where the two alignments are on different chromosomes and emits four classes based on strand × soft-clip combinations (one of `CT=3to3`, `5to5`, `3to5`, `5to3` per call). Unlike intra-chromosomal events (where `3to5` becomes DEL and `5to3` becomes DUP), all four classes survive as TRA records because the chromosome pair is unordered and every orientation combination is meaningful.

---

## 3. The four inter-chromosomal CT classes

Translocations have four geometrically distinct configurations, each represented by a different `CT` value:

| `CT` value | Geometry (chrA's break joined to chrB's break) | Common biological context |
|---|---|---|
| `3to3` | chrA's right side joined to chrB's right side (head-to-head) | Reciprocal translocation, one of the two breakpoints |
| `5to5` | chrA's left side joined to chrB's left side (tail-to-tail) | Reciprocal translocation, the partner breakpoint |
| `3to5` | chrA's right side joined to chrB's left side | Reciprocal translocation, alternate orientation |
| `5to3` | chrA's left side joined to chrB's right side | Reciprocal translocation, alternate orientation partner |

**A single reciprocal translocation produces TWO records** in the catalog (typically one `3to5` and one `5to3` at the same breakpoint pair, corresponding to the two chimeric chromosomes formed by the exchange). Counting unique reciprocal events therefore requires pairing records — usually by checking that two records share the same chromosome pair and have inverted breakpoint coordinates.

### How to identify reciprocal pairs

For every TRA record, check for a partner record matching:
- Same chromosome pair (`CHROM_A == CHROM_B'` and `CHR2_A == CHR2_B'` swapped)
- Breakpoint coordinates within a small window of each other
- Complementary CT classes (3to3 ↔ 5to5, or 3to5 ↔ 5to3)

If found → reciprocal translocation. If only one record exists → either a non-reciprocal break or a missed partner due to coverage/quality.

### Non-reciprocal breaks

Not all translocations are reciprocal. **Non-reciprocal translocations** (one segment moves to a new position without reciprocal exchange) appear as single TRA records without partners. These are typically associated with chromosome instability events (chromothripsis, breakage-fusion-bridge cycles) and are rare in germline contexts. In the catfish hybrid cohort, any non-reciprocal TRA call should be treated with extra scepticism — more likely a missed partner than a true non-reciprocal event.

---

## 4. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
PASS + INFO/PRECISE=1 + QUAL ≥ 300 + INFO/PE ≥ 3 + INFO/CHR2 ≠ CHROM
```

Source: `00_module4f_config.sh` lines 70–71 and `slurm/SLURM_A03_merge_genotype.sh` lines 123 + 197:

```bash
STRICT_REQUIRE_PASS=1
STRICT_REQUIRE_PRECISE=1
STRICT_MIN_QUAL=300
STRICT_MIN_PE=3

# In SLURM_A03 (two-stage filter):
# Stage 1: standard strict filter
STRICT_EXPR='INFO/PRECISE=1 && QUAL>=300 && INFO/PE>=3'
# Stage 2: TRA-specific filter
bcftools view -i 'INFO/CHR2!=CHROM' input.bcf
```

### Why these specific values

Identical to MODULE_4D INV and MODULE_4E BND (300 / 3 / PRECISE), and for the same reasons. See [MODULE_4D § 3](../MODULE_4D_delly_inv/WIKI.md#3-filter-thresholds-verified-from-code) for the full rationale on PE ≥ 3 vs higher.

The TRA-specific addition is the `INFO/CHR2 != CHROM` filter. Without it, MODULE_4F's catalog would be identical to MODULE_4E BND. With it, only true inter-chromosomal events remain.

---

## 5. Reading a TRA VCF record

<!-- TODO: replace synthetic example with a real DELLY TRA record from the production catalog -->

Synthetic example (a reciprocal translocation breakpoint on Scaffold02 ↔ Scaffold11):

```
Scaffold02  4521000  TRA00000018  T  T]Scaffold11:11203450]  340  PASS
PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv1.7.3;CHR2=Scaffold11;
END=11203450;PE=4;MAPQ=44;CT=3to3;CIPOS=-15,15;CIEND=-18,18;
SRMAPQ=58;INSLEN=0;HOMLEN=2;SR=3;SRQ=0.94715;
CONSENSUS=AGCTAGCAGTACGCATGCAGCTAACGAGGCAT...
GT:GL:GQ:FT:RR:RV:DR:DV
0/1:-12.45,-1.82,0:18:PASS:11:2:9:3
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS | Scaffold02:4,521,000 | Left side of the translocation breakpoint |
| ID | TRA00000018 | DELLY-assigned ID (note: SVTYPE is still BND in the VCF) |
| ALT | `T]Scaffold11:11203450]` | Bracket pattern → CT=3to3, mate at Scaffold11:11,203,450 |
| QUAL | 340 | Above 300 strict cutoff |
| `SVTYPE=BND` | string | DELLY emits all TRAs as SVTYPE=BND with CHR2 set |
| `CHR2=Scaffold11` | string | Mate chromosome — **different from CHROM → inter-chromosomal → TRA** |
| `END=11203450` | int | Mate position |
| `CT=3to3` | string | One of four inter-chromosomal CT classes |
| `PE=4`, `SR=3` | int | 4 paired-end + 3 split-read supporting (above 3/0 cutoffs) |
| `CIPOS=-15,15`, `CIEND=-18,18` | range | Slightly wider than INV examples (translocation context) |

Per-sample (one column shown):
- `GT=0/1` heterozygous translocation (one chromosome carries the rearrangement)
- `RR/RV = 11/2` — 2 alt junction reads
- `DR/DV = 9/3` — 3 alt discordant pairs
- Per-sample alt evidence = 2 + 3 = **5 reads**

### Identifying the partner record

For this TRA00000018 with `CT=3to3` connecting Scaffold02:4,521,000 to Scaffold11:11,203,450, look in the catalog for another record with:
- `CHROM=Scaffold11`, `POS≈11,203,450`, `CHR2=Scaffold02`, `END≈4,521,000`
- Complementary CT (5to5 typically pairs with 3to3 in the same reciprocal event)

If found, the two records describe the two chimeric chromosomes formed by the same reciprocal exchange. If only this record exists, the partner may have been filtered out (low quality), or it may genuinely be a non-reciprocal break.

---

## 6. What the output feeds downstream

```
MODULE_4F TRA catalog
  │
  ├─▶ Inversion candidate false-positive control   [CONFIRM: planned or coded?]
  │   - For each MODULE_5A inversion candidate:
  │     check whether any TRA record has a breakpoint within ±N kb of
  │     either inversion breakpoint
  │   - If yes → flag candidate as potentially translocation-related
  │     - might still be a real inversion embedded in a complex event
  │     - might be a translocation mis-identified as an inversion
  │   - Manuscript should flag any such candidate explicitly
  │
  ├─▶ Complex rearrangement detection
  │   - TRA + INV at the same breakpoint = likely inv-trans complex event
  │   - TRA + DEL at the same breakpoint = likely del-trans (chromothripsis-like)
  │   - These are rare in germline but worth flagging when found
  │
  └─▶ Reciprocal translocation pairing analysis
        - Programmatic detection of partner pairs (3to3+5to5 or 3to5+5to3)
        - Reports unique reciprocal events vs single-sided breaks
        - Single-sided breaks → either missed partners or non-reciprocal
          (latter is rare in germline)
```

For the broader plug-map, see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map).

---

## 7. Limitations specific to DELLY TRA at ~5× coverage

1. **TRA detection is sensitivity-bounded by both breakpoints having reads.** A reciprocal translocation has two breakpoints; both must produce enough discordant pairs and split reads to pass DELLY's clustering thresholds. At ~5×, frequently only one side passes, leaving an unpaired BND that may not be confidently classified as TRA. Some real translocations therefore appear as single-sided BNDs in MODULE_4E rather than as paired records in MODULE_4F.

2. **Inter-chromosomal calls have a higher false-positive rate than intra-chromosomal.** Mapping artifacts where reads map to repetitive regions on multiple chromosomes can produce spurious inter-chromosomal discordant signal. The strict cutoff (PRECISE + QUAL ≥ 300 + PE ≥ 3) removes most of these; remaining false positives typically have low MAPQ.

3. **Reciprocal partner pairing is not always trivial.** When the breakpoints are in repeat-rich regions, the assembled junction sequences from the two sides may not assemble cleanly enough for partner identification by sequence alone. Pairing relies on chromosome pair + breakpoint distance + complementary CT, which is mostly correct but can fail at SD-rich loci.

4. **Catfish hybrid TRA prevalence is biologically expected to be low.** The cohort represents an intra-species hybrid (Mac × Gar) with a stable karyotype. Reciprocal translocations are not expected to segregate in the cohort. Most strict-filter TRA calls likely reflect either (a) low-frequency private events, (b) rare assembly errors interpreted as translocations, or (c) artifacts of homologous-but-distant SDs across chromosomes. Each TRA call should be evaluated individually rather than treated as a true translocation by default.

5. **`delly filter -f germline` was designed for human cohorts.** The germline filter applies allele-frequency / Hardy-Weinberg-like criteria that may be inappropriate for the catfish cohort's structure (hatchery population, F1 hybrids violating HWE). For TRA calls specifically — where every call deserves manual inspection anyway — this is less of a concern than for DEL or DUP.

---

## 8. Manuscript-ready prose snippets

### Methods — DELLY TRA calling paragraph

> Translocations were identified by running DELLY2 v1.7.3 in BND mode (`delly call -t BND`) and filtering the resulting BND catalog to inter-chromosomal events (`INFO/CHR2 ≠ CHROM`). Per-sample discovery, cohort merging, regenotyping, subsetting to 81 unrelated individuals, and germline filtering followed the same workflow as MODULE_4D INV. The strict catalog required `FILTER = PASS`, `INFO/PRECISE = 1`, `QUAL ≥ 300`, `INFO/PE ≥ 3`, and the inter-chromosomal filter. Translocations have four possible orientation classes (`CT ∈ {3to3, 5to5, 3to5, 5to3}`), each representing a different geometric configuration of the joined chromosomes. Reciprocal translocations were identified by pairing records with complementary CT classes and matching breakpoint coordinates on swapped chromosome pairs.

### Results — number reporting placeholder

> The MODULE_4F TRA catalog contained [TODO N_TRA_STRICT] inter-chromosomal breakend records after strict filtering, of which [TODO N_RECIPROCAL] could be paired into [TODO N_PAIRED_EVENTS] reciprocal translocation events. The remaining [TODO N_UNPAIRED] records were single-sided breaks, likely representing translocations with missed partners due to coverage limitations rather than true non-reciprocal events. Cross-referencing with the MODULE_5A inversion candidate set identified [TODO N_TRA_NEAR_INV] inversion candidates with at least one TRA breakpoint within ±50 kb of either inversion breakpoint. Of these, [TODO N_FLAGGED] were flagged in the manuscript as potential complex events or translocation artifacts.

### Discussion — limitations footnote

> Translocation detection at ~5× coverage is limited by the requirement that both breakpoints of a reciprocal event independently pass evidence thresholds. Single-sided BND records (one breakpoint detected, partner not detected) are common at low coverage and cannot be reliably classified as reciprocal translocations from short-read data alone. The catfish hybrid cohort is biologically expected to have a low rate of segregating reciprocal translocations, so the small TRA catalog size is consistent with biology rather than necessarily reflecting reduced sensitivity.

---

## 9. Q&A

### Q: Why does TRA have four CT classes when intra-chromosomal events have only two for inversions?

For intra-chromosomal events, `selectInversions` produces only `3to3` and `5to5` (the two breakpoints of an inversion), and `selectDeletions` / `selectDuplications` produce `3to5` and `5to3` respectively (which are emitted as `<DEL>` and `<DUP>` rather than as BND). For inter-chromosomal events, all four CT classes are biologically possible — chromosomes can be joined in any of four orientation combinations, and there is no "deletion" or "duplication" reinterpretation because the segments are on different chromosomes. So the same junction.h logic produces 4 classes for TRA but only 2 for INV. See [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2E](../SV_CALLING_FRAMEWORK_WIKI.md#2e--selecttranslocations--bnd--tra) for the full dissection.

### Q: My inversion candidate has a TRA breakpoint nearby. Is it real?

Three possibilities, in rough order of likelihood:

1. **The TRA call is a false positive.** TRAs at ~5× coverage have a higher false-positive rate than other SV types. Check the TRA's MAPQ, evidence counts, and whether the breakpoint falls in a repeat-rich region. If suspicious, downweight or ignore it.
2. **The TRA is real, the inversion is real, and they are separate events that happen to have nearby breakpoints.** Possible but coincidence-driven; only worth considering if both calls are high-confidence and far enough apart that breakpoint-overlap is unlikely.
3. **The "inversion" is actually part of a complex translocation event.** This is the case worth flagging in the manuscript. An inv-trans complex event can be misidentified as a simple inversion by population-signal methods. If MODULE_4D shows weak INV evidence and MODULE_4F shows strong TRA evidence at overlapping coordinates, the complex-event interpretation is more parsimonious.

### Q: Why does MODULE_4F use `delly call -t BND` instead of `-t TRA`?

DELLY does not have a `-t TRA` mode — translocations are not a separately-typed SV class in DELLY. They are derived from BND records by post-hoc filtering on `CHR2 != CHROM`. The MODULE_4F config sets `DELLY_CALL_TYPE="BND"` for this reason, then applies the TRA-specific filter in SLURM_A03.

### Q: Why doesn't MODULE_4F share the BND discovery step with MODULE_4E?

It could in principle — both modules run identical `delly call -t BND` invocations and then filter complementary subsets. Currently they each run their own discovery, which is computationally wasteful (~2× compute on the BND step). The reason for the separation is operational simplicity — each module is independently runnable and debuggable. A future optimization could centralize BND discovery in a shared module and have 4E and 4F consume the same upstream output. Not a priority because the BND discovery is not the rate-limiting step in the overall pipeline.

### Q: What does CHR2 == CHROM mean in a BND record?

It means the breakpoint is intra-chromosomal — both ends of the rearrangement are on the same chromosome. These are the records MODULE_4E keeps (inversion-rescue substrate). MODULE_4F filters them out and keeps only `CHR2 != CHROM` (true inter-chromosomal events).

### Q: Is a single-sided TRA record (no reciprocal partner) real?

At ~5× coverage, almost certainly the partner is missing due to coverage limitations rather than a genuine non-reciprocal break. Non-reciprocal translocations are biologically rare in germline and usually associated with chromosome instability syndromes. For a catfish hybrid cohort, treat single-sided TRA records as **incomplete reciprocal events** rather than as genuine non-reciprocal translocations. This means the actual number of "translocation events" is closer to (paired + single-sided/2) than to the raw record count.

### Q: How does MODULE_4F interact with MODULE_4G Manta's translocation calls?

Currently they are not cross-referenced. Manta's BND catalog (post-conversion) contains translocations exclusively (intra-chromosomal inversion-orientation BNDs were converted to `<INV>`). DELLY TRA + Manta translocation concordance would be analogous to the DELLY×Manta INV concordance in MODULE_5A2 STEP05 but is not implemented. If translocation concordance becomes important for the manuscript, the implementation pattern from STEP05 could be adapted. For the current scope (TRA as false-positive control for inversions), DELLY-only is sufficient.

### Q: I see an extra `INFO/POS2` field in some records. What is it?

`POS2` is the position on the partner chromosome (`CHR2`). It duplicates the information already encoded in `INFO/END` for inter-chromosomal records but is sometimes emitted as a separate field for downstream tools that don't parse the bracket notation in ALT. For DELLY 1.7.3 output it should always equal `END`. Use whichever your downstream consumer prefers; both are correct.
