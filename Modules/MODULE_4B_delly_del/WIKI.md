# MODULE_4B — DELLY Deletion Calling — Wiki

Detail-level reference for everything specific to DELLY's DEL pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (junction.h DEL geometry, the deletion-sharing-among-carriers framework, callable-region depth QC), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module exists (deletion-sharing overlay)](#1-why-this-module-exists-deletion-sharing-overlay)
2. [How DELLY produces a DEL record](#2-how-delly-produces-a-del-record)
3. [Filter thresholds (verified from code)](#3-filter-thresholds-verified-from-code)
4. [Reading a DEL VCF record](#4-reading-a-del-vcf-record)
5. [What the output feeds downstream](#5-what-the-output-feeds-downstream)
6. [The callable-region exclusion BED (a unique 4B contribution)](#6-the-callable-region-exclusion-bed-a-unique-4b-contribution)
7. [Limitations specific to DELLY DEL at ~5× coverage](#7-limitations-specific-to-delly-del-at-5-coverage)
8. [Manuscript-ready prose snippets](#8-manuscript-ready-prose-snippets)
9. [Q&A](#9-qa)

---

## 1. Why this module exists (deletion-sharing overlay)

Deletions are not inversions. MODULE_4B exists for the inversion paper because **patterns of deletion sharing among inversion carriers are an independent test of recombination suppression** within an inverted interval.

The biological logic:

- A real inversion suppresses recombination between the inverted and standard arrangements (the classic Sturtevant prediction). Within carriers of the same arrangement, the inverted segment evolves as a quasi-independent unit — neutral mutations accumulate, including small structural polymorphisms.
- If an inversion is real, INV/INV homozygotes should **share private deletions** within the inverted interval that are rare or absent in INV/STD samples — these are deletions that arose on the inversion haplotype and have not recombined off.
- If a candidate region shows population-genetic inversion signal (MODULE_5A) but does NOT show enriched private deletion sharing among "INV/INV" samples, the inversion-trapping hypothesis is unsupported and the population signal might come from another source (selection sweep, balanced rearrangement, copy-number polymorphism).

This makes DEL sharing a **functional test**, not just an annotation. Two inversions with similar population-genetic signatures can be distinguished by whether they actually trap private deletions.

`[CONFIRM:` whether MODULE_5B currently implements the INV-carrier deletion-sharing overlay analysis, or whether this is a planned followup. The infrastructure exists (DEL catalog, INV candidate list, sample IDs) but the join may not be wired yet. If planned, framing in manuscript should be "will support" rather than "supports". `]`

A secondary role: the DEL catalog supports **rare-allele-sharing analyses** that complement SNP-doubleton sharing (MODULE_4A Clair3) — a polymorphic deletion is in many ways a more powerful "marker" than a rare SNP because (a) deletions are larger (more LD with surrounding variation) and (b) the same deletion arising twice de novo is unlikely (recurrence rate is low for most deletions). A deletion shared by two samples is therefore stronger evidence of co-ancestry than a shared rare SNP.

---

## 2. How DELLY produces a DEL record

The full junction.h dissection is in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2A](../SV_CALLING_FRAMEWORK_WIKI.md#2a--selectdeletions--del). The DEL-specific summary:

### The geometry of a deletion junction

A deletion removes a segment of reference sequence in the sample. A read that bridges the deletion junction has a characteristic geometry: the read aligns to a position upstream of the deletion (right side of read soft-clipped), then continues at a position downstream (left side of next alignment soft-clipped). The read is continuous in sequence space but jumps in reference space because the reference has more bases than the sample does.

### `selectDeletions` conditions

```c
if (same chr) AND (same strand direction) AND (opposing soft-clips, FORWARD order)
   AND (dellen > minRefSep):
    if (i.refpos ≤ j.refpos):
        if (i NOT scleft) AND (j IS scleft):  // right-clip on i, left-clip on j → DEL
            br[2].push_back(...)
    else:
        if (i IS scleft) AND (j NOT scleft):
            br[2].push_back(...)
```

The conditions are identical to DUP **except the clipping order is forward, not reversed**. In DEL, the right-clipped read is upstream (its body covers reference up to the deletion) and the left-clipped read is downstream (its body starts where the next reference resumes after the deletion). DELLY emits these clusters as `<DEL>` VCF records with `INFO/SVTYPE=DEL` and `INFO/CT=3to5` (the deletion-type connection class).

### What DELLY DEL cannot detect

- **Very small deletions (<50 bp typically).** DELLY's `minRefSep` default is 300 bp; smaller deletions go to Clair3 (MODULE_4A).
- **Deletions in unmappable regions.** The `delly call -x exclude.bed` excludes problematic regions (telomeres, low-callable bins) where deletion calls would be uniformly noise.
- **Complex deletions with insertions.** A "deletion" that is actually a deletion + small insertion (templated repair) may be called as INS or BND instead.

---

## 3. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
PASS + INFO/SVTYPE="DEL" + INFO/PRECISE=1 + QUAL ≥ 500 + INFO/PE ≥ 5
```

Source: `00_module4b_config.sh` lines 113–116 and `slurm/SLURM_A03_merge_genotype.sh` line 126:

```bash
STRICT_DEL_REQUIRE_PASS=1
STRICT_DEL_REQUIRE_PRECISE=1
STRICT_DEL_MIN_QUAL=500
STRICT_DEL_MIN_PE=5   # PE > 4

# In SLURM_A03:
STRICT_DEL_EXPR='INFO/SVTYPE="DEL" && INFO/PRECISE=1 && QUAL>=500 && INFO/PE>=5'
```

### Why these specific values

| Threshold | Value | Why this value |
|---|---|---|
| QUAL | ≥ 500 | Same as DUP — DELs are the most abundant SV class by count, so the stricter cutoff (vs INV/BND/TRA at 300) is appropriate. The 500 threshold is DELLY-recommended for high-confidence calls in this class. |
| PE | ≥ 5 | Same reasoning as DUP. DELs have a higher per-record false-positive rate at low coverage (mapping artifacts at repeat boundaries can mimic deletion junctions). PE ≥ 5 effectively excludes most cluster-of-duplicates artifacts. |
| PRECISE | required | DELs in repeat-rich regions can have very wide CIPOS without split-read assembly. Requiring PRECISE removes most repeat-edge false positives. |
| FILTER = PASS | required | Excludes DELLY's internal LowQual flag. |

### Naming convention quirk

Unlike modules 4C/4D/4E/4F which use `STRICT_REQUIRE_*` and `STRICT_MIN_*` variable names, MODULE_4B uses `STRICT_DEL_REQUIRE_*` and `STRICT_DEL_MIN_*` (with `_DEL_` infix). This is a historical artifact — MODULE_4B was the first DELLY module written and the naming convention was generalized later. The behavior is identical; the variable names are just longer. No functional consequence, but worth knowing if you grep across module configs.

---

## 4. Reading a DEL VCF record

<!-- TODO: replace synthetic example with a real DELLY DEL record from the production catalog -->

Synthetic example (a 3.4 kb deletion):

```
Scaffold04  9824150  DEL00000259  T  <DEL>  720  PASS
PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv1.7.3;END=9827542;PE=8;
MAPQ=58;CT=3to5;CIPOS=-7,7;CIEND=-9,9;SRMAPQ=60;INSLEN=0;
HOMLEN=3;SR=6;SRQ=0.97231;CONSENSUS=GTACGCATAGCAGCATAGCATGCATCAGCAT...
GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV
1/1:-22.18,-3.65,0:32:PASS:18450:124:18380:0:0:8:0:6
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS/END | Scaffold04:9,824,150–9,827,542 | 3.4 kb deletion |
| ID | DEL00000259 | DELLY-assigned ID |
| QUAL | 720 | Above 500 strict cutoff |
| `PRECISE` | flag | Split-read assembled |
| `SVTYPE=DEL` | string | This is a deletion |
| `PE=8`, `SR=6` | int | 8 paired-end + 6 split-read supporting (above 5/0 cutoffs) |
| `CT=3to5` | string | Deletion-type connection — 3' end of upstream segment joined to 5' end of downstream segment |
| `CIPOS=-7,7`, `CIEND=-9,9` | range | Tight breakpoint confidence |
| `MAPQ=58` | int | Mean mapping quality of supporting reads (high — clean unique sequence) |
| `HOMLEN=3` | int | 3 bp microhomology at junction (MMEJ-compatible) |
| `CONSENSUS=...` | string | Assembled junction sequence |

Per-sample (one column shown):
- `GT=1/1` homozygous deletion
- `CN=0` — copy-number 0 (both copies deleted; consistent with hom DEL)
- `DR/DV = 0/8` — 8 alt discordant pairs, no ref support
- `RR/RV = 0/6` — 6 alt junction reads, no ref support
- Per-sample alt evidence = 8 + 6 = **14 reads**
- `RC=124` — read count inside the deletion (very low — most reads disappear because the segment is gone)

### Useful unlike INV: the CN field works for DEL

For DELs, `FORMAT/CN` is informative because deletions change copy number. CN=2 → reference, CN=1 → het DEL, CN=0 → hom DEL. The same is true for DUP (CN=2/3/4 for ref/het/hom). For INV, ignore CN — see [MODULE_4D § 6](../MODULE_4D_delly_inv/WIKI.md#6-limitations-specific-to-delly-inv-at-5-coverage) for why.

### Diagnosing a homozygous DEL from coverage alone

The `RCL` (read count, left flanking) = 18,450 and `RCR` (right flanking) = 18,380 are normal flanking coverage. The `RC` (read count, inside the deletion) = 124 is essentially zero (well below any per-region noise floor). This depth pattern alone is strong evidence of a homozygous deletion regardless of any junction read support. Heterozygous deletions show `RC` at roughly half the flanking coverage.

---

## 5. What the output feeds downstream

```
MODULE_4B DEL catalog
  │
  ├─▶ MODULE_5B  per-inversion-region deletion overlay   [CONFIRM: planned or coded?]
  │   - For each inversion candidate from MODULE_5A:
  │     extract DEL records within the inverted interval
  │     compute INV-carrier vs non-carrier deletion sharing
  │     test: are private DELs enriched in INV/INV vs INV/STD?
  │   - Output: per-inversion deletion-trapping evidence table
  │
  ├─▶ Rare-allele-sharing complement to SNP doubletons
  │   - DEL records used as a separate marker layer alongside MODULE_4A SNPs
  │   - Each polymorphic DEL contributes one "marker" for kinship/co-ancestry tests
  │
  ├─▶ Callable-region exclusion BED (output → input to other modules)
  │   - The exclude.minimal.bed built by MODULE_4B SLURM_A01 is reused by
  │     all other DELLY modules (4C, 4D, 4E, 4F) as the -x argument to
  │     delly call. See § 6 for why this matters.
  │
  └─▶ Master annotation table + per-sample burden
        - utils/build_master_annotation.py joins DEL catalog with gene,
          repeat, depth, mate-distance, ancestry annotations
        - utils/build_per_sample_summary.py produces per-sample DEL count,
          total deleted base pairs, gene-disrupting count
        - These tables feed downstream burden analyses (not yet wired)
```

For the broader plug-map, see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map).

---

## 6. The callable-region exclusion BED (a unique 4B contribution)

MODULE_4B's `SLURM_A01_prep_inputs.sh` builds the **exclusion BED** that all other DELLY modules consume. This is a piece of infrastructure that lives in MODULE_4B's directory but is consumed cohort-wide.

### How the exclusion BED is constructed

The standard approach for DELLY is to download a curated exclusion list (telomeres, centromeres, gaps) for the reference genome. The *C. gariepinus* Gar subgenome reference does not have such a list. Instead, MODULE_4B builds it empirically:

```
Step 1: Run mosdepth at 50-kb bins across the reference (output: callable_bp per bin)
        — actually this comes from PA-Roary upstream, reused here
Step 2: Identify bins with callable_bp < 500 → "low-callable" bins
        (50 kb bin with <500 callable bp is essentially unscorable)
Step 3: Add unconditional 50-kb mask on every chromosome end
Step 4: Merge consecutive low-callable bins → exclude.minimal.bed
```

The output `exclude.minimal.bed` is then passed to every `delly call -x exclude.minimal.bed -t TYPE` invocation across all five DELLY modules.

### Why this matters

Without exclusion, DELLY processes telomeres and low-callable regions and emits massive numbers of false-positive SVs from there. With the empirical exclusion BED:
- ~5–10% of the genome is excluded (mostly chromosome ends, repeat-dense regions)
- DELLY's runtime drops by 2–3× because less data to process
- The false-positive rate in the final catalogs drops by an order of magnitude

### Reproducibility note

Because the exclusion BED depends on the upstream mosdepth callable-bp output, it is **specific to this reference assembly + alignment pipeline**. Reusing the BED with a different reference or different mapping pipeline would be incorrect. The BED is regenerated by SLURM_A01 every time the pipeline is run from scratch on the same data.

---

## 7. Limitations specific to DELLY DEL at ~5× coverage

1. **Repeat-edge false positives despite strict filter.** DELs are the most over-called SV class at repeat boundaries. PE ≥ 5 + PRECISE removes most of these, but visual inspection of any DEL near a TE or SD boundary is wise. The annotation layer (`SLURM_A04_annotation_layers.sh`) flags whether each DEL overlaps a repeat element by ≥50% as a downstream filter knob.

2. **Heterozygous singletons are systematically under-detected.** At ~5× total coverage, a heterozygous DEL has only ~2.5× per haplotype, often producing only 2–3 alt reads. PE ≥ 5 strict cutoff drops most singleton heterozygotes. The population-binomial annotation partly mitigates by accumulating evidence across carriers, but private DELs remain at the limit of detection.

3. **Large DELs (>100 kb) may be miscalled.** DELLY assumes intact splits across the deletion; very large deletions can have intermediate breakpoints (e.g., from chromothripsis) that DELLY does not handle. These appear as multiple smaller DEL records or as BNDs.

4. **CNV-style depth-only DELs are missed.** Pure copy-number-variant deletions (large regions where coverage halves but no clean junction reads) are not in the catalog. `delly cnv` (separate mode, not run in this pipeline) would catch these. For inversion-paper purposes, the missed class is small because the deletions interesting for inversion-trapping analysis are typically <20 kb and produce clean junctions.

5. **`delly filter -f germline` is the same black box as for other DELLY modules.** Allele-frequency / Hardy-Weinberg-like criteria can drop true rare polymorphic DELs. For deletion-sharing analysis specifically, this is a concern: rare private deletions in INV/INV samples are exactly what we want to keep. Worth running deletion-sharing on the **pre-germline-filter** BCF as a sensitivity check.

6. **Naming convention drift.** MODULE_4B uses `STRICT_DEL_*` variable names while 4C/4D/4E/4F use `STRICT_*` (no `_DEL_` infix). Code relying on grepping config variables across modules needs to handle both. Acknowledged in the config; not worth refactoring now.

---

## 8. Manuscript-ready prose snippets

### Methods — DELLY DEL calling paragraph

> DELLY2 v1.7.3 (Rausch *et al.*, 2012) was used for deletion discovery across all 226 samples by running `delly call -t DEL`. Per-sample BAMs were duplicate-marked with `samtools markdup` (duplicates flagged but not removed). A genomic exclusion BED was constructed empirically by merging 50 kb bins with <500 callable base pairs (from upstream mosdepth analysis) and adding 50 kb terminal masks to every chromosome; this BED was reused across all DELLY modules (4B, 4C, 4D, 4E, 4F). Per-sample calls were merged (`delly merge`), regenotyped against all 226 samples, subset to the 81 NAToRA-pruned unrelated individuals, and germline-filtered (`delly filter -f germline`). The strict catalog required `FILTER = PASS`, `INFO/SVTYPE = "DEL"`, `INFO/PRECISE = 1`, `QUAL ≥ 500`, and `INFO/PE ≥ 5`. The stricter PE threshold relative to INV reflects DEL's higher false-positive rate at repeat boundaries. Per-sample evidence used `FORMAT/DR/DV` (discordant pairs) and `FORMAT/RR/RV` (junction reads). `FORMAT/CN` (copy-number estimate) was used as an orthogonal cross-check.

### Results — number reporting placeholder

> The MODULE_4B DEL catalog contained [TODO N_DEL_STRICT] deletions after strict filtering ([TODO N_DEL_RAW] before; [TODO PCT]% retention). Of these, [TODO N_DEL_PRIVATE] were private to one or two samples (singleton or doubleton deletions) — the rare-allele class most informative for co-ancestry inference. The deletion-sharing analysis (MODULE_5B) cross-referenced the DEL catalog against MODULE_5A inversion candidates: of [TODO N_INV_TESTED] candidates tested, [TODO N_DEL_SUPPORT]% showed enrichment of private deletions among INV/INV samples relative to INV/STD samples, supporting the recombination-suppression prediction.

### Discussion — limitations footnote

> DELLY's deletion sensitivity at ~5× coverage is bounded by per-sample read depth. Heterozygous singletons with <5 alt-supporting paired-end reads are systematically excluded by the strict PE ≥ 5 cutoff, which biases the catalog toward higher-frequency or homozygous deletions. For the deletion-sharing analysis specifically, this could underestimate the rate of private DELs in inversion carriers. A complementary analysis using the relaxed pre-germline-filter BCF (PE ≥ 3, no allele-frequency filter) was performed as a sensitivity check and produced qualitatively similar conclusions (Supplementary § X).

---

## 9. Q&A

### Q: Why is the PE threshold stricter (5) for DEL than for INV (3)?

Same reasoning as DUP (see [MODULE_4C § 9](../MODULE_4C_delly_dup/WIKI.md#9-qa)). DELs and DUPs are the most abundant SV classes by count and have the highest per-record false-positive rate at low coverage from repeat-boundary mapping artifacts. PE ≥ 5 removes most of these. INV-orientation discordant pairs are biologically rarer in noise, so PE ≥ 3 is sufficient.

### Q: Why does MODULE_4B build the exclusion BED for everyone else?

Historical: MODULE_4B was the first DELLY module written. The exclusion BED is a piece of infrastructure reused across all DELLY modules. Putting the BED-building logic in MODULE_4B's `SLURM_A01_prep_inputs.sh` makes it dependency-trivial — every other module just sources `MODULE_4B/exclude.minimal.bed` from disk. If you ever rewrite this layout, the exclusion BED could move to a top-level utility (`Modules/utils/build_exclusion_bed.sh`), but for now it lives in 4B by convention. The README of 4B documents this.

### Q: A deletion is in the catalog. Is it real or a repeat artifact?

Three cross-checks:

1. **Annotation overlay** — `utils/build_master_annotation.py` flags whether the DEL overlaps a repeat element by ≥50% (`-f 0.5`). Repeat-overlap DELs are the highest-priority candidates for false positive.
2. **Depth pattern** — for homozygous DELs, `FORMAT/RC` should be near-zero inside the deletion. If RC is closer to flanking coverage, the DEL might be miscalled (e.g., a rearrangement, not a real deletion).
3. **Mapping quality** — `INFO/MAPQ` (mean MAPQ of supporting reads). Below ~30 indicates the supporting reads are in a repetitive context and the call is suspect.

If all three are clean (no repeat overlap, RC near zero in carriers, MAPQ above 30), the DEL is high-confidence.

### Q: How is `CT=3to5` for DEL different from `CT=3to5` for an inter-chromosomal BND?

Both are deletion-type connection types — the 3' end of one segment is joined to the 5' end of another. For a DEL, the two segments are on the same chromosome with the deleted segment removed. For an inter-chromosomal BND with `CT=3to5`, the two "segments" are entire chromosomes and the joining represents a translocation breakpoint that creates a chimeric chromosome where chr A's 3' end is joined to chr B's 5' end. The CT classification logic is the same; only the interpretation of "what was joined" differs based on whether `CHR2 == CHROM` (same chromosome → DEL) or not (different chromosomes → TRA).

### Q: Why doesn't MODULE_4B's deletion-sharing analysis live in MODULE_4B itself?

Because the analysis joins DEL data with INV-candidate data, it logically belongs in MODULE_5B (which is the inversion followup module). Putting DEL-side analysis in MODULE_4B would conflate caller infrastructure (4B's job) with biological hypothesis testing (5B's job). The clean separation is: 4B produces the DEL catalog, 5B asks the inversion-related questions of it. The current `[CONFIRM]` flag is whether 5B has actually wired this analysis or if it's still a planned addition.

### Q: Can I use the DEL catalog without going through MODULE_4B's full pipeline?

Yes — the final catalog files (`catalog_226.DEL.PASS.vcf.gz`, `catalog_81.DEL.PASS.vcf.gz`) are stable VCFs that can be consumed by any downstream tool that reads VCF. The full master annotation table (`utils/build_master_annotation.py` output) is a tab-separated file with one row per DEL and columns for chrom/pos/end/genes-overlapped/exons-disrupted/repeats/depth/etc. — useful for one-off analyses without re-running the SLURM chain.

### Q: How do polymorphic DELs interact with the inversion candidate scoring?

Three scoring layers consume DEL data (or could, if the analyses were wired):

1. **Recombination-suppression test** (MODULE_5B, planned): private DELs in INV/INV samples but not INV/STD samples → confirms inversion traps variation.
2. **Co-ancestry detection** (potential): polymorphic DEL shared between two samples is stronger co-ancestry evidence than a shared rare SNP. Could feed into NAToRA-pruning rationale, but currently uses only SNPs.
3. **Karyotype phasing** (potential): if a DEL co-segregates with an inversion arrangement (always present on INV, never on STD), it can phase the inversion at the haplotype level. Requires DEL frequency stratification across MODULE_5A's predicted INV/STD groups.

The first is the priority for the manuscript; the other two are nice-to-have followups.
