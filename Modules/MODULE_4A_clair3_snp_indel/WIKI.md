# MODULE_4A — Clair3 SNP & INDEL Discovery — Wiki

Detail-level reference for everything specific to Clair3-based small-variant calling. Pairs with `README.md` (high-level + how to run). Note: MODULE_4A is fundamentally different from the other MODULE_4* modules — it does not call structural variants and shares no upstream code with DELLY or Manta. It is grouped under MODULE_4 because it operates on the same input BAMs and feeds many of the same downstream consumers.

For shared cross-cutting content (the SV calling framework that does NOT apply to this module), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md). For the broader context of how small variants feed into the inversion analysis (founder packs, deleterious-variant scoring), this wiki is the central reference.

---

## Table of contents

1. [Why this module exists (per-sample small variants for downstream)](#1-why-this-module-exists-per-sample-small-variants-for-downstream)
2. [Why Clair3 and not ANGSD or GATK](#2-why-clair3-and-not-angsd-or-gatk)
3. [How Clair3 produces variant calls](#3-how-clair3-produces-variant-calls)
4. [The five-phase post-processing pipeline](#4-the-five-phase-post-processing-pipeline)
5. [The three-branch rescue system](#5-the-three-branch-rescue-system)
6. [Reading a Clair3 VCF record](#6-reading-a-clair3-vcf-record)
7. [What the output feeds downstream](#7-what-the-output-feeds-downstream)
8. [Limitations specific to Clair3 at ~5× coverage](#8-limitations-specific-to-clair3-at-5-coverage)
9. [Manuscript-ready prose snippets](#9-manuscript-ready-prose-snippets)
10. [Q&A](#10-qa)

---

## 1. Why this module exists (per-sample small variants for downstream)

MODULE_4A produces **per-sample SNP and indel calls with hard phased genotypes**. The inversion paper depends on this catalog for two analyses that ANGSD's biSNP-only output cannot support:

**Founder-pack fragment classification (MODULE_6).** The 100-INDEL sliding window framework defines fragment classes based on per-sample indel haplotype patterns. Each window of 100 indels in pipeline order forms a fragment; samples with similar indel patterns within a fragment cluster into "founder-like packs". This requires (a) per-sample hard genotypes (not likelihoods), (b) indels (not just SNPs), and (c) phased data so that fragments can be resolved per haplotype. ANGSD provides none of these.

**Deleterious-variant scoring (MODULE_CON).** Annotation tools like SnpEff (functional impact), SIFT4G (amino-acid substitution effect), and VESM (zero-shot LLM-based scoring) operate on individual variant calls with hard genotypes. The output is the BC0–BC4 (breeding concern) classification overlaid on inversion candidate regions to identify inversions harboring elevated deleterious burden — a specific manuscript claim about which inversions have functional consequences.

`[CONFIRM:` founder packs use indels only, not SNPs+indels — based on memory context "100-INDEL windows" suggests indels only, but MODULE_4A produces both. The actual MODULE_6 input parser would resolve this. `]`

A secondary use: per-sample hard genotypes from MODULE_4A complement ANGSD's likelihood-based output by allowing per-sample analyses (rare-allele sharing, identity-by-descent, founder-haplotype reconstruction) that need definitive genotypes rather than probability distributions.

---

## 2. Why Clair3 and not ANGSD or GATK

The pipeline runs three small-variant callers in parallel for different purposes:

| Caller | Output type | Used for | Why not the others |
|---|---|---|---|
| **ANGSD** (MODULE_2A) | Genotype likelihoods at biallelic SNPs | Population genetics (Q/Fst/PCA/admix) | Doesn't call indels; doesn't produce hard genotypes |
| **Clair3** (MODULE_4A, this) | Per-sample hard SNP+indel calls, phased | Per-sample analyses (founder packs, deleterious scoring, ROH) | Slower than ANGSD; not designed for cohort-level GL operations |
| **GATK HaplotypeCaller** | Per-sample hard SNP+indel calls | (NOT in production pipeline) | Slower than Clair3 at low coverage; deep-learning-free model is less sensitive to coverage variation |

Clair3 specifically is chosen because:

1. **Designed for both short-read and long-read calling** with consistent accuracy across coverage levels. At ~5×, Clair3's deep-learning model handles coverage-dependent uncertainty better than HMM-based callers.
2. **Outputs phased VCFs natively** when paired with WhatsHap. Phased indels are essential for the founder-pack fragment-class framework.
3. **Faster than HaplotypeCaller** for cohort-scale per-sample calling — important when calling 226 samples × 28 chromosomes.
4. **Released model checkpoints for non-human species** — Clair3 has Illumina-trained models that work reasonably on non-human reference genomes without requiring custom training.

The trade-off is that Clair3's tensor-based filtering occasionally produces marginal-confidence indel calls in repeat-rich regions that benefit from manual rescue (see § 5 — the three-branch rescue system).

---

## 3. How Clair3 produces variant calls

Clair3's algorithm is described in detail in Zheng *et al.* 2022 (Nature Computational Science). The core architecture relevant to MODULE_4A:

### The two-stage neural network

For each candidate variant site, Clair3 applies two sequentially-trained neural networks:

1. **Pileup network**: examines the read pileup at the candidate position, classifying as REF / SNP / INSERTION / DELETION / COMPLEX. Fast first pass.
2. **Full-alignment network**: for sites flagged as variant by the pileup network, examines the local alignment context (~70 bp window around the candidate) for confirmation. Slower but more accurate, especially for indels.

Both networks output a probability distribution over genotype classes; the consensus call is the highest-probability class with associated quality (`QUAL`) and confidence (`FT` filter flags).

### Per-sample, per-chromosome execution

MODULE_4A runs Clair3 in parallel as a SLURM array: 226 samples × ~28 chromosomes = up to ~6,300 individual Clair3 jobs. Each emits a per-sample-per-chr VCF that is later concatenated. The pipeline's `SLURM_D02_clair3_discovery.sh` orchestrates this with a per-sample manifest and a dispatch table for chromosome batches.

### WhatsHap phasing

After Clair3 emits unphased VCFs, **WhatsHap** is invoked per sample to phase the variants using read-based phasing (which informative read carries which alleles). The output is a VCF with `FORMAT/PS` (phase set) and `|`-delimited GT (e.g., `0|1` instead of `0/1`). Phased blocks are bounded by phase-set switches — informally, "regions where we can confidently follow one haplotype through the variant set".

For founder-pack analysis, phase blocks are critical: each fragment class is defined within phase-block boundaries, so the founder-pack framework operates per-haplotype rather than per-genotype.

---

## 4. The five-phase post-processing pipeline

The raw Clair3 + WhatsHap output is per-sample-per-chromosome phased VCFs. The post-processing pipeline transforms these into the final cohort-wide phased catalog through five sequential phases.

### Phase 1 — Per-sample preparation (SLURM array × 226)

- **STEP_P01** — parse and annotate VCF (add INFO fields, normalize representation)
- **STEP_P02** — local event block detection (group nearby variants into "events")
- **STEP_P02B** — phase-aware blocks using WhatsHap PS + read-pair rescue for ambiguous boundaries
- **STEP_P03** — strong single-sample rescue (see § 5)
- **STEP_P04** — export weak indel candidates for population-level review

### Phase 2 — Population-level aggregation (single job)

- **STEP_P05** — group weak candidates across samples (one indel weak in 1 sample, supported in 30 others, becomes a regenotype-target catalog)
- **STEP_P06** — prepare regenotype catalog (write the joint site list with all samples' weak indels)

### Phase 3 — Cohort regenotyping (SLURM array × 226)

- **STEP_P07A** — for each sample, regenotype against the shared catalog. Forces a hard call (`0/0`, `0/1`, or `1/1`) at every site in the catalog regardless of whether Clair3 originally called it in that sample.

### Phase 4 — Merge regenotype outputs (single job)

- **STEP_P07B** — merge all 226 regenotyped VCFs, produce the cohort-level summary statistics

### Phase 5 — Per-sample final processing (SLURM array × 226)

- **STEP_P08–P10** — six-class final classification per sample (combining original calls + regenotype results), compute per-sample variant burden, produce the final per-sample VCF and the cohort-level shared catalog

The five-phase structure exists because Clair3's per-sample confidence is uneven. Phase 1 cleans up obvious errors. Phase 2 identifies weak calls that might be real but were unconfident. Phase 3 explicitly tests these against the BAM data (regenotyping is more accurate than discovery). Phase 4 aggregates. Phase 5 produces the canonical output that downstream consumers read.

---

## 5. The three-branch rescue system

The "three-branch rescue" is a Clair3-specific intervention that addresses systematic Clair3 failure modes at ~5× coverage.

The three branches:

| Branch | Triggers when | What it does | Where in the pipeline |
|---|---|---|---|
| **Branch A — strong-evidence rescue** | Clair3 says LOW QUAL at a site, but the BAM shows ≥3 high-quality alt reads in the same sample | Override Clair3's filter, emit as PASS with rescue tag | STEP_P03 |
| **Branch B — phase-aware rescue** | A weak call falls within a strong WhatsHap phase block | Use the surrounding phase context to restore confidence | STEP_P02B |
| **Branch C — population-evidence rescue** | A weak call in sample X is strongly supported in many other samples | Regenotype sample X against the population catalog, force a hard call | STEP_P03–P07A |

Each branch operates on a different type of marginal-confidence call. Together they recover indels that Clair3's tensor model marked as `LowQual` despite real underlying evidence — a substantial fraction at ~5× coverage.

### Why this matters for the inversion paper

Indels near inversion breakpoints are disproportionately in the marginal-confidence class because:

1. **Breakpoints are often in repetitive sequence** → Clair3's confidence drops in repeats
2. **Inversion-adjacent regions have unusual coverage patterns** → Clair3's coverage-dependent calibration can be off
3. **Microhomology-mediated breakpoints leave small indels right at the junction** → these are exactly the indels we want for breakpoint validation but Clair3 often misses

The rescue system is therefore directly useful for inversion-adjacent variant analysis. The three-branch structure was added explicitly to recover breakpoint-proximal indels.

### Validation

Each rescued call carries an INFO tag indicating which branch rescued it. Downstream consumers (MODULE_6 founder packs, MODULE_CON deleterious scoring) can choose to include or exclude rescued calls based on use-case sensitivity requirements. For high-stringency analyses, exclude all rescues. For sensitivity-prioritized analyses (rare-allele sharing, breakpoint context characterization), include them.

---

## 6. Reading a Clair3 VCF record

<!-- TODO: replace synthetic example with a real Clair3 record from the production catalog -->

Synthetic example (a heterozygous SNP with phasing):

```
Scaffold03  4521234  .  C  T  35  PASS
.
GT:GQ:DP:AD:AF:PS
0|1:23:8:5,3:0.375:4521198
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS | Scaffold03:4,521,234 | Variant position |
| REF / ALT | C / T | Single-nucleotide variant |
| QUAL | 35 | Clair3's combined quality score (phred-like) |
| FILTER | PASS | Passed Clair3's confidence threshold |

Per-sample (one column shown):
- `GT=0|1` heterozygous, **phased** (the `|` is the key indicator vs `/` for unphased)
- `GQ=23` genotype quality
- `DP=8` total read depth at this site (low — characteristic of ~5× coverage)
- `AD=5,3` allelic depth: 5 ref-supporting reads, 3 alt-supporting
- `AF=0.375` allele frequency in this sample (alt count / total = 3/8)
- `PS=4521198` phase set ID — this variant is in the same phase block as the variant at 4,521,198 (so haplotype follows from there)

### Reading an indel record

<!-- TODO: replace with real indel example -->

```
Scaffold03  5102456  .  ATGC  A  28  PASS
INFO/RESCUE_BRANCH=C;INFO/REGEN_SOURCE=population
GT:GQ:DP:AD:AF:PS
1|1:32:6:0,6:1.0:5102220
```

This is a **3 bp deletion** (`ATGC` → `A` removes `TGC`) at Scaffold03:5,102,456. Per-sample:
- `GT=1|1` homozygous deletion, phased
- `INFO/RESCUE_BRANCH=C` — rescued via Branch C (population-evidence regenotyping)
- `INFO/REGEN_SOURCE=population` — the call came from STEP_P07A regenotyping, not original Clair3 discovery

For founder-pack analysis, this indel would be one entry in the 100-INDEL window framework. The phase information (`1|1`, `PS=5102220`) lets MODULE_6 assign it to a specific phase-block-defined fragment.

---

## 7. What the output feeds downstream

```
MODULE_4A Clair3 catalog (per-sample phased VCFs + cohort shared catalog)
  │
  ├─▶ MODULE_6 founder packs                       [CONFIRM: indels only or SNPs+indels?]
  │   - 100-INDEL sliding windows define fragments
  │   - Per-sample fragment patterns cluster into founder-like packs
  │   - Phase-block boundaries determine fragment boundaries
  │   - Output: per-sample fragment-class assignments, founder-pack rosters
  │
  ├─▶ MODULE_CON deleterious-variant scoring
  │   - SnpEff for functional annotation (synonymous / missense / frameshift etc)
  │   - SIFT4G for amino-acid substitution effect prediction
  │   - VESM for zero-shot LLM-based deleteriousness scoring (Dinh et al. 2026)
  │   - Splice-site annotation for frameshift / splice-site indels
  │   - Output: BC0–BC4 (breeding concern) classification per variant, per sample
  │
  ├─▶ MODULE_3 ROH detection (indirect)
  │   - ngsF-HMM uses biallelic SNPs (from ANGSD MODULE_2A primarily)
  │   - But Clair3 phased data can supplement ROH boundary refinement
  │
  └─▶ Manuscript Figure (per-inversion deleterious burden)
        - Inversion intervals overlaid with BC1–BC4 variant counts
        - Per-arrangement burden comparison (INV/INV vs INV/STD vs STD/STD)
        - Breakpoint-proximal vs central deleterious distribution
```

The MODULE_4A output is the **foundation for all per-sample analyses** that need hard genotypes. ANGSD output cannot substitute for any of these consumers.

---

## 8. Limitations specific to Clair3 at ~5× coverage

1. **Indel sensitivity is bounded by per-sample coverage.** At ~5×, a heterozygous indel typically has only 2–3 alt-supporting reads. Clair3's full-alignment network requires sufficient context to confidently call indels, and confidence drops sharply below ~6× per-sample. The three-branch rescue system partially mitigates this for population-supported indels, but private heterozygous indels remain at the limit of detection.

2. **Repeat-region indels are systematically under-called.** Clair3's training data has limited representation of indels in long tandem repeats, microsatellites, and near transposable elements. The catalog will be sparse in these regions even though biologically these are indel-rich. The PA-Roary / MODULE_3 pipelines handle some of this via mosdepth callable-region masking, but Clair3 itself cannot recover what it doesn't call.

3. **WhatsHap phase blocks are short at ~5× coverage.** Phase blocks are bounded by informative reads spanning consecutive variants. At low coverage, fewer reads span variant pairs, so phase blocks are typically <10 kb. Long-range phasing across an inversion (which would help validate the inversion structure) is generally not available. For inversion analysis, this means founder-pack fragments are smaller than ideal — a 100-INDEL window may span multiple phase blocks at low coverage.

4. **The model is trained on Illumina human data.** Clair3's released checkpoints are tuned for human Illumina sequencing. Catfish-specific calibration was not performed (would require a truth set, which doesn't exist for catfish). Per-sample concordance against PCR validation would help quantify but has not been done.

5. **Clair3 has no cohort-level joint calling.** Unlike GATK's joint genotyping or DeepVariant's GVCF-based aggregation, Clair3 calls each sample independently. The pipeline's regenotyping (STEP_P07A) partially mitigates by forcing per-sample genotypes at all population sites, but this is regenotyping after independent discovery, not true joint calling. Cohort-level coverage information that could rescue per-sample weak calls is not used by Clair3 directly.

6. **The per-sample classification (BC0–BC4) requires SnpEff annotation.** SnpEff databases are reference-specific. MODULE_CON's scoring pipeline depends on a SnpEff database built for the catfish hybrid reference; quality of annotation depends on the gene model quality (`fClaHyb_Gar_LG.from_CGAR.gff3_polished.sorted.gff3`).

---

## 9. Manuscript-ready prose snippets

### Methods — Clair3 calling paragraph

> Per-sample SNP and indel discovery was performed with Clair3 v1.0.x (Zheng *et al.*, 2022) on duplicate-marked, MAPQ-filtered BAMs (MODULE_1). Each sample was called per-chromosome in parallel using SLURM arrays (226 samples × 28 chromosomes). Clair3's two-network architecture (pileup network for first-pass classification, full-alignment network for confirmation) was used with the released Illumina-trained model checkpoint. Per-sample VCFs were phased with WhatsHap v2.x using read-based phasing. A five-phase post-processing pipeline was then applied: (i) per-sample VCF parsing and event-block detection; (ii) population-level aggregation of weak indel candidates; (iii) cohort regenotyping forcing hard calls at all population sites; (iv) merge of regenotyped per-sample VCFs; (v) six-class classification combining original and regenotyped calls. A three-branch rescue system recovered marginal-confidence calls that had strong single-sample, phase-aware, or population-level evidence (Methods § Three-branch rescue). The final cohort catalog comprised per-sample phased VCFs and a cohort-shared site list with rescue-source provenance tags.

### Results — number reporting placeholder

> Clair3 produced [TODO N_RAW_SNP] SNP and [TODO N_RAW_INDEL] indel candidate sites across the cohort in initial discovery. After the five-phase post-processing pipeline including the three-branch rescue, the final catalog contained [TODO N_FINAL_SNP] SNPs and [TODO N_FINAL_INDEL] indels. Of these, [TODO N_RESCUED]% were rescued from initial LowQual filtering by at least one branch of the rescue system, predominantly via Branch C (population-evidence regenotyping). The catalog feeds two downstream analyses: founder-pack fragment classification (MODULE_6) using 100-indel sliding windows over the phased indel catalog, and deleterious-variant scoring (MODULE_CON) producing the BC0–BC4 classification used for per-inversion burden analysis.

### Discussion — limitations footnote

> Clair3's sensitivity for heterozygous indels is bounded by per-sample read depth. At ~5× coverage, a substantial fraction of true heterozygous indels produce only 2–3 alt-supporting reads, below Clair3's full-alignment confidence threshold for direct calling. The three-branch rescue system recovers indels with cross-sample population support, but private heterozygous indels remain at the limit of detection. WhatsHap phase blocks are typically <10 kb at this coverage, limiting per-haplotype analyses to local fragments rather than chromosome-level haplotype reconstruction. These limitations are acceptable for the founder-pack and deleterious-variant analyses that consume the catalog — both operate at fragment scale rather than requiring chromosome-spanning continuity.

---

## 10. Q&A

### Q: Why call SNPs with Clair3 here when MODULE_2A already calls them with ANGSD?

The two catalogs serve different purposes. ANGSD outputs **genotype likelihoods** (per-site, per-sample probability distributions over `0/0`, `0/1`, `1/1`) that propagate uncertainty correctly through population-genetic estimators (Q, Fst, π). It does not output hard genotypes — every per-sample call remains a probability distribution. ANGSD also doesn't call indels.

Clair3 outputs **hard genotypes** at per-sample resolution. Hard genotypes are required for any downstream analysis that needs definitive per-sample calls: founder packs (which sample carries which fragment), deleterious-variant scoring (does this sample have this variant?), per-sample variant burden, individual-pair IBD detection.

Both are correct for their use cases. Mixing them up — using ANGSD likelihoods where hard genotypes are needed, or hard genotypes where likelihoods would propagate uncertainty better — introduces bias. The pipeline runs both for this reason.

### Q: What does the `INFO/RESCUE_BRANCH` tag mean?

It indicates which of the three rescue branches recovered the call. `A` = strong-evidence rescue (single sample's BAM showed clear evidence despite Clair3's LowQual flag). `B` = phase-aware rescue (the call falls in a strong WhatsHap phase block that supports it). `C` = population-evidence rescue (the call was originally weak in this sample but strongly supported in many other samples; regenotyping forced a hard call). Records without the tag were not rescued — they were called confidently by Clair3 in original discovery.

For high-stringency analyses (e.g., manuscript Table 1 of "highly confident variants"), exclude all rescue-tagged calls. For sensitivity-prioritized analyses (rare allele sharing, breakpoint validation), include them.

### Q: Why does the founder-pack analysis use indels and not SNPs?

`[CONFIRM:` this should be answered by reading MODULE_6's input parser. The memory context says "100-INDEL sliding windows" suggesting indels only, but MODULE_4A produces both. Possible reasons indels would be preferred: (a) indels mutate at slower rate per base than SNPs, so the same window covers more genealogical depth; (b) indels are less affected by purifying selection than coding-region SNPs; (c) phase information is more reliable for indels because each indel is a single-event marker. Need to verify against MODULE_6 source to give a definitive answer. `]`

### Q: How are Clair3 phase blocks related to inversion breakpoints?

Phase blocks are determined by WhatsHap's read-pair phasing — they end where no informative read spans consecutive variants. For inversion analysis, this means:

1. **Phase blocks rarely span an inversion breakpoint.** Reads bridging a true inversion junction are split-aligned (DELLY/Manta territory), not standard read pairs that WhatsHap uses for phasing. So phase blocks tend to terminate at or near inversion breakpoints.
2. **The pattern of phase-block termination can itself be an inversion signal.** Multiple phase blocks systematically ending at the same coordinate across samples (where they shouldn't, given coverage) is consistent with a structural variant disrupting normal read-pair geometry there.
3. **Founder-pack fragments are bounded by phase blocks.** This means a single "fragment" never crosses an inversion breakpoint. A founder pack defined within an inversion candidate region therefore characterizes haplotype identity within the inverted segment specifically.

### Q: What's the difference between MODULE_4A and Clair3's standard tutorial workflow?

Clair3's documented usage is: per-sample VCF discovery → (optionally) WhatsHap phasing → done. MODULE_4A adds the **five-phase post-processing pipeline** plus the **three-branch rescue system**, neither of which Clair3 ships with. These additions exist because at ~5× coverage Clair3's per-sample output has too many marginal-confidence calls to use directly for downstream analyses. The pipeline's job is to convert raw Clair3 + WhatsHap output into a stable, downstream-consumable catalog with rescue provenance tags. None of the additions modify Clair3's core algorithm; they all operate post-hoc on Clair3's output.

### Q: How is the cohort regenotyping (STEP_P07A) different from Clair3 discovery?

Clair3 discovery on a single sample asks: "given this BAM, where are the variants?" It applies confidence thresholds to the candidate variant sites it considers, and may discard marginal sites entirely.

Regenotyping in STEP_P07A asks: "given a known site list, what is the genotype of this sample at each site?" The site list is the population-aggregated catalog from STEP_P05–P06. For each site in that catalog, the regenotyper extracts the per-sample BAM evidence and forces a hard call (`0/0`, `0/1`, `1/1`, or `./.` if no evidence). This recovers calls that Clair3 was too uncertain to emit in the original discovery but that, given the population context, can be confidently genotyped per sample.

The result is a complete cohort-wide hard-genotype matrix at every population site — what downstream consumers need.

### Q: Can MODULE_4A output be used standalone without going through MODULE_5/6/CON?

Yes. The per-sample phased VCFs and the cohort shared catalog are stable VCF files conforming to the standard format. Any tool that reads VCF can consume them. The `INFO/RESCUE_BRANCH` and `INFO/REGEN_SOURCE` tags are MODULE_4A-specific additions; downstream tools will ignore them (no harm). For a one-off per-sample variant analysis, just point your tool at the per-sample VCF in `MODULE_4A/05_per_sample_final/sample_X.vcf.gz`.
