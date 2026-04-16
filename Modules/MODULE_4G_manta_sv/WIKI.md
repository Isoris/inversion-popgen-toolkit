# MODULE_4G — Manta SV Calling — Wiki

Detail-level reference for everything specific to Manta's pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (DELLY×Manta concordance, BND orphan rescue, population binomial annotation), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module exists (parallel SV evidence)](#1-why-this-module-exists-parallel-sv-evidence)
2. [How Manta produces SV records (assembly-based)](#2-how-manta-produces-sv-records-assembly-based)
3. [The convertInversion problem and per-sample workflow](#3-the-convertinversion-problem-and-per-sample-workflow)
4. [The six output catalogs (DEL, DUP, INV, BND, INS_small, INS_large)](#4-the-six-output-catalogs)
5. [Filter thresholds (verified from code)](#5-filter-thresholds-verified-from-code)
6. [Reading a Manta VCF record](#6-reading-a-manta-vcf-record)
7. [What the output feeds downstream](#7-what-the-output-feeds-downstream)
8. [Limitations specific to Manta at ~5× coverage](#8-limitations-specific-to-manta-at-5-coverage)
9. [Manuscript-ready prose snippets](#9-manuscript-ready-prose-snippets)
10. [Q&A](#10-qa)

---

## 1. Why this module exists (parallel SV evidence)

Manta is the **independent second SV caller** for every type that DELLY produces. The motivation is methodological: DELLY and Manta use fundamentally different algorithms, so concordance between them is strong evidence that a call is real, while caller-specific calls are a known signal of caller bias.

- DELLY is a **discordant-read clusterer**. It groups paired-end and split-read evidence by orientation (`junction.h` selection logic) and produces SV records from clusters that pass evidence-count thresholds.
- Manta is an **assembly-based caller**. It builds a graph of breakend associations from anomalous reads, then performs **local de-novo assembly** of the candidate junction sequence, scoring each assembly with a probabilistic diploid likelihood model.

Three concrete consequences for the inversion paper:

1. **Concordance validation.** Every DELLY INV call is cross-referenced against Manta INV (MODULE_5A2 STEP05). Concordant calls (50% reciprocal overlap) form the highest-confidence class.
2. **Orphan rescue across callers.** Manta's BND records (post-conversion catalog contains only inter-chromosomal translocations; intrachromosomal inversion-orientation BNDs were converted to `<INV>`) plus DELLY's BND records pool into a combined orphan pool. MODULE_5A2 STEP06 pairs `INV3+INV5` from Manta and `CT=3to3+CT=5to5` from DELLY across both callers.
3. **Sole insertion source.** DELLY INS was dropped from the production pipeline because the false-negative rate at ~5× coverage was unacceptable. Manta carries the entire insertion catalog (`INS_small` + `INS_large` — see § 4).

---

## 2. How Manta produces SV records (assembly-based)

The canonical reference is Chen *et al.* 2016 (Bioinformatics). The summary as it applies to this pipeline:

### Step 1 — Breakend graph construction

Manta scans each BAM for **anomalous read pairs** (unexpected insert size, unexpected orientation) and **soft-clipped reads**. Each anomaly contributes an edge to a graph where nodes are genomic intervals. Connected components in this graph become **candidate SV regions** for assembly.

Configurable thresholds (from `configManta.custom.ini`):
- `minEdgeObservations = 3` — at least 3 reads required to form a graph edge
- `minCandidateSpanningCount = 3` — at least 3 spanning reads to consider a candidate
- `minCandidateVariantSize = 50` — only consider SVs ≥50 bp (raised from default 8 to avoid redundancy with Clair3 small-indel pipeline; see `configManta.custom.ini` line 4)

### Step 2 — Local de-novo assembly

For each candidate SV, Manta extracts all reads supporting both sides of the breakend graph and runs a **local de Bruijn assembly** of the junction. The output is either:

- **A complete assembled contig** spanning the breakpoint → emits as `PRECISE` with full ALT sequence (for DEL/DUP/INV) or a complete `SVLEN` (for INS_small)
- **A flanking-only assembly** where left and right sides assembled but didn't merge → emits as `IMPRECISE` for DEL/DUP/INV, or as `INS_large` for insertions (with `LEFT_SVINSSEQ` / `RIGHT_SVINSSEQ` flanking sequences but no central content)

This is the crucial difference from DELLY: Manta's `PRECISE`/`IMPRECISE` flag specifically tracks whether the **assembly succeeded**, not whether split reads exist. A Manta `PRECISE` call has a real reconstructed junction sequence; `IMPRECISE` means evidence exists but the junction didn't assemble cleanly.

### Step 3 — Diploid likelihood scoring

Manta scores each assembled SV with a likelihood model that compares evidence under three genotype hypotheses (`0/0`, `0/1`, `1/1`) and emits a **per-sample QUAL** combined into a record-level QUAL. The minimum score for diploid output is set in `configManta.custom.ini`:

- `minScoredVariantSize = 50` — same as candidate threshold
- `minDiploidVariantScore = 10` — record QUAL floor before VCF output

### Per-sample evidence fields

Manta emits **per-sample** `FORMAT/PR` (paired-read, ref:alt) and `FORMAT/SR` (split-read, ref:alt) — distinct from DELLY which puts global counts in INFO and per-sample in DR/DV/RR/RV.

```
FORMAT  GT:FT:GQ:PL:PR:SR
sample  0/1:PASS:240:283,0,999:8,3:12,5
                         ^^^ ^^^
                         |   split-read 12 ref + 5 alt
                         paired-read 8 ref + 3 alt
```

Per-sample alt evidence = `PR_alt + SR_alt` = `3 + 5` = 8 reads.

---

## 3. The convertInversion problem and per-sample workflow

Manta represents inversions natively as **two BND records sharing an `EVENT` tag**, with `INFO/INV3` (left breakpoint, 3to3 equivalent) and `INFO/INV5` (right breakpoint, 5to5 equivalent) flags. The standard tool to merge these into single `<INV>` VCF records is `convertInversion.py` shipped with Manta.

**Critical operational detail:** `convertInversion.py` must run **per-sample**, **before** cohort merging — never on the merged VCF. The reason is a TypeError crash documented at the top of `slurm/SLURM_A03_merge_and_split.sh`:

> After `bcftools merge`, the VCF ID column contains compound IDs (semicolon-separated from multiple samples), e.g., `MantaBND:148628:0:1:0:0:0:0;MantaBND:59701:0:1:0:0:0:1;...`. `convertInversion.py`'s `scanVcf()` uses `vcfRec.vid` to build the `invMateDict`, and later does `invMateDict[mateId]["CIPOS"]`. With compound IDs, the mate lookup returns an empty string `""` instead of a dict → TypeError: `string indices must be integers, not str`.

The solution wired into MODULE_4G is the standard nf-core / CGAP workflow:

```
For each sample:
   diploidSV.vcf.gz  →  convertInversion.py  →  diploidSV.converted.vcf.gz
After all samples converted:
   bcftools merge → cohort_226.merged.vcf.gz
Subset to unrelated:
   bcftools view -S samples_unrelated_81.txt → cohort_81.merged.vcf.gz
Split by SVTYPE → six per-type catalogs
```

The per-sample converted VCFs live alongside the raw outputs; `utils/convertInversion_py3.py` is the Python-3-adapted version of the script (Manta ships only Python 2).

### Consequence for orphan rescue

After per-sample conversion, the **post-conversion BND catalog** (`catalog_226.BND.PASS.vcf.gz`) contains **only inter-chromosomal translocations** — every intrachromosomal `INV3+INV5` pair was already merged to `<INV>`. Inversions with only one detected junction (where `INV3` and `INV5` did NOT pair successfully within `convertInversion.py`'s logic) are **lost from the BND catalog** because they were dropped during conversion.

To recover these orphan junctions, MODULE_5A2 STEP06 reads from the **raw pre-conversion merged VCF** (`cohort_226.ALL.raw.vcf.gz`) and classifies orientation from the BND ALT field bracket pattern instead of from `INV3`/`INV5` flags:

| ALT pattern | Orientation |
|---|---|
| `]p]t` | `3to3` (INV3 equivalent) |
| `t[p[` | `5to5` (INV5 equivalent) |

These pre-conversion orphan BNDs are then pooled with DELLY's `CT=3to3` and `CT=5to5` orphans, paired greedily within 5 Mb on the same chromosome, and cross-referenced against the full INV catalog (±1 kb) to identify rescued inversions both INV-typers missed.

---

## 4. The six output catalogs

Manta does not have separate runs per SV type the way DELLY does (`-t DEL`, `-t INV`, etc.). One Manta run produces **one diploid VCF per sample** containing all SV types. After cohort merging and convertInversion, the cohort VCF is split into six per-type catalogs by `slurm/SLURM_A03_merge_and_split.sh`:

| Catalog | What it contains | Notes |
|---|---|---|
| `catalog_226.DEL.PASS.vcf.gz` | All `SVTYPE=DEL` records | Direct from Manta |
| `catalog_226.DUP.PASS.vcf.gz` | All `SVTYPE=DUP` records (including `SVTYPE=DUP:TANDEM`) | Direct from Manta |
| `catalog_226.INV.PASS.vcf.gz` | All `<INV>` records produced by per-sample convertInversion | Merged from INV3+INV5 BND pairs |
| `catalog_226.BND.PASS.vcf.gz` | **Only inter-chromosomal translocations** (post-conversion) | Intra-chromosomal BNDs already became INV |
| `catalog_226.INS_small.PASS.vcf.gz` | Insertions with **complete** ALT assembly (`SVLEN` present, `SVINSSEQ` populated) | Typically 50–200 bp |
| `catalog_226.INS_large.PASS.vcf.gz` | Insertions with **flanking-only** assembly (`LEFT_SVINSSEQ`+`RIGHT_SVINSSEQ`, no central sequence) | Repetitive or large insertions where assembly stalled |

The INS_small / INS_large split is a Manta-specific feature with no DELLY equivalent. It exists because `INS_large` calls have known coordinates and known flanks but unknown insertion content — useful for "there is something inserted here" analysis but unsuitable for any analysis requiring the inserted sequence (e.g., TE family classification, gene-disruption prediction).

---

## 5. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
FILTER = PASS + QUAL ≥ 20
```

Source: `00_module4g_config.sh` lines 99–100 and `slurm/SLURM_A03_merge_and_split.sh` lines 284, 287:

```bash
STRICT_REQUIRE_PASS=1
STRICT_MIN_QUAL=20

# In SLURM_A03:
bcftools view -f PASS -i "QUAL>=${STRICT_MIN_QUAL}" ...
```

### Why these specific values

| Threshold | Value | Why this value |
|---|---|---|
| QUAL | ≥ 20 | Manta's QUAL is a phred-scaled likelihood from the diploid scoring model, **not comparable to DELLY's QUAL**. QUAL ≥ 20 corresponds to ~99% posterior confidence under Manta's model. The Manta authors recommend QUAL ≥ 20 for diploid germline calls. |
| FILTER = PASS | required | Manta's per-record FILTER includes failure flags like `MinQUAL`, `MinGQ`, `MaxDepth`, `MaxMQ0Frac`, `Ploidy`, `MinSomaticScore` (somatic mode only). PASS means none of these triggered. |
| Per-sample evidence | none at strict layer | Per-sample `PR_alt + SR_alt` thresholds are applied downstream in MODULE_5A2 STEP01 (≥6 total alt evidence for INV candidates), not at the cohort filter. |

### Why no PRECISE requirement

Manta's `IMPRECISE` flag means **assembly failed**, not that no evidence exists. Filtering all IMPRECISE calls would discard:

- Most large inversions (≥1 Mb) where neither breakpoint assembles cleanly
- Most calls in repeat-rich regions
- A substantial fraction of true heterozygous calls with insufficient coverage to support assembly

At ~5× coverage this would be a very aggressive filter. MODULE_4G keeps IMPRECISE records and uses them; downstream consumers (MODULE_5A2) use CIPOS to widen evidence-extraction windows for IMPRECISE calls rather than filtering them out.

### Why no PE / SR threshold at the cohort layer

Manta's QUAL already incorporates PR and SR evidence weights into the likelihood — adding a separate evidence threshold would double-count. The only place per-sample evidence is filtered is at the candidate-extraction step (MODULE_5A2 STEP01: PR_alt + SR_alt ≥ 6 for INV candidates), which is a candidate-quality filter rather than a discovery filter.

---

## 6. Reading a Manta VCF record

Two example records — one INV (post-convertInversion) and one BND (post-conversion, inter-chromosomal).

### Example 6A — A converted INV record

```
Scaffold03  7600386  MantaINV:142:1:0:0:0:0  T  <INV>  430  PASS
END=9802651;SVTYPE=INV;SVLEN=2202265;CIPOS=-12,12;CIEND=-15,15;
EVENT=MantaINV:142:1:0:0:0;HOMLEN=4;HOMSEQ=GTAC;
SVINSLEN=0;BND_DEPTH=8;MATE_BND_DEPTH=7
GT:FT:GQ:PL:PR:SR
0/1:PASS:240:283,0,999:8,3:12,5
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS/END | Scaffold03:7600386–9802651 | Inversion span |
| ID | `MantaINV:142:1:0:0:0:0` | Manta-assigned ID (post-convertInversion format) |
| ALT | `<INV>` | Symbolic — produced by convertInversion from INV3+INV5 BND pair |
| QUAL | 430 | Diploid likelihood score (well above 20 cutoff) |
| `SVLEN=2202265` | int | 2.2 Mb inversion |
| `CIPOS=-12,12` / `CIEND=-15,15` | range | Tight breakpoint confidence (assembled junctions) |
| `EVENT=MantaINV:142:1:0:0:0` | string | Original BND pair event ID before conversion |
| `HOMLEN=4`, `HOMSEQ=GTAC` | int + string | 4 bp microhomology at junction (MMEJ-compatible) |
| `BND_DEPTH=8`, `MATE_BND_DEPTH=7` | int | Read depth at the two breakpoints |

Per-sample (one column shown):
- `GT=0/1` — heterozygous inversion
- `PR=8,3` — 8 ref-supporting paired reads, 3 alt-supporting
- `SR=12,5` — 12 ref-supporting split reads, 5 alt-supporting
- Per-sample alt evidence = 3 + 5 = **8 reads**

### Example 6B — A post-conversion BND (inter-chromosomal translocation)

```
Scaffold03  4521000  MantaBND:891:0:1:0:0:0:0  T  T]Scaffold17:11203450]  120  PASS
SVTYPE=BND;MATEID=MantaBND:891:0:1:0:0:0:1;CIPOS=-50,50;
BND_DEPTH=6;MATE_BND_DEPTH=5
GT:FT:GQ:PL:PR:SR
0/1:PASS:80:120,0,500:4,2:0,0
```

| Field | Value | Meaning |
|---|---|---|
| ALT | `T]Scaffold17:11203450]` | Bracket pattern `]p]t` = 3to3 orientation, mate at Scaffold17:11203450 |
| `MATEID` | `MantaBND:891:...:1` | The other half of this translocation's BND pair |
| CIPOS=-50,50 | range | Wider than INV example (no split-read assembly here) |
| PR=4,2 | counts | 2 alt paired reads (no split-read evidence — note SR=0,0) |

### The four ALT bracket patterns (for orphan rescue)

| ALT pattern | Mate orientation | DELLY CT equivalent |
|---|---|---|
| `]p]t` | mate position with `]` brackets, base after | `3to3` |
| `t[p[` | base before, mate position with `[` brackets | `5to5` |
| `t]p]` | base before, mate position with `]` brackets | `3to5` (deletion-type) |
| `[p[t` | mate position with `[` brackets, base after | `5to3` (duplication-type) |

MODULE_5A2 STEP06's orphan rescue parses these patterns from the **raw pre-conversion VCF** because, after convertInversion, all `]p]t` and `t[p[` intrachromosomal BNDs have been merged to `<INV>` and are no longer in the BND catalog.

---

## 7. What the output feeds downstream

```
MODULE_4G Manta catalogs (six)
  │
  ├─▶ MODULE_5A2 STEP01  candidate extraction (Manta INV branch)
  │   - QUAL ≥ 50, total alt evidence (PR_alt + SR_alt) ≥ 6
  │   - Size window 5 kb ≤ SVLEN ≤ 20 Mb
  │   - Output: unified candidate table with caller='manta'
  │
  ├─▶ MODULE_5A2 STEP05  DELLY×Manta concordance
  │   - 50% reciprocal overlap test against MODULE_4D DELLY INV
  │   - Six classes (concordant, manta_only, delly_only, two BND-rescue classes, orphan)
  │
  ├─▶ MODULE_5A2 STEP06  orphan BND rescue
  │   - Reads pre-conversion VCF for INV3/INV5 orphans
  │   - Pools with DELLY CT=3to3/5to5 BNDs across both callers
  │   - Pairs within 5 Mb, cross-references against INV catalog
  │
  ├─▶ MODULE_4G INS_small / INS_large
  │   - Sole insertion source (DELLY INS dropped)
  │   - INS_small: complete ALT, used for TE family classification
  │   - INS_large: flanks only, used for "insertion exists here" analysis
  │
  └─▶ Population binomial annotation
      - Applied to every PASS record across all six catalogs
      - See ../SV_CALLING_FRAMEWORK_WIKI.md
```

For the full cross-module plug-map (cheats, framework methods), see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map).

---

## 8. Limitations specific to Manta at ~5× coverage

1. **Assembly success rate drops sharply below 10× coverage.** At ~5×, a substantial fraction of Manta candidates are emitted as `IMPRECISE` because the local de-novo assembly couldn't complete. These calls have correct *evidence* but uncertain *coordinates* (CIPOS often ±100–300 bp). Strict downstream analyses must either tolerate wide CIPOS or supplement Manta IMPRECISE calls with DELLY split-read coordinates where concordant.

2. **`INS_large` insertions cannot be classified by content.** Without the central insertion sequence, you cannot run RepeatMasker / TEclass / SnpEff on these — they exist as positions, not as sequence variants. Useful for counting "how many insertion events", useless for "what kind of insertion."

3. **Manta DUP can split a single tandem dup into multiple records.** When a tandem duplication has multiple breakpoints in the assembly graph (common with imperfect repeats), Manta may emit 2–3 overlapping DUP records for what is biologically one event. MODULE_4G does not currently de-overlap these; downstream consumers must apply 50% reciprocal overlap merging.

4. **convertInversion fails on compound IDs.** Documented at length in § 3. The wired-in workaround (per-sample conversion before merge) is correct but means cohort-merge-then-convert workflows from other Manta tutorials will silently fail or produce wrong results. Worth flagging when adapting code from other repos.

5. **Manta translocation calls require both BNDs to exist.** A single-side-only translocation is not callable — Manta requires both partners of a reciprocal break to assemble. DELLY can call one-sided BNDs (orphans). Combined with DELLY in MODULE_5A2, this asymmetry actually helps: Manta provides confirmed reciprocals, DELLY provides single-sided rescue evidence.

6. **No cohort-level discovery.** Each Manta run sees one BAM at a time; cohort-level information (joint PR/SR distributions, joint assembly) is not exploited. DELLY's cohort merge does not help here either — Manta needs one-shot joint calling for full cohort sensitivity, which it does not support. Both callers underestimate population frequency for inversions present in many samples; the binomial annotation is the partial mitigation.

---

## 9. Manuscript-ready prose snippets

### Methods — Manta calling paragraph

> Manta v1.6.0 (Chen *et al.*, 2016) was run as an independent SV caller alongside DELLY2 on all 226 samples. Per-sample `configManta.py` was invoked with `minCandidateVariantSize = 50` (raised from default 8 to delineate the responsibility split with the Clair3 small-variant pipeline; SVs <50 bp were left to Clair3). Each sample's `runWorkflow.py` produced a per-sample `diploidSV.vcf.gz`. Inversions, which Manta natively represents as paired BND records with `INV3` / `INV5` orientation flags, were merged to `<INV>` records using a Python 3-adapted version of `convertInversion.py` **per-sample** before cohort merging — necessary because the script's BND mate-resolution logic fails on the compound IDs produced by `bcftools merge`. The 226 per-sample converted VCFs were merged with `bcftools merge`, subset to 81 unrelated individuals, and split into six per-type catalogs (DEL, DUP, INV, BND, INS_small, INS_large) where INS_small contained insertions with complete ALT assembly and INS_large contained insertions with flanking-only assembly. The strict catalog required `FILTER = PASS` and `QUAL ≥ 20`. Per-sample evidence used `FORMAT/PR` (paired-read ref:alt) and `FORMAT/SR` (split-read ref:alt) summed across carriers.

### Results — number reporting placeholder

> Manta produced [TODO N_RAW_MANTA] candidate SVs across all six types in the cohort merge. After PASS + QUAL ≥ 20 filtering, [TODO N_STRICT_MANTA] survived. Among inversions specifically, [TODO N_MANTA_INV] passed strict filtering. Cross-comparison with the DELLY INV catalog (50% reciprocal overlap, MODULE_5A2 STEP05) classified [TODO N_CONC] as concordant, [TODO N_MO] as Manta-only, and [TODO N_DO] as DELLY-only. Insertion calls [TODO N_INS_SMALL] (complete assembly) + [TODO N_INS_LARGE] (flanking-only) form the entire insertion catalog, as DELLY INS calling was excluded from the production pipeline due to ~5× coverage limitations (Methods).

### Discussion — limitations footnote

> Manta's local de-novo assembly approach is more sensitive than DELLY's clustering for complex junctions but more dependent on coverage for assembly success. At ~5× depth, a substantial fraction of Manta calls emit as IMPRECISE (assembly incomplete) — these retain correct evidence-based detection but with wider breakpoint coordinates. The dual-caller design (DELLY + Manta concordance, MODULE_5A2) was specifically motivated by this asymmetry: DELLY provides tighter coordinates when split-read clustering succeeds, Manta provides assembly-confirmed junction sequences when assembly succeeds. The two callers' systematic biases are partially uncorrelated, making concordance a stronger signal than either alone.

---

## 10. Q&A

### Q: Why was DELLY INS dropped but Manta INS kept?

DELLY's INS detection (`junction.h::selectInsertions` + `bridgeInsertions`) clusters discordant pairs where the reference footprint is small but the sequence footprint is large. At ~5× coverage, these clusters rarely accumulate enough reads to pass DELLY's internal thresholds — false-negative rate was unacceptable in pipeline validation. Manta builds insertion calls from local assembly graphs, which can produce a call from as few as 2 anomalous reads if they assemble cleanly. INS_large catches insertions where assembly only resolved the flanks; this class has no DELLY equivalent.

### Q: What does it mean if `convertInversion.py` "fails silently"?

If you run `convertInversion.py` on a merged-cohort VCF with compound IDs (`MantaBND:148628:0:1:0:0:0:0;MantaBND:59701:...`), one of two things happens:
1. Python crashes with `TypeError: string indices must be integers, not str` — this is the loud failure
2. In some script versions, the script silently skips records it can't resolve, producing an output VCF that **looks correct** but is missing inversions

The wired-in workflow (per-sample conversion before merge) avoids both. If you ever debug a low Manta INV count: check whether `convertInversion.py` was run pre-merge or post-merge. Pre-merge is correct; post-merge is the bug.

### Q: Why is QUAL = 20 enough for Manta but DELLY needs QUAL ≥ 300?

They are **different scoring systems** and the numbers are not comparable. DELLY's QUAL combines paired-end count, split-read count, and mapping quality into a non-phred-scaled aggregate score where higher is better but the units are arbitrary. Manta's QUAL is a **phred-scaled likelihood** from the diploid scoring model: QUAL ≥ 20 means posterior probability ≥ 99% of being a real variant under Manta's model. Don't try to convert between the two. The threshold for each caller is calibrated independently to that caller's score distribution.

### Q: My MODULE_4G INV catalog is missing some inversions that exist in the raw per-sample VCFs. Why?

Three places this can happen:

1. **convertInversion didn't pair them.** If the matching `INV3` and `INV5` BNDs are too far apart, on different scaffolds, or one was filtered out, convertInversion can't merge them. They stay as BNDs in the per-sample VCF but get dropped from the post-conversion BND catalog (which only keeps inter-chromosomal). Fix: orphan rescue (MODULE_5A2 STEP06).
2. **bcftools merge dropped them.** If a record has a per-sample `FT` (filter) that fails the merge logic, it can be excluded. Check `bcftools view -h` of the merged VCF for filter rules.
3. **Strict filter dropped them.** PASS + QUAL ≥ 20 is lenient but not zero. Check the raw merged VCF (`cohort_226.ALL.raw.vcf.gz`) and grep for the missing INV ID to see what FILTER/QUAL values it had.

### Q: What is the difference between INS_small and INS_large?

Both are insertions Manta detected but they differ in assembly completeness:

- **INS_small** has a complete `SVLEN` field and a populated ALT sequence (`SVINSSEQ` or the ALT column itself for short insertions). The full insertion content was assembled. Useful for any downstream analysis requiring the inserted sequence (TE family, gene impact, PCR primer design).
- **INS_large** has `LEFT_SVINSSEQ` and `RIGHT_SVINSSEQ` (the flanking sequences either side of the insertion) but no central content. Manta knows where the insertion is and what the ends look like, but can't tell you what's in the middle. Common for large repetitive insertions (full-length TE insertions, segmental duplications). Useful for counting, useless for content classification.

The split is not arbitrary — it's a real distinction that affects what downstream questions you can ask. Manuscript should report both counts separately and clearly state that INS_large is content-blind.

### Q: Manta's BND catalog only has translocations. Where did the inversion BNDs go?

They were converted to `<INV>` records by per-sample `convertInversion.py` (see § 3). Anything that was an `INV3`+`INV5` BND pair before conversion is now in the INV catalog, not the BND catalog. The post-conversion BND catalog (`catalog_226.BND.PASS.vcf.gz`) contains **only inter-chromosomal translocations** by construction. To find inversion-orientation BNDs that didn't successfully convert (orphans), read from the raw pre-conversion VCF (`cohort_226.ALL.raw.vcf.gz`) and parse the ALT bracket patterns directly — see § 6 for the four bracket patterns. This is what MODULE_5A2 STEP06 does.

### Q: Should I trust Manta-only INV calls (no DELLY support)?

Conditionally yes. Manta uses an entirely different algorithm so Manta-only calls reflect Manta's assembly-based sensitivity to junctions DELLY's clustering missed. But the false-positive rate for Manta-only calls is also higher than concordant calls. The principled handling is in MODULE_5A2 STEP02 — re-extract per-sample BAM evidence for every Manta-only call and let STEP03's statistical tests decide which Manta-only candidates have population-level support strong enough to qualify as seeds. Manta-only calls that pass STEP03 are real; Manta-only calls that fail STEP03 are likely caller-specific artifacts.
