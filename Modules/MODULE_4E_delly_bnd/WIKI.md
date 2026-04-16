# MODULE_4E — DELLY Breakend (BND) Calling — Wiki

Detail-level reference for everything specific to DELLY's BND pipeline. Pairs with `README.md` (high-level + how to run). For cross-cutting content (how DELLY's `junction.h` classifies SVs in general, CT semantics, the orphan rescue framework that uses this module's output), see [`../SV_CALLING_FRAMEWORK_WIKI.md`](../SV_CALLING_FRAMEWORK_WIKI.md).

---

## Table of contents

1. [Why this module exists (orphan inversion rescue)](#1-why-this-module-exists-orphan-inversion-rescue)
2. [What a BND record actually is](#2-what-a-bnd-record-actually-is)
3. [How DELLY produces a BND](#3-how-delly-produces-a-bnd)
4. [The four CT classes and what each one means](#4-the-four-ct-classes-and-what-each-one-means)
5. [Filter thresholds (verified from code)](#5-filter-thresholds-verified-from-code)
6. [Reading a BND VCF record](#6-reading-a-bnd-vcf-record)
7. [The orphan inversion rescue pathway (downstream)](#7-the-orphan-inversion-rescue-pathway-downstream)
8. [Limitations specific to DELLY BND at ~5× coverage](#8-limitations-specific-to-delly-bnd-at-5-coverage)
9. [Manuscript-ready prose snippets](#9-manuscript-ready-prose-snippets)
10. [Q&A](#10-qa)

---

## 1. Why this module exists (orphan inversion rescue)

MODULE_4E is the **second-most-important MODULE_4 module for the inversion paper**, after MODULE_4D INV itself. Its role: catch the inversions that DELLY's INV typer missed.

The dynamic: DELLY's `selectInversions` function (in `junction.h`) classifies split-read pairs into two buckets — `3to3` left-breakpoint junctions and `5to5` right-breakpoint junctions (see [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2C](../SV_CALLING_FRAMEWORK_WIKI.md#2c--selectinversions--inv-3to3--5to5)). Then `delly call -t INV` runs a **pairing step** that tries to match each `3to3` to a `5to5` at the same locus to produce one `<INV>` VCF record. **When pairing fails**, the unpaired junction is demoted to a BND record carrying its original `CT` tag.

Pairing fails for two main reasons at ~5× coverage:

1. **Asymmetric breakpoint quality.** One breakpoint falls in clean unique sequence (strong split-read evidence), the other falls in a repeat (weak/noisy evidence). DELLY emits the strong side as a BND with `CT=3to3` or `CT=5to5`; the weak side may not even reach BND status.
2. **Distance exceeds DELLY's pairing window.** Very large inversions (>5 Mb in this pipeline's empirical experience) can have BND pairs that DELLY does not recognize as a single event.

Either way, the inversion is **not in the INV catalog** but its evidence is **in the BND catalog**. MODULE_5A2 STEP06 reads the BND catalog, extracts records with inversion-orientation `CT` tags, pairs them across both DELLY and Manta, and recovers these orphan inversions for the final candidate set.

Without MODULE_4E, every inversion missed by INV-typing is invisible to the rest of the pipeline. With MODULE_4E, we have a documented rescue pathway with a quantifiable contribution to the final inversion count.

---

## 2. What a BND record actually is

`BND` is shorthand for **"breakend"** — a single end of a structural rearrangement junction, recorded independently of whether its partner end was found. It is the most general SV record type in the VCF specification: any other SV (DEL, DUP, INV, INS, TRA) can in principle be represented as a pair of BND records.

DELLY uses BND records as the **catch-all for junctions it cannot confidently classify as a complete typed SV**. Three things end up in the BND catalog:

- **Inversion orphans** (`CT=3to3` or `CT=5to5` that did not pair into an INV — what we care about for the inversion paper)
- **Translocation breakends** (anything where `INFO/CHR2 ≠ CHROM` — these are the substrate for MODULE_4F TRA, see § 4)
- **Deletion- or duplication-orientation orphans** (`CT=3to5` or `CT=5to3` that did not pair into a DEL or DUP)

The same VCF record format is used for all three. The **`CT` tag** is what discriminates them (see § 4).

---

## 3. How DELLY produces a BND

DELLY does not have a dedicated `selectBreakends` function in `junction.h`. Instead, BND records are an **emission decision** made later in the workflow:

1. `junction.h` selection functions (`selectInversions`, `selectTranslocations`, etc.) cluster split-read pairs by orientation into per-type buckets.
2. The clustering step demands enough reads at a candidate junction to pass DELLY's internal evidence threshold.
3. The pairing step (in `delly call -t INV` for inversions, `delly call -t DEL` for deletions, etc.) tries to match complementary junctions into typed SV records.
4. **Junctions that pass clustering but fail pairing** are emitted as BND records carrying the `CT` tag from the original `junction.h` selection.

When you call `delly call -t BND` directly (which is what MODULE_4E does), you get **all** breakend records — the ones that would have been emitted by other typers as orphans, plus all inter-chromosomal events that have nowhere else to live. The `INFO/CT` field tells you which selection function produced each record.

### Per-sample evidence fields

DELLY BND records carry the same `FORMAT` fields as INV/DEL/DUP:
- `DR/DV` — discordant **R**eference / **V**ariant paired-end pairs
- `RR/RV` — junction (split-)read **R**eference / **V**ariant counts

Per-sample alt evidence = `DV + RV`. This is what MODULE_5A2 STEP02 re-extracts when validating BND-rescue candidates.

---

## 4. The four CT classes and what each one means

DELLY's `INFO/CT` field is the single most important annotation on a BND record. It encodes the **orientation** of the junction by reporting which DNA ends were joined.

| `CT` value | Junction type | Where it usually came from | Why it might appear in BND catalog |
|---|---|---|---|
| `3to3` | tail-to-tail, both reads right-clipped | INV left breakpoint | INV typer did not pair it with a 5to5 |
| `5to5` | head-to-head, both reads left-clipped | INV right breakpoint | INV typer did not pair it with a 3to3 |
| `3to5` | normal (deletion-type) | DEL | DEL typer did not produce a typed call |
| `5to3` | tandem-duplication-type | DUP | DUP typer did not produce a typed call |

**Important:** all four CT values can also appear on **inter-chromosomal** BNDs (where `CHR2 ≠ CHROM`). On inter-chromosomal records the four classes mean four different reciprocal-translocation orientations rather than the SV-type meanings above. See the `selectTranslocations` dissection in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 2E](../SV_CALLING_FRAMEWORK_WIKI.md#2e--selecttranslocations--bnd--tra) for the inter-chromosomal interpretation.

### How MODULE_4E and MODULE_4F divide the work

MODULE_4F TRA filters the BND catalog to keep **only** records where `CHR2 ≠ CHROM` (true inter-chromosomal translocations). MODULE_4E BND keeps **everything else** — the intra-chromosomal orphans that are the inversion-rescue substrate.

In practice this means the MODULE_4E catalog is dominated by:
- `CT=3to3` and `CT=5to5` orphan inversion junctions (relevant to the inversion paper)
- `CT=3to5` and `CT=5to3` orphans where DEL or DUP typing failed (less interesting; usually noise from repeat regions)

The orphan-inversion extraction in MODULE_5A2 STEP06 filters specifically:

```bash
bcftools view -i 'INFO/CT="3to3" || INFO/CT="5to5"' \
  catalog_226.BND.vcf.gz
```

The remaining DEL/DUP-orientation orphans (`CT=3to5` and `CT=5to3`) are kept in the catalog but are not currently consumed by any downstream module — they sit available if you ever want to do orphan-DEL or orphan-DUP rescue analogous to the inversion rescue.

---

## 5. Filter thresholds (verified from code)

The strict catalog filter is:

```bash
PASS + INFO/PRECISE=1 + QUAL ≥ 300 + INFO/PE ≥ 3
```

Source: `00_module4e_config.sh` lines 69–70 and `slurm/SLURM_A03_merge_genotype.sh` line 123:

```bash
STRICT_REQUIRE_PASS=1
STRICT_REQUIRE_PRECISE=1
STRICT_MIN_QUAL=300
STRICT_MIN_PE=3

# In SLURM_A03:
STRICT_EXPR='INFO/PRECISE=1 && QUAL>=300 && INFO/PE>=3'
```

### Why these specific values

Identical to MODULE_4D INV (300 / 3 / PRECISE), and for the same reasons:
- DELLY's QUAL ≥ 300 is the documented high-confidence threshold
- PE ≥ 3 is the inclusive cutoff that retains heterozygous singletons at ~5× coverage; downstream stringency lives in MODULE_5A2 STEP02/STEP03 evidence re-extraction
- PRECISE is required because IMPRECISE BNDs have CIPOS often >100 bp and cannot be reliably paired into orphan-INV rescue candidates within MODULE_5A2 STEP06's ±1 kb matching window

The rationale parallels [MODULE_4D § 3](../MODULE_4D_delly_inv/WIKI.md#3-filter-thresholds-verified-from-code) — see there for the full discussion of why not 4 or 5 PE.

### One implementation note worth flagging

Unlike INV/DEL/DUP, BND records do not carry `SVTYPE=DEL/DUP/INV/...` — they all carry `SVTYPE=BND`. The strict expression therefore does NOT include an `SVTYPE` clause. All filtering happens on QUAL, PE, and PRECISE. This is correct but means the BND catalog mixes inversion-orientation, DEL-orientation, DUP-orientation, and inter-chromosomal records in one VCF — they are distinguished only by `INFO/CT` and by whether `INFO/CHR2 == CHROM`.

---

## 6. Reading a BND VCF record

<!-- TODO: replace synthetic example with a real DELLY BND record from the production catalog -->

Synthetic example (intra-chromosomal `CT=3to3` orphan inversion left-breakpoint):

```
Scaffold05  3142000  BND00000412  G  G]Scaffold05:5891234]  340  PASS
PRECISE;SVTYPE=BND;SVMETHOD=EMBL.DELLYv1.7.3;CHR2=Scaffold05;
END=5891234;PE=4;MAPQ=42;CT=3to3;CIPOS=-8,8;CIEND=-12,12;
SRMAPQ=58;INSLEN=0;HOMLEN=2;SR=3;SRQ=0.94217;
CONSENSUS=AGCTGAATGCTAACGAATTTGCATCAGTCC...
GT:GL:GQ:FT:RR:RV:DR:DV
0/1:-12.45,-1.82,0:18:PASS:11:3:9:4
```

| Field | Value | Meaning |
|---|---|---|
| CHROM/POS | Scaffold05:3,142,000 | Left side of the breakend |
| ID | BND00000412 | DELLY-assigned record ID |
| ALT | `G]Scaffold05:5891234]` | Bracket pattern `]p]t` → `CT=3to3` orientation; mate at Scaffold05:5,891,234 |
| QUAL | 340 | Above 300 strict cutoff |
| `PRECISE` | flag | Split-read assembled |
| `SVTYPE=BND` | string | Always BND for this catalog |
| `CHR2=Scaffold05` | string | Mate chromosome — same as CHROM → **intra-chromosomal** |
| `END=5891234` | int | Mate position |
| `PE=4`, `SR=3` | int | 4 paired-end + 3 split-read supporting; both above strict cutoffs |
| `CT=3to3` | string | Tail-to-tail junction → **inversion left breakpoint orientation** |
| `CIPOS=-8,8`, `CIEND=-12,12` | range | Tight breakpoint confidence |
| `HOMLEN=2` | int | 2 bp microhomology (MMEJ-compatible) |
| `CONSENSUS=...` | string | Assembled junction sequence — can be BLAT'd |

Per-sample (one column shown):
- `GT=0/1` heterozygous
- `RR/RV = 11/3` — 3 alt-supporting junction reads
- `DR/DV = 9/4` — 4 alt-supporting discordant pairs
- Per-sample alt evidence = 3 + 4 = **7 reads**

### What this record means

This BND is a candidate **left breakpoint of an inversion** spanning Scaffold05:3,142,000–5,891,234 (~2.75 Mb). It got demoted from INV to BND because the matching `CT=5to5` right breakpoint was either too weak or too noisy for DELLY's INV-pairing logic to recognize as a partner.

If MODULE_5A2 STEP06 finds a `CT=5to5` BND on Scaffold05 within roughly the right distance window, the two will be paired into an orphan inversion candidate. If it doesn't find one in the DELLY catalog, the search extends to the Manta raw pre-conversion VCF for `t[p[` BNDs (Manta's `5to5` equivalent).

### Identifying inter-chromosomal BNDs

Same record format, but `INFO/CHR2 ≠ CHROM`. These are **NOT** in the MODULE_4E catalog — they were filtered out and went to MODULE_4F TRA. If you ever see a BND record in MODULE_4E with mismatched chromosomes, the filter regressed.

---

## 7. The orphan inversion rescue pathway (downstream)

The complete flow:

```
MODULE_4E BND catalog
  │
  ├─▶ MODULE_5A2 STEP06  orphan extraction
  │   bcftools view -i 'INFO/CT="3to3" || INFO/CT="5to5"' catalog_226.BND.vcf.gz
  │   → DELLY CT-based orphan list
  │
  ├─▶ Pool with Manta orphans
  │   MODULE_4G raw pre-conversion VCF (cohort_226.ALL.raw.vcf.gz)
  │   → ALT bracket parsing: ]p]t = 3to3, t[p[ = 5to5
  │   → Manta INV3/INV5 orphan list
  │   → Combined orphan pool across both callers
  │
  ├─▶ Greedy 5-Mb pairing
  │   For each unpaired CT=3to3 / INV3 left-junction:
  │     find nearest CT=5to5 / INV5 right-junction on same chr within 5 Mb
  │     each junction used at most once (greedy nearest-match)
  │   → Paired BND candidates table
  │
  ├─▶ Cross-reference with INV catalog (±1 kb)
  │   - If paired BNDs match a DELLY or Manta INV call → already known, skip
  │   - If paired BNDs match nothing → ORPHAN RESCUE CANDIDATE
  │
  └─▶ Orphan rescue list joins MODULE_5A candidate set
        Tagged class: bnd_rescue_delly / bnd_rescue_manta / orphan_bnd_pair
        Treated as lower-confidence than concordant calls but kept as
        candidates for breakpoint validation in MODULE_5A2 STEP02
```

The final classification (six classes total) is documented in [`../SV_CALLING_FRAMEWORK_WIKI.md` § 3](../SV_CALLING_FRAMEWORK_WIKI.md#3-what-each-callermodule-contributes-plug-map). Orphan-rescue candidates are explicitly marked as such in the candidate table; downstream consumers (manuscript figures, MODULE_5B followup analyses) can choose to include or exclude them by class.

---

## 8. Limitations specific to DELLY BND at ~5× coverage

1. **One-sided BNDs are common at low coverage.** Both junctions of an inversion need enough reads to pass the BND emission threshold. At ~5×, frequently only the "easier" junction (cleaner sequence context) makes it through. The other side is invisible in both INV and BND catalogs. These inversions cannot be rescued by MODULE_5A2 STEP06 — they can only be discovered by population-signal methods (MODULE_5A) and validated by re-extracting BAM evidence at the inferred breakpoint.

2. **Repeat-rich breakpoints fail twice.** A breakpoint falling in a repeat fails INV-typing (which is why it became a BND in the first place) AND often fails BND emission (because the split reads can't form a clean cluster). Inversions with both breakpoints in repeats are systematically under-detected by both INV and BND pathways.

3. **The 5 Mb pairing window is a heuristic.** MODULE_5A2 STEP06's greedy pairing assumes the inversion is ≤5 Mb. Larger inversions will have orphan BNDs that don't get paired and stay as singletons. There are no known catfish inversions >5 Mb in this cohort but the assumption is an artifact of the chosen window, not of the data.

4. **CT=3to5 and CT=5to3 orphans are kept but unused.** These could in principle be DEL or DUP rescue candidates analogous to inversion rescue. No downstream module currently consumes them. If a future analysis cares about hidden DELs or DUPs, the data is already in the catalog.

5. **`delly filter -f germline` removes some real BNDs.** The germline filter applies allele-frequency / Hardy-Weinberg-like criteria that can be inappropriate for true singletons or rare-frequency private inversions. If a candidate disappears at this stage, the pre-filter BCF is the audit path; downstream modules sometimes deliberately read from the pre-filter BCF for sensitivity analyses.

---

## 9. Manuscript-ready prose snippets

### Methods — DELLY BND calling paragraph

> DELLY2 v1.7.3 was used to discover all breakend (BND) records across the 226 samples by running `delly call -t BND`. BND records carry `INFO/CT` (connection type) tags that encode the orientation of the breakpoint junction: `CT=3to3` and `CT=5to5` correspond to the two breakpoint orientations of an inversion (left and right respectively), while `CT=3to5` and `CT=5to3` correspond to deletion- and duplication-type junctions. BND records are emitted by DELLY when junction clustering succeeds but typed-SV pairing fails — i.e., when a single inversion breakpoint passes evidence thresholds but its complementary breakpoint does not. Per-sample discovery, cohort merging, regenotyping, subsetting to 81 unrelated samples, and germline filtering followed the same workflow as MODULE_4D INV. The strict catalog required `FILTER = PASS`, `INFO/PRECISE = 1`, `QUAL ≥ 300`, and `INFO/PE ≥ 3`. Inter-chromosomal breakends (`INFO/CHR2 ≠ CHROM`) were routed to MODULE_4F TRA; intra-chromosomal records remained in the BND catalog and served as the substrate for the orphan-inversion rescue pipeline (MODULE_5A2 STEP06).

### Results — number reporting placeholder

> The MODULE_4E BND catalog contained [TODO N_BND_STRICT] records after strict filtering, of which [TODO N_3TO3] carried `CT=3to3` and [TODO N_5TO5] carried `CT=5to5` (the two inversion-orientation classes). MODULE_5A2 STEP06 paired these with the corresponding Manta orphan BNDs (parsed from the raw pre-conversion VCF by ALT bracket pattern) and identified [TODO N_RESCUED] candidate inversions that had not been emitted as `<INV>` by either DELLY's or Manta's INV typer. Of these orphan-rescue candidates, [TODO N_VALIDATED] passed downstream BAM-evidence re-extraction (MODULE_5A2 STEP02) and statistical seed qualification (STEP03). The orphan rescue pathway therefore contributed [TODO PCT]% of the final inversion candidate set.

### Discussion — limitations footnote

> The orphan-rescue pathway is asymmetric: it can recover inversions where one breakpoint is detectable but the other is not in the same caller's INV-typing output. It cannot recover inversions where both breakpoints fall in unscorable regions (repeats, low-mappability, assembly gaps). This bias is partially mitigated by the dual-caller design (DELLY + Manta in MODULE_4G), since the two callers have different sensitivity profiles to repeat context, but inversions with both breakpoints in highly repetitive sequence remain at the limit of detection from short-read data alone.

---

## 10. Q&A

### Q: What is the practical difference between an INV record and a `CT=3to3 + CT=5to5` BND pair?

The INV record is what you get when DELLY's pairing step succeeds. The BND pair is what you get when pairing fails. Both represent the same biological event (an inversion). The downstream pipeline treats them differently for confidence-class purposes:

- `<INV>` records → directly enter MODULE_5A2 STEP01 candidate extraction
- BND pairs → orphan rescue in STEP06, then enter STEP01 with rescue-class tags

Orphan-rescue candidates are weighted as lower-confidence in downstream prose (manuscript reporting separates the two classes), but if they pass STEP02 BAM re-extraction and STEP03 statistical tests, they are real inversions on equal footing with the others.

### Q: Why does MODULE_4E keep `CT=3to5` and `CT=5to3` records if nothing uses them?

Disk is cheap, and the upstream extraction (`bcftools view -i 'INFO/SVTYPE="BND"'`) catches them by default. Filtering them out at MODULE_4E would lose information that might be useful for hypothetical future analyses (orphan DEL or DUP rescue). The MODULE_5A2 STEP06 inversion-rescue extraction then explicitly filters for `CT=3to3 || CT=5to5`, ignoring the others. Two-stage retention is the standard pattern: keep everything at the catalog level, filter at the consumer level.

### Q: An inversion is in MODULE_5A's population candidates but NEITHER MODULE_4D nor MODULE_4E shows anything at those coordinates. Real?

Probably not, but proceed with care. Three real possibilities to check:

1. **Both breakpoints in highly repetitive sequence** — split reads can't form clusters. Check whether the breakpoints overlap repeat annotation (RepeatMasker, BISER2 SDs).
2. **Inversion is too large for DELLY's pairing** — check whether the population-signal region span is >5 Mb (DELLY's heuristic limit).
3. **Population signal is not actually an inversion** — could be a recombination-suppressed region from balancing selection, a sweep, or a copy-number artifact. Compare against the heterozygosity track and against MODULE_4B DEL / MODULE_4C DUP catalogs at the same coordinates.

If all three are excluded, the candidate is a high-priority target for long-read validation rather than a confident inversion call.

### Q: How do I tell whether a BND record came from `selectInversions` or `selectTranslocations`?

Two checks:

1. **Same-chromosome test**: `INFO/CHR2 == CHROM`. If yes → intra-chromosomal, came from `selectInversions` (CT=3to3 or 5to5) or from DEL/DUP typing failure (CT=3to5 or 5to3). If no → inter-chromosomal, came from `selectTranslocations`. (Inter-chromosomal records should be in MODULE_4F TRA, not MODULE_4E.)
2. **CT class for intra-chromosomal**: `CT=3to3` or `CT=5to5` → from `selectInversions`. `CT=3to5` or `CT=5to3` → from DEL or DUP typing failure.

There is no `INFO` field that explicitly names the source function. The combination of `CHR2 == CHROM` plus `CT` value is sufficient for full classification.

### Q: I see `IMPRECISE` BNDs in the catalog. Should I filter them out?

The strict catalog already requires PRECISE, so IMPRECISE BNDs are not in `catalog_226.BND.PASS.vcf.gz`. If you see IMPRECISE BNDs in your output, you are reading from the **pre-strict-filter** BCF (e.g., `cohort_226.merged.vcf.gz`). For orphan rescue specifically, IMPRECISE BNDs are problematic because their CIPOS is typically wide (±100–500 bp), which makes the ±1 kb matching window in MODULE_5A2 STEP06 unreliable — you can match almost anything. Stick to PRECISE for rescue.

### Q: A `CT=3to3` BND has no nearby `CT=5to5` partner. Is the inversion lost?

Yes, in the current pipeline. The orphan rescue requires both junctions; a single-sided orphan cannot be paired. This is a limitation, not a bug. The rescue candidate set is therefore conservatively biased toward inversions where both breakpoints have at least minimal evidence. Single-sided orphans accumulate in the BND catalog but contribute no rescue candidates. They could in principle be cross-referenced with the population-signal regions (MODULE_5A) to ask "does any population candidate have a single-sided BND at one of its breakpoints", but no module currently does this. If this is biologically interesting, it is a possible future analysis.
