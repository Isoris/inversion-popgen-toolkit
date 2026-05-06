# S5 — SV interpretation rules (discipline doc)

**Status:** spec / discipline doc.
**Used by:** P5.1 Section 6 (SV evidence), P5.1 Section 9 (group
association), all SV-related copy in tooltips and the manuscript.

This is not a code spec. It's the set of rules that govern what
the atlas can claim and what it cannot, when SV records appear
alongside boundary evidence.

---

## The five hard rules

### Rule 1 — INV ≠ raw read evidence

> INV is a VCF-level event/call. Raw evidence is in BAM
> alignments / read pairs / split reads.

When the atlas displays an INV record, the chip / label / tooltip
must make clear it's a **caller's claim**, not raw read counts.
INV records can be wrong. They can be inferred from spurious
discordant pairs. They can collapse a real complex rearrangement
into "INV" because that's the closest VCF type the caller knew.

The atlas chip for an INV looks like:

```
INV · LG28:14823100-17712450 · manta · QUAL=247
```

Not:

```
INV at LG28:14823100-17712450 (breakpoint)
```

Difference: the second implies the breakpoint is at those exact
coordinates; the first reports the call without endorsing it.

### Rule 2 — BND ≠ raw read either

> BND is a VCF-level breakend variant call.

A BND record is just two coordinates with a relationship — it's
"caller thinks something is going on between A and B". The actual
read evidence (split reads, supplementary alignments, discordant
mates) lives in BAM. The atlas displaying a BND should not pretend
it has shown the user the underlying reads.

### Rule 3 — INV does not consume BND

> Do not assume INV consumes BND.

Some callers emit ONE inversion as one INV record. Others emit
TWO BND records linked by EVENT/MATEID/PARID. Some emit BOTH.

When you see an INV at the boundary, you cannot assume there's no
BND nearby — there might be, with independent evidence. When you
see two BNDs nearby, you cannot assume they get merged into one
INV — they might already be linked, or they might be independent
events.

The atlas treats them as **independent records** that happen to be
clustered by zone (left_boundary / right_boundary / body / etc.)
and by EVENT/MATEID linkage. Clustering is for display, not for
deduplication of evidence.

### Rule 4 — A BND list is not a list of "raw pieces" of an INV

> Do not assume BND list necessarily contains raw pieces that made INV.

Manta and similar callers emit BND records with their own
provenance. They are NOT "the reads that supported the INV call".
They are separate variant claims that happen to be nearby.

If the user wants to know "what raw evidence supported this INV",
the answer is **read the BAM with samtools view at the INV
coordinates** — not "look at the BND records nearby".

### Rule 5 — Statistical transitions ≠ molecular breakpoints

> Do not claim exact molecular breakpoint from dosage / theta /
> PCA alone. Those define statistical transition zones. Exact
> breakpoint requires split reads, long reads, assembly junction,
> or strong BND/INV evidence.

This is the most important rule for the manuscript. The atlas
displays many "boundary zone" panels: PCA edge slope, dosage edge
slope, GHSL transitions, theta drops. All of these are
**statistical signals from the cohort**, not molecular precision.

The phrasing matters:

| ❌ Not allowed | ✅ Allowed |
|---|---|
| "Breakpoint at 14,823,100" | "Statistical boundary zone: 14,820k–14,830k" |
| "INV starts at position X" | "Cohort-level transition zone at ~14.82 Mb" |
| "We identified the inversion boundary" | "We bounded the inversion to a ~10 kb statistical boundary zone" |
| "Inversion 14.82M–17.71M" | "Inversion candidate, statistical boundary zones at ~14.82M and ~17.71M" |

The boundary widths in the manuscript should reflect the
statistical resolution (typically 1-10 kb at our window size), not
fake single-bp precision.

## The four-tier evidence framework

Tiered claims from strongest to weakest:

### Tier 1 — Breakpoint-resolving evidence

- A split-read alignment with both halves landing on opposite sides
  of the proposed breakpoint, with high MAPQ.
- An assembly junction (long-read or PacBio) crossing the breakpoint.
- A BND record from an independent caller, both ends with
  independent read support.

A Tier-1 event lets the atlas claim "breakpoint resolved at
position X" — but even then, only at split-read resolution
(typically ±10 bp).

### Tier 2 — Boundary-supporting evidence

- DEL/DUP/INS calls inside the boundary zone with H_class
  enrichment.
- Multiple SV calls clustering at the boundary zone from different
  callers.
- An INV record without independent read-support confirmation.

Claim: "boundary zone supported by SV cluster" — not "breakpoint
at X".

### Tier 3 — Linked internal markers

- SV records inside the inversion body that show H_class
  enrichment.
- These are **cargo / passenger** markers — they ride with the
  inversion but don't define its breakpoints.
- Useful for diagnostic PCR markers (page 12) but not for boundary
  resolution.

Claim: "linked marker within inversion body" — not "additional
breakpoint" or "internal breakpoint".

### Tier 4 — Family / repeat / noise artifacts

- SV calls correlating with hatchery family rather than H_class.
- SV calls in repeat-flagged regions or with high MAPQ0
  surroundings.
- Singletons with no group association.

These get filtered out of the catalogue — or shown only with a
"flagged: family artifact" / "flagged: repeat noise" badge so the
user can decide.

## Clustering vs deduplication

When the atlas sees:

```
INV   LG28:14820000-17710000     manta
BND   LG28:14820123 → LG28:17710003  manta  (EVENT=manta_inv_1)
BND   LG28:17710003 → LG28:14820123  manta  (EVENT=manta_inv_1, MATEID linked)
```

…this is **one cluster, one event**. The two BNDs are explicitly
linked to each other and to the INV. The atlas should:

- Display them as one "evidence cluster" in the SV table.
- Count it as ONE Tier-2 event (or Tier-1 if there's separate read
  support), not three.
- Show the cluster's strongest single record as the headline,
  with the others as expandable subrecords.

But this:

```
INV     LG28:14820000-17710000     manta     QUAL=247
INV     LG28:14821500-17708900     delly     QUAL=190
```

…is **two independent claims** (different callers, different
coordinates). Treat them as two records. Their agreement
strengthens the claim ("both manta and delly call an inversion in
this region") but they're independent evidence.

## "Evidence cluster" defined

A cluster is a set of SV records that satisfy ALL of:

- Same chromosome.
- Same zone (left_boundary / right_boundary / body / left_flank /
  right_flank).
- Within a tolerance Δ of each other (default Δ = 200 bp).
- Linked by EVENT/MATEID/PARID OR within Δ of each other from the
  same caller.

A cluster gets ONE entry in the table 5 enrichment results (S1).
The cluster's representative record is the one with the highest
QUAL (or most explicit SV type).

## What the atlas chrome must communicate

Every SV display element (chip, table row, hover, manuscript
exporter) should make these distinctions visible:

1. **Caller name** — never display an SV without naming who called it.
2. **Tier** — chip colour or label encodes Tier 1-4.
3. **Zone** — left_boundary / right_boundary / body / flank.
4. **Pattern label** — from S1's pattern vocabulary.
5. **Cluster size** — when the row represents a cluster, say so
   ("INV cluster of 3 records" rather than "INV at X").

These five fields combine to give the user (and the manuscript
reader) a defensible read of what the SV record means.

## Manuscript copy templates

Approved phrasings for SV-related claims in the manuscript draft:

✅ "Manta and DELLY independently call an INV spanning the
candidate region (Tier 2 evidence)."

✅ "Statistical boundary zones at ~14.82 Mb (left) and ~17.71 Mb
(right) are supported by clustered SV calls from two callers."

✅ "Within the inversion body, three DEL records show H_class
enrichment (linked internal markers; Tier 3)."

✅ "Breakpoint resolution requires long-read or assembly evidence;
short-read SV calls bound the breakpoints to a ~10 kb statistical
zone."

❌ "Breakpoint at 14,823,100." (without long-read or split-read evidence)

❌ "INV records confirm the inversion boundaries."
(INV records are CALLS, not confirmation.)

❌ "SVs identify the breakpoints precisely."
(They identify a zone, not a precise bp.)

❌ "Read evidence shows..."
(unless we actually showed reads via samtools view; SV records
≠ reads.)

## Where this discipline lives in code

This doc itself doesn't add code. It governs:

- Empty-state copy in P5.1 Section 6 (SV evidence) and Section 9
  (group association) — already partially baked in.
- Tooltip text on every SV chip / table column.
- The pattern-label vocabulary in S1's Table 5.
- The atlas exporter (TSV / Markdown) labelling.
- Future manuscript-figure exports.

When in doubt, the rule is: **understate**. The atlas should err
on the side of statistical-language framing. Going from "boundary
zone supported by SV cluster" to "breakpoint at X" requires a
specific Tier-1 event — not a vague aggregate of Tier-2 evidence.
