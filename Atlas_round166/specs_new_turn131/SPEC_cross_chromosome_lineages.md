# SPEC — Cross-chromosome lineage aggregator

**Status**: drafted turn 130 final session. Not yet implemented.
Estimated ~1 turn after fish-trajectory Slice 1 ships (it has).

**Trigger** (recovered from chat `f74cf5d4`):
> *"That's exactly approach 1 or 2 — the chromosomal background is
> preserved across the whole chromosome regardless of where on the
> chromosome you look. The signal is genome-wide consistency of band
> membership."*
>
> Implied across multiple chats: a fish that's lineage A on LG12 is
> probably lineage A on LG28 too, because their underlying chromosome
> ancestry is the same individual.

The fish-trajectory lineage compute (turn 130 Slice 1) runs **per chromosome**.
A fish gets lineage_id 0 on LG12 and lineage_id 0 on LG28 — but those
"0"s are independent labels. The chromosome backgrounds may or may not
match. This spec builds the cross-chromosome aggregator that says:

> *"Fish 12, 47, 89, 134 share lineage A on LG12. Of those, 12, 47, 89
> share lineage A on LG28. Cohort-wide, these 3 fish form a coherent
> chromosome-wide lineage group."*

---

## 1. Why this matters

Hatchery cohorts have ~8 founder broodlines (Quentin's K=8 NGSadmix
context). Each broodline gives rise to a coherent chromosome ancestry
across the genome. If a fish carries broodline-3 chromosomes, its
fish-trajectory lineage call should be consistent at every chromosome
**modulo recombination + inversion-system polymorphism**.

So:
- Fish with chromosome-wide consistent lineage calls → strong founder
  inheritance signal.
- Fish whose lineage calls differ between chromosomes → either
  inter-broodline recombination (cross between two founders' children)
  or chromosome-specific inversion polymorphism.

Both are biologically informative.

## 2. Algorithm

### 2.1 Input

`state.lineageResultsByChrom: Map<chrom, lineageResult>` — populated as
the user opens chromosomes. Each result has `lineage_id_per_sample` of
length n_samples.

### 2.2 Pairwise chromosome alignment

For each pair (chrom_a, chrom_b):

1. Build a contingency table `T[lineage_a][lineage_b]` = count of fish
   with lineage_a on chrom_a AND lineage_b on chrom_b.
2. Run Hungarian assignment to find the best lineage_a → lineage_b
   permutation (maximizing diagonal).
3. Compute alignment quality = sum(diagonal) / n_samples.
4. Apply the permutation to chrom_b's labels so the "0" lineage on
   chrom_a is the same as the "0" lineage on chrom_b after alignment.

### 2.3 Genome-wide lineage call

After all pairs aligned to a reference chromosome (default: longest
chromosome, or LG28 because it's the validation chrom):

```
fish_genome_lineage[fish_i] = mode of lineage_id across all chromosomes
                              where the fish has a valid call
```

Plus per-fish stability score:

```
stability[fish_i] = (count of chromosomes where fish_i is in mode lineage)
                  / (count of chromosomes with valid call)
```

Stable fish (stability ≥ 0.8) → strong genome-wide lineage signal.
Unstable fish (stability < 0.5) → recombinant or inversion-driven
chrom-specific lineage shifts.

## 3. Surface

### 3.1 New atlas page or G-panel tab

Tab `genome lineages` in the G-panel (sibling to per-chrom `lineages` tab).
Shows:

- 226 × n_chrom heatmap, rows = fish, columns = chromosomes, color =
  aligned lineage_id
- Stability score per fish (right-side bar)
- Click a fish → highlight in the per-sample-lines panel of every
  chromosome
- Filter by stability threshold

### 3.2 Export

TSV: `fish_id | per_chrom_lineage_array | mode_lineage | stability`.

## 4. Implementation slices

### Slice 1 — pairwise alignment (~0.5 turn)
- `_alignLineagesAcrossChroms(resultA, resultB)` Hungarian-aligns labels
- `_alignmentQuality(table)` returns diagonal sum / n
- Cache result keyed by chrom pair

### Slice 2 — genome-wide lineage call (~0.5 turn)
- Reference-chrom selection (longest by default, or user-set)
- Apply alignment cascade
- Compute mode + stability

### Slice 3 — heatmap UI (~0.5 turn)
- Canvas-based 226 × n_chrom heatmap
- Stability strip
- Hover: fish_id + per-chrom calls

## 5. Open questions

1. **Reference chromosome**: longest, validation (LG28), or user pick?
   Default: longest, override available.
2. **Alignment chain order**: pairwise to reference vs sequential
   (LG01 → LG02 → ... → LG28)? Pairwise is simpler and parallel.
3. **What about chrom-specific real biology?** A fish that's HET at an
   LG28 inversion could legitimately have a different lineage at LG28
   than elsewhere. Stability < 1.0 isn't necessarily noise; it's a
   signal of inversion polymorphism. Surface it as such.

## 6. Tests

- Synthetic case: 226 fish across 3 chroms, all share same lineage
  pattern → stability = 1.0 for every fish.
- Synthetic case: half the fish flip lineage on chrom 3 → stability
  drops to 2/3 for those fish; remains 1.0 for the others.
- Hungarian alignment with K=4 lineages → recovers the right permutation
  on a planted contingency table.

## 7. Dependencies

- Fish-trajectory lineage compute Slice 1 (shipped turn 130).
- Multi-chrom JSON load orchestrator (`SPEC_multichrom_load_orchestrator.md`).
- G-panel scaffold (`SPEC_g_panel_unified_groups.md` Slice 1).
