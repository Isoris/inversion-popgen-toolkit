# SPEC OVERVIEW — Multi-species architecture

**Status**: forward-looking integration map. Read this before the individual specs to understand how the pieces compose.

**Sister specs** (one per topic, intentionally separable):
- `SPEC_hypothesis_registry_and_multispecies.md` — the hypothesis-and-evidence object model that consumes everything below
- `SPEC_comparative_te_breakpoint_fragility.md` — TE density across species as fragility proxy
- `SPEC_phylogenetic_tree_integration.md` — using the catfish tree to place events on branches
- `SPEC_busco_anchors.md` — single-copy BUSCO proteins as deep-divergence anchors

---

## What this overview is for

The atlas is growing from a single-species inversion finder into a **multi-species evolutionary workbench**. Five evidence layers are involved at different divergence depths and resolutions. None of them is sufficient alone; together they triangulate.

The naive design would put all five into one giant spec. That'd be 3000+ lines of unreadable mush. The honest design is one spec per layer + this overview that shows how they compose.

## The five evidence layers, stacked by divergence depth

```
DEPTH        LAYER                        RESOLUTION    SPARSITY    USE FOR
─────────────────────────────────────────────────────────────────────────────────
Within-Cgar  226-sample resequencing      bp-level       dense       polymorphism (the only layer that proves it)
Cgar↔Cmac    wfmash 1-to-1 PAF             bp-level       dense       2-species breakpoint catalogue (current page 16)
~10-30 Mya   miniprot whole-proteome       protein-level  medium      multi-species synteny anchors
~30-150 Mya  BUSCO single-copy             protein-level  sparse      chromosome-homology backbone, deep tree
all depths   Comparative TE-density        100kb-windows  dense       fragility proxy (proxies polymorphism in non-resequenced species)
all depths   Phylogenetic tree             topology       —           branch placement of events; ancestral state reconstruction
```

Each layer answers a different question. None of them by themselves answer the manuscript's main biological question, which is:

> **Are some Gar polymorphic inversion breakpoints sitting on ancient, repeat-rich chromosome architecture that's been reused across catfish evolution?**

This question requires combining: per-Gar polymorphism (resequencing) + per-species TE density (assemblies) + chromosome homology (anchors) + tree topology (branch placement). Hence the multi-layer architecture.

## Resolution + sparsity tradeoff

The further from Cgar in evolutionary distance, the sparser the signal:

| Comparison | DNA alignment (wfmash) | Miniprot | BUSCO | TE annotation |
|---|---|---|---|---|
| Cgar vs Cmac (~25 Mya) | dense, breakpoint-resolved | dense gene anchors | sparse but present | per-assembly TE |
| Cgar vs Pangasius (~80 Mya) | gappy, low signal | medium | reliable | per-assembly TE |
| Cgar vs Trichomycterus (~150 Mya) | mostly fails | sparse | reliable | per-assembly TE |

This is why we need all four. Wfmash for the close pair, miniprot for the medium-divergence species, BUSCO for the deep outgroups, TE annotations everywhere because they're available everywhere.

The tree ties it together by placing each evidence observation on a branch.

---

## What the atlas needs from each layer

Each layer ships as a standardized JSON that the atlas reads, with format `<layer_name>_v1.json`. All schemas defined in their respective specs.

| Layer | JSON key | Atlas page that primarily consumes | Status |
|---|---|---|---|
| 226-sample data | (already in precomp) | pages 1-12 | shipped |
| Cgar↔Cmac wfmash | `cs_breakpoints_v1` | page 16 | shipped |
| Multi-species miniprot | `miniprot_anchors_v1` | page 16, page 14 | spec'd, not built |
| BUSCO anchors | `busco_anchors_v1` | page 16, page 14 | spec'd, not built |
| Comparative TE-density | `comparative_te_breakpoint_fragility_v1` | boundaries page (page 4), page 16, page 14 | spec'd, not built |
| Phylogenetic tree | `phylo_tree_v1` | page 14, page 16 | spec'd, not built |
| (Future) cross-species synteny graph | `synteny_multispecies_v1` | page 14, page 16 | spec'd in hypothesis-registry spec, not built |

The hypothesis registry on page 14 is the **integrator** — it consumes all of the above and scores hypotheses against them.

---

## Build order

If everything were built today, the rational order would be:

1. **Hypothesis registry data model** (chat 115a) — needed before anything else, because it's the consumer
2. **Comparative TE-density** (multi-chat) — most useful evidence we can produce from existing data, no resequencing needed
3. **Phylogenetic tree integration** (chat 115h) — Quentin already has the tree; small wiring job
4. **BUSCO anchors** (chats 115k-m) — Quentin has the faa, ~3 chat days
5. **Miniprot anchors** (chats 115c-e) — needs miniprot pipeline output first
6. **Cross-species synteny graph** (catfish-synteny-toolkit) — needs the toolkit run on real data first
7. **Auto-suggest classification** (chat 115f) — depends on layers 1-6
8. **Hypothesis-scoring engine** (chat 115g) — depends on 7

In practice, work proceeds opportunistically — whatever Quentin has data for next, plus whatever atlas turns are small enough to fit in one chat.

---

## What this overview does NOT do

- Detailed schemas — those are in the per-topic specs
- Implementation plans — those are in the per-topic specs
- Manuscript wording — captured in each spec separately + the help page on atlas

This overview is just the **map**. The detail lives in the topical specs.

---

## Anti-overclaiming discipline (cross-cutting)

This applies to ALL evidence layers above and is a critical thread through every spec:

| Confirmed | Allowed comparative inference | Never |
|---|---|---|
| polymorphic inversion in Gar broodstock | homologous breakpoint is TE-rich | species X is polymorphic |
| 226-sample carrier counts | breakpoint boundary conserved across catfish assemblies | this inversion exists in species X |
| local-PCA bands + dXY | candidate recurrent rearrangement hotspot | structural polymorphism in non-resequenced species |
| | possible ancient breakpoint substrate | |
| | predicted polymorphic hotspot | |

Every JSON output should carry explicit `polymorphism_confirmed_in_species: true|false` and `polymorphism_unknown_in_species: true|false` fields per species. The atlas's render code should refuse to print "polymorphic" claims about species where the former is false.

---

## Practical motivation (the part that matters for manuscript reviewers)

We will not be resequencing populations of Pangasius from Vietnam, Silurus from China, Trichomycterus from Peru, Heterobranchus from Africa anytime soon. The money + time + permits are not coming. But assemblies + TE annotations DO exist for all of these species already.

So instead of saying "we don't know," we say:
- "We confirm polymorphism in Cgar."
- "We predict which homologous regions in other catfish are likely structurally fragile, based on TE-rich, conserved synteny boundaries."
- "Predictions are testable when population data becomes available for those species."

This is honest, useful, and publishable. It's also a shape that other inversion-genomics groups can adopt for their own non-model species — turning the catfish manuscript into a small methods contribution beyond just the catfish biology.

---

## Reading order for an implementing chat

A fresh chat coming in to work on any of this should read:

1. `HANDOFF_2026-05-01_session_summary.md` — what's been done
2. `ATLAS_REFERENCE_for_phase4b_synthesis.md` — atlas internals
3. **This overview** — to understand which layer they're working on
4. The relevant single-topic spec for that layer
5. The atlas's help-page section "Breakpoint classification framework" (open `Inversion_atlas.html`, click tab 16) — for the schema vocabulary

Then pick one chat-sized scope and ship it. Don't try to build multiple layers at once.
