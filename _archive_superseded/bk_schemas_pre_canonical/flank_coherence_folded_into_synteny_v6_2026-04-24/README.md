# flank_coherence — superseded

**Date folded:** 2026-04-24 (chat B continuation)
**Folded into:** `synteny_v6.schema.json`

## Why

`flank_coherence` was declared in an earlier pass as a standalone block for the
gene-family order coherence data produced by `cross_species_bridge.py`. By the
time v6 of that script was written (`cross_species_bridge_v6.py`), the same
three fields were being emitted as part of the `synteny_v6` block payload —
but `synteny_v6`'s `keys_extracted` list never picked them up. Net effect: the
flank coherence keys were being persisted into the structured JSON but never
flowing to `keys.tsv`, while the `flank_coherence` schema sat unused.

## What changed

- The 5 flank-coherence properties (`q3_left_flank_coherence_score`,
  `q3_right_flank_coherence_score`, `q3_left_flank_n_families`,
  `q3_right_flank_n_families`, `q3_flank_coherence_class`) were already
  declared in `synteny_v6.schema.json#properties` — no change there.
- Added 5 corresponding entries to `synteny_v6.schema.json#keys_extracted`
  so those keys now land in `keys.tsv`.
- Writer side (`cross_species_bridge_v6.py`): **unchanged**. The data
  payload already carries these fields.
- This file (`flank_coherence.schema.json`) archived here as a historical
  record; no live code reads it.

## No action required from consumers

Downstream readers of `q3_*_flank_*` keys (e.g., `compute_candidate_status.R`)
see the same keys in `keys.tsv` as before the change; the difference is that
now they actually *appear* there, where previously they were silent.
