# Archived: ancestry_bridge.R (pre-unified)

Moved here 2026-04-17 (chat 15).

## What this was

`inversion_modules/utils/ancestry_bridge.R` — a legacy bridge script that
exposed `merge_local_Q_into_invlikeness()` for the inversion precomp
step to merge per-window Q summaries into `<chr>.precomp.rds`.

## Why it's archived

The same API is now provided by `unified_ancestry/wrappers/instant_q.R`
with two advantages:

1. **K-parameterised**: `merge_local_Q_into_invlikeness(inv_like_dt, K = NULL)`
   lets the caller pick any cached K level, not just canonical.
2. **Sample-set-tagged cache paths**: files under
   `<LOCAL_Q_DIR>/K<NN>/<chr>.<sample_set>.local_Q_*.tsv.gz` can't collide
   between subset analyses.

## Where callers should look now

Anywhere `source("inversion_modules/utils/ancestry_bridge.R")` appeared,
use `source("utils/load_bridge.R")` instead. The bridge sources
`unified_ancestry/wrappers/instant_q.R` and calls `configure_instant_q()`
automatically, so `merge_local_Q_into_invlikeness()` lands in the global
namespace with the same signature.

The updated `STEP_C01a_precompute.R` (phase_2_discovery/2c_precomp) does
this. No other script in the flattened layout still sources the old bridge.

## Retention

Kept for historical reference and to enable sanity-comparison of output
between the old and new pipelines if any merge-column value ever disagrees.
Safe to delete after ≥1 full pipeline rerun confirms identical results.
