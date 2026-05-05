# SPEC — Live L2 envelope split (cut L2 → make new L2)

**Status**: drafted turn 128, not implemented.
**Trigger**: Quentin's request (turn 128 user message): *"I would like to
have spacebar to cut the interval when in candidate mode we can cut L2
intervals to make more L2. Its easier to merge."*

**What shipped in turn 128 (minimum-viable interpretation):**
spacebar in candidate mode now acts as a shortcut for the existing `C`
key — toggle an L3 cut at cursor inside the current draft. The L2 registry
is **not** mutated. This is keyboard ergonomics only.

**What this spec covers (deferred):** literally subdividing an L2
envelope at the cursor window into two new L2 envelopes, mutating
`state.data.l2_envelopes` live in-session. After the split, the user can
extend a draft into one half independently of the other half (the merge
workflow Quentin described).

---

## Why this is non-trivial

The L2 envelope registry is currently treated as immutable across a
session — it ships in the precomp JSON, gets indexed once at load time
(`state.windowToL2`, `state.l2NeighborsInL1`), and downstream code
assumes the indexes stay stable. Splitting an L2 mid-session breaks
several invariants:

### 1. Index stability across `state.candidateList`

Existing saved candidates carry `l2_indices: [int]`. If we splice a new
L2 into the middle of the array, every saved candidate's indices >= the
splice point shift by +1. Either:
- (a) Use insertion-at-end (new L2 gets index `N`) and store an explicit
  `parent_l1_id` + `start_w` so neighbor lookup re-derives correctly. This
  preserves saved-candidate indices but breaks any UI that assumes L2
  arrays are sorted by `start_w`.
- (b) Insert in-place (sorted by `start_w`) and update every existing
  candidate's `l2_indices`, `ref_l2`, etc. Many touch points.
- (c) Use a stable string-key registry instead of integer indices. Major
  refactor; rejected for this feature.

Option (a) is the cleanest if downstream code can be made order-agnostic.

### 2. Cached K-cluster results per L2

Each L2 has cached K-means cluster assignments (`getL2Cluster(l2idx)`)
that depend on the L2's window range. When we split, both halves need
fresh K-means runs. The cache invalidation is straightforward but the
re-cluster cost is non-trivial for narrow halves (n_windows < 5 may not
have enough signal to cluster meaningfully).

Open question: should narrow halves be marked unclusterable and rendered
greyed out, or should the split be rejected when either half would have
< some threshold of windows?

### 3. Sim-mat envelope overlay drawing

`drawSim()` paints the L2 envelope rectangles on the sim_mat heatmap
with greens. After split, the new envelopes need to render correctly.
This is a redraw call, not a structural change.

### 4. Catalogue / page-4 / page-8 / Save Session round-trip

- Catalogue (page 4): rows reference `l2_indices` — see point 1.
- Page 4 karyotype: per-sample regime breakdown is computed from
  `state.candidate.l2_indices` → propagates from point 1.
- Page 8 windows table: filters by L2 membership using
  `state.windowToL2` — needs rebuild after split.
- Save Session: `state.data.l2_envelopes` is part of the snapshot; if
  the in-session split adds new L2s, Save must include them and Load
  must restore them. The schema currently has no `_synthetic_l2: true`
  flag to distinguish precomp L2s from session-edits.

### 5. Persistence to the catalogue export TSV

The catalogue TSV columns include `l2_indices` (joined). Exported
catalogues from a session that split L2s would have indices that
don't correspond to any single precomp's L2 numbering. Either:
- Tag exported indices with a session UUID so loaders know to re-index.
- Re-emit the L2 envelopes alongside the catalogue export.
- Document that splits are session-local and don't persist.

---

## Suggested minimum-viable design (for a future spec turn)

1. **Trigger**: spacebar in candidate mode, cursor strictly inside a
   single L2 (and not at its edges) AND no draft yet promoted from
   that L2. Add a clear hotkey hint.

2. **Mutation**:
   - Compute the focal L2 from `state.windowToL2[cur]`.
   - Split into L2_left (start_w → cur) and L2_right (cur+1 → end_w).
   - Append both to `state.data.l2_envelopes` at the end (option (a)
     above), each with a synthetic `id` like `<orig_id>_L` and `_R`,
     and `_session_split: true`.
   - Mark the original L2 as `_session_split_replaced_by: [Lidx, Ridx]`
     and skip rendering it. (Alternative: actually delete + reindex. More
     work but cleaner.)
   - Recompute `state.windowToL2` for windows in the affected range.
   - Recompute `state.l2NeighborsInL1` (only the chain the original L2
     was in).
   - Invalidate the K-cluster cache for the original L2 and prime new
     entries for the two halves.

3. **UI affordance**:
   - On split, briefly flash a gold outline on the two new L2 panes in
     the L3 contingency carousel so the user sees what just happened.
   - Header info chip: "L2 #N split → #N_L (n_w=12) + #N_R (n_w=8)".

4. **Undo**:
   - At minimum, an `Undo split` button in the L3 toolbar that reverts
     the last split. Stack-based undo is overkill for v1.
   - Saved-candidate `l2_indices` referencing a session-split L2 stay
     valid via the replacement marker.

5. **Persistence**:
   - Save Session round-trips the `_session_split: true` envelopes.
   - Catalogue export includes a sidecar `l2_session_splits.json` so the
     catalogue's `l2_indices` resolve correctly on a fresh load.
   - The cluster-side R pipeline should NEVER re-emit a precomp that
     contains session-split envelopes — they're atlas-side only.

6. **Tests**:
   - Pure-function: a `_splitL2EnvelopeAtWindow(state, l2idx, cut_w)` that
     returns the new state without mutating. Heavy unit tests on this.
   - Integration: split + extend-draft into half + commit candidate, then
     verify the catalogue TSV round-trips.
   - Regression: a session with N splits, save → reload → N splits still
     present.

---

## Estimated effort

- Pure mutation function + invariants: ~1 turn.
- Wire to keyboard handler + UI affordances: ~0.5 turn.
- Round-trip persistence + catalogue export: ~1 turn.
- Tests: ~0.5 turn.
- **Total: ~3 turns** if Quentin signs off on the design choices above.

The current turn-128 spacebar binding is the cheapest interpretation that
gets Quentin a faster cut shortcut. If the existing L3 cut model (cuts
inside a draft, applied at commit time to split the candidate into
sub-intervals) already covers the "make merging easier" workflow, the
literal L2-split is unnecessary. Discuss with Quentin before
implementing.
