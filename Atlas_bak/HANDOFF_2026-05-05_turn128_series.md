# HANDOFF — turns 128 / 128c / 128d (continuous session)

**Date**: 2026-05-04 / 2026-05-05
**Session bundle**: `sv_evidence_bundle_FULL_v3.zip` (this archive)
**Atlas line count**: 63,075 lines
**Test status**: 1653 PASS / 0 FAIL across all suites

---

## 1. Reload procedure for the next chat

```bash
cd /home/claude && mkdir -p Atlas && cd Atlas
unzip -oq /mnt/user-data/uploads/sv_evidence_bundle_FULL_v3.zip
pip install --break-system-packages --quiet fastapi pydantic pyarrow httpx
```

Test runner pattern (use this exact form — output format `PASS: N FAIL: N`
is what the runner greps):

```bash
total_p=0; total_f=0
for f in tests/sv_evidence/test_step*.js tests/dosage_bridge/*.js tests/dosage_viewer/*.js \
         tests/test_turn120_scree_inset.js tests/test_turn128_*.js tests/test_turn128c_*.js tests/test_turn128d_*.js; do
    out=$(node "$f" 2>&1 | tail -3)
    p=$(echo "$out" | grep -oP "PASS: \K\d+")
    fails=$(echo "$out" | grep -oP "FAIL: \K\d+")
    total_p=$((total_p + ${p:-0})); total_f=$((total_f + ${fails:-0}))
done
for f in tests/producers/*.py tests/dosage_bridge/*.py tests/dosage_viewer/*.py; do
    out=$(python3 "$f" 2>&1 | tail -3)
    p=$(echo "$out" | grep -oP "PASS: \K\d+")
    fails=$(echo "$out" | grep -oP "FAIL: \K\d+")
    total_p=$((total_p + ${p:-0})); total_f=$((total_f + ${fails:-0}))
done
echo "Combined total:  PASS $total_p / FAIL $total_f"
```

Expected: **1653 PASS / 0 FAIL**.

---

## 2. What shipped in this session

### Turn 128 (initial bugs from Quentin)

| # | Fix | Tests |
|---|---|---|
| Pre-existing | scree spec hygiene + rescued orphan test | 33 |
| 1 | JS scripts badge in header (sibling of #schemaBadge) | 55 |
| 2 | Symmetric remove-candidate (page-2 button & page-4 ✕) | 29 |
| 3 | Window-mode cycler "lit" state + sidebar sync | 26 |
| 4 | PCA tracked-sample trail alpha attenuation on non-PC1 axes | 16 |
| 5 | Spacebar contextual cut shortcut (candidate mode) | 16 |

### Turn 128c

| # | Fix | Tests |
|---|---|---|
| 6 | Cycler → L3 compareUnit sync (was calling never-defined helper) | 28 |
| 7 | CTRL+- DPR refit for L3 mini-PCAs (rAF + matchMedia) | 21 |
| 8 | Spacebar play/pause unbind (4 added assertions to test #5) | (in #5) |

### Turn 128d

| # | Fix | Tests |
|---|---|---|
| 9 | "Help > Help" single-child stage collapse | 24 |
| 10 | Compact tracked-samples settings parity with fixed | 32 |
| 11 | L3 contingency default layout = leftright (was focal) | 10 |

**Combined: 290 new tests added across the session, all green.**

---

## 3. Specs delivered this session (NOT YET IMPLEMENTED in code)

### `specs_todo/SPEC_l2_envelope_live_split.md` (turn 128)

The deferred half of Quentin's spacebar request. The spacebar shortcut for
L3 cuts shipped (turn 128 #5), but the bigger feature — actually splitting
an L2 envelope into two new L2s on `state.data.l2_envelopes` — is deferred
because it breaks three architectural invariants:

1. Saved-candidate `l2_indices` index stability
2. K-cluster cache invalidation
3. Save-Session round-trip

Spec details design choices (synthetic L2 ids, replacement markers, sidecar
catalogue export, undo). **Estimated: ~3 turns** if Quentin signs off.

### `specs_todo/SPEC_l3_het_dosage_coloring.md` (turn 128d)

Het-rate coloring on L3 contingency mini-PCAs. **Slice 1 partially started**:
`state.l3HetColoring` slot added, L3 toolbar `[ ] het` checkbox markup added,
but the renderer body + cache + drawMiniPCA hook + color-bar legend were
NOT plumbed (turn ran long).

**To finish Slice 1 (~0.5 turn)**:
- Wire the checkbox change handler to mutate `state.l3HetColoring` and
  call `renderL3Panel()`.
- Implement `_computeHetRateForL2(l2idx)` returning `Float32Array[n_samples]`
  of het rates from `dosage_chunks`. Cache on `state.__hetRateCache`.
- Implement `_hetRateColor(rate)` returning RdBu CSS color (FIG_C08 ramp).
- In `drawMiniPCA`: when `state.l3HetColoring`, override dot fill with
  `_hetRateColor(hetRate[si])`. Halo stays K-cluster colored.
- Color-bar legend under each mini-PCA when toggle is on.
- Cache invalidation hook on `dosage_chunks` load.
- Tests against a synthetic chunk.

The spec also defines Slice 2 (tracked-samples always-on het halo, deferred)
and lists 4 open design questions.

### `specs_todo/SPEC_g_panel_unified_groups.md` (turn 128d, Quentin-flagged)

The spec Quentin asked to add to the next chat: a unified G-panel popup
with three tabs (karyotype / inheritance / manual). Answers his question
about whether karyotype + inheritance + manual groups belong together.

**Yes** for those three (they're all "phenotype groupings of the cohort").
**No** for family groups (kinship is a separate ontology — different
color mode, different palette, different semantics).

The reason karyotype groups didn't load into the existing manual-groups
system: manual groups are static snapshots, karyotype groups are live
derivations from the focal candidate. The fix is a new popup with a
live karyotype tab, not auto-populating the manual list.

**Slice 1 (~1 turn)**: popup scaffold + manual tab re-host (existing
manualGroups UI in roomier popup). Slices 2+3 add karyotype + inheritance
tabs. **Estimated total: ~3 turns across 4 slices.**

The spec explicitly lists 4 open design questions for Quentin to resolve
before Slice 1:
1. Sidebar button location (above Atlas tools / in lines header / both?)
2. Hotkey: `g` or `G`?
3. K=6 karyotype labels — show 6 unnamed bands or refuse?
4. Family groups — fourth tab as read-only directory, or keep separate?

#### Late-turn diagnosis: why the grouping panel is bugged today

Quentin flagged after the spec was written:
> "The grouping panel is bugged — it doesn't auto select or compute when
> we click it for the candidate. It doesn't find the 'groups' also
> shouldn't it be related or linked to the groups in candidate page as
> well? Or its a different grouping?"

**Diagnosis**: there's a real bug in the data pipeline that the new
G-panel will need to fix. The inheritance compute pipeline reads from a
DIFFERENT registry than the one the user actually populates:

| Registry | Type | Populated by | Read by |
|---|---|---|---|
| `state.candidateList` | Array | user actions: promote, lock-and-promote, draft confirm | catalogue, page 2, candidate strip, exports |
| `state.candidates` | Dict (id→cand) | Save-Session import ONLY | `_gatherActiveCandidatesForInheritance` (line 34362) → `runInheritanceCompute` → I·g pills |
| `state.candidates_detailed` | Dict | mode-switch detailed pipeline only | same gatherer, detailed branch |

So when the user saves 5 candidates interactively, `state.candidateList`
has 5 entries but `state.candidates` is still empty. The I·g compute
silently finds < 2 items and returns null. **The grouping panel is
reading from the wrong slot.**

The candidate-page bands and inheritance groups ARE the same biological
concept — both derive from `candidate.locked_labels` (the K-means
cluster assignment frozen at promote time). They're not different
groupings; they're the same groupings rendered through divergent reader
code paths. The page-2 band-composition cards work fine because they
read `state.candidate.locked_labels` directly from a single candidate;
the inheritance pipeline tries to traverse `state.candidates` and finds
it empty.

**Fix options for the G-panel spec to absorb**:

a. **Bridge approach** (cheap): make `_gatherActiveCandidatesForInheritance`
   read from `state.candidateList` instead of `state.candidates`. One
   function change. Risks: anywhere downstream that mutates
   `state.candidates` (Save-Session import) loses its effect. May need
   to merge into candidateList on import.

b. **Single source of truth refactor** (clean): converge on
   `state.candidateList` everywhere; deprecate `state.candidates` /
   `state.candidates_detailed` dicts; rebuild the `state.candidate /
   .candidate_detailed / .activeMode` dispatch on top of array filters.
   2–3 turns of careful work.

c. **Auto-sync** (pragmatic middle ground): on every mutation of
   `state.candidateList` (add / remove / restore from localStorage),
   rebuild `state.candidates` as a derived index. Slim helper
   `_rebuildCandidateRegistries()` called from the existing 4–5
   mutation entry points. ~1 turn.

**Recommended**: option (c) — small, low-risk, unblocks the I·g compute
+ the future G-panel tabs. The G-panel Slice 1 should NOT ship the
popup until this bridge is in place; otherwise the karyotype +
inheritance tabs will be empty-state on first load even with saved
candidates.

**Updated G-panel ordering for the next chat**:

- **Pre-Slice (~0.3 turn)**: ship `_rebuildCandidateRegistries()` bridge.
  Hook into the candidateList mutation points. Tests verify that adding
  a candidate to candidateList causes `runInheritanceCompute()` to find
  it. Single-candidate test ⇒ inheritance returns null (need 2+);
  two-candidate test ⇒ result has > 0 groups.
- **Slice 1**: G-panel scaffold + manual tab.
- **Slice 2**: karyotype tab (now will work because the bridge exists).
- **Slice 3**: inheritance tab (now will compute correctly).

Quentin's intuition that the panels are linked is **correct** — they
are the same groupings. The codebase has lost the bridge between the
user's saved-list and the computation registry along the way; the
G-panel is the natural surface that surfaces both views, but the
plumbing fix has to land first.

---

## 4. Open observations (in `OBSERVATIONS_TO_FIX.txt`)

Still unfixed and waiting:

- **Bug 1**: L3 arrow up/down promote/demote sub-panel highlighting +
  merge broken (highest-risk).
- **Bug 3**: page-1 PCA legend overlap with leftmost stripe — blocked
  on screenshot from Quentin (couldn't locate the rendering site).
- Boundaries page Mb mismatch (~14 vs 15.75)
- Top-bar hover gap (session/mode/data folders refold)
- L3 detailed click no-op
- Diversity atlas plot placement (page 1 samples)
- ST2 plot bottom margin
- ROH heatmap S8b last-row collapse
- Loess fit at zero coverage

---

## 5. Project conventions (refresher)

### Cohort discipline (CRITICAL)

Three cohorts must NEVER be conflated:
1. **F₁ hybrid** (*C. gariepinus* × *C. macrocephalus*) — genome assembly
   paper only.
2. **226-sample pure *C. gariepinus*** hatchery — current inversion work
   (THIS bundle). K clusters = hatchery broodline structure, NOT species
   admixture.
3. **C. macrocephalus** wild — future paper.

Quentin's full name: **Quentin Andres**. Never invent surname.

### File / edit conventions

- **JS-string footgun**: literal `</script>` inside JS string terminates
  inline `<script>` block. Project idiom: `'<' + '/script>'`.
- **Test output format**: must be `PASS: N   FAIL: N` for runner grep.
- **Test path resolve**: `path.resolve(__dirname, '..', 'Inversion_atlas.html')`
  not hardcoded `/home/claude/work/build/`.
- **Spec discipline**: when a spec lands as code, move from `specs_todo/`
  to `_archive/specs_done/` (NOT `specs_done/`) with a one-line
  shipped-in-turn-N header.

### File creation strategy

- For surgical edits to existing code, use `str_replace`.
- For new tests, use `create_file` with full paths in `tests/`.
- For new specs, write to `specs_todo/SPEC_<topic>.md`.

---

## 6. Key infrastructure constants (unchanged)

- Reference: `fClaHyb_Gar_LG.fa` (28 pseudochromosomes)
- KING / Manichaikul thresholds: 1st≥0.177, 2nd≥0.0884
- NAToRA-pruned 81 unrelated samples (NGSadmix K=8)
- `ngsRelate` theta in column 18
- `state.k` = 3, `purity_threshold` = 0.80, `mergeThr` = 0.85
- `minNGroup` = 5, `minNWin` = 5, `alpha` = 0.001
- Scrubber JSON layers at `scrubber/data/<chrom>/`

### State slots added/changed in this session

- `state.l3Layout`: default changed `'focal'` → `'leftright'` (turn 128d)
- `state.l3HetColoring`: NEW boolean, default `false` (turn 128d, Slice 1
  partial)
- localStorage keys: `pca_scrubber_v3.stepmode` (now persisted from cycler
  too, turn 128c); `pca_scrubber_v3.l3HetColoring` (planned, turn 128d
  Slice 1)

---

## 7. Suggested next-turn priorities

In rough order of value-per-turn:

1. **Candidate registry bridge** (`_rebuildCandidateRegistries`, ~0.3 turn).
   Unblocks the I·g pills, fixes the "grouping panel doesn't compute"
   bug Quentin flagged at turn end, and is a hard prerequisite for the
   G-panel karyotype + inheritance tabs.

2. **Finish Slice 1 of het-coloring** — checkbox handler + `_computeHetRateForL2`
   + `_hetRateColor` + `drawMiniPCA` hook + color-bar legend + tests.
   ~0.5 turn. Quentin will see immediate visual payoff.

3. **G-panel Slice 1** — popup scaffold + manual tab re-host. ~1 turn.
   Lays groundwork for Slices 2+3 once the registry bridge is in place.

4. **Bug 1** (L3 arrow promote/demote broken) — high-risk, needs
   careful state-machine reading. ~1 turn.

5. **Other observations** — boundaries Mb mismatch, top-bar hover gap,
   etc. — short fixes, mix in opportunistically.

---

## 8. File deltas (this session)

```
Inversion_atlas.html              62338 → 63075  (+737 lines, 11 fixes)
OBSERVATIONS_TO_FIX.txt           updated with 6 [FIXED] markers + 2 spec
                                  pointers

NEW tests:
  tests/test_turn128_js_scripts_badge.js              55
  tests/test_turn128_remove_candidate_symmetric.js    29
  tests/test_turn128_window_cycler_lit.js             26
  tests/test_turn128_pca_trail_alpha.js               16
  tests/test_turn128_spacebar_cut_in_candidate_mode.js 20
  tests/test_turn128c_cycler_l3_sync.js               28
  tests/test_turn128c_resize_dpr.js                   21
  tests/test_turn128d_help_single_child.js            24
  tests/test_turn128d_compact_tracked_settings_parity.js 32
  tests/test_turn128d_l3_default_leftright.js         10

RESCUED test (previously orphaned):
  tests/test_turn120_scree_inset.js                   33

NEW specs:
  specs_todo/SPEC_l2_envelope_live_split.md
  specs_todo/SPEC_l3_het_dosage_coloring.md
  specs_todo/SPEC_g_panel_unified_groups.md      ← Quentin-flagged for next chat

MOVED:
  specs_todo/SPEC_pca_scree_inset.md → _archive/specs_done/  (was already shipped turn 120)
```

---

## 9. Quentin's running ask list (from his messages)

Quentin's explicit feature requests during this session, in order received:

1. ✅ JS scripts badge in top bar (turn 128)
2. ✅ Make page-2 candidate-clear and page-4 ✕ behave the same (turn 128)
3. ✅ Spacebar to cut L2 in candidate mode (turn 128, partial — see L2-split spec)
4. ✅ Window cycler buttons should "light up" when active (turn 128)
5. ✅ PCA trail alpha when both axes are non-PC1 (turn 128)
6. ✅ Window 5/10 cycler should sync to L3 compareUnit (turn 128c)
7. ✅ Unbind spacebar from play/pause (turn 128c)
8. ✅ CTRL+- distortion of L3 mini-PCAs (turn 128c)
9. 📋 Het (dosage) coloring on L3 panels — SPEC + state slot + UI checkbox shipped, renderer next turn (turn 128d)
10. ✅ "Help > Help" redundancy — collapse to just "Help" (turn 128d)
11. ✅ Compact tracked-samples settings parity with fixed (turn 128d)
12. ✅ L3 default layout = +1L+1R, not Focal (turn 128d)
13. 📋 G-panel with karyotype + inheritance + manual tabs — spec written for next chat (turn 128d)
14. 📋 Grouping panel doesn't auto-compute / find groups for the candidate — diagnosed turn 128d as `state.candidateList` (array, user-populated) vs `state.candidates` (dict, only loaded via Save-Session) divergence; the inheritance compute reads the wrong slot. Bridge fix queued as Pre-Slice in the G-panel work. Confirmed: it IS the same grouping as the candidate page (both derive from `locked_labels`), the codebase just lost the bridge.

---

## 10. End-of-session note

Atlas shipped, tests green, two specs queued for the next chat per
Quentin's explicit choice (het-coloring Slice 1 finish + G-panel). The
backups in `Inversion_atlas.html.bak_*` were dropped before bundling
to keep the archive size sane (38MB saved). Re-establish in-session
backups in the next chat per the existing pattern (`bak_pre_<feature>`,
`bak_post_<feature>`).
