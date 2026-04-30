# Session audit — April 29-30, 2026

**Session length:** large multi-turn architectural session, continuing
from earlier ghsl/θπ atlas work via the prior compaction summary.

**Session theme:** four-stream architecture lock-in (was 3, briefly 5),
phase-4 resolution layer specification (CUSUM, dual-clustering, shelf QC),
and updated Q-regime detector spec with new metric families.

This audit is for the next session's Claude (and for future-Quentin).
It inventories what was built, what was decided, what wasn't built but
was planned, and which decisions to NOT re-litigate.

---

## 1. What was actually produced this session

### 1.1 Documents shipped to /mnt/user-data/outputs/

All written by this Claude during this session:

| File | Lines | Purpose |
|---|---|---|
| `SESSION_MANIFEST.md` | ~400 | Inventory of all artifacts (R scripts, atlas patches, SLURM launchers, prose docs) with done/pending status |
| `ARCHITECTURE_DECISIONS.md` | ~470 | 12 ADRs locking core architectural choices (3-stream → 4-stream framing, detector hierarchy, JSON shapes, etc.) |
| `JSON_CONTRACT.md` | ~570 | Per-layer schema spec for GHSL 8 layers + θπ 5 layers + cross-cutting fields |
| `HANDOFF.md` | ~180 | Short "if you only read one file" entry point with critical path |
| `ADR-13_q_and_ancestry_pages.md` | 87 | Locks the Q-page (ancestry-Q discovery) + Ancestry-page (phase-QC integration) split with asymmetric naming |
| `QR03_SPEC.md` (v1) | 301 | Initial Q-regime detector spec recovered from past chats |
| `specs/QR03_SPEC_v2.md` | 628 | Expanded Q-regime spec with 8 metric families, exp(H) priority, component-by-component matrices, hatchery angle |
| `specs/CUSUM_SPEC.md` | 464 | Per-sample CUSUM detector spec for GHSL (R02), θπ (T05), with shared kernel + DC06 concordance |

### 1.2 Atlas tree (in /home/claude/atlas_layout/Atlas/)

Current state of the canonical Atlas folder layout:

```
Atlas/
├── README.md
├── atlas.html  (pca_scrubber_v4_full_page3.html, four patches applied)
├── json/  (empty — populated by user from cluster outputs)
├── saved/  (empty — for user-saved candidate sets)
└── _scripts/
    ├── ARCHITECTURE_DECISIONS.md
    ├── ADR-13_q_and_ancestry_pages.md
    ├── HANDOFF.md
    ├── JSON_CONTRACT.md
    ├── QR03_SPEC.md  (v1)
    ├── README.md
    ├── SESSION_MANIFEST.md
    ├── sync_jsons.sh
    └── specs/
        ├── CUSUM_SPEC.md
        └── QR03_SPEC_v2.md
```

The full bundle as a tarball is `Atlas_complete_documented.tar.gz`
(648 KB, md5 `81685d8fdb6979132411b9dbff3a0d3e` for the previous bundle;
new bundle this session adds the v2 spec and CUSUM spec).

---

## 2. Architectural decisions locked this session

These are NOT to be re-litigated by the next session. If a future session
proposes reverting them, refer to the ADR or spec that locked them:

### 2.1 Four discovery streams, not five (ADR-1 update + QR03_SPEC_v2 §8)

Final: **dosage local PCA + GHSL phased haplotype divergence + θπ nucleotide
diversity + ancestry-Q regime**. Rare-MAF demoted to:

- `rare_sfs_pairwise.c` + `07_plot_rare_sfs_heatmap.R` → diversity atlas
  (separate planned product, NOT inversion atlas)
- `rare_inversion_score.R` → per-candidate annotation field on dosage
  candidates, NOT a 5th discovery stream

**Walk-back trigger:** rare-MAF has family-structure noise issues that
the four robust streams are designed to avoid. Don't add it back.

### 2.2 Q page + Ancestry page split (ADR-13)

Two new pages, deliberately asymmetric labels:

- **"Q" page** = per-stream discovery view (mirrors page 3 GHSL, page 12 θπ).
  Tracks-only until QR03 lands.
- **"Ancestry" page** = phase-QC integration with 11-track shelf, per-sample
  karyotype calls, cross-method concordance grid. Cross-stream view at the
  candidate region.

**NOT a tab toggle.** Two separate pages because they have different
scopes, different load dependencies, different audiences. The
"Ancestry" label reflects the page's PURPOSE (confronting ancestry
confounding), not its literal contents (which include Fst, θπ-by-genotype,
etc. — most tracks are non-ancestry).

### 2.3 CUSUM as per-sample resolution layer (CUSUM_SPEC §1, §11)

Three CUSUM applications, one shared kernel:

- `STEP_R02_ghsl_cusum.R` — GHSL rolling-divergence matrix
- `STEP_T05_theta_cusum.R` — θπ rolling matrix
- `STEP_SC01_emit_subcandidates.R` cell cliff-walker — per-cell mean tracks

Kernel lives in `phase_4_resolution/shared_lib/cusum_core.R`; scripts
differ in calibration only.

CUSUM produces per-carrier boundary distributions (5'/3' median + IQR +
spread_class) and per-candidate sub-system clusters (e.g. 45 carriers
clustered at 18.95 Mb + 9 carriers internal at 17.2 Mb suggesting nested
inversion). **These outputs feed the GHSL page and θπ page boundary panels
when zoomed to a candidate, plus the Ancestry page concordance grid.**

**Walk-back trigger:** "let's compute CUSUM in JS to skip cluster step" —
reject. Input matrices too big (~15 MB per chromosome per signal). JS
can re-CUSUM on user-defined sample subsets later as a v2 feature, but
default is cluster-side only.

**No 4th CUSUM for Q-regime per-sample.** Q-regime needs dispersion /
bimodality / per-sample deviation, not level shifts. Different problem.

### 2.4 Q-regime detector specifications (QR03_SPEC_v2)

**8 metric families** to compute in QR02:

1. Central tendency (dom_Q, Δ_ij)
2. Diversity (exp(H), ENA) ← **exp(H) is the cleanest single track**
3. Pairwise gaps (Δ12, Δ13, Δ14, Δ23, Δ24)
4. Cohort dispersion (cohort_var_Q[K])
5. Bimodality (dip_stat per K-component)
6. Component-by-component per-sample matrices ← **the LG28-killer**
7. Δ-track of Q vector from chromosomal background ← **GHSL-contrast analog**
8. Per-sample deviation from genome-wide-mean Q

**NGSadmix stays genome-wide K=8** on NAToRA-pruned 81 unrelated
individuals. Family-bias-resistance: inversions are tiny fraction of
genome, don't bias the genome-wide decomposition. **Don't re-run
per-chromosome to "fix" LG28.**

**Manuscript angle (QR03_SPEC_v2 §7):** first published inversion pipeline
to properly handle hatchery / breeding-program structure. Methodological
contribution for aquaculture genomics.

### 2.5 Folder reorganization to Atlas/ layout

Reorganize script (`reorganize_to_atlas_layout.sh`) tested on synthetic
messy folder. Idempotent, dry-run available. Renames atlas variants
(`pca_scrubber_v4_*.html`) to `atlas.html`. Renames per-chromosome JSONs
from `C_gar_LG##_phase*_*.json` to `LG##_<stream>.json` at the boundary.

---

## 3. What was specified but not built

All gated on real LG28 data from the cluster:

| Item | Status | Estimate when built |
|---|---|---|
| Replace 4 stub renderers (drawGhZ, drawGhLines, drawGhPCA, drawGhMDS) | spec done in JSON_CONTRACT, blocked on real data | ~2 days |
| STEP_TR_B v5 retrofit (sim_mat + MDS + sign-aligned loadings for θπ) | spec done, blocked on real data | ~2 days |
| Page 12 panel renderers (mirror page 3) | architecture clear, blocked on real data | ~2 days |
| Page-3-bis K-stripe heatmap view | sketched, blocked on real data | ~1 day |
| `export_q_to_json.R` skeleton | spec in QR03_SPEC_v2, blocked on QR02 actual output shape | ~1 day |
| `export_shelf_qc_to_json.R` | architecture in ADR-13, blocked on STEP_C01i decomposition output format | ~2 days |
| Atlas "Q" page scaffold | clones page 3/12, blocked on data | ~2 days |
| Atlas "Ancestry" page scaffold | per-candidate 11-track shelf, blocked on data | ~3-5 days |
| `STEP_R02_ghsl_cusum.R` + `STEP_T05_theta_cusum.R` + `cusum_core.R` | spec done, ~3 days cluster work | ~3 days |
| `STEP_DC06_cusum_concordance.R` | spec done, blocked on R02 + T05 | ~0.5 day |
| LG28 validation script for QR02 metric predictions | spec done, can run independent of full pipeline | ~1 afternoon |
| `STEP_QR02` + `STEP_QR03` + Q-regime atlas page | spec done, gated on validation script | ~7-8 days |

**Total speculative work specified:** ~25-30 days of cluster+atlas work
across 4 streams + phase-4 layer.

**Total work blocked on cluster output:** all of it. No code written this
session against real data.

---

## 4. Critical biological constraints (carry forward, never violate)

- **226 pure C. gariepinus hatchery cohort ONLY.** NOT F1 hybrid (assembly
  paper), NOT C. macrocephalus (future paper).
- **9× coverage**, hatchery F1 violates HWE → folded SFS, doSaf 5,
  doMajorMinor 1
- **K=8 NGSadmix** on NAToRA-pruned 81 unrelated individuals (genome-wide).
  Don't re-run per-chromosome.
- **K clusters represent hatchery broodline structure**, NOT geographic
  populations
- **GHSL grid 5kb fixed**; θπ at win10000.step2000 (~16,500 windows on
  LG28); dosage variable
- **Phase blocks from WhatsHap on short reads: 44–819 bp**
- **BAM manifest:** `pa_roary_results/00_manifests/sample_bam_minimap2_vs_P99TLENMAPQ30.tsv`
  (col3 = filtered BAM)
- **Rscript binary:** `/lustrefs/disk/project/lt200308-agbsci/13-programs/mambaforge/envs/assembly/bin/Rscript`
- **Conda env:** `assembly`
- **Working dir:** `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/`

---

## 5. Failure-mode patterns observed this session

These are patterns this Claude fell into. Future Claude should watch
for them:

**Pattern 1: Reflexive symmetry-imposition.**
I tried to make ancestry-Q a per-stream discovery page mirroring GHSL
and θπ before checking whether the discovery detector existed. Quentin
pushed back: ancestry IS phase QC right now; QR03 doesn't exist yet.
**Lesson:** asymmetry is OK when the upstream pipeline is asymmetric.

**Pattern 2: Inventing work to do.**
When Quentin said "continue" without specific direction, I drifted toward
producing more documents. Quentin called it "Claude inventing work to keep
busy." **Lesson:** when nothing has positive expected value pending real
data, stop and let the cluster run. Ask Quentin what he actually wants
rather than producing speculative work.

**Pattern 3: Over-correction toward "remove the optional thing for safety."**
Repeated through the session. Quentin's correct response: "include cheap
things by default." This became ADR-12. **Lesson:** if it costs almost
nothing to include and might be useful, include it.

**Pattern 4: Conflating different versions / scope drift on CUSUM count.**
The other-Claude conversation went 1→2→3→4→3 on how many CUSUMs there are.
This Claude initially missed that θπ and GHSL also have CUSUM (only saw
it for the Q stream initially). **Lesson:** when the user says "you forgot
X for Y," check past chats / images carefully. The initial CUSUM_SPEC
draft was incomplete because of this.

**Pattern 5: Spec verbosity drift.**
QR03_SPEC v1 was 301 lines, v2 is 628. Some growth was warranted (new
metric families) but some was decoration. **Lesson:** specs should be
read by future-Claude; if a section won't change behavior, cut it.

---

## 6. What the next session should do FIRST

In strict priority order:

### 6.1 Cluster work (Quentin's, not Claude's)
1. Run LG28 dry-run from `RUNBOOK_produce_phase2_jsons.md` Section 1
   (~1 hour wall clock)
2. Browser smoke test: drag-drop `LG28_ghsl.json` and `LG28_theta.json`
   into `Atlas/atlas.html`, verify no console errors, verify panels render
   (~5 min)
3. Send back evidence: screenshot, structural dump, or console error log

### 6.2 If cluster work succeeded
Build phase-2 atlas in priority order:
1. Replace 4 stub renderers on page 3 with real renderers (highest ROI)
2. Add page 12 panel renderers (mirror page 3)
3. Page-3-bis K-stripe heatmap view (lowest ROI of the page 3/12 work)

### 6.3 If page 3 / 12 stable, then phase-4 begins
1. Build `cusum_core.R` shared kernel (~1 day)
2. Build `STEP_R02_ghsl_cusum.R` + `STEP_T05_theta_cusum.R` (~2 days
   together; they share the kernel)
3. Build `STEP_DC06_cusum_concordance.R` (~0.5 day)
4. Build CUSUM rendering on GHSL + θπ pages (per-carrier boundary
   panels)
5. Build Ancestry page scaffold (concordance grid + 11-track stack)

### 6.4 In parallel (gated on validation, not on phase-2/4)
1. Run LG28 Q-regime validation script (one afternoon, see QR03_SPEC_v2 §5)
2. If predictions hold → build STEP_QR02 (~3 days) → STEP_QR03 (~1 day)
   → export_q_to_json.R → Q page

### 6.5 What NOT to start until later
- Manuscript methods sections (need real results)
- Synthetic test JSON (only if used before LG28 — Quentin was not
  interested)
- Diversity atlas (separate product, not blocking inversion work)

---

## 7. Tools / files / paths (quick reference)

### Local working dirs
- `/home/claude/atlas_layout/Atlas/` — current Atlas tree
- `/home/claude/atlas_session/atlas_v4_session_apr29/` — session work files
- `/mnt/user-data/outputs/` — files shipped to user this session

### Cluster paths (for reference, no access from Claude)
- Codebase: `/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/inversion_codebase_v8.5/`
- LANTA cluster, SLURM, account `lt200308`
- Conda env `assembly`
- Rscript binary path above

### Atlas state
- Base atlas: `pca_scrubber_v4.html` (md5 `37eb28da49c79bf8b67322fa6ea4617c`)
- Final 4-patch: `pca_scrubber_v4_full_page3.html` (md5 `c6c3ba48bb797884f57b944d2a31b01e`)
- Bundle (previous): md5 `81685d8fdb6979132411b9dbff3a0d3e`

### Past-chat references for archaeology
If next session needs to recover more context:
- "Nested subgroup extraction from PCA and theta profiles"
  (`92fef806-5a81-45f7-bbdc-eaf0e1df6791`, April 25 2026) — canonical
  Q-regime architecture
- "Live PCA visualization with contingency tables" (`819d8454-...`,
  April 28) — phase 4 dual-clustering
- "Insulation scoring and TAD boundary detection" (`ac219483-...`,
  April 10) — engine B / Q-vector idea origin
- "MODULE_2B audit" (`5b793a68-...`, April 10) — full pop-genomics
  pipeline S00-S18 breakdown

---

## 8. Open questions / deferred decisions

These are explicitly not resolved this session and are flagged in the
relevant specs for later resolution:

- **θπ sim_mat sizing decision** (option A coarse / B banded / C int8) —
  deferred in ADR-7. Decide on first real run.
- **Multi-scale Δ12 in QR02** — defer to QR02 implementation choice.
- **Whether to ship per-component matrices in the JSON** — currently NO
  (too big, ~120 MB per chromosome at K=8). Could change if user wants
  interactive component browsing in atlas.
- **JS-side re-CUSUM for user-defined sample subsets** — currently NO.
  Defer until cluster-on-demand pattern proves insufficient.
- **`enriched_components` field translation to breeding-line labels** —
  needs MODULE_2B's Hungarian-matched K↔line mapping. Translation at
  atlas render time, not at JSON emit time.
- **Whether the diversity atlas is a separate HTML file or a mode-switch
  in the same atlas** — both architectures documented (line 4004 of the
  scrubber HTML mentions mode-switch). User preference TBD.

---

## 9. Final state assessment

**What's solid:**
- Architecture is locked, documented, walk-back-protected
- All specs are internally consistent
- 13 ADRs cover the core decisions
- Two detailed specs (CUSUM, QR03_v2) cover the unbuilt phase-4 layer

**What's NOT solid:**
- **Zero code has been tested against real data.** ~3500 lines R + 1500
  lines Python written across the session, all blocked on cluster output.
- Several "fragile assumptions" flagged in SESSION_MANIFEST that need
  first-run verification:
  - JSON layer field names match between R exporter and atlas renderer
  - sim_mat upper-triangle packing math is correct
  - sign-alignment of MDS eigenvectors works as specified
  - K=3 GHSL palette displays correctly across browsers
  - reorganize script handles edge cases (symlinks, partial Atlas trees)

**Honest scope statement:** Quentin has a fully-documented architecture
for a 4-stream + phase-4-resolution inversion pipeline. The architecture
is sound and reviewer-defensible. Implementation is approximately 1
week of cluster work + 1 week of atlas work + ~1 week of phase-4
extension, gated on each other. Total real-time: ~3-4 weeks of focused
work, contingent on no biological surprises in LG28 validation.

**Critical path:** LG28 dry-run → page 3/12 atlas validation → phase-4
build → manuscript figures.

Sleep / cluster work / Clair3 finishing — whichever comes next.
The architecture will hold.
