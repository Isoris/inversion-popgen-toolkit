# DROP_README — chat A deliverable (2026-04-24)

Three independent items in one drop. Drag the whole tree into GitHub Desktop;
no files need to be removed first, no folder moves, nothing renamed at the
path level.

---

## 1. Data flow audit — `docs/DATA_FLOW_AUDIT.md` (new file, 25 KB)

Chat A's primary deliverable. Static-analysis audit of the schema →
`build_key_spec()` → `characterize_candidate.R` wiring. Identifies:

- **7 schemas missing the `block_type` field** — P1 silent-failure, fix first.
- **`frequency.schema.json` + `frequency.v2.schema.json` share `block_type="frequency"`** — ambiguous `load_schema()` lookup.
- **208 real phantom keys** (spec expects, nothing writes) across Q1–Q_QC_SHELF.
- **149 orphan keys** (schemas produce, spec doesn't register).
- **42 of 90 consumer reads in `characterize_candidate.R` (46%) target phantoms** — actively reading keys nothing writes.

The audit is the input to chat B. No code fixes in this drop — the audit
is diagnosis, chat B is treatment. Read the §7 "what this audit does NOT
cover" before trusting any number blindly — static analysis has real
limits and I documented them.

---

## 2. Pass 15b — `v5_blocks` naming rename (2 files, 16 touchpoints)

Pure rename as logged in the pass-15 handoff. No logic change.

| Old | New |
|---|---|
| `V5_BLOCKS` (bash var) | `EVIDENCE_BLOCKS` |
| `v5_blocks` (subdir path) | `evidence_blocks` |
| `--v5_blocks_dir` (CLI flag) | `--evidence_blocks_dir` |
| `v5_dir` (python local) | `evidence_dir` |
| `cand_v5` (python local) | `cand_evidence` |

Files touched:

- `inversion_modules/phase_8_evidence_biology/run_evidence_biology.sh` (7 hits)
- `inversion_modules/phase_9_classification/assign_structural_class_v7.py` (9 hits: 1 argparse + 1 var + 7 cand_v5 call sites)

**Intentionally NOT touched:** `docs/HANDOFF_2026-04-23_v7_wiring_plan.md`
still contains `--v5_blocks_dir` at line 361. That's a dated historical
handoff and per the pass-15 session-audit lesson, historical docs are left
as-is. If you want it updated for consistency, it's a 1-line sed.

No external callers — grepped the whole tree for `v5_blocks_dir` and the
only call site is `run_evidence_biology.sh` itself calling
`assign_structural_class_v7.py`, both of which are in this drop.

---

## 3. Three parse bug fixes (3 files, ~3 lines each)

Each was flagged in pass 15 verification. All pre-existing, unrelated to
the audit or the rename. Each fix confirmed against the original to show
the bug was real and the fix resolves it.

### 3.1 `Modules/MODULE_4A_clair3_snp_indel/steps/postprocess/STEP_P10_publication_figure.R:282`

R uses `%in%` for set membership, not SQL-style `in`.

```r
# before
if (nrow(weak) > 0 && "N_READS_SUPPORT_INDEL" in names(weak)) {

# after
if (nrow(weak) > 0 && "N_READS_SUPPORT_INDEL" %in% names(weak)) {
```

### 3.2 `inversion_modules/phase_5_qc_triage/R/q04_compose_plot.R:394`

R parser requires `if { ... } else { ... }` to be one expression.
The original split `if ... "SNPs/kb"` and `else sprintf(...)` into
two statements across a newline, which R reads as assignment + orphan
`else`. Wrapped in braces.

```r
# before
sd_unit <- if (sd_scale_kb == 1) "SNPs/kb"
           else sprintf("SNPs/%gkb", sd_scale_kb)

# after
sd_unit <- if (sd_scale_kb == 1) {
  "SNPs/kb"
} else {
  sprintf("SNPs/%gkb", sd_scale_kb)
}
```

### 3.3 `unified_ancestry/engines/hobs_hwe/scripts/01_build_subset_bamlists.sh:245`

The `%+%` operator was being used inside an R heredoc but defined
*outside* the closing quote — bash saw `'%+%' <- function(a, b) paste0(a, b)`
as stray tokens after the `Rscript -e "..."` call, which the `2>/dev/null || true`
silently swallowed. Moved the `%+%` definition inside the heredoc where
it actually belongs, removed the stray trailing tokens.

---

## Verification

Every file in this drop passed parse checks in the sandbox:

- `bash -n` on both `.sh` files: OK
- `python3 ast.parse` on the `.py`: OK
- `Rscript -e "parse(...)"` on both `.R` files: OK
- MD fence balance on the audit: OK (0 fences — prose-only doc)

The 3 pre-existing parse bugs in the original source tree were confirmed
to fail the same checks; the drop versions all pass.

---

## What's not in this drop

- **Chat B fixes.** The audit is the deliverable; chat B does the actual
  schema/spec repairs. Wait for you to review the audit before picking
  what to fix first — several findings might be static-analysis
  false-positives (see §7 of the audit).
- **`phase_10_followup/` anything.** Off my plate per your instruction.
- **Folder renames / phase_7 wrapper README / sub-folder letter
  schemes.** Left alone per your decision in chat A.

---

## Workflow

1. Drag the four files + one new doc into GitHub Desktop at the
   corresponding paths (the tree structure here mirrors the repo).
2. Commit. No `rm` needed anywhere.
3. Push. No HPC pull required yet — all 3 items are diff-only at the
   laptop level; LANTA pulls after you commit.
