# STANDARDIZE_FOLDER.md
# -----------------------------------------------------------------------------
# Super-prompt for standardizing one folder of the inversion-popgen-toolkit.
# Paste everything below the `PROMPT BEGINS` line into a new chat, then attach
# or paste the scripts from one folder. The LLM will return standardized
# versions of each.
#
# Usage pattern:
#   1. Pick one folder (e.g. phase_qc_shelf/, or Modules/MODULE_4B_delly_del/)
#   2. New chat. Paste this prompt.
#   3. Attach every script file in the folder (or paste them one by one).
#   4. Review the output per-file. If it looks right, save over the originals.
#      Commit to git so you can see the diff and revert individual files.
#   5. Move to the next folder.
#
# Do NOT send more than ~20 files at once. Accuracy drops on long batches.
# Prefer 5-15 files per chat session.
# -----------------------------------------------------------------------------

===================  PROMPT BEGINS  ===================

You are standardizing scripts in Quentin's inversion-popgen-toolkit for a PhD
thesis + manuscript. The toolkit has ~640 scripts across many folders and
languages (bash, R, Python, C). Your job is to apply a consistent header
format + runtime-log format to one folder at a time.

## Hard rules (read before each file)

1. **Do NOT change executable logic.** No rewriting, no refactoring, no
   "improvements", no style normalization of code, no variable renaming.
2. Allowed changes:
   - Header comment block (new or replaced)
   - Inline comments (added, or existing ones kept intact)
   - Logging calls (add new `qc_log`/`message`/`print` lines; normalize
     existing ones to the conventions below)
   - Whitespace between sections (add blank lines for readability only)
3. If uncertain whether a line is logic or decoration: do NOT touch it.
4. If a file is compiled source (.c, .cpp, .h): only update the header.
   Do NOT add runtime logging to compiled binaries unless they already
   have printf/fprintf statements you're just normalizing.
5. If a file is under 20 lines and mostly logic: leave it unchanged and
   note "too small for headerable surface area."

## Output format (use this structure PER FILE)

For each file I send you, respond with exactly this structure:

```
================  FILE 1/N: <filename>  ================

**DIFF SUMMARY**
- line <N>: <what changed>    (one line per change)
- no logic lines modified.

**RISK FLAGS**
- <anything the user should look at>  (or "none")

**UPDATED FILE**
```<lang>
<complete updated file, every line>
```

================  END FILE 1/N  ================
```

If I send 5 files, give 5 of these blocks, one after another. Between
blocks, include a brief "moving to file N+1" line. This helps me review.

## Header standard (required, identical in every file)

Every script gets a two-layer header: a terse banner of fields, then a
narrative explanation. Language-appropriate comment syntax (`#` for bash/R/
Python/awk; `//` for C/C++).

### Layer 1: terse banner

```
# =============================================================================
#  <FILENAME>
#  -------------------------------------------------------------------
#  <ONE-LINE PURPOSE, MAX 100 CHARS>
# =============================================================================
#
#  PURPOSE   <1-line purpose, can repeat the one above>
#  INPUTS    <file pattern/path>   <short description>
#            <second input if any>  <...>
#            or "none (top-level driver)"
#  OUTPUTS   <file pattern/path>   <short description>
#            or "none (emits to stdout only)"
#  CONFIG    <env var or config knob>  <description>
#            or "none beyond 00_config.sh"
#  USAGE     <one or two example invocations>
#  STATUS    STABLE | WIP | ORPHAN              LAST UPDATED   <ISO date>
#  CALLED BY <parent script(s)>  or  "??? — grep to confirm"
#
```

### Layer 2: narrative explanation (REQUIRED for scripts over 30 lines)

```
# -------------------------------------------------------------------
#  What this step is doing
# -------------------------------------------------------------------
#  <2-4 sentences, plain English, explains why this exists and what
#   biological/technical signal it produces or consumes. If you don't
#   have enough context from the code alone, leave a placeholder
#   comment like "??? — Q, please fill in: what is the biological
#   reading of this output" so the user can revisit.>
# =============================================================================
```

For files under 30 lines, Layer 2 is optional — skip it and go straight
to code after Layer 1.

## Runtime-log standard (bash only)

Every bash STEP_Q* script uses helpers defined in `00_config.sh`:
- `qc_log "msg"` (timestamped info)
- `qc_warn "msg"` (⚠️ prefix, greppable)
- `qc_err "msg"` (❌ prefix, greppable)
- `qc_die "msg"` (error + exit 1)
- `qc_banner_open <STEP> <UNIT>` (open separator + banner)
- `qc_banner_close <STEP> <UNIT> <seconds>` (close banner with timing)
- `qc_config_snapshot <STEP> VAR1 VAR2 ...` (print config at job start)
- `qc_preview_file <label> <path> [n_lines]` (emit head-N with ASCII rails)
- `qc_log_grep_tips` (print grep recipes at end of run)

Required placements in a bash step script:
1. After argument parsing, call `qc_config_snapshot "Q0X" <relevant vars>`
2. Inside each per-unit function (e.g. `run_one`):
   - Start with `t0=$(date +%s)` and `qc_banner_open "Q0X" "${unit}"`
   - Preview every input file with `qc_preview_file "input" "${path}" 3`
   - Preview every output file after it's written: `qc_preview_file "output" "${out}" 3`
   - End with `elapsed=$(( $(date +%s) - t0 ))` and `qc_banner_close ...`
3. At script end: `qc_log "🏁 Q0X DONE. ..."` then `qc_log_grep_tips`

Any pre-existing `qc_log` / `echo` / `printf` calls should be kept if they
convey useful information, but normalized to use the helpers above.
Remove bare `echo` in favor of `qc_log` for consistent timestamping.

## Runtime-log standard (R workers)

R scripts use `message()` (stderr) for logs. No shared helpers library,
but mimic the bash pattern:

- First substantive line after library loading:
  `message("[q0X] starting: ", Sys.time())`
- Print config vars the R worker uses:
  `message("[q0X] args: precomp=", PRECOMP, "  chrom=", CHROM, "  ...")`
- At each major stage:
  `message("[q0X] reading ", FILE_PATH)`
  `message("[q0X] writing ", OUT_PATH)`
- If writing a TSV/CSV/RDS, emit a 3-line preview of the output:
  ```r
  preview <- utils::head(written_dt, 3)
  message("[q0X] output preview:")
  for (ln in capture.output(print(preview))) message("    ", ln)
  ```
- At end: `message("[q0X] done")`

## Runtime-log standard (Python utilities)

Python scripts use `print(..., file=sys.stderr)`:

```python
import sys
def log(msg):
    print(f"[{sys.argv[0].rsplit('/', 1)[-1]}] {msg}", file=sys.stderr)
```

Same cadence as R: starting, args, reading, writing, output preview, done.

## What to infer vs mark ???

Infer freely from the file content:
- PURPOSE / INPUTS / OUTPUTS from what the code does
- CONFIG from env vars referenced (`${BIN_MB}`, `Sys.getenv("X")`, etc.)
- STATUS: default to STABLE for files that look complete and WIP for
  those with TODO comments, commented-out sections, or placeholder logic

Leave as ??? (and flag in RISK):
- CALLED BY — you can't know this from a single file. Write "??? — grep
  for <basename> in phase_*/*sh to confirm"
- LAST UPDATED — use today's date
- Anything that requires reading another file to understand

## Per-folder adjustments

Before starting a folder, tell me:
- Folder name and its role (e.g. "phase_qc_shelf is the shelf-diagnosis
  QC module")
- Any folder-specific conventions I should apply (e.g. "Q-prefix for
  steps" vs "STEP_<N>_<name>" vs no prefix)

If the folder has no step prefix convention (e.g. the launchers in
Modules/MODULE_4*/launchers/), use the file's basename as the STEP label:
`qc_config_snapshot "LAUNCH_batches_4ch" ...` etc.

## Before you start: sanity check

Count the files I sent you. State the count back to me in your first line:
"Received N files. Processing in order: <list>. Starting with file 1."

If any file is larger than 20 kb, warn me and ask whether to truncate or
process in parts.

If any file looks like it might be generated output (.tsv, .csv, .rds
with .R extension, .pdf) rather than source code, ask before processing.

## Example

If I send `STEP_Q02_beagle_uncertainty.sh`, your response should look like:

```
Received 1 file. Processing in order: STEP_Q02_beagle_uncertainty.sh.
Starting with file 1.

================  FILE 1/1: STEP_Q02_beagle_uncertainty.sh  ================

**DIFF SUMMARY**
- lines 1-35: replaced 3-line header with two-layer standardized header
- line 42: added `qc_config_snapshot "Q02" BEAGLE_DIR QC_TRACKS`
- line 58: added `qc_banner_open "Q02" "${chr}"` at start of run_one
- line 62: added `qc_preview_file "input" "${beagle}" 3` before read loop
- line 130: added `qc_preview_file "output" "${out}" 3` after write
- line 135: added timing wrapper (t0/elapsed) and qc_banner_close
- line 142: added `qc_log_grep_tips` at end
- no logic lines modified.

**RISK FLAGS**
- CALLED BY left as ??? — grep for STEP_Q02_beagle_uncertainty in
  phase_qc_shelf/*.sh to confirm (looks like run_all.sh and run_chrom.sh).
- Line 89 has `|| true` swallowing an error from the R worker — left as-is
  but flagging so you can decide if that's correct.

**UPDATED FILE**
```bash
...complete file here...
```

================  END FILE 1/1  ================
```

## Ready

Send me:
1. The folder name and one-line role.
2. The list of files you're sending (count).
3. The files themselves (paste or attach).

I'll process them one at a time in the order received and emit the
standardized versions.

===================  PROMPT ENDS  ===================
