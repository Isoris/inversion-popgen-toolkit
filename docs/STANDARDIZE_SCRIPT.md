# STANDARDIZE_SCRIPT.md — super prompt for header/echo standardization

## Purpose

A single prompt that can be given to an LLM (with one script file attached)
to produce a standardized version: upgraded header, clear stage-transition
echoes, no logic changes. Designed to be run one file at a time, reviewed,
committed — iterated until the output quality is reliable, then batched.

## How to use

1. Copy the prompt body below into a chat with an LLM.
2. Attach ONE script file. The LLM must respond with the full updated file
   and a clearly-marked diff summary at the top.
3. Commit the original before applying the update so `git diff` is clean.
4. Read the diff summary. If any non-comment/non-echo lines changed, reject
   the update.

---

## THE PROMPT (copy everything below the `---` line)

---

You are standardizing a single script file for Quentin's pipeline
(`inversion-popgen-toolkit`). Your ONLY task is to upgrade the header
comment block and add/normalize progress echoes. You must NOT change
program logic.

### Hard constraints

1. **NEVER change executable logic.** Do not rewrite code, rename variables,
   restructure control flow, modernize idioms, "improve" style, reorder
   arguments, or anything else that could change behavior.
2. Changes allowed: header comments, inline comments, stage echoes, blank
   lines for readability.
3. Changes forbidden: anything that alters what the program does.
4. If you are uncertain whether a line is logic or comment, DO NOT TOUCH IT.
5. If the script is very short (<20 lines) and logic-heavy, return it
   unchanged with a note saying "no headerable surface area."

### Required output format

Respond with THREE sections, in order:

**1. DIFF SUMMARY** (markdown bullet list)

For each change, one line: `line N: <what changed>` where "what changed"
is one of: `added header`, `added echo`, `added inline comment`,
`normalized whitespace`. If you touched anything else, say so.

**2. RISK FLAGS** (markdown bullet list, or "none")

Any lines you were unsure about. Anything that looked like a bug but
you left alone. Anything the user should look at.

**3. UPDATED FILE** (the complete file inside a fenced code block with
language tag)

No truncation, no "... unchanged ..." ellipses, no summaries. Emit every
line of the full file exactly as it should be saved.

### Required header format

Every script gets a header block in its language's comment syntax
(`#` for bash/R/python/awk, `//` or `/* */` for C/C++). The block must
contain EXACTLY these fields, in this order:

```
<divider line of = characters>
<filename> — <one-line purpose, max 100 chars>
<divider line of = characters>

PURPOSE
  <1-3 sentences of what this file does and why it exists. Written in
   plain English. No jargon the user won't remember in 6 months.>

INPUTS
  - <name/path>: <short description>
  - <name/path>: <short description>
  [or: "none (top-level driver)"]

OUTPUTS
  - <path/pattern>: <short description>
  [or: "none (emits to stdout only)"]

DEPENDENCIES
  - <binary/library/tool>: <version if known>
  [or: "none beyond base language runtime"]

USED BY
  - <script that calls this one, if any>
  [or: "top-level — invoked by user"]

STATUS: STABLE | WIP | ORPHAN
LAST UPDATED: <date from git log if available, otherwise leave as
              "YYYY-MM-DD" placeholder>

USAGE
  <example command line, typically one or two lines>

NOTES
  - <anything tricky, any gotcha, any TODO>
  [or omit this block entirely if nothing to say]

<divider line of = characters>
```

### Required echo conventions

Every meaningful stage transition gets a stderr log line. Conventions:

- **Bash**: `qc_log "STEP_NAME: <what I'm about to do>"` if `qc_log` is
  defined in the sourced config, otherwise `printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "<msg>" >&2`.
- **R**: `message("[<script_name>] <what I'm about to do>")`.
- **Python**: `print(f"[<script_name>] <msg>", file=sys.stderr)`.
- **C**: `fprintf(stderr, "[<binary_name>] <msg>\n");`

Rules:
- Emit at file start ("starting"), at each major stage ("reading X",
  "processing Y", "wrote Z"), and at exit ("done in T seconds" if you
  have timing, otherwise just "done").
- Do NOT add echoes inside tight loops (> 1000 iterations). Only at
  chunk/batch boundaries.
- If the script already has echoes, normalize them to the convention
  above but do not remove information. If it prints Indonesian, French,
  or other language mix, keep the language but make the format consistent.

### What to infer vs what to mark ???

You may infer:
- PURPOSE from the existing code behavior (what files it reads, what it
  writes, how it's structured).
- INPUTS / OUTPUTS from `fread` / `read_csv` / `open()` / `cat` / `>` calls.
- DEPENDENCIES from `library()` / `import` / `source` / `#include`.

You must mark ??? and flag in RISK section:
- USED BY — you cannot know this from a single file; leave as
  "??? — search with grep after"
- LAST UPDATED — leave as today's date in ISO format
- STATUS — default to `WIP` unless obviously stable; user will correct.
- Anything about intent that requires knowledge outside this file.

### Example diff summary

```
**DIFF SUMMARY**
- lines 1-28: replaced old 3-line header with standardized header block
- line 45: added `qc_log "STEP_Q07: reading BEAGLE"` before the read loop
- line 89: normalized existing echo "Done writing" to `qc_log "STEP_Q07: wrote $out"`
- no logic lines modified

**RISK FLAGS**
- line 137: saw `rm -f ${tmpfile}` inside a loop — left as-is, but user should
  verify `tmpfile` is set correctly (no quotes around it in a set -u context).
- USED BY field left as ??? — grep for "STEP_Q07_popstats" in inversion_modules/
  to find callers.

**UPDATED FILE**
...
```

---

## End of prompt

After you iterate with the LLM on 5 sample files and the header standard
feels right, freeze this prompt. Then batch through the rest with a
simple wrapper script that:

1. For each script in a file list:
2. Call the LLM with this prompt + the script contents
3. Save the updated file to a staging directory
4. `git diff` the staging version against original
5. Reject if any non-comment lines changed (safety check)
6. If clean, write to the real location and `git add`

The batch wrapper is itself a one-time script we write AFTER the prompt
is proven on examples. Don't automate before you've verified the output
quality on a handful.
