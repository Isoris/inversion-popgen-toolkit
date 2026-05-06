# P1.4 — `LDSplitReq` vs `FastLDReq` naming

## Risk: low
## Lines changed: ~10
## Depends on: nothing
## Verification: server starts cleanly without ImportError; LD endpoint round-trips

---

## What's wrong

Two folders coexist:

- `Atlas/server_turn1/` (the canonical server)
- `Atlas/server_turn11c_ld_fast/` (the fast LD spike that got merged in)

`popstats_server.py` imports:
```python
from ld_endpoint import LDSplitReq, handle_split_heatmap
```

But the file actually deployed at `server_turn1/ld_endpoint.py` was
copied from `server_turn11c_ld_fast/fast_ld_endpoint.py` which
defines `class FastLDReq(BaseModel)`, not `LDSplitReq`.

Temporary patch added:
```python
LDSplitReq = FastLDReq  # compat alias
```

This is a hidden mismatch. Future readers will see `LDSplitReq`
imported and won't realize it's actually `FastLDReq` underneath.

## Two clean options

### Option A — rename `FastLDReq` → `LDSplitReq` (recommended)

This is the option that hides the spike origin and makes the
production code self-consistent.

In `Atlas/server_turn1/ld_endpoint.py`:

```python
# before
class FastLDReq(BaseModel):
    chrom: str
    start_bp: int
    end_bp: int
    # ...

LDSplitReq = FastLDReq   # <-- delete this line

# after
class LDSplitReq(BaseModel):
    """Request shape for /api/ld/split-heatmap.
    Was previously named FastLDReq during the server_turn11c_ld_fast
    spike; renamed in P1.4 to match the existing import in
    popstats_server.py.
    """
    chrom: str
    start_bp: int
    end_bp: int
    # ...
```

Search the file for any internal reference to `FastLDReq` and update.

### Option B — keep alias but add comment

If something downstream (a script, a notebook, a test) imports
`FastLDReq` directly, do NOT rename. Instead:

```python
class LDSplitReq(BaseModel):
    """Canonical name for the request shape. Used by
    popstats_server.py via `from ld_endpoint import LDSplitReq`."""
    chrom: str
    start_bp: int
    end_bp: int
    # ...

# Backward-compat alias for the server_turn11c_ld_fast era. Some
# tests or scripts may still reference the old name. Don't drop
# this without grepping the workspace.
FastLDReq = LDSplitReq
```

Then `LDSplitReq` is the canonical name and `FastLDReq` is the
alias (the OPPOSITE of what's there now). Imports of either name
work. The hidden mismatch is gone — both are documented to be the
same class.

## Decide between A and B

Run this in the Atlas root:
```bash
grep -r "FastLDReq" --include="*.py" --include="*.js" .
```

- Zero hits outside `ld_endpoint.py` → Option A.
- Hits in tests, notebooks, scripts → Option B.

## Cleanup of `server_turn11c_ld_fast/`

Once `server_turn1/` is canonical, the parallel folder shouldn't
exist. Either:

- Delete `server_turn11c_ld_fast/` if it's already fully merged.
- Or move it to `Atlas/_audits/server_turn11c_ld_fast_archive/`
  with a README explaining what was merged where.

Don't leave two parallel servers — exactly the source of the
import confusion.

Suggested README for the archive:
```md
# server_turn11c_ld_fast — archived (P1.4)

This folder contained the fast LD heatmap spike. Its contents were
merged into server_turn1/ on <DATE>:

- fast_ld_endpoint.py → server_turn1/ld_endpoint.py
- lazy_windows_json.py → server_turn1/lazy_windows_json.py

The class FastLDReq was renamed to LDSplitReq during the merge.
See patches/P1_4_ld_naming.md.
```

## Verification

```bash
cd Atlas/server_turn1
mamba activate atlas-popstats
python -c "from popstats_server import app; print('ok')"
```

Should print `ok` with no ImportError.

```bash
python popstats_server.py --config popstats_server.local.yaml
# in another shell:
curl -s -X POST http://127.0.0.1:8766/api/ld/split-heatmap \
  -H 'content-type: application/json' \
  -d '{"chrom":"C_gar_LG28","start_bp":15000000,"end_bp":18000000}' \
  | python3 -m json.tool | head -5
```

Should return JSON, not ValidationError.
