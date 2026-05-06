# P1.5 — Server startup `POPSTATS_CONFIG` race

## Risk: low
## Lines changed: ~20
## Depends on: nothing
## Verification: `python popstats_server.py --config <path>` always works first time

---

## What's wrong

You hit this:
```
"server not initialized — set POPSTATS_CONFIG env var
 or run via python popstats_server.py --config <path>"
```

…sometimes, and then it worked after you set the env var manually.

Root cause: when the server is launched via
`python popstats_server.py --config X`, the `__main__` block parses
`--config`, sets context, then calls `uvicorn.run("popstats_server:app", ...)`.
But uvicorn's reload mode (or its worker spawning) can re-import the
module fresh, at which point the `--config` argv is gone and only
`POPSTATS_CONFIG` env var survives.

So the worker's import of `popstats_server` fails the context check,
and the app responds with the not-initialized error until you also
set `POPSTATS_CONFIG`.

## The fix — set env var BEFORE uvicorn imports the module

In `popstats_server.py`, anchor on the `if __name__ == "__main__":`
block.

### Before (buggy)

```python
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--host", default=None)
    p.add_argument("--port", default=None, type=int)
    args = p.parse_args()

    _init_from_config(args.config)  # populates _ctx in this process

    import uvicorn
    uvicorn.run("popstats_server:app",
                host=args.host or _ctx.host,
                port=args.port or _ctx.port,
                reload=False)
```

The bug: `_init_from_config` ran in the parent process, but uvicorn's
spawn re-imports `popstats_server` in the child without that context.

### After (race-free)

```python
if __name__ == "__main__":
    import argparse, os
    p = argparse.ArgumentParser()
    p.add_argument("--config", required=True)
    p.add_argument("--host", default=None)
    p.add_argument("--port", default=None, type=int)
    args = p.parse_args()

    # Set the env var BEFORE uvicorn imports the module. Workers and
    # reloads will all see this env var and re-initialize from the
    # same config file. This eliminates the parent-vs-child race.
    cfg_abs = os.path.abspath(args.config)
    os.environ["POPSTATS_CONFIG"] = cfg_abs

    # Initialize for the current process too (so /api/health works
    # immediately even before the first uvicorn import).
    _init_from_config(cfg_abs)

    # Resolve host/port from CLI override → ctx → defaults.
    host = args.host or _ctx.host or "127.0.0.1"
    port = args.port or _ctx.port or 8766

    import uvicorn
    uvicorn.run("popstats_server:app",
                host=host, port=port,
                reload=False)
```

### Add a lazy initializer at module top level

Anchor: top of `popstats_server.py`, just after the imports but
before any `@app.get` decorator runs.

```python
# turn 129 P1.5: lazy init — if POPSTATS_CONFIG is set in env (e.g.
# from a uvicorn worker re-import), initialize ctx automatically
# so endpoints don't return "server not initialized".
import os as _os
_AUTO_CFG = _os.environ.get("POPSTATS_CONFIG")
if _AUTO_CFG and _ctx is None:
    try:
        _init_from_config(_AUTO_CFG)
    except Exception as e:
        # Don't crash on import — let the not-initialized path handle.
        import sys as _sys
        print(f"[popstats_server] warning: auto-init from "
              f"POPSTATS_CONFIG={_AUTO_CFG} failed: {e}",
              file=_sys.stderr)
```

(The exact name `_ctx` may differ in your code — check what
`_init_from_config` populates. The pattern is: at module top, if
the env var is set and the global state isn't yet populated, try
to populate it.)

### Improve the not-initialized error message

If init still fails after both paths, the error message should tell
the user what to do. Find:

```python
return {"error": "server not initialized — ..."}
```

and change to:

```python
return {
    "error": "server not initialized",
    "hint": (
        "Set POPSTATS_CONFIG env var to the absolute path of "
        "popstats_server.local.yaml, or run via "
        "`python popstats_server.py --config <path>`. "
        "Current POPSTATS_CONFIG: " + str(os.environ.get('POPSTATS_CONFIG'))
    ),
}
```

This way the error self-documents.

## Verification

```bash
cd Atlas/server_turn1
mamba activate atlas-popstats

# Cold start (no env var):
unset POPSTATS_CONFIG
python popstats_server.py --config popstats_server.local.yaml &
sleep 2
curl -s http://127.0.0.1:8766/api/health | python3 -m json.tool
kill %1

# Should return ok=true, n_samples=226 etc.

# Restart with env var only (no --config flag):
export POPSTATS_CONFIG="$(pwd)/popstats_server.local.yaml"
python -c "import uvicorn; uvicorn.run('popstats_server:app', host='127.0.0.1', port=8766)" &
sleep 2
curl -s http://127.0.0.1:8766/api/health | python3 -m json.tool
kill %1

# Should also return ok=true.
```

Both invocation styles should work without
"server not initialized".

## Test (optional)

`server_turn1/test_startup_paths.py`:

```python
import os, subprocess, time, signal, json, requests

ROOT = os.path.dirname(__file__)
CFG = os.path.join(ROOT, 'popstats_server.local.yaml')

def _wait_health(url, timeout=5):
    for _ in range(timeout * 5):
        try:
            r = requests.get(url + '/api/health', timeout=0.5)
            if r.ok and r.json().get('ok'):
                return True
        except Exception:
            pass
        time.sleep(0.2)
    return False

def test_via_cli_flag():
    env = dict(os.environ); env.pop('POPSTATS_CONFIG', None)
    p = subprocess.Popen(
        ['python', os.path.join(ROOT, 'popstats_server.py'),
         '--config', CFG],
        env=env)
    try:
        assert _wait_health('http://127.0.0.1:8766')
    finally:
        p.send_signal(signal.SIGTERM); p.wait(timeout=3)

def test_via_env_var():
    env = dict(os.environ); env['POPSTATS_CONFIG'] = CFG
    p = subprocess.Popen(
        ['python', '-c',
         "import uvicorn; uvicorn.run('popstats_server:app', host='127.0.0.1', port=8766)"],
        env=env, cwd=ROOT)
    try:
        assert _wait_health('http://127.0.0.1:8766')
    finally:
        p.send_signal(signal.SIGTERM); p.wait(timeout=3)
```
