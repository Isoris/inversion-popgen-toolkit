# P1.3 — `/health` alias preservation

## Risk: trivial
## Lines changed: 0 in current code (already done) — adds a regression test
## Depends on: nothing
## Verification: both `curl http://127.0.0.1:8766/health` and `curl http://127.0.0.1:8766/api/health` return ok

---

## What's wrong

The atlas pings `/health` for the server-status badge. The server
originally only exposed `/api/health`. You patched
`popstats_server.py` around line 1138 to add a plain `/health`
alias that delegates to the canonical `/api/health` handler.

Risk: a future refactor of `popstats_server.py` could drop the
alias and silently break the badge again. Add a test so any drop
is caught.

## What's already in place

```python
# popstats_server.py around line 1138
@app.get("/api/health")
async def health():
    return {
        "ok": True,
        "service": "popstats_server",
        "n_samples": _ctx.n_samples,
        "beagle_dir": str(_ctx.beagle_dir),
        # ...
    }

@app.get("/health")
async def health_alias():
    return await health()
```

## Add regression test

Create `server_turn1/test_health_alias.py`:

```python
"""
Regression test for the /health -> /api/health alias.

The atlas's server-status badge probes /health. The server's canonical
endpoint is /api/health. If the alias is ever dropped, this test will
catch it.

Run with the server NOT running -- this is a static-import test, not
an integration test.
"""
from fastapi.testclient import TestClient
import os, sys

# Set config env so the app can initialize
os.environ.setdefault('POPSTATS_CONFIG',
    os.path.join(os.path.dirname(__file__), 'popstats_server.config.example.yaml'))

sys.path.insert(0, os.path.dirname(__file__))
from popstats_server import app

client = TestClient(app)

def test_canonical_health():
    r = client.get('/api/health')
    assert r.status_code == 200
    body = r.json()
    assert body.get('ok') is True
    assert body.get('service') == 'popstats_server'

def test_alias_health():
    r = client.get('/health')
    assert r.status_code == 200
    body = r.json()
    assert body.get('ok') is True

def test_both_return_same_shape():
    a = client.get('/api/health').json()
    b = client.get('/health').json()
    # Same keys (values may differ in cache stats etc.)
    assert set(a.keys()) == set(b.keys()), \
        f"alias diverged from canonical: a={set(a.keys())}, b={set(b.keys())}"

if __name__ == '__main__':
    test_canonical_health(); print('PASS canonical')
    test_alias_health();     print('PASS alias')
    test_both_return_same_shape(); print('PASS shape')
```

Run via:
```
cd Atlas/server_turn1
mamba activate atlas-popstats
pytest test_health_alias.py
# OR plain python:
python test_health_alias.py
```

## Document the alias

Add a one-line comment above the alias in `popstats_server.py`:

```python
# /health alias for the atlas server-status badge. The atlas pings
# /health (no /api prefix). If you drop this alias, the green dot in
# the yellow header folder will go grey. See P1_3 patch.
@app.get("/health")
async def health_alias():
    return await health()
```

## Verification

With server running:
```
curl -s http://127.0.0.1:8766/health        | python3 -m json.tool
curl -s http://127.0.0.1:8766/api/health    | python3 -m json.tool
```

Both should return identical JSON shapes with `"ok": true`.
