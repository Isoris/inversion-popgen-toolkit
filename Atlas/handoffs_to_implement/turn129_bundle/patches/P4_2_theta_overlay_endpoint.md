# P4.2 — Theta candidate overlay endpoint

## Risk: medium
## Lines changed: ~100 server + ~40 atlas
## Depends on: P4.1 (same pattern, same plumbing)
## Verification: candidate panel shows per-sample θπ deviation tracks

---

## What you said

> Theta pi page status: There are two θπ concepts.
>
> 1. Full θπ scrubber page needs:
>    - theta_pi_per_window
>    - theta_pi_local_pca
>    - theta_pi_envelopes
>
> 2. Candidate θπ overlay can use:
>    - het_bridge_C_gar_LG28.sample_tP_windows.tsv.gz
>    - population_mean_tP_windows.C_gar_LG28.tsv.gz
>    - .pestPG-derived per-sample windows
>
> Do not assume theta results are missing. User says they already
> exist but locating/wiring them is slow.

This patch covers concept #2 — the candidate-level θπ overlay. The
full scrubber page (#1) needs cluster-side R scripts that haven't
shipped yet (turn 121 listed those layers as still pending).

## Step 1: locate the actual files

You wrote a search command. Run it first:

```bash
cd /mnt/e/01-catfish_assembly_manuscript_CGA
find . -type f \( \
  -name "het_bridge_C_gar_LG28.sample_tP_windows.tsv.gz" -o \
  -name "population_mean_tP_windows.C_gar_LG28.tsv.gz" -o \
  -name "*theta*LG28*.tsv.gz" -o \
  -name "*tP*LG28*.tsv.gz" -o \
  -name "*pestPG*" \
\) 2>/dev/null | tee atlas_theta_LG28_file_candidates.txt
```

Send me the output before this patch lands. The exact paths and
file shapes determine the endpoint reader.

## Expected file shapes (working assumption)

`het_bridge_C_gar_LG28.sample_tP_windows.tsv.gz`:
```
window_id  chrom  start  end  sample_id  theta_pi  n_sites
1          C_gar_LG28  0  10000  Cgar_001  0.00321  103
1          C_gar_LG28  0  10000  Cgar_002  0.00298  103
...
```

`population_mean_tP_windows.C_gar_LG28.tsv.gz`:
```
window_id  chrom  start  end  pop_mean_theta_pi  n_samples
1          C_gar_LG28  0  10000  0.00305  226
...
```

(These are guesses based on the filenames. Confirm with `zcat | head`.)

## Server endpoint

```python
@app.get("/api/theta/candidate")
async def theta_candidate(
    chrom: str,
    start: int,
    end: int,
    sample_ids: Optional[str] = None,  # comma-separated
):
    """Per-sample θπ in windows overlapping [start, end] for the
    candidate panel. Returns:
      - per-sample windows: each sample × window's θπ
      - population mean per window
      - z-score per sample × window vs population
    """
    if _ctx is None:
        return _server_not_initialized()

    base = _ctx.base
    sample_path = os.path.join(base, 'het_bridge',
                               f'het_bridge_{chrom}.sample_tP_windows.tsv.gz')
    popmean_path = os.path.join(base, 'het_bridge',
                                f'population_mean_tP_windows.{chrom}.tsv.gz')
    # !!! TBD: confirm actual paths from find output !!!

    if not os.path.isfile(sample_path):
        return {"error": f"per-sample theta file not found: {sample_path}"}
    if not os.path.isfile(popmean_path):
        return {"error": f"population mean theta file not found: {popmean_path}"}

    # Filter by chrom + window overlap with [start, end]
    # Returns:
    #   {
    #     chrom, start_bp, end_bp,
    #     windows: [{window_id, start, end, pop_mean, pop_n}],
    #     samples: [sample_id, ...],
    #     theta_matrix: [[s0_w0, s0_w1, ...], ...],  // n_samples × n_windows
    #     z_matrix:     [[z0_w0, z0_w1, ...], ...],
    #   }

    # ... (implementation follows P4.1 pattern: gzip + csv + filter)
    pass
```

This is a stub. **Don't implement until the file shapes are
confirmed** — I don't want to write code that reads a TSV with
columns we have to guess at.

## Atlas wiring (sketch)

Pattern parallel to P4.1: install a virtual layer when server is
available.

```js
function _atlasInstallServerThetaBridge() {
  if (!atlasServer || atlasServer.status !== 'available') return;
  // Mark per_sample_theta_pi as "available via server".
  state.layersPresent = state.layersPresent || {};
  state.layersPresent.per_sample_theta_pi = 'server';
}
```

The candidate panel renderer can then probe
`state.layersPresent.per_sample_theta_pi === 'server'` and call
`/api/theta/candidate`.

## Decision needed

Before this patch lands, send me:

1. The output of the `find` command above.
2. A `zcat <one-of-the-files> | head -3` for both file types.

Those two outputs let me write a correct endpoint instead of
guessing.

## Spec, not implementation

This is a placeholder spec. Treat it as:

- ✅ Endpoint name and shape are decided.
- ⏳ Reader implementation pending file inspection.
- ⏳ Atlas wiring pending the renderer's actual schema (need to
  inspect `state.data.theta_pi_*` accesses in your atlas).

When you ship the `find` output and a head of the files, I'll
expand this into a complete patch like P4.1.
