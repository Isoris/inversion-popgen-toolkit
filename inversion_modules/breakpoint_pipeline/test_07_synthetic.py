#!/usr/bin/env python3
# =============================================================================
# test_07_synthetic.py — smoke test for script 07 without needing real BAMs.
#
# Generates a synthetic per-sample evidence TSV with a known breakpoint signal:
#   - 10 INV samples (strong breakpoint evidence)
#   - 10 HET samples (weaker evidence)
#   - 10 REF samples (noise only)
# Then invokes 07_breakpoint_evidence_pileup.py in mode A to render.
#
# Use this to verify the script runs on your system + renders correctly
# before you trust it on real data.
#
# Usage:
#   python test_07_synthetic.py [--keep]
# =============================================================================

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--keep', action='store_true',
                     help='Keep temp files after test for inspection')
    ap.add_argument('--outdir', default=None,
                     help='Output directory (default: temp)')
    args = ap.parse_args()

    outdir = Path(args.outdir) if args.outdir else Path(tempfile.mkdtemp(prefix='test07_'))
    outdir.mkdir(parents=True, exist_ok=True)
    print(f'[test07] working in {outdir}')

    np.random.seed(42)

    # 30 synthetic samples: 10 INV + 10 HET + 10 REF
    samples = (
        [('CGA_INV_' + str(i).zfill(3), 'INV') for i in range(10)] +
        [('CGA_HET_' + str(i).zfill(3), 'HET') for i in range(10)] +
        [('CGA_REF_' + str(i).zfill(3), 'REF') for i in range(10)]
    )
    pd.DataFrame(samples, columns=['sample', 'group']).to_csv(
        outdir / 'samples.tsv', sep='\t', index=False
    )

    # Breakpoints
    bp_L, bp_R = 15115243, 18005891

    # Synthesize evidence
    rows = []
    for sid, grp in samples:
        n_ev = {'INV': (4, 6), 'HET': (1, 3), 'REF': (0, 1)}[grp]
        for bp_pos, bp_side, partner in [(bp_L, 'L', bp_R), (bp_R, 'R', bp_L)]:
            n = np.random.randint(n_ev[0], n_ev[1] + 1)
            for _ in range(n):
                jitter = np.random.randint(-150, 150)
                rs = bp_pos + jitter - np.random.randint(100, 200)
                re = bp_pos + jitter + np.random.randint(50, 150)
                rows.append(dict(
                    sample=sid, bp_side=bp_side, read_type='split_read',
                    read_pos=rs, read_end=re, strand='+', mate_pos=partner
                ))
                rows.append(dict(
                    sample=sid, bp_side=bp_side, read_type='discordant',
                    read_pos=bp_pos + np.random.randint(-500, -50),
                    read_end=0, strand='+',
                    mate_pos=bp_pos + np.random.randint(50, 500),
                    orient=np.random.choice(['FF', 'RR'])
                ))
                rows.append(dict(
                    sample=sid, bp_side=bp_side, read_type='soft_clip',
                    read_pos=bp_pos + np.random.randint(-50, 50),
                    clip_len=np.random.randint(15, 40),
                    side=np.random.choice(['left', 'right'])
                ))

            # Coverage with drop at breakpoint for carriers
            for bin_off in range(-2000, 2001, 50):
                depth = np.random.poisson(9)
                if abs(bin_off) < 100 and grp in ('INV', 'HET'):
                    depth = max(0, depth - np.random.randint(2, 5))
                rows.append(dict(
                    sample=sid, bp_side=bp_side, read_type='coverage',
                    read_pos=bp_pos + bin_off, depth=depth
                ))

    ev_path = outdir / 'evidence.tsv.gz'
    pd.DataFrame(rows).to_csv(ev_path, sep='\t', index=False, compression='gzip')
    print(f'[test07] wrote {len(rows)} evidence rows to {ev_path}')

    # Invoke the plot script
    script_dir = Path(__file__).parent
    script = script_dir / '07_breakpoint_evidence_pileup.py'
    out_pdf = outdir / 'test_output.pdf'

    cmd = [
        sys.executable, str(script),
        '--evidence-tsv', str(ev_path),
        '--sample-list',  str(outdir / 'samples.tsv'),
        '--chrom',        'C_gar_LG28',
        '--bp-left',      str(bp_L),
        '--bp-right',     str(bp_R),
        '--window',       '2000',
        '--max-samples',  '30',
        '--anonymize',
        '--out',          str(out_pdf),
    ]
    print(f'[test07] running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print('STDERR:', result.stderr)
        sys.exit(1)

    if out_pdf.exists():
        print(f'[test07] OK — {out_pdf} created')
        print(f'[test07]      and {out_pdf.with_suffix(".png")} preview')
    else:
        print('[test07] FAIL — output not found')
        sys.exit(1)

    if not args.keep:
        print(f'[test07] to inspect outputs, re-run with --keep and --outdir <path>')


if __name__ == '__main__':
    main()
