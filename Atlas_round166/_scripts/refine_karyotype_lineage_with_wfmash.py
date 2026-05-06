#!/usr/bin/env python3
"""
refine_karyotype_lineage_with_wfmash.py

Take an existing karyotype_lineage_v1.json (from
build_karyotype_lineage_v1_json.py) and refine each interesting cell
(class != '1-1', or any class if --refine-all) with a higher-resolution
wfmash pass. Writes a refined JSON with these new fields on each
classes_by_species[species_id] entry:

  refined_by_wfmash: 'confirmed' | 'refuted' | 'refined' | 'failed'
  wfmash_class:      class string after refinement (when present)
  wfmash_targets:    [{chrom, start_bp, end_bp, strand}, ...] (when present)
  wfmash_n_blocks:   integer
  wfmash_pct_identity: float (median across blocks)

The atlas page 16b reads these fields automatically — no atlas changes
needed once this output is loaded.

Honest limitation:
  This script orchestrates wfmash invocation + result classification.
  It does NOT discover new orthology — it only re-evaluates whether
  mashmap's chromosome-scale class holds at higher resolution.

Usage:
  python3 refine_karyotype_lineage_with_wfmash.py \\
    --input karyotype_lineage_v1.json \\
    --species-fasta-map species_fasta_map.tsv \\
    --output karyotype_lineage_v1.refined.json \\
    --threads 16 \\
    [--refine-all] \\
    [--wfmash-segment 50000] \\
    [--wfmash-block 250000] \\
    [--wfmash-identity 90]

species_fasta_map.tsv columns (whitespace separated):
  species_id<TAB>fasta_path

Example row:
  Cgar    /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Gar_LG.fa
  Cmac    /project/lt200308-agbsci/01-catfish_assembly/02-annot/fClaHyb_Mac_LG.fa
  Tros    /project/lt200308-agbsci/02-TE_catfish/00-GENOMES/T_rosablanca.fa

The script extracts the relevant focal-chr region from the focal species
fasta and aligns to the target species fasta with wfmash, then classifies
the result based on number of distinct target chromosomes hit.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True,
                   help='Input karyotype_lineage_v1.json')
    p.add_argument('--species-fasta-map', required=True,
                   help='TSV mapping species_id -> fasta path')
    p.add_argument('--output', required=True,
                   help='Output refined JSON')
    p.add_argument('--threads', type=int, default=16,
                   help='Threads for wfmash')
    p.add_argument('--refine-all', action='store_true',
                   help='Also refine 1-1 cells (default: only 1-2, 1-3, 1-4+)')
    p.add_argument('--wfmash-segment', type=int, default=50000,
                   help='wfmash -s segment length (bp). Default: 50000')
    p.add_argument('--wfmash-block', type=int, default=250000,
                   help='wfmash -l block length (bp). Default: 250000')
    p.add_argument('--wfmash-identity', type=int, default=90,
                   help='wfmash --pi identity percent. Default: 90')
    p.add_argument('--wfmash-binary', default='wfmash',
                   help='wfmash binary name/path. Default: wfmash')
    p.add_argument('--samtools-binary', default='samtools',
                   help='samtools binary name/path. Default: samtools')
    p.add_argument('--workdir', default=None,
                   help='Working dir for intermediate files. Default: tempdir.')
    p.add_argument('--keep-workdir', action='store_true',
                   help='Don\'t delete the workdir on exit (debugging)')
    return p.parse_args()


def load_species_fasta_map(path):
    m = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            m[parts[0]] = parts[1]
    return m


def need_refinement(cell, refine_all):
    cls = cell.get('class')
    if not cls:
        return False
    if refine_all:
        return True
    return cls in ('1-2', '1-3', '1-4+')


def extract_focal_region(focal_fasta, focal_chr, out_fa, samtools_bin):
    """Use samtools faidx to extract the focal chr to a temp fasta."""
    cmd = [samtools_bin, 'faidx', focal_fasta, focal_chr]
    with open(out_fa, 'w') as f:
        try:
            subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f'[warn] faidx failed for {focal_chr}: {e.stderr.decode()}\n')
            return False
    if os.path.getsize(out_fa) == 0:
        return False
    return True


def run_wfmash(query_fa, ref_fa, out_paf, threads, segment, block, identity,
               wfmash_bin):
    """Run wfmash; returns True on success and non-empty output."""
    cmd = [
        wfmash_bin,
        ref_fa, query_fa,
        '-t', str(threads),
        '-s', str(segment),
        '-l', str(block),
        '--pi', str(identity),
        '-N',  # no splitting, one-to-one
    ]
    try:
        with open(out_paf, 'w') as f:
            subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return False
    if os.path.getsize(out_paf) == 0:
        return False
    return True


def classify_paf(paf_path):
    """
    Read PAF; return dict:
      { 'targets': {target_chrom -> { 'start_bp': min, 'end_bp': max, 'strand': '+'|'-',
                                      'n_blocks': count, 'identity_sum': float }},
        'n_blocks_total': int,
        'mean_identity': float,
        'class': '1-1'|'1-2'|'1-3'|'1-4+'|'no_alignment' }
    PAF cols: 1 query, 2 qlen, 3 qstart, 4 qend, 5 strand, 6 ref, 7 reflen,
              8 rstart, 9 rend, 10 matches, 11 alnlen, 12 mapq, then tags.
    """
    targets = defaultdict(lambda: {'start_bp': None, 'end_bp': None,
                                   'strands': set(), 'n_blocks': 0,
                                   'identity_sum': 0.0, 'identity_n': 0})
    n_blocks = 0
    if not os.path.isfile(paf_path) or os.path.getsize(paf_path) == 0:
        return {'targets': {}, 'n_blocks_total': 0, 'mean_identity': 0.0,
                'class': 'no_alignment'}
    with open(paf_path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 12:
                continue
            strand = parts[4]
            ref = parts[5]
            try:
                rstart = int(parts[7])
                rend = int(parts[8])
                matches = int(parts[9])
                alnlen = int(parts[10])
            except ValueError:
                continue
            t = targets[ref]
            t['n_blocks'] += 1
            t['strands'].add(strand)
            t['start_bp'] = rstart if t['start_bp'] is None else min(t['start_bp'], rstart)
            t['end_bp']   = rend   if t['end_bp']   is None else max(t['end_bp'], rend)
            if alnlen > 0:
                t['identity_sum'] += matches / alnlen
                t['identity_n']   += 1
            n_blocks += 1
    n_targets = len(targets)
    cls = '1-1' if n_targets == 1 else \
          '1-2' if n_targets == 2 else \
          '1-3' if n_targets == 3 else \
          '1-4+' if n_targets >= 4 else 'no_alignment'
    mean_id = 0.0
    n_id = 0
    out_targets = {}
    for ref, t in targets.items():
        out_targets[ref] = {
            'chrom': ref,
            'start_bp': t['start_bp'],
            'end_bp': t['end_bp'],
            'strand': '+' if '+' in t['strands'] and '-' not in t['strands']
                      else ('-' if '-' in t['strands'] and '+' not in t['strands']
                            else '?'),
            'n_blocks': t['n_blocks'],
        }
        mean_id += t['identity_sum']
        n_id += t['identity_n']
    if n_id > 0:
        mean_id = (mean_id / n_id) * 100.0
    return {
        'targets': out_targets,
        'n_blocks_total': n_blocks,
        'mean_identity': mean_id,
        'class': cls,
    }


def classify_refinement_outcome(mashmap_class, wfmash_result):
    """Compare mashmap class with wfmash class to label the refinement state."""
    if wfmash_result['class'] == 'no_alignment' or wfmash_result['n_blocks_total'] == 0:
        return 'failed'
    if wfmash_result['class'] == mashmap_class:
        return 'confirmed'
    # Different class — was the mashmap call refuted (collapsed) or refined (different)?
    # We treat any change as 'refined' unless the new class is strictly less
    # fragmented than the mashmap class (i.e. 1-2 → 1-1, 1-3 → 1-1, 1-3 → 1-2),
    # in which case it's 'refuted'.
    rank = {'1-1': 1, '1-2': 2, '1-3': 3, '1-4+': 4}
    rm = rank.get(mashmap_class, 9)
    rw = rank.get(wfmash_result['class'], 9)
    if rw < rm:
        return 'refuted'
    return 'refined'


def main():
    args = parse_args()
    with open(args.input) as f:
        kl = json.load(f)
    species_map = load_species_fasta_map(args.species_fasta_map)
    if 'per_focal_chr' not in kl and 'per_cgar_chr' in kl:
        # Legacy schema — promote
        kl['per_focal_chr'] = []
        for e in kl['per_cgar_chr']:
            ee = dict(e)
            if 'cgar_chr' in ee and 'focal_chr' not in ee:
                ee['focal_chr'] = ee.pop('cgar_chr')
            kl['per_focal_chr'].append(ee)
    focal_species = kl.get('focal_species')
    if not focal_species:
        sys.exit('[error] input JSON has no focal_species; cannot run wfmash refinement '
                 '(don\'t know which fasta is the query). Re-emit with build_karyotype_lineage_v1_json.py '
                 'using --query-genome.')
    if focal_species not in species_map:
        sys.exit(f'[error] focal_species {focal_species} missing from species-fasta-map')
    focal_fasta = species_map[focal_species]
    if not os.path.isfile(focal_fasta):
        sys.exit(f'[error] focal fasta not found: {focal_fasta}')

    workdir = args.workdir or tempfile.mkdtemp(prefix='wfmash_refine_')
    os.makedirs(workdir, exist_ok=True)
    print(f'[info] workdir: {workdir}', file=sys.stderr)

    n_cells_to_refine = 0
    for entry in kl.get('per_focal_chr', []):
        cls = entry.get('classes_by_species', {})
        for sp, cell in cls.items():
            if sp == focal_species:
                continue
            if need_refinement(cell, args.refine_all):
                n_cells_to_refine += 1
    print(f'[info] {n_cells_to_refine} cells will be refined', file=sys.stderr)

    n_done = 0
    for entry in kl.get('per_focal_chr', []):
        focal_chr = entry.get('focal_chr')
        if not focal_chr:
            continue
        # Strip 'C_gar_' prefix etc. so samtools faidx finds it. But also try
        # the literal name first.
        focal_query_fa = os.path.join(workdir, f'focal_{focal_chr.replace("|", "_")}.fa')
        if not os.path.isfile(focal_query_fa):
            ok = extract_focal_region(focal_fasta, focal_chr, focal_query_fa, args.samtools_binary)
            if not ok:
                # Try short name (strip prefix like C_gar_)
                short = re.sub(r'^[A-Z]_[a-z]+_', '', focal_chr)
                if short != focal_chr:
                    ok = extract_focal_region(focal_fasta, short, focal_query_fa,
                                              args.samtools_binary)
            if not ok:
                print(f'[warn] could not extract focal region {focal_chr}; skipping',
                      file=sys.stderr)
                continue
        cls = entry.get('classes_by_species', {})
        for sp, cell in cls.items():
            if sp == focal_species:
                # Self mapping: trivially confirmed
                cell['refined_by_wfmash'] = 'confirmed'
                cell['wfmash_n_blocks'] = 1
                continue
            if not need_refinement(cell, args.refine_all):
                continue
            if sp not in species_map:
                print(f'[warn] species {sp} missing from species-fasta-map; skipping',
                      file=sys.stderr)
                cell['refined_by_wfmash'] = 'failed'
                continue
            ref_fa = species_map[sp]
            if not os.path.isfile(ref_fa):
                print(f'[warn] ref fasta missing for {sp}: {ref_fa}; skipping',
                      file=sys.stderr)
                cell['refined_by_wfmash'] = 'failed'
                continue
            paf = os.path.join(workdir, f'{focal_chr}__to__{sp}.paf')
            ok = run_wfmash(focal_query_fa, ref_fa, paf,
                            args.threads, args.wfmash_segment,
                            args.wfmash_block, args.wfmash_identity,
                            args.wfmash_binary)
            if not ok:
                cell['refined_by_wfmash'] = 'failed'
                n_done += 1
                continue
            res = classify_paf(paf)
            outcome = classify_refinement_outcome(cell['class'], res)
            cell['refined_by_wfmash'] = outcome
            if outcome != 'failed':
                cell['wfmash_class'] = res['class']
                cell['wfmash_n_blocks'] = res['n_blocks_total']
                cell['wfmash_pct_identity'] = round(res['mean_identity'], 2)
                # Sort wfmash targets by start_bp for stable output
                tgs = list(res['targets'].values())
                tgs.sort(key=lambda t: (t.get('chrom') or '', t.get('start_bp') or 0))
                cell['wfmash_targets'] = tgs
            n_done += 1
            if n_done % 10 == 0:
                print(f'[info] {n_done}/{n_cells_to_refine} cells refined',
                      file=sys.stderr)

    # Stamp refinement params
    kl.setdefault('params', {})
    kl['params']['wfmash_refinement'] = {
        'method': 'wfmash',
        'segment_bp': args.wfmash_segment,
        'block_bp': args.wfmash_block,
        'identity_pct': args.wfmash_identity,
        'filter': 'one-to-one',
        'applied_to_classes': ['ALL'] if args.refine_all else ['1-2', '1-3', '1-4+'],
        'refined_at': datetime.now(timezone.utc).isoformat(),
    }
    with open(args.output, 'w') as f:
        json.dump(kl, f, indent=2)
    print(f'[ok] wrote {args.output}', file=sys.stderr)
    print(f'[ok] {n_done} cells refined', file=sys.stderr)

    if not args.keep_workdir:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == '__main__':
    main()
