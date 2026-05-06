#!/usr/bin/env python3
"""
build_karyotype_lineage_v1_json.py

Convert mashmap_fusion_scan_lines pipeline output → karyotype_lineage_v1.json
for the atlas page 16b karyotype context layer.

Source pipeline:
  /project/lt200308-agbsci/01-catfish_assembly/05_ancestral_karyotype/mashmap_fusion_scan_lines/
  - mashmap_precurl_ancestral_events.slurm  (runs all-vs-all mashmap)
  - get_best_refs_and_1to2.sh               (picks best refs + 1-2 candidates)

Inputs:
  --summaries-dir   directory of {Q}__to__{R}.summary.tsv files
                    columns: query_chr  n_targets  class  targets
  --query-genome    the focal query genome basename
  --species-map     TSV mapping genome basename → atlas species id
  --outgroup        ONE OR MORE outgroup genome basenames (comma-separated)
  --sister          OPTIONAL sister-group genome basenames (comma-separated)
  --output          output .json path

The atlas is species-agnostic: focal_species / sister_species /
outgroup_species declared in the JSON drive the polarization rules.
This script reads --query-genome + --species-map to determine the focal
species, --outgroup and --sister to determine the others.

Honest limitation:
  This is chromosome-scale evidence (1 Mb / 85% identity), NOT
  breakpoint-resolution. The atlas displays it as a separate panel
  from the wfmash breakpoint-scale layer.

Example invocation on LANTA after re-running the all-vs-all with
fClaHyb_Gar_LG.fa as a query:

  python3 build_karyotype_lineage_v1_json.py \\
    --summaries-dir /project/lt200308-agbsci/01-catfish_assembly/05_ancestral_karyotype/mashmap_fusion_scan_lines/summaries \\
    --query-genome fClaHyb_Gar_LG.fa \\
    --sister fClaHyb_Mac_LG.fa \\
    --outgroup T_rosablanca.fa,N_graeffei.fa \\
    --species-map species_map.tsv \\
    --output karyotype_lineage_v1.json
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--summaries-dir', required=True,
                   help='Directory containing {Q}__to__{R}.summary.tsv files')
    p.add_argument('--query-genome', required=True,
                   help='Query genome basename (must match summary filenames). '
                        'This is the focal species; its species id from --species-map '
                        'goes into focal_species in the output JSON.')
    p.add_argument('--outgroup', required=True,
                   help='Outgroup genome basename(s), comma-separated '
                        '(e.g. "T_rosablanca.fa" or "T_rosablanca.fa,N_graeffei.fa")')
    p.add_argument('--sister', default='',
                   help='Sister-group genome basename(s), comma-separated. '
                        'Drives sister_lineage_fission verdicts. Optional.')
    p.add_argument('--species-map', required=True,
                   help='TSV mapping genome basename → atlas species id')
    p.add_argument('--output', required=True,
                   help='Output JSON path')
    p.add_argument('--focal-chr-prefix', default='C_gar_',
                   help='Prefix to add to query_chr names. Default: C_gar_ '
                        '(set this to match the focal species\'s chromosome '
                        'naming, e.g. C_mac_ when focal is Cmac)')
    p.add_argument('--segment-bp', type=int, default=1000000,
                   help='mashmap segment length in bp (for params block)')
    p.add_argument('--identity-pct', type=int, default=85,
                   help='mashmap identity threshold percent (for params block)')
    return p.parse_args()


def tagify(name: str) -> str:
    """Match the sanitization in get_best_refs_and_1to2.sh."""
    return re.sub(r'[^A-Za-z0-9._-]', '_', name)


def load_species_map(path: str) -> dict:
    """TSV: genome_basename<TAB>atlas_species_id"""
    m = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            m[parts[0]] = parts[1]
    return m


def parse_summary_file(path: str) -> dict:
    """
    Parse a summary file:
      query_chr  n_targets  class  targets
    Returns: { query_chr: { 'class': '1-1'|'1-2'|..., 'targets': [..] } }
    """
    rows = {}
    if not os.path.isfile(path):
        return rows
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            qchr, n_targets, cls, targets_str = parts[0], parts[1], parts[2], parts[3]
            targets = [t.strip() for t in targets_str.split(',') if t.strip()]
            rows[qchr] = {'class': cls, 'targets': targets}
    return rows


def comma_list(s: str):
    if not s:
        return []
    return [x.strip() for x in s.split(',') if x.strip()]


def main():
    args = parse_args()
    if not os.path.isdir(args.summaries_dir):
        sys.exit(f'[error] summaries dir not found: {args.summaries_dir}')
    if not os.path.isfile(args.species_map):
        sys.exit(f'[error] species map not found: {args.species_map}')

    species_map = load_species_map(args.species_map)
    if args.query_genome not in species_map:
        sys.exit(f'[error] query genome {args.query_genome} missing from species map')

    qtag = tagify(args.query_genome)
    focal_species = species_map[args.query_genome]

    outgroup_genomes = comma_list(args.outgroup)
    outgroup_species_list = []
    for og in outgroup_genomes:
        if og in species_map:
            outgroup_species_list.append(species_map[og])
        else:
            print(f'[warn] outgroup {og} not in species map; skipping',
                  file=sys.stderr)

    sister_genomes = comma_list(args.sister)
    sister_species_list = []
    for sg in sister_genomes:
        if sg in species_map:
            sister_species_list.append(species_map[sg])
        else:
            print(f'[warn] sister {sg} not in species map; skipping',
                  file=sys.stderr)

    # Iterate over all summaries with our query as the query side
    per_chr_classes = defaultdict(dict)  # focal_chr → species_id → {class, targets}
    n_pairs = 0
    for fn in sorted(os.listdir(args.summaries_dir)):
        if not fn.endswith('.summary.tsv'):
            continue
        m = re.match(r'^(.+?)__to__(.+?)\.summary\.tsv$', fn)
        if not m:
            continue
        q_in_file, r_in_file = m.group(1), m.group(2)
        if q_in_file != qtag:
            continue
        ref_basename = None
        for genome_basename in species_map:
            if tagify(genome_basename) == r_in_file:
                ref_basename = genome_basename
                break
        if not ref_basename:
            print(f'[warn] could not resolve ref tag {r_in_file} to a known genome',
                  file=sys.stderr)
            continue
        ref_species = species_map[ref_basename]
        if ref_species == focal_species:
            continue
        sumpath = os.path.join(args.summaries_dir, fn)
        rows = parse_summary_file(sumpath)
        n_pairs += 1
        for qchr, entry in rows.items():
            if not qchr.startswith(args.focal_chr_prefix):
                focal_chr = args.focal_chr_prefix + qchr.lstrip(':|/')
            else:
                focal_chr = qchr
            if '|' in focal_chr:
                focal_chr = args.focal_chr_prefix + focal_chr.split('|')[-1]
            per_chr_classes[focal_chr][ref_species] = entry

    if n_pairs == 0:
        sys.exit(f'[error] no summaries matched query {args.query_genome} in {args.summaries_dir}')

    # Add the focal species itself as 1-1 to its own chromosomes
    for focal_chr in list(per_chr_classes.keys()):
        per_chr_classes[focal_chr][focal_species] = {
            'class': '1-1',
            'targets': [focal_chr],
        }

    per_focal_chr = []
    for focal_chr in sorted(per_chr_classes.keys()):
        per_focal_chr.append({
            'focal_chr': focal_chr,
            'classes_by_species': per_chr_classes[focal_chr],
        })

    output = {
        'tool': 'karyotype_lineage_v1',
        'schema_version': 2,
        'generated_at': datetime.now(timezone.utc).isoformat(),
        'params': {
            'method': 'mashmap_one_to_one',
            'segment_bp': args.segment_bp,
            'identity_pct': args.identity_pct,
            'filter': 'one-to-one',
            'source_pipeline': args.summaries_dir,
            'query_genome': args.query_genome,
        },
        'focal_species': focal_species,
        'sister_species': sister_species_list,
        'outgroup_species': outgroup_species_list,
        'per_focal_chr': per_focal_chr,
    }

    with open(args.output, 'w') as f:
        json.dump(output, f, indent=2)

    print(f'[ok] wrote {args.output}', file=sys.stderr)
    print(f'[ok] {len(per_focal_chr)} {focal_species} chromosomes covered',
          file=sys.stderr)
    print(f'[ok] {n_pairs} (query, ref) pairs processed', file=sys.stderr)
    print(f'[ok] focal: {focal_species}', file=sys.stderr)
    print(f'[ok] sister: {sister_species_list}', file=sys.stderr)
    print(f'[ok] outgroup: {outgroup_species_list}', file=sys.stderr)


if __name__ == '__main__':
    main()
