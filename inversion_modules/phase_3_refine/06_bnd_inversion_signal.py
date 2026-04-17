#!/usr/bin/env python3
"""
06_bnd_inversion_signal.py — v2: dual-caller BND inversion signal
=================================================================
Extracts inversion-orientation BND records from:
  - DELLY2 BND catalog: uses INFO/CT=3to3 or CT=5to5
  - Manta RAW pre-conversion VCF: uses ALT bracket patterns for orientation
    (NOT the post-conversion BND catalog, which has NO inversion BNDs left)

CRITICAL: After convertInversion_py3.py, ALL intrachromosomal inversion-
orientation BNDs are converted to <INV>. The remaining Manta BND catalog
contains ONLY inter-chromosomal translocations. We must read the RAW
pre-conversion merged VCF for Manta inversion BND evidence.

DELLY2 evidence fields: INFO/PE, INFO/SR, INFO/CT, INFO/MAPQ
Manta evidence fields:  FORMAT/PR:ref,alt, FORMAT/SR:ref,alt (aggregated)
"""
import argparse, gzip, os, sys
from collections import Counter, defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--delly_bnd_vcf", default="")
    p.add_argument("--manta_raw_vcf", default="",
                   help="Manta RAW pre-conversion merged VCF")
    p.add_argument("--delly_inv_vcf", default="")
    p.add_argument("--manta_inv_vcf", default="")
    p.add_argument("--validated_tests", default="")
    p.add_argument("--pair_window", type=int, default=5000000)
    p.add_argument("--bp_match_window", type=int, default=1000)
    p.add_argument("--outdir", required=True)
    # BUGFIX 2026-04-17 (chat 5, FIX 29 v2): registry wiring for BND
    # rescue. Each paired junction (orphan or already-matched) writes an
    # existence_layer_b_bnd_rescue block to the evidence registry,
    # materialising q7b_bnd_rescued / q7b_bnd_rescue_source / etc. keys
    # that phase 4a Layer B scoring reads alongside the primary
    # existence_layer_b block from C00.
    p.add_argument("--registries_root", default="",
                   help="Path to registries/ root. Writes one "
                        "existence_layer_b_bnd_rescue block per paired "
                        "junction. Leave empty for standalone mode.")
    p.add_argument("--candidate_map", default="",
                   help="Optional TSV mapping chrom:bp1-bp2 → candidate_id "
                        "for aligning rescue records with phase_4's "
                        "canonical candidate IDs.")
    return p.parse_args()


def _load_registry(registries_root):
    """Locate and import the Python registry API. Silent no-op if unavailable."""
    if not registries_root:
        return None
    api_path = os.path.join(registries_root, "api", "python")
    if not os.path.isdir(api_path):
        print(f"[registry] API not found at {api_path} — BND rescue blocks skipped")
        return None
    if api_path not in sys.path:
        sys.path.insert(0, api_path)
    try:
        from registry_loader import load_registry  # type: ignore
    except ImportError as e:
        print(f"[registry] import failed: {e}")
        return None
    try:
        reg = load_registry(registries_root)
    except Exception as e:
        print(f"[registry] load_registry failed: {e}")
        return None
    print(f"[registry] loaded from {registries_root}")
    return reg


def _load_candidate_map_bnd(path):
    """
    Load TSV mapping (chrom, bp1, bp2) → candidate_id. Mapping is matched
    by reciprocal proximity (±10 kb on each breakpoint) at lookup time, so
    the TSV just needs to carry the canonical candidate's coordinates.
    Format: candidate_id<TAB>chrom<TAB>bp1<TAB>bp2
    """
    m = []
    if not path or not os.path.exists(path):
        return m
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        want = {"candidate_id", "chrom", "bp1", "bp2"}
        if not want.issubset(header):
            print(f"[registry] candidate_map {path} missing cols {want - set(header)} — ignored")
            return m
        ic = header.index("candidate_id")
        ih = header.index("chrom")
        ib1 = header.index("bp1")
        ib2 = header.index("bp2")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ic, ih, ib1, ib2):
                continue
            try:
                m.append({"cid": p[ic], "chrom": p[ih],
                          "bp1": int(p[ib1]), "bp2": int(p[ib2])})
            except ValueError:
                continue
    print(f"[registry] candidate_map: {len(m)} entries")
    return m


def _match_cid(cand_map, chrom, bp1, bp2, tol=10000):
    """Find best candidate_id for a (chrom, bp1, bp2) triple; fallback synthesises one."""
    for e in cand_map:
        if e["chrom"] != chrom:
            continue
        if abs(e["bp1"] - bp1) <= tol and abs(e["bp2"] - bp2) <= tol:
            return e["cid"]
    # Synthetic ID for standalone runs. Matches the inv_id convention of
    # phase_3 STEP01 output so joins with STEP03's Layer D writes.
    return f"{chrom}:{bp1}-{bp2}:bnd_rescue"

def extract_delly_bnd_inv_signal(vcf_path):
    """Extract CT=3to3/5to5 BNDs from DELLY. Uses INFO/PE and INFO/SR."""
    records, ct_counts, n_total = [], Counter(), 0
    if not vcf_path or not os.path.exists(vcf_path): return records, ct_counts, n_total
    with (gzip.open if vcf_path.endswith('.gz') else open)(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            n_total += 1
            p = line.strip().split('\t')
            info = dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
            flags = set(kv for kv in p[7].split(';') if '=' not in kv)
            ct = info.get('CT', '.'); ct_counts[ct] += 1
            if ct not in ('3to3', '5to5'): continue
            # Count carriers
            n_het = sum(1 for g in p[9:] if g.split(':')[0] in ('0/1','0|1','1|0'))
            n_hom = sum(1 for g in p[9:] if g.split(':')[0] in ('1/1','1|1'))
            records.append({
                'source':'delly', 'bnd_id':p[2], 'chrom':p[0], 'pos':int(p[1]),
                'chr2':info.get('CHR2',p[0]), 'pos2':int(info.get('END',info.get('POS2',p[1]))),
                'ct':ct, 'junction_type': 'left_3to3' if ct=='3to3' else 'right_5to5',
                'qual':float(p[5]) if p[5]!='.' else 0,
                'pe':int(info.get('PE',0)), 'sr':int(info.get('SR',0)),
                'precise':1 if 'PRECISE' in flags else 0,
                'filter':p[6], 'n_carriers':n_het+n_hom,
            })
    return records, ct_counts, n_total

def classify_bnd_orientation(alt):
    """Classify BND orientation from VCF ALT bracket pattern."""
    if ']' in alt:
        idx = alt.index(']')
        if idx == 0: return '3to3'  # ]p]t
        elif alt.endswith(']'): return '3to5' if alt[0] not in '[]' else '3to3'
        else: return '3to3'
    elif '[' in alt:
        idx = alt.index('[')
        if idx == 0: return '5to3'  # [p[t
        elif alt.endswith('['): return '5to5'  # t[p[
        else: return '5to3'
    return 'unknown'

def parse_bnd_mate_pos(alt):
    for delim in [']','[']:
        if delim in alt:
            for part in alt.replace('[',']').split(']'):
                if ':' in part and not part.startswith('<'):
                    cp = part.split(':')
                    if len(cp)==2:
                        try: return cp[0], int(cp[1])
                        except: pass
    return '', 0

def extract_manta_bnd_inv_signal(vcf_path):
    """Extract inversion-orientation BNDs from Manta RAW pre-conversion VCF.
    Uses ALT bracket patterns + aggregates FORMAT/PR,SR across carriers."""
    records, orient_counts, n_total = [], Counter(), 0
    if not vcf_path or not os.path.exists(vcf_path):
        print(f"  WARNING: Manta raw VCF not found: {vcf_path}"); return records, orient_counts, n_total
    with (gzip.open if vcf_path.endswith('.gz') else open)(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            n_total += 1
            p = line.strip().split('\t')
            info = dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
            if info.get('SVTYPE','') != 'BND': continue
            chr2, pos2 = parse_bnd_mate_pos(p[4])
            if not chr2 or chr2 != p[0]:
                orient_counts['inter_chrom' if chr2 and chr2!=p[0] else 'no_mate'] += 1; continue
            orient = classify_bnd_orientation(p[4]); orient_counts[orient] += 1
            if orient not in ('3to3','5to5'): continue
            # Aggregate PR/SR across carriers
            fmt = p[8].split(':') if len(p)>8 else []
            pri = fmt.index('PR') if 'PR' in fmt else -1
            sri = fmt.index('SR') if 'SR' in fmt else -1
            gi = fmt.index('GT') if 'GT' in fmt else 0
            spr=ssr=ncar=0
            for sd in p[9:]:
                sv=sd.split(':'); gt=sv[gi] if gi<len(sv) else './.'
                if gt in ('./.','0/0','./0','0/.','.','.','0|0'): continue
                ncar+=1
                if pri>=0 and pri<len(sv):
                    try: spr+=int(sv[pri].split(',')[1])
                    except: pass
                if sri>=0 and sri<len(sv):
                    try: ssr+=int(sv[sri].split(',')[1])
                    except: pass
            records.append({
                'source':'manta_raw', 'bnd_id':p[2], 'chrom':p[0], 'pos':int(p[1]),
                'chr2':chr2, 'pos2':pos2, 'ct':orient,
                'junction_type': f"left_{orient}" if orient=='3to3' else f"right_{orient}",
                'qual':float(p[5]) if p[5]!='.' else 0,
                'pe':spr, 'sr':ssr, 'precise':0, 'filter':p[6], 'n_carriers':ncar,
            })
    return records, orient_counts, n_total

def pair_junctions(records, max_dist):
    pairs=[]; by_chr=defaultdict(lambda:{'3to3':[],'5to5':[]})
    for r in records:
        if r['chrom']==r.get('chr2',r['chrom']): by_chr[r['chrom']][r['ct']].append(r)
    for ch, cts in by_chr.items():
        left=sorted(cts['3to3'],key=lambda x:x['pos']); right=sorted(cts['5to5'],key=lambda x:x['pos'])
        used=set()
        for lb in left:
            best=None; bd=max_dist+1
            for i,rb in enumerate(right):
                if i in used: continue
                d=rb['pos']-lb['pos']
                if 0<d<=max_dist and d<bd: bd=d; best=(i,rb)
            if best:
                used.add(best[0]); rb=best[1]
                pairs.append({'chrom':ch,'bp1':lb['pos'],'bp2':rb['pos'],'size':rb['pos']-lb['pos'],
                              'left_id':lb['bnd_id'],'right_id':rb['bnd_id'],
                              'left_src':lb['source'],'right_src':rb['source'],
                              'left_pe':lb['pe'],'right_pe':rb['pe'],'left_sr':lb['sr'],'right_sr':rb['sr']})
    return pairs

def load_inv_pos(vcf_path):
    pos=[]
    if not vcf_path or not os.path.exists(vcf_path): return pos
    with (gzip.open if vcf_path.endswith('.gz') else open)(vcf_path,'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            p=line.strip().split('\t')
            info=dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
            pos.append({'id':p[2],'chrom':p[0],'start':int(p[1]),'end':int(info.get('END',p[1]))})
    return pos

def write_tsv(rows, path):
    if not rows: return
    cols=list(rows[0].keys())
    with open(path,'w') as f:
        f.write('\t'.join(cols)+'\n')
        for r in rows: f.write('\t'.join(str(r.get(c,'.')) for c in cols)+'\n')
    print(f"  Written: {path} ({len(rows)} rows)")

def main():
    args=parse_args(); os.makedirs(args.outdir,exist_ok=True)
    print("="*60+"\nBND INVERSION SIGNAL ANALYSIS (v2: dual-caller)\n"+"="*60)

    reg = _load_registry(args.registries_root)
    cand_map = _load_candidate_map_bnd(args.candidate_map)

    print("\n--- DELLY BNDs (CT field) ---")
    db,dct,nd=extract_delly_bnd_inv_signal(args.delly_bnd_vcf)
    print(f"  Total: {nd}, CT dist: {dict(dct)}, inv-signal: {len(db)}")

    print("\n--- Manta BNDs (from RAW pre-conversion VCF, bracket patterns) ---")
    mb,mor,nm=extract_manta_bnd_inv_signal(args.manta_raw_vcf)
    print(f"  Total records: {nm}, orient dist: {dict(mor)}, inv-signal: {len(mb)}")

    write_tsv(db, os.path.join(args.outdir,"bnd_inv_signal_delly.tsv"))
    write_tsv(mb, os.path.join(args.outdir,"bnd_inv_signal_manta.tsv"))

    print("\n--- Pairing junctions ---")
    allb=db+mb; pairs=pair_junctions(allb,args.pair_window)
    print(f"  Paired: {len(pairs)}")
    write_tsv(pairs, os.path.join(args.outdir,"bnd_inv_paired.tsv"))

    print("\n--- Cross-ref with INV catalogs ---")
    dip=load_inv_pos(args.delly_inv_vcf); mip=load_inv_pos(args.manta_inv_vcf)
    aip=dip+mip; print(f"  DELLY INV: {len(dip)}, Manta INV: {len(mip)}")
    xref=[]
    for b in allb:
        mt='.'; mtype='no_inv_match'
        for inv in aip:
            if b['chrom']!=inv['chrom']: continue
            if abs(b['pos']-inv['start'])<=args.bp_match_window:
                mt=inv['id']; mtype='matches_BP1'; break
            elif abs(b['pos']-inv['end'])<=args.bp_match_window:
                mt=inv['id']; mtype='matches_BP2'; break
        xref.append({**b,'matched_inv':mt,'match_type':mtype})
    write_tsv(xref, os.path.join(args.outdir,"bnd_inv_crossref.tsv"))
    orphans=[r for r in xref if r['match_type']=='no_inv_match']
    write_tsv(orphans, os.path.join(args.outdir,"bnd_inv_rescued_candidates.tsv"))

    # ── Registry writes (FIX 29 v2) ─────────────────────────────────────
    # One existence_layer_b_bnd_rescue block per *paired* junction. Pairs
    # with no INV-catalog match are true orphan rescues (bnd_rescued=True,
    # match_type=no_inv_match). Pairs matching an INV call are written as
    # confirmation (bnd_rescued=False, match_type=matches_BP1/BP2) so that
    # downstream 4a/4e can still see the BND-level supporting evidence.
    n_blocks_written = 0
    if reg is not None:
        # Build pair-level xref: for each pair, inherit the match_type from
        # whichever junction was checked (left or right, whichever matched).
        pair_xref_left = {r['bnd_id']: r for r in xref}
        for pair in pairs:
            l_meta = pair_xref_left.get(pair['left_id'],  {'match_type': 'no_inv_match', 'matched_inv': '.'})
            r_meta = pair_xref_left.get(pair['right_id'], {'match_type': 'no_inv_match', 'matched_inv': '.'})
            # A pair confirms an INV if EITHER junction matches the catalog
            if l_meta['match_type'] != 'no_inv_match':
                match_type = l_meta['match_type']; matched_inv = l_meta['matched_inv']
            elif r_meta['match_type'] != 'no_inv_match':
                match_type = r_meta['match_type']; matched_inv = r_meta['matched_inv']
            else:
                match_type = 'no_inv_match'; matched_inv = ''

            # Determine provenance
            l_src = pair['left_src']; r_src = pair['right_src']
            if l_src == 'delly' and r_src == 'delly':
                rescue_source = 'delly_only'
            elif l_src == 'manta_raw' and r_src == 'manta_raw':
                rescue_source = 'manta_raw_only'
            else:
                rescue_source = 'cross_caller'

            cid = _match_cid(cand_map, pair['chrom'], pair['bp1'], pair['bp2'])

            block_data = {
                "bnd_rescued": (match_type == 'no_inv_match'),
                "rescue_source": rescue_source,
                "bnd_pair_bp1": int(pair['bp1']),
                "bnd_pair_bp2": int(pair['bp2']),
                "bnd_pair_size_bp": int(pair['size']),
                "left_junction_id":  pair['left_id'],
                "right_junction_id": pair['right_id'],
                "left_junction_source":  l_src,
                "right_junction_source": r_src,
                "left_pe":  int(pair.get('left_pe',  0)),
                "right_pe": int(pair.get('right_pe', 0)),
                "left_sr":  int(pair.get('left_sr',  0)),
                "right_sr": int(pair.get('right_sr', 0)),
                "matched_inv_id": matched_inv if matched_inv != '.' else '',
                "match_type": match_type,
            }
            try:
                reg.evidence.write_block(
                    candidate_id=cid,
                    block_type="existence_layer_b_bnd_rescue",
                    data=block_data,
                    source_script="phase_3_refine/06_bnd_inversion_signal.py",
                )
                n_blocks_written += 1
            except Exception as e:
                print(f"  [registry] BND-rescue block write failed for {cid}: {e}")

        print(f"[registry] wrote {n_blocks_written} existence_layer_b_bnd_rescue blocks")

    # Stats
    mc=Counter(r['match_type'] for r in xref)
    matched=sum(v for k,v in mc.items() if k!='no_inv_match')
    print(f"  Matched: {matched}, Orphan rescue: {mc.get('no_inv_match',0)}")
    with open(os.path.join(args.outdir,"bnd_inv_statistics.tsv"),'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"delly_total_bnds\t{nd}\n"); f.write(f"manta_raw_total\t{nm}\n")
        for ct,n in sorted(dct.items()): f.write(f"delly_CT_{ct}\t{n}\n")
        for o,n in sorted(mor.items()): f.write(f"manta_orient_{o}\t{n}\n")
        f.write(f"delly_inv_signal\t{len(db)}\nmanta_inv_signal\t{len(mb)}\n")
        f.write(f"paired_inversions\t{len(pairs)}\nmatched_to_inv\t{matched}\n")
        f.write(f"orphan_rescue\t{mc.get('no_inv_match',0)}\n")
        if reg is not None:
            f.write(f"registry_blocks_written\t{n_blocks_written}\n")
    print(f"\n{'='*60}\nCOMPLETE\n{'='*60}")

if __name__=="__main__": main()
