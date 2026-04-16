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
    return p.parse_args()

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
    print(f"\n{'='*60}\nCOMPLETE\n{'='*60}")

if __name__=="__main__": main()
