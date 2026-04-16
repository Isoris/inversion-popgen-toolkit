#!/usr/bin/env python3
import gzip
import csv
from collections import defaultdict

BASE = "/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
INVDIR = f"{BASE}/inversion_localpca_v7"
OUTDIR = f"{INVDIR}/supplementary_tables"
CAND = f"{INVDIR}/06_mds_candidates/inversion_localpca.candidate_regions.tsv.gz"
OUT = f"{OUTDIR}/supp_candidate_master_table.tsv"

import os
os.makedirs(OUTDIR, exist_ok=True)

rows = []
with gzip.open(CAND, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for i, row in enumerate(reader, start=1):
        row["_genome_rank"] = i
        rows.append(row)

per_chr_count = defaultdict(int)
for row in rows:
    chrom = row["chrom"]
    per_chr_count[chrom] += 1
    row["_chr_rank"] = per_chr_count[chrom]

fieldnames = list(rows[0].keys()) if rows else []
# put ranks first
new_fields = ["_genome_rank", "_chr_rank"] + [x for x in fieldnames if x not in {"_genome_rank", "_chr_rank"}]

with open(OUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=new_fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: row.get(k, "") for k in new_fields})

print(f"[DONE] {OUT}")
print(f"[INFO] n_candidates={len(rows)}")
