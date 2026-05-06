#!/usr/bin/env python3
# =============================================================================
# 01_normalize_te/normalize_te_annotations.py
# =============================================================================
"""
Reads the TE manifest (one row per file) and emits, per species, a standardised
sorted BED file at:

    output/normalized_te/<species_id>.te.bed.gz

BED columns (9):
    chrom  start0  end  repeat_id  score  strand  repeat_class  repeat_family  source_file

Coordinate convention: 0-based half-open (BED). GFF3 input is 1-based inclusive
and is converted via gff_to_bed_coords().

What it understands:
  - edta_intact   (GFF3 with TE features + TSDs; we keep the TE features only,
                   and pass through target_site_duplication if present so the
                   downstream TSD layer can pick them up)
  - edta_anno     (GFF3 from EDTA TEanno; whole-element features)
  - gff3          (generic GFF3; keeps lines whose type or attributes look TE)
  - repeatmasker_out (.out fixed-column table; classic RM output)
  - bed           (passthrough if it already has at minimum chrom/start/end;
                   class/family parsed from a column if present)

Anything malformed is logged (file path, line number, reason) and skipped.

Usage:
    python normalize_te_annotations.py \
        --manifest config/te_file_manifest.tsv \
        --species  config/species_manifest.tsv \
        --out_dir  output/normalized_te \
        [--species_id C_gariepinus]   # restrict to one species
        [--keep_classes DNA,LTR,LINE,SINE,Helitron,TIR,MITE,satellite,low_complexity,unknown]
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    bucket_te_class,
    get_logger,
    gff_to_bed_coords,
    normalize_chrom,
    open_maybe_gz,
    read_tsv,
    write_run_manifest,
)


# -----------------------------------------------------------------------------
# parsers
# -----------------------------------------------------------------------------
GFF_ATTR_RE = re.compile(r"(\w+)=([^;]+)")


def _parse_gff_attrs(s: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in GFF_ATTR_RE.finditer(s)}


def parse_gff3(path: Path, source_label: str, log) -> list[tuple]:
    """Yield BED-row tuples for TE-like features in a GFF3."""
    out: list[tuple] = []
    n_skip = 0
    with open_maybe_gz(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                n_skip += 1
                continue
            chrom, _src, ftype, s1, e1, score, strand, _phase, attrs = parts[:9]
            try:
                s1i, e1i = int(s1), int(e1)
            except ValueError:
                n_skip += 1
                continue

            ftype_low = ftype.lower()
            attr = _parse_gff_attrs(attrs)
            cls_raw = attr.get("Classification") or attr.get("class") or ftype
            fam_raw = attr.get("Name") or attr.get("Family") or attr.get("ID") or ""
            tid     = attr.get("ID") or attr.get("Name") or f"{path.stem}_{ln}"

            # Skip purely structural lines we don't want as repeat intervals
            if ftype_low in ("region", "chromosome", "scaffold", "contig",
                             "gene", "mrna", "exon", "cds"):
                continue

            # Keep TE-like records and TSDs (TSDs are useful elsewhere; we tag class)
            keep = (
                "transposable_element" in ftype_low
                or "repeat_region"     in ftype_low
                or "ltr"               in ftype_low
                or "tir"               in ftype_low
                or "mite"              in ftype_low
                or "helitron"          in ftype_low
                or "line_element"      in ftype_low
                or "sine_element"      in ftype_low
                or "dna_transposon"    in ftype_low
                or "target_site_dup"   in ftype_low
                or "transposon"        in ftype_low
            )
            if not keep:
                continue

            bucket = bucket_te_class(cls_raw, fam_raw)
            s0, e0 = gff_to_bed_coords(s1i, e1i)
            out.append((
                normalize_chrom(chrom), s0, e0,
                str(tid), str(score) if score != "." else "0",
                strand if strand in ("+", "-") else ".",
                bucket, str(fam_raw) or "NA",
                source_label,
            ))
    if n_skip:
        log.info("  %s: skipped %d malformed lines", path, n_skip)
    return out


# RepeatMasker .out fixed-format columns (1-indexed in spec):
#  1=score  2=%div  3=%del  4=%ins  5=qname  6=qstart  7=qend  8=qleft  9=strand
#  10=match 11=class/family ...
RM_LINE_RE = re.compile(r"\s+")


def parse_repeatmasker_out(path: Path, source_label: str, log) -> list[tuple]:
    out: list[tuple] = []
    n_skip = 0
    with open_maybe_gz(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip():
                continue
            # skip header lines (non-numeric first token)
            tok0 = line.lstrip().split(None, 1)[0]
            if not tok0[0].isdigit():
                continue
            cols = RM_LINE_RE.split(line.strip())
            if len(cols) < 11:
                n_skip += 1
                continue
            try:
                score = cols[0]
                qname = cols[4]
                qs    = int(cols[5])
                qe    = int(cols[6])
                strd  = cols[8]
                match = cols[9]
                clsfam= cols[10]
            except (ValueError, IndexError):
                n_skip += 1
                continue

            cls_fam = clsfam.split("/")
            cls_raw = cls_fam[0] if cls_fam else clsfam
            fam_raw = cls_fam[1] if len(cls_fam) > 1 else "NA"
            bucket = bucket_te_class(cls_raw, fam_raw)
            # RM is 1-based inclusive
            s0, e0 = gff_to_bed_coords(qs, qe)
            strand_norm = "+" if strd == "+" else ("-" if strd in ("C", "-") else ".")
            out.append((
                normalize_chrom(qname), s0, e0,
                str(match), str(score), strand_norm,
                bucket, fam_raw, source_label,
            ))
    if n_skip:
        log.info("  %s: skipped %d malformed lines", path, n_skip)
    return out


def parse_bed(path: Path, source_label: str, log) -> list[tuple]:
    out: list[tuple] = []
    n_skip = 0
    with open_maybe_gz(path, "rt") as fh:
        for ln, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                n_skip += 1
                continue
            try:
                chrom = parts[0]
                s0    = int(parts[1])
                e0    = int(parts[2])
            except ValueError:
                n_skip += 1
                continue
            tid    = parts[3] if len(parts) > 3 else f"{path.stem}_{ln}"
            score  = parts[4] if len(parts) > 4 else "0"
            strand = parts[5] if len(parts) > 5 else "."
            cls_raw = parts[6] if len(parts) > 6 else ""
            fam_raw = parts[7] if len(parts) > 7 else "NA"
            bucket  = bucket_te_class(cls_raw, fam_raw)
            out.append((
                normalize_chrom(chrom), s0, e0,
                str(tid), str(score), strand,
                bucket, str(fam_raw), source_label,
            ))
    if n_skip:
        log.info("  %s: skipped %d malformed lines", path, n_skip)
    return out


PARSERS = {
    "edta_intact":      parse_gff3,
    "edta_anno":        parse_gff3,
    "gff3":             parse_gff3,
    "repeatmasker_out": parse_repeatmasker_out,
    "bed":              parse_bed,
}


# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--species", required=True, type=Path,
                    help="species manifest TSV (for assembly name validation)")
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--species_id", default=None,
                    help="restrict to one species_id")
    ap.add_argument("--keep_classes", default=None,
                    help="comma-list of TE classes to keep; default keeps all")
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("01_normalize_te", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species_rows = read_tsv(args.species,
                            required_cols=["species_id", "assembly_name"])
    known_species = {r["species_id"] for r in species_rows}

    manifest = read_tsv(args.manifest,
                        required_cols=["species_id", "file_path", "file_format"])

    # group files by species (priority order ascending = preferred first)
    by_sp: dict[str, list[dict]] = defaultdict(list)
    for r in manifest:
        if args.species_id and r["species_id"] != args.species_id:
            continue
        if r["species_id"] not in known_species:
            log.warning("species_id %s in TE manifest is not in species_manifest.tsv",
                        r["species_id"])
        try:
            r["_prio"] = int(r.get("priority") or 99)
        except ValueError:
            r["_prio"] = 99
        by_sp[r["species_id"]].append(r)

    keep_set = None
    if args.keep_classes:
        keep_set = {x.strip() for x in args.keep_classes.split(",") if x.strip()}

    for sp, rows in by_sp.items():
        rows.sort(key=lambda x: (x["_prio"], x["file_path"]))
        out_bed = args.out_dir / f"{sp}.te.bed.gz"
        log.info("normalising %d files for %s -> %s",
                 len(rows), sp, out_bed)

        all_recs: list[tuple] = []
        for r in rows:
            fmt = r["file_format"]
            parser = PARSERS.get(fmt)
            if parser is None:
                log.info("  skip (format=%s): %s", fmt, r["file_path"])
                continue
            try:
                recs = parser(Path(r["file_path"]), Path(r["file_path"]).name, log)
            except Exception as exc:
                log.error("  parser error %s on %s: %s",
                          fmt, r["file_path"], exc)
                continue
            log.info("  %s: %d records (%s)",
                     Path(r["file_path"]).name, len(recs), fmt)
            all_recs.extend(recs)

        if keep_set is not None:
            all_recs = [t for t in all_recs if t[6] in keep_set]

        # sort by chrom, start, end
        all_recs.sort(key=lambda t: (t[0], t[1], t[2]))

        with gzip.open(out_bed, "wt") as fh:
            for t in all_recs:
                fh.write("\t".join(str(x) for x in t) + "\n")

        write_run_manifest(out_bed,
                           species_id=sp,
                           n_input_files=len(rows),
                           n_records=len(all_recs))
        log.info("  -> %s (%d records)", out_bed, len(all_recs))

    log.info("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
