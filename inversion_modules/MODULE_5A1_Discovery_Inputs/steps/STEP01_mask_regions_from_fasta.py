#!/usr/bin/env python3
"""
STEP01_mask_regions_from_fasta.py

Extract BED intervals + per-contig stats for:
  1) hard-masked regions:            N (and n)
  2) soft-masked regions:            a/c/g/t
  3) normal uppercase regions:       A/C/G/T
  4) soft+hard masked regions:       a/c/g/t + N/n
  5) all-case ACGT regions:          A/C/G/T + a/c/g/t

Outputs (0-based, half-open BED):
  <prefix>.hardN.bed
  <prefix>.softacgt.bed
  <prefix>.normalACGT.bed
  <prefix>.masked_acgtN.bed
  <prefix>.allcaseACGTacgt.bed
  <prefix>.stats.tsv

Stats are reported per contig and for __TOTAL__.
For each mask track, the script reports:
  - total bp
  - percent of contig/genome
  - number of intervals
  - mean / median / min / max interval length
  - SD of interval lengths
  - N50 of interval lengths
"""

import argparse
import statistics
from dataclasses import dataclass, field
from typing import List, Optional, TextIO


@dataclass
class RunState:
    start: Optional[int] = None


@dataclass
class TrackStats:
    bp: int = 0
    n_intervals: int = 0
    max_len: int = 0
    lengths: List[int] = field(default_factory=list)

    def add_length(self, length: int) -> None:
        if length <= 0:
            return
        self.bp += length
        self.n_intervals += 1
        self.max_len = max(self.max_len, length)
        self.lengths.append(length)

    def merge_from(self, other: "TrackStats") -> None:
        self.bp += other.bp
        self.n_intervals += other.n_intervals
        self.max_len = max(self.max_len, other.max_len)
        self.lengths.extend(other.lengths)

    @staticmethod
    def _n50(lengths: List[int]) -> int:
        if not lengths:
            return 0
        total = sum(lengths)
        half = total / 2.0
        running = 0
        for x in sorted(lengths, reverse=True):
            running += x
            if running >= half:
                return x
        return 0

    def summary(self):
        if not self.lengths:
            return {
                "mean": 0.0,
                "median": 0.0,
                "min": 0,
                "max": 0,
                "sd": 0.0,
                "n50": 0,
            }

        mean_v = sum(self.lengths) / len(self.lengths)
        median_v = statistics.median(self.lengths)
        min_v = min(self.lengths)
        max_v = max(self.lengths)
        sd_v = statistics.stdev(self.lengths) if len(self.lengths) > 1 else 0.0
        n50_v = self._n50(self.lengths)

        return {
            "mean": mean_v,
            "median": median_v,
            "min": min_v,
            "max": max_v,
            "sd": sd_v,
            "n50": n50_v,
        }


def close_run(chrom: str, pos: int, state: RunState, out: TextIO, label: str, tstats: TrackStats):
    if state.start is None:
        return
    s = state.start
    e = pos
    if e > s:
        out.write(f"{chrom}\t{s}\t{e}\t{label}\n")
        tstats.add_length(e - s)
    state.start = None


def open_run(pos: int, state: RunState):
    if state.start is None:
        state.start = pos


def flush_contig(
    chrom,
    pos,
    hard_state,
    soft_state,
    norm_state,
    mask_state,
    allcase_state,
    hard_out,
    soft_out,
    norm_out,
    mask_out,
    allcase_out,
    hard_t,
    soft_t,
    norm_t,
    mask_t,
    allcase_t,
):
    close_run(chrom, pos, hard_state, hard_out, "hardN", hard_t)
    close_run(chrom, pos, soft_state, soft_out, "softacgt", soft_t)
    close_run(chrom, pos, norm_state, norm_out, "normalACGT", norm_t)
    close_run(chrom, pos, mask_state, mask_out, "masked_acgtN", mask_t)
    close_run(chrom, pos, allcase_state, allcase_out, "allcaseACGTacgt", allcase_t)


def fmt_float(x: float) -> str:
    return f"{x:.6f}"


def add_track_columns(prefix: str, header: List[str]) -> None:
    header.extend([
        f"bp_{prefix}",
        f"pct_{prefix}",
        f"nint_{prefix}",
        f"meanlen_{prefix}",
        f"medianlen_{prefix}",
        f"minlen_{prefix}",
        f"maxlen_{prefix}",
        f"sdlen_{prefix}",
        f"n50_{prefix}",
    ])


def append_track_values(length_bp: int, tstats: TrackStats, row: List[str]) -> None:
    pct = (100.0 * tstats.bp / length_bp) if length_bp else 0.0
    summ = tstats.summary()
    row.extend([
        str(tstats.bp),
        fmt_float(pct),
        str(tstats.n_intervals),
        fmt_float(summ["mean"]),
        fmt_float(float(summ["median"])),
        str(summ["min"]),
        str(summ["max"]),
        fmt_float(summ["sd"]),
        str(summ["n50"]),
    ])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--fasta", required=True, help="Input FASTA")
    ap.add_argument("-p", "--prefix", required=True, help="Output prefix")
    args = ap.parse_args()

    hard_bed_path = f"{args.prefix}.hardN.bed"
    soft_bed_path = f"{args.prefix}.softacgt.bed"
    norm_bed_path = f"{args.prefix}.normalACGT.bed"
    mask_bed_path = f"{args.prefix}.masked_acgtN.bed"
    allcase_bed_path = f"{args.prefix}.allcaseACGTacgt.bed"
    stats_path = f"{args.prefix}.stats.tsv"

    HARD = set("Nn")
    SOFT = set("acgt")
    NORM = set("ACGT")
    ALLCASE = set("ACGTacgt")

    with open(hard_bed_path, "w") as hard_out, \
         open(soft_bed_path, "w") as soft_out, \
         open(norm_bed_path, "w") as norm_out, \
         open(mask_bed_path, "w") as mask_out, \
         open(allcase_bed_path, "w") as allcase_out, \
         open(stats_path, "w") as stats_out:

        header = ["contig", "length_bp", "bp_other", "pct_other"]
        add_track_columns("hardN", header)
        add_track_columns("softacgt", header)
        add_track_columns("normalACGT", header)
        add_track_columns("masked_acgtN", header)
        add_track_columns("allcaseACGTacgt", header)
        stats_out.write("\t".join(header) + "\n")

        hard_state = RunState()
        soft_state = RunState()
        norm_state = RunState()
        mask_state = RunState()
        allcase_state = RunState()

        chrom = None
        pos = 0
        bp_hard = bp_soft = bp_norm = bp_mask = bp_allcase = bp_other = 0

        hard_t = TrackStats()
        soft_t = TrackStats()
        norm_t = TrackStats()
        mask_t = TrackStats()
        allcase_t = TrackStats()

        tot_len = 0
        tot_other = 0

        tot_hard_t = TrackStats()
        tot_soft_t = TrackStats()
        tot_norm_t = TrackStats()
        tot_mask_t = TrackStats()
        tot_allcase_t = TrackStats()

        def write_contig_stats(contig_name, length, bo, ht, st, nt, mt, at):
            if length == 0:
                return

            pct_other = (100.0 * bo / length) if length else 0.0
            row = [
                contig_name,
                str(length),
                str(bo),
                fmt_float(pct_other),
            ]
            append_track_values(length, ht, row)
            append_track_values(length, st, row)
            append_track_values(length, nt, row)
            append_track_values(length, mt, row)
            append_track_values(length, at, row)
            stats_out.write("\t".join(row) + "\n")

        with open(args.fasta, "r") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue

                if line.startswith(">"):
                    if chrom is not None:
                        flush_contig(
                            chrom,
                            pos,
                            hard_state,
                            soft_state,
                            norm_state,
                            mask_state,
                            allcase_state,
                            hard_out,
                            soft_out,
                            norm_out,
                            mask_out,
                            allcase_out,
                            hard_t,
                            soft_t,
                            norm_t,
                            mask_t,
                            allcase_t,
                        )
                        write_contig_stats(
                            chrom,
                            pos,
                            bp_other,
                            hard_t,
                            soft_t,
                            norm_t,
                            mask_t,
                            allcase_t,
                        )

                        tot_len += pos
                        tot_other += bp_other
                        tot_hard_t.merge_from(hard_t)
                        tot_soft_t.merge_from(soft_t)
                        tot_norm_t.merge_from(norm_t)
                        tot_mask_t.merge_from(mask_t)
                        tot_allcase_t.merge_from(allcase_t)

                    chrom = line[1:].split()[0]
                    pos = 0

                    hard_state = RunState()
                    soft_state = RunState()
                    norm_state = RunState()
                    mask_state = RunState()
                    allcase_state = RunState()

                    bp_hard = bp_soft = bp_norm = bp_mask = bp_allcase = bp_other = 0

                    hard_t = TrackStats()
                    soft_t = TrackStats()
                    norm_t = TrackStats()
                    mask_t = TrackStats()
                    allcase_t = TrackStats()
                    continue

                for c in line:
                    in_hard = c in HARD
                    in_soft = c in SOFT
                    in_norm = c in NORM
                    in_mask = in_hard or in_soft
                    in_allcase = c in ALLCASE

                    if in_hard:
                        bp_hard += 1
                    elif in_soft:
                        bp_soft += 1
                    elif in_norm:
                        bp_norm += 1
                    else:
                        bp_other += 1

                    if in_mask:
                        bp_mask += 1
                    if in_allcase:
                        bp_allcase += 1

                    if in_hard:
                        open_run(pos, hard_state)
                    else:
                        close_run(chrom, pos, hard_state, hard_out, "hardN", hard_t)

                    if in_soft:
                        open_run(pos, soft_state)
                    else:
                        close_run(chrom, pos, soft_state, soft_out, "softacgt", soft_t)

                    if in_norm:
                        open_run(pos, norm_state)
                    else:
                        close_run(chrom, pos, norm_state, norm_out, "normalACGT", norm_t)

                    if in_mask:
                        open_run(pos, mask_state)
                    else:
                        close_run(chrom, pos, mask_state, mask_out, "masked_acgtN", mask_t)

                    if in_allcase:
                        open_run(pos, allcase_state)
                    else:
                        close_run(chrom, pos, allcase_state, allcase_out, "allcaseACGTacgt", allcase_t)

                    pos += 1

        if chrom is not None:
            flush_contig(
                chrom,
                pos,
                hard_state,
                soft_state,
                norm_state,
                mask_state,
                allcase_state,
                hard_out,
                soft_out,
                norm_out,
                mask_out,
                allcase_out,
                hard_t,
                soft_t,
                norm_t,
                mask_t,
                allcase_t,
            )
            write_contig_stats(
                chrom,
                pos,
                bp_other,
                hard_t,
                soft_t,
                norm_t,
                mask_t,
                allcase_t,
            )

            tot_len += pos
            tot_other += bp_other
            tot_hard_t.merge_from(hard_t)
            tot_soft_t.merge_from(soft_t)
            tot_norm_t.merge_from(norm_t)
            tot_mask_t.merge_from(mask_t)
            tot_allcase_t.merge_from(allcase_t)

        pct_other_total = (100.0 * tot_other / tot_len) if tot_len else 0.0
        row = [
            "__TOTAL__",
            str(tot_len),
            str(tot_other),
            fmt_float(pct_other_total),
        ]
        append_track_values(tot_len, tot_hard_t, row)
        append_track_values(tot_len, tot_soft_t, row)
        append_track_values(tot_len, tot_norm_t, row)
        append_track_values(tot_len, tot_mask_t, row)
        append_track_values(tot_len, tot_allcase_t, row)
        stats_out.write("\t".join(row) + "\n")

    print("Wrote:")
    print(" ", hard_bed_path)
    print(" ", soft_bed_path)
    print(" ", norm_bed_path)
    print(" ", mask_bed_path)
    print(" ", allcase_bed_path)
    print(" ", stats_path)


if __name__ == "__main__":
    main()
