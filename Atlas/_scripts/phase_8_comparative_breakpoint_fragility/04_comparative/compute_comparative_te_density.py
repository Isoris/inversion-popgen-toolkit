#!/usr/bin/env python3
# =============================================================================
# 04_comparative/compute_comparative_te_density.py
# =============================================================================
"""
Driver script: for every non-focal species that has both a TE BED
(output/normalized_te/<sp>.te.bed.gz) and a comparative window file
(output/breakpoint_windows/<sp>.windows.bed), invoke 03_density
to produce per-candidate density JSONs under
output/density/comparative/.

This is a thin wrapper. If you've already run 03 for the focal species, this
just calls it once per non-focal species.

Usage:
    python compute_comparative_te_density.py \
        --species_manifest config/species_manifest.tsv \
        --windows_dir output/breakpoint_windows \
        --te_dir output/normalized_te \
        --out_dir output/density/comparative
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

_THIS = Path(__file__).resolve()
sys.path.insert(0, str(_THIS.parent.parent))
from common import (  # noqa: E402
    get_logger,
    read_tsv,
    write_run_manifest,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species_manifest", required=True, type=Path)
    ap.add_argument("--windows_dir", required=True, type=Path)
    ap.add_argument("--te_dir", required=True, type=Path)
    ap.add_argument("--out_dir", required=True, type=Path)
    ap.add_argument("--log_dir", type=Path, default=Path("output/logs"))
    args = ap.parse_args()

    log = get_logger("04_compute_comparative", args.log_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    species = read_tsv(args.species_manifest,
                       required_cols=["species_id", "fai_path", "is_focal"])

    density_script = (_THIS.parent.parent / "03_density"
                      / "compute_te_density.py").resolve()

    n_run = n_skip = 0
    for r in species:
        sp = r["species_id"]
        if str(r["is_focal"]).lower() == "true":
            log.info("skip focal species %s", sp)
            continue

        win = args.windows_dir / f"{sp}.windows.bed"
        te  = args.te_dir / f"{sp}.te.bed.gz"
        fai = Path(r["fai_path"])

        if not win.exists():
            log.info("no windows file for %s (no synteny mapping?); skipping", sp)
            n_skip += 1
            continue
        if not te.exists():
            log.info("no normalised TE BED for %s; skipping", sp)
            n_skip += 1
            continue
        if not fai.exists():
            log.warning("no .fai for %s at %s; skipping", sp, fai)
            n_skip += 1
            continue

        cmd = [
            sys.executable, str(density_script),
            "--windows", str(win),
            "--te_bed",  str(te),
            "--fai",     str(fai),
            "--species_id", sp,
            "--out_dir", str(args.out_dir),
            "--log_dir", str(args.log_dir),
        ]
        log.info("running density for %s: %s", sp, " ".join(cmd))
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            log.error("density step exited rc=%d for %s", rc, sp)
        else:
            n_run += 1

    write_run_manifest(args.out_dir / "_comparative.run.json",
                       n_species_run=n_run,
                       n_species_skipped=n_skip)
    log.info("done: %d species processed, %d skipped", n_run, n_skip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
