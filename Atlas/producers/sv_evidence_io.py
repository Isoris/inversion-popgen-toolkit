"""
producers/sv_evidence_io.py
============================
Common I/O for the per-candidate folder writers.

Layout convention (per Quentin's step-4 chat decision):

    data/
    └── <chrom>/
        └── candidates/
            ├── <cand_id>/
            │   ├── sv_genotype_counts.json        (STEP_SV_GT_AGG)
            │   ├── sv_evidence_combinations.json  (STEP_SV_EVID_COMB)
            │   └── sv_support_by_sample.json      (STEP_SV_SUPPORT)
            └── <cand_id>/...

The atlas-side drag-drop walks these folders via webkitGetAsEntry()
and dispatches each *.json by `format_version`. So the filenames inside
the folder are conventional but not strictly required — the schema
version drives routing.

Used by:
  - producers/STEP_SV_GT_AGG_aggregate_genotype_counts.py
  - producers/STEP_SV_EVID_COMB_emit_combinations.py
  - producers/STEP_SV_SUPPORT_emit_support_by_sample.py

This is intentionally a thin module — most of the heavy work is in
each STEP_*. The shared bits are folder layout + small validators
that mirror the atlas-side ones (so producer + consumer agree).
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Folder layout
# ---------------------------------------------------------------------------

@dataclass
class CandidateDir:
    """Resolved path to one candidate's per-candidate folder."""
    out_root: Path        # data/
    chrom:    str         # 'C_gar_LG28' or 'LG28' (caller's choice)
    cand_id:  str         # 'INV_LG28_003'

    @property
    def path(self) -> Path:
        return self.out_root / self.chrom / 'candidates' / self.cand_id

    def ensure(self) -> Path:
        self.path.mkdir(parents=True, exist_ok=True)
        return self.path

    def file(self, layer_name: str) -> Path:
        # layer_name is one of:
        #   'sv_genotype_counts'
        #   'sv_evidence_combinations'
        #   'sv_support_by_sample'
        #   'candidate'   (manifest, optional)
        return self.path / f'{layer_name}.json'


def write_layer(cand_dir: CandidateDir, layer_name: str, payload: dict) -> Path:
    """
    Write a single JSON layer into the per-candidate folder. Returns the
    written path. Validates `format_version` against `layer_name` so we
    don't silently put a combinations layer where gt_counts should be.
    """
    expected = {
        'sv_genotype_counts':       'sv_genotype_counts_v1',
        'sv_evidence_combinations': 'sv_evidence_combinations_v1',
        'sv_support_by_sample':     'sv_support_by_sample_v1',
        'candidate':                'candidate_manifest_v1',
    }.get(layer_name)
    if expected is None:
        raise ValueError(f"unknown layer name: {layer_name}")
    fv = payload.get('format_version')
    if fv != expected:
        raise ValueError(
            f"format_version mismatch for layer {layer_name!r}: "
            f"expected {expected!r}, got {fv!r}"
        )
    cand_dir.ensure()
    out = cand_dir.file(layer_name)
    with open(out, 'w') as f:
        # Compact JSON — atlas parses it fine, and per-candidate folders
        # add up across hundreds of candidates so size matters at scale.
        json.dump(payload, f, separators=(',', ':'))
    return out


def write_manifest(cand_dir: CandidateDir, **fields) -> Path:
    """
    Write a small candidate.json manifest summarising the candidate. Lets
    a reviewer answer "what is this folder?" without parsing every layer.
    Lightweight — boundaries, span, tier, validation status.
    """
    payload = {
        'format_version': 'candidate_manifest_v1',
        'candidate_id':   cand_dir.cand_id,
        'chrom':          cand_dir.chrom,
        **fields,
    }
    return write_layer(cand_dir, 'candidate', payload)


# ---------------------------------------------------------------------------
# Cohort metadata helpers
# ---------------------------------------------------------------------------

def load_locked_labels(path: str | Path) -> list[dict]:
    """
    Load the `state.candidateState[cid].locked_labels` JSON the atlas uses.
    Each entry: {'sample_id': str, 'label': 'HOMO_1' | 'HET' | 'HOMO_2'}.
    The atlas remaps HOMO_1/HOMO_2 to H1/H1 / H2/H2 conventionally.
    """
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and 'locked_labels' in data:
        return data['locked_labels']
    raise ValueError(f"unexpected locked_labels shape in {path}")


def label_to_group(label: str) -> str:
    """Map atlas karyotype label → group key used in JSON outputs."""
    return {'HOMO_1': 'H1/H1', 'HET': 'H1/H2', 'HOMO_2': 'H2/H2'}.get(label, label)


def split_samples_by_group(locks: list[dict]) -> dict[str, list[str]]:
    """{'H1/H1': [sample_ids], 'H1/H2': [...], 'H2/H2': [...]}"""
    out = {'H1/H1': [], 'H1/H2': [], 'H2/H2': []}
    for entry in locks:
        g = label_to_group(entry['label'])
        if g in out:
            out[g].append(entry['sample_id'])
    return out


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def print_summary(cand_dir: CandidateDir) -> None:
    """Print a small report for what's in the candidate folder."""
    p = cand_dir.path
    if not p.exists():
        print(f"  {cand_dir.cand_id}: folder does not exist yet")
        return
    print(f"  {cand_dir.cand_id}: {p}")
    for layer in ('sv_genotype_counts', 'sv_evidence_combinations',
                  'sv_support_by_sample', 'candidate'):
        f = cand_dir.file(layer)
        if f.exists():
            sz = f.stat().st_size
            print(f"    [✓] {layer}.json  ({sz:,} bytes)")
        else:
            print(f"    [ ] {layer}.json")
