#!/usr/bin/env python3
"""
registry_loader.py — Python mirror of registry_loader.R

Minimal API subset for Python pipeline scripts (STEP03_statistical_tests.py,
breakpoint_evidence_audit.py). Reads and writes the same filesystem layout
as the R library so the two can interoperate.

Layer 1 (atomic) only. Queries and compute that need Engine B stay in R.
Python scripts should:
  - read sample groups from sample_registry (R wrote them)
  - write Tier-2 JSON blocks via write_block() (R reads them back)
  - NOT duplicate R's query/compute layers

Usage:
    from registry_loader import load_registry
    reg = load_registry()

    # Read groups registered by C01i (R)
    hom_inv = reg.samples.get_group("inv_LG12_17_HOM_INV")

    # Write Layer D block
    reg.evidence.write_block(
        candidate_id="LG12_17",
        block_type="existence_layer_d",
        data={
            "fisher_or": 23.7,
            "fisher_p": 4.98e-06,
            "contingency_table": [[45, 8], [12, 161]],
            "armitage_z": 5.12,
            "armitage_p": 3.01e-07,
            "n_inv_with_support": 45,
            "n_inv_total": 53,
            "samples_inv_with_support": ["CGA009", "CGA045"],
        },
    )
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional


def _now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _read_tsv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def _write_tsv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


# =============================================================================
# Samples
# =============================================================================
class SamplesAPI:
    def __init__(self, sample_dir: Path):
        self.sample_dir = sample_dir
        self.master_file = sample_dir / "sample_master.tsv"
        self.groups_file = sample_dir / "sample_groups.tsv"
        self.groups_dir = sample_dir / "groups"
        self.groups_dir.mkdir(parents=True, exist_ok=True)

    def get_master(self) -> List[Dict[str, str]]:
        return _read_tsv(self.master_file)

    def get_group(self, group_id: str) -> List[str]:
        mf = self.groups_dir / f"{group_id}.txt"
        if not mf.exists():
            return []
        return [line.strip() for line in mf.read_text().splitlines() if line.strip()]

    def has_group(self, group_id: str) -> bool:
        return (self.groups_dir / f"{group_id}.txt").exists()

    def list_groups(self) -> List[Dict[str, str]]:
        return _read_tsv(self.groups_file)

    def get_carriers(self, candidate_id: str, status: str = "HOM_INV") -> List[str]:
        """Convenience: get samples for inv_<candidate_id>_<status>."""
        return self.get_group(f"inv_{candidate_id}_{status}")

    def add_group(
        self,
        group_id: str,
        sample_ids: List[str],
        chrom: str = ".",
        inv_id: str = ".",
        subgroup: str = ".",
        description: str = "",
        overwrite: bool = False,
    ) -> bool:
        # Load existing
        groups = _read_tsv(self.groups_file)
        exists = any(g.get("group_id") == group_id for g in groups)
        if exists and not overwrite:
            print(f"[registry] Group exists: {group_id} — skip")
            return False
        if exists:
            groups = [g for g in groups if g.get("group_id") != group_id]

        # Write members file
        mf = self.groups_dir / f"{group_id}.txt"
        mf.write_text("\n".join(sample_ids) + "\n")

        # Append row
        row = {
            "group_id": group_id,
            "n": str(len(sample_ids)),
            "chrom": chrom,
            "inv_id": inv_id,
            "subgroup": subgroup,
            "description": description,
            "members_file": str(mf),
            "created": _now_iso(),
        }
        groups.append(row)
        _write_tsv(
            self.groups_file,
            groups,
            ["group_id", "n", "chrom", "inv_id", "subgroup",
             "description", "members_file", "created"],
        )
        return True


# =============================================================================
# Evidence
# =============================================================================
class EvidenceAPI:
    def __init__(self, evid_dir: Path, schema_dir: Path):
        self.evid_dir = evid_dir
        self.percand_dir = evid_dir / "per_candidate"
        self.global_dir = evid_dir / "global"
        self.schema_dir = schema_dir
        self.schema_blocks = schema_dir / "structured_block_schemas"

        for d in (self.percand_dir, self.global_dir):
            d.mkdir(parents=True, exist_ok=True)

        self._schema_cache: Dict[str, dict] = {}

    def _cand_dir(self, cid: str) -> Path:
        d = self.percand_dir / cid
        (d / "raw").mkdir(parents=True, exist_ok=True)
        (d / "structured").mkdir(parents=True, exist_ok=True)
        return d

    def _keys_file(self, cid: str) -> Path:
        return self._cand_dir(cid) / "keys.tsv"

    def _load_schema(self, block_type: str) -> Optional[dict]:
        if block_type in self._schema_cache:
            return self._schema_cache[block_type]
        sf = self.schema_blocks / f"{block_type}.schema.json"
        if not sf.exists():
            # Try templated boundary schema
            if block_type.startswith("boundary_"):
                sf = self.schema_blocks / "boundary.schema.json"
                if not sf.exists():
                    return None
            else:
                return None
        try:
            schema = json.loads(sf.read_text())
        except Exception:
            return None
        self._schema_cache[block_type] = schema
        return schema

    def _extract_keys_from_schema(
        self, schema: dict, block_type: str, data: dict
    ) -> List[tuple]:
        """Return list of (key_name, value) pairs per keys_extracted directive.

        Supports dotted paths in `from` (e.g., "gds_sub_block.gds_gap") and
        {side} templating for boundary_{left|right}.
        """
        kex = schema.get("keys_extracted") or []
        side = None
        if block_type.startswith("boundary_"):
            side = block_type.split("_", 1)[1]  # "left" or "right"

        out = []
        for ke in kex:
            key = ke.get("key")
            frm = ke.get("from")
            if not key or not frm:
                continue
            if side and "{side}" in key:
                key = key.replace("{side}", side)

            # Walk dotted path
            val = data
            for part in frm.split("."):
                if isinstance(val, dict):
                    val = val.get(part)
                else:
                    val = None
                    break
            if val is None:
                continue
            out.append((key, val))
        return out

    def write_block(
        self,
        candidate_id: str,
        block_type: str,
        data: dict,
        source_script: Optional[str] = None,
    ) -> dict:
        """Write a Tier-2 JSON block and extract flat keys per schema."""
        schema = self._load_schema(block_type)
        d = self._cand_dir(candidate_id)
        block_path = d / "structured" / f"{block_type}.json"

        block = {
            "block_type": block_type,
            "candidate_id": candidate_id,
            "source_script": source_script or os.environ.get(
                "CURRENT_SCRIPT", "python_script"
            ),
            "computed_at": _now_iso(),
            "data": data,
        }

        # Schema validation
        if schema is not None:
            required = schema.get("required") or []
            missing = [r for r in required if r not in data]
            if missing:
                print(
                    f"[registry] block {block_type} missing required: {missing} "
                    f"— writing anyway with incomplete flag"
                )
                block["validation_status"] = "incomplete"
            else:
                block["validation_status"] = "validated"
        else:
            block["validation_status"] = "no_schema"

        # Write
        block_path.write_text(json.dumps(block, indent=2, default=str))

        # Extract flat keys
        n_extracted = 0
        if schema is not None:
            extracted = self._extract_keys_from_schema(schema, block_type, data)
            if extracted:
                kf = self._keys_file(candidate_id)
                existing = _read_tsv(kf)
                extracted_keys = {k for k, _ in extracted}
                existing = [r for r in existing if r.get("key") not in extracted_keys]
                for key, val in extracted:
                    existing.append({
                        "key": key,
                        "value": str(val),
                        "source_script": block["source_script"],
                        "scope": "candidate",
                        "candidate_id": candidate_id,
                        "block_type": block_type,
                        "timestamp": block["computed_at"],
                    })
                    n_extracted += 1
                _write_tsv(
                    kf, existing,
                    ["key", "value", "source_script", "scope",
                     "candidate_id", "block_type", "timestamp"],
                )

        print(f"[registry] wrote block {block_type} for {candidate_id} "
              f"({n_extracted} keys extracted, status={block['validation_status']})")
        return {"path": str(block_path), "n_keys": n_extracted,
                "status": block["validation_status"]}

    def read_block(self, candidate_id: str, block_type: str) -> Optional[dict]:
        bp = self._cand_dir(candidate_id) / "structured" / f"{block_type}.json"
        if not bp.exists():
            return None
        return json.loads(bp.read_text())

    def get_keys(self, candidate_id: str, key: Optional[str] = None) -> List[dict]:
        kf = self._keys_file(candidate_id)
        rows = _read_tsv(kf)
        if key:
            rows = [r for r in rows if r.get("key") == key]
        return rows

    def has_evidence(self, candidate_id: str, key: str) -> bool:
        return len(self.get_keys(candidate_id, key)) > 0

    def add_evidence(
        self, candidate_id: str, key: str, value: str,
        script: str = "python", file_path: str = ""
    ) -> bool:
        """Append/replace a flat key (v9 compat)."""
        kf = self._keys_file(candidate_id)
        rows = _read_tsv(kf)
        rows = [r for r in rows if r.get("key") != key]
        rows.append({
            "key": key, "value": str(value),
            "source_script": script, "scope": "candidate",
            "candidate_id": candidate_id, "block_type": "",
            "timestamp": _now_iso(),
        })
        _write_tsv(
            kf, rows,
            ["key", "value", "source_script", "scope",
             "candidate_id", "block_type", "timestamp"],
        )
        return True

    def list_candidates(self) -> List[str]:
        if not self.percand_dir.exists():
            return []
        return sorted(p.name for p in self.percand_dir.iterdir() if p.is_dir())


# =============================================================================
# Intervals (minimal — Python scripts don't usually need much from here)
# =============================================================================
class IntervalsAPI:
    def __init__(self, interval_dir: Path):
        self.interval_dir = interval_dir
        self.cand_file = interval_dir / "candidate_intervals.tsv"

    def get_candidate(self, cid: str) -> Optional[Dict[str, str]]:
        for r in _read_tsv(self.cand_file):
            if r.get("candidate_id") == cid:
                return r
        return None

    def get_candidate_boundaries(self, cid: str) -> Optional[Dict[str, int]]:
        r = self.get_candidate(cid)
        if r is None:
            return None
        return {"left_bp": int(r["start_bp"]), "right_bp": int(r["end_bp"])}


# =============================================================================
# Top-level loader
# =============================================================================
@dataclass
class Registry:
    root: Path
    samples: SamplesAPI
    intervals: IntervalsAPI
    evidence: EvidenceAPI


def load_registry(registries_root: Optional[str] = None) -> Registry:
    if registries_root is None:
        registries_root = os.environ.get("REGISTRIES")
        if not registries_root:
            base = os.environ.get("BASE", ".")
            cand = Path(base) / "inversion-popgen-toolkit" / "registries"
            if cand.exists():
                registries_root = str(cand)
            else:
                registries_root = str(Path(base) / "sample_registry_v10_fallback")

    root = Path(registries_root)
    data_dir = root / "data"
    schema_dir = root / "schemas"

    for sub in ("sample_registry", "interval_registry",
                "evidence_registry/per_candidate",
                "evidence_registry/global"):
        (data_dir / sub).mkdir(parents=True, exist_ok=True)

    print(f"[registry_loader.py] root = {root}")

    return Registry(
        root=root,
        samples=SamplesAPI(data_dir / "sample_registry"),
        intervals=IntervalsAPI(data_dir / "interval_registry"),
        evidence=EvidenceAPI(data_dir / "evidence_registry", schema_dir),
    )


if __name__ == "__main__":
    # Smoke test
    reg = load_registry()
    print("samples groups:", len(reg.samples.list_groups()))
    print("candidates:", len(reg.evidence.list_candidates()))
