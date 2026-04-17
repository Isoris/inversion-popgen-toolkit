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

    # ══════════════ chat 11 extensions ══════════════════════════════════
    # Atomic reads over sample_master + sample_groups + groups/<gid>.txt.
    # Mirrors R API method-for-method (same names, same semantics) so
    # Python consumers (STEP03, breakpoint_evidence_audit) can use the
    # same vocabulary as R pipeline code.

    def _read_members(self, group_id: str) -> List[str]:
        """Read a group's member file from the groups/ subdir."""
        mf = self.groups_dir / f"{group_id}.txt"
        if not mf.exists():
            return []
        return [line.strip() for line in mf.read_text().splitlines() if line.strip()]

    _fam_warned: bool = False

    def _check_family_col(self, rows: List[Dict[str, str]]) -> bool:
        """Return True if sample_master has a family_id column; warn-once if not."""
        if not rows:
            return False
        has = "family_id" in rows[0]
        if not has and not SamplesAPI._fam_warned:
            print("[registry] sample_master.tsv has no 'family_id' column — "
                  "family methods return None/[].")
            SamplesAPI._fam_warned = True
        return has

    def get_sample_metadata(self, sid: str) -> Optional[Dict[str, str]]:
        """Return the sample_master row for sid, or None."""
        for r in self.get_master():
            if r.get("sample_id") == sid:
                return r
        return None

    def get_sample_groups(self, sid: str) -> List[str]:
        """Return all group_ids containing sid."""
        groups = self.list_groups()
        hits = []
        for g in groups:
            gid = g.get("group_id")
            if gid and sid in self._read_members(gid):
                hits.append(gid)
        return hits

    def get_family(self, sid: str) -> Optional[str]:
        """Return family_id for sid, or None if no family_id column / unknown sample."""
        rows = self.get_master()
        if not self._check_family_col(rows):
            return None
        for r in rows:
            if r.get("sample_id") == sid:
                f = r.get("family_id")
                return f if f else None
        return None

    def list_families(self) -> List[str]:
        """Return sorted distinct non-empty family IDs."""
        rows = self.get_master()
        if not self._check_family_col(rows):
            return []
        fams = {r.get("family_id") for r in rows if r.get("family_id")}
        return sorted(f for f in fams if f)

    def get_family_members(self, fam_id: str) -> List[str]:
        """Return sample IDs belonging to fam_id."""
        rows = self.get_master()
        if not self._check_family_col(rows):
            return []
        return [r["sample_id"] for r in rows if r.get("family_id") == fam_id]

    def list_carriers_across_candidates(self, status: str = "HOM_INV") -> List[str]:
        """Return unique sample IDs that are `status` for any candidate.

        Uses the 'subgroup' column if present, else falls back to
        inv_<cid>_<status> name-suffix matching.
        """
        groups = self.list_groups()
        gids: List[str] = []
        for g in groups:
            gid = g.get("group_id", "")
            sub = g.get("subgroup")
            if sub:
                if sub == status:
                    gids.append(gid)
            else:
                if gid.endswith(f"_{status}"):
                    gids.append(gid)
        samples: set = set()
        for gid in gids:
            samples.update(self._read_members(gid))
        return sorted(samples)

    def compute_group_overlap(self, gid1: str, gid2: str) -> Dict[str, Any]:
        """Return overlap stats for two groups: intersection, union, jaccard, set-diffs."""
        s1 = set(self._read_members(gid1))
        s2 = set(self._read_members(gid2))
        inter = s1 & s2
        uni = s1 | s2
        return {
            "gid1": gid1, "gid2": gid2,
            "n1": len(s1), "n2": len(s2),
            "intersection": len(inter),
            "union": len(uni),
            "jaccard": (len(inter) / len(uni)) if uni else None,
            "only_in_1": sorted(s1 - s2),
            "only_in_2": sorted(s2 - s1),
            "shared": sorted(inter),
        }

    def list_recombinant_candidates(
        self, sid: str, include_subclasses: bool = True
    ) -> List[str]:
        """Return CIDs where sid is RECOMBINANT (optionally including GC/DCO subclasses)."""
        import re
        groups = self.get_sample_groups(sid)
        if include_subclasses:
            pat = re.compile(r"^inv_(.+)_RECOMBINANT(?:_GC|_DCO)?$")
        else:
            pat = re.compile(r"^inv_(.+)_RECOMBINANT$")
        cids: List[str] = []
        for g in groups:
            m = pat.match(g)
            if m:
                cids.append(m.group(1))
        # dedupe preserving order
        seen: set = set()
        return [c for c in cids if not (c in seen or seen.add(c))]

    def get_groups_for_candidate(self, cid: str) -> Dict[str, List[str]]:
        """Return all sample groups for candidate cid, keyed by status.

        Each slot is [] if that group isn't registered.
        """
        statuses = ["HOM_REF", "HET", "HOM_INV", "RECOMBINANT",
                    "RECOMBINANT_GC", "RECOMBINANT_DCO", "HOM_STD"]
        return {
            st: (self.get_group(f"inv_{cid}_{st}") if self.has_group(f"inv_{cid}_{st}")
                 else [])
            for st in statuses
        }

    def find_co_segregating_groups(
        self,
        min_jaccard: float = 0.9,
        group_prefix: Optional[str] = "inv_",
        min_size: int = 3,
    ) -> List[Dict[str, Any]]:
        """Return pairs of groups whose membership overlap >= min_jaccard.

        Restricts to groups matching `group_prefix` and size >= min_size.
        Returns a list of dicts with keys: gid1, gid2, n1, n2, intersection, jaccard.
        """
        groups = self.list_groups()
        gids = [g["group_id"] for g in groups if "group_id" in g]
        if group_prefix:
            gids = [g for g in gids if g.startswith(group_prefix)]
        # Preload members
        mem = {g: set(self._read_members(g)) for g in gids}
        gids = [g for g in gids if len(mem[g]) >= min_size]
        out: List[Dict[str, Any]] = []
        for i in range(len(gids) - 1):
            s1 = mem[gids[i]]
            for j in range(i + 1, len(gids)):
                s2 = mem[gids[j]]
                uni = s1 | s2
                if not uni:
                    continue
                inter = s1 & s2
                jac = len(inter) / len(uni)
                if jac >= min_jaccard:
                    out.append({
                        "gid1": gids[i], "gid2": gids[j],
                        "n1": len(s1), "n2": len(s2),
                        "intersection": len(inter),
                        "jaccard": round(jac, 4),
                    })
        return out


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
# Intervals
# =============================================================================
class IntervalsAPI:
    def __init__(self, interval_dir: Path):
        self.interval_dir = interval_dir
        self.cand_file = interval_dir / "candidate_intervals.tsv"

    def _read_cands(self) -> List[Dict[str, str]]:
        """Read candidate_intervals.tsv; empty list if missing."""
        return _read_tsv(self.cand_file)

    def get_candidate(self, cid: str) -> Optional[Dict[str, str]]:
        for r in self._read_cands():
            if r.get("candidate_id") == cid:
                return r
        return None

    def get_candidate_boundaries(self, cid: str) -> Optional[Dict[str, int]]:
        r = self.get_candidate(cid)
        if r is None:
            return None
        return {"left_bp": int(r["start_bp"]), "right_bp": int(r["end_bp"])}

    # ══════════════ chat 11 extensions ══════════════════════════════════
    # Atomic mirror of the R interval API. Semantics identical. parent_id
    # may be "", ".", or missing — all treated as "no parent".

    @staticmethod
    def _has_parent(row: Dict[str, str]) -> bool:
        p = row.get("parent_id", "")
        return bool(p) and p != "." and p != "NA"

    def get_children(self, cid: str) -> List[str]:
        """CIDs whose parent_id == cid."""
        return [r["candidate_id"] for r in self._read_cands()
                if self._has_parent(r) and r.get("parent_id") == cid]

    def get_parent(self, cid: str) -> Optional[str]:
        """Direct parent CID, or None."""
        r = self.get_candidate(cid)
        if r is None or not self._has_parent(r):
            return None
        return r["parent_id"]

    def get_ancestors(self, cid: str) -> List[str]:
        """Walk up parent_id chain, closest first. Cycle-safe."""
        out: List[str] = []
        seen = {cid}
        cur = cid
        while True:
            r = self.get_candidate(cur)
            if r is None or not self._has_parent(r):
                break
            p = r["parent_id"]
            if p in seen:
                break
            out.append(p)
            seen.add(p)
            cur = p
        return out

    def get_descendants(self, cid: str, depth: Optional[int] = None) -> List[str]:
        """BFS down parent_id relation. depth=None means full subtree; 1=direct children."""
        rows = self._read_cands()
        out: List[str] = []
        frontier = [cid]
        level = 0
        max_depth = depth if depth is not None else float("inf")
        out_set = {cid}
        while frontier and level < max_depth:
            next_f: List[str] = []
            for r in rows:
                if self._has_parent(r) and r.get("parent_id") in frontier:
                    ccid = r["candidate_id"]
                    if ccid not in out_set:
                        next_f.append(ccid)
                        out.append(ccid)
                        out_set.add(ccid)
            if not next_f:
                break
            frontier = next_f
            level += 1
        return out

    def get_overlapping(
        self,
        chrom: str,
        start_bp: int,
        end_bp: int,
        exclude_cid: Optional[str] = None,
    ) -> List[str]:
        """CIDs on chrom whose interval overlaps [start_bp, end_bp]."""
        out = []
        for r in self._read_cands():
            if r.get("chrom") != chrom:
                continue
            try:
                s = int(r["start_bp"])
                e = int(r["end_bp"])
            except (KeyError, ValueError):
                continue
            if s > end_bp or e < start_bp:
                continue
            if exclude_cid is not None and r.get("candidate_id") == exclude_cid:
                continue
            out.append(r["candidate_id"])
        return out

    def get_nested_within(self, cid: str) -> List[str]:
        """CIDs strictly contained in cid's interval (coordinate-based, ignores parent_id)."""
        me = self.get_candidate(cid)
        if me is None:
            return []
        try:
            ms, me_e = int(me["start_bp"]), int(me["end_bp"])
        except (KeyError, ValueError):
            return []
        out = []
        for r in self._read_cands():
            if r.get("candidate_id") == cid or r.get("chrom") != me.get("chrom"):
                continue
            try:
                s, e = int(r["start_bp"]), int(r["end_bp"])
            except (KeyError, ValueError):
                continue
            if s > ms and e < me_e:
                out.append(r["candidate_id"])
        return out

    def classify_relationship(self, cid1: str, cid2: str) -> Optional[str]:
        """equal | nested_1_in_2 | nested_2_in_1 | partial_overlap | disjoint | None."""
        r1 = self.get_candidate(cid1)
        r2 = self.get_candidate(cid2)
        if r1 is None or r2 is None:
            return None
        if r1.get("chrom") != r2.get("chrom"):
            return "disjoint"
        try:
            s1, e1 = int(r1["start_bp"]), int(r1["end_bp"])
            s2, e2 = int(r2["start_bp"]), int(r2["end_bp"])
        except (KeyError, ValueError):
            return None
        if s1 == s2 and e1 == e2:
            return "equal"
        if e1 < s2 or e2 < s1:
            return "disjoint"
        if s1 >= s2 and e1 <= e2 and (s1 > s2 or e1 < e2):
            return "nested_1_in_2"
        if s2 >= s1 and e2 <= e1 and (s2 > s1 or e2 < e1):
            return "nested_2_in_1"
        return "partial_overlap"

    def add_candidate(
        self,
        cid: str,
        chrom: str,
        start_bp: int,
        end_bp: int,
        scale: str = "inversion_block",
        parent_id: Optional[str] = None,
    ) -> bool:
        """Append a candidate. Skips if cid already present (returns False)."""
        rows = self._read_cands()
        if any(r.get("candidate_id") == cid for r in rows):
            return False
        size_kb = round((int(end_bp) - int(start_bp)) / 1000.0, 1)
        rows.append({
            "candidate_id": cid,
            "chrom": chrom,
            "start_bp": str(int(start_bp)),
            "end_bp": str(int(end_bp)),
            "size_kb": str(size_kb),
            "scale": scale,
            "parent_id": parent_id or "",
        })
        _write_tsv(
            self.cand_file, rows,
            ["candidate_id", "chrom", "start_bp", "end_bp",
             "size_kb", "scale", "parent_id"],
        )
        return True

    def update_candidate(self, cid: str, **updates) -> bool:
        """Partial update of a single candidate row. Returns False if cid not found."""
        rows = self._read_cands()
        if not rows:
            return False
        fieldnames = list(rows[0].keys())
        unknown = [k for k in updates if k not in fieldnames]
        if unknown:
            print(f"[registry] update_candidate: unknown columns skipped: {unknown}")
        idx = next((i for i, r in enumerate(rows) if r.get("candidate_id") == cid), None)
        if idx is None:
            print(f"[registry] update_candidate: unknown CID {cid}")
            return False
        for k, v in updates.items():
            if k in fieldnames:
                rows[idx][k] = "" if v is None else str(v)
        if "start_bp" in updates or "end_bp" in updates:
            try:
                s = int(rows[idx]["start_bp"])
                e = int(rows[idx]["end_bp"])
                rows[idx]["size_kb"] = str(round((e - s) / 1000.0, 1))
            except (KeyError, ValueError):
                pass
        _write_tsv(self.cand_file, rows, fieldnames)
        return True

    def bulk_add_candidates(self, rows: List[Dict[str, Any]]) -> int:
        """Batch insert. rows: list of dicts with required keys
        candidate_id, chrom, start_bp, end_bp; optional: scale, parent_id.
        Duplicates skipped silently. Returns count inserted.
        """
        required = {"candidate_id", "chrom", "start_bp", "end_bp"}
        existing = self._read_cands()
        fieldnames = (
            list(existing[0].keys()) if existing
            else ["candidate_id", "chrom", "start_bp", "end_bp",
                  "size_kb", "scale", "parent_id"]
        )
        existing_ids = {r["candidate_id"] for r in existing if "candidate_id" in r}
        n_ins = 0
        for r in rows:
            missing = required - set(r.keys())
            if missing:
                raise ValueError(
                    f"bulk_add_candidates: row missing required cols {missing}: {r}"
                )
            if r["candidate_id"] in existing_ids:
                continue
            s = int(r["start_bp"])
            e = int(r["end_bp"])
            new_row = {
                "candidate_id": r["candidate_id"],
                "chrom": r["chrom"],
                "start_bp": str(s),
                "end_bp": str(e),
                "size_kb": str(round((e - s) / 1000.0, 1)),
                "scale": r.get("scale", "inversion_block"),
                "parent_id": r.get("parent_id") or "",
            }
            # Fill in any extra columns the existing file has
            for fn in fieldnames:
                new_row.setdefault(fn, "")
            existing.append(new_row)
            existing_ids.add(r["candidate_id"])
            n_ins += 1
        if n_ins > 0:
            _write_tsv(self.cand_file, existing, fieldnames)
        return n_ins


# =============================================================================
# Top-level loader
# =============================================================================
class ResultsAPI:
    """Python binding for the results_registry (chat 16).

    Atomic writer and read-only query. Parity with the R API for the
    common path; advanced queries (overlap semantics, integrity_check)
    stay in R.

    See registries/DATABASE_DESIGN.md for the data model and
    registries/schemas/registry_schemas/result_row.schema.json for the
    manifest row schema.
    """

    MANIFEST_COLS = [
        "row_id", "kind",
        "chrom", "start_bp", "end_bp", "candidate_id",
        "group_1", "group_1_version",
        "group_2", "group_2_version",
        "stat", "K",
        "file",
        "n_rows", "n_cols", "sha256",
        "engine", "engine_version", "source_script",
        "run_id", "config_hash", "upstream_files",
        "timestamp",
    ]

    def __init__(self, results_dir: Path,
                  samples: Optional["SamplesAPI"] = None,
                  intervals: Optional["IntervalsAPI"] = None):
        self.results_dir = Path(results_dir)
        self.samples = samples
        self.intervals = intervals
        for sub in ("pairwise", "candidate", "interval"):
            (self.results_dir / sub).mkdir(parents=True, exist_ok=True)
        self.manifest_file = self.results_dir / "manifest.tsv"

    # ── helpers ──
    @staticmethod
    def _pair_name(g1: str, g2: str) -> str:
        a, b = sorted([str(g1), str(g2)])
        return f"{a}__vs__{b}"

    @staticmethod
    def _gen_uuid() -> str:
        import uuid
        return str(uuid.uuid4())

    def _get_group_version(self, group_id: str) -> str:
        if self.samples is None:
            return ""
        try:
            for g in self.samples.list_groups():
                if g.get("group_id") == group_id:
                    return str(g.get("created", ""))
        except Exception:
            pass
        return ""

    def _check_group_fk(self, group_id: str, context: str = "group") -> None:
        if not group_id:
            raise ValueError(
                f"[results_registry] {context} group_id empty — "
                "refusing to write")
        if self.samples is not None:
            if not self.samples.has_group(group_id):
                raise ValueError(
                    f"[results_registry] {context} group_id '{group_id}' not "
                    "registered in sample_registry")

    def _check_candidate_fk(self, cid: str) -> None:
        if not cid:
            raise ValueError("[results_registry] candidate_id empty")
        if self.intervals is not None:
            c = self.intervals.get_candidate(cid)
            if c is None:
                raise ValueError(
                    f"[results_registry] candidate_id '{cid}' not registered "
                    "in interval_registry")

    @staticmethod
    def _sha256_file(path: Path) -> str:
        import hashlib
        if not path.exists():
            return ""
        h = hashlib.sha256()
        try:
            with open(path, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    h.update(chunk)
            return h.hexdigest()
        except Exception:
            return ""

    def _append_row(self, row: Dict[str, Any]) -> None:
        # Fill missing cols with empty
        full = {c: row.get(c, "") for c in self.MANIFEST_COLS}
        need_header = not self.manifest_file.exists()
        mode = "a" if self.manifest_file.exists() else "w"
        with open(self.manifest_file, mode) as fh:
            if need_header:
                fh.write("\t".join(self.MANIFEST_COLS) + "\n")
            fh.write("\t".join(
                "" if full[c] is None else str(full[c])
                for c in self.MANIFEST_COLS) + "\n")

    # ── writers ──
    def put_pairwise(self, chrom: str, group_1: str, group_2: str,
                      stat: str, rows: List[Dict[str, Any]],
                      fieldnames: List[str],
                      source_script: str = "unknown",
                      engine: Optional[str] = None,
                      engine_version: Optional[str] = None,
                      run_id: Optional[str] = None,
                      config_hash: Optional[str] = None,
                      upstream_files: Optional[List[str]] = None) -> Path:
        self._check_group_fk(group_1, "group_1")
        self._check_group_fk(group_2, "group_2")
        outdir = self.results_dir / "pairwise" / stat / chrom
        outdir.mkdir(parents=True, exist_ok=True)
        pair = self._pair_name(group_1, group_2)
        # Python writes uncompressed .tsv; the consumer (R) can read either.
        out = outdir / f"{pair}.tsv"
        _write_tsv(out, rows, fieldnames)
        rel = out.relative_to(self.results_dir).as_posix()
        self._append_row({
            "row_id": self._gen_uuid(), "kind": "pairwise",
            "chrom": chrom,
            "group_1": group_1, "group_1_version": self._get_group_version(group_1),
            "group_2": group_2, "group_2_version": self._get_group_version(group_2),
            "stat": stat,
            "file": rel,
            "n_rows": len(rows),
            "n_cols": len(fieldnames),
            "sha256": self._sha256_file(out),
            "engine": engine or "",
            "engine_version": engine_version or "",
            "source_script": source_script,
            "run_id": run_id or os.environ.get("SLURM_JOB_ID", ""),
            "config_hash": config_hash or "",
            "upstream_files": ";".join(upstream_files) if upstream_files else "",
            "timestamp": _now_iso(),
        })
        return out

    def put_interval_summary(self, chrom: str, start_bp: int, end_bp: int,
                              group_1: str, stat: str,
                              rows: List[Dict[str, Any]],
                              fieldnames: List[str],
                              K: Optional[int] = None,
                              source_script: str = "unknown",
                              engine: Optional[str] = None,
                              engine_version: Optional[str] = None,
                              run_id: Optional[str] = None,
                              config_hash: Optional[str] = None,
                              upstream_files: Optional[List[str]] = None) -> Path:
        self._check_group_fk(group_1, "group_1")
        if stat not in ("delta12", "entropy", "ena", "ancestry_q_mean"):
            raise ValueError(
                f"[results_registry] stat '{stat}' not valid for interval_summary")
        outdir = self.results_dir / "interval" / chrom
        outdir.mkdir(parents=True, exist_ok=True)
        k_suffix = f".K{int(K):02d}" if K is not None else ""
        out = outdir / f"{int(start_bp)}_{int(end_bp)}.{group_1}.{stat}{k_suffix}.tsv"
        _write_tsv(out, rows, fieldnames)
        rel = out.relative_to(self.results_dir).as_posix()
        self._append_row({
            "row_id": self._gen_uuid(), "kind": "interval_summary",
            "chrom": chrom, "start_bp": int(start_bp), "end_bp": int(end_bp),
            "group_1": group_1, "group_1_version": self._get_group_version(group_1),
            "stat": stat, "K": "" if K is None else int(K),
            "file": rel,
            "n_rows": len(rows), "n_cols": len(fieldnames),
            "sha256": self._sha256_file(out),
            "engine": engine or "", "engine_version": engine_version or "",
            "source_script": source_script,
            "run_id": run_id or os.environ.get("SLURM_JOB_ID", ""),
            "config_hash": config_hash or "",
            "upstream_files": ";".join(upstream_files) if upstream_files else "",
            "timestamp": _now_iso(),
        })
        return out

    # ── readers ──
    def list_manifest(self, kind: Optional[str] = None) -> List[Dict[str, str]]:
        if not self.manifest_file.exists():
            return []
        rows = _read_tsv(self.manifest_file)
        if kind is not None:
            rows = [r for r in rows if r.get("kind") == kind]
        return rows

    def ask(self, where: Optional[Dict[str, Any]] = None,
             who: Optional[Any] = None,
             what: Optional[str] = None,
             kind: Optional[str] = None,
             K: Optional[int] = None) -> List[Dict[str, str]]:
        """Python-side read-only query. No overlap semantics — use R's
        reg$results$ask() for interval overlap queries."""
        rows = self.list_manifest(kind=kind)
        if what is not None:
            rows = [r for r in rows if r.get("stat") == what]
        if K is not None:
            rows = [r for r in rows if str(r.get("K", "")) == str(K)]
        if who is not None:
            wset = {who} if isinstance(who, str) else set(who)
            rows = [r for r in rows
                    if r.get("group_1") in wset or r.get("group_2") in wset]
        if where:
            if where.get("candidate_id"):
                rows = [r for r in rows
                        if r.get("candidate_id") == where["candidate_id"]]
            if where.get("chrom"):
                rows = [r for r in rows if r.get("chrom") == where["chrom"]]
        return rows


@dataclass
class Registry:
    root: Path
    samples: SamplesAPI
    intervals: IntervalsAPI
    evidence: EvidenceAPI
    results: ResultsAPI


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
                "evidence_registry/global",
                "results_registry"):
        (data_dir / sub).mkdir(parents=True, exist_ok=True)

    # Chat-16: auto-migrate chat-15 stats_cache/ → results_registry/
    legacy_stats = data_dir / "stats_cache"
    new_results = data_dir / "results_registry"
    if legacy_stats.exists() and not any(new_results.iterdir()):
        try:
            # Move contents of legacy into new_results
            for p in legacy_stats.iterdir():
                p.rename(new_results / p.name)
            legacy_stats.rmdir()
            print("[registry_loader.py] migrated stats_cache/ → results_registry/")
        except Exception as e:
            print(f"[registry_loader.py] could not migrate stats_cache/: {e}")

    print(f"[registry_loader.py] root = {root}")

    samples   = SamplesAPI(data_dir / "sample_registry")
    intervals = IntervalsAPI(data_dir / "interval_registry")
    evidence  = EvidenceAPI(data_dir / "evidence_registry", schema_dir)
    # Resolve results dir: RESULTS_REGISTRY_DIR env > STATS_CACHE_DIR env > default
    rdir_env = os.environ.get("RESULTS_REGISTRY_DIR") or \
                os.environ.get("STATS_CACHE_DIR")
    results = ResultsAPI(
        Path(rdir_env) if rdir_env else (data_dir / "results_registry"),
        samples=samples, intervals=intervals)

    return Registry(
        root=root,
        samples=samples,
        intervals=intervals,
        evidence=evidence,
        results=results,
    )


if __name__ == "__main__":
    # Smoke test
    reg = load_registry()
    print("samples groups:", len(reg.samples.list_groups()))
    print("candidates:", len(reg.evidence.list_candidates()))
    print("results rows:", len(reg.results.list_manifest()))
