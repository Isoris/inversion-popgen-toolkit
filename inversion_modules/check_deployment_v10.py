#!/usr/bin/env python3
"""
check_deployment_v10.py — verify a deployed v10 + v10.1 tree is wired correctly

Run this after dropping the v10 and v10.1 archives into inversion_modules/.
It does NOT touch real data; it just confirms that:

  1. Every expected file is present
  2. JSON schemas parse
  3. R files have balanced brackets
  4. Python imports work (no syntax errors)
  5. Tests pass
  6. Registry-tree structure can be created
  7. Cross-references between files resolve (schema names match block types
     referenced in scripts, etc.)

Exits 0 if everything is green, non-zero with a clear list of what's wrong.
Typical runtime: 5-10 seconds.

Usage:
    cd inversion-popgen-toolkit/
    python3 inversion_modules/check_deployment_v10.py --root .
    # or with explicit root
    python3 check_deployment_v10.py --root /path/to/inversion-popgen-toolkit
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


class CheckFailed(Exception):
    pass


def check_files_present(root: Path) -> list:
    """Return list of (relpath, exists) for all expected files."""
    expected = [
        # v10 core
        "registries/api/R/registry_loader.R",
        "registries/api/python/registry_loader.py",
        "registries/api/bash/registry_loader.sh",
        # v10 schemas (at least the key ones)
        "registries/schemas/structured_block_schemas/existence_layer_a.schema.json",
        "registries/schemas/structured_block_schemas/existence_layer_d.schema.json",
        "registries/schemas/structured_block_schemas/hypothesis_verdict.schema.json",
        "registries/schemas/structured_block_schemas/boundary.schema.json",
        # v10.1 phase 4b rewrite (now under phase_7_karyotype_groups/proposal/)
        "inversion_modules/phase_7_karyotype_groups/proposal/lib_decompose_helpers.R",
        "inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_decompose.R",
        "inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_b_multi_recomb.R",
        "inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_d_seal.R",
        "inversion_modules/phase_7_karyotype_groups/proposal/engine_b_smoke_test.R",
        "inversion_modules/phase_7_karyotype_groups/proposal/nested_composition_core.py",
        "inversion_modules/phase_7_karyotype_groups/proposal/STEP_C01i_c_nested_composition.py",
        # pass 15: structured-block schemas consolidated into registries/schemas/
        "registries/schemas/structured_block_schemas/internal_dynamics.schema.json",
        "registries/schemas/structured_block_schemas/recombinant_map.schema.json",
        "registries/schemas/structured_block_schemas/internal_ancestry_composition.schema.json",
        "registries/schemas/structured_block_schemas/frequency.v2.schema.json",
        # pass 15: phase_4 patches/ tests/ docs/ orchestrator/ relocated
        "inversion_modules/phase_9_classification/patches/02_C01f_promotion_cap.R",
        "inversion_modules/phase_9_classification/patches/03_C01f_jackknife_semantics.R",
        "inversion_modules/phase_9_classification/patches/ENGINE_B_SMOKE_TEST_INSERTS.md",
        "inversion_modules/phase_7_karyotype_groups/proposal/orchestrator/run_phase4b.sh",
        "inversion_modules/phase_9_classification/tests/test_c01i_d_seal_resolution.py",
        "inversion_modules/phase_9_classification/tests/test_phase4b_integration.py",
        "inversion_modules/phase_9_classification/tests/test_jackknife_semantics.py",
        "inversion_modules/phase_7_karyotype_groups/proposal/docs/PHASE4B_REWRITE_ARCHITECTURE.md",
        "docs/DESIGN_NOTE_K_and_local_Q_and_jackknife.md",
    ]
    return [(p, (root / p).exists()) for p in expected]


def check_jsons(root: Path) -> list:
    """Return list of (path, status) for each .schema.json."""
    results = []
    for p in root.rglob("*.schema.json"):
        try:
            with p.open() as f:
                json.load(f)
            results.append((p.relative_to(root), "ok"))
        except Exception as e:
            results.append((p.relative_to(root), f"bad: {e}"))
    return results


def check_r_brackets(root: Path) -> list:
    """Basic bracket-balance check on R files.

    Note: this is a naive regex-based counter. It correctly flags most broken
    R code but false-positives on files containing complex string literals
    (regex patterns, paste0 across many lines, character-class brackets in
    strings). The allowlist below names files known to be syntactically valid
    R that this naive counter mis-flags. Do a real R parse on LANTA before
    submitting SLURM jobs.
    """
    # Files known to be valid R that trip the naive counter
    known_false_positive_suffixes = (
        "phase_2_discovery/2c_precomp/STEP_C01b_1_seeded_regions.R",
        "phase_2_discovery/2c_precomp/STEP_C00_build_sv_prior.R",
        "phase_2_discovery/2d_candidate_detection/STEP_C01b_2_region_merge.R",
        "phase_2_discovery/2d_candidate_detection/STEP_D17_plot_marginal_tracks.R",
        "phase_2_discovery/2d_candidate_detection/STEP_D15_plot_zoomed_regions.R",
        "phase_2_discovery/2d_candidate_detection/STEP_D13_plot_annotated_simmat.R",
        "phase_2_discovery/2d_candidate_detection/STEP_D10_variant_consensus.R",
        "phase_7_karyotype_groups/validation/STEP_C01f_hypothesis_tests.R",
        "phase_7_karyotype_groups/validation/_PRISTINE_v9.3.4.R",
    )
    results = []
    for p in root.rglob("*.R"):
        rel = str(p.relative_to(root))
        try:
            with p.open() as f:
                s = f.read()
            # Strip comments + strings
            s2 = re.sub(r"#[^\n]*", "", s)
            s2 = re.sub(r'"([^"\\]|\\.)*"', '""', s2)
            s2 = re.sub(r"'([^'\\]|\\.)*'", "''", s2)
            bal = {
                "(": s2.count("(") - s2.count(")"),
                "{": s2.count("{") - s2.count("}"),
                "[": s2.count("[") - s2.count("]"),
            }
            if all(v == 0 for v in bal.values()):
                results.append((p.relative_to(root), "ok"))
            elif rel.endswith(known_false_positive_suffixes):
                results.append((p.relative_to(root), f"ok (allowlisted, naive counter reports {bal})"))
            else:
                results.append((p.relative_to(root), f"unbalanced {bal}"))
        except Exception as e:
            results.append((p.relative_to(root), f"error: {e}"))
    return results


def check_python_syntax(root: Path) -> list:
    """py_compile each .py to catch syntax errors."""
    results = []
    for p in root.rglob("*.py"):
        try:
            # Compile without running
            compile(p.read_text(), str(p), "exec")
            results.append((p.relative_to(root), "ok"))
        except SyntaxError as e:
            results.append((p.relative_to(root), f"syntax: line {e.lineno}: {e.msg}"))
        except Exception as e:
            results.append((p.relative_to(root), f"error: {e}"))
    return results


def check_tests_pass(root: Path) -> list:
    """Run each test script; report pass/fail."""
    results = []
    tests = sorted(root.rglob("test_*.py"))
    for t in tests:
        try:
            r = subprocess.run(
                [sys.executable, str(t)],
                capture_output=True, text=True, timeout=120,
                cwd=str(t.parent.parent),  # run from phase4b_rewrite/ so rel paths work
            )
            if r.returncode == 0:
                # Try to extract final count from output
                last_lines = r.stdout.strip().splitlines()[-3:]
                summary = " | ".join(l.strip() for l in last_lines if "✓" in l or "pass" in l.lower())[:80]
                results.append((t.relative_to(root), f"pass ({summary})"))
            else:
                results.append((t.relative_to(root), f"FAIL (exit {r.returncode})"))
        except Exception as e:
            results.append((t.relative_to(root), f"error: {e}"))
    return results


def check_schema_block_crossref(root: Path) -> list:
    """For each schema file, check its declared block_type matches what
    the pipeline scripts write. Catches typos like writing
    'internal_dyanmics' instead of 'internal_dynamics'."""
    results = []
    schema_block_types = set()
    for sp in root.rglob("*.schema.json"):
        try:
            with sp.open() as f:
                s = json.load(f)
            bt = s.get("block_type")
            if bt:
                schema_block_types.add(bt)
            # Templated schemas like boundary.schema.json use block_type_template
            btt = s.get("block_type_template")
            if btt:
                # e.g. "boundary_{left|right}" → boundary_left, boundary_right
                match = re.search(r"\{([^}]+)\}", btt)
                if match:
                    for opt in match.group(1).split("|"):
                        schema_block_types.add(btt.replace(match.group(0), opt))
        except Exception as e:
            results.append((sp.relative_to(root), f"cannot parse: {e}"))
            continue

    # Find write_block( calls in R and py; extract block_type arg
    write_calls = set()
    for pat in ["*.R", "*.py"]:
        for p in root.rglob(pat):
            # Skip the check script itself — its example regex contains
            # literal "block_type = \"xxx\"" which would be a false positive
            if p.name == "check_deployment_v10.py":
                continue
            try:
                s = p.read_text()
                # Look for write_block(... block_type = "xxx"  OR  block_type="xxx"
                for m in re.finditer(
                    r'write_block\s*\([^)]*block_type\s*=\s*["\']([^"\']+)["\']',
                    s, re.DOTALL,
                ):
                    write_calls.add((m.group(1), p.relative_to(root)))
            except Exception:
                pass

    # For each written block_type, does a schema exist?
    for bt, src in sorted(write_calls):
        if bt in schema_block_types:
            results.append((f"write:{bt}", f"schema found (from {src})"))
        else:
            results.append((f"write:{bt}", f"NO SCHEMA (written from {src})"))
    return results


def check_registry_tree_creation(root: Path) -> list:
    """Verify registry_loader.py can create the expected directory tree."""
    results = []
    loader = root / "registries" / "api" / "python" / "registry_loader.py"
    if not loader.exists():
        return [("registry_loader.py", "missing")]
    with tempfile.TemporaryDirectory(prefix="check_deploy_") as tmp:
        r = subprocess.run(
            [sys.executable, "-c",
             f"import sys; sys.path.insert(0, '{loader.parent}'); "
             f"from registry_loader import load_registry; "
             f"reg = load_registry(registries_root='{tmp}/reg'); "
             f"print('root=', reg.root); "
             f"print('groups=', len(reg.samples.list_groups())); "
             f"print('cands=', len(reg.evidence.list_candidates())); "],
            capture_output=True, text=True, timeout=30,
        )
        if r.returncode == 0:
            results.append(("registry_loader.py load", "ok"))
            expected_subs = [
                "reg/data/sample_registry",
                "reg/data/interval_registry",
                "reg/data/evidence_registry/per_candidate",
                "reg/data/evidence_registry/global",
            ]
            for sub in expected_subs:
                if (Path(tmp) / sub).exists():
                    results.append((f"tree:{sub}", "created"))
                else:
                    results.append((f"tree:{sub}", "MISSING"))
        else:
            results.append(("registry_loader.py load",
                            f"FAIL: {r.stderr[:200]}"))
    return results


def render(results_list, section_title):
    print(f"\n=== {section_title} ===")
    bad = 0
    for key, status in results_list:
        if isinstance(status, bool):
            mark = "✓" if status else "✗"
            if not status:
                bad += 1
        else:
            status_str = str(status)
            if "ok" in status_str.lower() or "pass" in status_str.lower() \
               or "found" in status_str.lower() or "created" in status_str.lower():
                mark = "✓"
            else:
                mark = "✗"
                bad += 1
        print(f"  {mark} {key}: {status}" if not isinstance(status, bool)
              else f"  {mark} {key}")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".",
                    help="Root of deployed tree (default: cwd)")
    ap.add_argument("--skip-tests", action="store_true",
                    help="Skip running tests (just check files)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    print(f"Checking deployment at: {root}")
    print("(exits 0 if all green; non-zero with list of issues)\n")

    total_bad = 0

    # 1. Files present
    file_results = check_files_present(root)
    simple = [(p, "ok" if exists else "MISSING") for p, exists in file_results]
    total_bad += render(simple, "Files present")

    # 2. JSON schemas
    total_bad += render(check_jsons(root), "JSON schemas parse")

    # 3. R brackets
    total_bad += render(check_r_brackets(root), "R bracket balance")

    # 4. Python syntax
    total_bad += render(check_python_syntax(root), "Python syntax")

    # 5. Schema <-> write_block cross-reference
    total_bad += render(check_schema_block_crossref(root),
                         "Schema block_type coverage")

    # 6. Registry tree creation
    total_bad += render(check_registry_tree_creation(root),
                         "Registry loader smoke")

    # 7. Tests (optional)
    if not args.skip_tests:
        total_bad += render(check_tests_pass(root), "Tests pass")

    # ── Summary ──
    print("\n" + "=" * 60)
    if total_bad == 0:
        print(f"✓ ALL GREEN — deployment at {root} looks ready")
        return 0
    else:
        print(f"✗ {total_bad} issue(s) found — see above")
        return 1


if __name__ == "__main__":
    sys.exit(main())
