# SCRIPT_CONTRACT.md — registry contract declared in every script header

Companion to `STANDARDIZE_SCRIPT.md`. The standardization prompt defines the
shape of the script-level header block (PURPOSE / INPUTS / OUTPUTS / DEPENDENCIES
/ USED BY / STATUS / USAGE / NOTES). This doc adds **one required optional
section** to that block for scripts that participate in the registry: the
REGISTRY_CONTRACT.

## Why

The bucket B audit (2026-04-24) found that half the "PRODUCES_BUT_NOT_WIRED"
schemas had the wrong producer named in their `source_script` field. The
actual data the schema wanted lived in a different script, or wasn't computed
at all. There was no cheap way to detect this from either side:

- Schemas don't know whether their claimed source actually emits their fields.
- Scripts don't know whether the registry expects anything from them.
- `docs/FIELD_MAPPING_PRODUCES.md` is a manual catalog that has to be hand-
  maintained against moving code.

A per-script REGISTRY_CONTRACT header block closes the loop. Each registry-
participating script says, in one place, which blocks it writes, which flat
keys that produces, and where the schemas live. A scanner
(`tools/scan_script_contracts.py`) cross-references contracts against schemas
and flags drift automatically.

## Where the block goes

Directly inside the script's header comment block, between `OUTPUTS` and
`DEPENDENCIES`. Only present for scripts that write to the registry.

## Required fields

```
REGISTRY_CONTRACT
  BLOCKS_WRITTEN:
    - <block_type>: <schema_path relative to repo root>
      keys: <comma-separated list of the flat keys this writer actually
             produces — must match, or be a subset of, the schema's
             keys_extracted>
      status: WIRED | PRODUCES_BUT_NOT_WIRED | BLOCKED_ON_<reason>
      [note: <free text, optional>]
  KEYS_IN:
    - <key>: <which script writes it, if known; or "— (from iv_dt)" etc.>
    [or "none" when the script reads no registry keys]
```

### BLOCKS_WRITTEN semantics

- List **every** `block_type` the script writes via `reg$evidence$write_block`
  or helper calls. A script that writes two different block types (like
  `register_C01d_keys` writing both `existence_layer_a` and `morphology`)
  lists both.
- `keys:` is the subset of the schema's `keys_extracted` that the writer
  actually populates with non-NA values in typical operation. NA-only fields
  are listed separately under a `keys_na:` sub-field when relevant.
- `status`:
  - `WIRED` — writer calls the registry and populates the listed keys.
  - `PRODUCES_BUT_NOT_WIRED` — writer computes the data but doesn't call
    the registry yet. Block contents should land via helper in a future pass.
  - `BLOCKED_ON_<reason>` — short reason tag. Current reasons seen in this
    codebase:
    - `BLOCKED_ON_NO_CANDIDATE_JOIN` — script runs before candidates exist.
    - `BLOCKED_ON_MISSING_PRODUCER_FIELDS` — schema asks for fields the
      upstream producer doesn't compute.
    - `BLOCKED_ON_MISSING_UPSTREAM_COMPUTATION` — a required computation
      doesn't live anywhere in the codebase yet.

### KEYS_IN semantics

Rarely useful today but included for symmetry. Scripts that read keys from
the registry (C01f reading C01i's seal verdicts, C01d reading hyp_verd from
C01f) would list them here. In practice most scripts read directly from
files or in-memory data.tables, not from the registry, so this is usually
`none`.

## Example — ship-ready script

```
REGISTRY_CONTRACT
  BLOCKS_WRITTEN:
    - synteny_v6: registries/schemas/structured_block_schemas/synteny_v6.schema.json
      keys: q5_bs_event_overlap, q5_bs_event_type,
            q5_tree_polarization_direction, q5_tree_polarization_confidence,
            q5_dollo_vs_tree_concordance, q5_conservation_class,
            q3_left_flank_coherence_score, q3_right_flank_coherence_score,
            q3_left_flank_n_families, q3_right_flank_n_families,
            q3_flank_coherence_class
      status: WIRED
      note: Flank-coherence keys folded in from superseded flank_coherence
            schema 2026-04-24; no writer change needed.
  KEYS_IN: none
```

## Example — blocked script

```
REGISTRY_CONTRACT
  BLOCKS_WRITTEN:
    - existence_layer_b: registries/schemas/structured_block_schemas/existence_layer_b.schema.json
      keys: (none — not currently written from this script)
      status: BLOCKED_ON_NO_CANDIDATE_JOIN
      note: C00 produces per-chromosome SV priors BEFORE candidates exist
            (see L1243-1245 of this file). Per-candidate existence_layer_b
            writes need to move to a consumer that runs after C01d creates
            candidate_ids — candidate consumers are C01d itself or a
            dedicated STEP_C00b_attribute_sv_per_candidate.R. Schema
            asks for 6 keys; iv_dt currently propagates only n_sv_hits
            (DRIFT → n_sv_calls). The remaining 5 need a fresh join of
            sv_prior_<chr>.rds × candidate regions.
  KEYS_IN: none
```

## Scanner

`tools/scan_script_contracts.py` walks a directory tree and:

1. For every `.R` / `.py` / `.sh` file, looks for a `REGISTRY_CONTRACT` block
   in the first 200 lines.
2. Parses out each `BLOCKS_WRITTEN` entry.
3. For each `(block_type, schema_path)` pair:
   - Verifies the schema file exists.
   - Parses `keys_extracted` from the schema.
   - Checks the contract's `keys:` list is a subset of the schema's
     extracted keys (or empty, for blocked/not-wired statuses).
   - Flags any orphan schema-declared keys not claimed by any contract.
4. Reports summary by status (WIRED / PRODUCES_BUT_NOT_WIRED / BLOCKED_ON_*).
5. Non-zero exit if any WIRED contract names a key that isn't in its
   schema's `keys_extracted` (drift from writer to schema).

## Adoption

Not all scripts need this at once. Priority order:

1. Every script named as `source_script` in a schema.
2. Every helper in `utils/registry_key_helpers.R` (the helper file itself gets
   one contract block; the listed functions don't need individual blocks).
3. Phase 4 scripts that call `reg$evidence$write_block` or helpers directly.

Scripts with no registry interaction do not need this section. The
standardized header (from `STANDARDIZE_SCRIPT.md`) still applies.

## Reviewing changes

When editing a schema's `keys_extracted`, the scanner will flag any WIRED
contract that named the removed/renamed key. When editing a writer, any
contract that no longer matches the code triggers the scanner. Either way,
the contract forces the edit to be deliberate rather than silent.
