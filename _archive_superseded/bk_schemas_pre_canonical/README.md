# bk_schemas_pre_canonical/

Archive of files superseded by chat-15 Finding BK (phase-4b schema
canonicalisation).

## internal_dynamics.schema.json.pre_BK

The pre-chat-12 flavor of the canonical internal_dynamics schema
that formerly lived at
`registries/schemas/structured_block_schemas/internal_dynamics.schema.json`.

Described the recombinant-carrying block (per-sample recombinants
array, cheat24 posterior, event_class enum of gene_conversion /
double_crossover / suspicious) that STEP_C01i_b wrote before the
chat-12 rewrite. Every `keys_extracted` entry referenced a field that
C01i_decompose does not write, so all 12 extractions silently skipped.

Superseded by the chat-15 rewrite which matches what
STEP_C01i_decompose.R actually writes today (silhouette / BIC /
flashlight / seed counts / PCA class counts).

## phase4b_drafts/

Draft schemas that formerly sat at
`inversion_modules/phase_4_postprocessing/schemas/` — design-time files
never read by `registries/api/R/registry_loader.R`. They lived parallel
to the canonical directory because the phase-4b rewrite was in flight
when they were authored.

- `internal_dynamics.schema.json` — draft was also decompose-flavor
  (closer to reality than the canonical pre-BK file) but its
  `keys_extracted` set was narrower and omitted `decomp_status` /
  `decomp_reason` which are the entry-point fields for downstream
  consumers to branch on emit-variant.
- `recombinant_map.schema.json` — draft `keys_extracted` named the
  pre-chat-12 fields (`n_gene_conversion`, `n_double_crossover`,
  `mean_posterior`, `mean_dco_prior`, `max_dco_prior`). The chat-12
  rewrite of C01i_b emits `n_recombinant_gc`, `n_recombinant_dco`,
  `n_disputed`, `n_ghsl_only`, `cheat24_version`, `combination_rule`,
  `gate_params.*` instead.
- `internal_ancestry_composition.schema.json` — the one draft whose
  `keys_extracted` more-or-less matched current C01i_c output. The
  canonical version keeps the same field definitions and adds
  `K_used` / `n_samples_analyzed` / `pct_diffuse_mixed` to the
  extracted set (12 keys vs draft's 9).

Not archived: `frequency.v2.schema.json` — out of BK scope; no canonical
`frequency.schema.json` exists yet either, and the chat-15 BK fix was
strictly about the three named blocks that C01i_b / C01i_c / C01i_decompose
actually write today.

## Not to be confused with

- `_archive_superseded/ghsl_v5_consumers_pre_v6_rewire/` — chat-14
  pre-patch copies of the four consumer scripts that were rewired from
  v5 to v6 GHSL.
