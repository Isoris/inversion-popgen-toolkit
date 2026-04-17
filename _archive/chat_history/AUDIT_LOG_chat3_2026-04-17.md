# AUDIT LOG — chat 3 of the 4a audit series (2026-04-17)

This chat was a *continuation* of the phase 2c / 4a audit. Earlier
chats in this session had already identified the 14 bugs and applied
fixes 1-14. This chat's contribution was:

- **Full parse-check of every modified R file** with the R interpreter
  (not just a syntactic grep). All 5 modified files parse cleanly.
  Full-repo parse sweep found 6 pre-existing parse errors, all in
  `phase_2_discovery/2c_precomp/patches/` (documented legacy reference
  material, never sourced at runtime — not my fault, not my fix).
- **Downstream column contract verification** — traced every column
  name that 4b/4c/4e/test_registry_sanity.py reads from
  `candidate_scores.tsv.gz` back to a column that C01d writes. All
  names match including v9.3.2 additions (`d11_boundary_concordance`,
  `d12_snake_concordance`, `snake_overlap`, `cheat25_status`).
- **One late-session find** that was NOT fixed this chat: C01f
  (`phase_4_postprocessing/4c_group_validation/STEP_C01f_hypothesis_tests.R`)
  has the same triangle_sample_composition read pattern that was fixed
  in C01e last chat. In C01f it's partially papered over by a
  registry-first fallback (line 354-368: `build_comp_with_fallback`),
  so it's not a silent wrong-output bug in the usual case — but the
  second fallback is the obsolete triangle file instead of on-the-fly
  k-means. Port the `compute_bands_for_candidate()` helper from C01e
  to C01f in the next session.

**Non-technical discussion this chat:**

- **Meta-concern raised by previous chats that I agreed with:** we are
  in an optimization trap where code polishing is legible progress but
  the manuscript is what actually needs finishing. Of the 14 fixes,
  ~4 were material (BUG 27, 5, 3, 9) and ~10 were documentation /
  stale-flag hygiene. That 4:10 ratio is the tell. Compute pipeline
  is running in parallel (Clair3 LG13+ tonight, LG17-LG27 over the
  next three days, C00 sv_prior rebuild queued for tomorrow). Use
  compute-wait time for the manuscript, not for watching SLURM.

- **On fossil detection** (Cheat 17 archived this session): Quentin
  correctly noted that in a 226-sample F₁ hybrid hatchery on a
  reasonable timescale, the pop-gen signature of a fossil inversion
  (elevated inv_likeness + no trimodality + no SV) is not
  distinguishable from active rare inversions whose HET samples got
  binned away by k-means, assembly artifacts with mild elevated sim,
  or segmental duplications. He doesn't have the sample count or the
  phylogenetic depth or an ancestral outgroup for unfolded SFS to make
  that call. Archiving the feature is the right call — it would be
  an untestable hypothesis in this system.

- **Tomorrow's plan agreed:**
  - 10 min: smoke test the 14 fixes on LG12 (the Tier-1 candidate at
    52.4-54.6 Mb) against existing precomp cache. Gate.
  - 30 min: submit C00 sv_prior rebuild on all 28 chr.
  - 0 min: Clair3 LG17-LG20 already scheduled.
  - Rest of the day: manuscript.
  - Skip the C01a full rerun unless a diff between old
    `sv_flashlight_<chr>.rds` and new `sv_prior_<chr>.rds` shows
    meaningfully different `inv_calls` tables. In practice, same VCFs
    → same output → skip the rerun → save a day.

**Tarball shipped at end of chat 2 of this series:**
`inversion-popgen-toolkit_phase2_phase4a_fixes_2026-04-17.tar`
(9.5 MB, 957 files)

**What chat 3 of this series added:** the HANDOFF_PROMPT_next_chat
file at repo root describing scope for chat 4 (phase 2e audit + phase
3 audit + finishing the C01f triangle fallback).

**Explicit non-scope for chat 4:**
- Don't touch phase 4b/4c/4d/4e internals beyond what's needed
- Don't touch C01a internals (science audit, separate project)
- Don't touch the flashlight → sv_prior rename of the C01d output
  column name (`cheat5_family_fst_ratio` stays legacy to avoid
  cascading changes to 4e registry keys)
- Don't slide into more code polishing if the code is not actively
  broken — push back if Quentin asks for more audits when the
  manuscript is the real bottleneck.
