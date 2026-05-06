# Specs deferred to later

Specs in this folder are NOT prioritized for current implementation. Reasons
are recorded below.

## SPEC_xpehh_track.md

**Status**: Deferred (not blocking — don't pull into next session unless
explicitly requested).

**Why deferred**:
- Requires phased VCF (BEAGLE rerun) on the 226-sample cohort
- Requires a reference cohort decision (HOM_REF vs HOM_INV haplotypes
  within the 226 — specced day-1 path — vs. an external reference)
- The other four specs in `specs_todo/` are higher value because they
  have no upstream blockers

**To re-prioritize**: move back to `specs_todo/`.
