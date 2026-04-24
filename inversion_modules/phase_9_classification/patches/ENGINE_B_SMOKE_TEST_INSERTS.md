# =============================================================================
# PATCH: insert Engine B smoke test at top of every script that calls
#        the unified ancestry dispatcher
# =============================================================================
# Apply to:
#   - STEP_C01a_snake1_precompute_*.R       BEFORE the main loop
#   - ancestry_bridge.R                      BEFORE --prepare does any work
#   - STEP_C01i_b_multi_recomb.R             IF phase-switch + Engine B used
#   - STEP_C01i_c_nested_composition.py      separate Python check
#   - cheat6_ancestry_jackknife_v934.R       BEFORE any Fst calls
#
# The check adds ~5s to startup. Skipping it saves 5s; running it on a
# broken Engine B saves ~6 hours of computed-then-discarded NAs.
#
# Override with SKIP_SMOKE_TEST=1 env var for dev work when you know
# the test will fail but want to run anyway.
# =============================================================================


# ─── FOR ANY R SCRIPT THAT SOURCES load_bridge.R ─────────────────────────────
#
# Insert immediately AFTER sourcing load_bridge.R, BEFORE any get_region_stats
# or get_Q call. Easiest marker: look for the line that sources load_bridge
# and add the smoke test right after.
#
# BEFORE:
#   source(bridge_path)  # or similar
#   # ... main work begins here ...
#
# AFTER:
#   source(bridge_path)
#
#   # ── Engine B self-test (fail fast, skip with SKIP_SMOKE_TEST=1) ──
#   if (!nzchar(Sys.getenv("SKIP_SMOKE_TEST", ""))) {
#     for (sp in c("utils/engine_b_smoke_test.R",
#                  "../utils/engine_b_smoke_test.R",
#                  file.path(Sys.getenv("BASE", ""),
#                            "inversion-popgen-toolkit/utils/engine_b_smoke_test.R"))) {
#       if (file.exists(sp)) { source(sp); break }
#     }
#     if (exists("run_engine_b_smoke_test", mode = "function")) {
#       ok <- run_engine_b_smoke_test()
#       if (!ok) {
#         stop("[", script_name, "] Engine B smoke test failed — aborting before ",
#              "long compute. Override with SKIP_SMOKE_TEST=1 to bypass.")
#       }
#     }
#   }
#   # ... main work begins here ...


# ─── FOR ancestry_bridge.R --prepare ─────────────────────────────────────────
#
# Insert at the top of the --prepare branch. Same code block as above. If
# the smoke test fails, ancestry_bridge should abort before running
# instant_q for 28 chromosomes.


# ─── FOR Python scripts (STEP_C01i_c_nested_composition.py) ──────────────────
#
# The Python wrapper doesn't call Engine B directly — it reads the cached
# .local_Q_samples.tsv.gz produced by ancestry_bridge. So the smoke test
# for Python is different: check that the cache exists and is non-empty.
#
# Insert early in main() BEFORE parent processing:
#
# def check_q_cache(cache_dir: str, sample_chr: str = "LG01") -> bool:
#     """5-second smoke test: can we read one chromosome's Q cache?"""
#     import gzip
#     if not cache_dir or not os.path.isdir(cache_dir):
#         print(f"[smoke_test] ✗ Q cache dir missing: {cache_dir}")
#         print("[smoke_test] Run ancestry_bridge.R --prepare first.")
#         return False
#     for ext in [".local_Q_samples.tsv.gz", ".local_Q_samples.tsv"]:
#         path = os.path.join(cache_dir, sample_chr + ext)
#         if os.path.isfile(path):
#             opener = gzip.open if path.endswith(".gz") else open
#             try:
#                 with opener(path, "rt") as f:
#                     header = f.readline()
#                     if "assigned_pop" not in header:
#                         print(f"[smoke_test] ✗ {path} missing assigned_pop column")
#                         return False
#                     first = f.readline()
#                     if not first:
#                         print(f"[smoke_test] ✗ {path} empty")
#                         return False
#                 print(f"[smoke_test] ✓ Q cache OK ({path})")
#                 return True
#             except Exception as e:
#                 print(f"[smoke_test] ✗ {path} read error: {e}")
#                 return False
#     print(f"[smoke_test] ✗ no .local_Q_samples file in {cache_dir}")
#     return False
#
#
# And in main():
#   if not os.environ.get("SKIP_SMOKE_TEST"):
#       if q_cache_available and not check_q_cache(q_cache, ...):
#           print("[nested_comp] Q cache smoke test failed — writing stubs")
#           q_cache_available = False   # fall through to stub-writing path


# ─── FOR SHELL ORCHESTRATOR run_phase4b.sh ──────────────────────────────────
#
# Before submitting long SLURM jobs, run a small interactive check:
#
# SELF_CHECK_OK=1
# if [[ -z "${SKIP_SMOKE_TEST:-}" ]]; then
#   echo "[run_phase4b] running Engine B smoke test..."
#   RUN_SMOKE_TEST=1 Rscript "${MODULE_DIR}/utils/engine_b_smoke_test.R" \
#     || { echo "[run_phase4b] smoke test failed; aborting"; SELF_CHECK_OK=0; }
# fi
# if [[ "${SELF_CHECK_OK}" -eq 0 ]]; then
#   echo ""
#   echo "[run_phase4b] Engine B is not ready. Common fixes:"
#   echo "  cd \${BASE}/unified_ancestry/src && make"
#   echo "  cd \${BASE}/unified_ancestry/engines/fst_dxy && make"
#   echo "  Rscript \${BASE}/utils/ancestry_bridge.R --prepare --K 8"
#   echo ""
#   echo "Or bypass with: SKIP_SMOKE_TEST=1 bash run_phase4b.sh"
#   exit 1
# fi


# ─── SUMMARY ────────────────────────────────────────────────────────────────
#
# Five seconds at startup, catches ~90% of "silent Engine B failure"
# scenarios. Override with SKIP_SMOKE_TEST=1 when you explicitly want to
# run without Engine B (e.g., on a dev node where the binaries aren't
# compiled yet and you're testing the R logic).
#
# This is defensive; it does NOT test whether the numbers are correct.
# It tests whether the plumbing works. For correctness you still need
# to eyeball the first few candidates or run a known-good reference.
