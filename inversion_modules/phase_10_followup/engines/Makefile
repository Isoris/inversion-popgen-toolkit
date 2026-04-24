CC      = gcc
CFLAGS  = -O3 -march=native -fopenmp -Wall
LDFLAGS = -lz -lm -fopenmp

# -----------------------------------------------------------------------------
# MODULE_5B_inversion_followup/engines/Makefile
#
# Builds the two C engines that feed the 5B inversion-followup figure scripts.
# These engines were moved out of unified_ancestry/engines/fst_dxy/ because
# they exist specifically to produce data for 5B R plots:
#
#   export_q_residual_dosage  →  plot_ld_panels.R (Panel 3: Q-corrected LD)
#   rare_sfs_pairwise         →  07_plot_rare_sfs_heatmap.R
#
# The candidate-agnostic region_popstats engine remains in
# unified_ancestry/engines/fst_dxy/.
# -----------------------------------------------------------------------------

all: export_q_residual_dosage rare_sfs_pairwise

export_q_residual_dosage: export_q_residual_dosage.c
	$(CC) $(CFLAGS) -o $@ $< -lz -lm

rare_sfs_pairwise: rare_sfs_pairwise.c
	$(CC) $(CFLAGS) -o $@ $< -lz -lm

clean:
	rm -f export_q_residual_dosage rare_sfs_pairwise

test: all
	@echo "Compile test passed: 2 binaries (export_q_residual_dosage, rare_sfs_pairwise)"

.PHONY: all clean test
