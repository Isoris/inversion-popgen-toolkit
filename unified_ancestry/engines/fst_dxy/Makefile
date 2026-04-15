CC      = gcc
CFLAGS  = -O3 -march=native -fopenmp -Wall
LDFLAGS = -lz -lm -fopenmp

all: region_popstats rare_sfs_pairwise export_q_residual_dosage

region_popstats: region_popstats.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

rare_sfs_pairwise: rare_sfs_pairwise.c
	$(CC) $(CFLAGS) -o $@ $< -lz -lm

export_q_residual_dosage: export_q_residual_dosage.c
	$(CC) $(CFLAGS) -o $@ $< -lz -lm

clean:
	rm -f region_popstats rare_sfs_pairwise export_q_residual_dosage

test: all
	@echo "Compile test passed: 3 binaries"

.PHONY: all clean test
