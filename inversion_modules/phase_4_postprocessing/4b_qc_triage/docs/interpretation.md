# Interpreting the Q04 diagnostic figure

Five tracks per chromosome, top to bottom:

| # | Track | Source | What it measures |
|---|-------|--------|------------------|
| 1 | **Robust Z** | precomp `dt$robust_z` | Local PCA outlier score relative to chromosome median |
| 2 | **SNPs / window** | precomp `dt$n_snps` | Number of polymorphic sites in the 100-SNP window (this is the sites that entered the PCA, not raw bp density) |
| 3 | **Uncertain fraction** | Q02 | Mean across samples of the fraction of sites where `max(P00, P01, P11) < 0.9`. Low = confident calls; high = caller lacked data |
| 4 | **Mean coverage (colored by CV across samples)** | Q03 | Blue line = uniform coverage across samples; red line = highly variable across samples |
| 5 | **# samples with low coverage** | Q03 | Count of samples whose coverage in this bin falls below 50% of their own chromosome-median |

## Decision tree for a Z plateau

Walk down the tracks at the shelf region.

**1. Is Z actually a plateau or does it have texture?**
   - Flat → proceed to track 2
   - Textured (variable, peaks/valleys within shelf) → **likely real signal, stop**

**2. Does the SNP count drop in the shelf?**
   - Yes, to <50% of chromosome median → **low-polymorphism region**. Weak local PCA, signal plausibly from a handful of strong SNPs. Not necessarily an artifact but interpretation is fragile.
   - No → proceed to track 3

**3. Does BEAGLE uncertainty rise in the shelf?**
   - Yes, 2×+ above baseline → **caller confidence is low**. Posteriors are flat; PCA picks up residual bias. This is a calling-quality artifact. Flag region, note in Methods.
   - No → proceed to track 4

**4. Does coverage drop AND CV stay low (blue line)?**
   - Coverage drops uniformly across samples (low CV = blue) → **mappability / repeat region**. Every sample affected the same way. The PCA signal comes from consistent read-mapping bias. Artifact.
   - Coverage drops but CV high (red) → **sample-specific drop**. Maybe a batch effect, maybe a subset of samples carrying a large deletion. Cross-check with track 5.
   - No coverage drop → proceed to track 5

**5. Does track 5 show many low-coverage samples (tall bars in shelf)?**
   - Yes, ~half the cohort → **real deletion polymorphism**. Half the fish have the region, half don't. This is a biological signal.
   - Yes, a small fraction of samples (<10%) → **minor variant**, probably not driving the Z signal.
   - No → signal must be from another mechanism.

## If tracks 2–5 are all clean and Z is still plateau

Then the shelf is real biology of some kind — introgression block, population substructure, or a moderate-frequency inversion with a uniform full-length signal. Move to the scrubber to identify which samples are clustering together.

## Console summary metrics

The Q04 script prints a shelf-vs-reference comparison:

```
Z flatness (shelf SD/|mean|): 0.04
```
- **< 0.10** → flat plateau, needs a mechanism (check tracks 2–5)
- **0.10–0.15** → borderline
- **> 0.15** → textured profile, likely real signal

```
n_snps mean       shelf=25  ref=95   ratio=0.26x
BEAGLE uncertain  shelf=0.18 ref=0.06  ratio=3.0x
coverage mean     shelf=4.5 ref=9.1   ratio=0.49x
coverage CV       shelf=0.08 ref=0.12 ratio=0.67x
```

This particular pattern (low SNPs + high uncertainty + low coverage + low CV) is the signature of a **repeat-heterochromatin mapping artifact**: every sample fails the same way, so the PCA finds the same null axis at every window. No samples are singled out; it's a uniform failure.

The alternative pattern (normal SNPs + normal uncertainty + normal coverage + high CV on some samples) would point to real biological variation.

## Writing this up

A clean Methods paragraph would read:

> Local-PCA Z plateaus were classified as artifact versus biological signal using
> four independent QC tracks: (i) SNP density per local-PCA window, (ii) BEAGLE
> posterior-uncertainty rate, (iii) per-bin mean coverage across samples, and
> (iv) per-sample coverage deviation. Regions with concurrent low SNP density
> (<50% of chromosome median), elevated uncertainty (>2× baseline), and uniform
> coverage depression (low cross-sample CV) were flagged as low-information
> mapping artifacts and excluded from the final inversion catalog.
