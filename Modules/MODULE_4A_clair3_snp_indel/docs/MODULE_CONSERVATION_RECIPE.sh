#!/usr/bin/env bash
# =============================================================================
# MODULE_CONSERVATION — Conservation & deleterious prediction for catfish
#
# This module builds TWO independent scores for the VEF framework:
#
# A. SIFT4G (protein-level)
#    - Predicts whether amino acid substitutions affect protein function
#    - Based on sequence homology across UniRef90
#    - Works on ANY species with a genome + gene annotation (GTF)
#    - Input: your catfish genome FASTA + GTF + VCF
#    - Output: per-variant SIFT score (0-1, <0.05 = deleterious)
#
# B. Conservation scores (genome-level)
#    - GERP++, phastCons, phyloP from multi-species alignment
#    - Requires: whole-genome alignment of ~15-30 fish species
#    - Pipeline: download genomes → Cactus/Progressive alignment → HAL →
#                hal2maf → GERP++ / phastCons / phyloP
#    - Output: per-base conservation scores (BED/bigWig)
#
# This script is a RECIPE — it generates the actual pipeline scripts
# for your HPC.  Read it, adjust paths, then run each step.
# =============================================================================

cat << 'RECIPE'

╔══════════════════════════════════════════════════════════════╗
║         CONSERVATION SCORING PIPELINE FOR CATFISH           ║
║                                                             ║
║  Part A: SIFT4G (easy, ~1 day compute)                     ║
║  Part B: Multi-species GERP/phastCons (hard, ~1 week)      ║
╚══════════════════════════════════════════════════════════════╝

RECIPE

# ─────────────────────────────────────────────────────────────
# PART A: SIFT4G — Protein-level deleterious prediction
# ─────────────────────────────────────────────────────────────
#
# SIFT4G works on ANY organism.  It needs:
#   1. Reference genome FASTA
#   2. Gene annotation GTF
#   3. UniRef90 protein database (downloaded once)
#   4. Your VCF(s) to annotate
#
# Steps:
#   A1. Install SIFT4G
#   A2. Download UniRef90
#   A3. Build SIFT4G database for your catfish genome
#   A4. Annotate your VCF
#
# Total time: ~12-24h for database build, minutes for annotation
# ─────────────────────────────────────────────────────────────

cat << 'SIFT_INSTALL'

# === A1. Install SIFT4G ===

# Clone repos
cd /scratch/lt200308-agbsci/13-programs
git clone --recursive https://github.com/rvaser/sift4g.git
cd sift4g && make -j 16
# Binary: sift4g/bin/sift4g

# Database builder
cd /scratch/lt200308-agbsci/13-programs
git clone https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB.git

# Annotator
wget https://github.com/pauline-ng/SIFT4G_Annotator/raw/master/SIFT4G_Annotator.jar \
     -O /scratch/lt200308-agbsci/13-programs/SIFT4G_Annotator.jar

SIFT_INSTALL


cat << 'SIFT_DB'

# === A2. Download UniRef90 ===
# ~30 GB download, needed once

UNIREF_DIR="/scratch/lt200308-agbsci/databases/uniref90"
mkdir -p "$UNIREF_DIR"
cd "$UNIREF_DIR"
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
# Index it
/scratch/lt200308-agbsci/13-programs/sift4g/bin/sift4g \
    --generate-db uniref90.fasta --out uniref90.sift4g_db

SIFT_DB


cat << 'SIFT_BUILD_DB'

# === A3. Build SIFT4G database for your catfish genome ===
# This creates pre-computed SIFT predictions for every possible
# amino acid substitution at every coding position.

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
SIFT_DIR="${BASE}/MODULE_CONSERVATION/sift4g"
mkdir -p "$SIFT_DIR"

# Create config file
cat > "${SIFT_DIR}/catfish_sift4g.config" << EOF
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard
MITO_GENETIC_CODE_TABLE=2
MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial

PARENT_DIR=${SIFT_DIR}
ORG=fClaHyb
ORG_VERSION=v1

#FASTA and GFF/GTF — your genome and annotation
GENOME_FASTA=${BASE}/00-samples/fClaHyb_Gar_LG.fa
GENE_ANNOTATION=${BASE}/00-samples/fClaHyb_Gar_LG.gtf

DBNAME=fClaHyb_v1

PROTEIN_DB=/scratch/lt200308-agbsci/databases/uniref90/uniref90.fasta
SIFT4G_PATH=/scratch/lt200308-agbsci/13-programs/sift4g/bin/sift4g

# Threads
SIFT4G_THREADS=32
EOF

# Build the database (submit as SLURM job)
# This takes ~12-24h depending on genome size
cat > "${SIFT_DIR}/build_sift4g_db.sh" << 'JOBEOF'
#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1 -n 32 --mem=64G
#SBATCH -t 2-00:00:00
#SBATCH -A lt200308
#SBATCH -J sift4g_build
set -euo pipefail

SIFT_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_CONSERVATION/sift4g"
SIFT_SCRIPTS="/scratch/lt200308-agbsci/13-programs/SIFT4G_Create_Genomic_DB"

cd "$SIFT_SCRIPTS"
perl make-SIFT-db-all.pl -config "${SIFT_DIR}/catfish_sift4g.config"
JOBEOF

echo "Submit: sbatch ${SIFT_DIR}/build_sift4g_db.sh"

SIFT_BUILD_DB


cat << 'SIFT_ANNOTATE'

# === A4. Annotate your VCF ===
# After database is built, annotate SNPs in seconds

BASE="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04"
SIFT_DIR="${BASE}/MODULE_CONSERVATION/sift4g"
SIFT_DB="${SIFT_DIR}/fClaHyb/v1"
VCF_DIR="${BASE}/MODULE_4A_SNP_INDEL50_Clair3/vcf"

# Annotate per-chromosome merged VCF
# (or use your per-sample VCFs)
for VCF in "${VCF_DIR}"/C_gar_LG01/*.vcf.gz; do
    SAMPLE=$(basename "$VCF" | cut -d. -f1)
    java -jar /scratch/lt200308-agbsci/13-programs/SIFT4G_Annotator.jar \
        -c -i "$VCF" \
        -d "$SIFT_DB" \
        -r "${SIFT_DIR}/results/${SAMPLE}" \
        -t
done

# Output: per-sample TSV with SIFT_SCORE and SIFT_PREDICTION columns
# SIFT_SCORE < 0.05 → DELETERIOUS
# SIFT_SCORE >= 0.05 → TOLERATED

SIFT_ANNOTATE


cat << 'SIFT_PARSE'

# === A5. Parse SIFT results into VEF format ===
# This Python snippet converts SIFT annotator output to the VEF ESM/SIFT input format

python3 << 'PYEOF'
import os, glob
sift_results = "/path/to/sift4g/results"
outfile = "/path/to/sift_scores_for_vef.tsv"

with open(outfile, 'w') as out:
    out.write("VAR_KEY\tSIFT_SCORE\tSIFT_PREDICTION\n")
    for tsv in glob.glob(os.path.join(sift_results, "*", "*SIFTpredictions.tsv")):
        with open(tsv) as f:
            header = f.readline()  # skip
            for line in f:
                p = line.strip().split('\t')
                if len(p) < 15: continue
                chrom, pos, ref, alt = p[0], p[1], p[3], p[4]
                sift_score = p[12] if len(p) > 12 else "."
                sift_pred  = p[13] if len(p) > 13 else "."
                vk = f"{chrom}:{pos}:{ref}:{alt}"
                out.write(f"{vk}\t{sift_score}\t{sift_pred}\n")
PYEOF

SIFT_PARSE


# ─────────────────────────────────────────────────────────────
# PART B: Multi-species alignment → GERP / phastCons / phyloP
# ─────────────────────────────────────────────────────────────
#
# This is the hard part.  Steps:
#   B1. Download ~20 fish genomes from NCBI
#   B2. Build species tree
#   B3. Run Cactus progressive alignment → HAL
#   B4. hal2maf → per-chromosome MAF
#   B5. Run GERP++ on MAF
#   B6. Run phastCons / phyloP on MAF (PHAST package)
#   B7. Convert to per-base BED/bigWig scores
#
# Practical time: ~3-7 days on HPC
# ─────────────────────────────────────────────────────────────


cat << 'SPECIES_LIST'

# === B1. Species selection ===
# We need ~15-25 genomes at varying evolutionary distances:
#   - Close: other Clarias (C. fuscus, C. magur, C. batrachus)
#   - Medium: other Siluriformes families
#   - Distant: other Ostariophysi + outgroups (zebrafish, medaka)
#
# Recommended species and NCBI accessions:

# CLOSE (Clariidae + close Siluriformes)
# 1.  Clarias fuscus          GCA_039881885.1  (Hong Kong catfish, chr-level)
# 2.  Clarias magur            GCA_023372715.1  (walking catfish)
# 3.  Ictalurus punctatus      GCF_001660625.3  (channel catfish, RefSeq)
# 4.  Ictalurus furcatus       GCA_022429615.1  (blue catfish, chr-level)
# 5.  Pelteobagrus fulvidraco  GCA_004120215.1  (yellow catfish)
# 6.  Silurus meridionalis     GCA_014805685.1  (southern catfish)
# 7.  Pangasianodon hypophthalmus GCA_009078355.1 (striped catfish)

# MEDIUM (other Siluriformes families)
# 8.  Malapterurus electricus  GCA_030060505.1  (electric catfish)
# 9.  Bagarius yarrelli        GCA_004120335.1  (giant devil catfish)
# 10. Kryptopterus vitreolus   GCA_041355785.1  (glass catfish)
# 11. Pterocryptis cochinchinensis GCA_054855095.1

# DISTANT (Ostariophysi outgroups)
# 12. Danio rerio              GCF_000002035.6  (zebrafish, gold standard)
# 13. Astyanax mexicanus       GCF_023375975.1  (Mexican tetra)
# 14. Pygocentrus nattereri    GCF_015220715.1  (red-bellied piranha)
# 15. Electrophorus electricus GCF_013358815.1  (electric eel)

# OUTGROUP (non-Ostariophysi teleosts)
# 16. Oryzias latipes          GCF_002234675.1  (medaka)
# 17. Oreochromis niloticus    GCF_001858045.2  (Nile tilapia)
# 18. Gasterosteus aculeatus   GCF_006229115.1  (stickleback)
# 19. Lepisosteus oculatus     GCF_000242695.1  (spotted gar, deep outgroup)

SPECIES_LIST


cat << 'DOWNLOAD_GENOMES'

# === B1b. Download genomes ===

GENOMES_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_CONSERVATION/genomes"
mkdir -p "$GENOMES_DIR"

# Example download script (repeat for each species):
download_genome() {
    local ACC="$1" NAME="$2"
    local URL="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${ACC}/download?include_annotation_type=GENOME_FASTA"
    mkdir -p "${GENOMES_DIR}/${NAME}"
    cd "${GENOMES_DIR}/${NAME}"
    datasets download genome accession "$ACC" --include genome
    unzip ncbi_dataset.zip
    # Softmask or not — Cactus prefers softmasked
    cp ncbi_dataset/data/${ACC}/*_genomic.fna "${NAME}.fa"
    samtools faidx "${NAME}.fa"
    echo "Downloaded: ${NAME}"
}

# Run for each species:
download_genome GCA_039881885.1 Clarias_fuscus
download_genome GCF_001660625.3 Ictalurus_punctatus
download_genome GCF_000002035.6 Danio_rerio
# ... etc for all ~19 species

DOWNLOAD_GENOMES


cat << 'CACTUS_ALIGNMENT'

# === B3. Progressive Cactus alignment ===
# Cactus produces a HAL file containing the whole-genome alignment
#
# Install: mamba install -c conda-forge -c bioconda cactus
# OR use Singularity: singularity pull docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.3

CONS_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_CONSERVATION"
GENOMES="${CONS_DIR}/genomes"

# Species tree (Newick format) — approximate from TimeTree or literature
# Adjust branch lengths and topology to match actual phylogeny
cat > "${CONS_DIR}/species_tree.nwk" << 'TREE'
((((((Clarias_gariepinus_x_macrocephalus:5,Clarias_fuscus:5):15,Clarias_magur:20):30,((Ictalurus_punctatus:25,Ictalurus_furcatus:25):15,Pelteobagrus_fulvidraco:40):10):20,(Silurus_meridionalis:55,Pangasianodon_hypophthalmus:55):15):30,(Malapterurus_electricus:80,((Bagarius_yarrelli:50,Pterocryptis_cochinchinensis:50):20,Kryptopterus_vitreolus:70):10):20):40,(((Danio_rerio:100,Astyanax_mexicanus:100):20,(Pygocentrus_nattereri:90,Electrophorus_electricus:90):30):30,((Oryzias_latipes:110,Oreochromis_niloticus:110):20,(Gasterosteus_aculeatus:120,Lepisosteus_oculatus:120):10):20):10);
TREE

# Seqfile mapping species names to FASTA paths
cat > "${CONS_DIR}/seqfile.txt" << 'SEQFILE'
Clarias_gariepinus_x_macrocephalus /path/to/fClaHyb_Gar_LG.fa
Clarias_fuscus /path/to/genomes/Clarias_fuscus/Clarias_fuscus.fa
Ictalurus_punctatus /path/to/genomes/Ictalurus_punctatus/Ictalurus_punctatus.fa
Danio_rerio /path/to/genomes/Danio_rerio/Danio_rerio.fa
# ... add all species
SEQFILE

# Run Cactus (this is the big job — 1-3 days with ~64 cores)
# Use --binariesMode local if not using Docker/Singularity
cactus \
    "${CONS_DIR}/cactus_jobstore" \
    "${CONS_DIR}/seqfile.txt" \
    "${CONS_DIR}/fish_alignment.hal" \
    --maxCores 64 \
    --maxMemory 256G

CACTUS_ALIGNMENT


cat << 'GERP_PHASTCONS'

# === B4-B6. Extract MAF and run GERP / phastCons ===

CONS_DIR="/scratch/lt200308-agbsci/Quentin_project_KEEP_2026-02-04/MODULE_CONSERVATION"
HAL="${CONS_DIR}/fish_alignment.hal"
REF_GENOME="Clarias_gariepinus_x_macrocephalus"

# B4. HAL → MAF (per chromosome)
mkdir -p "${CONS_DIR}/maf"
while read CHROM LEN REST; do
    hal2maf "$HAL" "${CONS_DIR}/maf/${CHROM}.maf" \
        --refGenome "$REF_GENOME" \
        --refSequence "$CHROM" \
        --noAncestors
done < "${CONS_DIR}/../00-samples/fClaHyb_Gar_LG.fa.fai"

# B5. GERP++ (install: mamba install -c bioconda gerp)
mkdir -p "${CONS_DIR}/gerp"
# GERP needs a neutral tree — estimate from 4-fold degenerate sites
# For now, use the species tree with branch lengths as proxy
for MAF in "${CONS_DIR}/maf"/*.maf; do
    CHROM=$(basename "$MAF" .maf)
    gerpcol -t "${CONS_DIR}/species_tree.nwk" -f "$MAF" \
        -e "$REF_GENOME" -j \
        -a > "${CONS_DIR}/gerp/${CHROM}.gerp.rates"
done

# Convert GERP to BED
python3 << 'PYEOF'
import os, glob
gerp_dir = os.environ.get("CONS_DIR", ".") + "/gerp"
for rates_file in glob.glob(os.path.join(gerp_dir, "*.gerp.rates")):
    chrom = os.path.basename(rates_file).replace(".gerp.rates", "")
    bed_out = rates_file.replace(".rates", ".bed")
    with open(rates_file) as f, open(bed_out, 'w') as out:
        for i, line in enumerate(f):
            p = line.strip().split('\t')
            if len(p) >= 2:
                # GERP: neutral_rate, RS_score
                rs = float(p[1])
                out.write(f"{chrom}\t{i}\t{i+1}\t{rs:.4f}\n")
PYEOF

# B6. phastCons / phyloP (install: mamba install -c bioconda phast)
mkdir -p "${CONS_DIR}/phast"

# Step 1: Fit neutral model from 4-fold degenerate sites
# Extract 4-fold sites from your annotation
# (this requires msa_view and phyloFit from PHAST)

# For each chromosome MAF:
for MAF in "${CONS_DIR}/maf"/*.maf; do
    CHROM=$(basename "$MAF" .maf)

    # phyloFit: estimate neutral model
    phyloFit "$MAF" --tree "${CONS_DIR}/species_tree.nwk" \
        --subst-mod REV --out-root "${CONS_DIR}/phast/${CHROM}.neutral"

    # phastCons: conservation scores
    phastCons "$MAF" "${CONS_DIR}/phast/${CHROM}.neutral.mod" \
        --most-conserved "${CONS_DIR}/phast/${CHROM}.most_conserved.bed" \
        --score > "${CONS_DIR}/phast/${CHROM}.phastcons.wig"

    # phyloP: per-base conservation p-values
    phyloP "${CONS_DIR}/phast/${CHROM}.neutral.mod" "$MAF" \
        --method LRT --mode CONACC \
        > "${CONS_DIR}/phast/${CHROM}.phylop.wig"
done

# Convert WIG to BED for easy intersection
for WIG in "${CONS_DIR}/phast"/*.phastcons.wig; do
    CHROM=$(basename "$WIG" .phastcons.wig)
    awk -v c="$CHROM" '/^[0-9]/{print c"\t"NR-1"\t"NR"\t"$1}' "$WIG" \
        > "${CONS_DIR}/phast/${CHROM}.phastcons.bed"
done

GERP_PHASTCONS


cat << 'INTEGRATE_VEF'

# === B7. Integrate conservation scores into VEF ===

# After SIFT, GERP, and phastCons are computed, merge them
# into a single score file for the VEF pipeline.

python3 << 'PYEOF'
"""
Merge conservation scores into VEF-compatible format.

Output: conservation_scores.tsv with columns:
  VAR_KEY  SIFT_SCORE  SIFT_PRED  GERP_RS  PHASTCONS  PHYLOP
"""
import os

# Paths — adjust to your setup
sift_file    = "sift_scores_for_vef.tsv"
gerp_dir     = "gerp/"
phastcons_dir = "phast/"

# Load SIFT
sift = {}
if os.path.isfile(sift_file):
    with open(sift_file) as f:
        f.readline()
        for line in f:
            p = line.strip().split('\t')
            sift[p[0]] = {"score": p[1], "pred": p[2]}

# Load GERP (indexed by chrom:pos)
gerp = {}
for bed in os.listdir(gerp_dir):
    if not bed.endswith(".gerp.bed"): continue
    with open(os.path.join(gerp_dir, bed)) as f:
        for line in f:
            p = line.strip().split('\t')
            key = f"{p[0]}:{int(p[1])+1}"  # convert to 1-based
            gerp[key] = p[3]

# Load phastCons
phastcons = {}
for bed in os.listdir(phastcons_dir):
    if not bed.endswith(".phastcons.bed"): continue
    with open(os.path.join(phastcons_dir, bed)) as f:
        for line in f:
            p = line.strip().split('\t')
            key = f"{p[0]}:{int(p[1])+1}"
            phastcons[key] = p[3]

# Write merged output
# This file can be passed to 20_annotate_variant_consequences.py via --esm_tsv
# (or we add a dedicated --conservation_tsv flag)
with open("conservation_scores.tsv", 'w') as out:
    out.write("VAR_KEY\tSIFT_SCORE\tSIFT_PRED\tGERP_RS\tPHASTCONS\n")
    # Iterate over all SIFT-scored variants
    for vk, s in sift.items():
        parts = vk.split(":")
        pos_key = f"{parts[0]}:{parts[1]}"
        g = gerp.get(pos_key, ".")
        pc = phastcons.get(pos_key, ".")
        out.write(f"{vk}\t{s['score']}\t{s['pred']}\t{g}\t{pc}\n")

print("Merged conservation scores written.")
PYEOF

INTEGRATE_VEF


cat << 'SUMMARY'

╔══════════════════════════════════════════════════════════════╗
║                    PRACTICAL SUMMARY                        ║
╠══════════════════════════════════════════════════════════════╣
║                                                             ║
║  EASY (do this week):                                       ║
║    SIFT4G — just needs your genome + GTF + UniRef90         ║
║    Time: ~1 day for DB build, minutes to annotate           ║
║    Result: SIFT_SCORE per missense variant                  ║
║                                                             ║
║  MEDIUM (do this month):                                    ║
║    Download 15-20 fish genomes from NCBI                    ║
║    Run Cactus progressive alignment                         ║
║    Time: ~3-5 days on 64 cores                              ║
║    Result: HAL whole-genome alignment                       ║
║                                                             ║
║  THEN (same week as Cactus):                                ║
║    GERP++ from MAF → per-base rejected substitution score   ║
║    phastCons from MAF → per-base conservation probability   ║
║    phyloP from MAF → per-base evolutionary rate test        ║
║    Time: ~1-2 days per tool                                 ║
║                                                             ║
║  INTEGRATION:                                               ║
║    Feed scores into VEF script 20/21:                       ║
║    - SIFT < 0.05 → add evidence code SU1                    ║
║    - GERP RS > 2 → add evidence code SU_GERP (new)          ║
║    - phastCons > 0.9 → add evidence code SU_CONSERVED       ║
║    - All three together → strong support for deleterious     ║
║                                                             ║
║  WHAT THIS GIVES YOU FOR THE PAPER:                         ║
║    "We computed per-site conservation across N fish          ║
║     genomes using GERP++ and phastCons, and predicted        ║
║     missense deleteriousness using SIFT4G, to classify       ║
║     catfish variants by breeding concern."                   ║
║                                                             ║
║  This is publishable in a fish genetics / aquaculture        ║
║  journal — nobody has done this for Clarias hybrids.        ║
╚══════════════════════════════════════════════════════════════╝

SUMMARY
