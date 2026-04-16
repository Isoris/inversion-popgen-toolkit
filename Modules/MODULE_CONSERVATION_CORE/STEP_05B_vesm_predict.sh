#!/usr/bin/env bash
#SBATCH --job-name=CONS_05B_vesm
#SBATCH --account=lt200308
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/05B_vesm_%j.out
#SBATCH --error=slurm_logs/05B_vesm_%j.err
set -euo pipefail

# =============================================================================
# STEP 05B — VESM: Protein Language Model Missense Effect Prediction
# =============================================================================
# Uses VESM (Dinh et al. 2026, Nature Methods) — co-distilled ESM family model
# that achieves state-of-the-art zero-shot variant effect prediction.
#
# NO alignment needed. NO species-specific training.
# Works on ANY protein sequence from ANY organism.
# This is the key advantage over SIFT4G: no UniRef90 download, no DB build.
#
# Method:
#   For each catfish protein, compute LLR (log-likelihood ratio) for all 19
#   possible amino acid substitutions at every position. More negative LLR =
#   more likely damaging. Then intersect with actual missense variants from
#   SnpEff to get per-variant VESM scores.
#
# Model: VESM_650M (MIT license, 650M params, ~4GB VRAM)
#   - Loaded from HuggingFace: ntranoslab/vesm
#   - Based on ESM2-650M backbone, co-distilled across ESM family
#   - Outperforms AlphaMissense on rare variants (no MAF circularity)
#
# Reference:
#   Dinh et al. (2026) "Compressing the collective knowledge of ESM into a
#   single protein language model" Nature Methods.
#   https://github.com/ntranoslab/vesm
#
# Inputs:
#   - Catfish protein FASTA (from 00_setup.sh)
#   - SnpEff missense variants (from 04_snpeff_annotate.sh)
#
# Outputs:
#   - vesm_all_proteins_llr.tsv.gz  (all possible missense LLR scores)
#   - vesm_variant_scores.tsv       (scores for actual observed variants)
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config/pipeline.config.sh"

VESM_DIR="${DIR_SIFT4G}/../05B_vesm"
mkdir -p "${VESM_DIR}" "${DIR_LOGS}"

echo "=== STEP 05B: VESM Missense Effect Prediction ==="

# ─────────────────────────────────────────
# A. Install VESM dependencies if needed
# ─────────────────────────────────────────
echo "=== A. Checking VESM dependencies ==="

python3 -c "import torch; import transformers; import huggingface_hub" 2>/dev/null || {
    echo "  Installing VESM dependencies..."
    pip install torch transformers huggingface_hub biopython tqdm --break-system-packages 2>/dev/null || \
    pip install torch transformers huggingface_hub biopython tqdm
}

# ─────────────────────────────────────────
# B. Run VESM scoring on catfish proteins
# ─────────────────────────────────────────
echo "=== B. Running VESM on catfish proteins ==="

python3 - "${REF_PEP}" "${DIR_SNPEFF}/snpeff_canonical_consequences.tsv" "${VESM_DIR}" <<'PYEOF'
"""
VESM missense variant effect prediction for catfish proteins.

1. Load VESM_650M model from HuggingFace
2. For each protein: compute LLR for all 19 substitutions at every position
3. Map actual SnpEff missense variants to VESM scores
4. Output: per-variant VESM LLR score + classification
"""
import sys, os, gzip
from collections import defaultdict

pep_fasta = sys.argv[1]
snpeff_tsv = sys.argv[2]
outdir = sys.argv[3]

# ── Load model ──
print("  Loading VESM_650M model...")
import torch
from huggingface_hub import hf_hub_download
from transformers import AutoTokenizer, EsmForMaskedLM
import torch.nn.functional as F

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"  Device: {device}")

# Download VESM weights
model_name = "VESM_650M"
base_model = 'facebook/esm2_t33_650M_UR50D'
local_dir = os.path.join(outdir, 'vesm_weights')
os.makedirs(local_dir, exist_ok=True)

vesm_weights = os.path.join(local_dir, f"{model_name}.pth")
if not os.path.exists(vesm_weights):
    print(f"  Downloading {model_name} weights from HuggingFace...")
    hf_hub_download(repo_id="ntranoslab/vesm", filename=f"{model_name}.pth", local_dir=local_dir)

# Load base ESM2-650M + VESM weights
print("  Loading model weights...")
tokenizer = AutoTokenizer.from_pretrained(base_model)
model = EsmForMaskedLM.from_pretrained(base_model).to(device)
model.load_state_dict(torch.load(vesm_weights, map_location=device), strict=False)
model.eval()

vocab = tokenizer.get_vocab()
# Standard amino acids in ESM vocabulary
AA_TOKENS = {aa: vocab[aa] for aa in 'ACDEFGHIKLMNPQRSTVWY' if aa in vocab}
print(f"  Model loaded. Vocabulary: {len(AA_TOKENS)} amino acids")

# ── Helper: compute LLR for a protein sequence ──
def compute_llrs(sequence, max_len=1022):
    """Compute log-likelihood ratios for all positions in a protein sequence.
    
    Returns: dict of (position, mut_aa) → LLR score
    LLR < 0 = mutant less likely than wildtype = potentially damaging
    More negative = more damaging
    """
    results = {}
    
    # Handle long proteins with sliding window (ESM2 max context = 1022)
    if len(sequence) <= max_len:
        windows = [(0, sequence)]
    else:
        # Sliding window with 50% overlap
        step = max_len // 2
        windows = []
        for start in range(0, len(sequence), step):
            end = min(start + max_len, len(sequence))
            windows.append((start, sequence[start:end]))
            if end == len(sequence):
                break
    
    for win_start, win_seq in windows:
        tokens = tokenizer([win_seq], return_tensors='pt').to(device)
        with torch.no_grad():
            outputs = model(**tokens)
        
        logits = outputs['logits'][0]  # (seq_len, vocab_size)
        input_ids = tokens['input_ids'][0]
        
        # Log-softmax → log probabilities
        log_probs = torch.log_softmax(logits, dim=-1)
        
        # For each position (skip special tokens at 0 and -1)
        for i in range(1, len(input_ids) - 1):
            abs_pos = win_start + i - 1  # 0-indexed position in full protein
            
            # Skip if already scored from a previous window
            if abs_pos in {pos for (pos, _) in results}:
                continue
            
            wt_token = input_ids[i].item()
            wt_logprob = log_probs[i, wt_token].item()
            
            for aa, aa_token in AA_TOKENS.items():
                if aa_token == wt_token:
                    continue  # Skip wildtype
                mut_logprob = log_probs[i, aa_token].item()
                llr = mut_logprob - wt_logprob  # negative = damaging
                results[(abs_pos, aa)] = llr
    
    return results

# ── Load protein sequences ──
print("  Loading protein sequences...")
proteins = {}  # protein_id → sequence
with open(pep_fasta) as f:
    name, seq = '', ''
    for line in f:
        if line.startswith('>'):
            if name and seq:
                # Clean: remove stop codons (*) and newlines
                proteins[name] = seq.replace('*', '').replace('.', '')
            name = line[1:].split()[0]
            seq = ''
        else:
            seq += line.strip()
    if name and seq:
        proteins[name] = seq.replace('*', '').replace('.', '')

print(f"  Proteins loaded: {len(proteins)}")

# ── Load missense variants from SnpEff ──
print("  Loading SnpEff missense variants...")
missense_variants = []  # list of (var_key, gene_id, transcript_id, aa_change)

if os.path.exists(snpeff_tsv):
    with open(snpeff_tsv) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            d = dict(zip(header, parts))
            ann = d.get('annotation', '').lower()
            if 'missense' not in ann:
                continue
            
            vk = d.get('var_key', '')
            gid = d.get('gene_id', '')
            tid = d.get('feature_id', '')
            hgvs_p = d.get('hgvs_p', '')
            
            if hgvs_p and vk:
                missense_variants.append((vk, gid, tid, hgvs_p))

print(f"  Missense variants: {len(missense_variants)}")

# ── Build transcript → protein mapping ──
# Try to map transcript IDs to protein sequences
# In gffread output, protein headers often match transcript IDs
tx_to_prot = {}
for pid, seq in proteins.items():
    # Try various ID formats
    tx_to_prot[pid] = seq
    # Also try without version suffix
    base_id = pid.split('.')[0]
    tx_to_prot[base_id] = seq

# ── Score all proteins ──
print("\n  Scoring proteins with VESM...")
all_llr_path = os.path.join(outdir, 'vesm_all_proteins_llr.tsv.gz')
n_scored = 0
n_total = len(proteins)

with gzip.open(all_llr_path, 'wt') as out:
    out.write('protein_id\tposition\twt_aa\tmut_aa\tvesm_llr\tvesm_class\n')
    
    for idx, (pid, seq) in enumerate(proteins.items()):
        if idx % 100 == 0:
            print(f"    [{idx+1}/{n_total}] Scoring {pid} ({len(seq)} aa)...")
        
        if len(seq) < 5:
            continue  # Skip very short sequences
        
        try:
            llrs = compute_llrs(seq)
        except Exception as e:
            print(f"    WARNING: Failed on {pid}: {e}")
            continue
        
        for (pos, mut_aa), llr in llrs.items():
            if pos >= len(seq):
                continue
            wt_aa = seq[pos]
            # Classification: LLR < -7 = likely_damaging, < -3 = possibly_damaging, else benign
            if llr < -7:
                vclass = 'likely_damaging'
            elif llr < -3:
                vclass = 'possibly_damaging'
            else:
                vclass = 'benign'
            
            out.write(f"{pid}\t{pos+1}\t{wt_aa}\t{mut_aa}\t{llr:.4f}\t{vclass}\n")
        
        n_scored += 1

print(f"\n  Scored {n_scored}/{n_total} proteins")
print(f"  All LLRs: {all_llr_path}")

# ── Map to actual variants ──
print("\n  Mapping VESM scores to observed missense variants...")

# Parse HGVS protein notation: p.Xxx123Yyy → (wt_aa, pos, mut_aa)
import re
AA3TO1 = {
    'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E',
    'Gly':'G','His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F',
    'Pro':'P','Ser':'S','Thr':'T','Trp':'W','Tyr':'Y','Val':'V',
    'Ter':'*','Stop':'*'
}

def parse_hgvs_p(hgvs):
    """Parse p.Met1Val → ('M', 1, 'V') or None."""
    m = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', hgvs)
    if not m:
        return None
    wt3, pos, mt3 = m.group(1), int(m.group(2)), m.group(3)
    wt1 = AA3TO1.get(wt3)
    mt1 = AA3TO1.get(mt3)
    if wt1 and mt1 and mt1 != '*':
        return (wt1, pos, mt1)
    return None

# Load all LLR scores into lookup
print("  Loading VESM scores for lookup...")
vesm_lookup = {}  # (protein_id, pos, mut_aa) → llr
with gzip.open(all_llr_path, 'rt') as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            pid, pos, wt, mut, llr = parts[0], int(parts[1]), parts[2], parts[3], float(parts[4])
            vesm_lookup[(pid, pos, mut)] = llr

print(f"  VESM lookup: {len(vesm_lookup)} scores")

# Map variants
var_scores_path = os.path.join(outdir, 'vesm_variant_scores.tsv')
n_mapped = 0
n_unmapped = 0

with open(var_scores_path, 'w') as out:
    out.write('var_key\tgene_id\ttranscript_id\thgvs_p\twt_aa\tposition\tmut_aa\t'
              'vesm_llr\tvesm_class\n')
    
    for vk, gid, tid, hgvs_p in missense_variants:
        parsed = parse_hgvs_p(hgvs_p)
        if not parsed:
            n_unmapped += 1
            continue
        
        wt_aa, pos, mut_aa = parsed
        
        # Try to find VESM score
        llr = None
        for pid_try in [tid, tid.split('.')[0], gid, gid.split('.')[0]]:
            key = (pid_try, pos, mut_aa)
            if key in vesm_lookup:
                llr = vesm_lookup[key]
                break
        
        if llr is not None:
            if llr < -7:
                vclass = 'likely_damaging'
            elif llr < -3:
                vclass = 'possibly_damaging'
            else:
                vclass = 'benign'
            out.write(f"{vk}\t{gid}\t{tid}\t{hgvs_p}\t{wt_aa}\t{pos}\t{mut_aa}\t"
                      f"{llr:.4f}\t{vclass}\n")
            n_mapped += 1
        else:
            n_unmapped += 1

print(f"\n  Mapped: {n_mapped} variants with VESM scores")
print(f"  Unmapped: {n_unmapped}")
print(f"  Output: {var_scores_path}")

# ── Summary statistics ──
print("\n  VESM score distribution for observed missense variants:")
scores = []
with open(var_scores_path) as f:
    f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            try:
                scores.append(float(parts[7]))
            except ValueError:
                pass

if scores:
    import statistics
    print(f"    N variants scored: {len(scores)}")
    print(f"    Mean LLR:  {statistics.mean(scores):.3f}")
    print(f"    Median LLR:{statistics.median(scores):.3f}")
    print(f"    Likely damaging (LLR < -7): {sum(1 for s in scores if s < -7)}")
    print(f"    Possibly damaging (-7 ≤ LLR < -3): {sum(1 for s in scores if -7 <= s < -3)}")
    print(f"    Benign (LLR ≥ -3): {sum(1 for s in scores if s >= -3)}")

PYEOF

echo "=== STEP 05B complete ==="
echo "  VESM all-protein LLRs: ${VESM_DIR}/vesm_all_proteins_llr.tsv.gz"
echo "  VESM variant scores:   ${VESM_DIR}/vesm_variant_scores.tsv"
echo ""
echo "  Feed into STEP14 via the vesm_variant_scores.tsv file."
echo "  VESM LLR < -7 = likely damaging, < -3 = possibly damaging"
