#!/usr/bin/env bash
# =============================================================================
# run_theta_all_scales.sh
# =============================================================================
# Run build_theta_supp_tables.py on all available theta scales.
#
# Produces per-scale output directories:
#   13_theta_pi_supp/win500000_step500000/   <- use this for the MS numbers
#   13_theta_pi_supp/win50000_step10000/     <- fine-scale (not for tests)
#   13_theta_pi_supp/win10000_step2000/      <- finer
#   13_theta_pi_supp/win5000_step1000/       <- finest (noisy)
#
# IMPORTANT: only the 500kb non-overlapping scale gives statistically valid
# quantile/percentile summaries (windows are independent). The sliding-window
# scales are for PLOTTING ONLY — do not use their outlier lists for statistics.
# =============================================================================

set -euo pipefail

SCRIPT="${1:-build_theta_supp_tables.py}"   # path to the python script
THETA_BASE="${2:-02_heterozygosity/03_theta}" # base dir containing pestPG files
SAMPLES="${3:-01_inputs_check/samples_qcpass.txt}"
OUT_BASE="${4:-13_theta_pi_supp}"

mkdir -p "$OUT_BASE"

# Scale -> (window, step, min_sites, is_main)
# is_main=1 means this scale is recommended for statistical summaries
declare -a SCALES=(
    # main scale — non-overlapping, used for manuscript numbers
    "500000:500000:1000:main:."
    # fine scales — sliding, for plotting only
    "50000:10000:500:fine:multiscale"
    "10000:2000:300:fine:multiscale"
    "5000:1000:200:fine:multiscale"
)

for entry in "${SCALES[@]}"; do
    IFS=':' read -r win step min_sites kind subdir <<< "$entry"
    scale_dir="$OUT_BASE/win${win}_step${step}"
    mkdir -p "$scale_dir"

    if [[ "$subdir" == "." ]]; then
        theta_dir="$THETA_BASE"
    else
        theta_dir="$THETA_BASE/$subdir"
    fi

    echo ""
    echo "============================================================"
    echo "Scale: win=$win step=$step  [$kind]  min_sites=$min_sites"
    echo "Input: $theta_dir"
    echo "Output: $scale_dir"
    echo "============================================================"

    if ! ls "$theta_dir"/*.win${win}.step${step}.pestPG >/dev/null 2>&1; then
        echo "  (no pestPG files at this scale — skipping)"
        continue
    fi

    python3 "$SCRIPT" \
        --theta-dir "$theta_dir" \
        --window "$win" --step "$step" \
        --out-dir "$scale_dir" \
        --samples "$SAMPLES" \
        --min-sites "$min_sites"

    # Write a README marking which scale is the statistical one
    cat > "$scale_dir/README.txt" <<EOF
Theta window scale: win=${win} step=${step}  [${kind}]

$(if [[ "$kind" == "main" ]]; then
    echo "This is the MAIN scale — non-overlapping windows."
    echo "Use this for:"
    echo "  - Manuscript Results paragraph numbers (mean, median, IQR of θ_π)"
    echo "  - Supplementary Table ST (per-chromosome summary)"
    echo "  - Outlier window list (statistically valid because windows are independent)"
else
    echo "This is a FINE-SCALE sliding window — for plotting only."
    echo "DO NOT use for:"
    echo "  - Quantile / percentile summaries (windows overlap, not independent)"
    echo "  - Kruskal-Wallis or any test assuming independent observations"
    echo "Use for:"
    echo "  - Fine-scale ideogram plots of θ_π"
    echo "  - Inspecting specific genomic regions (e.g. inversion candidates)"
fi)

Minimum sites per window filter: $min_sites
Generated: $(date -u +'%Y-%m-%dT%H:%M:%SZ')
EOF

done

echo ""
echo "============================================================"
echo "DONE. Output tree:"
echo "============================================================"
find "$OUT_BASE" -maxdepth 2 -type f | sort
