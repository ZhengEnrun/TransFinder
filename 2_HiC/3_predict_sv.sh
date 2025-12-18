#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Step 7: Call SVs from Hi-C / Micro-C using predictSV (NeoLoopFinder)
#
# Input:
#   - mcool from step 1 + 2:
#       2_get_hic_mcool/<sample>/<sample>_contact.mcool
#
# For each sample, we run predictSV using:
#   - 5 kb    matrix    (resolutions/5000)
#   - 10 kb   matrix    (resolutions/10000)
#   - 50 kb   matrix    (resolutions/50000)
#
# Output:
#   7_predictSV/<sample>/<sample>.predictsv.txt
###############################################################################

# Threads used internally by predictSV (via numexpr)
export NUMEXPR_MAX_THREADS=30

############################## CONFIG #########################################

# Micro-C samples (patients)
MICROC_SAMPLES=(PT1 PT2 PT3)

# In situ Hi-C samples (cell lines)
INSITU_SAMPLES=(LP1 KMS11 U266 MM1S RPMI8226)

# Project root (same as previous scripts)
ROOT_DIR="$(pwd)"

# Unified mcool directory (from 1_hic_to_mcool.sh + 2_cnv_and_correct.sh)
MCOOL_DIR="${ROOT_DIR}/2_get_hic_mcool"

# Output root directory for predictSV (step 7)
PREDICTSV_OUT="${ROOT_DIR}/7_predictSV"

# Genome build
GENOME="hg38"

# Resolutions used by predictSV
RES_5K=5000
RES_10K=10000
RES_50K=50000

# Output format from predictSV (fixed to NeoLoopFinder as in your code)
FORMAT="NeoLoopFinder"

# Balance type:
#   - Micro-C: use Raw matrix
#   - In situ Hi-C: use CNV-corrected matrix
BALANCE_TYPE_MICROC="Raw"
BALANCE_TYPE_INSITU="CNV"

###############################################################################
#                              FUNCTIONS                                      #
###############################################################################

run_predictsv_for_sample() {
    local sample="$1"
    local balance_type="$2"

    echo "========== predictSV for ${sample} (balance_type=${balance_type}) =========="

    # mcool paths at different resolutions
    local HIC_5K="${MCOOL_DIR}/${sample}/${sample}_contact.mcool::/resolutions/${RES_5K}"
    local HIC_10K="${MCOOL_DIR}/${sample}/${sample}_contact.mcool::/resolutions/${RES_10K}"
    local HIC_50K="${MCOOL_DIR}/${sample}/${sample}_contact.mcool::/resolutions/${RES_50K}"

    # Basic existence check
    if [[ ! -f "${MCOOL_DIR}/${sample}/${sample}_contact.mcool" ]]; then
        echo "[WARN] mcool not found for ${sample}: ${MCOOL_DIR}/${sample}/${sample}_contact.mcool" >&2
        return
    fi

    # Output directory: 7_predictSV/<sample>/
    local outdir="${PREDICTSV_OUT}/${sample}"
    mkdir -p "${outdir}"
    cd "${outdir}"

    # Output file name
    local outfile="${sample}.predictsv.txt"

    echo "[INFO] (${sample}) running predictSV..."
    predictSV \
        --hic-5k  "${HIC_5K}" \
        --hic-10k "${HIC_10K}" \
        --hic-50k "${HIC_50K}" \
        -O "${outfile}" \
        -g "${GENOME}" \
        --balance-type "${balance_type}" \
        --output-format "${FORMAT}" \
        --prob-cutoff-5k  0.8 \
        --prob-cutoff-10k 0.8 \
        --prob-cutoff-50k 0.99999

    cd "${ROOT_DIR}"

    echo "========== predictSV for ${sample} DONE =========="
}

###############################################################################
#                                   MAIN                                      #
###############################################################################

echo "########## Micro-C samples (PT1/PT2/PT3, balance=Raw) ##########"
for s in "${MICROC_SAMPLES[@]}"; do
    run_predictsv_for_sample "${s}" "${BALANCE_TYPE_MICROC}"
done

echo "########## In situ Hi-C samples (LP1/KMS11/U266/MM1S/RPMI8226, balance=CNV) ##########"
for s in "${INSITU_SAMPLES[@]}"; do
    run_predictsv_for_sample "${s}" "${BALANCE_TYPE_INSITU}"
done

echo "All predictSV jobs finished."
