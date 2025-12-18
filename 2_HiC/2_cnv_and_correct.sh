#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# CNV from Hi-C / Micro-C (mcool) + CNV-based normalization (correct-cnv)
#
# Steps:
#   1) calculate-cnv  : from .mcool to CNV profile (bedGraph)
#   2) segment-cnv    : segment CNV with proper ploidy
#   3) plot-cnv       : optional genome-wide CNV plot (for QC)
#   4) correct-cnv    : normalize mcool using segmented CNV
#
# Input mcool (from previous script 1_hic_to_mcool.sh):
#   2_get_hic_mcool/<sample>/<sample>_contact.mcool
#
# Outputs (renamed to 3/4/5/6* because 1/2 are used already):
#   3_calculate-cnv/<sample>/<res>/<sample>_<res>.CNV-profile.bedGraph
#   4_segment-cnv  /<sample>/<res>/<sample>_<res>.CNV-seg.bedGraph
#   5_plot-cnv     /<sample>/<res>/<sample>_<res>.CNV.genome-wide.png
#   6_correct-cnv  /<sample>/log  (correct-cnv modifies mcool in-place)
###############################################################################

############################## CONFIG #########################################

# In situ Hi-C samples (digested by MboI)
INSITU_SAMPLES=(LP1 KMS11 U266 MM1S RPMI8226)

# Micro-C samples (no restriction enzyme, use "uniform")
MICROC_SAMPLES=(PT1 PT2 PT3)

# Samples expected to be triploid
PLOIDY3_SAMPLES=(KMS11 LP1 RPMI8226)

# Samples expected to be diploid (in situ)
PLOIDY2_INSITU_SAMPLES=(MM1S U266)

# Micro-C patients: assume diploid (ploidy 2)
PLOIDY2_MICROC_SAMPLES=(PT1 PT2 PT3)

# Resolutions to use for CNV calling / correction
# You can add more (e.g. 5000 10000 25000) if needed
RESOLUTIONS=(50000)

# Project root
ROOT_DIR="$(pwd)"

# mcool directory (from previous Hi-C / Micro-C pipeline)
MCOOL_DIR="${ROOT_DIR}/2_get_hic_mcool"

# Genome / enzyme for calculate-cnv
GENOME="hg38"

# Enzyme for in situ Hi-C
INSITU_ENZYME="MboI"

# Enzyme for Micro-C (no restriction site, use "uniform")
MICROC_ENZYME="uniform"

###############################################################################
#                         HELPER FUNCTIONS                                    #
###############################################################################

is_in_array() {
    local elem="$1"; shift
    local arr=("$@")
    for a in "${arr[@]}"; do
        [[ "$a" == "$elem" ]] && return 0
    done
    return 1
}

get_ploidy() {
    local sample="$1"

    # Triploid samples
    if is_in_array "${sample}" "${PLOIDY3_SAMPLES[@]}"; then
        echo 3
        return
    fi

    # Diploid in situ samples
    if is_in_array "${sample}" "${PLOIDY2_INSITU_SAMPLES[@]}"; then
        echo 2
        return
    fi

    # Diploid Micro-C samples
    if is_in_array "${sample}" "${PLOIDY2_MICROC_SAMPLES[@]}"; then
        echo 2
        return
    fi

    # Default
    echo 2
}

get_enzyme() {
    local sample="$1"

    if is_in_array "${sample}" "${INSITU_SAMPLES[@]}"; then
        echo "${INSITU_ENZYME}"
        return
    fi

    if is_in_array "${sample}" "${MICROC_SAMPLES[@]}"; then
        echo "${MICROC_ENZYME}"
        return
    fi

    # Default: in situ enzyme
    echo "${INSITU_ENZYME}"
}

###############################################################################
#                           MAIN LOGIC                                        #
###############################################################################

# 1) calculate-cnv + segment-cnv + plot-cnv for all samples
ALL_SAMPLES=("${INSITU_SAMPLES[@]}" "${MICROC_SAMPLES[@]}")

for sample in "${ALL_SAMPLES[@]}"; do
    echo "================ CNV pipeline for ${sample} ================"

    # Determine enzyme and ploidy for this sample
    ENZ=$(get_enzyme "${sample}")
    PLOIDY=$(get_ploidy "${sample}")

    for r in "${RESOLUTIONS[@]}"; do
        echo "[INFO] (${sample}) resolution ${r} bp"

        # Path to mcool resolution
        MCOOL_PATH="${MCOOL_DIR}/${sample}/${sample}_contact.mcool::/resolutions/${r}"

        ########################
        # 1. calculate-cnv
        ########################
        echo "[INFO] (${sample}) calculate-cnv (enzyme=${ENZ})"
        CNV_OUT_DIR="${ROOT_DIR}/3_calculate-cnv/${sample}/${r}"
        mkdir -p "${CNV_OUT_DIR}"

        calculate-cnv \
            -H "${MCOOL_PATH}" \
            -g "${GENOME}" \
            -e "${ENZ}" \
            --output "${CNV_OUT_DIR}/${sample}_${r}.CNV-profile.bedGraph" \
            >> "${CNV_OUT_DIR}/log" 2>&1

        ########################
        # 2. segment-cnv
        ########################
        echo "[INFO] (${sample}) segment-cnv (ploidy=${PLOIDY})"
        SEG_OUT_DIR="${ROOT_DIR}/4_segment-cnv/${sample}/${r}"
        mkdir -p "${SEG_OUT_DIR}"

        segment-cnv \
            --cnv-file "${CNV_OUT_DIR}/${sample}_${r}.CNV-profile.bedGraph" \
            --binsize "${r}" \
            --ploidy "${PLOIDY}" \
            --output "${SEG_OUT_DIR}/${sample}_${r}.CNV-seg.bedGraph" \
            --nproc 10 \
            >> "${SEG_OUT_DIR}/log" 2>&1

        ########################
        # 3. plot-cnv (optional but useful for QC)
        ########################
        echo "[INFO] (${sample}) plot-cnv"
        PLOT_OUT_DIR="${ROOT_DIR}/5_plot-cnv/${sample}/${r}"
        mkdir -p "${PLOT_OUT_DIR}"

        plot-cnv \
            --cnv-profile "${CNV_OUT_DIR}/${sample}_${r}.CNV-profile.bedGraph" \
            --cnv-segment "${SEG_OUT_DIR}/${sample}_${r}.CNV-seg.bedGraph" \
            --output-figure-name "${PLOT_OUT_DIR}/${sample}_${r}.CNV.genome-wide.png" \
            --dot-size 0.5 \
            --dot-alpha 0.2 \
            --line-width 1 \
            --boundary-width 0.5 \
            --label-size 7 \
            --tick-label-size 6 \
            --clean-mode \
            >> "${PLOT_OUT_DIR}/log" 2>&1

    done

    echo "================ CNV pipeline for ${sample} DONE ================"
done

###############################################################################
# 2) correct-cnv: CNV-based normalization of mcool for all samples
###############################################################################

for sample in "${ALL_SAMPLES[@]}"; do
    echo "================ correct-cnv for ${sample} ================""

    CORR_OUT_DIR="${ROOT_DIR}/6_correct-cnv/${sample}"
    mkdir -p "${CORR_OUT_DIR}"

    for r in "${RESOLUTIONS[@]}"; do
        echo "[INFO] (${sample}) correct-cnv at ${r} bp"
        MCOOL_PATH="${MCOOL_DIR}/${sample}/${sample}_contact.mcool::/resolutions/${r}"
        SEG_FILE="${ROOT_DIR}/4_segment-cnv/${sample}/${r}/${sample}_${r}.CNV-seg.bedGraph"

        correct-cnv \
            -H "${MCOOL_PATH}" \
            --cnv-file "${SEG_FILE}" \
            --nproc 10 \
            -f \
            >> "${CORR_OUT_DIR}/log" 2>&1
    done

    echo "================ correct-cnv for ${sample} DONE ================"
done

echo "All CNV calculation, segmentation, plotting and correction finished."
