#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RNA-seq Step 2
#   Gene-level TPM/FPKM quantification per replicate using StringTie
#
# Input BAM:
#   2_hisat2_mapping/<sample>/<rep>/<sample>_<rep>.bam
#
# Outputs:
#   7_stringtie_tpm/<sample>/<sample>_<rep>.tpm   (gene table with Coverage/FPKM/TPM)
#   7_stringtie_tpm/<sample>/<sample>_<rep>.gtf   (transcript model for this replicate)
###############################################################################

ROOT_DIR="$(pwd)"

GTF="/mnt/f/zer/reference/hg38/gencode.v40.annotation.gtf"
THREADS=20

SAMPLES_3REP=(PT1 PT2 LP1 KMS11 U266 MM1S RPMI8226)
SAMPLES_2REP=(PT3)
ALL_SAMPLES=("${SAMPLES_3REP[@]}" "${SAMPLES_2REP[@]}")

###############################################################################
# Helper: return replicate names for a given sample
###############################################################################
get_repeats() {
    local s="$1"
    for x in "${SAMPLES_2REP[@]}"; do
        [[ "$x" == "$s" ]] && { echo "R1 R2"; return; }
    done
    echo "R1 R2 R3"
}

###############################################################################
# Main loop
###############################################################################
for sample in "${ALL_SAMPLES[@]}"; do
    REPEATS=$(get_repeats "${sample}")
    echo "================ StringTie: ${sample} ================"

    OUT_DIR="${ROOT_DIR}/7_stringtie_tpm/${sample}"
    mkdir -p "${OUT_DIR}"

    for rep in ${REPEATS}; do
        echo "[INFO] StringTie: ${sample} ${rep}"

        BAM="${ROOT_DIR}/2_hisat2_mapping/${sample}/${rep}/${sample}_${rep}.bam"

        cd "${OUT_DIR}"

        stringtie \
            -p "${THREADS}" \
            -G "${GTF}" \
            -o "${sample}_${rep}.gtf" \
            -e \
            -A "${sample}_${rep}.tpm" \
            -l "${sample}" \
            "${BAM}"

        cd "${ROOT_DIR}"
    done
done

echo "RNA-seq Step 2 (StringTie TPM/FPKM) finished."
