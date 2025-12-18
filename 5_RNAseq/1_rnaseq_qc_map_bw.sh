#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RNA-seq Step 1
#   1) Read trimming and QC with fastp
#   2) Alignment to hg38 with HISAT2 + sorting/indexing with samtools
#   3) BPM-normalized bigWig generation with bamCoverage
#
# Expected folder structure:
#   RNAseq/
#     rawdata/
#       PT1/PT1_mRNA_R1_1.fq.gz
#       PT1/PT1_mRNA_R1_2.fq.gz
#       PT1/PT1_mRNA_R2_1.fq.gz
#       PT1/PT1_mRNA_R2_2.fq.gz
#       ...
#       PT3/PT3_mRNA_R1_1.fq.gz
#       PT3/PT3_mRNA_R1_2.fq.gz
#       PT3/PT3_mRNA_R2_1.fq.gz
#       PT3/PT3_mRNA_R2_2.fq.gz
#
# Outputs:
#   1_cleandata/<sample>/<rep>/*.fastq.gz      (trimmed reads + fastp logs)
#   2_hisat2_mapping/<sample>/<rep>/*.bam      (sorted, indexed BAM)
#   3_bamCoverage/<sample>/<rep>/*.bw          (BPM-normalized bigWig)
###############################################################################

ROOT_DIR="$(pwd)"

# HISAT2 index prefix for hg38
HISAT2_INDEX="/mnt/f/zer/reference/hg38/index/hisat2/hg38"

# Threads
FASTP_THREADS=16
HISAT2_THREADS=20
SAMTOOLS_THREADS=10
BAMCOV_THREADS=20

# Samples:
#   3 replicates: R1, R2, R3
SAMPLES_3REP=(PT1 PT2 LP1 KMS11 U266 MM1S RPMI8226)
#   2 replicates: R1, R2
SAMPLES_2REP=(PT3)

ALL_SAMPLES=("${SAMPLES_3REP[@]}" "${SAMPLES_2REP[@]}")

###############################################################################
# Helper function: return replicate names for a given sample
###############################################################################
get_repeats() {
    local s="$1"
    # 2-replicate samples (PT3)
    for x in "${SAMPLES_2REP[@]}"; do
        [[ "$x" == "$s" ]] && { echo "R1 R2"; return; }
    done
    # Default: 3 replicates
    echo "R1 R2 R3"
}

###############################################################################
# Prepare output directories
###############################################################################
mkdir -p 1_cleandata 2_hisat2_mapping 3_bamCoverage

###############################################################################
# Main loop
###############################################################################
for sample in "${ALL_SAMPLES[@]}"; do
    echo "================ RNA-seq: ${sample} ================"
    REPEATS=$(get_repeats "${sample}")

    for rep in ${REPEATS}; do
        echo "[INFO] Processing ${sample} ${rep}"

        ######################################################################
        # 1) fastp: trimming and QC
        ######################################################################
        CLEAN_DIR="${ROOT_DIR}/1_cleandata/${sample}/${rep}"
        mkdir -p "${CLEAN_DIR}"

        RAW_R1="${ROOT_DIR}/rawdata/${sample}/${sample}_mRNA_${rep}_1.fq.gz"
        RAW_R2="${ROOT_DIR}/rawdata/${sample}/${sample}_mRNA_${rep}_2.fq.gz"

        OUT_R1="${CLEAN_DIR}/${sample}_mRNA_${rep}_1.fastq.gz"
        OUT_R2="${CLEAN_DIR}/${sample}_mRNA_${rep}_2.fastq.gz"

        echo "[INFO] fastp: ${RAW_R1}, ${RAW_R2}"
        fastp \
            -i "${RAW_R1}" \
            -o "${OUT_R1}" \
            -I "${RAW_R2}" \
            -O "${OUT_R2}" \
            -w "${FASTP_THREADS}" \
            --detect_adapter_for_pe \
            2> "${CLEAN_DIR}/${sample}_mRNA_${rep}.fastp.log"

        ######################################################################
        # 2) HISAT2 mapping + BAM sorting/indexing
        ######################################################################
        MAP_DIR="${ROOT_DIR}/2_hisat2_mapping/${sample}/${rep}"
        mkdir -p "${MAP_DIR}"
        cd "${MAP_DIR}"

        BAM="${sample}_${rep}.bam"

        echo "[INFO] HISAT2: ${sample} ${rep}"
        hisat2 \
            -p "${HISAT2_THREADS}" \
            -x "${HISAT2_INDEX}" \
            -1 "${OUT_R1}" \
            -2 "${OUT_R2}" \
            2> "${sample}_${rep}.hisat2.log" | \
        samtools view -bhS - | \
        samtools sort -@ "${SAMTOOLS_THREADS}" -o "${BAM}" -

        samtools index "${BAM}"

        ######################################################################
        # 3) bamCoverage: BPM-normalized bigWig
        ######################################################################
        BW_DIR="${ROOT_DIR}/3_bamCoverage/${sample}/${rep}"
        mkdir -p "${BW_DIR}"
        cd "${BW_DIR}"

        echo "[INFO] bamCoverage: ${sample} ${rep}"
        bamCoverage \
            -p "${BAMCOV_THREADS}" \
            --normalizeUsing BPM \
            -b "${MAP_DIR}/${BAM}" \
            -o "${sample}_${rep}.bw"

        cd "${ROOT_DIR}"
    done
done

echo "RNA-seq Step 1 (fastp, HISAT2, bamCoverage) finished."
