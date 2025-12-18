#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Hi-C / Micro-C pipeline: raw reads -> .mcool (balanced)
#
# Part 1: in situ Hi-C (5 samples) using HiC-Pro + hicpro2juicebox + hic2cool
# Part 2: Micro-C (3 samples) using bwa + pairtools + juicer_tools + hic2cool
#
# Final outputs:
#   2_get_hic_mcool/<sample>/<sample>_contact.mcool
# for all 8 samples (5 in situ + 3 Micro-C).
###############################################################################

############################## GLOBAL CONFIG ##################################

# In situ Hi-C samples (processed with HiC-Pro)
INSITU_SAMPLES=(LP1 KMS11 U266 MM1S RPMI8226)

# Micro-C samples (processed with bwa + pairtools)
MICROC_SAMPLES=(PT1 PT2 PT3)

# Current working directory (project root)
ROOT_DIR="$(pwd)"

######## Shared tools / reference (used by BOTH in situ Hi-C and Micro-C) ####

# Unified juicer tools path
JUICER_TOOLS="/mnt/d/linux/software/HiC-Pro-master/juicer_tools_1.22.01.jar"

# Unified chrom sizes file
CHROMSIZES="/mnt/d/linux/reference/hg38/hg38.chrom.sizes.txt"

######## In situ Hi-C specific paths ########

HICPRO_BIN="/mnt/d/linux/software/HiC-Pro-master/bin/HiC-Pro"
HICPRO2JB_SH="/mnt/d/linux/software/HiC-Pro-master/bin/utils/hicpro2juicebox.sh"
INSITU_RESTRICTION_BED="/mnt/d/linux/reference/hg38/hic/hg38_dpnii.bed"
HICPRO_CONFIG="/mnt/d/linux/reference/hg38/hic/config-hicpro.txt"

# Raw data directory for in situ Hi-C (each sample in 0_rawdata/<sample>)
INSITU_RAW_DIR="${ROOT_DIR}/0_rawdata"

# HiC-Pro output directory
HICPRO_OUT_DIR="${ROOT_DIR}/1_hicpro"

# Unified mcool directory for BOTH in situ Hi-C and Micro-C
MCOOL_DIR="${ROOT_DIR}/2_get_hic_mcool"

# Resolutions to balance for in situ Hi-C
INSITU_RESOLUTIONS=(5000 10000 25000 50000)

######## Micro-C specific paths ########

BWA_INDEX="/mnt/f/zer/reference/hg38/index/bwa/hg38.p13.fa"
MICROC_RAW_DIR="/mnt/f/zer/mm_patient/6_microC/0_cleandata"
MICROC_MAP_DIR="${ROOT_DIR}/1_bwa-pairtools"

###############################################################################
#                           FUNCTIONS                                         #
###############################################################################

run_insitu_sample() {
    local sample="$1"
    echo "========== [In situ Hi-C] ${sample} =========="

    ########################
    # 1) Run HiC-Pro
    ########################

    echo "[INFO] (${sample}) Running HiC-Pro..."
    "${HICPRO_BIN}" \
        -i "${INSITU_RAW_DIR}/${sample}" \
        -o "${HICPRO_OUT_DIR}/${sample}" \
        -c "${HICPRO_CONFIG}"

    ########################
    # 2) Convert HiC-Pro allValidPairs -> .hic -> .mcool
    ########################

    mkdir -p "${MCOOL_DIR}/${sample}"
    cd "${MCOOL_DIR}/${sample}"

    local valid_pairs="${HICPRO_OUT_DIR}/${sample}/allvalidpairs_results/${sample}.allValidPairs"

    echo "[INFO] (${sample}) allValidPairs -> .hic (hicpro2juicebox)..."
    "${HICPRO2JB_SH}" \
        -i "${valid_pairs}" \
        -g "${CHROMSIZES}" \
        -j "${JUICER_TOOLS}" \
        -r "${INSITU_RESTRICTION_BED}" \
        -t ./tmp \
        -o ./

    # This should generate: ${sample}.allValidPairs.hic
    echo "[INFO] (${sample}) .hic -> .mcool..."
    hic2cool convert "${sample}.allValidPairs.hic" "${sample}_contact.mcool" -r 0 -p 30
    echo "[INFO] (${sample}) hic2cool done."

    ########################
    # 3) Balance .mcool at selected resolutions
    ########################
    for r in "${INSITU_RESOLUTIONS[@]}"; do
        echo "[INFO] (${sample}) cooler balance at ${r} bp..."
        cooler balance "${sample}_contact.mcool::/resolutions/${r}" -p 30
    done

    cd "${ROOT_DIR}"
    echo "========== [In situ Hi-C] ${sample} DONE =========="
}


run_microc_sample() {
    local sample="$1"
    echo "========== [Micro-C] ${sample} =========="

    ########################
    # 1) Mapping + pairtools pipeline
    ########################

    mkdir -p "${MICROC_MAP_DIR}/${sample}"
    cd "${MICROC_MAP_DIR}/${sample}"
    mkdir -p temp

    local fq_path="${MICROC_RAW_DIR}/${sample}"

    echo "[INFO] (${sample}) Running bwa + pairtools + samtools..."

    bwa mem -5SP -T0 -t 12 "${BWA_INDEX}" \
        <(zcat "${fq_path}/${sample}_1_R1.fastq.gz" "${fq_path}/${sample}_2_R1.fastq.gz") \
        <(zcat "${fq_path}/${sample}_1_R2.fastq.gz" "${fq_path}/${sample}_2_R2.fastq.gz") | \
    pairtools parse \
        --min-mapq 40 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --nproc-in 12 \
        --nproc-out 12 \
        --chroms-path "${CHROMSIZES}" | \
    pairtools sort \
        --tmpdir="${MICROC_MAP_DIR}/${sample}/temp/" \
        --nproc 12 | \
    pairtools dedup \
        --backend cython \
        --nproc-in 12 \
        --nproc-out 12 \
        --mark-dups \
        --output-stats stats.txt | \
    pairtools split \
        --nproc-in 12 \
        --nproc-out 12 \
        --output-pairs "${sample}_mapped.pairs" \
        --output-sam - | \
    samtools view -bS -@ 12 | \
    samtools sort -@ 12 -o "${sample}_mapped.bam"

    samtools index "${sample}_mapped.bam"

    echo "[INFO] (${sample}) Running Micro-C QC..."
    python3 /mnt/f/zer/software/micro-C/Micro-C-main/get_qc.py \
        -p stats.txt > "${sample}_qc.log"

    cd "${ROOT_DIR}"

    ########################
    # 2) pairs -> .hic -> .mcool in 2_get_hic_mcool/
    ########################

    mkdir -p "${MCOOL_DIR}/${sample}"
    cd "${MCOOL_DIR}/${sample}"

    echo "[INFO] (${sample}) pairs -> .hic (juicer tools)..."
    java -Xmx48000m -Djava.awt.headless=true -jar "${JUICER_TOOLS}" pre \
        --threads 20 \
        "${MICROC_MAP_DIR}/${sample}/${sample}_mapped.pairs" \
        "${sample}_contact_map.hic" \
        "${CHROMSIZES}"

    echo "[INFO] (${sample}) .hic -> .mcool..."
    hic2cool convert "${sample}_contact_map.hic" "${sample}_contact.mcool" -r 0 -p 20
    echo "[INFO] (${sample}) hic2cool done."

    echo "[INFO] (${sample}) Balancing all resolutions in .mcool..."
    resolutions=$(cooler ls "${sample}_contact.mcool")
    for r in ${resolutions}; do
        echo "[INFO] (${sample}) cooler balance ${r}"
        cooler balance "${r}" -p 20
    done

    cd "${ROOT_DIR}"
    echo "========== [Micro-C] ${sample} DONE =========="
}

###############################################################################
#                                   MAIN                                      #
###############################################################################

echo "########## In situ Hi-C samples ##########"
for s in "${INSITU_SAMPLES[@]}"; do
    run_insitu_sample "${s}"
done

echo "########## Micro-C samples ##########"
for m in "${MICROC_SAMPLES[@]}"; do
    run_microc_sample "${m}"
done

echo "All Hi-C / Micro-C processing to .mcool finished."
