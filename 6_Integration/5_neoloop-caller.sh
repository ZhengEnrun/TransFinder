#!/usr/bin/env bash
set -euo pipefail

export NUMEXPR_MAX_THREADS=30

# SAMPLES=(KMS11 LP1 MM1S RPMI8226 U266 PT1 PT2 PT3)
SAMPLES=(PT3)

# 
HERE="$(pwd)"

MCOOL_ROOT="../2_HiC/2_get_hic_mcool"
ASSEMBLY_ROOT="${HERE}/4_complex_bnd"
OUT_ROOT="${HERE}/5_neoloop-caller"

mkdir -p "${OUT_ROOT}"

for sample in "${SAMPLES[@]}"; do
    echo "===== NeoLoop calling: ${sample} ====="

    ASSEMBLY_FILE="${ASSEMBLY_ROOT}/${sample}/${sample}.assemblies.txt"
    MCOOL_FILE="${MCOOL_ROOT}/${sample}/${sample}_contact.mcool"
    SAMPLE_OUT="${OUT_ROOT}/${sample}"

    if [[ ! -f "${ASSEMBLY_FILE}" ]]; then
        echo "[WARN] Missing assembly for ${sample}, skip."
        continue
    fi

    if [[ ! -f "${MCOOL_FILE}" ]]; then
        echo "[WARN] Missing mcool for ${sample}, skip."
        continue
    fi

    mkdir -p "${SAMPLE_OUT}"

    neoloop-caller \
        -O "${SAMPLE_OUT}/${sample}.neo-loops.txt" \
        --assembly "${ASSEMBLY_FILE}" \
        --balance-type CNV \
        --protocol insitu \
        --prob 0.95 \
        --nproc 30 \
        -H \
        "${MCOOL_FILE}::resolutions/25000" \
        "${MCOOL_FILE}::resolutions/10000" \
        "${MCOOL_FILE}::resolutions/5000"

    echo "===== ${sample} DONE ====="
done
