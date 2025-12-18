#!/usr/bin/env bash
set -euo pipefail

export NUMEXPR_MAX_THREADS=30

# SAMPLES=(KMS11 LP1 MM1S RPMI8226 U266 PT1 PT2 PT3)
SAMPLES=(PT3)

# pwd：TransFinder/6_Integration
HERE="$(pwd)"


MCOOL_ROOT="../2_HiC/2_get_hic_mcool"
BND_ROOT="./2_intersection"
OUT_ROOT="./4_complex_bnd"

mkdir -p "${OUT_ROOT}"

for sample in "${SAMPLES[@]}"; do
    echo "===== Processing ${sample} ====="

    BND_FILE="${BND_ROOT}/${sample}/${sample}_transfinder_bnd.tsv"
    SAMPLE_OUT="${OUT_ROOT}/${sample}"

    if [[ ! -f "${BND_FILE}" ]]; then
        echo "[WARN] Missing BND file for ${sample}, skip."
        continue
    fi

    mkdir -p "${SAMPLE_OUT}"

    assemble-complexSVs \
        -O "${SAMPLE_OUT}/${sample}" \
        -B "${BND_FILE}" \
        --balance-type CNV \
        --protocol insitu \
        --nproc 30 \
        -H \
        ${MCOOL_ROOT}/${sample}/${sample}_contact.mcool::resolutions/25000 \
        ${MCOOL_ROOT}/${sample}/${sample}_contact.mcool::resolutions/10000 \
        ${MCOOL_ROOT}/${sample}/${sample}_contact.mcool::resolutions/5000

    echo "===== ${sample} DONE ====="
done
