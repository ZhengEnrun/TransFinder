#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Step 7: BAM -> bigWig (RPGC) for pbmm2 outputs
#
# Script location:
#   1_SMRT-seq/4_run_bam2bw.sh
#
# Input BAMs:
#   1_SMRT-seq/1-pbmm2/*_hg38_chr1_22xym.bam
#
# Outputs:
#   1_SMRT-seq/7-bam2bw/<sample>.bw
#
# Requirements:
#   - samtools
#   - bamCoverage (deepTools)
###############################################################################

SMRT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PBMM2_DIR="${SMRT_DIR}/1-pbmm2"
OUT_DIR="${SMRT_DIR}/7-bam2bw"

THREADS=20
EFFECTIVE_GENOME_SIZE=2913022398   # hg38 RPGC effective genome size

mkdir -p "${OUT_DIR}"

shopt -s nullglob
BAMS=("${PBMM2_DIR}"/*_hg38_chr1_22xym.bam)

if [[ ${#BAMS[@]} -eq 0 ]]; then
  echo "[ERROR] No BAMs found in: ${PBMM2_DIR}"
  echo "        Expected pattern: *_hg38_chr1_22xym.bam"
  exit 1
fi

for bam in "${BAMS[@]}"; do
  base="$(basename "${bam}")"
  sample="${base%_hg38_chr1_22xym.bam}"
  bw="${OUT_DIR}/${sample}.bw"

  # Ensure BAM index exists
  if [[ ! -f "${bam}.bai" ]]; then
    echo "[INFO] samtools index: ${bam}"
    samtools index "${bam}"
  fi

  echo "[INFO] bamCoverage RPGC -> ${bw}"
  bamCoverage -p "${THREADS}" \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
    -b "${bam}" \
    -o "${bw}"
done

echo "[DONE] bigWigs written to: ${OUT_DIR}/*.bw"
