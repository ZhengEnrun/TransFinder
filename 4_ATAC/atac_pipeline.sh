#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# ATAC-seq pipeline (single replicate R1)
#
# Steps per sample:
#   1) Adapter trimming and QC (fastp + fastqc)
#   2) Alignment with bowtie2, filtering and sorting (samtools)
#   3) PCR duplicate removal (picard MarkDuplicates) and Tn5 shift
#   4) BigWig generation (bamCoverage)
#   5) Peak calling (macs2)
#
# Expected rawdata file pattern (paired-end):
#   rawdata/<SAMPLE>_1.fq.gz
#   rawdata/<SAMPLE>_2.fq.gz
#
# Example:
#   rawdata/PT1_ATAC_R1_1.fq.gz
#   rawdata/PT1_ATAC_R1_2.fq.gz
#   -> SAMPLE = PT1_ATAC_R1
#
# Outputs:
#   1_cleandata/   : cleaned fastq (fastp)
#   2_bam/         : aligned + dedup + shifted BAM
#   3_bw/          : BigWig tracks (RPGC normalized)
#   4_macs2/       : macs2 peaks
#   fastqc/        : fastqc reports on cleaned fastq
###############################################################################

############################## CONFIG #########################################

ROOT_DIR="$(pwd)"
REF_GENOME="hg38"

# bowtie2 index prefix (bowtie2 -x)
BOWTIE2_INDEX="/mnt/f/zer/hg38_chr1-22x/2-hic/2-hicpro/1-data/hg38/hg38"

# effective genome size for hg38 (for bamCoverage --normalizeUsing RPGC)
EFFECTIVE_GENOME_SIZE=2913022398

# Threads
FASTP_THREADS=16
FASTQC_THREADS=20
BOWTIE2_THREADS=24
SAMTOOLS_SORT_THREADS=10
BAMCOVERAGE_THREADS=24

###############################################################################
#                         PREPARE DIRECTORIES                                 #
###############################################################################

mkdir -p 1_cleandata 2_bam 3_bw 4_macs2 fastqc

###############################################################################
#                         DETECT SAMPLES (R1 ONLY)                            #
###############################################################################

# Take files ending with "R1_1.fq.gz" as sample identifiers
# Example: rawdata/PT1_ATAC_R1_1.fq.gz -> PT1_ATAC_R1
mapfile -t SAMPLES < <(ls rawdata/*R1_1.fq.gz 2>/dev/null \
  | xargs -n1 basename \
  | sed 's/_1\.fq\.gz$//' \
  | sort -u)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "[ERROR] No samples detected in rawdata/*R1_1.fq.gz" >&2
  exit 1
fi

echo "Detected ATAC-seq samples (R1 only):"
printf '  %s\n' "${SAMPLES[@]}"

echo "Pipeline for ATAC-seq (single replicate R1)"

###############################################################################
#                           MAIN LOOP                                         #
###############################################################################

for sample in "${SAMPLES[@]}"; do
  echo "================ Processing sample: ${sample} ================"

  ############################################################################
  # 1) Adapter trimming with fastp
  ############################################################################
  RAW_R1="rawdata/${sample}_1.fq.gz"
  RAW_R2="rawdata/${sample}_2.fq.gz"
  CLEAN_R1="1_cleandata/${sample}_1.fq.gz"
  CLEAN_R2="1_cleandata/${sample}_2.fq.gz"

  echo "[INFO] fastp: ${RAW_R1}, ${RAW_R2}"
  fastp \
    -i "${RAW_R1}" \
    -o "${CLEAN_R1}" \
    -I "${RAW_R2}" \
    -O "${CLEAN_R2}" \
    --detect_adapter_for_pe \
    -w "${FASTP_THREADS}" \
    2> "${sample}.fastp.log"

  ############################################################################
  # 2) fastqc on cleaned fastq files
  ############################################################################
  echo "[INFO] fastqc on cleaned fastq for ${sample}"
  fastqc "${CLEAN_R1}" "${CLEAN_R2}" -t "${FASTQC_THREADS}" -o fastqc

  ############################################################################
  # 3) Alignment with bowtie2, filter and sort BAM
  ############################################################################
  cd "${ROOT_DIR}/2_bam"

  BAM_SORT="${sample}.${REF_GENOME}.sort.bam"

  echo "[INFO] bowtie2 + samtools for ${sample}"
  bowtie2 \
    -p "${BOWTIE2_THREADS}" \
    --very-sensitive \
    -X 2000 \
    -x "${BOWTIE2_INDEX}" \
    -1 "${ROOT_DIR}/${CLEAN_R1}" \
    -2 "${ROOT_DIR}/${CLEAN_R2}" \
    2> "${sample}.bowtie.log" | \
  samtools view -bh -q 30 | \
  samtools sort -@ "${SAMTOOLS_SORT_THREADS}" -o "${BAM_SORT}" -

  samtools index "${BAM_SORT}"

  ############################################################################
  # 4) Remove PCR duplicates (picard MarkDuplicates)
  ############################################################################
  BAM_NODUP_TMP="${sample}.${REF_GENOME}.sort.noDup.bam.tmp"
  BAM_NODUP="${sample}.${REF_GENOME}.sort.noDup.bam"
  METRICS_LOG="${sample}.markduplicates.log"

  echo "[INFO] picard MarkDuplicates for ${sample}"
  picard MarkDuplicates \
    QUIET=true \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=LENIENT \
    INPUT="${BAM_SORT}" \
    OUTPUT="${BAM_NODUP_TMP}" \
    METRICS_FILE="${METRICS_LOG}"

  samtools view -bh "${BAM_NODUP_TMP}" | \
  samtools sort -@ "${SAMTOOLS_SORT_THREADS}" -o "${BAM_NODUP}" -

  samtools index "${BAM_NODUP}"
  rm -f "${BAM_NODUP_TMP}"

  ############################################################################
  # 5) Tn5 shift using alignmentSieve (--ATACshift)
  ############################################################################
  BAM_SHIFT_TMP="${sample}.${REF_GENOME}.sort.noDup.shift.bam"
  BAM_SHIFT_SORT="${sample}.${REF_GENOME}.sort.noDup.shift.sort.bam"

  echo "[INFO] alignmentSieve --ATACshift for ${sample}"
  alignmentSieve \
    --numberOfProcessors 20 \
    --ATACshift \
    --bam "${BAM_NODUP}" \
    -o "${BAM_SHIFT_TMP}"

  samtools view -bh "${BAM_SHIFT_TMP}" | \
  samtools sort -@ "${SAMTOOLS_SORT_THREADS}" -o "${BAM_SHIFT_SORT}" -

  samtools index "${BAM_SHIFT_SORT}"
  rm -f "${BAM_SHIFT_TMP}"

  ############################################################################
  # 6) BigWig generation (bamCoverage)
  ############################################################################
  cd "${ROOT_DIR}/3_bw"

  BW_FILE="${sample}.bw"
  echo "[INFO] bamCoverage for ${sample}"
  bamCoverage \
    -p "${BAMCOVERAGE_THREADS}" \
    --normalizeUsing RPGC \
    --effectiveGenomeSize "${EFFECTIVE_GENOME_SIZE}" \
    -b "${ROOT_DIR}/2_bam/${BAM_SHIFT_SORT}" \
    -o "${BW_FILE}"

  ############################################################################
  # 7) Peak calling with macs2
  ############################################################################
  cd "${ROOT_DIR}/4_macs2"

  echo "[INFO] macs2 callpeak for ${sample}"
  macs2 callpeak \
    -t "${ROOT_DIR}/2_bam/${BAM_SHIFT_SORT}" \
    -n "${sample}" \
    --shift -75 \
    --extsize 150 \
    --nomodel \
    -B \
    --SPMR \
    -g hs \
    -f BAMPE \
    -q 0.05 \
    --keep-dup all

  # Back to project root
  cd "${ROOT_DIR}"

  echo "================ Sample ${sample} DONE ================"
done

echo "[INFO] Running multiqc over fastqc directory..."
multiqc fastqc --force

