#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SMRT-seq inter-chromosomal translocation (BND) calling – CLR + HiFi
#
# This script:
#   1) Aligns CLR/HiFi PacBio reads to a reference genome using pbmm2
#   2) Calls structural variants with pbsv (discover + call)
#   3) Performs basic quality filtering on SV calls
#   4) Extracts BND (SVTYPE=BND) events
#   5) Keeps only inter-chromosomal BND events (chr1 != chr2), based on pbsv ID
#
# Outputs (per sample, in 4-filtersv/):
#   - <sample>_hg38_chr1_22x_50bp_filter.vcf          (filtered SV VCF, all SV types)
#   - <sample>_translocation.vcf                      (BND-only VCF, with header)
#   - <sample>_inter-translocation.vcf                (inter-chromosomal BND only, NO header)
#
# Requirements:
#   - pbmm2
#   - pbsv
#
# Author: (your name)
# Project: TransFinder
###############################################################################


############################ CONFIG ###########################################

# Reference genome FASTA
REF_FASTA="/mnt/f/zer/reference/hg38/fa/hg38chr1-22xym/hg38_chr1_22xym.p13.fa"

# Output root directory (default: current working directory)
OUT_DIR="${PWD}"

# CLR raw subreads BAM directory
CLR_RAW_DIR="/mnt/g/long_read/rawdata/longDNA_rawdata/Genome/subreads"

# HiFi raw BAM directory
HIFI_RAW_DIR="/mnt/g/DATA-BACKUP/mm-patient/pacbio"

# CLR sample IDs
CLR_SAMPLES=("LP1" "KMS11" "U266" "MM1S" "RPMI8226")

# HiFi sample IDs
HIFI_SAMPLES=("PT1" "PT2" "PT3")

# Subdirectories for intermediate outputs
PBMM2_DIR="${OUT_DIR}/1-pbmm2"
SVSIG_DIR="${OUT_DIR}/2-svsig"
CALLSV_DIR="${OUT_DIR}/3-callsv"
FILTER_DIR="${OUT_DIR}/4-filtersv"   # keep original naming

###############################################################################


mkdir -p "${PBMM2_DIR}" "${SVSIG_DIR}" "${CALLSV_DIR}" "${FILTER_DIR}"


run_clr_sample() {
    local sample="$1"
    echo "[CLR] Processing sample: ${sample}"

    # 1) Align CLR subreads to the reference genome
    pbmm2 align \
        "${REF_FASTA}" \
        "${CLR_RAW_DIR}/${sample}_subreads.bam" \
        "${PBMM2_DIR}/${sample}_hg38_chr1_22xym.bam" \
        --median-filter \
        --preset SUBREAD \
        --sort -j 24

    # 2) Generate SV signatures with pbsv (discover)
    pbsv discover \
        "${PBMM2_DIR}/${sample}_hg38_chr1_22xym.bam" \
        "${SVSIG_DIR}/${sample}_hg38_chr1_22xym.svsig.gz"

    # 3) Call SVs with pbsv (call)
    pbsv call -m 50 \
        "${REF_FASTA}" \
        "${SVSIG_DIR}/${sample}_hg38_chr1_22xym.svsig.gz" \
        "${CALLSV_DIR}/${sample}_hg38_chr1_22xym_50bp.vcf"

    # 4) Basic quality filtering:
    #    - keep header lines
    #    - remove chrM and chrY (optional)
    #    - keep only PASS
    #    - remove IMPRECISE calls
    local raw_vcf="${CALLSV_DIR}/${sample}_hg38_chr1_22xym_50bp.vcf"
    local filt_vcf="${FILTER_DIR}/${sample}_SVs_hg38.vcf"

    awk '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            # Skip chrM and chrY if present (adjust if needed)
            if ($1 == "chrM" || $1 == "chrY") next;

            # Only keep PASS
            if ($7 != "PASS") next;

            # Remove IMPRECISE calls
            if ($0 ~ /IMPRECISE/) next;

            print;
        }
    ' "${raw_vcf}" > "${filt_vcf}"

    # 5) Extract BND (SVTYPE=BND) from filtered VCF
    local bnd_vcf="${FILTER_DIR}/${sample}_translocation.vcf"

    awk '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            # INFO is in column 8
            split($8, a, ";")
            # Assuming SVTYPE is the first tag, e.g. SVTYPE=BND
            if (a[1] == "SVTYPE=BND")
                print
        }
    ' "${filt_vcf}" > "${bnd_vcf}"

    # 6) Keep only inter-chromosomal BND:
    #    - We parse chr2 from the pbsv ID and require chr1 != chr2
    #    - NO header lines are written to the final VCF
    local inter_vcf="${FILTER_DIR}/${sample}_inter-translocation.vcf"

    awk '
        BEGIN { OFS="\t" }
        !/^#/ {
            # VCF col1 = chr1
            chr1 = $1
            # VCF col3 = ID, e.g. pbsv.BND.chr16:79932630-chr22:22922727
            id   = $3

            # Split on "-" to get second part: chr2:pos2
            split(id, a, "-")
            split(a[2], b, ":")
            chr2 = b[1]

            if (chr1 != chr2)
                print
        }
    ' "${bnd_vcf}" > "${inter_vcf}"

    echo "[CLR] Inter-chromosomal BND VCF (no header) written: ${inter_vcf}"
}


run_hifi_sample() {
    local sample="$1"
    echo "[HiFi] Processing sample: ${sample}"

    # 1) Align HiFi reads to the reference genome
    pbmm2 align \
        "${REF_FASTA}" \
        "${HIFI_RAW_DIR}/${sample}/${sample}.bam" \
        "${PBMM2_DIR}/${sample}_hg38_chr1_22xym.bam" \
        --preset HIFI \
        --sort -j 20

    # 2) Generate SV signatures with pbsv (HiFi mode)
    pbsv discover --hifi \
        "${PBMM2_DIR}/${sample}_hg38_chr1_22xym.bam" \
        "${SVSIG_DIR}/${sample}_hg38_chr1_22xym.svsig.gz"

    # 3) Call SVs with pbsv (HiFi mode)
    pbsv call -m 50 --hifi \
        "${REF_FASTA}" \
        "${SVSIG_DIR}/${sample}_hg38_chr1_22xym.svsig.gz" \
        "${CALLSV_DIR}/${sample}_hg38_chr1_22xym_50bp.vcf"

    # 4) Basic quality filtering (HiFi):
    #    - keep header lines
    #    - keep only PASS
    #    - remove IMPRECISE calls
    local raw_vcf="${CALLSV_DIR}/${sample}_hg38_chr1_22xym_50bp.vcf"
    local filt_vcf="${FILTER_DIR}/${sample}_SVs_hg38.vcf"

    awk '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            # Skip chrM and chrY if present (adjust if needed)
            if ($1 == "chrM" || $1 == "chrY") next;

            # Only keep PASS
            if ($7 != "PASS") next;

            # Remove IMPRECISE calls
            if ($0 ~ /IMPRECISE/) next;

            print;
        }
    ' "${raw_vcf}" > "${filt_vcf}"

    # 5) Extract BND (SVTYPE=BND)
    local bnd_vcf="${FILTER_DIR}/${sample}_translocation.vcf"

    awk '
        BEGIN { OFS="\t" }
        /^#/ { print; next }
        {
            split($8, a, ";")
            if (a[1] == "SVTYPE=BND")
                print
        }
    ' "${filt_vcf}" > "${bnd_vcf}"

    # 6) Keep only inter-chromosomal BND, no header lines
    local inter_vcf="${FILTER_DIR}/${sample}_inter-translocation.vcf"

    awk '
        BEGIN { OFS="\t" }
        !/^#/ {
            chr1 = $1
            id   = $3

            split(id, a, "-")
            split(a[2], b, ":")
            chr2 = b[1]

            if (chr1 != chr2)
                print
        }
    ' "${bnd_vcf}" > "${inter_vcf}"

    echo "[HiFi] Inter-chromosomal BND VCF (no header) written: ${inter_vcf}"
}


############################ MAIN #############################################

echo "=== Calling inter-chromosomal translocations for CLR samples ==="
for s in "${CLR_SAMPLES[@]}"; do
    run_clr_sample "${s}"
done

echo "=== Calling inter-chromosomal translocations for HiFi samples ==="
for s in "${HIFI_SAMPLES[@]}"; do
    run_hifi_sample "${s}"
done

echo "All SMRT-seq inter-chromosomal BND calling completed."
echo "Final inter-BND VCFs (no header): ${FILTER_DIR}/*_inter-translocation.vcf"
