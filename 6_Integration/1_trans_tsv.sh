#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Step 1: Normalize Hi-C and long-read translocations to a common TSV / BED
#
# Script location:
#   /mnt/f/zer/TransFinder/6_Integration/1_trans_tsv.sh
#
# Project root is inferred as the parent of 6_Integration:
#   /mnt/f/zer/TransFinder
#
# Inputs
#   Hi-C SV predictions (5K_combined):
#     2_HiC/7_predictSV/<sample>/<name>.predictsv.txt.CNN_SVs.5K_combined.txt
#       - PT1: PT1.predictsv.txt.CNN_SVs.5K_combined.txt
#       - PT2: PT2.predictsv.txt.CNN_SVs.5K_combined.txt
#       - PT3: PT3.predictsv.txt.CNN_SVs.5K_combined.txt
#       - others: <sample>.predictsv.txt.CNN_SVs.5K_combined.txt
#
#   Long-read BND with strand (from SMRT-seq):
#     1_SMRT-seq/6-bnd_with_strand/<sample>_bnd_with_strand.bed
#       columns:
#         1: chrA
#         2: posA
#         3: chrB
#         4: posB
#         5: strandA (+/-)
#         6: strandB (+/-)
#
# Outputs (per sample) in:
#   6_Integration/1_trans_tsv/<sample>/
#
#   Hi-C:
#     <sample>_hic_translocation.bed
#         chrA  posA  chrB  posB  ori(++, +-, -+, --)  
#     <sample>_hic.tsv
#         source  id     sample  chrA  posA  chrB  posB  strandA  strandB
#
#   Long-read:
#     <sample>_longread_translocation_with_strand.bed
#         chrA  posA  strandA  chrB  posB  strandB
#     <sample>_longread_translocation.bed
#         chrA  posA  chrB  posB
#     <sample>_longread.tsv
#         source  id       sample  chrA  posA  chrB  posB  strandA  strandB
###############################################################################

# Directory of this script (6_Integration)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root = parent of 6_Integration
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "[INFO] Script dir : ${SCRIPT_DIR}"
echo "[INFO] Project dir: ${ROOT_DIR}"

# Sample list: 5 cell lines + 3 patients
samples=("KMS11" "LP1" "MM1S" "RPMI8226" "U266" "PT1" "PT2" "PT3")

# Output root directory (under 6_Integration)
out_root="${ROOT_DIR}/6_Integration/1_trans_tsv"
mkdir -p "${out_root}"

for sample in "${samples[@]}"; do
  echo ">>> Processing sample: ${sample}"

  # Per-sample subfolder
  sample_dir="${out_root}/${sample}"
  mkdir -p "${sample_dir}"
  cd "${sample_dir}"

  ############################
  # 1. Process Hi-C translocations
  ############################

  # Hi-C predictSV result path (5K_combined), sample-dependent naming
  case "${sample}" in
    PT1)
      hic_src="${ROOT_DIR}/2_HiC/7_predictSV/PT1/PT1.predictsv.txt.CNN_SVs.5K_combined.txt"
      ;;
    PT2)
      hic_src="${ROOT_DIR}/2_HiC/7_predictSV/PT2/PT2.predictsv.txt.CNN_SVs.5K_combined.txt"
      ;;
    PT3)
      hic_src="${ROOT_DIR}/2_HiC/7_predictSV/PT3/PT3.predictsv.txt.CNN_SVs.5K_combined.txt"
      ;;
    *)
      hic_src="${ROOT_DIR}/2_HiC/7_predictSV/${sample}/${sample}.predictsv.txt.CNN_SVs.5K_combined.txt"
      ;;
  esac

  if [[ ! -f "${hic_src}" ]]; then
    echo "    [WARN] Hi-C source file not found: ${hic_src}, skip Hi-C for this sample."
  else
    # 1.1 Extract inter-chromosomal SVs:
    #     Input columns (CNN_SVs):
    #       chr1 = $1, chr2 = $2, ori = $3 (++, +-, -+, --), pos1 = $4, pos2 = $5
    #     Output temp: chr1 pos1 chr2 pos2 ori
    awk -v OFS="\t" '
      $1 != $2 {
        print $1, $4, $2, $5, $3;
      }
    ' "${hic_src}" | sort -k1,1 -k2,2n > "${sample}_hic_translocation.bed.tmp"

    # 1.2 Normalize chromosome order: ensure chrA <= chrB (by numeric part),
    #     and adjust orientation accordingly.
    #
    # Input:  chr1 pos1 chr2 pos2 ori
    # Output: chrA posA chrB posB ori_after_sort
    awk -v OFS="\t" '
      {
        chr1=$1; pos1=$2;
        chr2=$3; pos2=$4;
        ori =$5;               # e.g. "+-", "++", etc.

        # split ori into two strands
        s1 = substr(ori,1,1);
        s2 = substr(ori,2,1);

        split(chr1,a,"hr");
        split(chr2,b,"hr");

        if (a[2] <= b[2]) {
          # order unchanged
          new_ori = s1 s2;
          print chr1, pos1, chr2, pos2, new_ori;
        } else {
          # swap both breakpoints and strands
          new_ori = s2 s1;
          print chr2, pos2, chr1, pos1, new_ori;
        }
      }
    ' "${sample}_hic_translocation.bed.tmp" > "${sample}_hic_translocation.bed"

    rm -f "${sample}_hic_translocation.bed.tmp"

    # 1.3 Convert to TSV:
    #   source, id, sample, chrA, posA, chrB, posB, strandA, strandB
    #   strandA/strandB 
    awk -v OFS="\t" -v s="${sample}" '
      {
        sv_id = "HIC_" NR;
        ori   = $5;
        sA = substr(ori,1,1);
        sB = substr(ori,2,1);
        print "hic", sv_id, s, $1, $2, $3, $4, sA, sB;
      }
    ' "${sample}_hic_translocation.bed" > "${sample}_hic.tsv"

    echo "    Hi-C BED : ${sample_dir}/${sample}_hic_translocation.bed"
    echo "    Hi-C TSV : ${sample_dir}/${sample}_hic.tsv"
  fi

  ############################
  # 2. Process long-read BND (with strand info)
  ############################

  # Long-read BND with strand (6 columns: chrA posA chrB posB strandA strandB)
  lr_src="${ROOT_DIR}/1_SMRT-seq/6-bnd_with_strand/${sample}_bnd_with_strand.bed"

  if [[ ! -f "${lr_src}" ]]; then
    echo "  [WARN] Long-read source file not found: ${lr_src}, skip long-read for this sample."
    continue
  fi

  # 2.1 Select inter-chromosomal events, keep strand info
  #     Temporary output:
  #       chr1, pos1, strand1, chr2, pos2, strand2
  awk -v OFS="\t" '
    NF >= 6 && $1 != $3 {
      chr1 = $1; pos1 = $2;
      chr2 = $3; pos2 = $4;
      str1 = $5; str2 = $6;
      print chr1, pos1, str1, chr2, pos2, str2;
    }
  ' "${lr_src}" > "${sample}_longread_translocation_with_strand.tmp"

  n_lines=$(wc -l < "${sample}_longread_translocation_with_strand.tmp")
  echo "    [INFO] ${sample}: long-read inter-chr events (raw) = ${n_lines}"

  # 2.2 Normalize chromosome order: ensure chrA <= chrB,
  #     swapping coordinates and strands if needed.
  #     Output (with strand):
  #       chrA, posA, strandA, chrB, posB, strandB
  awk -v OFS="\t" '
    NF >= 6 {
      chr1=$1; pos1=$2; str1=$3;
      chr2=$4; pos2=$5; str2=$6;

      split(chr1,a,"hr");
      split(chr2,b,"hr");

      if (a[2] <= b[2]) {
        print chr1, pos1, str1, chr2, pos2, str2;
      } else {
        print chr2, pos2, str2, chr1, pos1, str1;
      }
    }
  ' "${sample}_longread_translocation_with_strand.tmp" > "${sample}_longread_translocation_with_strand.bed"

  rm -f "${sample}_longread_translocation_with_strand.tmp"

  # 2.3 Generate 4-column BED in the same coordinate format as Hi-C:
  #     chrA, posA, chrB, posB
  cut -f1,2,4,5 "${sample}_longread_translocation_with_strand.bed" \
    > "${sample}_longread_translocation.bed"

  # 2.4 Generate TSV, keeping strand information:
  #   columns:
  #     source  id       sample  chrA  posA  chrB  posB  strandA  strandB
  awk -v OFS="\t" -v s="${sample}" '
    {
      sv_id = "LR_" NR;
      print "longread", sv_id, s, $1, $2, $4, $5, $3, $6;
    }
  ' "${sample}_longread_translocation_with_strand.bed" > "${sample}_longread.tsv"

  echo "    Long-read BED (4-cols)      : ${sample_dir}/${sample}_longread_translocation.bed"
  echo "    Long-read BED (with strand) : ${sample_dir}/${sample}_longread_translocation_with_strand.bed"
  echo "    Long-read TSV (with strand) : ${sample_dir}/${sample}_longread.tsv"

done  # sample loop

echo ">>> Step 1 finished."
echo "    All results are under: ${out_root}"
