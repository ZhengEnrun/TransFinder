#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Build unique BND pairs and BED-like coordinates from translocation VCFs
#
# For each sample, this script:
#   1) Reads <cell>_translocation.vcf (BND-only VCF)
#   2) Extracts the pbsv BND ID (e.g. pbsv.BND.chr1:pos1-chr2:pos2 → chr1:pos1-chr2:pos2)
#   3) Deduplicates reciprocal pairs (A-B vs B-A) using an inline Python script
#   4) Normalizes pair order by chromosome index and sorts
#   5) Outputs a 4-column TSV/BED-like file: chr1  pos1  chr2  pos2
#
# Inputs:
#   6-svtype_vcf/<cell>/<cell}_translocation.vcf
#
# Outputs (per cell, in ${OUT_DIR}/${cell}/):
#   1_bnd_ID                 # raw BND IDs: chrA:posA-chrB:posB
#   2_pair_bnd               # deduplicated pairs: chrA:posA-chrB:posB
#   3_pair_bnd.sort          # normalized and sorted pairs
#   <cell>_bnd.bed           # chr1  pos1  chr2  pos2
#
# Author: (your name)
# Project: TransFinder
###############################################################################


############################ CONFIG ###########################################

# Project root (adjust as needed)
BASE_DIR="$(pwd)"

# Directory containing per-sample translocation VCFs
TRANS_VCF_DIR="${BASE_DIR}/4-filtersv"

# Output directory for BND pairs
OUT_DIR="${BASE_DIR}/5-bnd_pairs"

# List of samples to process
SAMPLES=("LP1" "KMS11" "U266" "MM1S" "RPMI8226" "PT1" "PT2" "PT3")

###############################################################################


mkdir -p "${OUT_DIR}"

for cell in "${SAMPLES[@]}"; do
    echo "[BND] Processing sample: ${cell}"

    sample_vcf="${TRANS_VCF_DIR}/${cell}_inter-translocation.vcf"
    sample_out_dir="${OUT_DIR}/${cell}"

    if [[ ! -f "${sample_vcf}" ]]; then
        echo "  [WARN] VCF not found for ${cell}: ${sample_vcf}, skip."
        continue
    fi

    mkdir -p "${sample_out_dir}"
    cd "${sample_out_dir}"

    ###########################################################################
    # 1) Extract BND IDs from VCF
    #
    # pbsv ID format example:
    #   pbsv.BND.chr16:79932630-chr22:22922727
    # We keep only the third field after splitting by ".", i.e.:
    #   chr16:79932630-chr22:22922727
    ###########################################################################

    awk '
        BEGIN { OFS="\t" }
        !/^#/ {
            split($3, a, ".");
            print a[3];
        }
    ' "${sample_vcf}" > 1_bnd_ID

    ###########################################################################
    # 2) Deduplicate reciprocal BND pairs using inline Python
    #
    # Input:  1_bnd_ID      (one line per BND: "chrA:posA-chrB:posB")
    # Output: 2_pair_bnd    (deduplicated pairs, A-B/B-A kept once)
    ###########################################################################

    python - << 'PY'
import pandas as pd

# Read BND ID list
df = pd.read_csv("1_bnd_ID", header=None, names=["translocation"], sep="\t")

seen = set()
unique_pairs = []

for row in df["translocation"]:
    if not isinstance(row, str):
        continue
    parts = row.split("-")
    if len(parts) != 2:
        continue
    left, right = parts[0], parts[1]

    # Normalize pair so that (A-B) and (B-A) are treated as the same
    key = tuple(sorted([left, right]))
    if key not in seen:
        seen.add(key)
        unique_pairs.append(f"{key[0]}-{key[1]}")

# Write deduplicated pairs
out_df = pd.DataFrame({"translocation": unique_pairs})
out_df.to_csv("2_pair_bnd", index=False, header=False)
PY

    ###########################################################################
    # 3) Normalize chromosome order and sort pairs
    #
    # Each line in 2_pair_bnd: "chrA:posA-chrB:posB"
    # We enforce that the smaller chromosome index is in the first position.
    #
    # Trick:
    #   "chr16" → split by "hr" → ["c", "16"] → use 16 for numeric comparison.
    ###########################################################################

    awk -F"-" -v OFS="-" '
    {
        # $1 = chrA:posA, $2 = chrB:posB
        split($1, a, ":");
        split(a[1], b, "hr");  # b[2] = chromosome index (e.g. "16")

        split($2, c, ":");
        split(c[1], d, "hr");  # d[2] = chromosome index

        # Ensure smaller chromosome index comes first
        if (b[2] <= d[2]) {
            print $1, $2
        } else {
            print $2, $1
        }
    }
    ' 2_pair_bnd \
    | sort -t ':' -k1,1 \
    > 3_pair_bnd.sort

    ###########################################################################
    # 4) Convert "chrA:posA-chrB:posB" → "chrA  posA  chrB  posB"
    #
    # This is a 4-column BED-like file that can be used for downstream
    # integration with Hi-C, enhancer peaks, genes, etc.
    ###########################################################################

    sed 's/:/\t/g; s/-/\t/g' 3_pair_bnd.sort > "${cell}_bnd.bed"

    echo "  [BND] Finished sample ${cell}. Output: ${sample_out_dir}/${cell}_bnd.bed"
done

echo "All samples processed. BND pair files are in: ${OUT_DIR}"
