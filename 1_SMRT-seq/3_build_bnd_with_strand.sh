#!/usr/bin/env bash
set -euo pipefail

######################### CONFIG ####################################

# Path to longrange_vcf_to_bedpe.py (edit this)
VCF2BED_PATH="./"

# Directory containing filtered SV VCFs:
# Expected filename pattern: 4-filtersv/<sample>_SVs_hg38.vcf
SV_VCF_DIR="4-filtersv"

# Directory containing BND coordinate pairs:
# Expected file: 5-bnd_pairs/<sample>/<sample>_bnd.bed
BND_PAIR_DIR="5-bnd_pairs"

# Output directory (starts with 6- as requested)
OUT_DIR="6-bnd_with_strand"
mkdir -p "${OUT_DIR}"

# All samples
SAMPLES=(PT1 PT2 PT3 LP1 KMS11 U266 MM1S RPMI8226)

#####################################################################

for sample in "${SAMPLES[@]}"; do
    echo "================ ${sample} ================"

    # 1) Input SV VCF for this sample
    #    Example: 4-filtersv/PT1_SVs_hg38.vcf
    SV_VCF="${SV_VCF_DIR}/${sample}_SVs_hg38.vcf"

    if [[ ! -f "${SV_VCF}" ]]; then
        echo "[WARN] SV VCF not found for ${sample}: ${SV_VCF}" >&2
        continue
    fi

    # 2) Input BND coordinate pairs (4 columns: chr1 pos1 chr2 pos2)
    #    Example: 5-bnd_pairs/PT1/PT1_bnd.bed
    BND_IN="${BND_PAIR_DIR}/${sample}/${sample}_bnd.bed"
    if [[ ! -f "${BND_IN}" ]]; then
        echo "[WARN] BND pair file not found for ${sample}: ${BND_IN}" >&2
        continue
    fi

    # 3) Intermediate files:
    #    - SV_BEDPE: all SVs converted to bedpe (CSV)
    #    - BND_STRAND_SRC: only BND from bedpe, with strand info (TAB)
    SV_BEDPE="${OUT_DIR}/${sample}_sv.bedpe"                          # CSV
    BND_STRAND_SRC="${OUT_DIR}/${sample}_bnd_with_strand_source.bedpe" # TAB, 12 columns

    echo "[INFO] (${sample}) 1/3: Converting VCF to bedpe (all SV types)"
    "${VCF2BED_PATH}/longrange_vcf_to_bedpe.py" \
        -input "${SV_VCF}" \
        -out "${SV_BEDPE}"

    echo "[INFO] (${sample}) 2/3: Extracting BND with strand information from bedpe"

    # Input:  SV_BEDPE (CSV with header:
    #         chrom1,start1,stop1,chrom2,start2,stop2,variant_name,score,strand1,strand2,variant_type,split)
    # Output: BND_STRAND_SRC (TAB, 12 columns, with chr prefix added back)
    awk -F',' 'NR>1 && $11=="BND" {
        printf "chr%s\t%s\t%s\tchr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
               $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
    }' "${SV_BEDPE}" > "${BND_STRAND_SRC}"

    # 4) Final output: merge BND pairs with strand info
    #    Format: chr1 pos1 chr2 pos2 strand1 strand2 (TAB-separated)
    BND_OUT="${OUT_DIR}/${sample}_bnd_with_strand.bed"

    echo "[INFO] (${sample}) 3/3: Mapping strand info to BND pairs"

    awk -v FS="\t" -v OFS="\t" '
        # First pass: read BND_STRAND_SRC (12-column bedpe with strands)
        NR == FNR {
            # key = chr1:pos1:chr2:pos2
            key = $1 ":" $2 ":" $4 ":" $5
            strands[key] = $9 "\t" $10
            next
        }
        # Second pass: read BND_IN (4-column BND pairs: chr1 pos1 chr2 pos2)
        {
            key  = $1 ":" $2 ":" $3 ":" $4       # forward key
            key2 = $3 ":" $4 ":" $1 ":" $2       # reverse key

            if (key in strands) {
                # Forward match: keep strand order as is
                print $1, $2, $3, $4, strands[key]
            } else if (key2 in strands) {
                # Reverse match: BND is in opposite direction;
                # swap strand1 and strand2 when outputting
                split(strands[key2], a, "\t")
                print $1, $2, $3, $4, a[2], a[1]
            } else {
                # Not found (optional warning)
                print "WARNING: no strand found for " $1, $2, $3, $4 > "/dev/stderr"
            }
        }
    ' "${BND_STRAND_SRC}" "${BND_IN}" > "${BND_OUT}"

    echo "[OK] (${sample}) Final BND+strand file: ${BND_OUT}"
done

echo "All samples finished. Results are in: ${OUT_DIR}"
