#!/usr/bin/env bash

here=/mnt/f/zer/TransFinder/

# Promoter regions (1 kb around TSS)
promoter="${here}/6_Integration/coding_gene_promoters_1kb.bed"

# List of samples to process
SAMPLES=(LP1 KMS11 MM1S RPMI8226 U266 PT1 PT2 PT3)

OUT_ROOT="${here}/6_Integration/6_bnd-ep-loop-gene"
mkdir -p "${OUT_ROOT}"

for sample in "${SAMPLES[@]}"; do
    echo "===== Annotating translocation-induced neo E-P loops: ${sample} ====="

    # enhancer=${here}/3_CUTtag/4_macs2/${sample}_*_H3K27ac*_peaks.narrowPeak
        # --- resolve enhancer safely (critical) ---
    enhancer=$(ls ${here}/3_CUTtag/4_macs2/${sample}*H3K27ac*peaks.narrowPeak 2>/dev/null | head -n 1)
    if [[ -z "${enhancer}" ]]; then
        echo "[WARN] enhancer peak file not found for ${sample}. Skip."
        continue
    fi
    gene_tpm=${here}/5_RNAseq/7_stringtie_tpm/${sample}/${sample}_coding_genes_tpm_fpkm.csv
    assembleBND=${here}/6_Integration/4_complex_bnd/${sample}/${sample}.assemblies.txt
    neo_loop=${here}/6_Integration/5_neoloop-caller/${sample}/${sample}.neo-loops.txt

    sample_out="${OUT_ROOT}/${sample}"
    mkdir -p "${sample_out}"
    cd "${sample_out}"

    ############################
    # 1) Translocation-induced neo-loops
    ############################
    awk -v OFS="\t" '{split($7,a,","); if(a[3]==1) print $0 }' \
        "${neo_loop}" >1_neo-loop.tsv

    ############################
    # 2) Split loop anchors
    ############################
    cat -n 1_neo-loop.tsv | awk -v OFS="\t" '{print $2,$3,$4,$1}' >2_anchor_left.tsv
    cat -n 1_neo-loop.tsv | awk -v OFS="\t" '{print $5,$6,$7,$1}' >2_anchor_right.tsv

    ############################
    # 3) Promoter–Enhancer (P–E)
    ############################
    bedtools intersect -a 2_anchor_left.tsv -b "${promoter}" -wa -wb \
        | awk -v OFS="\t" '{print $4,$1,$2,$3,$8}' | sort -u >tmp1
    bedtools intersect -a 2_anchor_right.tsv -b ${enhancer} -wa \
        | awk -v OFS="\t" '{print $4,$1,$2,$3}' | sort -u >tmp2

    sort -k1,1 -k2,2n tmp1 -o tmp1
    sort -k1,1 -k2,2n tmp2 -o tmp2
    
    join -t $'\t' tmp1 tmp2 \
        | awk -v OFS="\t" '{print $2,$3,$4,$6,$7,$8,$5,$1}' >3_neo-p-e-loop.tsv

    ############################
    # 4) Enhancer–Promoter (E–P)
    ############################
    bedtools intersect -a 2_anchor_left.tsv -b ${enhancer} -wa \
        | awk -v OFS="\t" '{print $4,$1,$2,$3}' | sort -u >tmp3
    bedtools intersect -a 2_anchor_right.tsv -b "${promoter}" -wa -wb \
        | awk -v OFS="\t" '{print $4,$1,$2,$3,$8}' | sort -u >tmp4

        sort -k1,1 -k2,2n tmp3 -o tmp3
        sort -k1,1 -k2,2n tmp4 -o tmp4

    join -t $'\t' tmp3 tmp4 \
        | awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$1}' >3_neo-e-p-loop.tsv

    ############################
    # 5) Combine P–E and E–P
    ############################
    cat 3_neo-p-e-loop.tsv 3_neo-e-p-loop.tsv | sort -u >4_neo-ep-loop.tsv

    cat 4_neo-ep-loop.tsv|cut -f 1-6|grep -w -f - 1_neo-loop.tsv >5_neo-ep-loop-bnd.tsv

    ############################
    # 6) Genes with TPM ≥ 10
    ############################
    grep -w -f <(cut -f7 4_neo-ep-loop.tsv) "${gene_tpm}" \
        | awk -v OFS="\t" '{if($9>=10) print}' >6_gene.tsv

    ############################
    # 7) Assemble BNDs supporting these loops
    ############################
    cat 6_gene.tsv|cut -f 2|grep -w -f - 4_neo-ep-loop.tsv |cut -f 1-6|grep -w -f - 5_neo-ep-loop-bnd.tsv >5_2_neo-ep-loop_plot.tsv

    cat 5_2_neo-ep-loop_plot.tsv|cut -f 7|grep -oE '\b[A-C][0-9]+\b' -|sort -u|grep -w -f - ${assembleBND} >7_assembleBND.txt

    awk 'NR==FNR {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; if (key in gene) gene[key]=gene[key]","$7; else gene[key]=$7; next} 
     {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; if (key in gene) print $0, gene[key],$3-$2; else print $0, ".",$3-$2}' \
    4_neo-ep-loop.tsv 5_2_neo-ep-loop_plot.tsv >8_bnd_neo-ep-loop-gene.tsv


    ############################
    # Cleanup
    ############################
    rm -f tmp1 tmp2 tmp3 tmp4

done
