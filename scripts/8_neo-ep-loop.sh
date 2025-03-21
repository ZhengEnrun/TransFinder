#!/bin/bash
#annotate translocation-induced e-p neo-loops using rna and chip-seq/cut-tag for H3K27ac data

here=`pwd`
promoter=/mnt/d/linux/reference/hg38/hg38p13/coding_gene_promoters_1kb.bed

for sample in LP1 
do
  for bp in 120000
  do
    enhancer=/mnt/f/zer/hg38_chr1-22x/7_chipseq/${sample}/4_macs2/${sample}_ChIP_H3K27ac_peaks.narrowPeak
    gene_tpm=/mnt/f/zer/hg38_chr1-22x/0_RNA/7_stringtie_tpm/${sample}/${sample}_coding_genes_tpm_fpkm.csv
    assembleSV=${here}/6_assemble-complexSVs/longread/${sample}/${bp}/${sample}.assemblies.txt

    mkdir -p 8_neo-ep-loop/${sample}/${bp}
    cd 8_neo-ep-loop/${sample}/${bp}

    cat ${here}/7_neoloop-caller/${sample}/${bp}/${sample}.neo-loops.txt|awk -v OFS="\t" '{split($7,a,",");if(a[3]==1) print $0 }' >1_neo-loop.tsv
    cat -n 1_neo-loop.tsv |awk -v OFS="\t" '{print $2,$3,$4,$1}' >2_neo-loop_anchor_left.tsv
    cat -n 1_neo-loop.tsv |awk -v OFS="\t" '{print $5,$6,$7,$1}' >2_neo-loop_anchor_right.tsv
    #将neo-loop anchor与E-P loop做overlap
    #左边promoter，右边enhancer
    bedtools intersect -a 2_neo-loop_anchor_left.tsv -b ${promoter} -wa -wb|awk -v OFS="\t" '{print $4,$1,$2,$3,$8}'|sort -u|sort -k1,1 >tmp1
    bedtools intersect -a 2_neo-loop_anchor_right.tsv -b ${enhancer} -wa|awk -v OFS="\t" '{print $4,$1,$2,$3}'|sort -u|sort -k1,1 >tmp2
    join -t $'\t' tmp1 tmp2|awk -v OFS="\t" '{print $2,$3,$4,$6,$7,$8,$5,$1}' >3_neo-p-e-loop.tsv

    #左边enhancer，右边promoter
    bedtools intersect -a 2_neo-loop_anchor_left.tsv -b ${enhancer} -wa|awk -v OFS="\t" '{print $4,$1,$2,$3}'|sort -u|sort -k1,1 >tmp3
    bedtools intersect -a 2_neo-loop_anchor_right.tsv -b ${promoter} -wa -wb|awk -v OFS="\t" '{print $4,$1,$2,$3,$8}' |sort -u|sort -k1,1 >tmp4
    join -t $'\t' tmp3 tmp4|awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$1}' >3_neo-e-p-loop.tsv
    #将ep和pe的结果结合在一起
    cat 3_neo-p-e-loop.tsv 3_neo-e-p-loop.tsv |sort -u >4_neo-ep-loop.tsv
    #loop文件
    cat 4_neo-ep-loop.tsv|cut -f 1-6|grep -w -f - 1_neo-loop.tsv >5_neo-ep-SV.tsv
    #这些ep对应的基因tpm>10
    cat 4_neo-ep-loop.tsv|cut -f 7|grep -w -f - ${gene_tpm}|awk -v OFS="\t" '{if($9>=10) print $0}' >6_neo-ep-diffgene.tsv
    #差异基因的loop所对应的assemble sv
    cat 6_neo-ep-diffgene.tsv|cut -f 2|grep -w -f - 4_neo-ep-loop.tsv |cut -f 1-6|grep -w -f - 5_neo-ep-SV.tsv >5_2_diffgene_neo-ep-SV.tsv

    cat 5_2_diffgene_neo-ep-SV.tsv|cut -f 7|grep -oE '\b[A-C][0-9]+\b' -|sort -u|grep -w -f - ${assembleSV} >7_assembleSV.txt

    awk 'NR==FNR {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; if (key in gene) gene[key]=gene[key]","$7; else gene[key]=$7; next} 
     {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; if (key in gene) print $0, gene[key],$3-$2; else print $0, ".",$3-$2}' \
    4_neo-ep-loop.tsv 5_2_diffgene_neo-ep-SV.tsv >5_3_diffgene_neo-ep-SV.tsv

    cat 6_neo-ep-diffgene.tsv|cut -f 2 >gene_list_promoters_1kb
    rm tmp1 tmp2 tmp3 tmp4
    cd ${here}
  done
done



