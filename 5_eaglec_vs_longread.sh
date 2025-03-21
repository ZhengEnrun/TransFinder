#!/bin/bash
#compare translocation between hic and long-read sequencing data

here=`pwd`

chr_size=/mnt/d/linux/reference/hg38/hg38p13r/hg38_chr1_22xym.chrom.sizes.txt

for sample in LP1 
do

  mkdir -p ${here}/5_compareSV/1_bed/${sample}

  for method in CNV
  do

    cd ${here}/5_compareSV/1_bed/${sample}

    cat ${here}/4_predictSV/NeoLoopFinder/${sample}/${sample}.predictsv.txt.CNN_SVs.5K_combined.txt| \
    awk -v OFS="\t" '{if($1!=$2) print$1,$4,$2,$5,$3}'|sort -k1,1 -k2,2n >${sample}_${method}_translocation.bed.tmp


    cat ${sample}_${method}_translocation.bed.tmp |awk -v OFS="\t" '{split($1,a,"hr");split($3,b,"hr");if(a[2]<=b[2]) {print $1,$2,$3,$4,$5} else{print $3,$4,$1,$2,$5}}' >${sample}_${method}_translocation.bed
    rm ${sample}_${method}_translocation.bed.tmp

    cat -n ${sample}_${method}_translocation.bed|awk -v OFS="\t" -v cell=${sample} '{print $2,$3,$3+1,$1cell,$6}' >${sample}_${method}_translocation_left.bed
    cat -n ${sample}_${method}_translocation.bed|awk -v OFS="\t" -v cell=${sample} '{print $4,$5,$5+1,$1cell,$6}' >${sample}_${method}_translocation_right.bed


    cat /mnt/f/zer/hg38_chr1-22x/1-pbsv/8_sv_range_bed/${sample}/${sample}_bnd.bed|awk -v OFS="\t" '{if($1!=$3) print $0}' >${sample}_longread_translocation.bed

    cat -n ${sample}_longread_translocation.bed|awk -v OFS="\t" -v cell=${sample} '{print $2,$3,$3+1,$1cell}' >${sample}_longread_translocation_left.bed
    cat -n ${sample}_longread_translocation.bed|awk -v OFS="\t" -v cell=${sample} '{print $4,$5,$5+1,$1cell}' >${sample}_longread_translocation_right.bed


    for bp in 120000  
    do
      mkdir -p ${here}/5_compareSV/2_vs/${sample}/${method}/${bp}
      cd ${here}/5_compareSV/2_vs/${sample}/${method}/${bp}

        for k in left right
        do

          bedtools slop -i ${here}/5_compareSV/1_bed/${sample}/${sample}_${method}_translocation_${k}.bed -g ${chr_size} \
          -b ${bp} |awk -v OFS="\t" '{print $1,$2,$3-1,$4,$5}' > 1_${sample}_${method}_translocation_${k}_slop${bp}.bed

          bedtools slop -i ${here}/5_compareSV/1_bed/${sample}/${sample}_longread_translocation_${k}.bed -g ${chr_size} \
          -b ${bp} |awk -v OFS="\t" '{print $1,$2,$3-1,$4}' > 1_${sample}_longread_translocation_${k}_slop${bp}.bed

          bedtools intersect -a 1_${sample}_longread_translocation_${k}_slop${bp}.bed -b 1_${sample}_${method}_translocation_${k}_slop${bp}.bed -wa -wb \
          |sort -u |awk -v OFS="\t" '{print $4,$1,$2,$3,$5,$6,$7,$8,$9}'|sort > 2_${method}_longread_${k}_combined.set.tmp
        done


      join -t $'\t' 2_${method}_longread_left_combined.set.tmp 2_${method}_longread_right_combined.set.tmp | \
      awk -v OFS="\t" '{if($8==$16) print $0}'|sort -u >2_${method}_eaglec_longread_combined.set.join

      cat 2_${method}_eaglec_longread_combined.set.join|awk -v OFS="\t" -v r="${bp}" '{print $5,$6+r,$13,$14+r,$9}'|sort -u >3_2methods_eaglec_bnd.bed


      cat 2_${method}_eaglec_longread_combined.set.join|awk -v OFS="\t" -v r="${bp}" '{print $2,$3+r,$10,$11+r,$9}'|sort -u >3_2methods_sv_longread_bnd.bed


      cd ${here}
    done
  done
done



