#!/bin/bash

#assemble-complexSVs 

here=`pwd`
export NUMEXPR_MAX_THREADS=30
for sample in LP1 
do

  for r in 120000
  do
    high_confidence_bnd=${here}/5_compareSV/2_vs/${sample}/CNV/${r}/3_2methods_sv_longread_bnd.bed
    mkdir -p 6_assemble-complexSVs/longread/${sample}/${r}
    cd 6_assemble-complexSVs/longread/${sample}/${r}

    ribbon_bnd_bedpe=/mnt/f/zer/hg38_chr1-22x/1-pbsv/10_vcf2bed_ribbon/${sample}/${sample}_bnd_whith_strand.bedpe
    #取三代和hic都找到的易位，用三代的断点,三代的方向
    cat ${high_confidence_bnd}|sed 's/\t/:/1'|sed 's/\t/-/1'|sed 's/\t/:/1'|cut -f 1|grep -f - ${ribbon_bnd_bedpe} >2methods_bnd_longread.bed.tmp
    cat 2methods_bnd_longread.bed.tmp|awk -v FS="\t" -v OFS="\t" '{print $1,$4,$9$10,$2,$5,"translocation"}' >${sample}_2methods_bnd_longread.bed.tmp1
    #取三代和hic都找到的易位，用三代的断点,hic的方向
    cat ${high_confidence_bnd}|awk -v FS="\t" -v OFS="\t" '{print $1,$3,$5,$2,$4,"translocation"}' >${sample}_2methods_bnd_longread.bed.tmp2

    cat ${sample}_2methods_bnd_longread.bed.tmp1 ${sample}_2methods_bnd_longread.bed.tmp2|sort -k1,1 -k2,2 |sort -u >${sample}_2methods_bnd_longread_both_direction.bed
    rm 2methods_bnd_longread.bed.tmp ${sample}_2methods_bnd_longread.bed.tmp1 ${sample}_2methods_bnd_longread.bed.tmp2

    python ${here}/6.2_change_bnd_position.py -i ${sample}_2methods_bnd_longread_both_direction.bed -o ${sample}_2methods_bnd_longread.bed


    assemble-complexSVs -O ${sample} \
                        -B ${sample}_2methods_bnd_longread.bed \
                        --balance-type CNV --protocol insitu --nproc 30 \
                        -H ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/25000 \
                        ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/10000 \
                        ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/5000

    cd ${here}
  done
done

