#!/bin/bash
#可视化eaglec识别出来的sv
mkdir -p 1-egalec/3-plot
for i in RPMI8226 U266
do
  mkdir -p 1-egalec/3-plot/${i}/inter
  mkdir -p 1-egalec/3-plot/${i}/intra
  cp ../4-longDNA/3-11cell_combined_set/8-eaglec/${i}/${i}_translocation.tsv 1-egalec/3-plot/${i}/${i}_translocation.tsv
  cat 1-egalec/3-plot/${i}/${i}_translocation.tsv|cut -f 1,3|sort -u|awk -v OFS="\t" '{if($1==$2) {print $0}}' >1-egalec/3-plot/${i}/intra/bnd.txt
  cat 1-egalec/3-plot/${i}/${i}_translocation.tsv|cut -f 1,3|sort -u|awk -v OFS="\t" '{if($1!=$2) {print $0}}' >1-egalec/3-plot/${i}/inter/bnd.txt
  #translocation分为inter和intra
  for j in intra inter
  #循环每一行
  do
    cat 1-egalec/3-plot/${i}/${j}/bnd.txt |while read line
    do
      trans1=`echo $line | awk '{print $1}'`
      trans2=`echo $line | awk '{print $2}'`
      for k in 100000 250000 500000 1000000
      do
        mkdir -p 1-egalec/3-plot/${i}/${j}/${k}
        plot-interSVs --cool-uri 1-egalec/1-hic2mcool/${i}/${i}_mboi.mcool::resolutions/${k}  \
        --full-sv-file ../4-longDNA/3-11cell_combined_set/8-eaglec/${i}/${i}_eaglec.tsv \
        --output-figure-name 1-egalec/3-plot/${i}/${j}/${k}/${trans1}-${trans2}.png -C ${trans1} ${trans2} --balance-type ICE --dpi 800
      done
    done
  done
done




