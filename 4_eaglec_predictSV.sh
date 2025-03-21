#!/bin/bash

#predict SVs from hic data using EagleC
export NUMEXPR_MAX_THREADS=30

here=`pwd`
for sample in LP1
do
  for f in NeoLoopFinder 
  do
    mkdir -p 4_predictSV/${sample}/${f}
    cd 4_predictSV/${sample}
    
    predictSV --hic-5k ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/5000 \
                    --hic-10k ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/10000 \
                    --hic-50k ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/50000 \
                    -O ${f}/${sample}.predictsv -g hg38 --balance-type CNV --output-format ${f} \
                    --prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999 
    cd ${here}
  done
done
