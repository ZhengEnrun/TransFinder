#!/bin/bash

#call cnv and normalize hic maps
export NUMEXPR_MAX_THREADS=30

here=`pwd`

for sample in LP1 
do

  for r in 5000 10000 25000 50000
  do
    mkdir -p 3_correct-cnv/${sample}/${r}

    echo ${r} "calculate-cnv"
    calculate-cnv -H ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/${r} -g hg38 \
                  -e MboI --output ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-profile.bedGraph 

    echo ${r} "segment-cnv"
    segment-cnv --cnv-file ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-profile.bedGraph --binsize ${r} \
                --ploidy 3 --output ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-seg.bedGraph --nproc 30 
    
    echo ${r} "plot-cnv"
    plot-cnv --cnv-profile ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-profile.bedGraph \
            --cnv-segment ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-seg.bedGraph \
            --output-figure-name ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV.genome-wide.png \
            --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 \
            --label-size 7 --tick-label-size 6 --clean-mode 

    echo ${r} "correct-cnv"
    correct-cnv -H ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/${r} \
            --cnv-file ${here}/3_correct-cnv/${sample}/${r}/${sample}_${r}.CNV-seg.bedGraph --nproc 30 -f 

  done
done
