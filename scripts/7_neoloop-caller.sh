#!/bin/bash
#call translocation-induced-loop
export NUMEXPR_MAX_THREADS=30
here=`pwd`
for sample in LP1 
  for r in 120000
  do
    mkdir -p 7_neoloop-caller/${sample}/${r}
    cd 7_neoloop-caller
    neoloop-caller -O ${sample}/${r}/${sample}.neo-loops.txt \
                  --assembly ${here}/6_assemble-complexSVs/longread/${sample}/${r}/${sample}.assemblies.txt \
                  --balance-type CNV --protocol insitu --prob 0.95 --nproc 30 \
                    -H ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/25000 \
                    ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/10000 \
                    ${here}/2_get_hic_mcool/${sample}/${sample}_contact.mcool::resolutions/5000

    cd ${here}
  done
done
