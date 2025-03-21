#!/bin/bash
here=`pwd`
#convert hic data format
for sample in LP1 
do
    mkdir -p 2_get_hic_mcool/${sample}
    cd 2_get_hic_mcool/${sample}
    
    /mnt/d/linux/software/HiC-Pro-master/bin/utils/hicpro2juicebox.sh \
    -i ${here}/1_hicpro/${sample}/allvalidpairs_results/${sample}.allValidPairs \
	-g /mnt/d/linux/reference/hg38/hg38p13r/hg38_chr1_22xym.chrom.sizes.txt \
	-j /mnt/d/linux/software/HiC-Pro-master/juicer_tools_1.22.01.jar \
	-r /mnt/d/linux/reference/hg38/hg38p13r/hic/hg38_mboi.bed \
	-t ./tmp \
	-o ./ 


    hic2cool convert ${sample}.allValidPairs.hic ${sample}_contact.mcool -r 0 -p 30

    echo "hic2mcool done!!"

    for r in 5000 10000 25000 50000
    do
        echo "cooler_balance" ${sample} ${r}
        cooler balance ${sample}_contact.mcool::/resolutions/${r} -p 30
    done
    cd ${here}
done

