#!/bin/bash
for sample in LP1
do
    /mnt/d/linux/software/HiC-Pro-master/bin/HiC-Pro -i 0_rawdata/${sample} \
    -o 1_hicpro/${sample} -c /mnt/d/linux/reference/hg38/hg38p13r/hic/config-hicpro.txt 
done

