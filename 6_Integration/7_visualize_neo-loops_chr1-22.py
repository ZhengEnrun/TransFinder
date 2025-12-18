from neoloop.visualize.core import *
import cooler
import os
import pandas as pd


cells = ['PT3']

resolutions = ['5000','10000','25000']

for cell in cells:
    for r in resolutions:
        os.system("mkdir -p 7_visualized-neo-ep-loop" )
        clr = cooler.Cooler('/mnt/f/zer/TransFinder/2_HiC/2_get_hic_mcool/%s/%s_contact.mcool::resolutions/%s' % (cell, cell,r))
        neoloop = "./6_bnd-ep-loop-gene/%s/5_2_neo-ep-loop_plot.tsv"% (cell)
        
        cancer_atac ='/mnt/f/zer/TransFinder/4_ATAC/3_bw/%s.bw'%(cell)
        cancer_H3K27ac_chip = "/mnt/f/zer/TransFinder/3_CUTtag/3_bw/%s.bw"% (cell)
        cancer_RNAseq = "/mnt/f/zer/TransFinder/5_RNAseq/3_bamCoverage/%s/R1/%s_R1.bw"% (cell,cell)
        cancer_longread = "/mnt/f/zer/TransFinder/1_SMRT-seq/7-bam2bw/%s.bw"% (cell)
        
        genelist = []


        #pt3
        assembly = 'C0	translocation,1,178385884,-,22,39280931,+	1,179975000	22,38925000'
        ac = assembly.split('\t')[0]

        print(assembly)


        vis = Triangle(clr, assembly, n_rows=8, figsize=(8, 6),
                        track_partition=[8, 0.4, 0.4, 0.8, 0.8, 0.8, 0.8, 0.5], correct='sweight', span=1000000,
                        slopes={(0, 0): 1, (0, 1): 0.3, (1, 1): 1})
   
        vis.matrix_plot(vmin=0,vmax=0.01)
        vis.plot_chromosome_bounds(linewidth=2)
        if open(neoloop, 'r').read().find(ac) > 0:
            vis.plot_loops(neoloop, face_color='none', marker_size=50, cluster=False, filter_by_res=True, onlyneo=True)

            vis.plot_genes(release=106, filter_=['CBX7','APOBEC3C'], label_aligns={'CBX7': 'right'}, fontsize=8)
            
            vis.plot_signal('ATAC', cancer_atac, label_size=8, data_range_size=9,max_value=100,color='#d77800')
            vis.plot_signal('H3K27ac', cancer_H3K27ac_chip, label_size=8, data_range_size=9, max_value=50, color='#6A3D9A')
            vis.plot_signal('RNA-seq', cancer_RNAseq, label_size=8, data_range_size=9, max_value=10, color='#E31A1C')
            vis.plot_signal('SMRT-seq', cancer_longread, label_size=8, data_range_size=9, max_value=5, color='#808080')

            vis.plot_chromosome_bar(name_size=10, coord_size=8)

            vis.outfig("./7_visualized-neo-ep-loop/%s_%s_%s.pdf"% (cell, ac,r), dpi=300)

