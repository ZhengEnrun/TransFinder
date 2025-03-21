from neoloop.visualize.core import *
import cooler


clr = cooler.Cooler('PT3_contact.mcool::resolutions/5000')
neoloop = "5_3_diffgene_neo-ep-SV.tsv"

cancer_atac ='PT3_ATACseq.bw'
cancer_H3K27ac_chip = 'PT3_H3K27ac_CUTTag.bw'
cancer_RNAseq = 'PT3_RNAseq.bw'
cancer_longread = 'PT3_SMRTseq.bw'

assembly = 'C0	translocation,1,178385884,-,22,39280931,+	1,179975000	22,38925000'


ac = assembly.split('\t')[0]

print(assembly)


vis = Triangle(clr, assembly, n_rows=7, figsize=(8, 6),
                track_partition=[8, 0.4, 0.8, 0.8, 0.8, 0.8, 0.5], correct='sweight', span=1000000,
                slopes={(0, 0): 1, (0, 1): 0.3, (1, 1): 1})

vis.matrix_plot(vmin=0,vmax=0.01)
vis.plot_chromosome_bounds(linewidth=2)
if open(neoloop, 'r').read().find(ac) > 0:
    vis.plot_loops(neoloop, face_color='none', marker_size=50, cluster=False, filter_by_res=True, onlyneo=True)

    vis.plot_genes(release=106, filter_=['CBX7','APOBEC3C'], label_aligns={'CBX7': 'right'}, fontsize=8)

    
    vis.plot_signal('ATAC', cancer_atac, label_size=8, data_range_size=9,max_value=100,color='#d77800')
    vis.plot_signal('H3K27ac', cancer_H3K27ac_chip, label_size=8, data_range_size=9, max_value=50, color='#6A3D9A')
    vis.plot_signal('RNA-seq', cancer_RNAseq, label_size=8, data_range_size=9, max_value=10, color='#E31A1C')
    vis.plot_signal('SMRTseq', cancer_longread, label_size=8, data_range_size=9, max_value=5, color='#808080')

    vis.plot_chromosome_bar(name_size=10, coord_size=8)

    vis.outfig("PT3_chr1chr22_CBX7_APOBEC3C.pdf")

