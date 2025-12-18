# TransFinder

**TransFinder** is a multi-omics analytical pipeline for identifying high-confidence chromosomal translocations and uncovering their 3D regulatory consequences, with a focus on translocation-induced enhancer–promoter (E–P) neo-loops and candidate oncogenes in cancer.

By integrating long-read sequencing (SMRT-seq), Hi-C, CUT&Tag/ChIP-seq, and RNA-seq, TransFinder enables accurate translocation detection and systematic reconstruction of post-translocation chromatin architecture.

---

## Key features

- **High-confidence translocation detection**  
  Translocations are retained only when supported by both SMRT-seq and Hi-C, achieving ~98% experimental validation.

- **Precise breakpoint and orientation resolution**  
  Long-read sequencing enables accurate definition of translocation breakpoints and fragment orientation, which is essential for reconstructing rearranged chromatin topology.

- **Direct linkage to regulatory rewiring**  
  CNV-aware Hi-C reconstruction combined with epigenomic and transcriptomic annotation allows direct identification of translocation-induced neo-loops and candidate oncogenes.

---

## Overview of the workflow

1. Structural variants are detected independently from Hi-C contact matrices and SMRT-seq data.  
2. High-confidence chromosomal translocations are defined by cross-validation between the two platforms.  
3. Breakpoints and fragment orientation are refined using long-read sequencing.  
4. CNV-normalized Hi-C matrices are locally reconstructed to predict translocation-induced neo-loops.  
5. Enhancer–promoter neo-loops and candidate oncogenes are annotated using CUT&Tag/ChIP-seq and RNA-seq data.  
6. Translocation-induced neo-loops and associated regulatory features are visualized using NeoLoopFinder, integrating Hi-C contact maps with H3K27ac CUT&Tag/ChIP-seq, RNA-seq, and long-read sequencing signals.

---

## Software requirements

### Long-read sequencing and structural variant analysis
- pbmm2 v1.13.1  
- pbsv v2.9.0  

### Hi-C / Micro-C analysis
- Trimmomatic v0.39  
- HiC-Pro v3.1.0  
- cooler v0.9.3  
- hic2cool v0.8.3  
- pairtools v1.0.3  
- BWA-MEM v0.7.17  
- EagleC v0.1.8  

### CUT&Tag / ChIP-seq
- fastp v0.23.2  
- Bowtie2 v2.2.5  
- Picard v2.27.5  
- MACS2 v2.2.7.1  

### RNA-seq
- fastp v0.23.2  
- HISAT2 v2.2.1  
- StringTie v2.2.3  

### Integrated neo-loop analysis and visualization
- NeoLoopFinder v0.4.3  
- bedtools v2.30.0  

### Programming environment
- Python ≥ 3.7  
- R ≥ 4.0  
