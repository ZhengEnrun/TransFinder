#!/usr/bin/env python

import pandas as pd
import glob
import os

###############################################################################
# RNA-seq Step 3
#   For each sample:
#     1) Merge all StringTie *.tpm files (all replicates)
#     2) Compute mean Coverage, FPKM and TPM per gene
#     3) Save full gene table
#     4) Filter to coding genes using a coding-gene BED file (4th column = gene ID)
#
# Inputs:
#   7_stringtie_tpm/<sample>/*.tpm
#
# Outputs:
#   7_stringtie_tpm/<sample>/<sample>_average_tpm_fpkm_final.tsv
#   7_stringtie_tpm/<sample>/<sample>_coding_genes_tpm_fpkm.tsv
###############################################################################

# Project root (assumed to be current directory)
ROOT_DIR = os.path.abspath(".")

# All RNA-seq samples
samples = ['PT1', 'PT2', 'PT3', 'LP1', 'KMS11', 'U266', 'MM1S', 'RPMI8226']

# BED file containing coding genes (4th column must be gene ID)
coding_genes_bed = "/mnt/d/linux/reference/hg38/hg38p13/gencode.v40_coding_genes.bed"

for sample in samples:
    tpm_dir = os.path.join(ROOT_DIR, "7_stringtie_tpm", sample)
    file_pattern = os.path.join(tpm_dir, "*.tpm")

    file_paths = glob.glob(file_pattern)
    if not file_paths:
        print(f"[WARN] No .tpm files found for sample {sample}, skip.")
        continue

    print(f"[INFO] Sample {sample}: found {len(file_paths)} TPM files.")

    dfs = []
    for fp in file_paths:
        df = pd.read_csv(fp, sep="\t")
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)

    # Group by gene information, average Coverage/FPKM/TPM across replicates
    group_cols = ['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End']
    value_cols = ['Coverage', 'FPKM', 'TPM']

    result = (
        merged
        .groupby(group_cols, as_index=False)[value_cols]
        .mean()
    )

    # Save full gene table
    avg_out = os.path.join(tpm_dir, f"{sample}_average_tpm_fpkm_final.tsv")
    result.to_csv(avg_out, sep="\t", index=False)
    print(f"[INFO] Averaged TPM/FPKM written to {avg_out}")

    # Extract coding gene IDs from BED
    coding_ids_tmp = os.path.join(tpm_dir, f"{sample}_coding_ids.tmp")
    os.system(f"cut -f 4 {coding_genes_bed} | sort -u > {coding_ids_tmp}")

    # Filter averaged table by coding gene IDs (Gene ID column)
    coding_out = os.path.join(tpm_dir, f"{sample}_coding_genes_tpm_fpkm.tsv")
    cmd = f"grep -w -F -f {coding_ids_tmp} {avg_out} > {coding_out}"
    os.system(cmd)
    os.remove(coding_ids_tmp)

    print(f"[INFO] Coding genes TPM/FPKM written to {coding_out}")

print("All RNA-seq TPM/FPKM averaging and coding-gene filtering done.")
