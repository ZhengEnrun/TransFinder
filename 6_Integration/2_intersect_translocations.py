#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pipeline per sample (after Step1 normalization):

  1) Merge LRS translocations (ignore orientation), delta=500bp, dual-end joint clustering.
     - DOES NOT use strand for clustering.
     - Keeps Step1 order chrA <= chrB in output.
     Outputs into:
       6_Integration/2_intersection/<sample>/
         <sample>_longread_merged.tsv
         <sample>_longread_merge_map.tsv

  2) Intersect Hi-C vs merged LRS with D_BETWEEN=240kb.
     - No Hi-C expansion.
     - No Hi-C internal dedup/merge.
     Scheme A counts (ground truth = merged LRS events).
     Outputs into:
       6_Integration/2_intersection/<sample>/
         <sample>_intersection_hic.tsv
         <sample>_intersection_longread.tsv
         <sample>_orientation_conflicts.tsv
         <sample>_shared_components.tsv
         <sample>_exact_event_summary.txt
         <sample>_exact_shared_longread_breakpoints.bed

  3) Combined summary across all samples:
       6_Integration/2_intersection/all_samples_exact_event_summary.tsv

Inputs (per sample):
  6_Integration/1_trans_tsv/<sample>/<sample>_hic.tsv
  6_Integration/1_trans_tsv/<sample>/<sample>_longread.tsv

TSV columns (both):
  source  id  sample  chrA  posA  chrB  posB  strandA  strandB

Confirmed rules:
  - LRS merge:
      delta = 500 bp
      same ordered chr-pair (chrA, chrB) [Step1 already enforced]
      |posA_i - posA_j| <= 500 AND |posB_i - posB_j| <= 500
      strand ignored
      representative = closest to FLOAT medians
      tie-breaker = larger posB, then larger posA, then smaller id
  - Hi-C vs LRS overlap:
      tolerance = 240 kb
      direct dual-end matching (order already standardized by Step1)
      strand ignored for matching
  - Scheme A counts:
      Shared = number of merged LRS events supported by >=1 Hi-C call
"""

import os
import sys
from collections import defaultdict, deque
from statistics import median

# thresholds
DELTA_LR   = 500       # LRS internal merge threshold (bp), ignore orientation
D_BETWEEN  = 240000    # Hi-C vs LRS matching tolerance (bp)

SAMPLES = ["KMS11","LP1","MM1S","RPMI8226","U266","PT1","PT2","PT3"]


###############################################################################
# Utils
###############################################################################

def get_root_dir():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, ".."))

def load_tsv(path, expected_source=None):
    svs = []
    if not os.path.exists(path):
        return svs

    with open(path) as f:
        first = f.readline()
        if not first:
            return svs

        
        first_cols = first.rstrip("\n").split("\t")
        has_header = (len(first_cols) >= 1 and first_cols[0].lower() == "source")

        
        lines_iter = f
        if not has_header:
            lines_iter = [first] + list(f)

        for line in lines_iter:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue

            source, sv_id, sample, chrA, posA, chrB, posB, sA, sB = cols[:9]
            if expected_source and source != expected_source:
                continue
            try:
                posA = int(posA); posB = int(posB)
            except ValueError:
                continue

            svs.append({
                "source": source, "id": sv_id, "sample": sample,
                "chrA": chrA, "posA": posA,
                "chrB": chrB, "posB": posB,
                "strandA": sA, "strandB": sB
            })

    return svs


def group_by_chrpair_ordered(svs):
    # Step1 already guaranteed chrA<=chrB for both hic & lr
    d = defaultdict(list)
    for sv in svs:
        d[(sv["chrA"], sv["chrB"])].append(sv)
    return d


###############################################################################
# Step 1: Merge LRS (ignore orientation)
###############################################################################

def can_merge(repA, repB, sv, delta=DELTA_LR):
    return abs(sv["posA"] - repA) <= delta and abs(sv["posB"] - repB) <= delta

def pick_lr_representative(cluster):
    """
    Representative = closest to FLOAT medians.
    Tie-breaker: larger posB, then larger posA, then smaller id.
    Keeps Step1 order chrA<=chrB.
    """
    medA = median([x["posA"] for x in cluster])
    medB = median([x["posB"] for x in cluster])

    def key_fn(x):
        dist = abs(x["posA"] - medA) + abs(x["posB"] - medB)
        return (dist, -x["posB"], -x["posA"], x["id"])

    return min(cluster, key=key_fn)

def merge_lrs_ignore_orientation(lrs_svs, delta=DELTA_LR):
    """
    Merge LRS SVs into clusters per (chrA, chrB), ignoring strands.
    Returns list of clusters (each cluster is list of SV dicts).
    """
    clusters_all = []
    groups = group_by_chrpair_ordered(lrs_svs)

    for key, lst in groups.items():
        lst.sort(key=lambda x: (x["posA"], x["posB"]))
        cur = []
        repA = repB = None

        for sv in lst:
            if not cur:
                cur = [sv]
                repA, repB = sv["posA"], sv["posB"]
                continue

            if can_merge(repA, repB, sv, delta):
                cur.append(sv)
                repA = median([x["posA"] for x in cur])
                repB = median([x["posB"] for x in cur])
            else:
                clusters_all.append(cur)
                cur = [sv]
                repA, repB = sv["posA"], sv["posB"]

        if cur:
            clusters_all.append(cur)

    return clusters_all

def write_merged_lrs(sample, out_dir, lrs_svs):
    """
    Merge LRS (ignore orientation) and write outputs into out_dir:
      - <sample>_longread_merged.tsv
      - <sample>_longread_merge_map.tsv
    """
    clusters = merge_lrs_ignore_orientation(lrs_svs, DELTA_LR)

    os.makedirs(out_dir, exist_ok=True)

    merged_path = os.path.join(out_dir, f"{sample}_longread_merged.tsv")
    map_path    = os.path.join(out_dir, f"{sample}_longread_merge_map.tsv")

    with open(merged_path, "w") as out_m, open(map_path, "w") as out_map:
        out_m.write("\t".join([
            "source","id","sample","chrA","posA","chrB","posB","strandA","strandB"
        ]) + "\n")
        out_map.write("\t".join(["merged_id","original_id"]) + "\n")

        for i, cluster in enumerate(clusters, 1):
            rep = pick_lr_representative(cluster)
            mid = f"LRM_{i}"

            out_m.write("\t".join(map(str, [
                "longread", mid, rep["sample"],
                rep["chrA"], rep["posA"], rep["chrB"], rep["posB"],
                rep["strandA"], rep["strandB"]
            ])) + "\n")

            for sv in cluster:
                out_map.write("\t".join([mid, sv["id"]]) + "\n")

    return merged_path, clusters


###############################################################################
# Step 2: Hi-C vs merged LRS intersection (Scheme A)
###############################################################################

def write_transfinder_bnd(sample, lr_merged_svs, lr_ids_with_match, out_dir):
    """
    Generate bidirectional breakend file for downstream analysis.

    Output:
      <sample>_transfinder_bnd.tsv

    Each translocation is expanded into two records:
      chrA -> chrB
      chrB -> chrA

    Breakpoint orientation is defined by LRS.
    """
    out_path = os.path.join(out_dir, f"{sample}_transfinder_bnd.tsv")

    with open(out_path, "w") as out:
        for sv in lr_merged_svs:
            if sv["id"] not in lr_ids_with_match:
                continue

            chrA, posA, sA = sv["chrA"], sv["posA"], sv["strandA"]
            chrB, posB, sB = sv["chrB"], sv["posB"], sv["strandB"]

            # A -> B
            out.write("\t".join(map(str, [
                chrA, chrB, f"{sA}{sB}", posA, posB, "translocation"
            ])) + "\n")

            # B -> A
            out.write("\t".join(map(str, [
                chrB, chrA, f"{sB}{sA}", posB, posA, "translocation"
            ])) + "\n")

    return out_path


def intersect_hic_vs_lr(sample, hic_svs, lr_merged_svs, out_dir):
    hic_by_id = {sv["id"]: sv for sv in hic_svs}
    lr_by_id  = {sv["id"]: sv for sv in lr_merged_svs}

    hic_groups = group_by_chrpair_ordered(hic_svs)
    lr_groups  = group_by_chrpair_ordered(lr_merged_svs)

    matched_pairs = []
    hic_ids_with_match = set()
    lr_ids_with_match  = set()
    orientation_conflicts = []

    for key in hic_groups.keys() & lr_groups.keys():
        for h in hic_groups[key]:
            for l in lr_groups[key]:
                if (abs(h["posA"] - l["posA"]) <= D_BETWEEN and
                    abs(h["posB"] - l["posB"]) <= D_BETWEEN):

                    matched_pairs.append((h["id"], l["id"]))
                    hic_ids_with_match.add(h["id"])
                    lr_ids_with_match.add(l["id"])

                    if (h["strandA"], h["strandB"]) != (l["strandA"], l["strandB"]):
                        orientation_conflicts.append({
                            "sample": sample,
                            "chrA": h["chrA"], "posA": h["posA"],
                            "chrB": h["chrB"], "posB": h["posB"],
                            "hic_id": h["id"], "hic_strandA": h["strandA"], "hic_strandB": h["strandB"],
                            "lr_id": l["id"],  "lr_strandA": l["strandA"],  "lr_strandB": l["strandB"],
                        })

    # Build shared components (annotation only)
    adj = defaultdict(set)
    for h_id, l_id in matched_pairs:
        adj[("H", h_id)].add(("L", l_id))
        adj[("L", l_id)].add(("H", h_id))

    components = []
    visited = set()
    for node in adj.keys():
        if node in visited:
            continue
        q = deque([node])
        visited.add(node)
        comp_h, comp_l = set(), set()
        while q:
            side, cid = q.popleft()
            (comp_h if side == "H" else comp_l).add(cid)
            for nb in adj[(side, cid)]:
                if nb not in visited:
                    visited.add(nb)
                    q.append(nb)
        components.append({"hic": comp_h, "lr": comp_l})

    # ===== Scheme A counts =====
    all_hic_ids = set(hic_by_id.keys())
    all_lr_ids  = set(lr_by_id.keys())

    N_HiC_total = len(all_hic_ids)
    N_LR_total  = len(all_lr_ids)

    N_Shared = len(lr_ids_with_match)  # intersection on merged LRS events
    N_LongRead_only = N_LR_total - N_Shared
    N_HiC_only      = N_HiC_total - len(hic_ids_with_match)

    os.makedirs(out_dir, exist_ok=True)

    # 1) intersection hic
    with open(os.path.join(out_dir, f"{sample}_intersection_hic.tsv"), "w") as out:
        out.write("\t".join([
            "source","id","sample","chrA","posA","chrB","posB","strandA","strandB"
        ]) + "\n")
        for sv in hic_svs:
            if sv["id"] in hic_ids_with_match:
                out.write("\t".join(map(str, [
                    sv["source"], sv["id"], sv["sample"],
                    sv["chrA"], sv["posA"], sv["chrB"], sv["posB"],
                    sv["strandA"], sv["strandB"]
                ])) + "\n")

    # 2) intersection longread
    with open(os.path.join(out_dir, f"{sample}_intersection_longread.tsv"), "w") as out:
        out.write("\t".join([
            "source","id","sample","chrA","posA","chrB","posB","strandA","strandB"
        ]) + "\n")
        for sv in lr_merged_svs:
            if sv["id"] in lr_ids_with_match:
                out.write("\t".join(map(str, [
                    sv["source"], sv["id"], sv["sample"],
                    sv["chrA"], sv["posA"], sv["chrB"], sv["posB"],
                    sv["strandA"], sv["strandB"]
                ])) + "\n")

    # 3) orientation conflicts
    with open(os.path.join(out_dir, f"{sample}_orientation_conflicts.tsv"), "w") as out:
        out.write("\t".join([
            "sample","chrA","posA","chrB","posB",
            "hic_id","hic_strandA","hic_strandB",
            "lr_id","lr_strandA","lr_strandB"
        ]) + "\n")
        for r in orientation_conflicts:
            out.write("\t".join(map(str, [
                r["sample"], r["chrA"], r["posA"], r["chrB"], r["posB"],
                r["hic_id"], r["hic_strandA"], r["hic_strandB"],
                r["lr_id"],  r["lr_strandA"],  r["lr_strandB"]
            ])) + "\n")

    # 4) shared components (optional)
    with open(os.path.join(out_dir, f"{sample}_shared_components.tsv"), "w") as out:
        out.write("\t".join(["sample","component_id","n_hic","n_lr","hic_ids","lr_ids"]) + "\n")
        cid = 0
        for c in components:
            if not (c["hic"] and c["lr"]):
                continue
            out.write("\t".join(map(str, [
                sample, f"C{cid}", len(c["hic"]), len(c["lr"]),
                ",".join(sorted(c["hic"])), ",".join(sorted(c["lr"]))
            ])) + "\n")
            cid += 1

    # 5) summary
    summary_path = os.path.join(out_dir, f"{sample}_exact_event_summary.txt")
    with open(summary_path, "w") as out:
        out.write("\t".join([
            "sample","DELTA_LR","D_BETWEEN",
            "N_HiC_total","N_LR_total",
            "N_HiC_only","N_LongRead_only","N_Shared"
        ]) + "\n")
        out.write("\t".join(map(str, [
            sample, DELTA_LR, D_BETWEEN,
            N_HiC_total, N_LR_total,
            N_HiC_only, N_LongRead_only, N_Shared
        ])) + "\n")

    # 6) shared LR breakpoints bed per shared component
    bed_path = os.path.join(out_dir, f"{sample}_exact_shared_longread_breakpoints.bed")
    with open(bed_path, "w") as out:
        eid = 0
        for c in components:
            if not (c["hic"] and c["lr"]):
                continue
            for lr_id in c["lr"]:
                sv = lr_by_id[lr_id]
                out.write("\t".join(map(str, [
                    sv["chrA"], sv["posA"], sv["chrB"], sv["posB"],
                    sv["id"], f"EVENT_{eid}"
                ])) + "\n")
            eid += 1

    # 7) TransFinder bidirectional BND for downstream analysis
    write_transfinder_bnd(
        sample,
        lr_merged_svs,
        lr_ids_with_match,
        out_dir
    )

    return {
        "sample": sample,
        "N_HiC_total": N_HiC_total,
        "N_LR_total": N_LR_total,
        "N_HiC_only": N_HiC_only,
        "N_LR_only": N_LongRead_only,
        "N_Shared": N_Shared
    }


###############################################################################
# Main
###############################################################################

def main():
    root_dir = get_root_dir()
    all_summary = []

    sys.stderr.write(
        f"[INFO] LRS merge(ignore orientation) DELTA={DELTA_LR} bp; "
        f"HiC-LRS match D_BETWEEN={D_BETWEEN} bp; Scheme A counts.\n"
    )

    for sample in SAMPLES:
        in_dir = os.path.join(root_dir, "6_Integration", "1_trans_tsv", sample)
        out_dir = os.path.join(root_dir, "6_Integration", "2_intersection", sample)

        hic_path = os.path.join(in_dir, f"{sample}_hic.tsv")
        lr_path  = os.path.join(in_dir, f"{sample}_longread.tsv")

        hic_svs = load_tsv(hic_path, expected_source="hic")
        lrs_svs = load_tsv(lr_path,  expected_source="longread")

        if not hic_svs or not lrs_svs:
            sys.stderr.write(f"[WARN] {sample}: missing hic or longread TSV, skip.\n")
            continue

        # Step 1 merge LRS -> write merged to out_dir
        merged_lr_path, _clusters = write_merged_lrs(sample, out_dir, lrs_svs)
        lr_merged_svs = load_tsv(merged_lr_path, expected_source="longread")

        # Step 2 intersect
        summ = intersect_hic_vs_lr(sample, hic_svs, lr_merged_svs, out_dir)
        all_summary.append(summ)

    # Combined summary across samples
    comb_path = os.path.join(root_dir, "6_Integration", "2_intersection",
                             "all_samples_exact_event_summary.tsv")
    os.makedirs(os.path.dirname(comb_path), exist_ok=True)
    with open(comb_path, "w") as out:
        out.write("\t".join([
            "sample","N_HiC_total","N_LR_total","N_HiC_only","N_LR_only","N_Shared"
        ]) + "\n")
        for r in all_summary:
            out.write("\t".join(map(str, [
                r["sample"], r["N_HiC_total"], r["N_LR_total"],
                r["N_HiC_only"], r["N_LR_only"], r["N_Shared"]
            ])) + "\n")

    sys.stderr.write(f"[INFO] Wrote combined summary: {comb_path}\n")


if __name__ == "__main__":
    main()
