#!/usr/bin/env Rscript

setwd("F:/zer/TransFinder/6_Integration")
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# CNS-style palette
CNS_COLORS <- c(
  "Orientation matched"    = "#E64B35",
  "Orientation mismatched" = "#4DBBD5",
  "Matched"                = "#E64B35",
  "Mismatched"             = "#4DBBD5",
  "LRS_only"               = "#00A087",
  "Shared"                 = "#3C5488",
  "HiC_only"               = "#F39B7F"
)

# =========================
# Config
# =========================
samples <- c("KMS11","LP1","MM1S","RPMI8226","U266","PT1","PT2","PT3")

# 
sample_order <- c("KMS11","LP1","MM1S","RPMI8226","U266","PT1","PT2","PT3")

wd <- getwd()
root_dir <- if (basename(wd) == "6_Integration") dirname(wd) else wd

dir_intersection <- file.path(root_dir, "6_Integration", "2_intersection")
dir_trans_tsv    <- file.path(root_dir, "6_Integration", "1_trans_tsv")
dir_plot         <- file.path(root_dir, "6_Integration", "3_plot")
if (!dir.exists(dir_plot)) dir.create(dir_plot, recursive = TRUE)

DELTA_LR  <- 500
D_BETWEEN <- 240000

# =========================
# Blacklist (remove PT3 event globally)
# =========================
EXCLUDE <- list(
  PT3 = list(
    lr_ids  = c("LRM_2"),
    hic_ids = c("HIC_2")
  )
)

# =========================
# Helper: safe fread header/no-header
# =========================
fread_maybe_noheader <- function(path, colnames_expected) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path, header = TRUE)
  if (!(names(dt)[1] %in% colnames_expected)) {
    dt <- fread(path, header = FALSE)
    setnames(dt, colnames_expected)
  }
  return(dt)
}

apply_blacklist <- function(dt, s, id_col="id", type=c("lr","hic")) {
  type <- match.arg(type)
  if (s %in% names(EXCLUDE)) {
    bl <- EXCLUDE[[s]]
    if (type=="lr" && !is.null(bl$lr_ids)) {
      dt <- dt[!(get(id_col) %in% bl$lr_ids)]
    }
    if (type=="hic" && !is.null(bl$hic_ids)) {
      dt <- dt[!(get(id_col) %in% bl$hic_ids)]
    }
  }
  return(dt)
}

# =========================
# Recompute per-sample counts AFTER blacklist
# (Scheme A: truth = merged LRS events)
# Also keep matched HIC-LR pairs for orientation analysis
# Also write per-event match detail (ALL merged LRS events)
# =========================
count_list <- list()
shared_pairs_all <- list()

hic_cols <- c("source","id","sample","chrA","posA","chrB","posB","strandA","strandB")

for (s in samples) {
  
  hic_path  <- file.path(dir_trans_tsv, s, sprintf("%s_hic.tsv", s))
  lr_m_path <- file.path(dir_intersection, s, sprintf("%s_longread_merged.tsv", s))
  
  hic_dt <- fread_maybe_noheader(hic_path, hic_cols)
  lr_dt  <- fread_maybe_noheader(lr_m_path, hic_cols)
  
  if (is.null(lr_dt) || nrow(lr_dt)==0) next
  if (is.null(hic_dt)) hic_dt <- data.table()
  
  # NOTE: do NOT apply blacklist here.
  # lr_dt  <- apply_blacklist(lr_dt, s, type="lr")
  # hic_dt <- apply_blacklist(hic_dt, s, type="hic")
  
  lr_dt[, posA := as.integer(posA)]
  lr_dt[, posB := as.integer(posB)]
  if (nrow(hic_dt)>0) {
    hic_dt[, posA := as.integer(posA)]
    hic_dt[, posB := as.integer(posB)]
  }
  
  N_LR_total  <- nrow(lr_dt)
  N_HiC_total <- nrow(hic_dt)
  
  per_event_dt <- NULL
  
  if (N_HiC_total == 0) {
    N_Shared   <- 0
    N_LR_only  <- N_LR_total
    N_HiC_only <- 0
    
    per_event_dt <- lr_dt[, .(
      sample = s,
      lr_id = id,
      chrA, posA, strandA,
      chrB, posB, strandB,
      is_shared = FALSE
    )]
    
    fwrite(
      per_event_dt,
      file.path(dir_plot, sprintf("%s_LRS_match_detail.tsv", s)),
      sep="\t"
    )
    
    count_list[[length(count_list)+1]] <- data.table(
      sample=s,
      N_HiC_total=N_HiC_total,
      N_LR_total=N_LR_total,
      N_HiC_only=N_HiC_only,
      N_LR_only=N_LR_only,
      N_Shared=N_Shared
    )
    next
  }
  
  # match within same chr-pair + dual-end 240kb tolerance
  hic_keys <- unique(hic_dt[, .(chrA, chrB)])
  lr_keys  <- unique(lr_dt[,  .(chrA, chrB)])
  keys <- merge(hic_keys, lr_keys, by=c("chrA","chrB"))
  
  hic_with_match <- character()
  lr_with_match  <- character()
  mp_list <- list()
  
  if (nrow(keys)>0) {
    for (i in 1:nrow(keys)) {
      ka <- keys$chrA[i]; kb <- keys$chrB[i]
      hic_sub <- hic_dt[chrA==ka & chrB==kb]
      lr_sub  <- lr_dt[ chrA==ka & chrB==kb]
      if (nrow(hic_sub)==0 || nrow(lr_sub)==0) next
      
      for (j in 1:nrow(lr_sub)) {
        l <- lr_sub[j]
        m <- hic_sub[
          abs(posA - l$posA) <= D_BETWEEN &
            abs(posB - l$posB) <= D_BETWEEN
        ]
        if (nrow(m)==0) next
        
        lr_with_match  <- c(lr_with_match, l$id)
        hic_with_match <- c(hic_with_match, m$id)
        
        # ===== store PAIR-LEVEL info with breakpoints on BOTH sides =====
        mp_list[[length(mp_list)+1]] <- data.table(
          sample = s,
          
          # LRS side (one event)
          lr_id = l$id,
          lr_chrA = l$chrA, lr_posA = l$posA, lr_strandA = l$strandA,
          lr_chrB = l$chrB, lr_posB = l$posB, lr_strandB = l$strandB,
          
          # Hi-C side (can be multiple rows)
          hic_id = m$id,
          hic_chrA = m$chrA, hic_posA = m$posA, hic_strandA = m$strandA,
          hic_chrB = m$chrB, hic_posB = m$posB, hic_strandB = m$strandB
        )
      }
    }
  }
  
  lr_with_match  <- unique(lr_with_match)
  hic_with_match <- unique(hic_with_match)
  
  # ---------- per-event match detail (ALL merged LRS events) ----------
  per_event_dt <- lr_dt[, .(
    sample = s,
    lr_id = id,
    chrA, posA, strandA,
    chrB, posB, strandB,
    is_shared = id %in% lr_with_match
  )]
  
  fwrite(
    per_event_dt,
    file.path(dir_plot, sprintf("%s_LRS_match_detail.tsv", s)),
    sep="\t"
  )
  
  N_Shared   <- length(lr_with_match)                 # Scheme A
  N_LR_only  <- N_LR_total - N_Shared
  N_HiC_only <- N_HiC_total - length(hic_with_match)
  
  count_list[[length(count_list)+1]] <- data.table(
    sample=s,
    N_HiC_total=N_HiC_total,
    N_LR_total=N_LR_total,
    N_HiC_only=N_HiC_only,
    N_LR_only=N_LR_only,
    N_Shared=N_Shared
  )
  
  if (length(mp_list)>0) {
    shared_pairs_all[[length(shared_pairs_all)+1]] <- rbindlist(mp_list)
  }
}

sum_dt <- rbindlist(count_list, use.names=TRUE, fill=TRUE)
shared_pairs <- rbindlist(shared_pairs_all, use.names=TRUE, fill=TRUE)

# write a NEW blacklist-aware summary
fwrite(sum_dt, file.path(dir_plot, "all_samples_exact_event_summary.tsv"), sep="\t")

# ---------- combine per-sample per-event files into one ----------
detail_files <- file.path(dir_plot, sprintf("%s_LRS_match_detail.tsv", samples))
detail_files <- detail_files[file.exists(detail_files)]
if (length(detail_files) > 0) {
  all_detail <- rbindlist(lapply(detail_files, fread), use.names=TRUE, fill=TRUE)
  fwrite(all_detail,
         file.path(dir_plot, "all_samples_LRS_match_detail.tsv"),
         sep="\t")
}

# =========================
# 1) Stacked bar for 8 samples (blacklist aware)
# =========================
plot_dt <- melt(
  sum_dt[, .(sample, HiC_only=N_HiC_only, Shared=N_Shared, LRS_only=N_LR_only)],
  id.vars = "sample",
  variable.name = "category",
  value.name = "count"
)
plot_dt[, category := factor(category, levels=c("LRS_only","Shared","HiC_only"))]
plot_dt[, sample := factor(sample, levels = sample_order)]  

p_bar <- ggplot(plot_dt, aes(x=sample, y=count, fill=category)) +
  geom_col(width=0.7) +
  theme_classic() +   
  theme(panel.grid = element_blank(),  
        axis.line.x = element_line(),   
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1))+ 
  labs(x="", y="Number of inter-translocations",
       title="Hi-C vs LRS translocations") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  scale_fill_manual(values = CNS_COLORS)

p_bar

ggsave(file.path(dir_plot, "stacked_bar_all_samples.pdf"), p_bar, width=8, height=5)
ggsave(file.path(dir_plot, "stacked_bar_all_samples.png"), p_bar, width=8, height=5, dpi=300)

# =========================
# 2) Per-sample Venn diagrams (blacklist aware)
# =========================
venn_ok <- requireNamespace("VennDiagram", quietly=TRUE)

for (s in samples) {
  
  s_summ <- sum_dt[sample==s]
  if (nrow(s_summ)==0) next
  
  LRS_only <- s_summ$N_LR_only
  HiC_only <- s_summ$N_HiC_only
  Shared   <- s_summ$N_Shared
  
  out_pdf <- file.path(dir_plot, sprintf("venn_%s.pdf", s))
  
  if (venn_ok) {
    library(VennDiagram)
    pdf(out_pdf, width=5, height=5)
    grid::grid.newpage()
    draw.pairwise.venn(
      area1 = LRS_only + Shared,
      area2 = HiC_only + Shared,
      cross.area = Shared,
      category = c("Long-read (merged)", "Hi-C"),
      fill = c("lightblue", "salmon"),
      lty = "blank",
      cex = 1.2,
      cat.cex = 1.2,
      cat.pos = c(-20, 20),
      euler.d = FALSE,
      scaled = FALSE
    )
    dev.off()
  } else {
    p_txt <- ggplot() +
      theme_void() +
      annotate("text", x=0, y=0.2,
               label=sprintf("%s\nLRS only: %d\nHi-C only: %d\nShared: %d",
                             s, LRS_only, HiC_only, Shared),
               size=6)
    ggsave(out_pdf, p_txt, width=4, height=3)
  }
}

# =========================
# 2b) Global Venn diagram (NO blacklist)
# overlap should be 43
# =========================
venn_ok <- requireNamespace("VennDiagram", quietly=TRUE)
if (venn_ok) {
  library(VennDiagram)
  
  global_HiC_total <- sum(sum_dt$N_HiC_total)
  global_LRS_total <- sum(sum_dt$N_LR_total)
  global_shared    <- sum(sum_dt$N_Shared)  
  
  out_pdf <- file.path(dir_plot, "venn_all_samples.pdf")
  out_png <- file.path(dir_plot, "venn_all_samples.png")
  
  pdf(out_pdf, width=5, height=5)
  grid::grid.newpage()
  draw.pairwise.venn(
    area1      = global_LRS_total,
    area2      = global_HiC_total,
    cross.area = global_shared,
    category   = c("Long-read (merged)", "Hi-C"),
    fill       = c("lightblue", "salmon"),
    lty        = "blank",
    cex        = 1.5,
    cat.cex    = 1.2,
    cat.pos    = c(-20, 20),
    euler.d    = FALSE,
    scaled     = FALSE
  )
  dev.off()
  
  png(out_png, width=1600, height=1600, res=300)
  grid::grid.newpage()
  draw.pairwise.venn(
    area1      = global_LRS_total,
    area2      = global_HiC_total,
    cross.area = global_shared,
    category   = c("Long-read (merged)", "Hi-C"),
    fill       = c("lightblue", "salmon"),
    lty        = "blank",
    cex        = 1.5,
    cat.cex    = 1.2,
    cat.pos    = c(-20, 20),
    euler.d    = FALSE,
    scaled     = FALSE
  )
  dev.off()
  
  message(sprintf("[INFO] Global overlap (NO blacklist) = %d (expect 43).", global_shared))
}

# =========================
# 3) Orientation match among shared LRS events
#    (EVENT-level only)

# =========================
if (nrow(shared_pairs)==0) {
  message("[WARN] No shared LRS events found for orientation plots.")
} else {

  
  shared_pairs <- shared_pairs[!(sample=="PT3" & lr_id %in% EXCLUDE$PT3$lr_ids)]
  shared_pairs <- shared_pairs[!(sample=="PT3" & hic_id %in% EXCLUDE$PT3$hic_ids)]
  
  
  # pair-level orientation agreement (internal only, used to derive event-level)
  shared_pairs[, orient_same := (hic_strandA==lr_strandA & hic_strandB==lr_strandB)]
  
  # ---------- EVENT-level orientation ----------
  orient_dt <- shared_pairs[, .(
    orientation_match = any(orient_same),
    n_hic_support = uniqueN(hic_id)
  ), by=.(sample, lr_id)]
  
  # per-event TSV (event-level only)
  fwrite(orient_dt,
         file.path(dir_plot, "orientation_match_per_event.tsv"),
         sep="\t")
  
  # also split to matched / mismatched EVENTS with breakpoint info
  detail_path <- file.path(dir_plot, "all_samples_LRS_match_detail.tsv")
  all_detail <- fread(detail_path)
  lr_shared_detail <- all_detail[is_shared == TRUE]
  
  orient_with_bp <- merge(
    lr_shared_detail,
    orient_dt,
    by.x = c("sample","lr_id"),
    by.y = c("sample","lr_id"),
    all.x = TRUE
  )
  orient_with_bp[is.na(orientation_match), orientation_match := FALSE]
  
  fwrite(
    orient_with_bp[orientation_match==TRUE,
                   .(sample, lr_id, chrA, posA, strandA, chrB, posB, strandB, n_hic_support)],
    file.path(dir_plot, "orientation_match_events_with_breakpoints.tsv"),
    sep="\t"
  )
  
  # ---------- mismatch EVENTS with BOTH LRS + Hi-C breakpoints ----------

  mismatch_lr <- orient_with_bp[orientation_match==FALSE,
                                .(sample, lr_id, chrA, posA, strandA, chrB, posB, strandB, n_hic_support)]
  

  mismatch_pairs <- shared_pairs[lr_id %in% mismatch_lr$lr_id]
  

  hic_agg <- mismatch_pairs[, .(
    hic_id_list       = paste(unique(hic_id), collapse=";"),
    hic_chrA_list     = paste(unique(hic_chrA), collapse=";"),
    hic_posA_list     = paste(unique(hic_posA), collapse=";"),
    hic_strandA_list  = paste(unique(hic_strandA), collapse=";"),
    hic_chrB_list     = paste(unique(hic_chrB), collapse=";"),
    hic_posB_list     = paste(unique(hic_posB), collapse=";"),
    hic_strandB_list  = paste(unique(hic_strandB), collapse=";")
  ), by=.(sample, lr_id)]
  

  mismatch_full <- merge(mismatch_lr, hic_agg, by=c("sample","lr_id"), all.x=TRUE)
  
  fwrite(
    mismatch_full,
    file.path(dir_plot, "orientation_mismatch_events_with_breakpoints.tsv"),
    sep="\t"
  )
  
  # ---------- global counts (EVENT-level) ----------
  orient_sum <- orient_dt[, .N, by=orientation_match]
  orient_sum[, label := ifelse(orientation_match, "Orientation matched", "Orientation mismatched")]
  
  fwrite(orient_sum[, .(label,N)],
         file.path(dir_plot, "orientation_match_counts.tsv"),
         sep="\t")
  
  # global PIE (EVENT-level)
  orient_sum[, prop := N / sum(N)]
  orient_sum[, pct  := sprintf("%.1f%%", prop * 100)]
  orient_sum[, y_pos := cumsum(prop) - prop/2]
  
  p_orient_global <- ggplot(orient_sum, aes(x="", y=prop, fill=label)) +
    geom_col(width=1, color="white") +
    coord_polar(theta="y") +
    theme_void(base_size=12) +
    labs(title="Orientation agreement in shared translocations") +
    geom_text(aes(y=y_pos, label=paste0(label, "\n", N, " (", pct, ")")), size=4) +
    scale_fill_manual(values = CNS_COLORS)+labs(fill = NULL)
  
  p_orient_global
  
  
  
  ggsave(file.path(dir_plot, "orientation_match_shared_LRS.pdf"),
         p_orient_global, width=6, height=5)
  ggsave(file.path(dir_plot, "orientation_match_shared_LRS.png"),
         p_orient_global, width=6, height=5, dpi=300)
  
  # ---------- per-sample counts (EVENT-level) ----------
  per_sample_sum <- orient_dt[, .N, by=.(sample, orientation_match)]
  per_sample_sum[, label := ifelse(orientation_match, "Matched", "Mismatched")]
  
  fwrite(per_sample_sum[, .(sample,label,N)],
         file.path(dir_plot, "orientation_match_per_sample.tsv"),
         sep="\t")
  

  per_sample_sum[, sample := factor(sample, levels = sample_order)]
  
  p_orient_sample <- ggplot(per_sample_sum, aes(x=sample, y=N, fill=label)) +
    geom_col(position="stack", width=0.7) +
    theme_classic() +  
    theme(panel.grid = element_blank(),  
          axis.line.x = element_line(),  
          axis.line.y = element_line(),
          axis.text.x = element_text(angle = 45, hjust = 1))+
    labs(x="", y="Number of high-confidence translocations",
         title="Orientation agreement per sample") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_manual(values = CNS_COLORS) +
    scale_y_continuous(
      breaks = function(x) seq(0, 12, by = 2),
      limits = function(x) c(0, 12)
    )
  p_orient_sample
  

  
  ggsave(file.path(dir_plot, "orientation_match_per_sample.pdf"),
         p_orient_sample, width=8, height=5)
  ggsave(file.path(dir_plot, "orientation_match_per_sample.png"),
         p_orient_sample, width=8, height=5, dpi=300)
}


library(data.table)
library(gridExtra)
library(grid)

setwd("F:/zer/TransFinder/6_Integration")
dir_plot <- "3_plot"


mismatch <- fread(file.path(dir_plot, "orientation_mismatch_events_with_breakpoints.tsv"))


table_fig <- mismatch[, .(
  Sample = sample,
  `LRS translocation` = sprintf(
    "%s:%d (%s) - %s:%d (%s)",
    chrA, posA, strandA,
    chrB, posB, strandB
  ),
  `Hi-C translocation` = sprintf(
    "%s:%s (%s) - %s:%s (%s)",
    hic_chrA_list, hic_posA_list, hic_strandA_list,
    hic_chrB_list, hic_posB_list, hic_strandB_list
  ),
  `LRS ori`  = paste0(strandA, strandB),
  `Hi-C ori` = paste0(hic_strandA_list, hic_strandB_list)
)]


pdf(file.path(dir_plot, "orientation_mismatch_table.pdf"), width=11, height=3)
grid.newpage()
grid.table(table_fig)
dev.off()

message("PDF: 3_plot/orientation_mismatch_table.pdf")
