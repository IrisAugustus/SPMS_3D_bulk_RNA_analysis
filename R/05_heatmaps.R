# Pathway DEG heatmaps (MS / PMS / SPMS)--------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(pheatmap)
})
dir.create(out("05_heatmaps"), showWarnings = FALSE)

res       <- load_rds("res_list")[[MAIN_COMPARISON]]
expr      <- load_rds("log2_norm_counts")   # log2(DESeq2 normalized counts + 1)
gene_sets <- load_rds("gene_sets")
pathways  <- c(gene_sets$theme, gene_sets$signaling)

ann_col    <- data.frame(Group = GROUPS, row.names = SAMPLES)
ann_colors <- list(Group = GROUP_COLORS)
breaks     <- seq(-2, 2, length.out = 101)   # fixed color range of the row Z-score: -2 to 2
max_genes  <- 200

draw_heatmap <- function(pw, exclude, file) {
  genes <- intersect(pathways[[pw]], rownames(res))
  d <- res[genes, ]
  d <- d[is_deg(d) & !rownames(d) %in% exclude, ]
  if (nrow(d) < 2) { cat(pw, ": fewer than 2 DEGs, skipped\n"); return(NULL) }

  # Order by log2FC, descending: upregulated genes on top, downregulated below
  d <- d[order(-d$log2FoldChange), ]
  up   <- rownames(d)[d$log2FoldChange > 0]
  down <- rownames(d)[d$log2FoldChange < 0]
  # If there are too many genes, take about half from each direction (upregulated: largest FC; downregulated: first in the ordering)
  if (length(up) + length(down) > max_genes) {
    n_up <- min(length(up), max_genes %/% 2)
    up   <- head(up, n_up)
    down <- head(down, max_genes - n_up)
  }
  mat <- expr[c(up, down), SAMPLES, drop = FALSE]

  n <- nrow(mat)
  pdf(file, width = 7, height = max(5, n * 0.16 + 2))
  pheatmap(mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = ann_col, annotation_colors = ann_colors,
           color = colorRampPalette(c("#2878B5", "white", "#C82423"))(100), breaks = breaks,
           legend_breaks = -2:2, border_color = NA, angle_col = 45,
           gaps_row = if (length(up) > 0 && length(down) > 0) length(up) else NULL,
           fontsize_row = if (n > 80) 6 else if (n > 40) 8 else 10,
           main = sprintf("%s DEGs (%s; up %d, down %d)", gsub("_", "-", pw),
                          gsub("_", " ", MAIN_COMPARISON), length(up), length(down)))
  dev.off()
  data.frame(Pathway = pw, Up = length(up), Down = length(down), Genes = paste(c(up, down), collapse = ";"))
}

heatmap_log <- data.frame()
for (pw in names(pathways)) {
  exclude <- HEATMAP_EXCLUDE[[pw]]
  stem <- out("05_heatmaps", paste0("heatmap_DEG_", pw, "_", MAIN_COMPARISON))
  heatmap_log <- rbind(heatmap_log, draw_heatmap(pw, exclude, paste0(stem, ".pdf")))
  if (length(exclude) > 0) draw_heatmap(pw, character(0), paste0(stem, "_noExclusion.pdf"))
}
heatmap_log$Excluded <- sapply(heatmap_log$Pathway, function(pw) paste(HEATMAP_EXCLUDE[[pw]], collapse = ";"))
write.csv(heatmap_log, out("05_heatmaps", "heatmap_gene_counts.csv"), row.names = FALSE)
print(heatmap_log[, c("Pathway", "Up", "Down", "Excluded")])

# Statistics of the excluded genes (to document the exclusion rationale in the paper)
if (length(HEATMAP_EXCLUDE) > 0) {
  excl <- do.call(rbind, lapply(names(HEATMAP_EXCLUDE), function(pw) {
    g <- HEATMAP_EXCLUDE[[pw]]
    data.frame(Pathway = pw, Gene = g,
               baseMean = res[g, "baseMean"], log2FoldChange = res[g, "log2FoldChange"], padj = res[g, "padj"],
               Is_DEG = g %in% rownames(res)[is_deg(res)])
  }))
  write.csv(excl, out("05_heatmaps", "heatmap_excluded_genes.csv"), row.names = FALSE)
}
