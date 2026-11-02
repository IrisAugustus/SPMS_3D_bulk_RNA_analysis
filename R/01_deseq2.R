# Preprocessing + DESeq2 differential expression + PCA + volcano plots -------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})
dir.create(out("01_DEG"), showWarnings = FALSE)

# Read the count matrix ------------------------------------------------------------
geneCount <- read.delim(COUNT_FILE, check.names = FALSE)
stopifnot(all(SAMPLES %in% colnames(geneCount)))

allcounts <- geneCount[, SAMPLES]
gene_name <- geneCount[["Gene name"]]

# Drop rows without a gene name; keep the first entry of duplicated names
keep <- !is.na(gene_name) & gene_name != "" & !duplicated(gene_name)
cat("Genes in raw matrix:", nrow(allcounts),
    "| without gene name:", sum(is.na(gene_name) | gene_name == ""),
    "| duplicated names:", sum(duplicated(gene_name) & gene_name != ""), "\n")
allcounts <- allcounts[keep, ]
rownames(allcounts) <- gene_name[keep]
allcounts <- round(as.matrix(allcounts))   # round RSEM expected counts to integers
cat("Genes kept:", nrow(allcounts), "\n")

# DESeq2 ---------------------------------------------------------------------
coldata <- data.frame(row.names = SAMPLES, condition = GROUPS)
dds <- DESeqDataSetFromMatrix(countData = allcounts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

res_list <- lapply(COMPARISONS, function(cmp) {
  res <- results(dds, contrast = c("condition", cmp[1], cmp[2]), alpha = PADJ_CUTOFF)
  res <- as.data.frame(res)
  res[!is.na(res$padj), ]
})

# Result tables and DEG summary -----------------------------------------------------------
deg_summary <- data.frame()
for (cmp in names(res_list)) {
  res <- res_list[[cmp]]
  res_out <- res %>% rownames_to_column("SYMBOL") %>% arrange(padj)
  write.csv(res_out, out("01_DEG", paste0("DESeq2_all_", cmp, ".csv")), row.names = FALSE)

  deg <- res_out[is_deg(res_out), ]
  write.csv(filter(deg, log2FoldChange > 0), out("01_DEG", paste0("DEG_up_", cmp, ".csv")), row.names = FALSE)
  write.csv(filter(deg, log2FoldChange < 0), out("01_DEG", paste0("DEG_down_", cmp, ".csv")), row.names = FALSE)

  deg_summary <- rbind(deg_summary, data.frame(
    Comparison = cmp,
    Up   = sum(deg$log2FoldChange > 0),
    Down = sum(deg$log2FoldChange < 0),
    Total = nrow(deg)
  ))
}
deg_summary$Criteria <- sprintf("padj < %s & |log2FC| > %s", PADJ_CUTOFF, LFC_CUTOFF)
write.csv(deg_summary, out("01_DEG", "DEG_summary.csv"), row.names = FALSE)
print(deg_summary)

# PCA (vst-transformed)----------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"), 1)
p_pca <- ggplot(pca, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = name), size = 3, show.legend = FALSE) +
  scale_color_manual(values = GROUP_COLORS) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), color = "Group") +
  theme_bw()
ggsave(out("01_DEG", "PCA_vst.pdf"), p_pca, width = 6, height = 5)

# Volcano plots ---------------------------------------------------------------------
for (cmp in names(res_list)) {
  vd <- res_list[[cmp]] %>%
    rownames_to_column("SYMBOL") %>%
    mutate(Significance = case_when(
      is_deg(.) & log2FoldChange > 0 ~ "Up",
      is_deg(.) & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    ))
  n_up <- sum(vd$Significance == "Up"); n_down <- sum(vd$Significance == "Down")
  top <- vd %>% filter(Significance != "NS") %>% arrange(padj) %>% head(15)

  p <- ggplot(vd, aes(log2FoldChange, -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_text_repel(data = top, aes(label = SYMBOL), size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = c(Up = "#C82423", Down = "#2878B5", NS = "grey70"),
                       labels = c(Up = paste0("Up (", n_up, ")"), Down = paste0("Down (", n_down, ")"), NS = "NS")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", alpha = 0.5) +
    labs(title = gsub("_", " ", cmp), x = "log2 fold change", y = "-log10 adjusted P") +
    theme_bw()
  ggsave(out("01_DEG", paste0("Volcano_", cmp, ".pdf")), p, width = 7, height = 6)
}

# Save for the following steps ---------------------------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
save_rds(dds, "dds")
save_rds(res_list, "res_list")
save_rds(log2(norm_counts + 1), "log2_norm_counts")
