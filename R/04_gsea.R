# GSEA: four theme gene sets + osteogenesis-related signaling pathways ------------------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(fgsea)
  library(clusterProfiler)
  library(GseaVis)
  library(ggplot2)
})
dir.create(out("04_GSEA"), showWarnings = FALSE)

res_list  <- load_rds("res_list")
gene_sets <- load_rds("gene_sets")
pathways  <- c(gene_sets$theme, gene_sets$signaling)

make_ranks <- function(res) {
  r <- if (GSEA_RANK_BY == "stat") res$stat else res$log2FoldChange
  names(r) <- rownames(res)
  sort(r[!is.na(r)], decreasing = TRUE)
}

gsea_all <- list()
for (cmp in c("SPMS_vs_MS", "SPMS_vs_PMS")) {
  ranks <- make_ranks(res_list[[cmp]])

  # fgsea (multilevel): used to report NES / P / padj
  set.seed(SEED)
  fg <- fgsea(pathways = pathways, stats = ranks, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE)
  fg <- fg[order(fg$pval), ]
  fg_df <- as.data.frame(fg)
  fg_df$leadingEdge <- sapply(fg_df$leadingEdge, paste, collapse = ";")
  fg_df$Comparison <- cmp
  write.csv(fg_df, out("04_GSEA", paste0("fgsea_", cmp, ".csv")), row.names = FALSE)
  gsea_all[[cmp]] <- fg_df
  cat("\n==", cmp, "==\n"); print(fg_df[, c("pathway", "size", "NES", "pval", "padj")])

  # clusterProfiler::GSEA: used only for GseaVis plots
  term2gene <- data.frame(term = rep(names(pathways), lengths(pathways)), gene = unlist(pathways))
  set.seed(SEED)
  gs <- GSEA(geneList = ranks, TERM2GENE = term2gene, minGSSize = GSEA_MIN_SIZE,
             maxGSSize = GSEA_MAX_SIZE, pvalueCutoff = 1, eps = 0, seed = TRUE, verbose = FALSE)

  for (pw in intersect(names(pathways), gs@result$ID)) {
    p <- gseaNb(object = gs, geneSetID = pw, subPlot = 3, lineSize = 1,
                curveCol = c("#76BA99", "#EB4747", "#996699"), segCol = "black",
                addPval = TRUE, pvalX = 0.65, pvalY = 0.85, pvalSize = 4,
                termWidth = 35, rankSeq = 5000, htHeight = 0.3)
    ggsave(out("04_GSEA", paste0("gseaNb_", pw, "_", cmp, ".pdf")), p, width = 8, height = 6)
  }
}

summary_df <- do.call(rbind, gsea_all)[, c("Comparison", "pathway", "size", "NES", "pval", "padj")]
write.csv(summary_df, out("04_GSEA", "GSEA_summary.csv"), row.names = FALSE)
save_rds(gsea_all, "gsea_all")

# Sensitivity to the ranking method -----------------------------------------------------------
# The results are sensitive to the ranking: the original analysis also ranked low-expression genes removed by DESeq2 independent filtering (padj = NA),
# whose log2FC values are very noisy. Results for three ranking schemes are given side by side so that one can be chosen and justified in the paper.
suppressPackageStartupMessages(library(DESeq2))
dds <- load_rds("dds")
schemes <- list(
  log2FC_all_genes      = function(r) setNames(r$log2FoldChange, rownames(r)),
  log2FC_filtered_genes = function(r) { r <- r[!is.na(r$padj), ]; setNames(r$log2FoldChange, rownames(r)) },
  stat_filtered_genes   = function(r) { r <- r[!is.na(r$padj), ]; setNames(r$stat, rownames(r)) }
)
sens <- do.call(rbind, lapply(c("SPMS_vs_MS", "SPMS_vs_PMS"), function(cmp) {
  r <- as.data.frame(results(dds, contrast = c("condition", COMPARISONS[[cmp]])))
  do.call(rbind, lapply(names(schemes), function(sc) {
    ranks <- schemes[[sc]](r)
    ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
    set.seed(SEED)
    fg <- fgsea(pathways = pathways, stats = ranks, minSize = GSEA_MIN_SIZE, maxSize = GSEA_MAX_SIZE)
    data.frame(Comparison = cmp, Ranking = sc, pathway = fg$pathway, NES = fg$NES, pval = fg$pval, padj = fg$padj)
  }))
}))
sens_wide <- tidyr::pivot_wider(sens, id_cols = c(Comparison, pathway), names_from = Ranking,
                                values_from = c(NES, pval, padj), names_glue = "{Ranking}_{.value}")
write.csv(sens_wide, out("04_GSEA", "GSEA_ranking_sensitivity.csv"), row.names = FALSE)
