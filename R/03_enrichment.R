# GO / KEGG over-representation analysis (ORA) + hypergeometric test of the theme gene sets -------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Rn.eg.db)
  library(ggplot2)
})
dir.create(out("03_enrichment"), showWarnings = FALSE)

res_list  <- load_rds("res_list")
gene_sets <- load_rds("gene_sets")

# The KEGG annotation is downloaded once and cached (results/kegg_cache/; the first run takes about 2-3 minutes).
# the KEGG REST link file is large and often exceeds the default 60 s download timeout of R, which gives "cannot read from connection" and looks like a network failure; enrichKEGG also re-downloads on every call.
# The modification time of the cache file is the KEGG access date (cite it in the Methods); delete the file to refresh KEGG.
options(timeout = 1800)
kegg_file <- file.path(OUT_DIR, "kegg_cache", "rno_KEGG.rds")
dir.create(dirname(kegg_file), showWarnings = FALSE)
if (!file.exists(kegg_file)) {
  kegg_db <- tryCatch(download_KEGG("rno"),
                      error = function(e) { message("KEGG download failed: ", conditionMessage(e)); NULL })
  if (!is.null(kegg_db)) saveRDS(kegg_db, kegg_file)
} else {
  kegg_db <- readRDS(kegg_file)
}
if (!is.null(kegg_db)) cat("KEGG annotation date:", format(file.mtime(kegg_file), "%Y-%m-%d"), "\n")

save_plot <- function(p, file, w = 9, h = 7) ggsave(out("03_enrichment", file), p, width = w, height = h)

for (cmp in c("SPMS_vs_MS", "SPMS_vs_PMS")) {
  res <- res_list[[cmp]]
  universe <- rownames(res)   # background: all genes that passed DESeq2 independent filtering
  deg <- list(
    All  = rownames(res)[is_deg(res)],
    Up   = rownames(res)[is_deg(res) & res$log2FoldChange > 0],
    Down = rownames(res)[is_deg(res) & res$log2FoldChange < 0]
  )

  for (dir in names(deg)) {
    tag <- paste0(cmp, "_", dir)

    # GO
    for (ont in c("BP", "CC", "MF")) {
      ego <- enrichGO(gene = deg[[dir]], universe = universe, OrgDb = org.Rn.eg.db,
                      keyType = "SYMBOL", ont = ont, pAdjustMethod = "BH",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (is.null(ego) || nrow(ego) == 0) next
      write.csv(as.data.frame(ego), out("03_enrichment", paste0("GO_", ont, "_", tag, ".csv")), row.names = FALSE)
      save_plot(dotplot(ego, showCategory = 15) + ggtitle(paste("GO", ont, "-", gsub("_", " ", tag))),
                paste0("GO_", ont, "_", tag, "_dotplot.pdf"))
    }

    # KEGG: hypergeometric test on the cached annotation; results are identical to enrichKEGG (verified)
    ekegg <- NULL
    if (!is.null(kegg_db)) {
      entrez <- bitr(deg[[dir]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)$ENTREZID
      bg     <- bitr(universe,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)$ENTREZID
      ekegg <- enricher(gene = entrez, universe = bg, TERM2GENE = kegg_db$KEGGPATHID2EXTID,
                        TERM2NAME = kegg_db$KEGGPATHID2NAME, pvalueCutoff = 0.05, qvalueCutoff = 0.2)
    }
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      ekegg <- setReadable(ekegg, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
      write.csv(as.data.frame(ekegg), out("03_enrichment", paste0("KEGG_", tag, ".csv")), row.names = FALSE)
      save_plot(dotplot(ekegg, showCategory = 15) + ggtitle(paste("KEGG -", gsub("_", " ", tag))),
                paste0("KEGG_", tag, "_dotplot.pdf"))
    }
  }

  # Are the theme gene sets enriched among DEGs? (hypergeometric test)
  # The original script mistakenly used the number of DEGs in the pathway as m; here m is the number of pathway genes in the background
  N <- length(universe); k <- length(deg$All)
  hyper <- do.call(rbind, lapply(names(gene_sets$theme), function(nm) {
    pw <- intersect(gene_sets$theme[[nm]], universe)
    hit <- intersect(pw, deg$All)
    data.frame(Pathway = nm, Genes_in_background = length(pw), DEG_in_pathway = length(hit),
               Up = sum(res[hit, "log2FoldChange"] > 0), Down = sum(res[hit, "log2FoldChange"] < 0),
               Expected = round(k * length(pw) / N, 1),
               P = phyper(length(hit) - 1, length(pw), N - length(pw), k, lower.tail = FALSE))
  }))
  hyper$Padj <- p.adjust(hyper$P, "BH")
  write.csv(hyper, out("03_enrichment", paste0("ThemeSets_hypergeometric_", cmp, ".csv")), row.names = FALSE)
  cat("\n==", cmp, "theme gene set hypergeometric test ==\n"); print(hyper)
}
