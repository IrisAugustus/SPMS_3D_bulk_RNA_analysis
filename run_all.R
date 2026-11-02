steps <- c(
  "01_deseq2.R",
  "02_gene_sets.R",
  "03_enrichment.R",
  "04_gsea.R",
  "05_heatmaps.R",
  "06_pathway_DEGs_PPI.R"
)

for (s in steps) {
  cat("\n================ ", s, " ================\n")
  source(file.path("R", s), local = new.env())
}

writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
cat("\n Code running finished, see in results/\n")
