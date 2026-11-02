# Export pathway DEGs + STRING PPI network (for Cytoscape)-----------------------------
source("R/00_config.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(STRINGdb)
})
dir.create(out("06_PPI"), showWarnings = FALSE)

res       <- load_rds("res_list")[[MAIN_COMPARISON]]
gene_sets <- load_rds("gene_sets")
pathways  <- c(gene_sets$theme, gene_sets$signaling)

# DEGs of each pathway (some steps of the original scripts used the unadjusted pvalue; padj is used throughout here)-------------
deg_genes <- rownames(res)[is_deg(res)]
pathway_deg <- do.call(rbind, lapply(names(pathways), function(pw) {
  g <- intersect(pathways[[pw]], deg_genes)
  if (length(g) == 0) return(NULL)
  data.frame(Pathway = pw, Gene = g, log2FoldChange = res[g, "log2FoldChange"], padj = res[g, "padj"])
})) %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(Pathway, desc(log2FoldChange))
write.csv(pathway_deg, out("06_PPI", paste0("pathway_DEGs_long_", MAIN_COMPARISON, ".csv")), row.names = FALSE)

# One row per gene, listing all pathways it belongs to (a gene can belong to several pathways)
node_attr <- pathway_deg %>%
  group_by(Gene) %>%
  summarise(log2FoldChange = first(log2FoldChange), padj = first(padj),
            Direction = first(Direction), Pathways = paste(Pathway, collapse = ";"),
            Primary_pathway = Pathway[which.min(match(Pathway, names(pathways)))],
            .groups = "drop")
write.csv(node_attr, out("06_PPI", paste0("pathway_DEGs_node_attributes_", MAIN_COMPARISON, ".csv")), row.names = FALSE)

deg_count <- pathway_deg %>% count(Pathway, Direction) %>% tidyr::pivot_wider(names_from = Direction, values_from = n, values_fill = 0)
print(deg_count)

# STRING (the first run needs network access to download the database)---------------------------------------------
dir.create(file.path(OUT_DIR, "string_cache"), showWarnings = FALSE)
options(timeout = 1800)   # the protein.links file is about 90 MB and would time out at the default 60 s
string_db <- tryCatch(
  STRINGdb$new(version = STRING_VERSION, species = STRING_SPECIES,
               score_threshold = STRING_SCORE, input_directory = file.path(OUT_DIR, "string_cache")),
  error = function(e) { message("STRINGdb initialization failed: ", conditionMessage(e)); NULL })

export_network <- function(genes, prefix) {
  mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
  cat(prefix, ": ", nrow(mapped), "/", length(genes), " genes mapped to STRING\n", sep = "")
  edges <- string_db$get_interactions(mapped$STRING_id) %>%
    distinct(from, to, .keep_all = TRUE) %>%
    left_join(mapped, by = c("from" = "STRING_id")) %>% rename(from_gene = gene) %>%
    left_join(mapped, by = c("to" = "STRING_id"))   %>% rename(to_gene = gene) %>%
    select(from_gene, to_gene, combined_score)
  write.table(edges, out("06_PPI", paste0(prefix, "_edges.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  edges:", nrow(edges), "\n")
}

if (!is.null(string_db)) {
  # Network of all pathway DEGs
  export_network(unique(pathway_deg$Gene), paste0("PPI_pathway_DEGs_", MAIN_COMPARISON))

  # Network of the leading-edge genes of the osteogenesis signaling pathways (GSEA)
  gsea <- load_rds("gsea_all")[[MAIN_COMPARISON]]
  le <- unique(unlist(strsplit(gsea$leadingEdge[gsea$pathway %in% names(gene_sets$signaling)], ";")))
  export_network(le, paste0("PPI_signaling_leadingEdge_", MAIN_COMPARISON))
}
