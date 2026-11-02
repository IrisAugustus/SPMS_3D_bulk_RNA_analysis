# Global configuration ------------------------------------------------------------------
# Paths and thresholds shared by all scripts are set here so that results stay consistent

# Paths (relative to the project root; run the scripts from the project root)
if (!file.exists(file.path("R", "00_config.R")))
  stop("Set the working directory to the project root (the folder containing R/ and data/) first. Current directory: ", getwd())
OUT_DIR <- "results"
RDS_DIR <- file.path(OUT_DIR, "rds")

COUNT_FILE    <- file.path("data", "gene.count.matrix.annotation.xls")  # raw count matrix (not included in the repository; place it in data/)
GENESET4_FILE <- file.path("data", "all_acquired_genes.csv")            # four theme gene sets (keyword search of MSigDB C5)

# Samples (columns 5-13 of the matrix; SPMS_Zhu1-3 are not used)
SAMPLES <- c("MS1", "MS2", "MS3", "PMS1", "PMS2", "PMS3", "SPMS1", "SPMS2", "SPMS3")
GROUPS  <- factor(rep(c("MS", "PMS", "SPMS"), each = 3), levels = c("MS", "PMS", "SPMS"))
GROUP_COLORS <- c("MS" = "#2E802E", "PMS" = "#FFA500", "SPMS" = "#E91919")

# Comparisons (treatment, control)
COMPARISONS <- list(
  SPMS_vs_MS  = c("SPMS", "MS"),
  SPMS_vs_PMS = c("SPMS", "PMS"),
  PMS_vs_MS   = c("PMS", "MS")
)
MAIN_COMPARISON <- "SPMS_vs_MS"   # comparison used for heatmaps and PPI

# DEG thresholds (shared by all steps)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# GSEA
GSEA_RANK_BY  <- "log2FC"   # "log2FC" (as in the original analysis) or "stat" (DESeq2 Wald statistic)
GSEA_MIN_SIZE <- 3
GSEA_MAX_SIZE <- 3000       # originally 1000, which dropped Cell-ECM / Cell-cell
SEED          <- 123

# Genes manually excluded from the heatmaps (set per pathway).
# Matches the two final figures (deseq2/heatmap_Osteogenesis_all_genes_MS_PMS_SPMS.pdf,
# deseq2/heatmap_cell_ecm_all_genes_MS_PMS_SPMS.pdf); the exclusion rationale must be stated in the paper.
# Statistics of the excluded genes are written to results/05_heatmaps/heatmap_excluded_genes.csv,
# and a version without exclusion (*_noExclusion.pdf) is also saved for comparison. Set to list() to exclude nothing.
HEATMAP_EXCLUDE <- list(
  Osteogenesis = c("Tacr1", "Lrp3"),
  Cell_ECM     = c("Cd200", "Cxadr", "Flrt2", "Nox1", "Tgm2", "Epha1", "Kif26b",
                   "Ripor2", "Il1a", "Il23r", "Cd72")
)

# STRING
STRING_VERSION <- "11.5"
STRING_SPECIES <- 10116     # Rattus norvegicus
STRING_SCORE   <- 400

dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
options(stringsAsFactors = FALSE)

out <- function(...) file.path(OUT_DIR, ...)
load_rds <- function(name) readRDS(file.path(RDS_DIR, paste0(name, ".rds")))
save_rds <- function(obj, name) saveRDS(obj, file.path(RDS_DIR, paste0(name, ".rds")))

is_deg <- function(res) !is.na(res$padj) & res$padj < PADJ_CUTOFF & abs(res$log2FoldChange) > LFC_CUTOFF
