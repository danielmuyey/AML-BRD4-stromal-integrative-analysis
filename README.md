############################################################
# AML BRD4-Sensitive Stromal Remodeling Analysis Pipeline
# ALL-IN-ONE R SCRIPT
#
# Integrative transcriptomics, differential expression, PCA,
# directional reversal, GO/KEGG enrichment, STRING/PPI,
# candidate prioritization, GSE107490 human validation,
# TCGA-LAML exploratory analysis, and figure generation.
#
# Author: Daniel Muteb Muyey
# Repository: AML-BRD4-stromal-integrative-analysis
# Version: 1.0.0
#
# Public datasets:
#   GSE148625  - murine AML-associated Nestin+ BM MSC discovery
#   GSE101449  - human MSC ARV-825 / BRD4 perturbation
#   GSE107490  - human AML/MDS MSC validation
#   TCGA-LAML  - exploratory clinical/molecular context
#
# IMPORTANT:
# This script preserves the conservative reproducibility design.
# It does not invent missing raw expression matrices, metadata,
# orthologue mappings, STRING settings, or the exact upstream
# TCGA nine-gene score-construction implementation.
############################################################

options(stringsAsFactors = FALSE)
set.seed(2026)

dir.create("results", showWarnings = FALSE)
dir.create("results/DEGs", recursive = TRUE, showWarnings = FALSE)
dir.create("results/reversal", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/TCGA", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)



####################################################################
# MODULE: 01_load_data.R
####################################################################
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

read_if_exists <- function(path) {
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}

# Public datasets used in the manuscript
# Raw/processed source files should be downloaded from GEO and placed in data/.
gse148625_expr <- read_if_exists("data/GSE148625_expression.csv")
gse148625_meta <- read_if_exists("data/GSE148625_metadata.csv")

gse101449_expr <- read_if_exists("data/GSE101449_expression.csv")
gse101449_meta <- read_if_exists("data/GSE101449_metadata.csv")

gse107490_expr <- read_if_exists("data/GSE107490_expression.csv")
gse107490_meta <- read_if_exists("data/GSE107490_metadata.csv")

# Verified/processed result tables, when supplied
deg_results <- read_if_exists("data/deg_results.csv")
reversal_candidates <- read_if_exists("data/Table_3_reversal_candidates.csv")
enrichment_results <- read_if_exists("data/Table_4_GO_KEGG_enrichment.csv")
ppi_nodes <- read_if_exists("data/ppi_nodes.csv")
ppi_edges <- read_if_exists("data/ppi_edges.csv")
validation_table <- read_if_exists("data/Table_6_GSE107490_validation.csv")

# Completed exploratory TCGA-LAML outputs
tcga_survival <- read_if_exists("data/TCGA_LAML_stromal_score_survival_table.csv")
tcga_cox <- read_if_exists("data/TCGA_LAML_stromal_score_cox_results.csv")
tcga_mutations <- read_if_exists("data/TCGA_LAML_stromal_score_mutation_associations.csv")


####################################################################
# MODULE: 02_preprocess.R
####################################################################
suppressPackageStartupMessages({library(dplyr); library(tibble)})

harmonize_gene_symbols <- function(df, gene_col = "gene") {
  stopifnot(gene_col %in% names(df))
  df[[gene_col]] <- toupper(trimws(as.character(df[[gene_col]])))
  df %>% filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "")
}

# Use only when appropriate for the actual source expression scale.
log2_if_needed <- function(x, pseudocount = 1) {
  x <- as.matrix(x)
  if (max(x, na.rm = TRUE) > 50) log2(x + pseudocount) else x
}

prepare_expression <- function(expr_df, gene_col = "gene") {
  expr_df <- harmonize_gene_symbols(expr_df, gene_col)
  genes <- expr_df[[gene_col]]
  mat <- as.matrix(expr_df[, setdiff(names(expr_df), gene_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat <- log2_if_needed(mat)
  rownames(mat) <- genes
  mat
}


####################################################################
# MODULE: 03_DE_analysis.R
####################################################################
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
})

run_limma_DE <- function(expr, group, contrast_string, dataset_name) {
  group <- factor(group)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(expr, design)
  cont <- makeContrasts(contrasts = contrast_string, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, cont))

  out <- topTable(fit2, number = Inf, sort.by = "P") %>%
    rownames_to_column("gene") %>%
    transmute(
      dataset = dataset_name,
      contrast = contrast_string,
      gene = toupper(gene),
      log2FC = logFC,
      pvalue = P.Value,
      FDR = adj.P.Val
    )
  out
}

# Examples once expression matrices + metadata are present:
# g148 <- prepare_expression(gse148625_expr)
# de_aml_ctrl   <- run_limma_DE(g148, gse148625_meta$group, "AML-Control", "GSE148625")
# de_early_ctrl <- run_limma_DE(g148, gse148625_meta$group, "Early_AML-Control", "GSE148625")
# de_late_ctrl  <- run_limma_DE(g148, gse148625_meta$group, "Late_AML-Control", "GSE148625")
#
# g101 <- prepare_expression(gse101449_expr)
# de_arv_dmso <- run_limma_DE(g101, gse101449_meta$group, "ARV825-DMSO", "GSE101449")
#
# write_csv(bind_rows(de_aml_ctrl, de_early_ctrl, de_late_ctrl, de_arv_dmso),
#           "results/DEGs/deg_results.csv")


####################################################################
# MODULE: 04_PCA_analysis.R
####################################################################
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})

run_PCA <- function(expr, metadata, group_col, title, outfile) {
  pca <- prcomp(t(expr), scale. = TRUE)
  var <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  d <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2]
  ) %>% left_join(metadata, by = "sample")

  p <- ggplot(d, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
    geom_point(size = 3) +
    labs(title = title,
         x = sprintf("PC1 %.1f%%", var[1]),
         y = sprintf("PC2 %.1f%%", var[2]),
         color = NULL) +
    theme_classic(base_size = 12)

  ggsave(outfile, p, width = 6, height = 5, dpi = 600)
  list(pca = pca, scores = d, plot = p)
}


####################################################################
# MODULE: 05_directional_reversal.R
####################################################################
suppressPackageStartupMessages({library(dplyr); library(readr)})

directional_reversal <- function(aml_de, arv_de,
                                 aml_fdr = 0.05,
                                 arv_fdr = 0.05,
                                 min_abs_log2fc = 0) {
  a <- aml_de %>%
    transmute(gene = toupper(gene),
              AML_log2FC = log2FC,
              AML_FDR = FDR)

  b <- arv_de %>%
    transmute(gene = toupper(gene),
              ARV825_log2FC = log2FC,
              ARV825_FDR = FDR)

  inner_join(a, b, by = "gene") %>%
    mutate(
      direction_class = case_when(
        AML_log2FC > 0 & ARV825_log2FC < 0 ~ "AML-up / ARV-825-down",
        AML_log2FC < 0 & ARV825_log2FC > 0 ~ "AML-down / ARV-825-up",
        TRUE ~ "concordant/other"
      ),
      reversal_score = abs(AML_log2FC) + abs(ARV825_log2FC),
      relaxed_reversal = direction_class != "concordant/other",
      strict_reversal =
        relaxed_reversal &
        AML_FDR < aml_fdr & ARV825_FDR < arv_fdr &
        abs(AML_log2FC) >= min_abs_log2fc &
        abs(ARV825_log2FC) >= min_abs_log2fc
    ) %>%
    arrange(desc(reversal_score))
}

# IMPORTANT:
# Use the exact thresholds documented in the final manuscript.
# Do not change thresholds merely to reproduce a desired candidate count.


####################################################################
# MODULE: 06_GO_KEGG_enrichment.R
####################################################################
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

run_GO_KEGG <- function(symbols, prefix = "reversal") {
  conv <- bitr(unique(symbols),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
  entrez <- unique(conv$ENTREZID)

  go <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID", ont = "ALL",
                 pAdjustMethod = "BH", readable = TRUE)

  kegg <- enrichKEGG(gene = entrez, organism = "hsa",
                     pAdjustMethod = "BH")

  dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
  write_csv(as.data.frame(go),
            paste0("results/enrichment/", prefix, "_GO.csv"))
  write_csv(as.data.frame(kegg),
            paste0("results/enrichment/", prefix, "_KEGG.csv"))

  list(GO = go, KEGG = kegg)
}


####################################################################
# MODULE: 07_STRING_PPI.R
####################################################################
suppressPackageStartupMessages({
  library(STRINGdb)
  library(dplyr)
  library(readr)
  library(igraph)
})

run_STRING_PPI <- function(symbols, species = 9606, score_threshold = 400) {
  string_db <- STRINGdb$new(version = "11.5",
                            species = species,
                            score_threshold = score_threshold)

  genes <- data.frame(gene = unique(toupper(symbols)))
  mapped <- string_db$map(genes, "gene", removeUnmappedRows = TRUE)
  edges <- string_db$get_interactions(mapped$STRING_id)

  g <- graph_from_data_frame(edges[,c("from","to")], directed = FALSE)
  hubs <- data.frame(
    STRING_id = names(degree(g)),
    degree = as.numeric(degree(g))
  ) %>% arrange(desc(degree))

  dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
  write_csv(edges, "results/networks/STRING_edges.csv")
  write_csv(hubs, "results/networks/STRING_hubs.csv")
  list(mapped = mapped, edges = edges, hubs = hubs, graph = g)
}


####################################################################
# MODULE: 08_candidate_prioritization.R
####################################################################
suppressPackageStartupMessages({library(dplyr); library(readr)})

priority_genes <- c(
  "CXCL12","IGFBP3","PLOD2","VCAN","COL3A1",
  "SULF1","SOX9","BGLAP","NUDT1"
)

prioritize_candidates <- function(reversal_table, ppi_hubs = NULL) {
  x <- reversal_table %>%
    mutate(gene = toupper(gene),
           manuscript_candidate = gene %in% priority_genes)

  if (!is.null(ppi_hubs) && "gene" %in% names(ppi_hubs)) {
    x <- left_join(x, ppi_hubs, by = "gene")
  }
  x %>% arrange(desc(manuscript_candidate), desc(reversal_score))
}


####################################################################
# MODULE: 09_GSE107490_validation.R
####################################################################
suppressPackageStartupMessages({library(dplyr); library(readr)})

priority_genes <- c(
  "CXCL12","IGFBP3","PLOD2","VCAN","COL3A1",
  "SULF1","SOX9","BGLAP","NUDT1"
)

extract_validation <- function(de_table) {
  de_table %>%
    mutate(gene = toupper(gene)) %>%
    filter(gene %in% priority_genes) %>%
    arrange(match(gene, priority_genes))
}

# Manuscript interpretation is intentionally candidate-level and cautious:
# directional support was partial/context-dependent rather than universal.


####################################################################
# MODULE: 10_TCGA_LAML_analysis.R
####################################################################
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

d <- read_csv("data/TCGA_LAML_stromal_score_survival_table.csv",
              show_col_types = FALSE)

cat("TCGA-LAML rows:", nrow(d), "\n")
cat("Columns:\n"); print(names(d))

# The verified study output contains 132 patients, 80 deaths,
# 52 censored observations, and a median score cut point of -0.1430.
#
# Because archived output column names may vary, map these four variables
# explicitly before executing the models:
#
# time_col  <- "OS_time"
# event_col <- "OS_event"
# score_col <- "stromal_score"
# group_col <- "score_group"
#
# dd <- d %>% transmute(
#   time = .data[[time_col]],
#   event = .data[[event_col]],
#   stromal_score = .data[[score_col]],
#   score_group = factor(.data[[group_col]]),
#   age = age
# )
#
# km <- survfit(Surv(time, event) ~ score_group, data = dd)
# print(survdiff(Surv(time, event) ~ score_group, data = dd))
#
# cox_group <- coxph(Surv(time,event) ~ score_group, data=dd)
# cox_cont  <- coxph(Surv(time,event) ~ stromal_score, data=dd)
# cox_age   <- coxph(Surv(time,event) ~ stromal_score + age, data=dd)
#
# PH diagnostics should be retained when the exact generating analysis is rerun:
# print(cox.zph(cox_group))
# print(cox.zph(cox_cont))
# print(cox.zph(cox_age))
#
# IMPORTANT:
# This repository does not invent the upstream nine-gene score equation.
# Archive the exact original score-construction code when available.


####################################################################
# MODULE: 11_figure_generation.R
####################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

dir.create("figures", showWarnings = FALSE)

save_pub <- function(plot, filename, width = 12, height = 8) {
  ggsave(paste0("figures/", filename, ".png"),
         plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(paste0("figures/", filename, ".pdf"),
         plot, width = width, height = height, bg = "white")
}

# Figure 3 example: reversal candidates
rev_path <- "data/Table_3_reversal_candidates.csv"
if (file.exists(rev_path)) {
  rev <- read_csv(rev_path, show_col_types = FALSE)

  if (all(c("AML_log2FC","ARV825_log2FC") %in% names(rev))) {
    p3a <- ggplot(rev, aes(AML_log2FC, ARV825_log2FC)) +
      geom_hline(yintercept=0, linewidth=.3) +
      geom_vline(xintercept=0, linewidth=.3) +
      geom_point(alpha=.65) +
      labs(x="AML vs control log2FC",
           y="ARV-825 vs DMSO log2FC",
           title="Opposing AML and ARV-825 effects") +
      theme_classic(base_size=12)
    save_pub(p3a, "Figure_03_reversal_scatter", 7, 6)
  }
}

# Figure 4 example: enrichment
enr_path <- "data/Table_4_GO_KEGG_enrichment.csv"
if (file.exists(enr_path)) {
  enr <- read_csv(enr_path, show_col_types = FALSE)
  # Map exact archived column names here if needed.
}

# Figure 8 is generated from completed TCGA-LAML outputs only.
# BeatAML, scRNA-seq and ligand-receptor analyses are not generated here
# because they are not completed findings in the final manuscript.


####################################################################
# MODULE: 12_generate_all_figures.R
####################################################################
# Central figure-generation entry point.
# Each figure should be produced from archived processed tables.
source("R/11_figure_generation.R")

cat("Figure-generation workflow completed.\n")
cat("Expected final manuscript figures:\n")
cat("Figure 1  - study design / analytical workflow\n")
cat("Figure 2  - transcriptomic structure and differential expression\n")
cat("Figure 3  - directional reversal and candidate prioritization\n")
cat("Figure 4  - GO/KEGG functional enrichment\n")
cat("Figure 5  - STRING/PPI network\n")
cat("Figure 6  - integrated candidate-gene patterns\n")
cat("Figure 7  - GSE107490 human AML/MDS MSC validation\n")
cat("Figure 8  - exploratory TCGA-LAML clinical/molecular context\n")
cat("Graphical abstract - integrated study summary\n")

####################################################################
# MASTER EXECUTION NOTES
####################################################################

cat("\n============================================================\n")
cat("AML BRD4 stromal integrative analysis pipeline loaded.\n")
cat("Datasets: GSE148625, GSE101449, GSE107490, TCGA-LAML\n")
cat("============================================================\n")

cat("\nPrioritized manuscript genes:\n")
print(c("CXCL12","IGFBP3","PLOD2","VCAN","COL3A1",
        "SULF1","SOX9","BGLAP","NUDT1"))

cat("\nFinal manuscript figure map:\n")
cat("Figure 1: Study design and analytical workflow\n")
cat("Figure 2: Transcriptomic structure and differential expression\n")
cat("Figure 3: Directional reversal and candidate prioritization\n")
cat("Figure 4: GO/KEGG functional enrichment\n")
cat("Figure 5: STRING/PPI network analysis\n")
cat("Figure 6: Integrated candidate-gene patterns\n")
cat("Figure 7: GSE107490 human AML/MDS MSC validation\n")
cat("Figure 8: TCGA-LAML exploratory clinical/molecular context\n")
cat("Graphical abstract: Integrated study summary\n")

cat("\nReproducibility note:\n")
cat("Before public Zenodo deposition, archive the exact source/processed\n")
cat("inputs used for each reported analysis and the original TCGA score-\n")
cat("construction code when available. Proposed BeatAML, single-cell and\n")
cat("ligand-receptor analyses are not treated as completed results.\n")

cat("\nSession information:\n")
print(sessionInfo())
