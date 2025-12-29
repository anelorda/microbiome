library(ggplot2)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(vegan)
library(dplyr)
library(scatterplot3d)
library(tidyr)
library(tibble)
library(scales)

# ============================================================
# DATA LOADING
# ============================================================

meta <- read.csv("soil/metadata_filtered.csv", row.names = 1)

meta$salinity <- factor(
  meta$salinity,
  levels = c("non_saline", "moderately_saline", "strongly_saline"),
  labels = c("NS", "MS", "HS")
)

meta_ids <- rownames(meta)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

prepare_abundance <- function(abund_raw, meta_ids) {
  abund_sub <- abund_raw[rownames(abund_raw) %in% meta_ids, , drop = FALSE]
  abund_sub[] <- lapply(abund_sub, as.numeric)
  abund_sub <- abund_sub[, colSums(abund_sub, na.rm = TRUE) > 0, drop = FALSE]
  abund_sub[meta_ids[meta_ids %in% rownames(abund_sub)], , drop = FALSE]
}

calc_relative_abundance <- function(abund_mat) {
  sweep(abund_mat, 1, rowSums(abund_mat), "/") * 100
}

plot_relative_abundance <- function(abund_mat, meta, level_name, top_n = 20, outfile) {
  
  rel_abund <- calc_relative_abundance(abund_mat)
  
  taxa_totals <- colSums(abund_mat)
  top_taxa <- names(sort(taxa_totals, decreasing = TRUE))[1:min(top_n, ncol(abund_mat))]
  other_taxa <- setdiff(colnames(rel_abund), top_taxa)
  
  rel_abund <- rel_abund[, top_taxa, drop = FALSE]
  if (length(other_taxa) > 0) {
    rel_abund$Other <- rowSums(calc_relative_abundance(abund_mat)[, other_taxa, drop = FALSE])
  }
  
  df <- as.data.frame(rel_abund) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "Taxon", values_to = "Abundance") %>%
    left_join(
      data.frame(SampleID = rownames(meta), Salinity = meta$salinity),
      by = "SampleID"
    )
  
  cols <- hue_pal(l = 65, c = 100)(length(unique(df$Taxon)))
  names(cols) <- unique(df$Taxon)
  if ("Other" %in% names(cols)) cols["Other"] <- "grey80"
  
  p <- ggplot(df, aes(SampleID, Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_wrap(~Salinity, scales = "free_x") +
    scale_fill_manual(values = cols) +
    theme_bw() +
    labs(
      title = paste("Relative abundance at", level_name, "level"),
      y = "Relative abundance (%)",
      x = "Sample"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.35, "cm"),
      legend.title = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  
  
  
  
  ggsave(outfile, p, width = 14, height = 7, dpi = 300)
}


plot_abundance_heatmap <- function(abund_mat, meta, level_name = "Taxa",
                                   salinity_col = "salinity",
                                   top_n = 50,
                                   filename = NULL) {
  
  # --- Select top N taxa by total abundance ---
  taxa_totals <- colSums(abund_mat, na.rm = TRUE)
  top_taxa <- names(sort(taxa_totals, decreasing = TRUE))[1:min(top_n, length(taxa_totals))]
  abund_top <- abund_mat[, top_taxa, drop = FALSE]
  
  cat("Plotting top", ncol(abund_top), level_name, "\n")
  
  # --- Log10 transform ---
  log_counts <- log10(abund_top + 1)
  
  # --- Order samples by salinity ---
  salinity_order <- c("NS", "MS", "HS")
  meta_ordered <- meta[order(factor(meta[[salinity_col]], levels = salinity_order)), , drop = FALSE]
  log_counts <- log_counts[rownames(meta_ordered), , drop = FALSE]
  
  # --- Annotation ---
  annotation_col <- data.frame(
    Salinity = factor(meta_ordered[[salinity_col]], levels = salinity_order),
    row.names = rownames(meta_ordered)
  )
  
  annotation_colors <- list(
    Salinity = c(NS = "#1f78b4", MS = "#33a02c", HS = "#e31a1c")
  )
  
  # --- Plot ---
  pheatmap(
    t(log_counts),                     # taxa = rows
    scale = "row",                     # Z-score per taxon
    clustering_distance_rows = "correlation",
    clustering_method = "complete",
    cluster_cols = FALSE,              # keep salinity order
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    main = paste0(level_name, " (top ", top_n, ") abundance heatmap (log10)"),
    filename = filename,
    width = 14,
    height = 10
  )
}



# ============================================================
# ALPHA + BETA DIVERSITY (SPECIES ONLY: level-7)
# ============================================================

abund_species <- prepare_abundance(
  read.csv("soil/level-7.csv", row.names = 1),
  meta_ids
)

# Alpha diversity
alpha_df <- data.frame(
  SampleID = rownames(abund_species),
  Chao1 = estimateR(abund_species)["S.chao1", ],
  Shannon = diversity(abund_species),
  Observed = specnumber(abund_species),
  Salinity = meta$salinity
)

ggsave(
  "soil/Fig1_Alpha_diversity_salinity.pdf",
  ggplot(alpha_df, aes(Salinity, Shannon, color = Salinity)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    theme_bw(),
  width = 8, height = 6
)

# Beta diversity
bray <- vegdist(abund_species, method = "bray")

pheatmap(
  as.matrix(bray),
  annotation_row = data.frame(Salinity = meta$salinity),
  annotation_col = data.frame(Salinity = meta$salinity),
  annotation_colors = list(Salinity = c(NS="#1f78b4",MS="#33a02c",HS="#e31a1c")),
  main = "Bray–Curtis dissimilarity",
  filename = "soil/Fig_BrayCurtis_heatmap.pdf",
  width = 10, height = 9
)

# ============================================================
# AUTOMATIC ANALYSIS FOR ALL LEVELS
# ============================================================

level_files <- list.files("soil", pattern = "^level-[0-9]+\\.csv$", full.names = TRUE)

level_names <- c(
  "1" = "Kingdom",
  "2" = "Phylum",
  "3" = "Class",
  "4" = "Order",
  "5" = "Family",
  "6" = "Genus",
  "7" = "Species"
)

for (f in level_files) {
  
  lvl <- gsub(".*level-([0-9]+)\\.csv", "\\1", f)
  lvl_name <- level_names[lvl]
  
  cat("Processing level", lvl, "-", lvl_name, "\n")
  
  abund <- prepare_abundance(read.csv(f, row.names = 1), meta_ids)
  
  plot_relative_abundance(
    abund, meta, lvl_name,
    outfile = paste0("soil/Fig_", lvl_name, "_relative_abundance.pdf")
  )
  
  plot_abundance_heatmap(
    abund, meta, lvl_name,
    filename = paste0("soil/Fig_", lvl_name, "_heatmap.pdf")
  )
}


cat("\n=== ALL LEVELS ANALYSIS COMPLETE ===\n")
