# ============================================================
# LIBRARIES
# ============================================================

library(ggplot2)
library(tidyverse)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(scatterplot3d)

# ============================================================
# USER SETTINGS
# ============================================================

remove_samples <- c("MetaG-Soil-6", "MetaG-Soil-7")

# ============================================================
# DATA LOADING
# ============================================================

# ---------- Metadata ----------
meta <- read.csv(
  "soil/metadata_filtered.csv",
  row.names = 1,
  stringsAsFactors = FALSE
)

# Remove unwanted samples from metadata
meta <- meta[!rownames(meta) %in% remove_samples, ]

# ---------- OTU table ----------
otu_raw <- read.delim(
  "soil/Soil-1-47/otu/otu_table.tsv",
  sep = "\t",
  comment.char = "#",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# First column = OTU IDs
rownames(otu_raw) <- otu_raw[, 1]
otu_raw <- otu_raw[, -1]

# Convert to numeric
otu_raw[] <- lapply(otu_raw, as.numeric)

# Transpose: samples as rows
otu_raw <- t(otu_raw)

# Remove unwanted samples from OTU table
otu_raw <- otu_raw[!rownames(otu_raw) %in% remove_samples, , drop = FALSE]

# ============================================================
# MATCH OTU TABLE AND METADATA
# ============================================================

# Keep only shared samples
common_samples <- intersect(rownames(meta), rownames(otu_raw))

meta <- meta[common_samples, , drop = FALSE]
otu_raw <- otu_raw[common_samples, , drop = FALSE]

# Ensure identical order
stopifnot(all(rownames(meta) == rownames(otu_raw)))

# ============================================================
# METADATA PREPARATION
# ============================================================

meta$salinity <- factor(
  meta$salinity,
  levels = c("non_saline", "moderately_saline", "strongly_saline"),
  labels = c("NS", "MS", "HS")
)

# ============================================================
# OTU ABUNDANCE MATRIX
# ============================================================

# Remove zero-abundance OTUs
abund_otu <- otu_raw[, colSums(otu_raw) > 0, drop = FALSE]

# Sanity checks
stopifnot(nrow(abund_otu) > 0, ncol(abund_otu) > 0)

# ============================================================
# RAREFACTION (JOURNAL-RECOMMENDED)
# ============================================================

set.seed(123)
abund_otu <- rrarefy(abund_otu, min(rowSums(abund_otu)))

# ============================================================
# ALPHA DIVERSITY (OTU-BASED)
# ============================================================

cat("\n=== ALPHA DIVERSITY (OTU LEVEL) ===\n")

Observed <- specnumber(abund_otu)
Shannon  <- diversity(abund_otu, index = "shannon")
Chao1    <- estimateR(abund_otu)["S.chao1", ]
Goods    <- apply(abund_otu, 1, function(x) 1 - sum(x == 1) / sum(x))

alpha_df <- data.frame(
  SampleID = rownames(abund_otu),
  Chao1 = Chao1,
  Shannon = Shannon,
  Observed = Observed,
  Goods = Goods,
  Salinity = meta$salinity
)

alpha_long <- alpha_df %>%
  pivot_longer(
    cols = c("Chao1", "Shannon", "Observed", "Goods"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  mutate(
    Index = factor(
      Index,
      levels = c("Chao1", "Shannon", "Observed", "Goods"),
      labels = c(
        "(a) Chao1 richness",
        "(b) Shannon's diversity",
        "(c) Observed OTUs",
        "(d) Good's coverage"
      )
    )
  )

p_alpha <- ggplot(alpha_long, aes(x = Salinity, y = Value, color = Salinity)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~Index, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(x = "", y = NULL) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

print(p_alpha)

ggsave(
  "soil/Fig1_Alpha_diversity_OTU_no_Soil6_7.pdf",
  p_alpha,
  width = 10,
  height = 7,
  dpi = 300
)

# Statistics
cat("\nKruskal-Wallis tests:\n")
print(kruskal.test(Chao1 ~ Salinity, data = alpha_df))
print(kruskal.test(Shannon ~ Salinity, data = alpha_df))
print(kruskal.test(Observed ~ Salinity, data = alpha_df))
print(kruskal.test(Goods ~ Salinity, data = alpha_df))

# ============================================================
# BETA DIVERSITY (OTU-BASED)
# ============================================================

cat("\n=== BETA DIVERSITY (OTU LEVEL) ===\n")

bray_dist <- vegdist(abund_otu, method = "bray")

# ---------- ANOSIM ----------
anosim_res <- anosim(bray_dist, grouping = meta$salinity, permutations = 999)
anosim_R <- round(anosim_res$statistic, 3)
anosim_p <- anosim_res$signif

print(anosim_res)

pdf("soil/Fig2a_ANOSIM_OTU_no_Soil6_7.pdf", width = 7, height = 6)
plot(anosim_res, main = "ANOSIM based on Bray–Curtis distance (OTU)")
dev.off()

# ============================================================
# PCoA
# ============================================================

pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 3)

pcoa_df <- data.frame(
  SampleID = rownames(meta),
  PC1 = pcoa_res$points[,1],
  PC2 = pcoa_res$points[,2],
  PC3 = pcoa_res$points[,3],
  Salinity = meta$salinity
)

var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)

# ---------- 2D PCoA ----------
p_pcoa2d <- ggplot(pcoa_df, aes(PC1, PC2, color = Salinity)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_bw() +
  labs(
    x = paste0("PCoA1 (", var_exp[1], "%)"),
    y = paste0("PCoA2 (", var_exp[2], "%)"),
    title = "(b) PCoA (Bray–Curtis, OTU)"
  ) +
  annotate(
    "text",
    x = min(pcoa_df$PC1),
    y = max(pcoa_df$PC2),
    hjust = 0,
    vjust = 1,
    size = 4,
    label = paste0("ANOSIM R = ", anosim_R, "\nP = ", anosim_p)
  )

print(p_pcoa2d)

ggsave(
  "soil/Fig2b_PCoA_2D_OTU_no_Soil6_7.pdf",
  p_pcoa2d,
  width = 7,
  height = 6
)

# ============================================================
# 3D PCoA
# ============================================================

cols <- c("NS" = "forestgreen", "MS" = "orange", "HS" = "red")

pdf("soil/Fig2c_PCoA_3D_OTU_no_Soil6_7.pdf", width = 7, height = 7)
scatterplot3d(
  pcoa_df$PC1,
  pcoa_df$PC2,
  pcoa_df$PC3,
  color = cols[pcoa_df$Salinity],
  pch = 19,
  xlab = paste0("PCoA1 (", var_exp[1], "%)"),
  ylab = paste0("PCoA2 (", var_exp[2], "%)"),
  zlab = paste0("PCoA3 (", var_exp[3], "%)"),
  main = "(c) PCoA 3D (Bray–Curtis, OTU)"
)
dev.off()
