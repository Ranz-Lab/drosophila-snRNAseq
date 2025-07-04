# Testis Gene Expression Correlation Analysis
# Hariyani et al., 2025

library(Seurat)
library(ggplot2)
library(tidyr)
library(forcats)
library(ggcorrplot)

# Load testis Seurat object
testis <- readRDS("~/FINAL-testis-nounknown.rds")

# 1A. Correlation matrix across cell types (combined samples)
av.exp <- AverageExpression(testis)$RNA
cor.exp <- cor(av.exp)
cor.df <- as.data.frame(cor.exp)
cor.df$x <- rownames(cor.df)

desired_order <- c(
  "Epithelial Cells", "Cyst Cells", "Hub Cells", 
  "GSC / Early Spermatogonia", "Late Spermatogonia", 
  "Early Spermatocytes", "Mid Spermatocytes", 
  "Late Spermatocytes", "Maturing Primary Spermatocytes", "Spermatids"
)

cor_long <- pivot_longer(cor.df, -x, names_to = "y", values_to = "correlation")
ggplot(cor_long, aes(factor(x, levels = desired_order), factor(y, levels = desired_order), fill = correlation)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Testis: Cell Type Correlation Heatmap", x = NULL, y = NULL)

# 1B. Per-strain correlation matrices
get_correlation_plot <- function(seurat_obj, title_text) {
  avg_exp <- AverageExpression(seurat_obj)$RNA
  cor_mat <- cor(avg_exp)
  cor_df <- as.data.frame(cor_mat)
  cor_df$x <- rownames(cor_df)

  cor_long <- pivot_longer(cor_df, -x, names_to = "y", values_to = "correlation")
  
  ggplot(cor_long, aes(factor(x, levels = desired_order), factor(y, levels = desired_order), fill = correlation)) +
    geom_tile() +
    scale_fill_gradient(low = "lightblue2", high = "maroon4") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

ISO1 <- subset(testis, subset = strain == "ISO1")
A4   <- subset(testis, subset = strain == "A4")
w501 <- subset(testis, subset = strain == "w501")

get_correlation_plot(ISO1, "ISO1: Cell Type Correlation")
get_correlation_plot(A4, "A4: Cell Type Correlation")
get_correlation_plot(w501, "w501: Cell Type Correlation")

# 1C. Cross-strain expression correlation matrix (filtered orthologs)
genes_ISO1 <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")$gene_name
genes_A4   <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")$gene_name
genes_w501 <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")$gene_name

common_genes <- Reduce(intersect, list(genes_ISO1, genes_A4, genes_w501))
testis_orthologs <- subset(testis, features = common_genes)

av.exp <- AverageExpression(testis_orthologs, add.ident = "strain")$RNA

# Optional outlier filtering
upper_bounds <- apply(av.exp, 2, quantile, 0.95, na.rm = TRUE)

filtered_cor <- function(x, y, x_thresh, y_thresh) {
  idx <- (x >= 0.01 & x <= x_thresh) & (y >= 0.01 & y <= y_thresh)
  if (sum(idx) > 1) cor(x[idx], y[idx], method = "pearson") else NA
}

# Correlation matrix after outlier removal
cor_mat <- matrix(NA, ncol = ncol(av.exp), nrow = ncol(av.exp))
colnames(cor_mat) <- rownames(cor_mat) <- colnames(av.exp)

for (i in 1:ncol(av.exp)) {
  for (j in 1:ncol(av.exp)) {
    cor_mat[i, j] <- filtered_cor(av.exp[, i], av.exp[, j], upper_bounds[i], upper_bounds[j])
  }
}

# Rename & reorder
cor_mat <- as.data.frame(cor_mat)
cor_mat <- cor_mat[order(sub(".*_", "", rownames(cor_mat))), order(sub(".*_", "", colnames(cor_mat)))]
rownames(cor_mat) <- sub("-.*", "", rownames(cor_mat))
colnames(cor_mat) <- sub("-.*", "", colnames(cor_mat))

ggcorrplot(cor_mat, lab = FALSE) +
  scale_fill_gradient(low = "lightblue", high = "maroon4", limits = c(0.1, 1))

# Gene count table for each correlation
gene_counts <- matrix(NA, ncol = ncol(av.exp), nrow = ncol(av.exp))
colnames(gene_counts) <- rownames(gene_counts) <- colnames(av.exp)

for (i in 1:ncol(av.exp)) {
  for (j in 1:ncol(av.exp)) {
    x <- av.exp[, i]
    y <- av.exp[, j]
    valid_idx <- (x >= 0.01 & x <= upper_bounds[i]) & (y >= 0.01 & y <= upper_bounds[j])
    gene_counts[i, j] <- sum(valid_idx, na.rm = TRUE)
  }
}

gene_counts_df <- as.data.frame(gene_counts)
rownames(gene_counts_df) <- rownames(cor_mat)
colnames(gene_counts_df) <- colnames(cor_mat)
write.csv(gene_counts_df, "testis_expression_correlation-gene_count_table.csv")
