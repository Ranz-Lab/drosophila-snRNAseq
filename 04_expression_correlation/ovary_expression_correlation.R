# Ovary Gene Expression Correlation Analysis
# Hariyani et al., 2025

library(Seurat)
library(tidyr)
library(ggplot2)
library(ggcorrplot)

# Load ovary Seurat object
ovary <- readRDS("~/FINAL-ovary-nounknown.rds")

# 1A. Correlation matrix of gene expression (all ovary samples combined)
av.exp <- AverageExpression(ovary)$RNA
cor.exp <- cor(as.matrix(av.exp))

desired_order <- c(
  "Germline - GSC / Germarium 1 and 2a", "Germline - Germarium 2a-2b",
  "Germline - Germarium 2b-3", "Stalk & Polar Cells", 
  "Follicle Stem Cells (FSC) / preFCs", "Early Follicle Cells", 
  "Mitotic Follicle Cells Stage 1-5", "Post-Mitotic Follicle Cells Stage 6",
  "Vitellogenic MBFCs Stage 7", "Vitellogenic MBFCs Stage 8",
  "Vitellogenic MBFCs Stage 9-10A",
  "Choriogenic MBFCs Stage 12", "Choriogenic MBFCs Stage 14",
  "Terminal Corpus Luteum Cells", "Stretch Cells", "Ovarian Sheath Muscle", "Oviduct"
)

cor.df <- tidyr::gather(as.data.frame(cor.exp), y, correlation)
cor.df$x <- rep(rownames(cor.exp), ncol(cor.exp))

ggplot(cor.df, aes(factor(x, levels = desired_order), factor(y, levels = desired_order), fill = correlation)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Ovary Cell-Type Correlation Heatmap", x = "Cell Type", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# 1B. Per-strain correlation heatmaps
A4   <- subset(ovary, subset = strain == "A4")
ISO1 <- subset(ovary, subset = strain == "ISO1")
w501 <- subset(ovary, subset = strain == "w501")

corr_mat <- function(seurat_obj, title_suffix) {
  av.exp <- AverageExpression(seurat_obj)$RNA
  cor.exp <- cor(as.matrix(av.exp))
  
  cor.df <- tidyr::gather(as.data.frame(cor.exp), y, correlation)
  cor.df$x <- rep(rownames(cor.exp), ncol(cor.exp))
  
  ggplot(cor.df, aes(factor(x, levels = desired_order), factor(y, levels = desired_order), fill = correlation)) +
    geom_tile() +
    scale_fill_gradient(low = "lightblue2", high = "maroon4") +
    labs(title = paste("Ovary Correlation Heatmap -", title_suffix), x = "Cell Type", y = "Cell Type") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

corr_mat(ISO1, "ISO1")
corr_mat(A4, "A4")
corr_mat(w501, "w501")

# 1C. Cross-strain correlation after filtering orthologs
gene_lists <- list(
  A4   = read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")$gene_name,
  ISO1 = read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")$gene_name,
  w501 = read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")$gene_name
)
common_genes <- Reduce(intersect, gene_lists)

ovary_orthologs <- subset(ovary, features = common_genes)
av.exp <- AverageExpression(ovary_orthologs, add.ident = "strain")$RNA

# Compute 95th percentile bounds per cell type
upper_bounds <- apply(av.exp, 2, function(x) quantile(x, 0.95, na.rm = TRUE))

# Filtered Pearson correlation
filtered_cor <- function(x, y, x_max, y_max) {
  idx <- (x >= 0.01 & x <= x_max) & (y >= 0.01 & y <= y_max)
  if (sum(idx) > 1) cor(x[idx], y[idx], method = "pearson") else NA
}

# Compute pairwise correlation
cor.exp <- matrix(NA, ncol = ncol(av.exp), nrow = ncol(av.exp))
colnames(cor.exp) <- colnames(av.exp)
rownames(cor.exp) <- colnames(av.exp)

for (i in seq_len(ncol(av.exp))) {
  for (j in seq_len(ncol(av.exp))) {
    cor.exp[i, j] <- filtered_cor(av.exp[, i], av.exp[, j], upper_bounds[i], upper_bounds[j])
  }
}

# Remove TCLC, Stretch Cells
remove_rows <- grep("^(Terminal Corpus Luteum Cells|Stretch Cells)", rownames(cor.exp))
if (length(remove_rows) > 0) cor.exp <- cor.exp[-remove_rows, -remove_rows]

# Final heatmap
ggcorrplot(cor.exp, lab = FALSE) +
  scale_fill_gradient(low = "lightblue", high = "maroon4", limits = c(0.75, 1))

# Count number of genes used per pairwise comparison
gene_count_matrix <- matrix(NA, ncol = ncol(av.exp), nrow = ncol(av.exp))
colnames(gene_count_matrix) <- colnames(av.exp)
rownames(gene_count_matrix) <- colnames(av.exp)

for (i in seq_len(ncol(av.exp))) {
  for (j in seq_len(ncol(av.exp))) {
    x <- av.exp[, i]; y <- av.exp[, j]
    valid <- (x >= 0.01 & x <= upper_bounds[i]) & (y >= 0.01 & y <= upper_bounds[j])
    gene_count_matrix[i, j] <- sum(valid, na.rm = TRUE)
  }
}

gene_count_df <- as.data.frame(gene_count_matrix)
if (length(remove_rows) > 0) gene_count_df <- gene_count_df[-remove_rows, -remove_rows]

write.csv(gene_count_df, "ovary_expression_correlation-gene_count_table.csv", row.names = TRUE)
