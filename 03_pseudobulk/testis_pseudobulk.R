# Load libraries
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(magrittr)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load testis Seurat object
testis <- readRDS("~/FINAL-testis-nounknown.rds")

# Subset to only keep orthologs
A4_genes   <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")
ISO1_genes <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")
w501_genes <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")
common_gene_names <- Reduce(intersect, list(A4_genes$gene_name, ISO1_genes$gene_name, w501_genes$gene_name))
testis <- subset(testis, features = common_gene_names)

# Extract counts and metadata
counts <- GetAssayData(object = testis, slot = "counts", assay = "RNA")
metadata <- testis@meta.data
metadata$cluster_id <- factor(testis@active.ident)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Filter lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]

# Aggregate counts per cluster-sample pair
groups <- colData(sce)[, c("cluster_id", "sample")]
groups$sample <- factor(groups$sample)
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum")

# Split by cluster
splitf <- sapply(stringr::str_split(rownames(pb), "_", n = 2), `[`, 1)
pb <- split.data.frame(pb, factor(splitf)) %>%
  lapply(function(u) set_colnames(t(u), gsub("^[^_]+_", "", rownames(u))))

# Rename sample columns
sample_rename_map <- c(
  "A4sample4_testis"    = "A4_sample1t",
  "A4sample10_testis"   = "A4_sample2t",
  "ISO1sample5_testis"  = "ISO1_sample3t",
  "ISO1sample9_testis"  = "ISO1_sample4t",
  "w501sample11"        = "w501_sample5t",
  "w501sample12"        = "w501_sample6t"
)
pb <- lapply(pb, function(mat) {
  colnames(mat) <- sample_rename_map[colnames(mat)]
  mat
})

# Merge matrices into single pseudobulk matrix
pb_list <- lapply(names(pb), function(ct) {
  mat <- pb[[ct]]
  colnames(mat) <- paste(ct, colnames(mat), sep = "__")
  mat
})
combined_pb <- do.call(cbind, pb_list)

# Create metadata
metadata_combined <- data.frame(
  cell_type = sapply(strsplit(colnames(combined_pb), "__"), `[`, 1),
  sample    = sapply(strsplit(colnames(combined_pb), "__"), `[`, 2)
)
rownames(metadata_combined) <- colnames(combined_pb)

# Normalize via rlog
dds <- DESeqDataSetFromMatrix(combined_pb, colData = metadata_combined, design = ~1)
dds <- dds[rowSums(counts(dds)) > 0, ]
vsd <- rlog(dds)

# PCA
pca <- prcomp(t(assay(vsd)))
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  cell_type = metadata_combined$cell_type,
  sample = metadata_combined$sample
)

# Define shape and color mappings
sample_shapes <- c(
  "A4_sample1t"    = 16, "A4_sample2t"    = 1,
  "ISO1_sample3t"  = 17, "ISO1_sample4t"  = 2,
  "w501_sample5t"  = 15, "w501_sample6t"  = 0
)
celltype_colors <- c(
  "Epithelial Cells"               = "#C83658", 
  "Hub Cells"                      = "#B85617",
  "Cyst Cells"                     = "#997900", 
  "GSC / Early Spermatogonia"     = "#6D8300",  
  "Late Spermatogonia"            = "#018E1E", 
  "Early Spermatocytes"           = "#00946C",  
  "Mid Spermatocytes"             = "#019399",  
  "Late Spermatocytes"            = "#008DBD",
  "Maturing Primary Spermatocytes"= "#4458CD",
  "Spermatids"                     = "#A03DC3"
)

# Plot PCA
pca_df$cell_type <- factor(pca_df$cell_type, levels = names(celltype_colors))
pca_df <- na.omit(pca_df)
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cell_type, shape = sample)) +
  geom_point(size = 4) +
  scale_color_manual(values = celltype_colors) +
  scale_shape_manual(values = sample_shapes) +
  labs(x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

ggsave("Testis_Pseudobulk_PCA_CellType.png", plot = p, width = 6, height = 6, dpi = 300)
print(p)
