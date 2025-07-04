# Load libraries
library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(magrittr)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(DESeq2)
library(ggrepel)

# Load Seurat object
ovary <- readRDS("~/FINAL-ovary-nounknown.rds")

# Subset to only keep orthologs

# Import gene names for all 3 strains
A4_genes <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")
ISO1_genes <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")
w501_genes <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes_orthologs.csv")

# Identify genes with data for all three strains
common_gene_names <- Reduce(intersect, list(A4_genes$gene_name, ISO1_genes$gene_name, w501_genes$gene_name))

# Subset the Seurat object to retain only these genes
ovary <- subset(ovary, features=common_gene_names)

# Extract counts and metadata
counts <- GetAssayData(object = ovary, slot = "counts", assay = "RNA")
metadata <- ovary@meta.data
metadata$cluster_id <- factor(ovary@active.ident)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Filter out lowly expressed genes
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
  "A4sample1_ovary"   = "A4_sample1o",
  "A4sample6_ovary"   = "A4_sample2o",
  "ISO1sample2_ovary" = "ISO1_sample3o",
  "ISO1sample7_ovary" = "ISO1_sample4o",
  "w501sample3"       = "w501_sample5o",
  "w501sample8"       = "w501_sample6o"
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
  "A4_sample1o" = 16, "A4_sample2o" = 1,
  "ISO1_sample3o" = 17, "ISO1_sample4o" = 2,
  "w501_sample5o" = 15, "w501_sample6o" = 0
)
celltype_colors <- c(
  "Germline - GSC / Germarium 1 and 2a" = "#CC3555", 
  "Germline - Germarium 2a-2b" ="#C1472F",
  "Germline - Germarium 2b-3" ="#B85C00",
  "Stalk & Polar Cells" = "#A27511",
  "Follicle Stem Cells (FSC) / preFCs" ="#847D00",
  "Early Follicle Cells" = "#698401",
  "Mitotic Follicle Cells Stage 1-5" = "#358A0F",
  "Post-Mitotic Follicle Cells Stage 6" = "#178D4D",
  "Vitellogenic MBFCs Stage 7" = "#009466",
  "Vitellogenic MBFCs Stage 8" = "#009684",
  "Vitellogenic MBFCs Stage 9-10A" = "#00929E",
  "Choriogenic MBFCs Stage 12" = "#0090B3",
  "Choriogenic MBFCs Stage 14" = "#0088C7",
  "Stretch Cells" = "#8044CB",
  "Terminal Corpus Luteum Cells" = "#3E5ACB",
  "Ovarian Sheath Muscle" = "#AA3DB8",
  "Oviduct" = "#C032A8"
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

ggsave("Ovary_Pseudobulk_PCA_CellType.png", plot = p, width = 6, height = 6, dpi = 300)
print(p)
