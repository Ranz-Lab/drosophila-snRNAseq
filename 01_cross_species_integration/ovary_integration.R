## Hariyani et al., 2025
## R Script - Ovary Filtering, Normalization, Doublet Removal, Clustering, Integration

## Single-nucleus analysis of ovary - Drosophila melanogaster & Drosophila simulans
## ISO1, A4, w501

# Load libraries
# Must use Seurat v5, R version 4.3 or later

library(Seurat)
library(glmGamPoi)
library(dplyr)
library(Matrix)
library(ggplot2)
library(dittoSeq)
library(plotly)
library(SeuratWrappers)
library(clustree)
library(DoubletFinder)

# Function to create and filter Seurat object and perform QC
stringent_seurat <- function(strain_folder, ovary_samples) {
  mat <- readMM(paste0(strain_folder, "count_matrix.mtx"))
  cell_meta <- read.delim(paste0(strain_folder, "cell_metadata.csv"), stringsAsFactor = FALSE, sep = ",")
  genes <- read.delim(paste0(strain_folder, "all_genes.csv"), stringsAsFactor = FALSE, sep = ",")

  rownames(cell_meta) <- cell_meta$bc_wells
  genes$gene_name <- make.unique(genes$gene_name, sep = "_dup")
  colnames(mat) <- genes$gene_name
  rownames(mat) <- rownames(cell_meta)
  mat_t <- t(mat)
  mat_t <- mat_t[(rownames(mat_t) != ""),]

  strain <- CreateSeuratObject(mat_t, min.features = 200, min.cells = 3, meta.data = cell_meta)
  strain[["percent.mt"]] <- PercentageFeatureSet(strain, pattern = "^mt:")
  ovary <- subset(strain, subset = sample %in% ovary_samples)
  ovary <- subset(ovary, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & nCount_RNA <= 100000 & nFeature_RNA < 7500 & percent.mt <= 5)

  vln2 <- VlnPlot(ovary, features = c("nFeature_RNA"), group.by = "sample") + scale_y_continuous(limits = c(500,7500))
  vln4 <- VlnPlot(ovary, features = c("percent.mt"), group.by = "sample") + scale_y_continuous(limits = c(0,5))
  vln6 <- VlnPlot(ovary, features = c("nCount_RNA"), group.by = "sample") + scale_y_continuous(limits = c(0,20000))
  ggsave(paste0(strain_folder, "vln2_s.png"), plot = vln2, width = 7, height = 4)
  ggsave(paste0(strain_folder, "vln4_s.png"), plot = vln4, width = 7, height = 4)
  ggsave(paste0(strain_folder, "vln6_s.png"), plot = vln6, width = 7, height = 4)

  plot3 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
  plot4 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  ovary_corr <- CombinePlots(plots = list(plot3, plot4))
  ggsave(paste0(strain_folder, "ovary_corr_s.jpeg"), plot = ovary_corr, width = 10, height = 6)

  print("Ovary QC Summary:")
  print("Genes")
  print(summary(ovary$nFeature_RNA))
  print("Transcripts")
  print(summary(ovary$nCount_RNA))

  return(list(ovary = ovary))
}

# Sample definitions
A4_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/"
A4_ovary_samples <- c("A4_sample1o", "A4_sample2o")

ISO1_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/"
ISO1_ovary_samples <- c("ISO1_sample3o", "ISO1_sample4o")

w501_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/"
w501_ovary_samples <- c("w501_sample5o", "w501_sample6o")

# Create objects
A4_ovary <- stringent_seurat(A4_DGE_folder, A4_ovary_samples)$ovary
ISO1_ovary <- stringent_seurat(ISO1_DGE_folder, ISO1_ovary_samples)$ovary
w501_ovary <- stringent_seurat(w501_DGE_folder, w501_ovary_samples)$ovary

# Normalize and reduce
normalize_seurat <- function(strain, strain_folder) {
  strain <- NormalizeData(strain)
  strain <- FindVariableFeatures(strain)
  strain <- ScaleData(strain)
  strain <- RunPCA(strain)
  strain <- FindNeighbors(strain, dims = 1:15)
  strain <- FindClusters(strain, resolution = 1, cluster.name = "unintegrated_clusters")
  strain <- RunUMAP(strain, dims=1:15, reduction.name = "umap.unintegrated.normalized")
  ggsave(paste0(strain_folder, "umap.ovary.unintegrated.normalized.pdf"), plot = DimPlot(strain, reduction = "umap.unintegrated.normalized", group.by = "seurat_clusters"), width = 10, height = 10)
  return(strain)
}

A4_ovary <- normalize_seurat(A4_ovary, A4_DGE_folder)
ISO1_ovary <- normalize_seurat(ISO1_ovary, ISO1_DGE_folder)
w501_ovary <- normalize_seurat(w501_ovary, w501_DGE_folder)

saveRDS(A4_ovary, file = "./r_A4_ovary_normalized.rds")
saveRDS(ISO1_ovary, file = "./r_ISO1_ovary_normalized.rds")
saveRDS(w501_ovary, file = "./r_w501_ovary_normalized.rds")

# Doublet detection
identify_pk <- function(strain) {
  sweep.res.list <- paramSweep(strain, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  find.pK(sweep.stats)
}

doublet_estimate <- function(strain, strain_folder, pk) {
  annotations <- strain@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.03*nrow(strain@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  strain <- doubletFinder(strain, PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  relevant_columns <- grep("^DF.classifications_", colnames(strain@meta.data), value = TRUE)
  ggsave(paste0(strain_folder, "umap.ovary.doublets.pdf"), plot = DimPlot(strain, group.by = relevant_columns, reduction = "umap.unintegrated.normalized"), width = 10, height = 10)
  return(strain)
}

A4_pK <- identify_pk(A4_ovary)
A4_ovary <- doublet_estimate(A4_ovary, A4_DGE_folder, 0.26)
ISO1_pK <- identify_pk(ISO1_ovary)
ISO1_ovary <- doublet_estimate(ISO1_ovary, ISO1_DGE_folder, 0.23)
w501_pK <- identify_pk(w501_ovary)
w501_ovary <- doublet_estimate(w501_ovary, w501_DGE_folder, 0.25)

A4_ovary <- subset(A4_ovary, subset = DF.classifications_0.25_0.26_XXX == "Singlet")
ISO1_ovary <- subset(ISO1_ovary, subset = DF.classifications_0.25_0.23_XXX == "Singlet")
w501_ovary <- subset(w501_ovary, subset = DF.classifications_0.25_0.25_XXX == "Singlet")

saveRDS(A4_ovary, file = "./r_A4_ovary_NoDoublets.rds")
saveRDS(ISO1_ovary, file = "./r_ISO1_ovary_NoDoublets.rds")
saveRDS(w501_ovary, file = "./r_w501_ovary_NoDoublets.rds")

# Merge and integrate
ovary <- merge(x = ISO1_ovary, y = list(A4_ovary, w501_ovary))
ovary <- NormalizeData(ovary)
ovary <- FindVariableFeatures(ovary)
ovary <- ScaleData(ovary)
ovary <- RunPCA(ovary)

ovary <- IntegrateLayers(
  object = ovary, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony', assay = 'RNA',
  verbose = TRUE, normalization.method = "logNormalization")

ovary <- JoinLayers(ovary)

ElbowPlot(ovary, ndims=50)

ovary <- FindNeighbors(ovary, reduction = "harmony", dims = 1:15)
ovary <- FindClusters(ovary, resolution = 1.2, cluster.name = "harmony_clusters")
ovary <- RunUMAP(ovary, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
ovary <- RunTSNE(ovary, reduction = "harmony", dims = 1:20, tsne.method = "Rtsne", dim.embed = 2, reduction.name = "tsne.harmony")

saveRDS(ovary, file = "./ovary_normalized_nodoublets_integrated_unannotated_umap_tsne.rds")
