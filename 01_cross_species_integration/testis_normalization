## Hariyani et al., 2025
## R Script - Testis Filtering, Normalization, Doublet Removal, Clustering, Integration

## Single-nucleus analysis of testis - Drosophila melanogaster & Drosophila simulans
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

# Function to create and filter stringent Seurat object and perform QC
relaxed_seurat <- function(strain_folder, testis_samples) {
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
  testis <- subset(strain, subset = sample %in% testis_samples)
  testis <- subset(testis, subset = nFeature_RNA >= 200 & nCount_RNA >= 300 & nCount_RNA <= 100000 & nFeature_RNA < 7500 & percent.mt <= 5)

  vln2 <- VlnPlot(testis, features = c("nFeature_RNA"), group.by = "sample") + scale_y_continuous(limits = c(500,7500))
  vln4 <- VlnPlot(testis, features = c("percent.mt"), group.by = "sample") + scale_y_continuous(limits = c(0,5))
  vln6 <- VlnPlot(testis, features = c("nCount_RNA"), group.by = "sample") + scale_y_continuous(limits = c(0,20000))
  ggsave(paste0(strain_folder, "vln2_s.png"), plot = vln2, width = 7, height = 4)
  ggsave(paste0(strain_folder, "vln4_s.png"), plot = vln4, width = 7, height = 4)
  ggsave(paste0(strain_folder, "vln6_s.png"), plot = vln6, width = 7, height = 4)

  plot3 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample")
  plot4 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  testis_corr <- CombinePlots(plots = list(plot3, plot4))
  ggsave(paste0(strain_folder, "testis_corr_s.jpeg"), plot = testis_corr, width = 10, height = 6)

  print("Testis QC Summary:")
  print("Genes")
  print(summary(testis$nFeature_RNA))
  print("Transcripts")
  print(summary(testis$nCount_RNA))

  return(list(testis = testis))
}

# Define strain paths and samples
A4_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/"
A4_testis_samples <- c("A4_sample1t", "A4_sample2t")
ISO1_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/"
ISO1_testis_samples <- c("ISO1_sample3t", "ISO1_sample4t")
w501_DGE_folder <- "/Users/Imtiyaz/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/"
w501_testis_samples <- c("w501_sample5t", "w501_sample6t")

# Create Seurat objects
A4_testis <- relaxed_seurat(A4_DGE_folder, A4_testis_samples)$testis
ISO1_testis <- relaxed_seurat(ISO1_DGE_folder, ISO1_testis_samples)$testis
w501_testis <- relaxed_seurat(w501_DGE_folder, w501_testis_samples)$testis

# Normalize and reduce
normalize_seurat <- function(strain, strain_folder) {
  strain <- NormalizeData(strain)
  strain <- FindVariableFeatures(strain)
  strain <- ScaleData(strain)
  strain <- RunPCA(strain)
  strain <- FindNeighbors(strain, dims = 1:15)
  strain <- FindClusters(strain, resolution = 1, cluster.name = "unintegrated_clusters")
  strain <- RunUMAP(strain, dims=1:15, reduction.name = "umap.unintegrated.normalized")
  ggsave(paste0(strain_folder, "umap.testis.unintegrated.normalized.pdf"), plot = DimPlot(strain, reduction = "umap.unintegrated.normalized", group.by = "seurat_clusters"), width = 10, height = 10)
  return(strain)
}

A4_testis <- normalize_seurat(A4_testis, A4_DGE_folder)
ISO1_testis <- normalize_seurat(ISO1_testis, ISO1_DGE_folder)
w501_testis <- normalize_seurat(w501_testis, w501_DGE_folder)

saveRDS(A4_testis, file = "./r_A4_testis_normalized.rds")
saveRDS(ISO1_testis, file = "./r_ISO1_testis_normalized.rds")
saveRDS(w501_testis, file = "./r_w501_testis_normalized.rds")

# Doublet Detection
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
  ggsave(paste0(strain_folder, "umap.testis.doublets.pdf"), plot = DimPlot(strain, group.by = relevant_columns, reduction = "umap.unintegrated.normalized"), width = 10, height = 10)
  return(strain)
}

A4_pK <- identify_pk(A4_testis)
A4_testis <- doublet_estimate(A4_testis, A4_DGE_folder, 0.26)
ISO1_pK <- identify_pk(ISO1_testis)
ISO1_testis <- doublet_estimate(ISO1_testis, ISO1_DGE_folder, 0.23)
w501_pK <- identify_pk(w501_testis)
w501_testis <- doublet_estimate(w501_testis, w501_DGE_folder, 0.25)

A4_testis <- subset(A4_testis, subset = DF.classifications_0.25_0.26_566 == "Singlet")
ISO1_testis <- subset(ISO1_testis, subset = DF.classifications_0.25_0.23_371 == "Singlet")
w501_testis <- subset(w501_testis, subset = DF.classifications_0.25_0.25_340 == "Singlet")

saveRDS(A4_testis, file = "./r_A4_testis_NoDoublets.rds")
saveRDS(ISO1_testis, file = "./r_ISO1_testis_NoDoublets.rds")
saveRDS(w501_testis, file = "./r_w501_testis_NoDoublets.rds")

# Final merge and Harmony integration
A4_testis <- NormalizeData(A4_testis)
ISO1_testis <- NormalizeData(ISO1_testis)
w501_testis <- NormalizeData(w501_testis)

testis <- merge(x = ISO1_testis, y = list(A4_testis, w501_testis))
testis <- NormalizeData(testis)
testis <- FindVariableFeatures(testis)
testis <- ScaleData(testis)
testis <- RunPCA(testis)

testis <- IntegrateLayers(
  object = testis, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = 'harmony', assay = 'RNA',
  verbose = TRUE, normalization.method = "logNormalization")

testis <- JoinLayers(testis)

ElbowPlot(testis, ndims=50)

testis <- FindNeighbors(testis, reduction = "harmony", dims = 1:15)
testis <- FindClusters(testis, resolution = 1.2, cluster.name = "harmony_clusters")
testis <- RunUMAP(testis, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
testis <- RunTSNE(testis, reduction = "harmony", dims = 1:20, tsne.method = "Rtsne", dim.embed = 2, reduction.name = "tsne.harmony")

saveRDS(testis, file = "./testis_normalized_nodoublets_integrated_unannotated_umap_tsne.rds")
