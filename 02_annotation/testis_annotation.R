## Hariyani et al., 2025
## 02_annotation - Testis annotation and visualization

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load integrated Seurat object
testis <- readRDS("testis_normalized_nodoublets_integrated_unannotated_umap_tsne.rds")

# Marker gene expression visualization
dotplot1 <- DotPlot(testis, features = c("His2Av","bam","aub","vas","nos","stg"), cols = c("blue","red")) #Spermatogonia
dotplot2 <- DotPlot(testis, features = c("bol","comr","topi","Rbp4","fzo","twe"), cols = c("blue","red")) #Spermatocytes
dotplot3 <- DotPlot(testis, features = c("Hsp23","Fas3","CadN", "Socs36E", "hh","rdo","eya","Egfr"), cols = c("blue","red")) #Somatic cells
dotplot4 <- DotPlot(testis, features = c("Dpy-30L2","sunz","whip","boly","soti","p-cup"), cols = c("blue","red")) #Spermatids

# Combine dotplots
combined_dotplot <- dotplot1 + dotplot2 + dotplot3 + dotplot4
print(combined_dotplot)

# Rename clusters based on marker expression
new.cluster.ids <- c(
  "Early Spermatocytes", "Early Spermatocytes", "Early Spermatocytes", "Early Spermatocytes",
  "Mid Spermatocytes", "Spermatids", "Mid Spermatocytes", "Late Spermatogonia", "Cyst Cells",
  "Late Spermatocytes", "Maturing Primary Spermatocytes", "Late Spermatogonia",
  "Early Spermatocytes", "Spermatids", "Epithelial Cells", "GSC / Early Spermatogonia",
  "Early Spermatocytes", "Mid Spermatocytes", "Early Spermatocytes",
  "Maturing Primary Spermatocytes", "Maturing Primary Spermatocytes", "Hub Cells",
  "Late Spermatogonia", "Maturing Primary Spermatocytes", "Unannotated"
)
names(new.cluster.ids) <- levels(testis)
testis <- RenameIdents(testis, new.cluster.ids)

# Define and apply preferred cluster order
my_levels <- c(
  "Epithelial Cells", "Hub Cells", "Cyst Cells", "GSC / Early Spermatogonia",
  "Late Spermatogonia", "Early Spermatocytes", "Mid Spermatocytes",
  "Late Spermatocytes", "Maturing Primary Spermatocytes", "Spermatids", "Unannotated"
)
Idents(testis) <- factor(Idents(testis), levels = my_levels)

# UMAP plot
umap_plot <- DimPlot(
  testis, reduction = "umap.harmony", label = TRUE,
  label.size = 2, repel = TRUE, label.box = TRUE
)
umap_plot + xlab("UMAP_1") + ylab("UMAP_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

# t-SNE plot
tsne_plot <- DimPlot(
  testis, reduction = "tsne.harmony", label = TRUE,
  label.size = 2, repel = TRUE, label.box = TRUE
)
tsne_plot + xlab("tSNE_1") + ylab("tSNE_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

# Save full object with annotations
saveRDS(testis, file = "./testis.rds")

# Remove unannotated clusters for downstream analyses
testis <- subset(testis, idents = "Unannotated", invert = TRUE)
saveRDS(testis, file = "./FINAL-testis-nounknown.rds")
