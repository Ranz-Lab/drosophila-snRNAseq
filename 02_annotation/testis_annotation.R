## Hariyani et al., 2025
## 02_annotation - Testis annotation, visualization and trajectory inference

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(monocle3)

# Load integrated Seurat object
testis <- readRDS("testis_normalized_nodoublets_integrated_unannotated_umap_tsne.rds")

# ---------------------------------------
# Marker gene visualization (DotPlots)
# ---------------------------------------

# Marker gene expression visualization
#Spermatogonia
dotplot1 <- DotPlot(testis, features = c("His2Av","bam","aub","vas","nos","stg"), split.by="strain",cols=c("#C83658","#3A8A00","#0088C5")) 

#Spermatocytes
dotplot2 <- DotPlot(testis, features = c( "bol","comr","topi","Rbp4","fzo","twe"),cols=c("#C83658","#3A8A00","#0088C5"), split.by="strain") 

#Somatic cells
dotplot3 <- DotPlot(testis, features = c("Hsp23","Fas3","CadN", "Socs36E", "hh","rdo","eya","Egfr"),cols=c("#C83658","#3A8A00","#0088C5"), split.by="strain") 

#Spermatids
dotplot4 <- DotPlot(testis, features = c("Dpy-30L2","sunz","whip","boly","soti","p-cup"),cols=c("#C83658","#3A8A00","#0088C5"), split.by="strain") 

# Combine dotplots
combined_dotplot <- dotplot1 + dotplot2 + dotplot3 + dotplot4
print(combined_dotplot)

# ---------------------------------------
# Manual annotation
# ---------------------------------------

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

# ---------------------------------------
# Dimensionality reduction plots
# ---------------------------------------

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

# ---------------------------------------
# Trajectory analysis - Monocle3
# ---------------------------------------

# Set UMAP for Monocle3
testis[["UMAP"]] <- testis[["umap.harmony"]]

# Convert to Monocle3 object
cds <- as.cell_data_set(testis)
cds <- cluster_cells(cds, resolution = 1e-5, k = 10)

# Plot pre-trajectory cluster and partition structure
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE) +
  plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

# Assign cluster identities from Seurat
cds@clusters@listData[["UMAP"]][["clusters"]] <- testis@active.ident

# Assign UMAP coordinates
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- testis@reductions$umap@cell.embeddings

# Visualize clusters before trajectory
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = FALSE,
           group_label_size = 5) + theme(legend.position = "right")

# Learn trajectory graph
cds <- learn_graph(cds, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 20))

# Plot trajectory
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = FALSE,
           label_branch_points = TRUE, label_roots = TRUE, label_leaves = FALSE,
           group_label_size = 5)

# Order cells in pseudotime
cds <- order_cells(cds)

# Plot pseudotime
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE,
           label_branch_points = FALSE, label_roots = FALSE, label_leaves = FALSE,
           trajectory_graph_color = "red")

# Save pseudotime
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Final pseudotime boxplot
ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = celltype, fill = celltype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12)
  ) +
  labs(x = "Pseudotime")
