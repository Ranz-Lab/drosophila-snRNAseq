## Hariyani et al., 2025
## 02_annotation - Ovary annotation and visualization

# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load integrated Seurat object
ovary <- readRDS("ovary_normalized_nodoublets_integrated_unannotated_umap_tsne.rds")

# Marker gene visualization

# Germline Cells; ~8-cell stage
dotplot1 <- DotPlot(ovary, features = c("corolla", "bam", "osk", "orb", "gro"),
        split.by = "strain", cols = c("#C83658", "#3A8A00", "#0088C5"))

# Stalk & Polar Cells; FSC and preFCs
dotplot2 <- DotPlot(ovary, features = c("drl", "cas", "upd1", "upd3", "Wnt4", "eya"),
        split.by = "strain", cols = c("#C83658", "#3A8A00", "#0088C5"))

# Main Body Follicle Cells
dotplot3 <- DotPlot(ovary, features = c("Vha16-1", "yellow-g", "yellow-g2", "Yp1", "Sox14", "senju", "Acer", "CG31926", "Femcoat"),
        split.by = "strain", cols = c("#C83658", "#3A8A00", "#0088C5"))

# Ovarian Sheath Muscle + Oviduct
dotplot4 <- DotPlot(ovary, features = c("Act57B", "Mhc", "Mlc2", "Mlp60A", "Mlp84B", "abd-A", "vir-1"),
        split.by = "strain", cols = c("#C83658", "#3A8A00", "#0088C5"))

# Combine dotplots
combined_dotplot <- dotplot1 + dotplot2 + dotplot3 + dotplot4
print(combined_dotplot)

# Cluster annotation
new.cluster.ids <- c(
  "Germline - GSC / Germarium 1 and 2a", "Germline - Germarium 2a-2b", "Germline - Germarium 2a-2b",
  "Germline - Germarium 2a-2b", "Germline - Germarium 2a-2b", "Germline - Germarium 2b-3",
  "Stalk & Polar Cells", "Follicle Stem Cells (FSC) / preFCs", "Early Follicle Cells",
  "Mitotic Follicle Cells Stage 1-5", "Post-Mitotic Follicle Cells Stage 6", "Post-Mitotic Follicle Cells Stage 6",
  "Post-Mitotic Follicle Cells Stage 6", "Post-Mitotic Follicle Cells Stage 6", "Post-Mitotic Follicle Cells Stage 6",
  "Post-Mitotic Follicle Cells Stage 6", "Vitellogenic MBFCs Stage 7", "Vitellogenic MBFCs Stage 7",
  "Vitellogenic MBFCs Stage 8", "Vitellogenic MBFCs Stage 9-10A", "Vitellogenic MBFCs Stage 9-10A",
  "Vitellogenic MBFCs Stage 9-10A", "Choriogenic MBFCs Stage 12", "Choriogenic MBFCs Stage 12",
  "Choriogenic MBFCs Stage 14", "Terminal Corpus Luteum Cells", "Stretch Cells",
  "Ovarian Sheath Muscle", "Oviduct", "Unannotated"
)
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)

# Set preferred cluster order
my_levels <- new.cluster.ids
Idents(ovary) <- factor(Idents(ovary), levels = my_levels)

# UMAP plot
umap_plot <- DimPlot(
  ovary, reduction = "umap.harmony", label = TRUE,
  label.size = 2, repel = TRUE, label.box = TRUE
)
umap_plot + xlab("UMAP_1") + ylab("UMAP_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

# t-SNE plot
tsne_plot <- DimPlot(
  ovary, reduction = "tsne.harmony", label = TRUE,
  label.size = 2, repel = TRUE, label.box = TRUE
)
tsne_plot + xlab("tSNE_1") + ylab("tSNE_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5)))

# Save full object with annotations
saveRDS(ovary, file = "./ovary.rds")

# Remove unannotated clusters for downstream analysis
ovary <- subset(ovary, idents = "Unannotated", invert = TRUE)
saveRDS(ovary, file = "./FINAL-ovary.rds")
