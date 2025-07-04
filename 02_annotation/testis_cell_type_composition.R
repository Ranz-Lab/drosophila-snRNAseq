## Hariyani et al., 2025
## 02_annotation - Testis cell type composition

# Load libraries
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(dplyr)
library(SCpubr)

# Load object
testis <- readRDS("FINAL-testis-nounknown.rds")

# ---------------------------------------
# Assign active identity to metadata
# ---------------------------------------
testis@meta.data$celltype <- testis@active.ident

# ---------------------------------------
# Quick composition visualization
# ---------------------------------------
dittoBarPlot(
  object = testis,
  var = "ident",
  group.by = "strain"
) + labs(x = "Strain", y = "Proportion of Cells")

# ---------------------------------------
# Calculate proportions
# ---------------------------------------
proportion_data <- testis@meta.data %>%
  group_by(celltype, strain) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(celltype) %>%
  mutate(proportion = n / sum(n))

# ---------------------------------------
# Color palette
# ---------------------------------------
cbPalette <- c("lightblue", "lightblue1", "lightblue3", "steelblue1", "steelblue2",
               "#ffb6db", "lightpink", "#ff6db6", "#CC79A7",
               "#F0E442", "#ffff6d", "khaki1", "khaki3", "khaki4",
               "#009292", "#999999")

# ---------------------------------------
# Raw count barplot
# ---------------------------------------
cell_type_composition <- ggplot(testis@meta.data, aes(x = strain, fill = celltype)) +
  geom_bar() +
  scale_fill_manual(values = cbPalette) +
  geom_text(
    aes(label = ..count..),
    stat = "count",
    size = 3,
    position = position_stack(vjust = 0.5)
  )
cell_type_composition

# ---------------------------------------
# Proportion barplot
# ---------------------------------------
cell_type_proportion <- ggplot(proportion_data, aes(x = celltype, y = proportion, fill = strain)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cbPalette) +
  geom_text(
    aes(label = scales::percent(proportion)),
    position = position_stack(vjust = 0.5),
    size = 3
  )
cell_type_proportion

# ---------------------------------------
# SCpubr visualizations
# ---------------------------------------
proportion1_scpubr <- SCpubr::do_BarPlot(
  testis,
  group.by = "celltype",
  split.by = "strain",
  position = "fill",
  flip = FALSE
)

proportion2_scpubr <- SCpubr::do_BarPlot(
  testis,
  group.by = "strain",
  split.by = "celltype",
  position = "fill",
  flip = FALSE
)

# ---------------------------------------
# Save output
# ---------------------------------------
write.csv(proportion_data, "cell_type_composition_proportion_testis.csv", row.names = FALSE)
