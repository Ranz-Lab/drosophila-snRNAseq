## Hariyani et al., 2025
## 02_annotation - Ovary cell type composition

# Load libraries
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(dplyr)

# Load object
ovary <- readRDS("FINAL-ovary.rds")

# ---------------------------------------
# Assign active identity to metadata
# ---------------------------------------
ovary@meta.data$celltype <- ovary@active.ident

# ---------------------------------------
# Quick composition visualization
# ---------------------------------------
dittoBarPlot(
  object = ovary,
  var = "ident",
  group.by = "strain"
) + labs(x = "Strain", y = "Proportion of Cells")

# ---------------------------------------
# Calculate proportions
# ---------------------------------------
proportion_data <- ovary@meta.data %>%
  group_by(strain, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(strain) %>%
  mutate(proportion = n / sum(n))

# ---------------------------------------
# Color palette
# ---------------------------------------
cbPalette <- c(
  "lightblue", "lightblue1", "lightblue3", "lightblue4", "steelblue1", "steelblue2", "steelblue3",
  "#ffb6db", "#ff6db6", "#CC79A7", "#D55E00", "#E69F00", "#F0E442", "#ffff6d",
  "khaki1", "khaki3", "#006ddb", "#999999", "#56B4E9", "#6db6ff", "#b6dbff", "#009292", "#009E73",
  "springgreen4", "springgreen3", "springgreen2", "springgreen", "#24ff24", "lawngreen",
  "indianred", "indianred1", "indianred2", "indianred3", "indianred4",
  "#920000", "#924900", "#db6d00", "#999999"
)

# ---------------------------------------
# Raw count barplot
# ---------------------------------------
cell_type_composition <- ggplot(ovary@meta.data, aes(x = strain, fill = celltype)) +
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
  scale_fill_manual(values = c("steelblue1", "indianred1", "lightgreen")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
cell_type_proportion

# ---------------------------------------
# Save output
# ---------------------------------------
write.csv(proportion_data, "cell_type_composition_proportion_ovary.csv", row.names = FALSE)
