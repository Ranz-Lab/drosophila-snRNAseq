# Differential Expression UpSet Plot - Testis
# Hariyani et al., 2025

# Load libraries
library(Seurat)
library(MAST)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

# Load DE results
combined_results <- read.delim("combined_results.txt")

# Extract unique cell types and comparisons
cell_types <- unique(combined_results$celltype)
comparisons <- unique(combined_results$comparison)

# Define ordered cell types for plotting
set_columns <- c(
  "Spermatids", "Maturing Primary Spermatocytes", "Late Spermatocytes",
  "Mid Spermatocytes", "Early Spermatocytes", "Late Spermatogonia",
  "Early Spermatogonia", "Cyst Cells", "Hub Cells", "Epithelial Cells"
)

# Iterate over all comparisons to generate UpSet plots
for (comparison in comparisons) {
  message("Processing: ", comparison)
  
  # Build gene list for current comparison
  DE_results_list <- list()
  for (celltype in cell_types) {
    subset_df <- combined_results %>%
      filter(celltype == !!celltype, comparison == !!comparison)
    
    up_genes   <- subset_df$gene[subset_df$avg_log2FC > 1  & subset_df$adj_p_val < 0.01]
    down_genes <- subset_df$gene[subset_df$avg_log2FC < -1 & subset_df$adj_p_val < 0.01]
    
    all_genes <- unique(c(up_genes, down_genes))
    DE_results_list[[celltype]] <- all_genes
  }

  # Build long-form data frame of DEGs
  combined_genes <- data.frame()
  for (celltype in names(DE_results_list)) {
    genes <- DE_results_list[[celltype]]
    if (length(genes) > 0) {
      combined_genes <- rbind(combined_genes, data.frame(
        gene = genes,
        celltype = celltype,
        regulation = 'DE'
      ))
    }
  }

  # Convert to wide format
  wide_genes <- combined_genes %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = celltype, values_from = presence, values_fill = list(presence = 0))

  # Generate and save UpSet plot
  plot <- upset(
    wide_genes,
    set_columns,
    n_intersections = 30,
    min_size = 1,
    width_ratio = 0.1,
    sort_sets = FALSE,
    base_annotations = list(
      'Intersection size' = intersection_size(counts = FALSE, mapping = aes(fill = 'bars_color')) +
        scale_fill_manual(values = c('bars_color' = 'lightblue'), guide = 'none') +
        stat_summary(fun = sum, geom = 'text', position = position_stack(vjust = 0.5), aes(label = after_stat(y))) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = 'grey')
        ) +
        ylim(c(0, 650))
    ),
    set_sizes = upset_set_size() +
      ylab('Number of DEGs') +
      theme(axis.text.x = element_text(angle = 90)) +
      theme(axis.ticks.x = element_line()) +
      ylim(c(1570, 0))
  )

  ggsave(paste0("UpSet_", gsub(" ", "_", comparison), ".pdf"), plot, width = 10, height = 8)
}
