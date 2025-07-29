# Differential Expression Analysis Ovary
# Hariyani et al., 2025

# Load libraries
library(Seurat)
library(MAST)
library(dplyr)

# Load processed testis Seurat object
testis <- readRDS("~/FINAL-ovary-nounknown.rds")

# Load gene lists for both D. melanogaster strains
ISO1_genes <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")
A4_genes   <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")
w501_genes   <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")

# Identify shared gene names
common_gene_names <- Reduce(intersect, list(ISO1_genes$gene_name, A4_genes$gene_name, w501_genes$gene_name))

# Subset Seurat object to retain only shared genes
testis <- subset(testis, features = common_gene_names)

# Create composite label for cell type and strain
testis$celltype.strain <- paste(testis$celltype, testis$strain, sep = "_")

# Initialize storage for differential expression results
DE_results_list <- list()
unique_celltypes <- unique(testis$celltype)

# List of all pairwise comparisons
comparisons <- list(c("ISO1", "A4"), c("ISO1", "w501"), c("A4", "w501"))

# Perform DE for each pairwise comparison per cell type
for (comp in comparisons) {
  strain1 <- comp[1]
  strain2 <- comp[2]
  comp_label <- paste0(strain1, "_vs_", strain2)
  
  for (celltype in unique_celltypes) {
    ident_1 <- paste(celltype, strain1, sep = "_")
    ident_2 <- paste(celltype, strain2, sep = "_")
    
    Idents(testis) <- "celltype.strain"
    DE_results <- FindMarkers(testis, ident.1 = ident_1, ident.2 = ident_2, min.pct = 0, logfc.threshold = 0, test.use = "MAST")
    DE_results_list[[paste(celltype, comp_label, sep = "_")]] <- DE_results
  }
}

# Classify genes by direction of regulation
upregulated_genes_list   <- list()
downregulated_genes_list <- list()
ns_genes_list            <- list()

for (celltype in names(DE_results_list)) {
  results <- DE_results_list[[celltype]]
  
  up_genes   <- rownames(results[results$avg_log2FC > 1  & results$p_val_adj < 0.01, ])
  down_genes <- rownames(results[results$avg_log2FC < -1 & results$p_val_adj < 0.01, ])
  ns_genes   <- rownames(results[!(results$p_val_adj < 0.01 & abs(results$avg_log2FC) > 1), ])
  
  upregulated_genes_list[[celltype]]   <- up_genes
  downregulated_genes_list[[celltype]] <- down_genes
  ns_genes_list[[celltype]]            <- ns_genes
}

# Combine into a single dataframe
combined_genes <- data.frame()

combine_genes <- function(gene_list, celltype, regulation) {
  genes <- unlist(gene_list[[celltype]])
  df <- DE_results_list[[celltype]][genes, c("p_val_adj", "avg_log2FC", "pct.1", "pct.2")]
  df$gene       <- rownames(df)
  df$celltype   <- celltype
  df$regulation <- regulation
  return(df[, c("gene", "celltype", "regulation", "p_val_adj", "avg_log2FC", "pct.1", "pct.2")])
}

for (ct in names(upregulated_genes_list)) {
  combined_genes <- rbind(combined_genes, combine_genes(upregulated_genes_list, ct, "Up"))
  combined_genes <- rbind(combined_genes, combine_genes(downregulated_genes_list, ct, "Down"))
  combined_genes <- rbind(combined_genes, combine_genes(ns_genes_list, ct, "NS"))
}

# Export DE tables per cell type
filter_and_write <- function(celltype, filename) {
  subset_df <- combined_genes %>% 
    filter(celltype == !!celltype) %>%
    mutate(`min.pct>=0.01` = pct.1 >= 0.01 | pct.2 >= 0.01)
  write.table(subset_df, filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Define output filenames for each cell type
celltypes <- c("Germline - GSC / Germarium 1 and 2a","Germline - Germarium 2a-2b","Germline - Germarium 2b-3","Stalk & Polar Cells",
                 "Follicle Stem Cells (FSC) / preFCs","Early Follicle Cells","Mitotic Follicle Cells Stage 1-5",
                 "Post-Mitotic Follicle Cells Stage 6",
                 "Vitellogenic MBFCs Stage 7","Vitellogenic MBFCs Stage 8","Vitellogenic MBFCs Stage 9-10A",
                 "Choriogenic MBFCs Stage 12",
                 "Choriogenic MBFCs Stage 14","Stretch Cells","Ovarian Sheath Muscle","Oviduct"

)

# Get all keys from DE_results_list (formatted as "celltype_comparison")
all_comparisons <- names(DE_results_list)

# Create output filenames for each
filenames <- paste0("DGE_", gsub(" ", "-", gsub("/", "_", all_comparisons)), ".pval.txt")

# Write tables using updated names
for (i in seq_along(all_comparisons)) {
  filter_and_write(all_comparisons[i], filenames[i])
}
