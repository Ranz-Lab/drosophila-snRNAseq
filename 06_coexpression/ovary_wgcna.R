# install additional packages:
#BiocManager::install(c("WGCNA", "igraph", "devtools", "GeneOverlap", "ggrepel", "UCell"))
#Sys.unsetenv("GITHUB_PAT")
#devtools::install_github("NightingaleHealth/ggforestplot")
#devtools::install_github('smorabit/hdWGCNA', ref='dev')

# load packages
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# Load ovary dataset
ovary <- readRDS("~/FINAL-ovary.rds")

# rename GSC & FSC
current_idents <- ovary@active.ident
levels(current_idents)[levels(current_idents) == "Follicle Stem Cells (FSC) / preFCs"] <- "Follicle Stem Cells"
levels(current_idents)[levels(current_idents) == "Germline - GSC / Germarium 1 and 2a"] <- "Germline - Germarium 1 and 2a"
ovary@active.ident <- current_idents
ovary@meta.data[["celltype"]] <- as.character(ovary@meta.data[["celltype"]])
ovary@meta.data[["celltype"]][ovary@meta.data[["celltype"]] == "Follicle Stem Cells (FSC) / preFCs"] <- "Follicle Stem Cells"
ovary@meta.data[["celltype"]][ovary@meta.data[["celltype"]] == "Germline - GSC / Germarium 1 and 2a"] <- "Germline - Germarium 1 and 2a"

# Import gene names for all 3 strains
ISO1_genes <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")
A4_genes <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")
w501_genes <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")

# Identify genes with data for all three species
common_gene_names <- Reduce(intersect, list(ISO1_genes$gene_name, A4_genes$gene_name, w501_genes$gene_name))

# Subset the Seurat object to retain only these genes
ovary <- subset(ovary, features=common_gene_names)

# Create a species column based on the strain information
#ovary$species <- ifelse(ovary$strain %in% c("A4", "ISO1"), "Dmel", 
                         #ifelse(ovary$strain == "w501", "Dsim", "Unknown"))
seurat_obj <- ovary

# plot UMAP to check
#p <- DimPlot(seurat_obj, group.by='celltype', label=TRUE) +
# umap_theme() + ggtitle('ovary of Drosophila melanogaster and simulans') + NoLegend()

#p

# select genes
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "realrun-1_ovary" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype", "strain"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 15, # maximum number of shared cells between two metacells
  ident.group = 'celltype', # set the Idents of the metacell seurat object
  min_cells = 50
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

# List of cell types
set_columns <- c("Germline - Germarium 1 and 2a","Germline - Germarium 2a-2b","Germline - Germarium 2b-3","Stalk & Polar Cells","Follicle Stem Cells","Early Follicle Cells",
                 "Mitotic Follicle Cells Stage 1-5","Post Mitotic Follicle Cells Stage 6","Vitellogenic MBFCs Stage 7",
                 "Vitellogenic MBFCs Stage 8","Vitellogenic MBFCs Stage 9 10A","Choriogenic MBFCs Stage 12","Choriogenic MBFCs Stage 14",
                 "Terminal Corpus Luteum Cells","Stretch Cells","Ovarian Sheath Muscle","Oviduct")

# Loop through each cell type
for (cell_type in set_columns) {
  
  temp_obj <- seurat_obj
  
  # Setup WGCNA for the specific cell type
  temp_obj <- SetDatExpr(
    temp_obj,
    group_name = cell_type,  # use the current cell type
    group.by = 'celltype',   # the metadata column containing the cell type info
    assay = 'RNA',           # using RNA assay
    slot = 'data'            # using normalized data
  )
  
  # Test different soft powers
  temp_obj <- TestSoftPowers(
    temp_obj,
    networkType = 'signed'  # you can also use "unsigned" or "signed hybrid"
  )
  
  # Plot soft power results
  plot_list <- PlotSoftPowers(temp_obj)
  
  # Plot_list is a list, combine them into one plot
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  
  # Save the combined plot
  ggsave(paste0("SoftPowers_ovary_", cell_type, ".png"), combined_plot, width = 10, height = 6)
  
  # Construct coexpression network for the current cell type
  temp_obj <- ConstructNetwork(
    temp_obj,
    tom_name = paste0(cell_type, '_network'), # dynamically name the TOM file
    deepSplit = 4,
    minModuleSize = 25,
    overwrite_tom = TRUE
  )
  
  # Plot dendrogram
  PlotDendrogram(temp_obj, main = paste0(cell_type, ' ovary hdWGCNA Dendrogram'))
  #ggsave(paste0("Dendrogram_ovary_", cell_type, ".png"), dendro_plot, width=10, height=6)
  
  # Compute module eigengenes (MEs)
  temp_obj <- ModuleEigengenes(
    temp_obj,
    group.by.vars="strain"  # Using strain as the grouping variable
  )
  
  # Harmonized module eigengenes (hMEs)
  hMEs <- GetMEs(temp_obj)
  
  # Module eigengenes (MEs)
  MEs <- GetMEs(temp_obj, harmonized = FALSE)
  
  write.csv(hMEs, file = paste0("hMEs_ovary_", cell_type, ".csv"))
  write.csv(MEs, file = paste0("MEs_ovary_", cell_type, ".csv"))
  
  # Compute eigengene-based connectivity (kME)
  temp_obj <- ModuleConnectivity(
    temp_obj,
    group.by = 'celltype', 
    group_name = cell_type
  )
  
  # Reset module names with the cell type
  temp_obj <- ResetModuleNames(
    temp_obj,
    new_name = paste0(cell_type, "-M")
  )
  
  # Plot genes ranked by kME for each module
  PlotKMEs(temp_obj, ncol = 3)
  #ggsave(paste0("kME_plot_ovary_", cell_type, ".png"), p, width=12, height=8)
  
  # get the module assignment table:
  modules <- GetModules(temp_obj) %>% subset(module != 'grey')
  write.csv(modules, file = paste0("modules_ovary_", cell_type, ".csv"))  
  
  # get hub genes
  hub_df <- GetHubGenes(temp_obj, n_hubs = 10)
  write.csv(hub_df, file = paste0("hubgenes_ovary_", cell_type, ".csv"))
  
  #Differential module eigenvalue analysis
  groupa <- temp_obj@meta.data %>% subset(celltype == cell_type & strain == "ISO1") %>% rownames
  groupb <- temp_obj@meta.data %>% subset(celltype == cell_type & strain == "A4") %>% rownames
  groupc <- temp_obj@meta.data %>% subset(celltype == cell_type & strain == "w501") %>% rownames
  
  #ISO1 vs A4
  DMEs_intra <- FindDMEs(
    temp_obj,
    barcodes1 = groupa,
    barcodes2 = groupb,
    test.use='wilcox',
    wgcna_name='realrun-1_ovary'
  )
  
  #IS01 vs w501
  DMEs_inter1 <- FindDMEs(
    temp_obj,
    barcodes1 = groupa,
    barcodes2 = groupc,
    test.use='wilcox',
    wgcna_name='realrun-1_ovary'
  )
  
  #A4 vs w501
  DMEs_inter2 <- FindDMEs(
    temp_obj,
    barcodes1 = groupb,
    barcodes2 = groupc,
    test.use='wilcox',
    wgcna_name='realrun-1_ovary'
  )
  
  #Save DMEs
  write.csv(DMEs_intra, file = paste0("DMEs_ISO1vsA4_ovary_", cell_type, ".csv"))
  write.csv(DMEs_inter1, file = paste0("DMEs_ISO1vsw501_ovary_", cell_type, ".csv"))
  write.csv(DMEs_inter2, file = paste0("DMEs_A4vsw501_ovary_", cell_type, ".csv"))
  
  #Plot
  PlotDMEsLollipop(
    temp_obj, 
    DMEs_intra, 
    wgcna_name="realrun-1_ovary",
    pvalue = "p_val_adj"
  )
  
  PlotDMEsLollipop(
    temp_obj, 
    DMEs_inter1, 
    wgcna_name="realrun-1_ovary", 
    pvalue = "p_val_adj"
  )
  
  PlotDMEsLollipop(
    temp_obj, 
    DMEs_inter2, 
    wgcna_name="realrun-1_ovary", 
    pvalue = "p_val_adj"
  )
  
  #Save plots
  #ggsave(paste0("DMEs_ISO1vsA4_ovary_", cell_type, ".png"), dmeplot1, width=12, height=8)  
  #ggsave(paste0("DMEs_ISO1vsw501_ovary_", cell_type, ".png"), dmeplot2, width=12, height=8)
  #ggsave(paste0("DMEs_A4vsw501_ovary_", cell_type, ".png"), dmeplot3, width=12, height=8)
  
  # Save the object for each cell type
  saveRDS(temp_obj, file = paste0('hdWGCNA_ovary_', cell_type, '_realrun-1.rds'))
  
}

#Dmel vs Dsim

# List of cell types
set_columns <- c("Ovarian Sheath Muscle","Oviduct")

# Loop through each cell type
for (cell_type in set_columns) {
  
  rdsFile <- paste0('hdWGCNA_ovary_', cell_type, '_realrun-1.rds')
  
  temp_obj <- readRDS(rdsFile)
  
  dmel <- temp_obj@meta.data %>% subset(celltype == cell_type & species == "Dmel") %>% rownames
  dsim <- temp_obj@meta.data %>% subset(celltype == cell_type & species == "Dsim") %>% rownames
  
  DMEs_inter <- FindDMEs(
    temp_obj,
    barcodes1 = dmel,
    barcodes2 = dsim,
    test.use='wilcox',
    min.pct=0.01,
    wgcna_name='realrun-1_ovary'
  )
  
  DMEfigure <- PlotDMEsLollipop(
    temp_obj, 
    DMEs_inter, 
    wgcna_name="realrun-1_ovary", 
    pvalue = "p_val_adj"
  )
  
  ggsave(paste0("DMEs_DmelvsDsim_ovary_", cell_type, ".jpeg"),
         plot = last_plot(),
         device = NULL,
         path = NULL,
         scale = 1,
         width = NA,
         height = NA,
         units = c("in", "cm", "mm", "px"),
         dpi = 300,
         limitsize = TRUE,
         bg = NULL)
  
  write.csv(DMEs_inter, file = paste0("DMEs_DmelvsDsim_ovary_", cell_type, ".csv"))
}


#Enrichment Analysis
# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# enrichment packages
library(clusterProfiler) 
library(org.Dm.eg.db)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

set_columns <- c("Vitellogenic MBFCs Stage 9-10A")

for (cell_type in set_columns) {
  rdsFile <- paste0('hdWGCNA_ovary_', cell_type, '_realrun-1.rds')
  
  # load the seurat object
  seurat_obj <- readRDS(rdsFile)

  # Run GO enrichment analysis using the modified function
  seurat_obj <- RunGOClusterProfilertest(
    seurat_obj, 
    organism = "org.Dm.eg.db", # Specify organism database; adjust if needed
    max_genes = 500            # Number of genes per module to test; use max_genes = Inf for all genes
  )

  # retrieve the output table
  enrich_df <- GetEnrichrTable(seurat_obj)
  
  # Save the enrichment table as a .csv file
  output_table_file <- paste0("hdWGCNA_ovary_GOenrichment_", cell_type, ".csv")
  write.csv(enrich_df, output_table_file, row.names = FALSE)

  # make GO term plots:
  ClusterProfilerDotPlot(
    seurat_obj,
    outdir = paste0('hdWGCNA_ovary_clusterprofiler_plots', cell_type), # name of output directory
    n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
    plot_size = c(8,12), # width, height of the output .pdfs
    logscale=TRUE # do you want to show the enrichment as a log scale?
  )

  # comparison dotplot across all modules within cell type
  dotplot <- ClusterProfilerComparisonPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    ontology = "BP", # this has to be one of the lists we used above!!!
    n_terms=2 # number of terms for each module
  )

  #save the figure
  ggsave(paste0("DmelvsDsim_ovary_GOenrichment_dotplot_", cell_type, ".jpeg"),
         plot = last_plot(),
         device = NULL,
         path = NULL,
         scale = 1,
         width = NA,
         height = NA,
         units = c("in", "cm", "mm", "px"),
         dpi = 300,
         limitsize = TRUE,
         bg = NULL)
}

##########################
# Top 10 Hub Genes GO Enrichment

# Load required libraries
library(clusterProfiler)
library(org.Dm.eg.db)  
library(dplyr)

# Simple GO enrichment function
SimpleGOEnrichment <- function(gene_symbols, organism_db = org.Dm.eg.db, p_cutoff = 0.05, q_cutoff = 0.2) {
  # Convert gene symbols to Entrez IDs
  mapped_genes <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism_db)
  
  if (nrow(mapped_genes) == 0) {
    stop("No genes could be mapped to Entrez IDs.")
  }
  
  # Perform GO enrichment
  ego <- enrichGO(
    gene          = mapped_genes$ENTREZID,
    OrgDb         = organism_db,
    keyType       = "ENTREZID",
    ont           = "ALL",  # options: "BP", "CC", "MF", or "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff  = p_cutoff,
    qvalueCutoff  = q_cutoff
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message("No significant GO enrichment results found.")
    return(NULL)
  }
  
  return(as.data.frame(ego))
}

genes <- read.table("~/Desktop/top1%-hub-genes-ovary.txt", quote="\"", comment.char="")
gene_list <- as.character(genes$V1)
go_results <- SimpleGOEnrichment(gene_list)

write.csv(go_results, "ovary_hub_top1%_go_results.csv")

genes <- read.table("~/Desktop/top10-hub-genes-ovary.txt", quote="\"", comment.char="")
gene_list <- as.character(genes$V1)
go_results <- SimpleGOEnrichment(gene_list)

write.csv(go_results, "ovary_hub_top10_go_results.csv")

##########################
#EnrichR - OLD CODE
# enrichr databases to test
dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

# enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2018", # this has to be one of the lists we used above!!!
  n_terms=2 # number of terms for each module
)
