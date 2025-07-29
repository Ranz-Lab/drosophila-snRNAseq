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

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the testis snRNA-seq dataset
testis <- readRDS("~/FINAL-testis-nounknown.rds")

# Import gene names for all 3 strains
ISO1_genes <- read.csv("~/Desktop/Parse_analysis_v3/ISO1/DGE_unfiltered/all_genes.csv")
A4_genes <- read.csv("~/Desktop/Parse_analysis_v3/A4/DGE_unfiltered/all_genes.csv")
w501_genes <- read.csv("~/Desktop/Parse_analysis_v3/w501/DGE_unfiltered/all_genes.csv")

# Identify genes with data for all three species
common_gene_names <- Reduce(intersect, list(ISO1_genes$gene_name, A4_genes$gene_name, w501_genes$gene_name))

# Subset the Seurat object to retain only these genes
testis <- subset(testis, features=common_gene_names)

# Create a species column based on the strain information
testis$species <- ifelse(testis$strain %in% c("A4", "ISO1"), "Dmel", 
                         ifelse(testis$strain == "w501", "Dsim", "Unknown"))
seurat_obj <- testis

# plot UMAP to check
p <- DimPlot(seurat_obj, group.by='celltype', label=TRUE) +
  umap_theme() + ggtitle('Testis of Drosophila melanogaster and simulans') + NoLegend()

p

# select genes
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "run-4" # the name of the hdWGCNA experiment
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

# perform 
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Late Spermatocytes", # the name of the group of interest in the group.by column
  group.by='celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct coexpression network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'Late Spermatocytes3', # name of the topoligical overlap matrix written to disk
  deepSplit = 4,
  minModuleSize = 25
)

# plot coexpression
PlotDendrogram(seurat_obj, main='Late Spermatocytes hdWGCNA Dendrogram')

# topoligcal overlap matrix - square matrix of genes by genes, where each value is the topoligcal overlap between the genes
TOM <- GetTOM(seurat_obj)

# compute all MEs in the full single-cell dataset  to summarize the gene expression profile of an entire co-expression module
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="strain"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype', group_name = 'Late Spermatocytes'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Late Spermatocytes-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=3)

p

# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA_testis_latespermatocytes-run3.rds')

##### VISUALIZATION #######

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction="umap.harmony",
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction="umap.harmony",
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)

# plot module correlagram
ModuleCorrelogram(seurat_obj)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'celltype')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

# individual module network plots
ModuleNetworkPlot(
  seurat_obj,
  outdir = 'ModuleNetworks'
)

ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks2', # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[1:2] # only select 2 modules
)

# umap 
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
)

g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)

##### DIFFERENTIAL EIGENVALUE #######

#Intraspecific

groupa <- seurat_obj@meta.data %>% subset(celltype == 'Late Spermatocytes' & strain == "ISO1") %>% rownames
groupb <- seurat_obj@meta.data %>% subset(celltype == 'Late Spermatocytes' & strain == "A4") %>% rownames

head(groupa)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = groupa,
  barcodes2 = groupb,
  test.use='wilcox',
  wgcna_name='run-3'
)

head(DMEs)

PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='run-3', 
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'run-3'
)


#Interspecific

seurat_obj$species <- ifelse(seurat_obj$strain %in% c("A4", "ISO1"), "Dmel", 
                             ifelse(seurat_obj$strain == "w501", "Dsim", "Unknown"))

group1 <- seurat_obj@meta.data %>% subset(celltype == 'Late Spermatocytes' & species == "Dmel") %>% rownames
group2 <- seurat_obj@meta.data %>% subset(celltype == 'Late Spermatocytes' & species == "Dsim") %>% rownames

head(group1)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='run-3'
)

head(DMEs)

PlotDMEsLollipop(
  seurat_obj, 
  DMEs, 
  wgcna_name='run-3', 
  pvalue = "p_val_adj"
)

PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'run-3'
)
