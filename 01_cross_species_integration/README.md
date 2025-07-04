# 01_cross_species_integration

This directory contains scripts and metadata for processing and integrating single-nucleus RNA-seq data from the testis and ovary of three *Drosophila* strains—*D. melanogaster* (ISO1, A4) and *D. simulans* (w501). This stage harmonizes gene identifiers across species, constructs filtered Seurat objects for each tissue and strain, removes doublets, performs normalization, and integrates the datasets using Harmony.

---

## Cross-Species Gene Mapping

To enable cross-species comparisons, *D. simulans* gene identifiers were converted to their *D. melanogaster* orthologs using:

- `mel-sim_orthologs.txt`: a two-column tab-delimited file with `Dmel_gene_id <tab> Dsim_gene_id`.

The mapping was performed using the script:

- `replace_gene_ortholog.py`: replaces gene IDs in `all_genes_w501.csv` with their *D. melanogaster* orthologs.

**Important note:** Gene IDs without known orthologs were retained to preserve species-specific expression profiles.

---

## Parse Digital Gene Expression Files Used

For each sample (testis and ovary), we downloaded the following unfiltered DGE files after running the Parse Biosciences pipeline:

- `count_matrix.mtx`
- `cell_metadata.csv`
- `all_genes.csv`

For *D. simulans* (w501), `all_genes.csv` was replaced by the ortholog-adjusted file `all_genes_orthologs.csv`.

Each sample was used to construct a separate Seurat object, which was later merged by tissue (testis or ovary).

---

## Seurat Object Construction and Filtering

Seurat v5 was used with R 4.3+ to generate objects for each tissue and strain via custom functions, `relaxed_seurat()` and `stringent_seurat()` for the testis and ovary samples respectively. Separate Seurat objects were built for each of the 12 samples.

### QC filtering criteria:

- `nFeature_RNA ≥ 200`
- `nCount_RNA ≥ 300 & ≤ 100000` (testis) and `nCount_RNA ≥ 500 & ≤ 100000` (ovary)
- `nFeature_RNA < 7500`
- `percent.mt ≤ 5%`

Violin plots and scatter plots (feature counts, mitochondrial content) were visualized for all datasets.

---

## Doublet Detection and Removal

DoubletFinder was used to detect and remove doublets from each tissue-strain combination. Key steps:

- Identification of optimal `pK` via parameter sweep.
- Expected doublet rate set to 3%.
- Homotypic proportion estimated from unintegrated clusters.
- Singlets were extracted based on `DF.classifications_*` and retained for downstream analysis.

---

## Normalization

After filtering and doublet removal, NormalizeData was applied individually to each Seurat object (testis and ovary for all three strains). PCA, variable feature selection, and scaling were also performed.

---

## Merging and Integration with Harmony

### Merging:
For each tissue, the Seurat objects corresponding to the three strains were merged into a single object using merge(), while preserving RNA expression layers via Seurat v5’s multilayer architecture.

### Integration:
Batch effects across strains were removed using `HarmonyIntegration` while preserving cross-species heterogeneity.
