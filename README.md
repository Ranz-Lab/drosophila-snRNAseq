# drosophila-snRNAseq

This repository contains all code used in the analysis of single-nucleus RNA-seq data from testis and ovary tissue of three *Drosophila* strains. The preprint for this study is available on [*bioRxiv*](https://www.biorxiv.org/content/10.1101/2025.07.04.663238v1).

## 📁 Repository Structure

| Level | Folder | What it contains |
|-------|--------|------------------|
| **0** | **`00_preprocessing/`** | Parse analysis pipeline: demultiplexing, reference preparation, alignment, and raw expression output |
|       |   ├── `00_demultiplexing/` 
|       |   ├── `01_prepare_references/` 
|       |   └── `02_parse_pipeline/` 
| **1** | **`01_cross_species_integration/`** | Ortholog mapping (w501 → Dmel), QC, doublet removal, normalization, Harmony integration |
| **2** | **`02_annotation/`** | Marker-gene expression and annotation, trajectory inference, cell-type composition analysis |
| **3** | **`03_pseudobulk/`** | Pseudobulk aggregation by cell type and sample |
| **4** | **`04_expression_correlation/`** | Expression correlation across strains / cell types |
| **5** | **`05_DE/`** | Differential expression analysis between species |
| **6** | **`06_coexpression/`** | Co-expression network analysis and divergence between species |

---

## 📦 Data Access

Raw FASTQ files are available under NCBI BioProject [PRJNA1273697](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1273697).
