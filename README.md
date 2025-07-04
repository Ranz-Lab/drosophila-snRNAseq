# drosophila-snRNAseq

This repository contains all code used in the analysis of single-nucleus RNA-seq data from testis and ovary tissue of three Drosophila strains. The preprint for this study is available on *bioRxiv*.

## 📁 Repository Structure

| Level | Folder | What it contains |
|-------|--------|------------------|
| **0** | **`00_preprocessing/`** | Data-handling steps that run **before** downstream R analysis  |
|       |   ├── `00_demultiplexing/` | well layout (`sample_sheet.csv`) & SLURM script to split raw FASTQs by **strain** |
|       |   ├── `01_prepare_references/` | Liftoff annotation transfer, protein-coding filters, and `split-pipe --mode mkref` scripts |
|       |   └── `02_parse_pipeline/` | SLURM jobs for Parse align (`01_parse_align_*`) and sublibrary combiner (`02_parse_combine_*`) |
| **1** | **`01_cross_species_integration/`** | Ortholog mapping (w501 → Dmel), QC, doublet removal, normalization, Harmony integration |
| **2** | **`02_annotation/`** | Marker-gene dot plots, manual cluster renaming, Monocle3 trajectories, cell-type composition |
| **3** | **`03_pseudobulk/`** | Pseudobulk aggregation and QC (empty placeholder for now) |
| **4** | **`04_expression_correlation/`** | Expression correlation across strains / cell types (placeholder) |
| **5** | **`05_DE/`** | Differential-expression analyses (placeholder) |
| **6** | **`06_coexpression/`** | Co-expression network inference (placeholder) |
| **7** | **`07_evolutionary_analysis/`** | Gene age, tissue-specificity (τ), Ka/Ks, etc. (placeholder) |

> *Empty analysis folders already contain a short “Work-in-progress” README so the tree is clear.*

---

## 📁 Repository Structure

- `00_preprocessing/` — Reference preparation and Parse analysis pipeline  
- `01_cross_species_integration/` — Ortholog mapping, normalization and integration
- `02_annotation/` — Cell type annotation, trajectory inference, and cell type composition analysis  
- `03_pseudobulk/` — Aggregation of gene counts by cell type and sample  
- `04_expression_correlation/` — Expression correlation across strains and cell types  
- `05_DE/` — Differential expression analysis between species 
- `06_coexpression/` — Coexpression network inference and comparison  
- `07_evolutionary_analysis/` — Gene age, tau, and evolutionary rate analyses


## 📦 Data Access

Raw FASTQ files are available under NCBI BioProject [PRJNA1273697](https://...).
