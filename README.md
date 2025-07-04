# drosophila-snRNAseq

This repository contains all code used in the analysis of single-nucleus RNA-seq data from testis and ovary tissue of three Drosophila strains. The preprint for this study is available on *bioRxiv*.

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
