# 03_pseudobulk

This directory contains scripts for pseudobulk analysis of single-nucleus RNA-seq data from *Drosophila* testis and ovary. This stage aggregates gene expression counts by cell type and sample to enable bulk-like comparisons across strains and plot a PCA.

---

###5. `testis_pseudobulk.R`  
###6. `ovary_pseudobulk.R`  
Generate pseudobulk expression matrices by aggregating counts across cell type–sample combinations, followed by PCA.

**Uses:**
- `aggregate.Matrix()` to build pseudobulk matrices  
- `DESeq2::rlog()` for variance-stabilizing transformation  
- `ggplot2` for PCA plots colored by cell type and shaped by sample  

**Outputs:**
- PCA plots (`.png`) of pseudobulk expression profiles by cell type  
- Intermediate `.csv` tables for metadata inspection (optional)
