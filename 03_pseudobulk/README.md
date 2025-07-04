# 03_pseudobulk

This directory contains scripts for pseudobulk analysis of single-nucleus RNA-seq data from *Drosophila* testis and ovary. Gene expression counts are aggregated by cell type and sample to enable bulk-like comparisons across strains, followed by PCA for visualization.

---

### 5. `testis_pseudobulk.R`  
### 6. `ovary_pseudobulk.R`  
Aggregate counts across cell type–sample combinations and perform PCA on normalized expression values.

**Uses:**
- `aggregate.Matrix()` to generate pseudobulk expression matrices  
- `DESeq2::rlog()` for variance-stabilizing transformation  
- `ggplot2` for PCA visualization, colored by cell type and shaped by sample  

**Output:**
- PCA plots (`.png`) of cell-type-level pseudobulk profiles
