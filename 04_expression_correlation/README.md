# 04_expression_correlation

This directory contains scripts to assess transcriptomic similarity across cell types and strains based on average gene expression profiles in the testis and ovary.

---

### `testis_expression_correlation.R`  
### `ovary_expression_correlation.R`  
Calculate correlation matrices of average expression values to quantify cross-strain similarity for each annotated cell type. Filtering is used to reduce the impact of outliers, and cell types with low gene coverage are excluded.

**Uses:**
- `AverageExpression()` to compute mean expression per cell type and strain  
- Pearson correlation coefficients, computed with and without outlier filtering  
- `ggcorrplot` and `ggplot2` for heatmap visualization  

**Outputs:**
- Correlation heatmaps (`.png`) showing cell type similarity within and across strains  
- Gene count tables used in each pairwise comparison (`.csv`)
