# 02_annotation

This directory contains scripts for cell type annotation and composition analysis of single-nucleus RNA-seq data from *Drosophila* testis and ovary. This stage defines biologically meaningful clusters, visualizes marker gene expression, and confirms cell type progression using trajectory inference. All steps were performed on Harmony-integrated Seurat objects from `01_cross_species_integration`.

---

## Overview

For both tissues (testis and ovary), we:

- Visualized canonical marker genes across clusters using `DotPlot`
- Assigned cell type labels to Harmony-integrated clusters
- Reordered cell types to reflect expected developmental progression
- Confirmed inferred trajectories using Monocle3 pseudotime analysis
- Quantified cell type proportions across strains

---

## Scripts Included

### `testis_annotation.R`

- Loads the integrated Seurat object (`testis.FINAL.collapsed.rds`)
- Visualizes canonical markers (e.g., *bam*, *His2Av*, *Rbp4*, *Fas3*, etc.)
- Assigns cluster identities to known testis cell types, including:
  - Early germline: GSCs, spermatogonia, early spermatocytes
  - Late germline: mid/late spermatocytes, spermatids
  - Somatic: hub, cyst, epithelial cells
- Excludes "Unannotated" clusters
- Saves final annotated object ("FINAL-testis-nounknown.rds")

### `ovary_annotation.R`

- Loads the integrated Seurat object (`ovary.annotated.v3.SCT.rds`)
- Visualizes canonical markers (e.g., *corolla*, *orb*, *cas*, *Sox14*, *Act57B*)
- Assigns cluster identities for germline, follicle, corpus luteum, and muscle/oviduct cells
- Excludes "Unannotated" and ambiguous MBFC clusters with widespread expression
- Saves final annotated object ("FINAL-ovary-nounknown.rds")

---

## Trajectory Inference and Validation

To validate the annotated cell types and confirm developmental order:

- Monocle3 was used to infer pseudotime trajectories from Harmony-reduced Seurat objects
- UMAP coordinates and cluster labels were transferred to Monocle3 objects
- Trajectory graphs were learned and pseudotime was estimated per cell
- Boxplots of pseudotime across annotated cell types were generated to confirm:
  - Germline development from stem cells to mature cells
  - Somatic lineage structure, particularly in the ovary

This served as an internal consistency check to support the biological accuracy of annotations.

---

## Cell Type Composition

### `testis_celltype_composition.R`
### `ovary_celltype_composition.R`

- Quantify and visualize the distribution of annotated cell types across strains
- Uses:
  - `dittoBarPlot` for quick proportion visualization
  - Manual barplots of cell counts and proportions using `ggplot2`
  - `SCpubr::do_BarPlot()` for stylized summary figures
- Outputs per-strain proportions as `.csv` files for downstream statistical comparison
