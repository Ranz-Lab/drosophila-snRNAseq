# 01_cross_species_integration

This folder documents the process of harmonizing gene annotations across species, constructing strain-specific Seurat objects, and integrating testis and ovary single-nucleus RNA-seq data from *Drosophila melanogaster* and *Drosophila simulans* for comparative analysis.

---

## Overview

We analyzed testis and ovary nuclei from three Drosophila strains:
- *D. melanogaster* (ISO1 and A4)
- *D. simulans* (w501)

Because gene annotations differ between species, we first converted *D. simulans* gene IDs to their *D. melanogaster* orthologs using a curated mapping file. These standardized gene IDs enabled meaningful cross-species integration and downstream analyses. Importantly, any *D. simulans*-specific genes or those lacking ortholog mappings were retained.

---

## Gene ID Harmonization

### Files:
- `mel-sim_orthologs.txt`:  
  A two-column tab-delimited file containing ortholog mappings: Used to convert *D. simulans* gene IDs into *D. melanogaster* orthologs.

- `replace_gene_ortholog.py`:  
  A Python script that replaces matching *D. simulans* gene names in `all_genes_w501.csv` with *D. melanogaster* orthologs and retains genes without orthologs (species-specific or unannotated)

---

## Parse Biosciences DGE Input Files

From the Parse Biosciences pipeline, we downloaded the following unfiltered files for each sample:
- `count_matrix.mtx`
- `all_genes.csv`
- `cell_metadata.csv`

For *D. simulans* (w501), `all_genes.csv` was replaced with the ortholog-converted `all_genes_orthologs+original.csv` before processing.

---

## Seurat Object Construction

Each strain’s Seurat object was created and processed using a standardized pipeline:

### Key steps:
1. Read count matrix, gene list, and metadata
2. Apply stringent QC filters:
 - `nFeature_RNA ≥ 200`
 - `nCount_RNA between 300 and 100000`
 - `nFeature_RNA < 7500`
 - `percent.mt ≤ 5`
3. Save violin plots for basic QC metrics
4. Construct Seurat object with normalized and transposed count matrix

---

## Doublet Detection and Removal

Doublets were identified using the **DoubletFinder** package:

- Each strain was analyzed independently
- pK values were estimated using param sweeps
- A 3% doublet rate was assumed
- Homotypic doublet proportions were modeled using pre-cluster annotations
- Only singlet cells were retained for downstream steps

---

## Normalization and Individual UMAPs

- Each strain was normalized using **SCTransform (v2)**  
- PCA, clustering, and UMAP were computed and saved for each strain
- UMAPs were used for pre-integration quality inspection

---

## Integration with Harmony

- All three testis datasets (ISO1, A4, w501) were merged using Seurat v5’s layer-aware `merge()`
- Datasets were normalized, variable features selected, and scaled
- PCA was computed followed by batch correction using **Harmony**
- Final dimensionality reductions included:
- UMAP (`umap.harmony`)
- tSNE (`tsne.harmony`)

---

## Output

Final integrated object: This file contains the integrated testis dataset across species and strains, ready for annotation, trajectory inference, and cross-species comparison.

