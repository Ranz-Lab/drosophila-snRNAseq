# 02_parse_pipeline

This directory contains SLURM batch scripts for running the Parse Biosciences split-pipe workflow on single-nucleus RNA-seq data from testis and ovary samples of three *Drosophila* strains: *D. melanogaster* (ISO1, A4) and *D. simulans* (w501).

The pipeline includes demultiplexing, alignment, and combining sub-libraries for each strain.

---

## Overview of Workflow

1. **Demultiplexing**  
   `01_demultiplex_by_strain.sub`
   Demultiplexing was performed by **strain**, using Parse Biosciences' in-house script `fastq_sep_groups.py`. The sublibraries were processed individually and wells were grouped into three strains according to the layout defined in `01_parse_input_prep/sample_sheet.csv`.
   See `01_demultiplex_by_strain.sub` for the full script.

2. **Alignment**  
   `02_parse_align_<strain>.sub`  
   - Aligns each strain separately to strain-specific Parse reference directories prepared in `00_prepare_references`.
   - Employs the Parse Biosciences pipeline (split-pipe v1.1.2) using `split-pipe --mode all`.

3. **Combining Sublibraries**  
   `03_parse_combine_<strain>.sub`  
   - Merges all sub-libraries within each strain into a single combined dataset.
   - Outputs include unfiltered gene expression matrices, cell metadata, and gene annotation tables.

---

## Input Files

- `parfile.txt`: Parameter file specifying experiment-level settings for Parse demultiplexing and alignment.
- `sample_sheet_<strain>.csv`: Strain-specific sample sheets defining sublibraries and sample barcodes.

---

## Notes

- Demultiplexing and alignment are performed separately for each strain to enable species-specific alignment and filtering.
- The Parse WT v2 (100k nuclei) kit was used for this experiment. As recommended by Parse Biosciences, the `--kit WT_mega` argument was used during demultiplexing and alignment to ensure compatibility with the updated barcode structure in v2 kits.
- The output from this stage is used as input for cross-species Seurat object construction in `01_cross_species_integration/`.
