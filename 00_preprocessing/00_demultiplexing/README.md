# 00_demultiplexing

This folder documents the sublibrary structure and includes scripts for demultiplexing raw FASTQ files.

## Sample and Sublibrary Structure

The experiment was sequenced across 8 sublibraries (S1–S8). Each sublibrary contained all 48 wells corresponding to 12 biological samples (2 replicates × 2 tissues × 3 strains) after split-pool combinatorial barcoding. Each biological sample was distributed across 4 wells.

## Sample Sheet

The file `sample_sheet.csv` contains:
- A unique sample name (`Sample_Name`)
- Strain, tissue type, and replicate number
- The 4-well group (`Wells`) corresponding to each sample

## Demultiplexing

For analysis purposes, demultiplexing was performed by **strain** using Parse Biosciences' in-house script `fastq_sep_groups.py`, based on the layout defined in `sample_sheet.csv`. Each sublibrary was processed individually. See `demultiplex_by_strain.sub` for the full SLURM script.

## NCBI Submission Format

Raw FASTQ files were demultiplexed by **sample** and uploaded to NCBI in the following format:
- 12 samples × 8 sublibraries = 96 FASTQ files for each R1 and R2
- File naming format: `split-fq-S1_A4_1o_R1.fastq.gz`, etc.

Each file corresponds to a sublibrary–sample combination for one read direction. To reproduce the analysis, group FASTQ files by strain for each sublibrary using the well assignments provided in `sample_sheet.csv`.
