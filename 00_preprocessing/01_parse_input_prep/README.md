# 00_parse_input_prep

This folder contains metadata describing the structure of the raw input data before Parse Biosciences demultiplexing.

## Sample and Sublibrary Structure

The experiment was sequenced across 8 sublibraries (S1–S8). Each sublibrary contained all 48 wells corresponding to 12 biological samples (2 replicates × 2 tissues × 3 strains). Each biological sample was distributed across 4 wells.

## Sample Sheet

The file `sample_sheet.csv` contains:
- A unique sample name (`Sample_Name`)
- Strain, tissue type, and replicate number
- The 4-well group (`Wells`) corresponding to each sample

## NCBI Submission Format

Raw FASTQ files were demultiplexed by **sample** and uploaded to NCBI in the format:
- 12 samples × 8 sublibraries = 96 FASTQ files for each R1 and R2
- File Format: `split-fq-S1_A4_1o_R1.fastq.gz`, etc.

Each file corresponds to a sublibrary–sample combination for one read direction.

## Demultiplexing Strategy

For analysis purposes, demultiplexing was performed by **strain**, using `fastq_sep_groups.py`, grouping wells according to the layout shown in `sample_sheet.csv`. The script used for this can be found in `02_parse_pipeline/01_demultiplex_by_strain.sub`.

To reproduce the analysis, combine FASTQs from the appropriate sublibraries by sample and strain using the well groupings defined in `sample_sheet.csv`.
