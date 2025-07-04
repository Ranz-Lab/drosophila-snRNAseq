# Parse Input Preparation

This folder contains metadata describing the structure of the raw input data before Parse Biosciences demultiplexing.

## Sample and Sublibrary Structure

The experiment was sequenced across 8 sublibraries (S1–S8). Each sublibrary contained multiple wells assigned to different biological samples, grouped by strain and tissue type.

## Sample Sheet

The file `sample_sheet.csv` contains well positions and corresponding sample IDs used during sequencing.

## NCBI Submission Format

Raw FASTQs were uploaded to NCBI with the structure:
- 12 samples × 8 sublibraries
- 96 FASTQ files for each R1 and R2
- Format: `split-fq-S1_A4_1o_R1.fastq.gz`, etc.

## Demultiplexing Strategy

Demultiplexing was performed by **strain**, using `fastq_sep_groups.py`, and groups were defined based on well layouts (see `02_parse_pipeline/01_demultiplex_by_strain.sh`).

If reproducing this analysis, users should first reconstruct the strain-level inputs by combining sublibrary FASTQs using the group layout described in `parse_strain_groupings.tsv`.
