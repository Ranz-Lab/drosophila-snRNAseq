# 00_prepare_references

This directory contains SLURM batch scripts for preparing genome references and annotations for the three *Drosophila* strains used in the study: *D. melanogaster* ISO1, *D. melanogaster* A4, and *D. simulans* w501. These references are used as input for the Parse Biosciences split-pipe alignment workflow.

---

## Overview of Workflow

For each strain, we:

1. **Download and prepare genome FASTA and GTF files**  
2. **Filter annotations to include only protein-coding genes**
3. **Build Parse-compatible references using `split-pipe --mode mkref`**
4. **Use Liftoff to lift ISO1 annotations onto A4**

---

## Scripts Included

### 1. `mkref_splitpipe_ISO1.sub`
- Prepares the ISO1 reference using the original BDGP6.46 genome and GTF.
- Filters for protein-coding genes using `mkgtf_filter_protein_coding.sub`
- Builds a reference with the Parse `mkref` tool.

### 2. `mkref_splitpipe_A4.sub`
- Uses Liftoff to map ISO1 GTF annotations onto the A4 genome.
- Filters the lifted GTF to retain only protein-coding genes.
- Builds the A4 Parse reference using the lifted annotations.

### 3. `mkref_splitpipe_w501.sub`
- Uses the w501 assembly and accompanying annotation GTF.
- Filters for protein-coding genes.
- Builds the Parse reference for *D. simulans* w501.

### 4. `liftoff_annotation_A4.sub`
- Performs annotation transfer using [Liftoff](https://github.com/agshumate/Liftoff), mapping gene models from ISO1 to A4.
- Ensures gene IDs are consistent across *D. melanogaster* strains.

### 5. `mkgtf_filter_protein_coding.sub`
- Filters any GTF file to retain only `gene_type == "protein_coding"` entries.
- Used as a preprocessing step before `mkref`.

---

## Output

Each script produces a strain-specific reference directory compatible with the Parse split-pipe workflow. These references are later used in demultiplexing and alignment steps under `00_processing/01_demultiplex/`.
