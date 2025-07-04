# 01_prepare_references

This directory contains SLURM batch scripts for preparing genome references and annotations for the three *Drosophila* strains used in the study: *D. melanogaster* ISO1, *D. melanogaster* A4, and *D. simulans* w501. These references are used as input for the Parse Biosciences split-pipe alignment workflow.

---

## Workflow Overview

For each strain, we:

1. **Download and prepare genome FASTA and GTF files**  
2. **Filter annotations to include only protein-coding genes**
3. **Build Parse-compatible references using `split-pipe --mode mkref`**
4. **Use Liftoff to lift ISO1 annotations onto A4**

---

## Script Descriptions

- `liftoff_annotation_A4.sub`: Maps ISO1 GTF annotations to the A4 genome using [Liftoff](https://github.com/agshumate/Liftoff)
- `mkgtf_filter_protein_coding.sub`: Filters any GTF file to retain `gene_biotype == "protein_coding"` entries using Cell Ranger
- `mkref_splitpipe_ISO1.sub`: Builds the Parse reference for ISO1 from BDGP6.46
- `mkref_splitpipe_A4.sub`: Uses Liftoff annotations and builds the A4 reference
- `mkref_splitpipe_w501.sub`: Builds the Parse reference for *D. simulans* w501
---

## Output

Each script produces a strain-specific reference directory compatible with the Parse split-pipe workflow. These references are later used in alignment steps under `00_processing/02_parse_pipeline/`.
