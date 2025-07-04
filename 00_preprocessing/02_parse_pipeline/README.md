# Parse Demultiplexing and Alignment

This folder contains scripts used to demultiplex and align Parse Biosciences WT_v2 single-cell RNA-seq data across three Drosophila strains.

## Demultiplexing

Demultiplexing was performed by **strain**, using `fastq_sep_groups.py`. The sublibraries were processed individually and wells were grouped into three strains according to the layout defined in `01_parse_input_prep/sample_sheet.csv`.

See `01_demultiplex_by_strain.sh` for the full script.

## Alignment

Each strain was aligned separately to its own strain-specific reference genome using `split-pipe --mode all`. The following scripts are included:

- `02_parse_align_A4.sh`
- `02_parse_align_ISO1.sh`
- `02_parse_align_w501.sh`

Each script defines the appropriate input FASTQ directory (`--fq_dir`), reference index directory (`--ref_dir`), and output path.

Parse reference indices were generated earlier in the pipeline using `split-pipe --mode mkref`.

## Kit Version Note

The Parse WT v2 (100k cell) kit was used for this experiment. As recommended by Parse Biosciences, the `--kit WT_mega` argument was used during demultiplexing and alignment to ensure compatibility with the updated barcode structure in v2 kits. This was critical for accurate parsing of sublibrary and well barcodes.
