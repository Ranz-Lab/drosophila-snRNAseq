# Parse Demultiplexing and Alignment

This folder contains scripts used to demultiplex and align Parse Biosciences WT single-nucleus RNA-seq data.

## Demultiplexing

Demultiplexing was performed by strain using `fastq_sep_groups.py`. See `01_demultiplex_by_strain.sh`.

## Alignment

Each strain was aligned to its own reference genome using `split-pipe --mode all`. The scripts:

- `02_parse_align_A4.sh`
- `02_parse_align_ISO1.sh`
- `02_parse_align_w501.sh`

Each script specifies the input FASTQ directory (`--fq_dir`), the appropriate reference index (`--ref_dir`), and the output directory.

Parse reference indices were generated using `split-pipe --mode mkref` earlier in the pipeline.
