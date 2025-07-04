#!/bin/bash

#SBATCH --job-name=demux_parse
#SBATCH -A JRANZ_LAB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --error=%x.%A.err
#SBATCH --output=%x.%A.out

# Load Parse environment
module load anaconda
source activate spipe

# Set paths
SCRIPTPATH="/pub/ihariyan/ranz_lab/parse_analysis/ParseBiosciences-Pipeline.1.1.2/fastq_sep_groups.py"
FQ_DIR="/dfs5/bio/ihariyan/ranz_lab/Parse/combined/"
OUT_DIR="/pub/ihariyan/ranz_lab/data/single-cell/parse/data/"

# Define well groupings per strain
GROUPS="--group ISO1 B5-B8,C9-C12,A5-A8,C1-C4 \
        --group A4 B1-B4,D1-D4,A1-A4,B9-B12 \
        --group w501 D5-D8,D9-D12,A9-A12,C5-C8"

# Demultiplex all 8 sublibraries
for i in {1..8}; do
  fq1="${FQ_DIR}S${i}_S${i}_L001_R1_001.fastq.gz"
  fq2="${FQ_DIR}S${i}_S${i}_L001_R2_001.fastq.gz"
  out_subdir="${OUT_DIR}split-fq-S${i}"

  echo "Demultiplexing sublibrary S${i}..."
  python "$SCRIPTPATH" \
    --kit WT_mega \
    --fq1 "$fq1" \
    --fq2 "$fq2" \
    --opath "$out_subdir" \
    $GROUPS \
    --kit_score_skip
done
