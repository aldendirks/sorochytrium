#!/bin/bash

# Load modules
module load Bioinformatics
module load R

# Define variables
illumina_dir="$project_dir"/seqs/genome/illumina_trimmed
gscope_dir="$project_dir"/seqs/genome/genomescope

# Run jellyfish
# According to GenomeScope GitHub, also use the -C flag
cd $gscope_dir
jellyfish count -C -m 23 -s 1000000000 -t 32 $illumina_dir/*fastq -o reads.jf
jellyfish histo -t 32 reads.jf > reads.histo

# Run genomescope
# Information: https://github.com/schatzlab/genomescope
Rscript genomescope.R reads.histo 23 151 .
