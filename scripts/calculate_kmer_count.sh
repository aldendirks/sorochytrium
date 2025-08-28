#!/bin/bash

#SBATCH --job-name=kmer
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20gb 
#SBATCH --time=1:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bbtools

# Run bbtools
illumina_dir="$project_dir"/seqs/genome/illumina_trimmed
output_dir="$project_dir"/seqs/genome/kmers
mkdir $output_dir
cd $output_dir
kmercountexact.sh in=$illumina_dir/filtered_1P.fastq.gz in2=$illumina_dir/filtered_2P.fastq.gz khist=hist.txt k=23 histcolumns=2 gchist=t peaks=peaks.txt
