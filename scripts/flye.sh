#!/bin/bash

#SBATCH --job-name=flye
#SBATCH --account=tyjames1
#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50gb
#SBATCH --time=24:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load gcc
module load spades

# Assemble genomes with flye
flye --nano-hq $project_dir/seqs/genome/raw/SRR26060970.fastq.gz -o $project_dir/seqs/genome/flye/hq -t 16
