#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=25gb
#SBATCH --time=1:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load blast-plus

# Set variables
cores="$SLURM_NTASKS"
query=$1
db=$2
out=$3

# BLAST
# outfmt 6 is necessary because outfmt 7 outputs 4 comment lines for EVERY query (so a few million unnecessary lines)
blastn -query $query -db $db -outfmt '6' -evalue 1e-5 -out $out -num_threads $cores