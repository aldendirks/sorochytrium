#!/bin/bash

#SBATCH --job-name=quast_busco
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=5gb 
#SBATCH --time=12:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python

# Set variables
genome="$project_dir"/seqs/genome/hapog_2/hapog_results/hapog.fasta
quast_dir="$project_dir"/seqs/genome/quast
mkdir -p $quast_dir
cd $quast_dir

# Run quast
python quast.py -o hapog_2 --fungus --rna-finding --split-scaffolds $genome


