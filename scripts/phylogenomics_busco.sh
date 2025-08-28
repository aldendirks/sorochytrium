#!/bin/bash

#SBATCH --job-name=busco
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=5gb 
#SBATCH --time=4:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load busco

# Set variables
sample=$1
name=$(basename ${sample%.fasta})

# Run busco
cd $project_dir/phylo/phylogenomics/busco
busco -m genome -i "$sample" -o "$name" -l $project_dir/phylo/phylogenomics/busco/busco_downloads/lineages/fungi_odb10/ --offline