#!/bin/bash

#SBATCH --job-name=extractITS
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=10gb 
#SBATCH --time=10:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bedtools2
module load blast-plus
module load python

# Set variables
reference="$project_dir/seqs/genome/funannotate/Sormil/fun_out/annotate_results/Sorochytrium_milnesiophthora_S1865.1.scaffolds.fa"
output_dir="$project_dir/its/extract"
mkdir -p $output_dir
cd $output_dir

# Run software
python extractITS.py -i $reference -o $output_dir -which all -name Sorochytrium
