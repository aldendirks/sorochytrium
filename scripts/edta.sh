#!/bin/bash

#SBATCH --job-name=edta
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=5gb
#SBATCH --time=1:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load singularity

# Set variables
cores="$SLURM_NTASKS"
output_dir="$project_dir/seqs/genome/edta"
genome="$project_dir/seqs/genome/funannotate/Sormil/fun_out/annotate_results/Sorochytrium_milnesiophthora_S1865.1.scaffolds.fa"
mkdir -p $output_dir
cd $output_dir

# Run edta
# Issue resolved here: https://github.com/oushujun/EDTA/issues/568#issuecomment-2913239008
unset -f which
singularity exec $EDTA_PATH/EDTA.sif EDTA.pl \
--genome $genome --species others --force 1 --anno 1 --evaluate 1 --sensitive 1 -t $cores 
