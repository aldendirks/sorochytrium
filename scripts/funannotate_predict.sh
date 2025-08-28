#!/bin/bash

#SBATCH --job-name=fun_predict
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=50gb
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load funannotate

# Set variables
cores="$SLURM_NTASKS"
genome="$project_dir"/seqs/genome/hapog_2/hapog_results/hapog.fasta
species="Sorochytrium milnesiophthora"
isolate="S1865.1"
mkdir $project_dir/seqs/genome/funannotate
cd $project_dir/seqs/genome/funannotate

# Prep genome and mask repetitive regions
cp $genome ./Sormil.fasta
funannotate clean -i Sormil.fasta -o Sormil_cleaned.fasta
funannotate sort -i Sormil_cleaned.fasta -o Sormil_cleaned_sorted.fasta
funannotate mask -i Sormil_cleaned_sorted.fasta -o Sormil_funmasked.fasta --cpus $cores
rm Sormil_cleaned.fasta Sormil_cleaned_sorted.fasta

# Predict genes
funannotate predict -i Sormil_funmasked.fasta -o fun_out \
--species "$species" --isolate "$isolate" \
--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
--cpus $cores \
--optimize_augustus