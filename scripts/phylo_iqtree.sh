#!/bin/bash

#SBATCH --job-name=iqtree2
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics

# Set variables
cores="$SLURM_NTASKS"
iqtree_dir="$project_dir/phylo/rdna/iqtree"
aln="seqs_formatted_combined_aligned_trimmed_formatted.fasta"
mkdir -p $iqtree_dir/alignments
cd $iqtree_dir

# Individual rDNA region phylogenies
for dir in ssu its lsu
do
    cp $project_dir/phylo/rdna/$dir/$aln $iqtree_dir/alignments/alignment_$dir.fasta
    iqtree3 -s $iqtree_dir/alignments/alignment_$dir.fasta  --alrt 1000 -B 1000 --prefix $dir -bb 1000 -T $cores -redo
done

# Concatenated phylogeny
iqtree3 -p $iqtree_dir/alignments --alrt 1000 -B 1000 --prefix concat -T $cores -redo
