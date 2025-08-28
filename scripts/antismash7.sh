#!/bin/bash

#SBATCH --job-name=antismash
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=20gb
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Set variables
cores="$SLURM_NTASKS"
species="Sorochytrium milnesiophthora"
isolate="S1865.1"
genome_dir="$project_dir"/seqs/genome
gff3="$(ls $genome_dir/funannotate/Sormil/fun_out/predict_results/*.gff3)"
fasta="$(ls $genome_dir/funannotate/Sormil/*funmasked.fasta)"
mkdir -p $genome_dir/antismash7
cd $genome_dir/antismash7

# Run antismash
antismash --taxon fungi --cassis --clusterhmmer --tigrfam --asf --cc-mibig --cb-general --cb-subclusters --cb-knownclusters --pfam2go --rre --smcog-trees --tfbs \
--output-dir . --cpus $cores --genefinding-gff3 $gff3 $fasta

# Convert files
# Convert GenBank output files to protein FASTA files
for file in *.gbk
do
    header="$species $isolate ${file%.region00*.gbk}"
    header="${header// /_}"
    python "$project_dir"/scripts/gbk_converter.py $file $header
done
file_name="$species $isolate"
file_name="${file_name// /_}"
cat *region*.fasta > "$file_name"_antismash.fasta
