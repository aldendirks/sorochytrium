#!/bin/bash

#SBATCH --job-name=bbmap
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=12gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,end,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bbtools
module load bedtools2

# Set variables
soro_genome="$project_dir/seqs/genome/funannotate/Sormil/fun_out/annotate_results/Sorochytrium_milnesiophthora_S1865.1.scaffolds.fa"
busco_dir="$project_dir/phylo/phylogenomics/busco/Sorochytrium_milnesiophthora_CZEUM2018/run_fungi_odb10"
illumina_dir="$project_dir/seqs/genome/illumina_trimmed"
R1="$illumina_dir/filtered_1P.fastq.gz"
R2="$illumina_dir/filtered_2P.fastq.gz"
output_dir="$project_dir/its/coverage"
cd $output_dir

# Get control sequences
# These are single-copy genes to serve as a coverage reference.
counter=1
ls $busco_dir/busco_sequences/single_copy_busco_sequences/*.faa | sort -R | tail -10 | while read file; do
    cp $file ./control_${counter}.faa
    awk 'sub(/^>/, "")' control_${counter}.faa | awk ' OFS="\t" { gsub(":","\t");gsub("-","\t");print } ' > control_${counter}.bed
    scaffold=$(cat control_${counter}.bed | cut -f 1)
    start=$(cat control_${counter}.bed | cut -f 2)
    end=$(cat control_${counter}.bed | cut -f 3)
    if [ "$start" -gt "$end" ]; then
        printf "$scaffold\t$end\t$start\n" > control_${counter}.bed
    fi
    bedtools getfasta -fi $soro_genome -fo control_${counter}.fasta -bed control_${counter}.bed
    counter=$((counter+1))
done

# Get coverage
# Genome
bbmap.sh in=$R1 in2=$R2 ref=$soro_genome \
nodisk covstats=genome_covstats.txt covhist=genome_histogram.txt \
basecov=genome_basecov.txt bincov=genome_bincov.txt 2> genome_log.txt
# BUSCO controls
for i in {1..10}
do
    bbmap.sh in=$R1 in2=$R2 ref=control_${i}.fasta \
    nodisk covstats=control_${i}_covstats.txt covhist=control_${i}_histogram.txt \
    basecov=control_${i}_basecov.txt bincov=control_${i}_bincov.txt 2> control_${i}_log.txt
done
# rDNA regions
cp $project_dir/data/soro_*.fasta .
for region in full SSU ITS1 5.8S ITS2 ITS LSU
do
    bbmap.sh in=$R1 in2=$R2 ref=$region.fasta \
    nodisk covstats=${region}_covstats.txt covhist=${region}_histogram.txt \
    basecov=${region}_basecov.txt bincov=${region}_bincov.txt 2> ${region}_log.txt
done

# Summarize findings
# Average coverage for the single-copy genes was calculated in Excel (218.3) and ratios were also calculated relative to that number
printf "subject\tcoverage\tratio\n" > summary.txt
genome_coverage=$(cat genome_log.txt | grep "Average coverage:" | grep -Eo '[0-9]*\.[0-9]*')
for LOG in *log.txt
do
    subject=$(basename ${LOG%_log.txt})
    coverage=$(cat $LOG | grep "Average coverage:" | grep -Eo '[0-9]*\.[0-9]*')
    ratio=$(echo "scale=2 ; $coverage / $genome_coverage" | bc)
    printf "$subject\t$coverage\t$ratio\n" >> summary.txt
done
