#!/bin/bash

#SBATCH --job-name=fastqc_trimmomatic
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=2:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load fastqc
module load trimmomatic

# Make directories
mkdir -p $project_dir/seqs/genome/fastqc/raw $project_dir/seqs/genome/fastqc/trimmed $project_dir/seqs/genome/illumina_trimmed

# Run fastqc on raw reads
fastqc -o $project_dir/seqs/genome/fastqc/raw $project_dir/seqs/genome/raw/SRR26060971*

# Run trimmomatic
cd $project_dir/seqs/genome/raw
TrimmomaticPE \
-trimlog $project_dir/seqs/genome/illumina_trimmed/trimmomatic.log \
-basein $project_dir/seqs/genome/raw/SRR26060971_1.fastq.gz \
-baseout $project_dir/seqs/genome/illumina_trimmed/filtered.fastq.gz \
ILLUMINACLIP:$project_dir/data/TruSeq3-PE-adapters.fa:3:30:10 SLIDINGWINDOW:4:26 LEADING:33 TRAILING:24 MINLEN:50

# Run fastqc on trimmed reads
fastqc -o $project_dir/seqs/genome/fastqc/trimmed $project_dir/seqs/genome/illumina_trimmed/*
