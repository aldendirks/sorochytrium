#!/bin/bash

#SBATCH --job-name=barb
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=12gb
#SBATCH --time=2:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load bbtools

# Set variables
barb_dir="$project_dir/seqs/barbacenia"
barb_dir_trimmed="$project_dir/seqs/barbacenia/trimmed"
output_dir="$project_dir/meta/barb"
map_dir="$output_dir/bbmap"
cd $output_dir

# Map to extract 5.8S
# Map paired reads, separate output
bbmap.sh in1=$barb_dir_trimmed/SRR8113281_trimmed_1P.fastq \
in2=$barb_dir_trimmed/SRR8113281_trimmed_2P.fastq \
out1=$map_dir/SRR8113281_trimmed_1P_mapped.fastq \
out2=$map_dir/SRR8113281_trimmed_2P_mapped.fastq \
mappedonly=t \
ref=$project_dir/data/soro_5.8S.fasta nodisk
### Map unpaired reads
bbmap.sh in=$barb_dir_trimmed/SRR8113281_trimmed_1U.fastq \
out=$map_dir/SRR8113281_trimmed_1U_mapped.fastq \
mappedonly=t \
ref=$project_dir/data/soro_5.8S.fasta nodisk
bbmap.sh in=$barb_dir_trimmed/SRR8113281_trimmed_2U.fastq \
out=$map_dir/SRR8113281_trimmed_2U_mapped.fastq \
mappedonly=t \
ref=$project_dir/data/soro_5.8S.fasta nodisk