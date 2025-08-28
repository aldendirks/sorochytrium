#!/bin/bash

#SBATCH --job-name=polish-phase
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --mem=50gb
#SBATCH --time=1:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python
module load bwa
module load samtools

# Set variables
genome_dir="$project_dir"/seqs/genome
assembly="$genome_dir"/flye/hq/assembly.fasta
R1="$genome_dir"/illumina_trimmed/filtered_1P.fastq.gz
R2="$genome_dir"/illumina_trimmed/filtered_2P.fastq.gz

# Run nextpolish
# Information: https://github.com/Nextomics/NextPolish?tab=readme-ov-file
mkdir -p "$genome_dir"/nextpolish
cd "$genome_dir"/nextpolish
ls $R1 $R2 > sgs.fofn
cat > run.cfg << EOL
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = $assembly #genome file
genome_size = auto
workdir = ./01_rundir
polish_options = -p {multithread_jobs}
[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa
EOL
nextPolish run.cfg

# Run Hapo-G (round 1)
# Information: https://github.com/institut-de-genomique/HAPO-G
# The output directory cannot already exist
np_genome="$genome_dir"/nextpolish/01_rundir/genome.nextpolish.fasta
python hapog.py --genome $np_genome --pe1 $R1 --pe2 $R2 -o $genome_dir/hapog_1 -t 36 -u

# Run Hapo-G (round 2)
hp_genome="$genome_dir"/hapog_1/hapog_results/hapog.fasta
python hapog.py --genome $hp_genome --pe1 $R1 --pe2 $R2 -o $genome_dir/hapog_2 -t 36 -u
