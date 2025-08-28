#!/bin/bash

#SBATCH --job-name=guppy
#SBATCH --account=tyjames1
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=4g
#SBATCH --time=12:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python
module load R
module load guppygpu
module load cuda

# Make directories
mkdir -p $project_dir/seqs/rdna/fast5 $project_dir/seqs/rdna/guppy $project_dir/seqs/rdna/summary

# Determine configuration file
# Apart from the FAST5 input files, Guppy requires a configuration file based on the flowcell and library preparation kit.
# The individual flowcell and kit can be specified or the corresponding config file (with the .cfg prefix) can be specified.
# This script assumes the flowcell FLO-MIN106 and the kit SQK-LSK109 were used. This corresponds to the configuration file dna_r9.4.1_450bps_hac.cfg.
# To list all supported flowcell and library kits type, run:
guppy_basecaller --print_workflows

# Run Guppy
# Relevant options:
#   --config = configuration file based on flowcell and library preparation kit
#   --flowcell = flow cell name (only required if config file not specified)
#   --kit = kit name (only required if config file not specified)
#   --trim_adapters = trim the adapters from the sequences in the output files
#   --trim_strategy = trimming strategy to apply: 'dna' or 'rna' (or 'none' to disable trimming)
#   --min_qscore = minimum acceptable qscore for a read to be filtered into the PASS folder
guppy_basecaller -i $project_dir/seqs/fast5 -s $project_dir/seqs/guppy/simplex_calls \
--config dna_r9.4.1_450bps_sup.cfg -x 'auto' --trim_adapters --trim_strategy dna

# Find duplex reads
duplex_tools pairs_from_summary $project_dir/seqs/guppy/simplex_calls/sequencing_summary.txt $project_dir/seqs/guppy/simplex_calls/pairs
duplex_tools filter_pairs $project_dir/seqs/guppy/simplex_calls/pairs/pair_ids.txt $project_dir/seqs/guppy/simplex_calls/pass

# Run Guppy duplex basecaller
guppy_basecaller_duplex -i $project_dir/seqs/fast5 -r -s $project_dir/seqs/guppy/duplex_calls \
--config dna_r9.4.1_450bps_sup.cfg -x 'auto' \
--duplex_pairing_mode from_pair_list --duplex_pairing_file $project_dir/seqs/guppy/simplex_calls/pairs/pair_ids_filtered.txt

# Show read count
cat $project_dir/seqs/guppy/simplex_calls/pass/*runid*.fastq > $project_dir/seqs/guppy/simplex_calls/pass/basecall.fastq
cat $project_dir/seqs/guppy/simplex_calls/pass/basecall.fastq | wc -l | awk '{print $1/4}' > $project_dir/seqs/rdna/summary/simplex_reads.txt
cat $project_dir/seqs/guppy/duplex_calls/pass/*runid*.fastq > $project_dir/seqs/guppy/duplex_calls/pass/basecall.fastq
cat $project_dir/seqs/guppy/duplex_calls/pass/basecall.fastq | wc -l | awk '{print $1/4}' > $project_dir/seqs/rdna/summary/duplex_reads.txt
