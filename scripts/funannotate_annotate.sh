#!/bin/bash

#SBATCH --job-name=fun_annotate
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=50gb
#SBATCH --time=4:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Set variables
cores="$SLURM_NTASKS"
species="Sorochytrium milnesiophthora"
isolate="S1865.1"
gbk=$project_dir/seqs/genome/antismash7/Sormil_funmasked.gbk
sbt=$project_dir/seqs/genome/funannotate/Sormil/template.sbt
cd $project_dir/seqs/genome/funannotate/Sormil

# Run interproscan
funannotate iprscan -i fun_out -m local --cpus $(($cores-2))

# Run funannotate
# If the interproscan file is in fun_out (as it shouild be), do not specify the file name in the command. It will find it automatically. 
# If you do specify it, it deletes the file for some reason. 
funannotate annotate -i fun_out --antismash $gbk --cpus $cores --tmpdir $TMPDIR --sbt $sbt
