#!/bin/bash

#SBATCH --job-name=lichen_metagenomics
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50gb
#SBATCH --time=12:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load blast-plus
module load python
module load seqtk
module load spades
module load sratoolkit
module load trimmomatic

# Set variables
cores="$SLURM_NTASKS"
reference="$project_dir/seqs/genome/funannotate/Sormil/fun_out/annotate_results/Sorochytrium_milnesiophthora_S1865.1.scaffolds.fa"
lichen_dir="$project_dir/seqs/parmelia"
lichen_dir_raw="$project_dir/seqs/parmelia/raw"
lichen_dir_trimmed="$project_dir/seqs/parmelia/trimmed"
output_dir="$project_dir/meta/parmelia"
db_dir="$output_dir/blastdb"
blast_dir="$output_dir/blast-out"
fastq_dir="$output_dir/fastq"
mkdir -p $output_dir $db_dir $blast_dir $fastq_dir
cd $output_dir

# Trim reads
mkdir -p $lichen_dir_trimmed
cd $lichen_dir_trimmed
TrimmomaticPE -trimlog $lichen_dir_trimmed/trimmomatic.log \
-basein $lichen_dir_raw/SRR34964407_1.fastq.gz \
-baseout $lichen_dir_trimmed/filtered.fastq.gz \
ILLUMINACLIP:$project_dir/data/TruSeq3-PE-adapters.fa:3:30:10 SLIDINGWINDOW:4:26 LEADING:33 TRAILING:24 MINLEN:50

# Make Sorochytrium BLAST database
makeblastdb -in $reference -dbtype nucl -out $db_dir/Soro_db

# BLAST reads and subset
# This loop BLASTs the two paired read files and the two unpaired read files
trimmed_fastq="$lichen_dir_trimmed/"*".fastq"
for fastq in $trimmed_fastq
do
    # Get file basename without file extension
    name="$(basename ${fastq%.fastq})"

    # Convert fastq to fasta
    seqtk seq -a $fastq > $blast_dir/tmp.fasta 

    # BLAST the trimmed reads against the reference genome
    blastn -query $blast_dir/tmp.fasta -db $db_dir/Soro_db -outfmt '6' -evalue 1e-5 \
    -out $blast_dir/${name}_blast.out -num_threads $cores

    # Get read names
    cat $blast_dir/${name}_blast.out | cut -f 1 > $blast_dir/${name}_blast.txt

    # Subset fastq
    seqtk subseq $fastq $blast_dir/${name}_blast.txt > $fastq_dir/${name}_subset.fastq

    # Remove temporary file
    rm $blast_dir/tmp.fasta
done

# Match up paired reads

# Make temporary file
> $fastq_dir/tmp.txt

# Add all fasta headers to temporary file
# Sometimes "@" is a quality symbol so need to add more information to grep
for fastq in $fastq_dir/*P_subset.fastq
do
    cat $fastq | grep "@SRR" | cut -d " " -f 1 | sed 's/@//' | sort | uniq >> $fastq_dir/tmp.txt
done

# Get list of unique headers
cat $fastq_dir/tmp.txt | sort | uniq > $fastq_dir/all_reads.txt

# Reextract sequences from fastq files 
seqtk subseq $(echo $trimmed_fastq | tr " " "\n" | grep "1P") $fastq_dir/all_reads.txt > $fastq_dir/Sample_8120-AD-1-filtered_Soro_1P_matched.fastq
seqtk subseq $(echo $trimmed_fastq | tr " " "\n" | grep "2P") $fastq_dir/all_reads.txt > $fastq_dir/Sample_8120-AD-1-filtered_Soro_2P_matched.fastq

# Remove files
rm $fastq_dir/tmp.txt $fastq_dir/all_reads.txt

# Compress files
gzip --fast $fastq_dir/*.fastq

# Assemble with spades
mkdir -p "$output_dir/spades"
spades.py --meta -o "$output_dir/spades" \
-1 "$fastq_dir/"*"1P_matched.fastq.gz" \
-2 "$fastq_dir/"*"2P_matched.fastq.gz" \
-s "$fastq_dir/"*"1U_subset.fastq.gz" \
-s "$fastq_dir/"*"2U_subset.fastq.gz"

# Run QUAST
mkdir -p "$output_dir/quast/" 
python quast.py -o $output_dir/quast/ --fungus --rna-finding --split-scaffolds $output_dir/spades/scaffolds.fasta
