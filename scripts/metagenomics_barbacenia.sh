#!/bin/bash

# Load modules
module load Bioinformatics
module load bbtools
module load blast-plus
module load seqtk
module load spades
module load sratoolkit
module load trimmomatic

# Set variables
reference="$project_dir/seqs/genome/funannotate/Sormil/fun_out/annotate_results/Sorochytrium_milnesiophthora_S1865.1.scaffolds.fa"
barb_dir="$project_dir/seqs/barbacenia"
barb_dir_raw="$project_dir/seqs/barbacenia/raw"
barb_dir_trimmed="$project_dir/seqs/barbacenia/trimmed"
output_dir="$project_dir/meta/barb"
db_dir="$output_dir/db" # Directory for BLAST databases
blast_dir="$output_dir/blast" # Directory for BLAST output
map_dir="$output_dir/bbmap"
seq_dir="$output_dir/seqs" # Directory for Barbacenia sequences
spades_dir="$output_dir/spades"
mkdir -p $output_dir $db_dir $blast_dir $map_dir $seq_dir/subset $spades_dir
cd $output_dir

# Trim SRA reads
mkdir -p $barb_dir_trimmed
cd $barb_dir_trimmed
TrimmomaticPE -trimlog $barb_dir_trimmed/trimmomatic.log \
-basein $barb_dir_raw/SRR8113281_1.fastq.gz \
-baseout $barb_dir_trimmed/filtered.fastq.gz \
ILLUMINACLIP:$project_dir/data/TruSeq3-PE-adapters.fa:3:30:10 SLIDINGWINDOW:4:26 LEADING:33 TRAILING:24 MINLEN:50

# Convert FASTQ to FASTA
for fastq in $barb_dir_trimmed/*.fastq
do
    name="$(basename ${fastq%.fastq})" # Get file basename without file extension
    seqtk seq -a $fastq > $barb_dir_trimmed/$name.fasta 
done

# Make BLAST database for querying
makeblastdb -in $reference -dbtype nucl -out $db_dir/soro_genome_db
makeblastdb -in $project_dir/data/soro_full.fasta -dbtype nucl -out $db_dir/soro_rdna_db
makeblastdb -in $project_dir/data/soro_SSU.fasta -dbtype nucl -out $db_dir/soro_ssu_db
makeblastdb -in $project_dir/data/soro_5.8s.fasta -dbtype nucl -out $db_dir/soro_58s_db
makeblastdb -in $project_dir/data/soro_LSU.fasta -dbtype nucl -out $db_dir/soro_lsu_db

# BLAST Barbacenia reads against the Sorochytrium databases
# Start a slurm job for each BLAST task (20 in total)
for db in genome rdna ssu 58s lsu
do
    for fasta in $barb_dir_trimmed/*.fasta
    do
        name="$(basename ${fasta%.fasta})"
        sbatch $project_dir/scripts/metagenomics_barbacenia_blast.sh $fasta $db_dir/soro_${db}_db $blast_dir/${name}_${db}_blast.out
    done
done

# Subset FASTQ files
# Based on BLAST output, subset trimmed reads and make files containing only matches
# This takes a little while
for db in genome rdna ssu 58s lsu
do
    mkdir -p $seq_dir/subset/$db
    for fasta in $barb_dir_trimmed/*.fasta
    do
        name="$(basename ${fasta%.fasta})"
        cat $blast_dir/${name}_${db}_blast.out | cut -f 1 | grep -v "#" > $blast_dir/${name}_${db}_blast.txt
        seqtk subseq $barb_dir_trimmed/$name.fastq $blast_dir/${name}_${db}_blast.txt > $seq_dir/subset/$db/${name}_subset.fastq
    done
done

# Match up paired reads
# For each BLAST db type, make an empty file for paired sequence names
# For each paired sequence file, get all the names and add them to this file
# Sort them and get the unique names
# Then subset the trimmed paired reads again
for db in genome rdna ssu 58s lsu
do
    >$seq_dir/subset/$db/paired_seqs.txt
    for fastq in $seq_dir/subset/$db/*P_subset.fastq
    do
        cat $fastq | grep "@SRR" | cut -d " " -f 1 | sed 's/@//' | sort | uniq >> $seq_dir/subset/$db/paired_seqs.txt
    done
    cat $seq_dir/subset/$db/paired_seqs.txt | sort | uniq > $seq_dir/subset/$db/paired_seqs_uniq.txt # Get list of unique read names
    for fastq in $seq_dir/trimmed/*P.fastq
    do
        name="$(basename ${fastq%.fastq})"
        seqtk subseq $fastq $seq_dir/subset/$db/paired_seqs_uniq.txt > $seq_dir/subset/$db/${name}_matched.fastq
    done
done

# Assemble reads
for db in genome rdna ssu 58s lsu
do
    mkdir -p $spades_dir/$db
    gzip --fast $seq_dir/subset/$db/*.fastq
    spades.py --meta -o "$spades_dir/$db" \
    -1 "$seq_dir/subset/$db/"*"1P_matched.fastq.gz" \
    -2 "$seq_dir/subset/$db/"*"2P_matched.fastq.gz" \
    -s "$seq_dir/subset/$db/"*"1U_subset.fastq.gz" \
    -s "$seq_dir/subset/$db/"*"2U_subset.fastq.gz"
done

# BLAST Sorochytrium against assembled scaffolds to identify putativehits
for db in genome rdna ssu 58s lsu
do
    makeblastdb -in $spades_dir/$db/scaffolds.fasta -dbtype nucl -out $spades_dir/$db/scaffolds_db
    > tmp.txt
    for ref in full SSU 5.8S LSU
    do
        blastn -query $project_dir/data/soro_$ref.fasta -db $spades_dir/$db/scaffolds_db -outfmt '6' -evalue 1e-5 -out $spades_dir/$db/scaffolds_${ref}_blast.out
        awk -F '\t' -v threshold1="90" -v threshold2="100" \
        '$3 >= threshold1 && $4 >= threshold2 {print $0}' \
        $spades_dir/$db/scaffolds_${ref}_blast.out > $spades_dir/$db/scaffolds_${ref}_blast_select.out 
        cat $spades_dir/$db/scaffolds_${ref}_blast_select.out | cut -f 2 >> tmp.txt
    done
    seqtk subseq $spades_dir/$db/scaffolds.fasta <(cat tmp.txt | sort | uniq) > $spades_dir/$db/scaffolds_select.fasta
    rm tmp.txt
done

# Pull out manually selected sequences
>$spades_dir/seqs_select.fasta
seqtk subseq $spades_dir/genome/scaffolds.fasta <(echo "NODE_772_length_284_cov_2.222707") >> $spades_dir/seqs_select.fasta
seqtk subseq $spades_dir/ssu/scaffolds.fasta <(echo "NODE_26_length_248_cov_0.989637") >> $spades_dir/seqs_select.fasta
seqtk subseq $spades_dir/lsu/scaffolds.fasta <(echo "NODE_10_length_508_cov_1.860927"; echo "NODE_47_length_268_cov_1.845070") >> $spades_dir/seqs_select.fasta
seqtk subseq $spades_dir/rdna/scaffolds.fasta <(echo "NODE_66_length_268_cov_1.845070"; echo "NODE_92_length_248_cov_0.989637") >> $spades_dir/seqs_select.fasta

# Map to extract 5.8S
sbatch $project_dir/scripts/metagenomics_barbacenia_map.sh

# Assemble reads
gzip $map_dir/*.fastq
mkdir -p $spades_dir/mapped
spades.py -o "$spades_dir/mapped" \
-1 "$map_dir/SRR8113281_trimmed_1P_mapped.fastq.gz" \
-2 "$map_dir/SRR8113281_trimmed_2P_mapped.fastq.gz" \
-s "$map_dir/SRR8113281_trimmed_1U_mapped.fastq.gz" \
-s "$map_dir/SRR8113281_trimmed_2U_mapped.fastq.gz" \
--careful -k 21
