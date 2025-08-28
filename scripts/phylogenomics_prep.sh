#!/bin/bash

# Set variables
list="$project_dir/data/genome_accessions.txt"
genomes_dir="$project_dir/seqs/phylogenomics"
mkdir -p $genomes_dir

# Format list
# The sample file needs to end in a new line in order for the last sample to be read in the while loop
dos2unix $list
function file_ends_with_newline() {
    [[ $(tail -c1 "$1" | wc -l) -gt 0 ]]
}
if ! file_ends_with_newline $list; then
    echo "" >> $list
fi

# Download genomes
cd $genomes_dir
while IFS=$'\t' read -r -a samples_array
do
    name="$(echo ${samples_array[1]} | tr " " "_" | sed 's/[_[:space:]]*$//')"
    strain="${samples_array[2]}"
    isolate="${samples_array[3]}"
    accession="${samples_array[5]}"
    if [ "$strain" == "NA" ]; then
        id=$isolate
    else
        id=$strain
    fi
    datasets download genome accession $accession --include genome
    unzip ncbi_dataset.zip
    mv $genomes_dir/ncbi_dataset/data/$accession/*.fna $genomes_dir/"$name"_"$id".fasta
    rm -r $genomes_dir/ncbi_dataset $genomes_dir/ncbi_dataset.zip $genomes_dir/README.md $genomes_dir/md5sum.txt
done < <(tail -n +2 $list) # skip first line (column headers)
