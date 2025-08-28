#!/bin/bash

# Set variables
list="$project_dir/data/rdna_accessions.txt"
SSU_dir="$project_dir/phylo/rdna/ssu"
ITS_dir="$project_dir/phylo/rdna/its"
LSU_dir="$project_dir/phylo/rdna/lsu"
longread_dir="$project_dir/phylo/rdna/long-read"
concat_dir="$project_dir/phylo/rdna/concatenated"

# Format list
# The sample file needs to end in a new line in order for the last sample to be read in the while loop
dos2unix $list
function file_ends_with_newline() {
    [[ $(tail -c1 "$1" | wc -l) -gt 0 ]]
}
if ! file_ends_with_newline $list
then
    echo "" >> $list
fi

# Organize sequences

# Make empty files to hold GenBank accessions and header information
mkdir $SSU_dir $ITS_dir $LSU_dir $longread_dir $concat_dir
for dir in $SSU_dir $ITS_dir $LSU_dir $longread_dir
do
    >$dir/seqs.txt
    >$dir/headers.txt
done

#Populate files with GenBank accessions and header information
while IFS=$'\t' read -r -a list_array
do
    if [ "${list_array[9]}" != "NA" ]; then
        seq="${list_array[9]}"
        header="${list_array[13]}"
        echo "$seq" >> $longread_dir/seqs.txt
        echo "$header" >> $longread_dir/headers.txt
    else
        if [ "${list_array[6]}" != "NA" ]; then
            seq="${list_array[6]}"
            header="${list_array[13]}"
            echo "$seq" >> $SSU_dir/seqs.txt
            echo "$header" >> $SSU_dir/headers.txt
        fi
        if [ "${list_array[7]}" != "NA" ]; then
            seq="${list_array[7]}"
            header="${list_array[13]}"
            echo "$seq" >> $ITS_dir/seqs.txt
            echo "$header" >> $ITS_dir/headers.txt
        fi
        if [ "${list_array[8]}" != "NA" ]; then
            seq="${list_array[8]}"
            header="${list_array[13]}"
            echo "$seq" >> $LSU_dir/seqs.txt
            echo "$header" >> $LSU_dir/headers.txt
        fi
    fi
done < <(tail -n +2 $list) # skip first line (column headers)

# Fetch sequences
for dir in $SSU_dir $ITS_dir $LSU_dir $longread_dir
do
    efetch -input "$dir/seqs.txt" -db nuccore -format fasta > "$dir/seqs.fasta"
    perl $project_dir/scripts/rename_multifasta_headers_according_to_list.pl "$dir/headers.txt" "$dir/seqs.fasta" > "$dir/seqs_formatted.fasta"
done

# Run ITSx on ITS sequences and long reads
mkdir $ITS_dir/itsx $longread_dir/itsx
ITSx -i $ITS_dir/seqs_formatted.fasta -o $ITS_dir/itsx/its --graphical F --save_regions all
ITSx -i $longread_dir/seqs_formatted.fasta -o $longread_dir/itsx/long-read --graphical F --save_regions all

# Check to see if output matches input
# Sorochytrium milnesiophthora and Sanchytriomycota taxa are missing from long-read 5.8S file
for dir in $ITS_dir $longread_dir
do
    cat $dir/seqs_formatted.fasta | grep ">" | wc -l
    for file in $dir/itsx/*.fasta
    do
        echo $file
        cat $file | grep ">" | wc -l
    done
done

# Reformat headers
sed -E 's/(^>.*)\|F\|.*/\1/' $ITS_dir/itsx/its.5_8S.fasta > $ITS_dir/itsx/its.5_8S_formatted.fasta
sed -E 's/(^>.*)\|F\|.*/\1/' $longread_dir/itsx/long-read.SSU.fasta > $longread_dir/itsx/long-read.SSU.formatted.fasta
sed -E 's/(^>.*)\|F\|.*/\1/' $longread_dir/itsx/long-read.5_8S.fasta > $longread_dir/itsx/long-read.5_8S.formatted.fasta
sed -E 's/(^>.*)\|F\|.*/\1/' $longread_dir/itsx/long-read.LSU.fasta > $longread_dir/itsx/long-read.LSU.formatted.fasta

# Combine ITSx output for long reads with single rDNA regions
for dir in $SSU_dir $ITS_dir $LSU_dir
do
    >$dir/seqs_formatted_combined.fasta
done
cat $SSU_dir/seqs_formatted.fasta >> $SSU_dir/seqs_formatted_combined.fasta
cat $longread_dir/itsx/long-read.SSU.formatted.fasta >> $SSU_dir/seqs_formatted_combined.fasta
cat $ITS_dir/itsx/its.5_8S_formatted.fasta >> $ITS_dir/seqs_formatted_combined.fasta
cat $longread_dir/itsx/long-read.5_8S.formatted.fasta >> $ITS_dir/seqs_formatted_combined.fasta 
cat $LSU_dir/seqs_formatted.fasta >> $LSU_dir/seqs_formatted_combined.fasta
cat $longread_dir/itsx/long-read.LSU.formatted.fasta >> $LSU_dir/seqs_formatted_combined.fasta

# Manual adjustments
# Add back in Sorochytrium and Sanchytriomycota 5.8S sequences
cat $project_dir/data/undetected_5_8S.fasta >> $ITS_dir/seqs_formatted_combined.fasta 
# Add Matteo Vecchi rock pool SSU ASVs
cat $project_dir/data/ASVs_SSU.fasta >> $SSU_dir/seqs_formatted_combined.fasta

# Align fasta files
fasta="seqs_formatted_combined.fasta"
for file in $SSU_dir/$fasta $ITS_dir/$fasta $LSU_dir/$fasta
do
    muscle -in $file -out ${file%.fasta}_aligned.fasta
done

# Trim alignment
fasta="seqs_formatted_combined_aligned.fasta"
for file in $SSU_dir/$fasta $ITS_dir/$fasta $LSU_dir/$fasta
do
    trimal -in $file -out ${file%.fasta}_trimmed.fasta -fasta -automated1 -htmlout ${file%seqs_formatted_combined_aligned.fasta}trimal.html
    sed 's/ [0-9]* bp//' ${file%.fasta}_trimmed.fasta > ${file%.fasta}_trimmed_formatted.fasta
done
