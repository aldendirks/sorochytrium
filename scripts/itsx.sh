#!/bin/bash

#SBATCH --job-name=itsx
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=40gb
#SBATCH --time=1-12:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Set variables
output="$project_dir/its"
dbs="$output/dbs"
itsonedb="$output/itsonedb"
rod="$output/rod"
unite="$output/unite"
mkdir -p $dbs $itsonedb $rod $unite

# Get datasets
cd $dbs
# ITSoneDB (transfer locally from https://itsonedb.cloud.ba.infn.it/index.jsp#)
# ROD (unzipped file is 386 MB)
git clone https://github.com/krabberod/ROD
gunzip $dbs/ROD/ROD_v1.2_operon_variants.fasta.gz
mv $dbs/ROD/ROD_v1.2_operon_variants.fasta $dbs
rm -r $dbs/ROD
# UNITE (unzipped file is 2.1 GB)
curl https://s3.hpc.ut.ee/plutof-public/original/909ba708-b58a-4b81-a17d-7c79f6e4b2a4.gz \
--output $dbs/UNITE_public_all_2025-02-19.fasta.gz
gunzip $dbs/UNITE_public_all_2025-02-19.fasta.gz

# Run ITSx
# ITSx fails on itsonedb because there are no conserved regions for border identification, just ITS1
cd $rod
rodfa="$dbs/ROD_v1.2_operon_variants.fasta"
ITSx -i $rodfa -o rod --cpu 24 --multi_thread T --graphical F --save_regions all --temp $tmp
cd $unite
unitefa="$dbs/UNITE_public_all_2025-02-19.fasta"
ITSx -i $unitefa -o unite --cpu 24 --multi_thread T --graphical F --save_regions all --temp $tmp

# Parse itsonedb data
cd $itsonedb
for file in $dbs/ITSoneDB_ITS1_ENA.fa $dbs/ITSoneDB_ITS1_HMM.fa $dbs/ITSoneDB_ITS1_ENAandHMM.fa $dbs/ITSoneDB_rep_seq_and_flanking_r144.fasta
do
    # Start
    echo "Processing $(basename $file)"

    # Set type variable based on file name
    if [ "$file" == "$dbs/ITSoneDB_ITS1_ENA.fa" ]; then
        type="ENA"
    elif [ "$file" == "$dbs/ITSoneDB_ITS1_HMM.fa" ]; then
        type="HMM"
    elif [ "$file" == "$dbs/ITSoneDB_ITS1_ENAandHMM.fa" ]; then
        type="both"
    else
        type="rep"
    fi
    
    # Get information 
    awk -F  '|' -v var="$type" 'BEGIN{printf "File\tAccession\tSpecies\tLength\n"} /^>/ {print var "\t" $1 "\t" $3 "\t" $4}' $file |\
    sed 's/>//g;s/ITS1 located by .* annotation, //g;s/\([0-9]*\)bp/\1/g' > $itsonedb/ITSoneDB_lengths_$type.txt

    # Sanity check: check to see if any lines are missing sequence length
    while read file accession species length 
    do
        if [ -z "$length" ]; then
            echo "$accession is missing sequence length"
        fi
    done < $itsonedb/ITSoneDB_lengths_$type.txt

    # Look at longest sequences
    # $'\t' is "ANSI-C Quoting", the dollar sign and single quotes tell bash to understand the special symbol
    sort -t $'\t' -nrk4 $itsonedb/ITSoneDB_lengths_$type.txt | head
    echo
done

# Parse ROD data
cd $rod
for region in ITS1 5_8S ITS2 full 
do
    # Start
    echo "Processing $region"

    # Using | or ( as delimiters, print first, second, third, and last columns of FASTA headers
    awk -F '[|(]' 'BEGIN{printf "Accession\tSequence\tTaxonomy\tLength\n"} /^>/ {print $1 "\t" $2 "\t" $3 "\t" $NF}' $rod/rod.$region.fasta |\
    sed 's/>//g;s/ bp)//g;s/ on main strand//g' > $rod/ROD_lengths_$region.txt

    # Sanity check: check to see if any lines are missing sequence length
    while read accession taxonomy length 
    do
        if [ -z "$length" ]; then
            echo "$accession is missing sequence length"
        fi
    done < $rod/ROD_lengths_$region.txt

    # Look at longest sequences
    cat $rod/ROD_lengths_$region.txt | cut -f 4 | sort -nr | head
    echo
done

# Parse UNITE data
cd $unite
for region in ITS1 5_8S ITS2 full 
do
    # Start
    echo "Processing $region"

    # Using | or ( as delimiters, print first, second, and last columns of FASTA headers (information in fifth or sixth column, but always last)
    awk -F '[|(]' 'BEGIN{printf "Accession\tTaxonomy\tLength\n"} /^>/ {print $1 "\t" $2 "\t" $NF}' $unite/unite.$region.fasta |\
    sed 's/>//g;s/ bp)//g;s/ on main strand//g' > $unite/UNITE_lengths_$region.txt

    # Sanity check: check to see if any lines are missing sequence length
    while read accession taxonomy length 
    do
        if [ -z "$length" ]; then
            echo "$accession is missing sequence length"
        fi
    done < $unite/UNITE_lengths_$region.txt

    # Look at longest sequences
    cat $unite/UNITE_lengths_$region.txt | cut -f 3 | sort -nr | head
    echo
done
