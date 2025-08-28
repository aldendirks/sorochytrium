#!/bin/bash

#SBATCH --job-name=infernal
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=10gb
#SBATCH --time=48:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics

# Set variables
cores="$SLURM_NTASKS"
dir="$project_dir/its"
dbs="$dir/dbs"
itsonedb="$dir/itsonedb"
rod="$dir/rod"
unite="$dir/unite"

# Get list of longest sequences for ROD and process with Infernal
cd $rod
mkdir -p $rod/infernal_2000/fasta $rod/infernal_2000/cmscan
for region in ITS1 ITS2 full
do
    # Start
    echo "Processing $region"

    # Sort and get longest sequences
    sort -nrk4 $rod/ROD_lengths_$region.txt | head -n2000 > $rod/infernal_2000/$region-top.txt

    # Pull out sequence and run Infernal cmscan
    while read -r accession sequence taxonomy length
    do
        id="$(echo $sequence | cut -d "/" -f 1)"
        segment="$(echo $sequence | cut -d "/" -f 2)"
        echo "$accession ${id}_${segment}"
        if [[ -s $rod/infernal_2000/fasta/$region-$accession-${id}_${segment}.fasta && -s $rod/infernal_2000/cmscan/$region-$accession-${id}_${segment}.out ]]; then
            echo "The FASTA and Infernal output files for '$accession ${id}_${segment}' exist and are not empty."
        else
            awk -v search="${sequence}" -f $scripts/fasta_search.awk $dbs/ROD_v1.2_operon_variants.fasta > $rod/infernal_2000/fasta/$region-$accession-${id}_${segment}.fasta
            cmscan --tblout $rod/infernal_2000/cmscan/$region-$accession-${id}_${segment}.out --fmt 3 --noali --cpu $cores $RFAM $rod/infernal_2000/fasta/$region-$accession-${id}_${segment}.fasta
        fi
    done < $rod/infernal_2000/$region-top.txt
    echo 
done

# Get sequence lengths from Infernal output
$project_dir/scripts/parse_cmscan.sh $rod/infernal_2000/cmscan/*.out > $rod/infernal_2000/infernal_lengths.tsv

# Get list of longest sequences for UNITE and process with infernal
cd $unite
mkdir -p $unite/infernal_2000/fasta $unite/infernal_2000/cmscan
for region in ITS1 ITS2 full
do
    # Start
    echo "Processing $region"

    # Sort and get longest sequences
    sort -nrk3 $unite/UNITE_lengths_$region.txt | head -n2000 > $unite/infernal_2000/$region-top.txt

    # Pull out sequence and run Infernal cmscan
    while read -r accession taxonomy length
    do
        echo $accession
        if [[ -s $unite/infernal_2000/fasta/$accession.fasta && -s $unite/infernal_2000/cmscan/$accession.out ]]; then
            echo "The FASTA and Infernal output files for '$accession' exist and are not empty."
        else
            awk -v search="${accession}\\\\|" -f $scripts/fasta_search.awk $dbs/UNITE_public_all_2025-02-19.fasta | sed 's/|.*//' > $unite/infernal_2000/fasta/$accession.fasta
            cmscan --tblout $unite/infernal_2000/cmscan/$accession.out --fmt 3 --noali --cpu $cores $RFAM $unite/infernal_2000/fasta/$accession.fasta
        fi
    done < $unite/infernal_2000/$region-top.txt
    echo 
done

# Get sequence lengths from Infernal output
$scripts/parse_cmscan.sh $unite/infernal_2000/cmscan/*.out > $unite/infernal_2000/infernal_lengths.tsv

# Get list of shortest sequences for UNITE and process with infernal
cd $unite
mkdir -p $unite/infernal_short/fasta $unite/infernal_short/cmscan
for region in full
do
    # Start
    echo "Processing $region"

    # Sort and get longest sequences
    sort -nrk3 $unite/UNITE_lengths_$region.txt | tail -n2000 > $unite/infernal_short/$region-bottom.txt

    # Pull out sequence and run Infernal cmscan
    while read -r accession taxonomy length
    do
        echo $accession
        if [[ -s $unite/infernal_short/fasta/$accession.fasta && -s $unite/infernal_short/cmscan/$accession.out ]]; then
            echo "The FASTA and Infernal output files for '$accession' exist and are not empty."
        else
            awk -v search="${accession}\\\\|" -f $scripts/fasta_search.awk $dbs/UNITE_public_all_2025-02-19.fasta | sed 's/|.*//' > $unite/infernal_short/fasta/$accession.fasta
            cmscan --tblout $unite/infernal_short/cmscan/$accession.out --fmt 3 --noali --cpu $cores $RFAM $unite/infernal_short/fasta/$accession.fasta
        fi
    done < $unite/infernal_short/$region-bottom.txt
    echo 
done

# Get sequence lengths from Infernal output
$scripts/parse_cmscan.sh $unite/infernal_short/cmscan/*.out > $unite/infernal_short/infernal_lengths.tsv























# # Combine data for ITS1
# printf "number\taccession\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tits1\tits2\tfull\tfive\tits1_infernal\n" > "$dir"/data/its1_infernal.tsv
# while read -r number accession kingdom phylum class order family genus species its1 its2 full five
# do
#     file_formatted="$(grep -v "#" "$dir"/its1/$accession.out | sed 's/  */ /g')"
#     ssu_end="$(grep "SSU_rRNA_eukarya" <(echo "$file_formatted") | cut -f 9 -d " ")"
#     five_start="$(grep "5_8S_rRNA" <(echo "$file_formatted") | cut -f 8 -d " " | head -n 1)"
#     its1_infernal=$((five_start-ssu_end))
#     printf "$number\t$accession\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\t$its1\t$its2\t$full\t$five\t$its1_infernal\n" >> "$dir"/data/its1_infernal.tsv
# done < <(tail -n +2 "$dir/data/topits1.tsv")
# cat "$dir"/data/its1_infernal.tsv | cut -f 14 | head

# # Combine data for ITS2
# printf "number\taccession\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tits1\tits2\tfull\tfive\tits2_infernal\n" > "$dir"/data/its2_infernal.tsv
# while read -r number accession kingdom phylum class order family genus species its1 its2 full five
# do
#     file_formatted="$(grep -v "#" "$dir"/its2/$accession.out | sed 's/  */ /g')"
#     lsu_start="$(grep "LSU_rRNA_eukarya" <(echo "$file_formatted") | cut -f 8 -d " " | head -n 1)"
#     five_end="$(grep "5_8S_rRNA" <(echo "$file_formatted") | cut -f 9 -d " " | head -n 1)"
#     its2_infernal=$((lsu_start-five_end))
#     printf "$number\t$accession\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\t$its1\t$its2\t$full\t$five\t$its2_infernal\n" >> "$dir"/data/its2_infernal.tsv
# done < <(tail -n +2 "$dir/data/topits2.tsv")
# cat "$dir"/data/its2_infernal.tsv | cut -f 14 | head

# # Combine data for ITS2 (501-2000)
# printf "number\taccession\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tits1\tits2\tfull\tfive\tits2_infernal\n" > "$dir"/data/its2_2_infernal.tsv
# while read -r number accession kingdom phylum class order family genus species its1 its2 full five
# do
#     file_formatted="$(grep -v "#" "$dir"/its2_2/$accession.out | sed 's/  */ /g')"
#     lsu_start="$(grep "LSU_rRNA_eukarya" <(echo "$file_formatted") | cut -f 8 -d " " | head -n 1)"
#     five_end="$(grep "5_8S_rRNA" <(echo "$file_formatted") | cut -f 9 -d " " | head -n 1)"
#     its2_infernal=$((lsu_start-five_end))
#     printf "$number\t$accession\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\t$its1\t$its2\t$full\t$five\t$its2_infernal\n" >> "$dir"/data/its2_2_infernal.tsv
# done < <(tail -n +2 "$dir/data/topits2_2.tsv")
# cat "$dir"/data/its2_2_infernal.tsv | cut -f 14 | head