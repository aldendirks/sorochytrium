#!/bin/bash

# Usage: ./batch_parse_cmscan.sh /path/to/*.out > combined_output.tsv

# Print header
echo -e "File\tITS1\t5.8S\tITS2\tFull_ITS\tNotes"

# Loop through each .out file given as argument
for infile in "$@"; do

    # Get clean filename (no path, no suffix)
    filename=$(basename "$infile" .out)

    # Extract relevant hits
    hits=$(awk '!/^#/ && NF > 0' "$infile")

    # Get sequence length
    fasta_length=$(echo "$hits" | awk '{ print $19 }' | uniq)

    # Function to get max interval from Infernal cmscan output
    get_max_interval() {
        model="$1"
        echo "$hits" | awk -v model="$model" '
        $1 == model {
            start = $8; end = $9;
            if (start > end) { temp = start; start = end; end = temp; }
            if (min_start == "" || start < min_start) {
                min_start = start;
            }
            if (max_end == "" || end > max_end) {
                max_end = end;
            }
            count++;
        }
        END {
            if (count == 0) {
                print "NA", "NA", "NA";
            } else {
                note = "ok";
                if (count > 1) {
                    note = "duplicate";
                }
                print min_start, max_end, note;
            }
        }'
    }

    # Function specifically for 5.8S when duplicated
    # Get the values that are between SSU and LSU
    get_5_8S_interval() {
        echo "$hits" | awk '
        $1 == "LSU_rRNA_eukarya" {
            startLSU = $8; endLSU = $9;
        }
        $1 == "SSU_rRNA_eukarya" {
            startSSU = $8; endSSU = $9;
        }
        $1 == "5_8S_rRNA" {
            start = $8; end = $9;
            if (start > endSSU && end < startLSU) {
                start58S = start;
                end58S = end;
            }
        }
        END {
            print start58S, end58S;
        }'
    }

    # Get positions
    read SSU_start SSU_end SSU_note <<< $(get_max_interval "SSU_rRNA_eukarya")
    read S58_start S58_end S58_note <<< $(get_max_interval "5_8S_rRNA")
    read LSU_start LSU_end LSU_note <<< $(get_max_interval "LSU_rRNA_eukarya")

    # Compile notes
    note_list=()
    [[ $SSU_note == "duplicate" ]] && note_list+=("SSU_dup")
    [[ $S58_note == "duplicate" ]] && note_list+=("5.8S_dup")
    [[ $LSU_note == "duplicate" ]] && note_list+=("LSU_dup")

    # Initialize
    ITS1="NA"; S58="NA"; ITS2="NA"; FULL="NA"; notes="NA"

    # Get positions
    if [[ $S58_note == "duplicate" ]]; then
        read S58_start S58_end <<< $(get_5_8S_interval)
    fi

    # Calculate
    if [[ $SSU_end != "NA" && $S58_start != "NA" ]]; then
        ITS1=$(( S58_start - SSU_end ))
    elif [[ $SSU_end == "NA" && $S58_start != "NA" ]]; then
        ITS1=$(( S58_start - 1 )) && note_list+=("SSU-missing")
    fi

    if [[ $S58_start != "NA" && $S58_end != "NA" ]]; then
        S58=$(( S58_end - S58_start + 1 ))
    fi

    if [[ $S58_end != "NA" && $LSU_start != "NA" ]]; then
        ITS2=$(( LSU_start - S58_end ))
    elif [[ $S58_end != "NA" && $LSU_start == "NA" ]]; then
        ITS2=$(( fasta_length - S58_end + 1 )) && note_list+=("LSU-missing")
    fi

    if [[ $SSU_end != "NA" && $LSU_start != "NA" ]]; then
        FULL=$(( LSU_start - SSU_end ))
    fi

    # Export notes
    if [[ ${#note_list[@]} -gt 0 ]]; then
        notes=$(IFS=','; echo "${note_list[*]}")
    fi

    # Print result
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$filename" "$ITS1" "$S58" "$ITS2" "$FULL" "$notes"
done
