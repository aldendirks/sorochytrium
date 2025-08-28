#!/bin/bash

#SBATCH --job-name=phylogenomics
#SBATCH --account=tyjames1
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=10gb 
#SBATCH --time=2:00:00
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/home/adirks/logs/%x-%j.log
#SBATCH --error=/home/adirks/logs/%x-%j-E.log

# Load modules
module load Bioinformatics
module load python
module load muscle
module load fasttree

# Set variables
cores="$SLURM_NTASKS"
phylo_dir="$project_dir/phylo/phylogenomics"
busco_dir="$phylo_dir/busco"
runs_dir="$phylo_dir/runs"
results_dir="$phylo_dir/results"
list="$project_dir/data/genome_accessions.txt"
mkdir -p $runs_dir $results_dir

# Soft link BUSCO runs
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
    label="$name"_"$id"
    run="$busco_dir/$label/run_fungi_odb10"
    ln -s $run $runs_dir/run_$label
done < $list

# Run BUSCO phylogenomics pipeline
# Concatenation
#   --stop_early generates the files that are used by IQTREE2 without continuing on to phylogenetic inference
python BUSCO_phylogenomics.py -d $runs_dir -o $results_dir/supermatrix_30 --supermatrix --threads $cores --percent_single_copy 30 --stop_early
python BUSCO_phylogenomics.py -d $runs_dir -o $results_dir/supermatrix_75 --supermatrix --threads $cores --percent_single_copy 75 --stop_early
python BUSCO_phylogenomics.py -d $runs_dir -o $results_dir/supermatrix_90 --supermatrix --threads $cores --percent_single_copy 90 --stop_early
# Coalescence
python BUSCO_phylogenomics.py -d $runs_dir -o $results_dir/supertree --supertree --threads $cores
# This exact formatting (no spaces between -D and quotation marks, space after =, etc.) was found to be necessary for some reason
java -D"java.library.path= $ASTRAL_LIB_PATH" -jar $ASTRAL_JAR_PATH \
-i $results_dir/supertree/ALL.trees -o $results_dir/supertree/astral.tree

# IQTREE
for cutoff in 30 75 90
do
    # Set up directories
    iqtree_analysis_dir="$results_dir/supermatrix_$cutoff"
    aln_dir="$results_dir/supermatrix_$cutoff/trimmed_alignments"
    cd $iqtree_analysis_dir

    # Infer a concatenation-based species tree with 1000 ultrafast bootstrap and an edge-linked partition model
    iqtree3 -p $aln_dir -B 1000 -mset WAG,LG,JTT --prefix concat -T $cores

    # Test phylogenetic assumptions
    #   stationarity: amino acid frequencies remain constant over time
    #   homogeneity: substitution rates remain constant over time
    iqtree3 -p $aln_dir --symtest-remove-bad -B 1000 -mset WAG,LG,JTT --prefix good_partitions -T $cores
    iqtree3 -p $aln_dir --symtest-remove-good -B 1000 -mset WAG,LG,JTT --prefix bad_partitions -T $cores

    # Compare trees from all partitions and from "good" partitions only
    cat concat.treefile good_partitions.treefile > species_trees.treefile
    iqtree3 -p concat.best_scheme.nex -z species_trees.treefile -n 0 -zb 10000 -zw -au --prefix topotest_all -T $cores
    iqtree3 -p good_partitions.best_scheme.nex -z species_trees.treefile -n 0 -zb 10000 -zw -au --prefix topotest_pruned_good -T $cores
    iqtree3 -p bad_partitions.best_scheme.nex -z species_trees.treefile -n 0 -zb 10000 -zw -au --prefix topotest_pruned_bad -T $cores

    # Infer the locus trees
    iqtree3 -S concat.best_scheme.nex --prefix loci -T $cores

    # Compute gene concordance factors and site concordance factors
    iqtree3 -te concat.treefile -p concat.best_scheme.nex --gcf loci.treefile --df-tree --cf-verbose --scf 100 --prefix concord -T $cores

    # Infer the locus trees for good partitions
    iqtree3 -S good_partitions.best_scheme.nex --prefix loci_good_partitions -T $cores

    # Compute gene concordance factors and site concordance factors for good partitions tree
    iqtree3 -te good_partitions.treefile -p good_partitions.best_scheme.nex --gcf loci_good_partitions.treefile --df-tree --cf-verbose --scf 100 --prefix concord_good_partitions -T $cores
done