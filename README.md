# Rediscovery of _Sorochytrium milnesiophthora_

<!-- vscode-markdown-toc -->
* 1. [Setup](#Setup)
	* 1.1. [Software](#Software)
	* 1.2. [Set up directory structure](#Setupdirectorystructure)
* 2. [Full-length rDNA sequencing](#Full-lengthrDNAsequencing)
	* 2.1. [Basecalling](#Basecalling)
	* 2.2. [Quality control](#Qualitycontrol)
	* 2.3. [Merge the simplex and duplex reads](#Mergethesimplexandduplexreads)
	* 2.4. [Demultiplexing](#Demultiplexing)
	* 2.5. [Create consensus sequences](#Createconsensussequences)
	* 2.6. [Aggregate consensus sequences](#Aggregateconsensussequences)
* 3. [Genome assembly](#Genomeassembly)
	* 3.1. [Acquire raw data](#Acquirerawdata)
		* 3.1.1. [Download raw sequencing reads for whole genome assembly](#Downloadrawsequencingreadsforwholegenomeassembly)
	* 3.2. [Assemble genome](#Assemblegenome)
	* 3.3. [Evaluate genome](#Evaluategenome)
	* 3.4. [Annotate genome](#Annotategenome)
* 4. [Multi-locus phylogenetics](#Multi-locusphylogenetics)
* 5. [ Phylogenomics](#Phylogenomics)
	* 5.1. [Phylogenomics analysis](#Phylogenomicsanalysis)
	* 5.2. [Visualize phylogenomics trees](#Visualizephylogenomicstrees)
* 6. [ITS analyses](#ITSanalyses)
	* 6.1. [Compare ITS lengths](#CompareITSlengths)
		* 6.1.1. [Download datasets and run ITSx](#DownloaddatasetsandrunITSx)
		* 6.1.2. [Run Infernal](#RunInfernal)
		* 6.1.3. [Visualize data](#Visualizedata)
	* 6.2. [ITS intragenomic variation](#ITSintragenomicvariation)
		* 6.2.1. [Calculate coverage to estimate rDNA copy number](#CalculatecoveragetoestimaterDNAcopynumber)
		* 6.2.2. [Extract rDNA sequences from the genome](#ExtractrDNAsequencesfromthegenome)
	* 6.3. [Global distribution of Sorochytrium](#GlobaldistributionofSorochytrium)
		* 6.3.1. [Extract _Sorochytrium_ from lichen metagenome](#Extract_Sorochytrium_fromlichenmetagenome)
		* 6.3.2. [Extract Sorochytrium from *Barbacenia* metagenome](#ExtractSorochytriumfromBarbaceniametagenome)
* 7. [Supplemental](#Supplemental)
	* 7.1. [Installing software](#Installingsoftware)
	* 7.2. [Installing NR database for AvP](#InstallingNRdatabaseforAvP)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Setup'></a>Setup

Throughout this document, code in brackets should be replaced with your own commands. 

###  1.1. <a name='Software'></a>Software

The following software need to be installed and findable in `$PATH`.

**Published software used in this pipeline**
* `antismash` 7.1.0
* `ASTRAL-III` 5.15.5 
* `AvP` 1.0.2
* `BBTools` 38.96
* `Bedtools` 2.30.0
* `BLAST+` 2.12.0
* `BUSCO` 5.3.0
* [`BUSCO_phylogenomics.py`](https://github.com/jamiemcg/BUSCO_phylogenomics)
* `DIAMOND` 2.0.12
* `Duplex Tools` 0.3.1
* `EDTA` 2.2.2
* `Entrez Direct`
* [`extractITS.py`](https://github.com/fantin-mesny/Extract-ITS-sequences-from-a-fungal-genome)
* `FastTree` 2.1.10
* `Flye` 2.9.2
* `Funannotate` 1.8.14 and dependencies
* `GenomeScope`
* `Guppy` 6.4.6
* `Hapo-G` 1.3.6 
* `ITSx` 1.1.3
* `IQ-TREE` 3.0.1
* `Jellyfish` 2.3.0
* `MiniBar` 0.21 
* `MinIONQC` 1.4.1
* `MUSCLE` 3.8.31
* `NCBI Datasets` 16.40.1
* `NextPolish` 1.4.1
* `NGSpeciesID` 0.1.3 
* `pycoQC` 2.5.2 
* `Python` 3.10.4
* `QUAST` 5.2.0
* `R` 4.2.0
* `R` packages `data.table`, `futile.logger`, `ggplot2`, `ggtree`, `optparse`, `plyr`, `readr`, `reshape2`, `scales`, `tidyverse`, `viridis`, `yaml`
* `Singularity` 4.1.3
* `Seqtk` 1.3
* `SPAdes` 3.15.5
* `sratoolkit` 3.1.1
* `trimAl` 1.2rev59

###  1.2. <a name='Setupdirectorystructure'></a>Set up directory structure

The directory containing this `README.md` file will be the project directory. Change to the parent directory and run the code below to clone the GitHub directory `sorochytrium`. 

```
git clone https://github.com/aldendirks/sorochytrium.git
export project_dir="$PWD/sorochytrium"
cd $project_dir
mkdir -p data figures its meta phylo results scripts seqs
```

##  2. <a name='Full-lengthrDNAsequencing'></a>Full-length rDNA sequencing

This code is modified primarily from the [Nanopore amplicon pipeline](https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3) developed by Steve Russell with contributions from previous work done by Rabern Simmons, Katelyn McKindles, and Michelle Orozco-Quime. 

###  2.1. <a name='Basecalling'></a>Basecalling

Basecalling is the process of converting the raw data in FAST5 format (nanopore electrical signals) to FASTQ files containing nucleotide sequences and quality information. *Sorochytrium* was only one sample in a multiplexed run. The raw FAST5 files are not public but can be made available upon request. 

> **NOTE:** In FASTQ files, [nucleotide quality scores](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScores_swBS.htm) are represented by [unique symbols](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm#).

A number of basecalling software exist, but the current "official" one is Guppy. Guppy basecalling is the most computationally intensive part of this pipeline. Running Guppy with CPUs is possible but takes way more time and resources than running it with GPUs. 

```
sbatch $project_dir/scripts/guppy.sh
```

This script first converts the FAST5 files to FASTQ format. Then, [Duplex Tools](https://github.com/nanoporetech/duplex-tools) is used to identify duplex reads. Around 20-40% of the reads are duplex reads, composed of multiple different sequences fused together. These are filtered and rerun with `guppy_basecaller_duplex`. 

###  2.2. <a name='Qualitycontrol'></a>Quality control

MinIONQC produces high-quality graphs from the Guppy summary files. 

```
Rscript $project_dir/scripts/MinIONQC.R -i $project_dir/seqs/rdna/guppy/simplex_calls/sequencing_summary.txt -o $project_dir/seqs/rdna/QC_reports
Rscript $project_dir/scripts/MinIONQC.R -i $project_dir/seqs/rdna/guppy/duplex_calls/sequencing_summary.txt -o $project_dir/seqs/rdna/QC_reports
```

###  2.3. <a name='Mergethesimplexandduplexreads'></a>Merge the simplex and duplex reads

```
mkdir $project_dir/seqs/rdna/guppy/merged_calls
{ sed 's/ /\n/' $project_dir/seqs/rdna/guppy/simplex_calls/pairs/pair_ids_filtered.txt" | seqkit grep -v -f - "$project_dir/seqs/rdna/guppy/simplex_calls/pass/basecall.fastq" ; cat $project_dir/seqs/rdna/guppy/duplex_calls/pass/*.fastq ; } > $project_dir/seqs/rdna/guppy/merged_calls/merged.fastq
cat $project_dir/seqs/rdna/guppy/merged_calls/merged.fastq | wc -l | awk '{print $1/4}' > $project_dir/seqs/rdna/summary/merged_reads.txt
```

###  2.4. <a name='Demultiplexing'></a>Demultiplexing

Use `minibar.py` to pull individual sample files from the merged data file with reference to a metadata file. 

```
mkdir $project_dir/seqs/rdna/NGSpeciesID
cd $project_dir/seqs/rdna/NGSpeciesID
$project_dir/scripts/minibar.py -F $project_dir/data/metadata.txt $project_dir/seqs/rdna/guppy/merged_calls/merged.fastq
```

###  2.5. <a name='Createconsensussequences'></a>Create consensus sequences

```
conda activate NGSpeciesID
cd $project_dir/seqs/rdna/NGSpeciesID
for file in *.fastq
do
    bn=`basename $file .fastq`
    NGSpeciesID --ont --consensus --sample_size 1000 --m 5000 --s 2000 --medaka --primer_file $project_dir/data/primers.txt --fastq $file --outfolder ${bn}
done

# Modify the command for files that fall outside the default range, e.g. Sorochytrium (large rDNA)
for file in sample_s1865.1*.fastq
do
    bn=`basename $file .fastq`
    NGSpeciesID --ont --consensus --sample_size 1000 --m 10000 --s 2000 --medaka --primer_file $project_dir/data/primers.txt --fastq $file --outfolder ${bn}
done
conda deactivate
```

###  2.6. <a name='Aggregateconsensussequences'></a>Aggregate consensus sequences

```
> $project_dir/seqs/rdna/summary/consensus_seqs.fasta
cd $project_dir/seqs/rdna/NGSpeciesID
for file in *.fastq
do
    bn=`basename $file .fastq`
    for file in $(ls $bn/consensus_reference_*.fasta)
    do
        sed 's%^>\(.*\)%>'$bn'_\1%I' $file >> $project_dir/seqs/rdna/summary/consensus_seqs.fasta
    done
done
cd $project_dir/seqs/rdna/summary
```

##  3. <a name='Genomeassembly'></a>Genome assembly

###  3.1. <a name='Acquirerawdata'></a>Acquire raw data

####  3.1.1. <a name='Downloadrawsequencingreadsforwholegenomeassembly'></a>Download raw sequencing reads for whole genome assembly

The Nanopore and Illumina raw reads are available in NCBI SRA under BioProject [PRJNA1015119](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1015119).

* Nanopore: run [SRR26060970](https://www.ncbi.nlm.nih.gov/sra/SRX21777098[accn])
* Illumina: run [SRR26060971](https://www.ncbi.nlm.nih.gov/sra/SRX21777097[accn])

```
mkdir -p $project_dir/seqs/genome/raw
cd $project_dir/seqs/genome/raw
for run in SRR26060970 SRR26060971
do
    prefetch $run
    fasterq-dump --gzip $run
done
```

###  3.2. <a name='Assemblegenome'></a>Assemble genome

Assemble the Nanopore reads. 

```
sbatch $project_dir/scripts/flye.sh
```

Filter and trim the Illumina reads. 

```
sbatch $project_dir/scripts/check_read_quality_and_trim.sh
```

Polish and phase the Nanopore assembly with the filtered and trimmed Illumina reads.

```
sbatch $project_dir/scripts/polish_and_phase.sh
```

###  3.3. <a name='Evaluategenome'></a>Evaluate genome

Evaluate genome assembly with quast. 

```
sbatch $project_dir/scripts/quast.sh
```

Run GenomeScope to estimate genome size and heterozygosity. 

```
sbatch $project_dir/scripts/genomescope.sh
```

Create a kmer histogram for visualization in R.

```
sbatch $project_dir/scripts/calculate_kmer_count.sh
cd "$project_dir"/seqs/genome/kmers
Rscript $project_dir/scripts/kmers_histogram.R
```

###  3.4. <a name='Annotategenome'></a>Annotate genome

Annotate the genome with Funannotate. 

```
sbatch $project_dir/scripts/funannotate_predict.sh
sbatch $project_dir/scripts/antismash7.sh
sbatch $project_dir/scripts/funannotate_annotate.sh
```

Identify transposable elements with EDTA. 

```
sbatch $project_dir/scripts/edta.sh
```

##  4. <a name='Multi-locusphylogenetics'></a>Multi-locus phylogenetics

Prepare files for phylogenetic analysis.

```
bash $project_dir/scripts/phylo_prep.sh
```

Run IQ-TREE. 

```
sbatch $project_dir/scripts/phylo_iqtree.sh
```

##  5. <a name='Phylogenomics'></a> Phylogenomics

###  5.1. <a name='Phylogenomicsanalysis'></a>Phylogenomics analysis

Download genomes from GenBank.

```
bash phylogenomics_prep.sh
```

Run BUSCO on genomes to identify conserved single-copy orthologs.

```
mkdir -p $project_dir/phylo/phylogenomics/busco
cd $project_dir/phylo/phylogenomics/busco
busco --download fungi_odb10
for fasta in $project_dir/seqs/phylogenomics/*.fasta
do
	sbatch $project_dir/scripts/phylogenomics_busco.sh $fasta
done
```

The [BUSCO phylogenomics pipeline](https://github.com/jamiemcg/BUSCO_phylogenomics) is used for coalescence and to prepare the dataset for IQTREE concatenation analysis. You will need to define the path to the ASTRAL library files and ASTRAL jar file, e.g., "/home/adirks/apps/ASTRAL/Astral/lib/" and "/home/adirks/apps/ASTRAL/Astral/astral.5.15.5.jar".

```
# Consider running the following code to find the paths if available, although find can take a while
#find / -name "astral.5.15.5.jar"
export ASTRAL_LIB_PATH="/path/to/Astral/lib/"
export ASTRAL_JAR_PATH="/path/to/astral.5.15.5.jar"
```

```
sbatch $project_dir/scripts/phylogenomics_trees.sh
```

###  5.2. <a name='Visualizephylogenomicstrees'></a>Visualize phylogenomics trees

```
Rscript $project_dir/scripts/phylogenomics_visualize.R
```

##  6. <a name='ITSanalyses'></a>ITS analyses

###  6.1. <a name='CompareITSlengths'></a>Compare ITS lengths

####  6.1.1. <a name='DownloaddatasetsandrunITSx'></a>Download datasets and run ITSx

Download the UNITE and ROD ITS datasets. ITSoneDB needs to be downloaded locally. ITSx identifies the SSU, 5.8S, and LSU bounds and reports the lengths of these different regions for the ITSoneDB, ROD, and UNITE datasets. The script below runs ITSx and parses the output files to extract length data. 

```
sbatch $project_dir/scripts/itsx.sh
```

####  6.1.2. <a name='RunInfernal'></a>Run Infernal

Process the longest sequences with Infernal to check the ITSx results. Make sure to designate the path to your locally downloaded Rfam database. 

```
export RFAM="/path/to/Rfam.cm"
sbatch $project_dir/scripts/infernal.sh
```

####  6.1.3. <a name='Visualizedata'></a>Visualize data

Open the document `$project_dir/scripts/its_lengths.R` in RSudio and work through the code. 

###  6.2. <a name='ITSintragenomicvariation'></a>ITS intragenomic variation

####  6.2.1. <a name='CalculatecoveragetoestimaterDNAcopynumber'></a>Calculate coverage to estimate rDNA copy number

BBMap is used to map the trimmed Illumina reads to the assembled genome and determine coverage. Coverage is calculated for the genome, 10 single-copy orthologs (to identify the coverage signature of a single-copy gene), and the rDNA sequences. 

```
sbatch $project_dir/scripts/calculate_coverage.sh
```

####  6.2.2. <a name='ExtractrDNAsequencesfromthegenome'></a>Extract rDNA sequences from the genome

```
sbatch $project_dir/scripts/extract_rdna.sh
```

###  6.3. <a name='GlobaldistributionofSorochytrium'></a>Global distribution of Sorochytrium

####  6.3.1. <a name='Extract_Sorochytrium_fromlichenmetagenome'></a>Extract _Sorochytrium_ from lichen metagenome

Download *Parmelia saxatilis* lichen Illumina reads (metagenome reads of sample from which Sorochytrium milnesiophthora was isolated).

```
lichen_dir="$project_dir/seqs/parmelia"
mkdir -p $lichen_dir/raw
cd $lichen_dir/raw
prefetch SRR34964407
fasterq-dump --gzip SRR34964407
```

Trim Illumina reads, BLAST against _Sorochytrium_ genome, subset matching reads, and assemble into a Sorochytrium metagenome-derived assembly. 

```
sbatch $project_dir/scripts/metagenomics_parmelia.sh
```

####  6.3.2. <a name='ExtractSorochytriumfromBarbaceniametagenome'></a>Extract Sorochytrium from *Barbacenia* metagenome

Pebblescout showed *Sorochytrium milnesiophthora* 5.8S to have a single hit to a rhizosphere metagenome belonging to *Barbacenia macrantha* from Brazil. Download the raw reads. 

```
barb_dir="$project_dir/seqs/parmelia"
mkdir -p $barb_dir/raw
cd $barb_dir/raw
prefetch SRR34964407
fasterq-dump --gzip SRR34964407
```

Extract Sorochytrium from the rhizosphere metagenome. 

```
sbatch $project_dir/scripts/metagenomics_barbacenia.sh
```

##  7. <a name='Supplemental'></a>Supplemental

###  7.1. <a name='Installingsoftware'></a>Installing software

This code is formatted for University of Michigan Advanced Genomics Core (AGC).

Install code for the Nanopore amplicon pipeline. 

```
module load Bioinformatics
module load python
pip install --user duplex_tools
```

> **NOTE:** Modules are preinstalled software on the AGC computing cluster. Whenever possible, use a module because the software already works. Installing your own software locally is doable but is trickier the more dependencies involved. Search for the existence of a module with `module spider [name]`. You need to load the python module to use pip, a python package management command. The flag `--user` is necessary so that the software is installed in your home directory rather than the default higher level directory, which requires admin permissions.

> **NOTE:** The default version for a module may change over time, which can impact the pipeline. For example, if you install `duplex_tools` at one point in time with the most recent (default) python module set to version 3.10.4, later on the default version may change to a more recent one like 3.12.1. `duplex_tools` would need to be reinstalled with this newer version of python or the corresponding version of 'python' would need to be specified in the module command (e.g., `module load python/3.10.4`) 

R packages need to be installed to use MinIONQC in the Nanopore amplicon pipeline. Launch R with `R` after loading the R module and in the R prompt install the following packages. More information on installing R packages on AGC can be found [here](https://arc.umich.edu/software/r/).

```
module load R
R
install.packages(c("data.table","futile.logger","ggplot2","optparse","plyr","readr","reshape2","scales","viridis","yaml"))
quit()
```

Install NGSpeciesID with Conda (this takes a little while).

```
conda create -n NGSpeciesID python=3.6 pip
conda activate NGSpeciesID
conda install --yes -c conda-forge -c bioconda medaka==0.11.5 openblas==0.3.3 spoa racon minimap2
pip install --user NGSpeciesID
#NGSpeciesID --help
```