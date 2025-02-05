# Rediscovery of _Sorochytrium milnesiophthora_

<!-- vscode-markdown-toc -->
* 1. [Setup](#Setup)
	* 1.1. [Software](#Software)
* 2. [Full-length rDNA sequencing](#Full-lengthrDNAsequencing)
* 3. [Genome assembly](#Genomeassembly)
* 4. [Multi-locus phylogenetics](#Multi-locusphylogenetics)
* 5. [Whole genome phylogenetics (phylogenomics)](#Wholegenomephylogeneticsphylogenomics)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

##  1. <a name='Setup'></a>Setup

Throughout this document, code in brackets should be replaced with your own commands. 

###  1.1. <a name='Software'></a>Software

The following software need to be installed and findable in $PATH.

**Published software used in this pipeline**
* 

**Other software used in this pipeline**
* 

###  1.2. <a name='Setupdirectorystructure'></a>Set up directory structure

The directory containing this README.md file will be the project directory. Change to the parent directory and run the code below to clone the GitHub directory `discinaceae_phylogenomics`. 

```
git clone https://github.com/aldendirks/sorochytrium.git
export PROJECT_DIR="$PWD/sorochytrium"
cd $PROJECT_DIR
mkdir -p data summary
```

Make sure the list of samples is unix compatible and ends with a space.


> **NOTE:** Modules are preinstalled software on the AGC computing cluster. Whenever possible, use a module because the software already works. Installing your own software locally is doable but is trickier the more dependencies involved. Search for the existence of a module with `module spider [name]`. You need to load the python module to use pip, a python package management command. The flag `--user` is necessary so that the software is installed in your home directory rather than the default higher level directory, which requires admin permissions.

> **NOTE:** The default version for a module may change over time, which can impact the pipeline. For example, if you install `duplex_tools` at one point in time with the most recent (default) python module set to version 3.10.4, later on the default version may change to a more recent one like 3.12.1. `duplex_tools` would need to be reinstalled with this newer version of python or the corresponding version of 'python' would need to be specified in the module command (e.g., `module load python/3.10.4`) 

R packages need to be installed to use MinIONQC. Launch R with `R` after loading the R module and in the R prompt install the following packages. More information on installing R packages on AGC can be found [here](https://arc.umich.edu/software/r/).


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

##  2. <a name='Full-lengthrDNAsequencing'></a>Full-length rDNA sequencing

This code is modified primarily from the [Nanopore amplicon pipeline](https://www.protocols.io/view/primary-data-analysis-basecalling-demultiplexing-a-dm6gpbm88lzp/v3) developed by Steve Russell with contributions from previous work done by Rabern Simmons, Katelyn McKindles, and Michelle Orozco-Quime. 

##  3. <a name='Genomeassembly'></a>Genome assembly

##  4. <a name='Multi-locusphylogenetics'></a>Multi-locus phylogenetics

##  5. <a name='Wholegenomephylogeneticsphylogenomics'></a>Whole genome phylogenetics (phylogenomics)