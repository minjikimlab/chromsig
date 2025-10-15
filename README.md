# Chrom-Sig: de-noising 1-dimensional genomic profiles by signal processing methods

## Features
* Chrom-Sig takes in Single-End or Paired-End read .bam or .bed files generated from ChIP-seq, ATAC-seq, snATAC-seq, or CUT&RUN experiments, among others 
* As an output, Chrom-Sig produces "pass" reads .bed file and .bedgraph files, along with peaks called by SICER on "pass" reads

## Installation
This program requires Python 3.9 to be able to run.
Start by installing Chrom-Sig using the following command:
```
git clone https://github.com/minjikimlab/chromsig.git
```
Navigating to the directory .../chromsig, do the following:
Install pyBedGraph and bedtools using the following commands:
### bedtools
```
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools-2.31.0.tar.gz
$ tar -zxvf bedtools-2.31.0.tar.gz
$ cd bedtools2
$ make
```
### pyBedGraph
```
git clone https://github.com/TheJacksonLaboratory/pyBedGraph.git
python -m venv pyBedGraph
source pyBedGraph/bin/activate
```
### Required Modules
Install all required packages listed in [requirements.txt](https://github.com/minjikimlab/chromsig/blob/main/v0.0.3/requirements.txt) using the following command:

`pip install -r requirements.txt`

If you are using conda, use the following command:

`conda install --file requirements.txt`

Install the chromsig package using the following command:

`git clone https://github.com/minjikimlab/chromsig.git`

## Running chromsig

Navigate into the chromsig directory, and create a folder called Data_Directory. In this folder, upload your input .bam or .bed file, and your reference genome sizes file (ex: [hg38.chrom.sizes](https://github.com/minjikimlab/chromsig/blob/main/v0.0.2/hg38.chrom.sizes)).
Return to the chromsig directory, and in you terminal enter the following command:

`sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed <BAM/BED> --dir <PATH_TO_DATA_DIR> --r <GENOME> --fdr <FDR_THRESH> --num <NUM_OF_PSEUDO_SAMPLES> --type <READ_TYPE> --lib <LIBRARY> --plot <TRUE>`

The following is a description of each of the arguments in the command:

ARGUMENTS:

    --bamorbed      Input bam or bed file containing your dataset.
    --dir           Directory containing the bam/bed file, as well as your reference genome sizes file (ex: /Users/chromsig/Data/).
    --r             Name of reference genome (ex: hg38).
    --fdr           FDR threshold for denoising (ex: 0.1).
    --num           Number of pseudo samples to be generated for denoising.
    --type          Type of read in your dataset, Single End or Paired End (ex: SE or PE).

OPTIONS:
      
    --lib_name      Name of dataset library (ex: For GM12878_ATAC-seq_ENCFF415FEC.bam, this would be ENCFF415FEC).
    --plot          Produce histogram plots of distances between read fragments (only for Paired End), type 'true' if you want plots, skip this option if you do not.
            
EXAMPLE:

    sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed GM12878_ATAC-seq_ENCFF415FEC.bam --dir /nfs/turbo/umms-minjilab/njgupta/chromsig/ATAC-seq_PE/ --r hg38 --fdr 0.1 --num 5000 --type PE

## Testing chromsig


## Citation
"Chrom-Sig: de-noising 1-dimensional genomic profiles by signal processing methods" by Nandita J. Gupta, Zachary Apell, and Minji Kim. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.12.670000v1) (2025), 670000.

For questions or bug reports, contact Nandita (njgupta@umich.edu) or visit the "Issues" page.
