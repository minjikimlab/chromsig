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
Navigate to the directory .../chromsig/v0.0.3.

Install all required packages listed in [requirements.txt](https://github.com/minjikimlab/chromsig/blob/main/v0.0.3/requirements.txt) using the following command:

`pip install -r requirements.txt`

If you are using conda, use the following command:

`conda install --file requirements.txt`

In each of the bash files (files ending in .sh), modify the parameters
```
#SBATCH --account=minjilab0 
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/overall_program_output/slurm-%j.out
#SBATCH --error=/nfs/turbo/umms-minjilab/njgupta/chromsig/overall_program_error/slurm-%j.out
#SBATCH --mail-user=njgupta@umich.edu
```
to match the account using which you are submitting the job, the output directory for the job report, the output directory for any error reports, and the email address you want the job report to be mailed to.

## Running chromsig
In the same directory (.../chromsig/v0.0.3), create a folder called Data_Directory. In this folder, upload your input .bam or .bed file, and your reference genome sizes file (ex: [hg38.chrom.sizes](https://github.com/minjikimlab/chromsig/blob/main/v0.0.3/hg38.chrom.sizes)).
Return to .../chromsig/v0.0.3, and in your terminal enter the following command:

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
      
    --lib           Name of dataset library (ex: For GM12878_ATAC-seq_ENCFF415FEC.bam, this would be ENCFF415FEC).
    --plot          Produce histogram plots of distances between read fragments (only for Paired End), type 'true' if you want plots, skip this option if you do not.
            
EXAMPLE:

    sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed GM12878_ATAC-seq_ENCFF415FEC.bam --dir /nfs/turbo/umms-minjilab/njgupta/chromsig/ATAC-seq_PE/ --r hg38 --fdr 0.1 --num 5000 --type PE

## Testing chromsig
Two test files have been included in the .../chromsig directory -- the bed file contains test data from chr21 of the GM12878_ATAC-seq_ENCFF415FEC.bam file, and the hg38_test.chrom.sizes file contains the size of chr21 needed to run this data. Move both of these files to your Data Directory folder .../chromsig/v0.0.3/Data_Directory. Navigate back to .../chromsig/v0.0.3, and run the following command:
```
sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed chromsig_test_file.bed --dir <PATH>/chromsig/v0.0.3/Data_Directory/ --r hg38_test --fdr 0.1 --num 10 --type PE
```

## Citation
"Chrom-Sig: de-noising 1-dimensional genomic profiles by signal processing methods" by Nandita J. Gupta, Zachary Apell, and Minji Kim. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.12.670000v1) (2025), 670000.

For questions or bug reports, contact Nandita (njgupta@umich.edu) or visit the "Issues" page.
