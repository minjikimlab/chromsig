# Chrom-Sig: de-noising 1-dimensional genomic profiles by signal processing methods

All processed data generated using chromsig v0.0.2 and v0.0.3 can be found at [this link](https://www.dropbox.com/scl/fo/dri58rnvghvbswzks5ox6/AFz6TKN7NbiASQMMI3HHo8c?rlkey=jel3nvo6azx5t89w76jxa5cuj&st=mcj5touo&dl=0).

## Features
* Chrom-Sig takes in Single-End or Paired-End read .bam or .bed files generated from ChIP-seq, ATAC-seq, snATAC-seq, or CUT&RUN experiments, among others 
* As an output, Chrom-Sig produces "pass" reads .bed file and .bedgraph files, along with peaks called by SICER on "pass" reads

## Installation
This program requires Python 3.9 to be able to run.
Start by installing Chrom-Sig using the following command:
```
git clone https://github.com/minjikimlab/chromsig.git
```
Navigating to the directory `chromsig`, do the following:

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
Navigate to the directory `chromsig/v0.0.4`.

Install all required packages listed in [requirements.txt](https://github.com/minjikimlab/chromsig/blob/main/v0.0.4/requirements.txt) using the following command:

`pip install -r requirements.txt`

If you are using conda, use the following command:

`conda install --file requirements.txt`

In each of the bash files (files ending in .sh), modify the parameters
```
#SBATCH --account=<ACCOUNT_NAME>
#SBATCH --output=<OUTPUT_FILE_PATH>/slurm-%j.out
#SBATCH --error=<ERROR_FILE_PATH>/slurm-%j.out
#SBATCH --mail-user=<USER_EMAIL>
```
to match the account using which you are submitting the job, the output directory for the job report, the output directory for any error reports, and the email address you want the job report to be mailed to.

## Running Chrom-Sig
In the same directory `chromsig/v0.0.4`, create a folder called Data_Directory. In this folder, upload your input .bam or .bed file, and your reference genome sizes file (ex: [hg38.chrom.sizes](https://github.com/minjikimlab/chromsig/blob/main/v0.0.4/hg38.chrom.sizes)).
Return to `chromsig/v0.0.4`, and in your terminal enter the following command:

`sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed <BAM/BED> --dir <PATH_TO_DATA_DIR> --r <GENOME> --fdr <FDR_THRESH> --num <NUM_OF_PSEUDO_SAMPLES> --type <READ_TYPE> --cov <COVERAGE> --lib <LIBRARY> --plot <TRUE> --modload <TRUE>`

The following is a description of each of the arguments in the command:

ARGUMENTS:

    --bamorbed      Input bam or bed file containing your dataset.
    --dir           Directory containing the bam/bed file, as well as your reference genome sizes file (ex: /Users/chromsig/Data/).
    --r             Name of reference genome (ex: hg38).
    --fdr           FDR threshold for denoising (ex: 0.1).
    --num           Number of pseudo samples to be generated for denoising.
    --type          Type of read in your dataset, Single End or Paired End (ex: SE or PE).
    --cov           Minimum read coverage required for a chromsome to be denoised (ex.: To set minimum coverage for a chromsome to be 5%, enter --cov 0.05)

OPTIONS:
      
    --lib           Name of dataset library (ex: For GM12878_ATAC-seq_ENCFF415FEC.bam, this would be ENCFF415FEC).
    --plot          Produce histogram plots of distances between read fragments (only for Paired End), type 'true' if you want plots, skip this option if you do not.
    --modload       Load python 3.9, bedtools, and samtools; type 'true' if you want to use module load to load these, skip this option if you do not.
            
EXAMPLE:

    sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed GM12878_ATAC-seq_ENCFF415FEC.bam --dir <FILE_PATH>/chromsig/ATAC-seq_PE/ --r hg38 --fdr 0.1 --num 5000 --type PE --cov 0.05 --modload true

## Recommended parameters
We recommend allocating about an hour of runtime and 25 GB of memory per 50 million paired-end reads. The results presented in the manuscript have been produced by running Chrom-Sig with --num 5000 --fdr 0.1 or --num 5000 --fdr 0.2.  

## Testing Chrom-Sig
Two test files have been included in the `chromsig` directory -- the bed file contains test data from chr18 of the GM12878_CTCF_ChIP-seq_ENCFF453EJM.bam file, and the hg38_test.chrom.sizes file contains the size of chr18 needed to run this data. Move both of these files to your Data Directory folder `chromsig/v0.0.4/Data_Directory`. Navigate back to `chromsig/v0.0.4`, and run the following command:
```
sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed chromsig_test_file.bed --dir <PATH>/chromsig/v0.0.4/Data_Directory/ --r hg38_test --fdr 0.1 --num 10 --type SE --modload true
```
A folder Chromsig_Test_Results.zip is included in [this dropbox folder](https://www.dropbox.com/scl/fo/dri58rnvghvbswzks5ox6/AFz6TKN7NbiASQMMI3HHo8c?rlkey=jel3nvo6azx5t89w76jxa5cuj&st=mcj5touo&dl=0) -- the Chromsig_Test_Results.zip folder contains the results of running the test files.

## Citation
"Chrom-Sig: de-noising 1-dimensional genomic profiles by signal processing methods" by Nandita J. Gupta, Zachary Apell, and Minji Kim. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.08.12.670000v1) (2025), 670000.

For questions or bug reports, contact Nandita (njgupta@umich.edu) or visit the "Issues" page.
