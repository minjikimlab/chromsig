#!/bin/bash
#SBATCH --job-name=chromHMM
#SBATCH --output=output/chromHMM_AW_%A_%a.out
#SBATCH --error=error/chromHMM_AW_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=10G
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --profile=Task
#SBATCH --array=0-1

source ~/.bashrc
conda activate chromsig

module load Bioinformatics
module load bedtools2

OUTDIR="chromHMM_annotation"
mkdir -p "$OUTDIR"

# annotation paths for each cell line
GM12878_ANNOT='BASIC_chromHMM_GM12878_hg38_Broad_UCSC.bed'
K562_ANNOT='BASIC_chromHMM_K562_hg38_Broad_UCSC.txt'

# path to peak files
# this directory has K562 RNA Pol2 ChIP-seq and GM12878 ATAC-seq data for Chrom-Sig SICER peaks and original SICER peaks
K562_FILES=(chromHMM_data/K562*SICER.bedgraph)
ATAC_FILES=(chromHMM_data/*ATAC*SICER.bedgraph)
# do this for either directory: K562 or ATAC
#F="${ATAC_FILES[$SLURM_ARRAY_TASK_ID]}"
#FNAME=$(basename "$F")

# AtacWorks peak path: only the chrom,start,end cols needed so bed file works here
AW_FILES=(data/AtacWorks_sorted/*GM12878*ATAC*.bed)

F="${AW_FILES[$SLURM_ARRAY_TASK_ID]}"
FNAME=$(basename "$F")

# annotate with chromHMM
bedtools intersect -wao -a <(cut -f1-3 "$F") -b "$GM12878_ANNOT" \
| cut -f1,2,3,7,8 \
| awk -v OFS="\t" '{ print $1";"$2";"$3";"$4, $5 }' \
| awk '{soma[$1]+=$2} END {for (item in soma) print item, soma[item]}' \
| sort \
| awk -v OFS="\t" '{print $1";"$2}' \
| awk -v OFS="\t" -F';' '{print $1,$2,$3,$4,$5}' \
> "$OUTDIR/${FNAME%.bed}_chromhmmanno.bed"
