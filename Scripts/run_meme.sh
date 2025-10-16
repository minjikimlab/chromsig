#!/bin/bash
#SBATCH --job-name=ctcf_meme
#SBATCH --output=output/meme/meme_%A_%a.out
#SBATCH --error=error/meme/meme_%A_%a.err
#SBATCH --time=05:00:00
#SBATCH --mem=10G
#SBATCH --account=minjilab0
#SBATCH --partition=standard
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --profile=Task
#SBATCH --array=1-12

source ~/.bashrc
conda activate chromsig

module load Bioinformatics
module load meme/5.5.5
module load bedtools2


AW_dir='/nfs/turbo/umms-minjilab/zapell/chrom-sig/data/AtacWorks_sorted'
CS_dir='/nfs/turbo/umms-minjilab/zapell/chrom-sig/data/CS_sorted'

files=("$AW_dir/AW_GM12878_CTCF_ChIP1_ENCFF355CYX_sorted.bed"
       "$AW_dir/AW_GM12878_CTCF_ChIP1_ENCFF453EJM_sorted.bed"
       "$AW_dir/AW_GM12878_CTCF_ChIP2_ENCFF500UFQ_sorted.bed"
       "$AW_dir/AW_GM12878_CTCF_ChIP2_ENCFF668ROZ_sorted.bed"
       "$AW_dir/AW_GM12878_CTCF_CUTRUN1_DR4_sorted.bed"
       "$AW_dir/AW_GM12878_CTCF_CUTRUN2_IB4_sorted.bed"
       "$CS_dir/GM12878_CTCF_ChIP1_ENCFF355CYX_sorted.bed"
       "$CS_dir/GM12878_CTCF_ChIP1_ENCFF453EJM_sorted.bed"
       "$CS_dir/GM12878_CTCF_ChIP2_ENCFF500UFQ_sorted.bed"
       "$CS_dir/GM12878_CTCF_ChIP2_ENCFF668ROZ_sorted.bed"
       "$CS_dir/GM12878_CTCF_CUTRUN1_DR4_sorted.bed"
       "$CS_dir/GM12878_CTCF_CUTRUN1_IB4_sorted.bed"
)


BG_FILE="${files[$SLURM_ARRAY_TASK_ID-1]}"
echo "Processing file: $BG_FILE"
MEME_DIR="/nfs/turbo/umms-minjilab/zapell/chrom-sig/meme_output"

mkdir -p ${MEME_DIR}
sample=$(basename "${BG_FILE}" _sorted.bed)
MEME_OUT_DIR="${MEME_DIR}/${sample}"

GENOME_FA="/nfs/turbo/umms-minjilab/ryan_lab/pipeline/ChIPseq_ATACseq_pipelines/reference_genomes/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
BED_FILE="${BG_FILE%_sorted.bed}.bed"
FASTA_FILE="${BG_FILE%_sorted.bed}.fasta"

# If BED file doesn't exist, create it
if [ ! -f "$BED_FILE" ]; then
    cut -f1-3 "$BG_FILE" > "$BED_FILE"
fi

# If FASTA file doesn't exist, create it
if [ ! -f "$FASTA_FILE" ]; then
    bedtools getfasta -fi "$GENOME_FA" -bed "$BED_FILE" -fo "$FASTA_FILE"
fi

meme-chip -oc "$MEME_OUT_DIR" -time 240 -ccut 100 -dna -order 2 -minw 6 -maxw 15 \
 -db data/motif_dbs/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
 -meme-mod zoops -meme-nmotifs 5 -meme-searchsize 100000 -streme-pvt 0.05 -streme-align center \
 -streme-totallength 4000000 -centrimo-score 5.0 -centrimo-ethresh 10.0 \
  "$FASTA_FILE"
