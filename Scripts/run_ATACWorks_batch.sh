#!/bin/bash
#SBATCH --job-name=atacWorks
#SBATCH --output=output/atacWorks_%A_%a.out
#SBATCH --error=error/atacWorks_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --mem=60G
#SBATCH --account=rjhryan0
#SBATCH --partition=gpu
#SBATCH --mail-user=zapell@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --profile=Task
#SBATCH --gres=gpu:1
#SBATCH --array=0-13   # <-- Set this to 0-(N-1) where N is the number of files

module load cuda

source ~/.bashrc
conda activate atacWorks

PRIM_DIR='/nfs/turbo/umms-minjilab/zapell/chrom-sig/AtacWorks'
cd $PRIM_DIR
MODEL_DIR="${PRIM_DIR}/models"

# download pretrained model - only needed once
# wget -P $MODEL_DIR --content-disposition 'https://api.ngc.nvidia.com/v2/models/org/nvidia/atac_bulk_lowcov_20m_50m/0.3/files?redirect=true&path=models/model.pth.tar' -O $MODEL_DIR/model.pth.tar
# low qual data
wget -P $MODEL_DIR --content-disposition 'https://api.ngc.nvidia.com/v2/models/org/nvidia/atac_bulk_lowqual_20m_20m/0.3/files?redirect=true&path=models/model.pth.tar' -O $MODEL_DIR/low_qual_model.pth.tar

ATAC_DIR='/nfs/turbo/umms-minjilab/njgupta/chromsig/ATAC-seq_PE'
CHIP_DIR='/nfs/turbo/umms-minjilab/njgupta/chromsig/ChIP-seq_SE'
CUTRUN_DIR='/nfs/turbo/umms-minjilab/njgupta/chromsig/CUTandRUN_PE'

# atac reps, ctcf rep pairs, yy1 rep pairs, cutrun ctcf reps
FILES=(
  "$ATAC_DIR/GM12878_ATAC-seq_ENCFF646NWY/GM12878_ATAC-seq_ENCFF646NWY_bed2bg.std.bedgraph"
  "$ATAC_DIR/GM12878_ATAC-seq_ENCFF415FEC/GM12878_ATAC-seq_ENCFF415FEC_bed2bg.std.bedgraph"
  "$ATAC_DIR/K562_ATAC-seq_ENCFF512VEZ/K562_ATAC-seq_ENCFF512VEZ_bed2bg.std.bedgraph"
  "$ATAC_DIR/K562_ATAC-seq_ENCFF987XOV/K562_ATAC-seq_ENCFF987XOV_bed2bg.std.bedgraph"

  "$CHIP_DIR/GM12878_CTCF_ChIP-seq_ENCFF355CYX/GM12878_CTCF_ChIP-seq_ENCFF355CYX_bed2bg.std.bedgraph"
  "$CHIP_DIR/GM12878_CTCF_ChIP-seq_ENCFF453EJM/GM12878_CTCF_ChIP-seq_ENCFF453EJM_bed2bg.std.bedgraph"

  "$CHIP_DIR/GM12878_CTCF_ChIP-seq_ENCFF500UFQ/GM12878_CTCF_ChIP-seq_ENCFF500UFQ_bed2bg.std.bedgraph"
  "$CHIP_DIR/GM12878_CTCF_ChIP-seq_ENCFF668ROZ/GM12878_CTCF_ChIP-seq_ENCFF668ROZ_bed2bg.std.bedgraph"

  "$CHIP_DIR/GM12878_YY1_ChIP-seq_ENCFF481OXC/GM12878_YY1_ChIP-seq_ENCFF481OXC_bed2bg.std.bedgraph"
  "$CHIP_DIR/GM12878_YY1_ChIP-seq_ENCFF676OVC/GM12878_YY1_ChIP-seq_ENCFF676OVC_bed2bg.std.bedgraph"

  "$CHIP_DIR/GM12878_YY1_ChIP-seq_ENCFF301QON/GM12878_YY1_ChIP-seq_ENCFF301QON_bed2bg.std.bedgraph"
  "$CHIP_DIR/GM12878_YY1_ChIP-seq_ENCFF493HLV/GM12878_YY1_ChIP-seq_ENCFF493HLV_bed2bg.std.bedgraph"

  "$CUTRUN_DIR/GM12878_CUTandRUN_CTCF_4DNFI2G71DR4/GM12878_CUTandRUN_CTCF_4DNFI2G71DR4_bed2bg.std.bedgraph"
  "$CUTRUN_DIR/GM12878_CUTandRUN_CTCF_4DNFI9U71IB4/GM12878_CUTandRUN_CTCF_4DNFI9U71IB4_bed2bg.std.bedgraph"

)

REF_GENOME="/nfs/turbo/umms-minjilab/zapell/chrom-sig/hg38.chrom.sizes"
CHROM_SIZES="/nfs/turbo/umms-minjilab/ryan_lab/chromsig_package/Data_Directory/hg38.chrom.sizes"

BW_DIR="/nfs/turbo/umms-minjilab/zapell/chrom-sig/AtacWorks/bigwig"
mkdir -p "$BW_DIR"

# Select the file for this array task
f="${FILES[$SLURM_ARRAY_TASK_ID]}"
BASE=$(basename "$f" .bedgraph)
B2=$(basename "$BASE" .std)

echo "Processing file: $f"
echo "BASE: $BASE"
echo "Current working directory: $(pwd)"

# Uncomment and run if needed
bedGraphToBigWig "$f" "$CHROM_SIZES" "$BW_DIR/${BASE}.bw"

# Uncomment to run atacworks denoise
atacworks denoise \
  --noisybw "$BW_DIR/${BASE}.bw" \
  --genome $REF_GENOME \
  --weights_path ./models/low_qual_model.pth.tar \
  --out_home "./" \
  --exp_name "atacworks_lq_model_denoise_$BASE" \
  --distributed \
  --num_workers 0

# Find the output directory (timestamped or _latest)
OUTDIR=$(ls -d atacworks_lq_model_denoise_${BASE}_* | sort | tail -n 1)

# Run peaksummary.py if output files exist
if [[ -d "$OUTDIR" ]]; then
    python ./scripts/peaksummary.py \
      --peakbw "${OUTDIR}/${B2}_infer.peaks.bw" \
      --trackbw "${OUTDIR}/${B2}_infer.track.bw" \
      --prefix "${B2}_infer.peak_calls" \
      --out_dir "${OUTDIR}" \
      --minlen 20
else
    echo "Output directory $OUTDIR not found for $BASE"
fi