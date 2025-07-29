#!/bin/bash
#SBATCH --job-name=freq_enrich_python
#SBATCH --account=minjilab0
#SBATCH --array=0-23 #24 cores for one bam file
#SBATCH --mem=100g
#SBATCH --time=00-01:00:00
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/denoise_output/slurm-%A_%a.out
#SBATCH --error=/nfs/turbo/umms-minjilab/njgupta/chromsig/denoise_outerr/slurm-%A_%a.out
#SBATCH --mail-user=njgupta@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --profile=Task

## The help message:
function usage
{
    echo -e "usage:
    sbatch denoise.sh --bg std_bg --bed std_bed --dir directory --lib lib_name --r ref --fdr fdr --num samp_size --type str_type --plot plot_hist"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bg )          shift
                                std_bg=$1
                                ;;
        -s | --bed )         shift
                                std_bed=$1
                                ;;
        -s | --dir )         shift
                                directory=$1
                                ;;
        -s | --lib )         shift
                                lib_name=$1
                                ;;
        -s | --r )           shift
                                ref=$1
                                ;;
        -s | --fdr )         shift
                                fdr=$1
                                ;;
        -s | --num )         shift
                                samp_size=$1
                                ;;
        -s | --type )        shift
                                str_type=$1
                                ;;
        -s | --plot )        shift
                                plot_hist=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

cd ${directory}

chroms=$(cut -f1 ${ref}.chrom.sizes)

chromtoremove="chrY"

chromslst=()
while read -r line; do chromslst+=("$line"); done <<< "$chroms"
finalchroms=()
for item in "${chromslst[@]}"; do
  # If the current item is not the one to remove, add it to the filtered list
  if [[ "$item" != "$chromtoremove" ]]; then
    finalchroms+=("$item")
  fi
done

#chromslst=(${chroms//$'\n'/ })
#IFS=$' ' read -r -a chromslst2 <<< "$chromslst"

if [ $SLURM_ARRAY_TASK_ID < ${#finalchroms[@]} ]
then
  current_chrom=${finalchroms[$SLURM_ARRAY_TASK_ID]}
  echo $current_chrom

  cd ..
  python chromsig_Freq_Enrich_Test.py "${lib_name} ${ref} ${fdr} ${current_chrom} ${samp_size} ${std_bg} ${directory} ${std_bed} ${str_type} ${plot_hist}"
fi