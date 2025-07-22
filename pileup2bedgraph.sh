#!/bin/bash -l
#SBATCH --job-name=pileup_to_bedgraph
#SBATCH --account=minjilab0
#SBATCH --array=0-23 #24 cores for one bam file
#SBATCH --mem=20g

## The help message:
function usage
{
    echo -e "usage:
    sbatch pileup2bedgraph.sh --od output_dir --dir directory --lib lib_name --r ref --fdr fdr --num samp_size --type str_type"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --od )          shift
                                output_dir=$1
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
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

cd ${output_dir}

chroms=$(cut -f1 ${directory}${ref}.chrom.sizes)

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
echo ${finalchroms[@]}
echo ${finalchroms[23]}

#chromslst=(${chroms//$'\n'/ })
#IFS=$' ' read -r -a chromslst2 <<< "$chromslst"
current_chrom=${finalchroms[$SLURM_ARRAY_TASK_ID]}

file_name_pass="${lib_name}_${current_chrom}_FDR_${fdr}_pseudoGEM_${samp_size}${str_type}_pass_pileup.bed"
file_name_fail="${lib_name}_${current_chrom}_FDR_${fdr}_pseudoGEM_${samp_size}${str_type}_fail_pileup.bed"

bedtools genomecov -i ${file_name_pass} -g "${directory}/${ref}.chrom.sizes" -bg > "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph"
bedtools genomecov -i ${file_name_fail} -g "${directory}/${ref}.chrom.sizes" -bg > "$(basename ${file_name_fail} .bed)_bed2bg.bedgraph"