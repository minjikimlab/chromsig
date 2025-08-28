#!/bin/bash -l
#SBATCH --job-name=pileup_to_bedgraph
#SBATCH --account=minjilab0
#SBATCH --array=0-23 #24 cores for one bam file
#SBATCH --mem=20g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/pileup2bedgraph_output/slurm-%j.out

## The help message:
function usage
{
    echo -e "usage:
    sbatch pileup2bedgraph.sh --od output_dir --dir directory --lib lib_name --r ref --fdr fdr --num samp_size --type str_type --bed std_bed"
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
        -s | --bed )         shift
                                std_bed=$1
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

coverage_arr=$(cut -f1 ${directory}${lib_name}_final_coverage.txt)
covlst=()
while read -r line; do covlst+=("$line"); done <<< "$coverage_arr"

chroms=$(cut -f1 ${directory}${ref}.chrom.sizes)

chromtoremove="chrY"

chromslst=()
while read -r line; do chromslst+=("$line"); done <<< "$chroms"
finalchroms=()
for item in "${!chromslst[@]}"; do
  echo "Current chrom: ${chromslst[$item]}"
  # If the current item is not the one to remove, add it to the filtered list
  if [[ "${chromslst[$item]}" != "$chromtoremove" ]]
  then
    if grep -qFw "${chromslst[$item]}" "${directory}${std_bed}"
    then
      echo ${covlst[$item]}
      num_to_compare=$(echo ${covlst[$item]} | cut -d ' ' -f2)
      echo "Coverage: $num_to_compare"
      if (( $(echo "$num_to_compare > 0.20" | bc -l) )); then
        finalchroms+=("${chromslst[$item]}")
      fi
    fi
  fi
done

echo "Chroms being processed:"
echo ${finalchroms[@]}


#chromslst=(${chroms//$'\n'/ })
#IFS=$' ' read -r -a chromslst2 <<< "$chromslst"
if (( $SLURM_ARRAY_TASK_ID < ${#finalchroms[@]} )); then
  current_chrom=${finalchroms[$SLURM_ARRAY_TASK_ID]}

  file_name_pass="${lib_name}_${current_chrom}_FDR_${fdr}_pseudoRead_${samp_size}${str_type}_pass_pileup.bed"
  file_name_fail="${lib_name}_${current_chrom}_FDR_${fdr}_pseudoRead_${samp_size}${str_type}_fail_pileup.bed"

  bedtools genomecov -i ${file_name_pass} -g "${directory}/${ref}.chrom.sizes" -bg > "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph"
  bedtools genomecov -i ${file_name_fail} -g "${directory}/${ref}.chrom.sizes" -bg > "$(basename ${file_name_fail} .bed)_bed2bg.bedgraph"
fi