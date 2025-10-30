#!/bin/bash -l
#SBATCH --job-name=total_bed
#SBATCH --account=minjilab0
#SBATCH --array=0-1

## The help message:
function usage
{
    echo -e "usage:
    sbatch totalbed.sh --od output_dir --lib lib_name --fdr fdr --num samp_size"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --od )          shift
                                output_dir=$1
                                ;;
        -s | --lib )         shift
                                lib_name=$1
                                ;;
        -s | --fdr )         shift
                                fdr=$1
                                ;;
        -s | --num )         shift
                                samp_size=$1
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

# total all chromosomes into whole genome files. To reduce time, this totals the pass files in parallel with the fail files.

bed_ids=("pass" "fail")
current_id=${bed_ids[$SLURM_ARRAY_TASK_ID]}

cat *_${current_id}_pileup.bed > "${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_${current_id}_pileup.bed";

cat *_${current_id}_pileup_bed2bg.bedgraph > "${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_${current_id}_pileup_bed2bg.bedgraph"