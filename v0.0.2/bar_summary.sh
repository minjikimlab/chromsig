#!/bin/bash
#SBATCH --job-name=bar_summary
#SBATCH --account=minjilab0
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/bar_summ/slurm-%j.out

function usage
{
    echo -e "usage:
    sbatch bar_summary.sh --dir directory"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --dir )         shift
                                directory=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

summ_folder="${directory}/Summary"

python barplot_summ.py "${summ_folder}";
