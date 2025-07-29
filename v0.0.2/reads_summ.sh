#!/bin/bash
#SBATCH --job-name=read_summary
#SBATCH --account=minjilab0
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/reads_summ/slurm-%j.out

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

#summ_folder="${directory}/Summary"

for f in ${directory}/*_sorted_bam2bed.bed; do echo "Number of Uniquely Mapped Reads is: $(wc -l ${f})"; done; 

python barplot_summ.py "${directory}/*_sorted_bam2bed.bed";
