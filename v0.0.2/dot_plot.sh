#!/bin/bash -l
#SBATCH --job-name=dot_plot
#SBATCH --account=minjilab0
#SBATCH --time=00-04:00:00
#SBATCH --mem=30g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/dot_plot_output/slurm-%j.out

ml python/3.9
cd ..
source pyBedGraph/bin/activate
cd chromsig

## The help message:
function usage
{
    echo -e "usage:
    sbatch dotplot.sh --od output_dir --pref pass_score_pref --bpref bam_prefix --lib lib_name"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

#echo "python venn_diag.py '${filea} ${fileb}'"
python dot_plot.py "ATAC-seq_PE ChIP-seq_PE CUTandRUN_PE snATAC-seq_PE"