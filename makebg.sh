#!/bin/bash
#SBATCH --job-name=stdbed_to_bg
#SBATCH --account=minjilab0
#SBATCH --mem=80g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/makebg_output/slurm-%j.out

function usage
{
    echo -e "usage:
    sbatch makebg.sh --bp bam_prefix --stdbed std_bed --dir directory --ref ref"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bp )          shift
                                bam_prefix=$1
                                ;;
        -s | --stdbed )      shift
                                std_bed=$1
                                ;;
        -s | --dir )         shift
                                directory=$1
                                ;;
        -s | --ref )         shift
                                ref=$1
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

pwd

echo "Standard Bedgraph file being made: ${bam_prefix}_bed2bg.std.bedgraph";
echo "bedtools genomecov -i ${std_bed} -g "${ref}.chrom.sizes" -bg > "${bam_prefix}_bed2bg.std.bedgraph"";
bedtools genomecov -i ${std_bed} -g "../${ref}.chrom.sizes" -bg > "${bam_prefix}_bed2bg.std.bedgraph";
echo "Standard Bedgraph file made"

pwd