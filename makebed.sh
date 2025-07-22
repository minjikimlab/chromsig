#!/bin/bash
#SBATCH --job-name=bam_to_bed
#SBATCH --account=minjilab0

function usage
{
    echo -e "usage:
    sbatch makebed.sh --bp bam_prefix --bam bam_file --dir directory"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bp )          shift
                                bam_prefix=$1
                                ;;
        -s | --bam )         shift
                                bam_file=$1
                                ;;
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

cd ${directory}

echo "Bed file being made: ${bam_prefix}_bam2bed.bed";
bedtools bamtobed -i ${bam_file} > "${bam_prefix}_bam2bed.bed";
echo "Bed file made"