#!/bin/bash
#SBATCH --job-name=bed_to_stdbed
#SBATCH --account=minjilab0

function usage
{
    echo -e "usage:
    sbatch makestdbed.sh --bp bam_prefix --bed bed_file --dir directory --ref ref"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bp )          shift
                                bam_prefix=$1
                                ;;
        -s | --bed )         shift
                                bed_file=$1
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

echo "Standard Bed file being made: ${bam_prefix}_bam2bed.std.bed";
python Read_Bed_GL.py "${directory} ${bed_file} ${ref}.chrom.sizes";
echo "Standard Bed file made"
