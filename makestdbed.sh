#!/bin/bash
#SBATCH --job-name=bed_to_stdbed
#SBATCH --account=minjilab0
#SBATCH --time=00-24:00:00
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/makestdbed_output/slurm-%j.out

function usage
{
    echo -e "DESCRIPTION:
    This script sorts the .bed file, and then runs it through the python program Read_Bed.py to remove any non-standard chromosomes.
    USAGE:
    sbatch makestdbed.sh --bp bam_prefix --bed bed_file --dir directory --ref ref
    
        ARGUMENTS:
            --bp    Name of the bam file with .bam extension removed (ex: For GM12878_ChIP-seq_RAD51_ENCFF725MWO.bam this would be GM12878_ChIP-seq_RAD51_ENCFF725MWO).
            --bed   Name of the bed file.
            --dir   Directory in which the bam file is stored.
            --ref   Name of the reference genome (ex: hg38).
        
        EXAMPLE:
            sbatch makestdbed.sh --bp GM12878_ChIP-seq_RAD51_ENCFF725MWO --bam GM12878_ChIP-seq_RAD51_ENCFF725MWO_bam2bed.bed --dir /Users/chromsig/Data/bam_file_name/ --ref hg38"
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

echo "Sorting Bed"
sort -k 1,1 -k2,2n ${directory}${bam_prefix}_bam2bed.bed > ${directory}${bam_prefix}_sorted_bam2bed.bed
echo "Bed Sorted"
sorted_bed="${bam_prefix}_sorted_bam2bed.bed"
echo "Standard Bed file being made: ${bam_prefix}_bam2bed.std.bed";
python Read_Bed.py "${directory} ${sorted_bed} ${ref}.chrom.sizes";
echo "Standard Bed file made"