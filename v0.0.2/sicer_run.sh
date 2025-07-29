#!/bin/bash -l
#SBATCH --job-name=sicer
#SBATCH --account=minjilab0
#SBATCH --mem=30g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/sicer_output/slurm-%j.out

## The help message:
function usage
{
    echo -e "usage:
    sbatch totalbed.sh --od output_dir --ref ref --bedf bed"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --od )          shift
                                output_dir=$1
                                ;;
        -s | --ref )         shift
                                ref=$1
                                ;;
        -s | --bedf )        shift
                                bed=$1
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
pwd
echo "SICER RUN for ${bed}"

echo "sicer -t ${bed} -s ${ref}"

sicer -t ${bed} -s ${ref}