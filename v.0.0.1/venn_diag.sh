#!/bin/bash -l
#SBATCH --job-name=venn_diag
#SBATCH --account=minjilab0
#SBATCH --time=00-04:00:00
#SBATCH --mem=30g

## The help message:
function usage
{
    echo -e "usage:
    sbatch totalbed.sh --od output_dir --pref pass_score_pref --bpref bam_prefix"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --od )          shift
                                output_dir=$1
                                ;;
        -s | --pref )        shift
                                pass_score_pref=$1
                                ;;
        -s | --bpref )       shift
                                bam_prefix=$1
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
cd ..
directory=$(pwd)

#seta=$(awk '{print $1, $2, $3}' "${bam_prefix}_bam2bed-W200-G600.scoreisland")
#setb=$(awk '{print $1, $2, $3}' "${output_dir}${pass_score_pref}-W200-G600.scoreisland")

#setastr="${seta[*]}"
#setbstr="${setb[*]}"

filea="${directory}/${bam_prefix}_bam2bed-W200-G600.scoreisland"
fileb="${output_dir}${pass_score_pref}-W200-G600.scoreisland"

mkdir "${directory}/venn_diagram"
cp ${filea} "${directory}/venn_diagram/"
cp ${fileb} "${directory}/venn_diagram/"

intervene venn -i ${directory}/venn_diagram/*.scoreisland --output "${directory}/venn_diagram"


#cd ..
#echo "python venn_diag.py '${filea} ${fileb}'"
#python venn_diag.py "${filea} ${fileb}"
