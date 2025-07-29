#!/bin/bash -l
#SBATCH --job-name=venn_diag
#SBATCH --account=minjilab0
#SBATCH --time=00-04:00:00
#SBATCH --mem=30g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/venn_diag_output/slurm-%j.out

## The help message:
function usage
{
    echo -e "usage:
    sbatch venn_diag.sh --od output_dir --pref pass_score_pref --bpref bam_prefix --lib lib_name"
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
        -s | --lib )         shift
                                lib_name=$1
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

filea="${directory}/${bam_prefix}_bam2bed.std-W200-G600.scoreisland"
fileb="${output_dir}${pass_score_pref}-W200-G600.scoreisland"

if [ ! -d "${directory}/venn_diagram" ]
then
    mkdir "${directory}/venn_diagram"
    cp ${filea} "${directory}/venn_diagram/"
    cp ${fileb} "${directory}/venn_diagram/"
fi

#intervene venn -i ${directory}/venn_diagram/*.scoreisland --output "${directory}/venn_diagram" --names="Pass Pileup Peaks","Total Bedgraph Peaks" --colors=g,r --figtype png > "${directory}/venn_diagram/${lib_name}_venn.png"

# Finding box plot

# Only A
bedtools intersect -a ${filea} -b ${fileb} -v > ${directory}/venn_diagram/only_a.scoreisland

# Only B
bedtools intersect -a ${fileb} -b ${filea} -v > ${directory}/venn_diagram/only_b.scoreisland

# Intersect A and B
bedtools intersect -a ${filea} -b ${fileb} > ${directory}/venn_diagram/a_intersect_b.scoreisland

onlya="${directory}/venn_diagram/only_a.scoreisland"
onlyb="${directory}/venn_diagram/only_b.scoreisland"
a_and_b="${directory}/venn_diagram/a_intersect_b.scoreisland"

cd ..
#echo "python venn_diag.py '${filea} ${fileb}'"
python venn_diag.py "${onlya} ${onlyb} ${a_and_b} ${directory}/venn_diagram/ ${lib_name}"