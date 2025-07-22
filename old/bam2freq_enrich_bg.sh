#!/bin/bash 

## The help message:
function usage
{
    echo -e "usage:
    bash bam2freq_enrich_bg.sh --bam bam_file --dir directory --lib lib_name --r ref --fdr fdr --num samp_size --type str_type"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bam )         shift
                                bam_file=$1
                                ;;
        -s | --dir )         shift
                                directory=$1
                                ;;
        -s | --lib )         shift
                                lib_name=$1
                                ;;
        -s | --r )           shift
                                ref=$1
                                ;;
        -s | --fdr )         shift
                                fdr=$1
                                ;;
        -s | --num )         shift
                                samp_size=$1
                                ;;
        -s | --type )        shift
                                str_type=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

#module load gcc/4.9.2
#module load bedtools/2.27.0

cd ${directory}
#source ${bam_file}

bam_prefix=$(basename ${bam_file} .bam)

# convert bam to bed
if [ ! -f "${bam_prefix}_bam2bed.bed" ]
then
    echo "Bed file being made: ${bam_prefix}_bam2bed.bed";
    bedtools bamtobed -i ${bam_file} > "${bam_prefix}_bam2bed.bed";
    echo "Bed file made"
fi
bed_file="${bam_prefix}_bam2bed.bed"

# convert bed to standard bed
if [ ! -f "${bam_prefix}_bam2bed.std.bed" ]
then
    echo "Standard Bed file being made: ${bam_prefix}_bam2bed.std.bed";
    cd ..;
    pwd;
    NB_ARGS="${directory} ${bed_file} ${ref}.chrom.sizes" jupyter nbconvert --to notebook --execute Read_Bed.ipynb;
    echo "Standard Bed file made"
fi
std_bed="${bam_prefix}_bam2bed.std.bed"

cd ${directory}
pwd

# convert bed to bedgraph
if [ ! -f "${bam_prefix}_bed2bg.std.bedgraph" ]
then
    echo "Standard Bedgraph file being made: ${bam_prefix}_bed2bg.std.bedgraph";
    bedtools genomecov -i ${std_bed} -g "${ref}.chrom.sizes" -bg > "${bam_prefix}_bed2bg.std.bedgraph";
    echo "Standard Bedgraph file made"
fi
std_bg="${bam_prefix}_bed2bg.std.bedgraph"

# run jupyter notebook
cd ..
pwd
NB_ARGS="${lib_name} ${ref} ${fdr} $(cut -f1 ${directory}${ref}.chrom.sizes) ${samp_size} ${std_bg} ${directory} ${bed_file} ${str_type}" jupyter nbconvert --to notebook --execute Freq_Enrich_Sigtest.ipynb

output_dir="${directory}${lib_name}_EnrichTest_FDR_${fdr}/"
cd ${output_dir}

pwd

# write out total result (all chromosomes)
#if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_enrichTest_master.txt" ]
#then
#    cat $(sed '1d' *_master.txt) > "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_enrichTest_master.txt"
#fi


# output pileups to bedgraph files
if [ ! -f "*.bedgraph" ]
then
    for file in *.bed; do
        file_pref=$(basename ${file} .bed);
        echo ${file_pref};
        bedtools genomecov -i ${file} -g "${directory}/${ref}.chrom.sizes" -bg > "${file_pref}_bed2bg.bedgraph";
    done
fi

if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_pass_pileup.bed" ]
then
    cat *_pass_pileup.bed > "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_pass_pileup.bed";
    total_pass="${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_pass_pileup.bed"
fi

if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_fail_pileup.bed" ]
then
    cat *_fail_pileup.bed > "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_fail_pileup.bed";
    total_fail="${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_fail_pileup.bed"
fi

if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_fail_pileup_bed2bg.bedgraph" ]
then
    cat *_fail_pileup_bed2bg.bedgraph > "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_fail_pileup_bed2bg.bedgraph"
fi

if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_pass_pileup_bed2bg.bedgraph" ]
then
    cat *_pass_pileup_bed2bg.bedgraph > "${lib_name}_total_FDR_${fdr}_pseudoGEM_${samp_size}_pass_pileup_bed2bg.bedgraph"
fi

pass_lines=$(wc -l < ${total_pass})
fail_lines=$(wc -l < ${total_fail})
percentage=100
echo "Pass Percentage: "
echo "scale=4 ; $pass_lines / ($pass_lines + $fail_lines) * $percentage" | bc
echo "\nFail Percentage: "
echo "scale=4 ; $fail_lines / ($pass_lines + $fail_lines) * $percentage" | bc
