#!/bin/bash
#SBATCH --account=minjilab0
#SBATCH --time=00-24:00:00
#SBATCH --mem=30g
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/overall_program_output/slurm-%j.out

chromsig_dir=$(pwd)
cd ..
if [ ! -d "pyBedGraph" ]
then
    git clone https://github.com/TheJacksonLaboratory/pyBedGraph.git
    python -m venv pyBedGraph
fi

source pyBedGraph/bin/activate

## The help message:
function usage
{
    echo -e "DESCRIPTION:
    The purpose of this command is to denoise a 1-dimensional chromatin profiles. This data must be in single-end or paired-end format. The command takes a .bam or .bed file as input, and outputs 
    denoised .bed, .bedgraph, and .scoreisland files for the full genome (.scoreisland is the result of the SICER peak-calling algorithm).
    USAGE:
    sbatch multijob_bam2freq_enrich_bg.sh --bamorbed input_file --dir data_directory --r ref --fdr fdr --num samp_size --type str_type --lib lib_name --plot true --modload true

        ARGUMENTS:
            --bamorbed      Input bam or bed file containing your dataset.
            --dir           Directory containing the bam/bed file, as well as your reference genome sizes file (ex: /Users/chromsig/Data/).
            --r             Name of reference genome (ex: hg38).
            --fdr           FDR threshold for denoising (ex: 0.1).
            --num           Number of pseudo samples to be generated for denoising.
            --type          Type of read in your dataset, Single End or Paired End (ex: SE or PE).

        OPTIONS:
            --lib           Name of dataset library (ex: For GM12878_ATAC-seq_ENCFF415FEC.bam, this would be ENCFF415FEC).
            --plot          Produce histogram plots of distances between read fragments (only for Paired End), type 'true' if you want plots, skip this option if you do not.
            --modload       Load python 3.9, bedtools, and samtools; type 'true' if you want to use module load to load these, skip this option if you do not.
            
        EXAMPLE:
            sbatch multiscript_bam2freq_enrich_bg.sh --bamorbed GM12878_ATAC-seq_ENCFF415FEC.bam --dir /nfs/turbo/umms-minjilab/njgupta/chromsig/ATAC-seq_PE/ --r hg38 --fdr 0.1 --num 5000 --type PE --modload true"
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --bamorbed )    shift
                                input_file=$1
                                ;;
        -s | --dir )         shift
                                data_directory=$1
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
        -s | --lib )         shift
                                lib_name=${1:-false}
                                ;;
        -s | --plot )        shift
                                plot_hist=${1:-false}
                                ;;
        -s | --modload )     shift
                                mod_load=${1:-false}
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# formating time for runtime result

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
}

# if user selected modload, then load necessary modules

if [[ -v mod_load ]]
then
    ml python/3.9
    module load Bioinformatics bedtools2
    module load samtools
fi


cd ${data_directory}

# determine input file type, format name for later use in script

if [[ "${input_file}" == *".bam"* ]]
then
    echo "Input file is bam"
    bam_prefix=$(basename ${input_file} .bam)
elif [[ "${input_file}" == *".bed"* ]]
then
    echo "Input file is bed"
    if [[ "${input_file}" == *"_bam2bed.bed"* ]]
    then
        bam_prefix=$(basename ${input_file} _bam2bed.bed)
    else
        bam_prefix=$(basename ${input_file} .bed)
        mv "${input_file}" "${bam_prefix}_bam2bed.bed"
        input_file="${bam_prefix}_bam2bed.bed"
        ls
    fi
fi

echo "Input file prefix:"
echo ${bam_prefix}

if [[ ! -v lib_name ]]
then
    lib_name=`echo ${bam_prefix} | rev | cut -d '_' -f 1 | rev`
fi

echo ${lib_name}

# move all relevant files associated with input to new directory to store all results for current input

if [ ! -d "${bam_prefix}" ]
then
    mkdir "${bam_prefix}"
fi

if [ ! -f "${data_directory}${bam_prefix}/${input_file}" ]
then
    mv "${input_file}" "${bam_prefix}"
fi

if [ ! -f "${data_directory}${bam_prefix}/${ref}.chrom.sizes" ]
then
    echo "Copying sizes file"
    cp "${ref}.chrom.sizes" "${bam_prefix}"
fi


directory="${data_directory}${bam_prefix}/"
echo ${directory}

cd ${directory}
# convert bam to bed
if [ ! -f "${bam_prefix}_bam2bed.bed" ]
then
    echo "In bam to bed"
    cd ${chromsig_dir}
    jid1=$(sbatch makebed.sh --bp ${bam_prefix} --bam ${input_file} --dir ${directory} --str ${str_type} | awk '{print $4}')
fi
bed_file="${bam_prefix}_bam2bed.bed"

# convert bed to standard bed
cd ${directory}

echo "${bam_prefix}_bam2bed.std.bed"
if [ ! -f "${bam_prefix}_bam2bed.std.bed" ]
then
    echo "In bed to std bed"
    pwd
    if [ ! -f "${bam_prefix}_bam2bed.bed" ]
    then
        cd ${chromsig_dir}
        jid2=$(sbatch --dependency=afterok:$jid1 makestdbed.sh --bp ${bam_prefix} --bed ${bed_file} --dir ${directory} --ref ${ref} | awk '{print $4}')
    else
        cd ${chromsig_dir}
        echo ${chromsig_dir}
        pwd
        jid2=$(sbatch makestdbed.sh --bp ${bam_prefix} --bed ${bed_file} --dir ${directory} --ref ${ref} | awk '{print $4}')
    fi
fi
std_bed="${bam_prefix}_bam2bed.std.bed"
echo ${std_bed}

cd ${directory}
pwd

# convert standard bed to bedgraph
if [ ! -f "${bam_prefix}_bed2bg.std.bedgraph" ]
then
    echo "In bed to bg"
    if [ ! -f "${bam_prefix}_bam2bed.std.bed" ]
    then
        cd ${chromsig_dir}
        jid3=$(sbatch --dependency=afterok:$jid2 makebg.sh --bp ${bam_prefix} --stdbed ${std_bed} --dir ${directory} --ref ${ref} | awk '{print $4}')
    else
        cd ${chromsig_dir}
        jid3=$(sbatch makebg.sh --bp ${bam_prefix} --stdbed ${std_bed} --dir ${directory} --ref ${ref} | awk '{print $4}')
    fi
fi
std_bg="${bam_prefix}_bed2bg.std.bedgraph"


# set names for denoising algorithm output files
cd ${directory}
pwd

output_dir="${directory}${lib_name}_EnrichTest_FDR_${fdr}_pseudoRead_${samp_size}/"
echo ${output_dir}

total_pass="${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_pass_pileup.bed"
total_fail="${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_fail_pileup.bed"

file_name_pass="${lib_name}_chr1_FDR_${fdr}_pseudoRead_${samp_size}${str_type}_pass_pileup.bed"
file_name_fail="${lib_name}_chr1_FDR_${fdr}_pseudoRead_${samp_size}${str_type}_fail_pileup.bed"

# if user selected plot, then plot histograms of read length distribution for each chromosome as part of denoising algorithm
if [[ ! -v plot_hist ]]
then
    plot_hist="False"
else
    plot_hist="True"
fi

# if the denoising results directory has not been created (first time running the program), then make directory and proceed through algorithm from beginning
if [ ! -d ${output_dir} ]
then
    mkdir ${output_dir}
    # run denoising algorithm
    if [ ! -f "${bam_prefix}_bed2bg.std.bedgraph" ]
    then
        cd ${chromsig_dir}
        jid4=$(sbatch --dependency=afterok:$jid3 denoise.sh --cd ${chromsig_dir} --bg ${std_bg} --bed ${std_bed} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --plot ${plot_hist} | awk '{print $4}')
    else
        cd ${chromsig_dir}
        jid4=$(sbatch denoise.sh --cd ${chromsig_dir} --bg ${std_bg} --bed ${std_bed} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --plot ${plot_hist} | awk '{print $4}')
    fi
    cd ${output_dir}
    # convert all output bed files from denoising algorithm into bedgraph files for easier visualization
    if [ ! -f "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph" ]
    then
        #sbatch pileup2bedgraph.sh --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type}
        if [ ! -f ${file_name_pass} ]
        then
            echo "Creating pass bedgraphs after waiting for denoise.sh"
            cd ${chromsig_dir}
            pwd
            jid5=$(sbatch --dependency=afterok:$jid4 pileup2bedgraph.sh --od ${output_dir} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --bed ${std_bed} | awk '{print $4}')
        else
            echo "Creating pass bedgraphs right away"
            cd ${chromsig_dir}
            jid5=$(sbatch pileup2bedgraph.sh --od ${output_dir} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --bed ${std_bed} | awk '{print $4}')
        fi
    fi
    cd ${output_dir}
    # once all bedgraphs have been made, store all chromosomes results together in file for whole genome.The total bedgraph can then be used to visualize the results for the whole genome.
    if [ ! -f ${total_pass} ]
    then
        if [ ! -f "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph" ]
        then
            cd ${chromsig_dir}
            jid6=$(sbatch --dependency=afterok:$jid5 totalbed.sh --od ${output_dir} --lib ${lib_name} --fdr ${fdr} --num ${samp_size} | awk '{print $4}')
        else
            cd ${chromsig_dir}
            jid6=$(sbatch totalbed.sh --od ${output_dir} --lib ${lib_name} --fdr ${fdr} --num ${samp_size} | awk '{print $4}')
        fi
    fi
    cd ${output_dir}
    pass_score_pref=$(basename ${total_pass} .bed)
    # run sicer peak-calling algorithm on both original bed file, and denoised whole genome bed file
    if [ ! -f "${pass_score_pref}-W200-G600.scoreisland" ]
    then
        if [ ! -f ${total_pass} ]
        then
            #sicer -t ${total_pass} -s ${ref}
            cd ${chromsig_dir}
            jid7=$(sbatch --dependency=afterok:$jid6 sicer_run.sh --od ${output_dir} --ref ${ref} --bedf ${total_pass} | awk '{print $4}')
        else
            jid7=$(sbatch sicer_run.sh --od ${output_dir} --ref ${ref} --bedf ${total_pass} | awk '{print $4}')
        fi
    fi
    cd ${directory}
    if [ ! -f "${bam_prefix}_bam2bed.std-W200-G600.scoreisland" ]
    then
        cd ${output_dir}
        if [ ! -f "${pass_score_pref}-W200-G600.scoreisland" ]
        then
            #sicer -t ${total_pass} -s ${ref}
            cd ${chromsig_dir}
            pwd
            echo "Creating total score island file"
            jid8=$(sbatch --dependency=afterok:$jid7 sicer_run.sh --od ${directory} --ref ${ref} --bedf ${std_bed} | awk '{print $4}')
        else
            cd ${chromsig_dir}
            echo "Creating total score island file (bed already made)"
            jid8=$(sbatch sicer_run.sh --od ${directory} --ref ${ref} --bedf ${std_bed} | awk '{print $4}')
        fi
    fi

    # total all file results -- this will produce percentage and runtime stats for the run.
    cd ${output_dir}
    WAIT_INTERVAL=60

    until [ -f "$total_fail" ]; do
        echo "File not found yet, waiting..."
        sleep "$WAIT_INTERVAL"
    done
    cd ${directory}
    pass_lines=$(wc -l < ${output_dir}${total_pass})
    fail_lines=$(wc -l < ${output_dir}${total_fail})
    percentage=100
    echo "Read Type: "
    echo "${str_type}"
    echo "Number of Samples: "
    echo "${samp_size}"
    echo "Denoise Runtime: "
    cd ${output_dir}
    runtime_lines_array=()
    for file in *_logFile.txt; do runtime_lines_array+=( $(grep "Total: " "$file") ); done
    echo ${runtime_lines_array[@]}
    runtime_array=()
    for element in "${runtime_lines_array[@]}"; do substring=$(echo "$element" | cut -d' ' -f2); float_value=$(echo "$substring + 0" | bc -l); runtime_array+=("$float_value"); done
    echo ${runtime_array[@]}
    max_value=$(printf "%s\n" "${runtime_array[@]}" | awk 'BEGIN { max = -infinity } { if ($1 > max) max = $1 } END { print max }')
    echo "${max_value}"
    echo "Pass Percentage: "
    echo "scale=4 ; $pass_lines / ($pass_lines + $fail_lines) * $percentage" | bc
    echo "Pass Reads: "
    echo "${pass_lines}"
    echo "Fail Percentage: "
    echo "scale=4 ; $fail_lines / ($pass_lines + $fail_lines) * $percentage" | bc
    echo "Fail Reads: "
    echo "${fail_lines}"
    echo "Total Reads: "
    echo "$(wc -l < ${std_bed})"
else
    # this is for the condition where a user may be running the algorithm starting from the middle (ex.: ran once before, so have all the starting files, but decided to rerun it). In this case, the denoise output directory will have already been created.
    cd ${output_dir}
    # if the pass results of denoising are not present in the denoise directory, then run denoising.
    if [ ! -f ${file_name_pass} ]
    then
        echo "${file_name_pass} does not exist"
        cd ..
        if [ ! -f "${bam_prefix}_bed2bg.std.bedgraph" ]
        then
            cd ${chromsig_dir}
            jid4=$(sbatch --dependency=afterok:$jid3 denoise.sh --cd ${chromsig_dir} --bg ${std_bg} --bed ${std_bed} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --plot ${plot_hist} | awk '{print $4}')
        else
            echo "Bedgraph already made, line 250"
            cd ${chromsig_dir}
            jid4=$(sbatch denoise.sh --cd ${chromsig_dir} --bg ${std_bg} --bed ${std_bed} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --plot ${plot_hist} | awk '{print $4}')
        fi
    else
        echo "${file_name_pass} already exists"
    fi
    cd ${output_dir}
    # if the bedgraph results of denoising are not present in the denoise directory, then run bed to bedgraph conversion.
    if [ ! -f "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph" ]
    then
        echo "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph does not exist"
        if [ ! -f ${file_name_pass} ]
        then
            cd ${chromsig_dir}
            jid5=$(sbatch --dependency=afterok:$jid4 pileup2bedgraph.sh --od ${output_dir} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --bed ${std_bed} | awk '{print $4}')
        else
            echo "Denoise already run, line 263"
            cd ${chromsig_dir}
            jid5=$(sbatch pileup2bedgraph.sh --od ${output_dir} --dir ${directory} --lib ${lib_name} --r ${ref} --fdr ${fdr} --num ${samp_size} --type ${str_type} --bed ${std_bed} | awk '{print $4}')
        fi
    else
        echo "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph already exists"
    fi
    cd ${output_dir}
    # If the whole genome file results are not present in the denoise directory, then total denoising results.
    if [ ! -f "${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_pass_pileup.bed" ]
    then
        echo "${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_pass_pileup.bed does not exist"
        if [ ! -f "$(basename ${file_name_pass} .bed)_bed2bg.bedgraph" ]
        then
            cd ${chromsig_dir}
            jid6=$(sbatch --dependency=afterok:$jid5 totalbed.sh --od ${output_dir} --lib ${lib_name} --fdr ${fdr} --num ${samp_size} | awk '{print $4}')
        else
            echo "Pass bedgraphs already made, line 277"
            cd ${chromsig_dir}
            jid6=$(sbatch totalbed.sh --od ${output_dir} --lib ${lib_name} --fdr ${fdr} --num ${samp_size} | awk '{print $4}')
        fi
    else
        echo "${lib_name}_total_FDR_${fdr}_pseudoRead_${samp_size}_pass_pileup.bed already exists"
    fi
    cd ${output_dir}
    pass_score_pref=$(basename ${total_pass} .bed)
    # If the scoreisland results are not present in the denoise directory, then run sicer on the denoising results.
    if [ ! -f "${pass_score_pref}-W200-G600.scoreisland" ]
    then
        echo "${pass_score_pref}-W200-G600.scoreisland does not exist"
        if [ ! -f ${total_pass} ]
        then
            #sicer -t ${total_pass} -s ${ref}
            cd ${chromsig_dir}
            jid7=$(sbatch --dependency=afterok:$jid6 sicer_run.sh --od ${output_dir} --ref ${ref} --bedf ${total_pass} | awk '{print $4}')
        else
            echo "Total pass already made, line 293"
            cd ${chromsig_dir}
            jid7=$(sbatch sicer_run.sh --od ${output_dir} --ref ${ref} --bedf ${total_pass} | awk '{print $4}')
        fi
    else
        echo "${pass_score_pref}-W200-G600.scoreisland already exists"
    fi
    # if the scoreisland results have not been created for the original bed file, then run sicer on the original file.
    cd ${directory}
    if [ ! -f "${bam_prefix}_bam2bed-W200-G600.scoreisland" ]
    then
        echo "${bam_prefix}_bam2bed-W200-G600.scoreisland does not exist"
        cd ${output_dir}
        if [ ! -f "${pass_score_pref}-W200-G600.scoreisland" ]
        then
            #sicer -t ${total_pass} -s ${ref}
            cd ${chromsig_dir}
            pwd
            echo "Creating total score island file"
            jid8=$(sbatch --dependency=afterok:$jid7 sicer_run.sh --od ${directory} --ref ${ref} --bedf ${std_bed} | awk '{print $4}')
        else
            cd ${chromsig_dir}
            echo "Pass scoreisland already made, line 311"
            jid8=$(sbatch sicer_run.sh --od ${directory} --ref ${ref} --bedf ${std_bed} | awk '{print $4}')
        fi
    else
        echo "${bam_prefix}_bam2bed-W200-G600.scoreisland already exists"
    fi
    
    # total all file results -- this will produce percentage and runtime stats for the run.
    cd ${output_dir}
    WAIT_INTERVAL=60

    until [ -f "$total_fail" ]; do
        echo "File not found yet, waiting..."
        sleep "$WAIT_INTERVAL"
    done

    cd ${directory}
    pass_lines=$(wc -l < ${output_dir}${total_pass})
    fail_lines=$(wc -l < ${output_dir}${total_fail})
    percentage=100
    echo "Read Type: "
    echo "${str_type}"
    echo "Number of Samples: "
    echo "${samp_size}"
    echo "Denoise Runtime: "
    cd ${output_dir}
    runtime_lines_array=()
    for file in *_logFile.txt; do runtime_lines_array+=( $(grep "Total: " "$file") ); done
    echo ${runtime_lines_array[@]}
    runtime_array=()
    for element in "${runtime_lines_array[@]}"; do substring=$(echo "$element" | cut -d' ' -f2); float_value=$(echo "$substring + 0" | bc -l); runtime_array+=("$float_value"); done
    echo ${runtime_array[@]}
    max_value=$(printf "%s\n" "${runtime_array[@]}" | awk 'BEGIN { max = -infinity } { if ($1 > max) max = $1 } END { print max }')
    echo "${max_value}"
    echo "Pass Percentage: "
    echo "scale=4 ; $pass_lines / ($pass_lines + $fail_lines) * $percentage" | bc
    echo "Pass Reads: "
    echo "${pass_lines}"
    echo "Fail Percentage: "
    echo "scale=4 ; $fail_lines / ($pass_lines + $fail_lines) * $percentage" | bc
    echo "Fail Reads: "
    echo "${fail_lines}"
    echo "Total Reads: "
    echo "$(wc -l < ${std_bed})"
fi

# report end of script.

cd ${chromsig_dir}

echo "Script completed in $(format_time $SECONDS)"

echo "[$(date)] DONE"
