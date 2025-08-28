#!/bin/bash
#SBATCH --job-name=freq_enrich_python
#SBATCH --account=minjilab0
#SBATCH --array=0-23 #24 cores for one bam file
#SBATCH --mem=50g
#SBATCH --time=00-24:00:00
#SBATCH --output=/nfs/turbo/umms-minjilab/njgupta/chromsig/denoise_output/slurm-%A_%a.out
#SBATCH --error=/nfs/turbo/umms-minjilab/njgupta/chromsig/denoise_outerr/slurm-%A_%a.out
#SBATCH --mail-user=njgupta@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --profile=Task

## The help message:
function usage
{
    echo -e "DESCRIPTION:
    This script runs the python program which performs the denoising process on the input dataset. It processes one chromosome at a time, performing denoising on all 
    chromosomes in parallel, and produces a logfile, master file, and pass and fail pileup beds for each chromosome.
    USAGE:
    sbatch denoise.sh --cd chromsig_dir --bg std_bg --bed std_bed --dir directory --lib lib_name --r ref --fdr fdr --num samp_size --type str_type --plot plot_hist
    
      ARGUMENTS:
        --cd    Directory of overall (multiscript_bam2freq_enrich_bg.sh) program.
        --bg    Standardized bedgraph file (standardized means all non-standard chromosomes have been removed, only chromosomes in .chrom.sizes file remain).
        --bed   Standardized bed file (standardized means all non-standard chromosomes have been removed, only chromosomes in .chrom.sizes file remain).
        --dir   Directory containing .bed and .bedgraph files (ex: /Users/chromsig/Data/input_file_name/).
        --lib   Name of dataset library.
        --r     Name of reference genome (ex: hg38).
        --fdr   FDR threshold value for denoising process.
        --num   Number of pseudo samples to be generated for denoising process.
        --type  Type of read in your dataset, Single End or Paired End (ex: SE or PE).
        --plot  True or False for whether to generate histogram plots for distances between read fragments.

    The output master.txt file for each chromosome has the following columns:
    READ_ID:      The ID of the read, taken from the fourth column of the input bed file. Paired fragments share an ID.
    Start Coord:  The starting index of the observed read along the chromosome.
    End Coord:    The ending index of the observed read along the chromosome.
    Category:     The source of the read (ex: For reads taken from the bed file, this would be 'Orig').
    Obs:          The enrichment value of the observed read.
    rawpval1:     The raw P-value calculated by comparing the observed read's enrichment to the enrichments of the corresponding read's pseudo samples.
    adjpval1:     The raw P-value adjusted based on the Benjamin-Hochberg FDR threshold.
    decis1:       Determination for whether the read in question passed or failed (ex: 'PASS' or 'FAIL')."
}

## Parse arguments from the command line
while [ "$1" != "" ]; do
    case $1 in
        -s | --cd )          shift
                                chromsig_dir=$1
                                ;;
        -s | --bg )          shift
                                std_bg=$1
                                ;;
        -s | --bed )         shift
                                std_bed=$1
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
        -s | --plot )        shift
                                plot_hist=$1
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

# determining amount of coverage

bedtools merge -i ${std_bed} > "${lib_name}_merge_result.bed"

bedtools summary -i "${lib_name}_merge_result.bed" -g hg38.chrom.sizes | column -t > "${lib_name}_coverage_summary.txt"


awk 'NR<=1 {next} { 
    if ($2 != 0) { 
        print $1, ($4 / $2) 
    } else { 
        print 0 
    } 
}' "${lib_name}_coverage_summary.txt" > "${lib_name}_final_coverage.txt"

coverage_arr=$(cut -f1 ${lib_name}_final_coverage.txt)
covlst=()
while read -r line; do covlst+=("$line"); done <<< "$coverage_arr"

chroms=$(cut -f1 ${ref}.chrom.sizes)

chromtoremove="chrY"

chromslst=()
while read -r line; do chromslst+=("$line"); done <<< "$chroms"
finalchroms=()
for item in "${!chromslst[@]}"; do
  echo "Current chrom: ${chromslst[$item]}"
  # If the current item is not the one to remove, add it to the filtered list
  if [[ "${chromslst[$item]}" != "$chromtoremove" ]]
  then
    if grep -qFw "${chromslst[$item]}" "$std_bed"
    then
      echo ${covlst[$item]}
      num_to_compare=$(echo ${covlst[$item]} | cut -d ' ' -f2)
      echo "Coverage: $num_to_compare"
      if (( $(echo "$num_to_compare > 0.20" | bc -l) )); then
        finalchroms+=("${chromslst[$item]}")
      fi
    fi
  fi
done

echo "Chroms being processed:"
echo ${finalchroms[@]}

#chromslst=(${chroms//$'\n'/ })
#IFS=$' ' read -r -a chromslst2 <<< "$chromslst"
if (( $SLURM_ARRAY_TASK_ID < ${#finalchroms[@]} )); then
  current_chrom=${finalchroms[$SLURM_ARRAY_TASK_ID]}
  echo $current_chrom

  cd ${chromsig_dir}
  echo ${std_bed}
  python chromsig_Freq_Enrich_Test.py "${lib_name} ${ref} ${fdr} ${current_chrom} ${samp_size} ${std_bg} ${directory} ${std_bed} ${str_type} ${plot_hist}"
fi
