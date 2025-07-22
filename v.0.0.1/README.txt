Download the full folder
In the folder, creat a directory to store your bam file, and the sizes file corresponding to your reference genome
Install the required packages outlined in requirements.txt:
pip install -r v.0.0.1/requirements.txt

Run the following command to run the main program:
sbatch multiscript_bam2freq_enrich_bg.sh --bam <bam_file> --dir <directory_with_bam_file> --r <reference_genome> --fdr <fdr_value> --num <number_of_pseudo_samples> --type <Single_End_or_Paired_End_ex_SE_or_PE>
