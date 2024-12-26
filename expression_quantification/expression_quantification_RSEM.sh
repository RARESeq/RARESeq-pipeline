input_file=$1
genome_folder=$2
output_folder=$3
n_threads=$4
mkdir -p $output_folder
 
echo $input_file
filename=$(basename $input_file)
output_file="${filename/.Aligned.toTranscriptome.out.bam/}"
echo $output_file

rsem-calculate-expression -p $n_threads --bam --paired-end $input_file $genome_folder $output_folder/$output_file
