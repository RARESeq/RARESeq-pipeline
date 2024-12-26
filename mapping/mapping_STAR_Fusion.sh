input_folder=$1
genome_folder=$2
output_folder=$3
nthreads=$4

mkdir -p $output_folder

#this command resets the shared memory. If STAR is interrupted before the genome is loaded, the subsequent runs will hang
#$star_cmd --genomeLoad Remove --genomeDir $genome_folder

for read1_file in $input_folder/*_R1.*
do
	echo $read1_file
	read2_file="${read1_file/_R1./_R2.}"
	echo $read2_file 
	tmp=$(basename $read1_file)
	output_file="${tmp/_R1.fastq/}"
	output_file="${output_file/.gz/}"
	echo $output_file

	STAR-Fusion --genome_lib_dir $genome_folder \
	--left_fq $read1_file  --right_fq $read2_file \
	--min_FFPM 0 \
	--CPU $nthreads --output_dir $output_folder/$output_file \
	--FusionInspector inspect
done
