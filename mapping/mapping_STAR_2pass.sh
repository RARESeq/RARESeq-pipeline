star_cmd=$1
input_folder=$2
second_pass_folder=$3
second_pass_genome=$4
nthreads=$5

mkdir -p $second_pass_folder
mkdir -p $second_pass_genome

#this command resets the shared memory. If STAR is interrupted before the genome is loaded, the subsequent runs will hang
#$star_cmd --genomeLoad Remove --genomeDir $first_pass_genome

for read1_file in $input_folder/*_R1.*
do
	echo $read1_file
	read2_file="${read1_file/_R1./_R2.}"
	echo $read2_file 
	tmp=$(basename $read1_file)
	output_file="${tmp/_R1.fastq/.}"
	output_file="${output_file/.gz/}"
	echo $output_file

	command=cat
	if [[ "$read1_file" == *.gz ]] 
	then
		echo "Using zcat"
		command=zcat
	fi

	$star_cmd --genomeDir $second_pass_genome \
	--readFilesIn $read1_file $read2_file --genomeLoad LoadAndKeep --readFilesCommand $command --runThreadN $nthreads --outFileNamePrefix $second_pass_folder/$output_file \
	--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 6000000000 \
	--quantMode TranscriptomeSAM \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
	--outFilterScoreMin 120 --outFilterMatchNmin 120 \
	--alignEndsProtrude 20 ConcordantPair
done