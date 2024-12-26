star_cmd=$1
input_folder=$2
first_pass_folder=$3
second_pass_genome=$4
genome_gtf=$5
genome_fa=$6
nthreads=$7

mkdir -p $first_pass_folder
mkdir -p $second_pass_genome

#this command resets the shared memory. If STAR is interrupted before the genome is loaded, the subsequent runs will hang
#$star_cmd --genomeLoad Remove --genomeDir $first_pass_genome

cat $first_pass_folder/*SJ.out.tab | awk '($5 > 0 && $7 > 5 && $6==0)' | cut -f1-6 | sort | uniq  > $second_pass_genome/SJ.filt.tab

$star_cmd --runThreadN $nthreads --runMode genomeGenerate --outFileNamePrefix $second_pass_genome/out \
	--genomeDir $second_pass_genome --genomeFastaFiles $genome_fa \
	--sjdbGTFfile $genome_gtf \
	--sjdbFileChrStartEnd $second_pass_genome/SJ.filt.tab
 