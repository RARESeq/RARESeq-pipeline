input_file=$1
read_list=$2
output_dir=$3

mkdir -p $output_dir

file_name=$(basename $read_list)
tag=${file_name/.sorted.dualindex-deduped.sorted.bam.read_ids/}
echo $tag

(samtools view -H $input_file; samtools view $input_file |
	awk '{
	
	if(FNR==NR){ 
		#print $1; 
		a[$1]=0; 
		next;
	} 
	if($1 in a){print}
	}' $read_list -) | \
	samtools view -Sbh - -o $output_dir/$tag.Aligned.toTranscriptome.out.bam
