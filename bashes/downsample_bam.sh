input_file=$1
read_list=$2
output_dir=$3
count=$4

mkdir -p $output_dir

file_name=$(basename $input_file)
tag=${file_name/.bam/}
tag=${tag/.Aligned.sortedByCoord.out/}
tag=${tag/.AlignedByCoord.out/}
tag=${tag/.sorted/}
tag=${tag/.Aligned.toTranscriptome.out/}
echo $tag

if [ -f $output_dir/$tag"__"$count.bam ]; then
	echo "File $output_dir/$tag__$count.bam exists."
	exit 1
fi

#samtools view $input_file | cut -f1 |  sort | uniq | sort -R > $output_dir/$tag.txt

#nreads=$(cat $read_list/$tag.txt | wc -l)
#echo $nreads
#echo $top
(samtools view -H $input_file; samtools view $input_file |
	awk -v top="$count" '{if(FNR==NR){ 
		if(NR<=top)
		{
			#print NR, top, $1;
			a[$1]=0;
		}
		next;
	} 
	if($1 in a){print}
	}' $read_list -) | \
	samtools view -Sbh - -o $output_dir/$tag"__"$count.bam

samtools sort $output_dir/$tag"__"$count.bam -o $output_dir/$tag"__"$count.sorted.bam
samtools index $output_dir/$tag"__"$count.sorted.bam
rm $output_dir/$tag"__"$count.bam
