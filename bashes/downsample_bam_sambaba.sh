sambaba_cmd=$1
input_file=$2
count=$3
threads=$4
output_dir=$5

mkdir -p $output_dir

file_name=$(basename $input_file)
tag=${file_name/.bam/}
tag=${tag/.Aligned.sortedByCoord.out/}
tag=${tag/.AlignedByCoord.out/}
tag=${tag/.sorted/}
tag=${tag/.filt.sorted/}
tag=${tag/.Aligned.toTranscriptome.out/}
echo $tag

if [ -f $output_dir/$tag"__"$count.bam ]; then
	echo "File $output_dir/$tag__$count.bam exists."
	exit 1
fi

if [ -z $DOCKERRUN ]
then
fpath=${pipeline_dir}
else
fpath=""
fi

#samtools view $input_file | cut -f1 |  sort | uniq | sort -R > $output_dir/$tag.txt

#nreads=$(cat $read_list/$tag.txt | wc -l)
#echo $nreads
#echo $top

if [ ! -f $input_file.bai ]; then
	samtools index $input_file
fi


## see also: http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html
FACTOR=$(samtools idxstats $input_file | cut -f3 | awk -v COUNT=$count 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
  then 
  echo 'Requested number of reads ('$count') exceeds total read count in' $input_file'. Factor:' $FACTOR'.Using all the reads in the file.'
  FACTOR=1
fi

${fpath}${sambaba_cmd} view -s $FACTOR -t $threads --subsampling-seed=123 -f bam -l 5 $input_file > $output_dir/$tag"__"$count.bam

samtools sort $output_dir/$tag"__"$count.bam -o $output_dir/$tag"__"$count.sorted.bam
samtools index $output_dir/$tag"__"$count.sorted.bam

rm $output_dir/$tag"__"$count.bam





