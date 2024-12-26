input_file=$1

tag=${input_file/.bam/}
## see also: http://crazyhottommy.blogspot.com/2016/05/downsampling-for-bam-files-to-certain.html

if [ ! -f $input_file.bai ]; then
	samtools index $input_file
fi

samtools idxstats $input_file | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {print total}' &>$tag.nr_reads
