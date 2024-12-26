file=$1
output_dir=$2
mkdir -p $output_dir
file_name=$(basename $file)
tag=${file_name/.bam/}
echo $tag

samtools view $file | cut -f1 |  sort | uniq | sort -R > $output_dir/$tag.read_ids
cat $output_dir/$tag.read_ids | wc -l >$output_dir/$tag.count

