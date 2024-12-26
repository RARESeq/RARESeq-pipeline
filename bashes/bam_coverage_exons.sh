input_file=$1
exons=$2
output_dir=$3

mkdir -p $output_dir

file_name=$(basename $input_file)
tag=${file_name/.bam/}

samtools depth -a -b $exons $input_file | awk 'BEGIN{OFS="\t"}{a[$3]++;}END{for(k in a){print k, a[k];}}' >$output_dir/$tag.exons_cov

#bedtools coverage -abam $input_file  -b $exons -split -sorted -d |\
#awk 'BEGIN{OFS="\t"}{a[$5]++;}END{for(k in a){print k, a[k];}}' >$output_dir/$tag.exons_cov
