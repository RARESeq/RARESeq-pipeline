file=$1
fa_index=$2
output_folder=$3

mkdir -p $output_folder

filename=$(basename $file)
tag=${filename/.bam/}

folder=$(dirname $0)

echo "Preparing file $filename for deduping"

samtools sort -n -@ 1 -o $output_folder/$tag.int.bam $file 

python $folder/dedupe_prepare_bam.py $output_folder/$tag.int.bam $fa_index $output_folder 

samtools sort -@ 1 -o $output_folder/$tag.sorted.bam $output_folder/$tag.int.bam
samtools index $output_folder/$tag.sorted.bam
rm $output_folder/$tag.int.bam
	



