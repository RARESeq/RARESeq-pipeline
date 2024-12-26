file=$1
output_folder=$2
mkdir -p $output_folder

# Extract the base name and tag from the input file
file_name=$(basename $file)
tag=${file_name/.Aligned.sortedByCoord.out.bam/}

# Sort the BAM file by read name
samtools sort -n $file -o $output_folder/$tag.sortedByName.bam

# Extract the header
samtools view -H $output_folder/$tag.sortedByName.bam > $output_folder/$tag.header.txt

# Extract uniquely mapping pairs without saving read names in a file
samtools view -q 255 $output_folder/$tag.sortedByName.bam |\
awk '{
  if ($1 == prev_name) {
    if (count == 1) {
      print prev_line;
      print;
    }
  } else {
    count = 0;
  }
  prev_name = $1;
  prev_line = $0;
  count++;
}' > $output_folder/$tag.uniquePairs.txt

# Combine the header and unique pairs into a BAM file
cat $output_folder/$tag.header.txt $output_folder/$tag.uniquePairs.txt |\
samtools view -Sbh - -o $output_folder/$tag.int.bam

# Sort the BAM file by coordinate
samtools sort $output_folder/$tag.int.bam -o $output_folder/$tag.bam

# Index the resulting BAM file
samtools index $output_folder/$tag.bam

# Remove temporary files
rm -f $output_folder/$tag.int.bam $output_folder/$tag.sortedByName.bam $output_folder/$tag.header.txt $output_folder/$tag.uniquePairs.txt
