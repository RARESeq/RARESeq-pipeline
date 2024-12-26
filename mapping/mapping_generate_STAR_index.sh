star_cmd=$1
rsem_ref_command=$2
genome_dir=$3
fasta_file=$4
gtf_file=$5
threads=$6

mkdir -p $genome_dir
mkdir -p $genome_dir/RSEM

$star_cmd --runThreadN $threads --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $fasta_file --sjdbGTFfile $gtf_file --outFileNamePrefix $genome_dir/out

$rsem_ref_command --gtf $gtf_file -p $threads $fasta_file $genome_dir/RSEM
  

