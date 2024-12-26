
rnaseqc_command=$1
gtf=$2
bam=$3
output_dir=$4 

if [ -z $DOCKERRUN ]
then
fpath=${pipeline_dir}
else
fpath=""
fi

${fpath}${rnaseqc_command} $gtf $bam $output_dir



